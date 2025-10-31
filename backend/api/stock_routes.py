"""
Stock Management API Routes

Provides RESTful endpoints for managing custom stock databases
for retrosynthesis planning.
"""

from fastapi import APIRouter, HTTPException, UploadFile, File
from pydantic import BaseModel, Field
from typing import List, Dict, Optional
import logging

from adapters.aizynthfinder.stock_manager import StockManager

router = APIRouter(prefix="/api/stock", tags=["stock"])
logger = logging.getLogger(__name__)


class CompoundRequest(BaseModel):
    """Request model for adding/removing compounds"""
    smiles: str = Field(..., description="SMILES string of the compound")
    name: Optional[str] = Field(None, description="Name of the compound")


class StockInfo(BaseModel):
    """Stock database information"""
    total_compounds: int
    csv_path: str
    hdf5_path: str
    csv_exists: bool
    hdf5_exists: bool


class CompoundInfo(BaseModel):
    """Compound information"""
    smiles: str
    name: str


@router.get("/info", response_model=StockInfo)
async def get_stock_info():
    """
    Get information about the current stock database.

    Returns:
        Stock database statistics and file paths
    """
    try:
        manager = StockManager()
        info = manager.get_stock_info()

        return StockInfo(
            total_compounds=info["total_compounds"],
            csv_path=info["csv_path"],
            hdf5_path=info["hdf5_path"],
            csv_exists=info["csv_exists"],
            hdf5_exists=info["hdf5_exists"]
        )

    except Exception as e:
        logger.error(f"Failed to get stock info: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/compounds", response_model=List[CompoundInfo])
async def list_compounds(limit: int = 100, offset: int = 0):
    """
    List compounds in the stock database.

    Args:
        limit: Maximum number of compounds to return
        offset: Number of compounds to skip

    Returns:
        List of compounds with SMILES and names
    """
    try:
        manager = StockManager()
        info = manager.get_stock_info()
        compounds = info["compounds"]

        # Paginate
        items = list(compounds.items())[offset:offset + limit]

        return [
            CompoundInfo(smiles=smiles, name=name)
            for smiles, name in items
        ]

    except Exception as e:
        logger.error(f"Failed to list compounds: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/compounds/add")
async def add_compound(request: CompoundRequest):
    """
    Add a compound to the stock database.

    Args:
        request: Compound SMILES and name

    Returns:
        Success message
    """
    try:
        manager = StockManager()
        name = request.name or "unnamed"

        success = manager.add_compound(request.smiles, name)

        if not success:
            raise HTTPException(status_code=400, detail="Invalid SMILES or failed to add compound")

        return {
            "message": f"Successfully added {name}",
            "smiles": request.smiles,
            "name": name
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to add compound: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/compounds/remove")
async def remove_compound(request: CompoundRequest):
    """
    Remove a compound from the stock database.

    Args:
        request: Compound SMILES

    Returns:
        Success message
    """
    try:
        manager = StockManager()

        success = manager.remove_compound(request.smiles)

        if not success:
            raise HTTPException(status_code=404, detail="Compound not found or invalid SMILES")

        return {
            "message": f"Successfully removed compound",
            "smiles": request.smiles
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to remove compound: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/initialize")
async def initialize_default_stock():
    """
    Initialize default lab stock database with common reagents.

    Returns:
        Success message and compound count
    """
    try:
        manager = StockManager()

        success = manager.initialize_default_stock()

        if not success:
            raise HTTPException(status_code=500, detail="Failed to initialize stock database")

        info = manager.get_stock_info()

        return {
            "message": "Default stock database initialized",
            "total_compounds": info["total_compounds"],
            "csv_path": info["csv_path"],
            "hdf5_path": info["hdf5_path"]
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to initialize stock: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/update")
async def update_stock_database():
    """
    Update HDF5 stock database from CSV file.

    Regenerates the HDF5 file used by AiZynthFinder from the editable CSV file.

    Returns:
        Success message
    """
    try:
        manager = StockManager()

        success = manager.update_stock_database()

        if not success:
            raise HTTPException(status_code=500, detail="Failed to update stock database")

        info = manager.get_stock_info()

        return {
            "message": "Stock database updated from CSV",
            "total_compounds": info["total_compounds"]
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to update stock: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/upload-csv")
async def upload_stock_csv(file: UploadFile = File(...)):
    """
    Upload a CSV file with stock compounds.

    CSV format: smiles,name (one compound per line)

    Args:
        file: CSV file upload

    Returns:
        Success message and compound count
    """
    try:
        # Read uploaded file
        contents = await file.read()

        # Save to custom_stock.csv
        manager = StockManager()
        csv_path = manager.stock_csv

        with open(csv_path, 'wb') as f:
            f.write(contents)

        # Update HDF5 from new CSV
        success = manager.update_stock_database()

        if not success:
            raise HTTPException(status_code=500, detail="Failed to process uploaded CSV")

        info = manager.get_stock_info()

        return {
            "message": "CSV uploaded and stock database updated",
            "total_compounds": info["total_compounds"],
            "filename": file.filename
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to upload CSV: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/download-csv")
async def download_stock_csv():
    """
    Download the current stock CSV file.

    Returns:
        CSV file with all stock compounds
    """
    try:
        from fastapi.responses import FileResponse

        manager = StockManager()
        csv_path = manager.stock_csv

        if not csv_path.exists():
            raise HTTPException(status_code=404, detail="Stock CSV file not found")

        return FileResponse(
            path=str(csv_path),
            filename="custom_stock.csv",
            media_type="text/csv"
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Failed to download CSV: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/validate-smiles")
async def validate_smiles(request: CompoundRequest):
    """
    Validate a SMILES string without adding it to the database.

    Args:
        request: SMILES string to validate

    Returns:
        Validation result and canonical form
    """
    try:
        manager = StockManager()

        canonical = manager.canonicalize_smiles(request.smiles)

        if canonical:
            return {
                "valid": True,
                "original": request.smiles,
                "canonical": canonical
            }
        else:
            return {
                "valid": False,
                "original": request.smiles,
                "error": "Invalid SMILES string"
            }

    except Exception as e:
        logger.error(f"Failed to validate SMILES: {e}")
        raise HTTPException(status_code=500, detail=str(e))
