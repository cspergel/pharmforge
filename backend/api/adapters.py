"""
Adapters API endpoint
Allows querying and testing adapters
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Any, Optional
import logging

from backend.core.adapters.protocol import registry

logger = logging.getLogger(__name__)
router = APIRouter()


class AdapterTestRequest(BaseModel):
    """Request model for testing an adapter"""
    adapter_name: str
    smiles: str


class AdapterTestResponse(BaseModel):
    """Response model for adapter test"""
    success: bool
    data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


@router.get("/adapters")
async def list_adapters():
    """
    List all registered adapters

    Returns:
        List of adapter metadata
    """
    try:
        adapters = registry.get_all_metadata()
        return {
            "adapters": adapters,
            "count": len(adapters)
        }
    except Exception as e:
        logger.error(f"Error listing adapters: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/adapters/{adapter_name}")
async def get_adapter_info(adapter_name: str):
    """
    Get information about a specific adapter

    Args:
        adapter_name: Name of the adapter

    Returns:
        Adapter metadata
    """
    adapter = registry.get(adapter_name)

    if not adapter:
        raise HTTPException(
            status_code=404,
            detail=f"Adapter '{adapter_name}' not found"
        )

    return adapter.get_metadata()


@router.post("/adapters/test", response_model=AdapterTestResponse)
async def test_adapter(request: AdapterTestRequest):
    """
    Test an adapter with a SMILES string

    Args:
        request: AdapterTestRequest with adapter_name and smiles

    Returns:
        AdapterTestResponse with results
    """
    adapter = registry.get(request.adapter_name)

    if not adapter:
        raise HTTPException(
            status_code=404,
            detail=f"Adapter '{request.adapter_name}' not found"
        )

    try:
        logger.info(f"Testing adapter {request.adapter_name} with SMILES: {request.smiles}")

        # Execute the adapter
        result = await adapter.execute(request.smiles)

        return AdapterTestResponse(
            success=result.success,
            data=result.data,
            error=result.error,
            metadata=result.metadata
        )

    except Exception as e:
        logger.error(f"Error testing adapter {request.adapter_name}: {e}")
        raise HTTPException(status_code=500, detail=str(e))
