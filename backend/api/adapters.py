"""
Adapters API Router - Backend API endpoints for PharmForge adapter system

Provides endpoints to:
- List all available adapters with metadata
- Get detailed info for specific adapters
- Execute single adapters on input data
- Execute multiple adapters in batch mode
"""

from fastapi import APIRouter, HTTPException, status
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional
import logging
import time
import asyncio

from backend.core.adapters.protocol import registry, AdapterResult

logger = logging.getLogger(__name__)
router = APIRouter()


# ============================================================================
# Pydantic Models for Request/Response
# ============================================================================


class AdapterMetadata(BaseModel):
    """Metadata for a single adapter"""
    name: str
    display_name: str
    adapter_type: str
    category: str
    status: str
    description: str
    version: str
    enabled: bool
    benchmarks: Optional[Dict[str, Any]] = None
    config: Optional[Dict[str, Any]] = None


class AdapterListResponse(BaseModel):
    """Response for GET /api/adapters/list"""
    adapters: List[AdapterMetadata]
    total_count: int
    categories: Dict[str, int]


class AdapterInfoResponse(BaseModel):
    """Response for GET /api/adapters/{adapter_name}/info"""
    name: str
    display_name: str
    adapter_type: str
    category: str
    description: str
    version: str
    enabled: bool
    required_inputs: List[str]
    optional_params: Dict[str, Any]
    output_schema: Dict[str, Any]
    benchmarks: Optional[Dict[str, Any]] = None
    config: Optional[Dict[str, Any]] = None


class ExecuteAdapterRequest(BaseModel):
    """Request for POST /api/adapters/{adapter_name}/execute"""
    input_data: Dict[str, Any] = Field(..., description="Input data for adapter (e.g., {smiles: 'CCO'})")
    params: Optional[Dict[str, Any]] = Field(default={}, description="Optional adapter-specific parameters")
    use_cache: bool = Field(default=True, description="Whether to use cache for this execution")


class ExecuteAdapterResponse(BaseModel):
    """Response for POST /api/adapters/{adapter_name}/execute"""
    success: bool
    data: Optional[Any] = None
    error: Optional[str] = None
    execution_time: float
    cache_hit: bool = False
    metadata: Optional[Dict[str, Any]] = None


class BatchExecuteRequest(BaseModel):
    """Request for POST /api/adapters/batch"""
    input_data: Dict[str, Any] = Field(..., description="Input data to pass to all adapters")
    adapter_names: List[str] = Field(..., description="List of adapter names to execute")
    params: Optional[Dict[str, Any]] = Field(default={}, description="Parameters to pass to all adapters")
    use_cache: bool = Field(default=True, description="Whether to use cache")


class BatchAdapterResult(BaseModel):
    """Result for a single adapter in batch execution"""
    adapter: str
    success: bool
    data: Optional[Any] = None
    error: Optional[str] = None
    execution_time: float
    cache_hit: bool = False


class BatchExecuteResponse(BaseModel):
    """Response for POST /api/adapters/batch"""
    results: List[BatchAdapterResult]
    total_time: float
    successful_count: int
    failed_count: int


# ============================================================================
# Helper Functions
# ============================================================================


def _categorize_adapter(adapter_name: str, adapter_type: str) -> str:
    """
    Categorize adapter based on name and type

    Args:
        adapter_name: Name of the adapter
        adapter_type: Type of adapter (api, local, ml)

    Returns:
        Category string
    """
    name_lower = adapter_name.lower()

    # Chemical databases
    if any(x in name_lower for x in ['pubchem', 'chembl', 'zinc', 'drugcentral', 'surechembl', 'coconut']):
        return "Chemical Databases"

    # Target & Disease
    if any(x in name_lower for x in ['opentargets', 'disgenet', 'uniprot', 'string', 'biogrid', 'intact']):
        return "Target & Disease"

    # Protein structures
    if any(x in name_lower for x in ['alphafold', 'pdb', 'swiss', 'rcsb', 'sabdab', 'immunebuilder']):
        return "Protein Structure"

    # Docking & Binding
    if any(x in name_lower for x in ['vina', 'gnina', 'diffdock', 'bindingdb', 'openmm']):
        return "Docking & Binding"

    # ADMET & Properties
    if any(x in name_lower for x in ['admet', 'pkcsm', 'rdkit', 'swisstarget', 'targetnet', 'tox21', 'comptox', 'xtb']):
        return "ADMET & Properties"

    # Retrosynthesis & De Novo
    if any(x in name_lower for x in ['aizynthfinder', 'retrosynthesis', 'reinvent', 'molgan', 'denovo', 'llm', 'ord']):
        return "Retrosynthesis & De Novo"

    # Clinical & Safety
    if any(x in name_lower for x in ['clinical', 'faers', 'fda']):
        return "Clinical & Safety"

    # Literature & Patents
    if any(x in name_lower for x in ['pubmed', 'europepmc', 'patent', 'lens', 'google']):
        return "Literature & Patents"

    # Genomics & Expression
    if any(x in name_lower for x in ['geo', 'gtex', 'reactome', 'kegg']):
        return "Genomics & Expression"

    # Metabolomics
    if any(x in name_lower for x in ['hmdb', 'metabol']):
        return "Metabolomics"

    # RNA Databases
    if any(x in name_lower for x in ['rna', 'mirna']):
        return "RNA Databases"

    return "Other"


def _generate_display_name(adapter_name: str) -> str:
    """Generate a human-readable display name from adapter name"""
    # Handle special cases
    special_names = {
        'pubchem': 'PubChem',
        'chembl': 'ChEMBL',
        'rdkit': 'RDKit',
        'admet': 'ADMET',
        'ai': 'AI',
        'alphafold': 'AlphaFold',
        'rcsb_pdb': 'RCSB PDB',
        'pdb': 'PDB',
        'gnina': 'GNINA',
        'vina': 'AutoDock Vina',
        'opentargets': 'Open Targets',
        'disgenet': 'DisGeNET',
        'uniprot': 'UniProt',
        'drugcentral': 'DrugCentral',
        'zinc': 'ZINC',
        'surechembl': 'SureChEMBL',
        'fda_faers': 'FDA FAERS',
        'gtex': 'GTEx',
        'geo': 'GEO',
        'kegg': 'KEGG',
        'pubmed': 'PubMed',
        'europepmc': 'Europe PMC',
        'diffdock': 'DiffDock',
        'swissmodel': 'SWISS-MODEL',
        'swisstarget': 'SwissTargetPrediction',
        'pkcsm': 'pkCSM',
        'aizynthfinder': 'AiZynthFinder',
        'llm_retrosynthesis': 'LLM Retrosynthesis',
        'targetnet': 'TargetNet',
        'openmm': 'OpenMM',
        'molgan': 'MolGAN',
        'reinvent': 'REINVENT',
        'coconut': 'COCONUT',
        'hmdb': 'HMDB',
        'tox21': 'Tox21',
        'comptox': 'CompTox',
        'sabdab': 'SAbDab',
        'immunebuilder': 'ImmuneBuilder',
        'rnacentral': 'RNAcentral',
        'ord': 'ORD',
        'intact': 'IntAct',
        'xtb': 'xTB',
    }

    if adapter_name in special_names:
        return special_names[adapter_name]

    # Default: capitalize and replace underscores
    return adapter_name.replace('_', ' ').title()


def _generate_description(adapter_name: str, category: str) -> str:
    """Generate a description for the adapter"""
    descriptions = {
        'pubchem': 'Fetch molecular properties and compound information from PubChem database',
        'chembl': 'Query bioactivity data and drug-target interactions from ChEMBL',
        'rdkit': 'Calculate molecular descriptors and properties using RDKit',
        'admet_ai': 'Predict ADMET properties using AI models',
        'vina': 'Perform molecular docking with AutoDock Vina',
        'aizynthfinder': 'AI-powered retrosynthetic route planning',
        'alphafold': 'Access AlphaFold protein structure predictions',
        'opentargets': 'Query disease-target associations from Open Targets',
        'disgenet': 'Access gene-disease associations from DisGeNET',
        'uniprot': 'Retrieve protein information and sequences from UniProt',
        'drugcentral': 'Query approved drug information and targets',
        'zinc': 'Search ZINC database for purchasable compounds',
        'gnina': 'CNN-enhanced molecular docking with GNINA',
        'bindingdb': 'Query experimental binding affinity data',
        'surechembl': 'Search patent chemistry from SureChEMBL',
        'rcsb_pdb': 'Access experimental protein structures from RCSB PDB',
        'swissmodel': 'Access homology models from SWISS-MODEL',
        'diffdock': 'Deep learning-based molecular docking',
        'pkcsm': 'Predict pharmacokinetic properties',
        'fda_faers': 'Query FDA adverse event reports',
        'pubmed': 'Search biomedical literature from PubMed',
        'clinicaltrials': 'Query clinical trial data',
        'gtex': 'Access gene expression data from GTEx',
        'kegg': 'Query KEGG pathway database',
        'reactome': 'Access Reactome pathway data',
        'coconut': 'Search 400k+ natural products from COCONUT database',
        'hmdb': 'Query human metabolome data and metabolite information',
        'tox21': 'Predict toxicity using Tox21 assay models',
        'comptox': 'Access EPA CompTox chemistry and toxicity data',
        'sabdab': 'Query antibody structures from SAbDab database',
        'immunebuilder': 'Fast antibody structure prediction (Fab, scFv, VHH)',
        'rnacentral': 'Search RNA sequences for RNA-targeting therapeutics',
        'ord': 'Query chemical reactions from Open Reaction Database',
        'intact': 'Access experimental protein-protein interactions',
        'xtb': 'Fast quantum mechanical calculations (GFN methods)',
    }

    return descriptions.get(adapter_name, f"{category} adapter for drug discovery workflows")


def _get_adapter_metadata(adapter) -> AdapterMetadata:
    """Extract metadata from an adapter instance"""
    adapter_name = adapter.name
    adapter_type = adapter.adapter_type
    category = _categorize_adapter(adapter_name, adapter_type)

    return AdapterMetadata(
        name=adapter_name,
        display_name=_generate_display_name(adapter_name),
        adapter_type=adapter_type,
        category=category,
        status="active" if adapter.enabled else "inactive",
        description=_generate_description(adapter_name, category),
        version=getattr(adapter, 'version', '1.0.0'),
        enabled=adapter.enabled,
        benchmarks=None,  # TODO: Add benchmark data if available
        config=adapter.config if hasattr(adapter, 'config') else None
    )


def _get_detailed_adapter_info(adapter) -> AdapterInfoResponse:
    """Get detailed information about an adapter"""
    metadata = _get_adapter_metadata(adapter)

    # Extract required inputs and optional params from adapter
    # This is a simplified version - you may want to enhance this
    required_inputs = ["smiles"]  # Default, most adapters take SMILES

    # Categorize by adapter type for more specific requirements
    if metadata.category == "Target & Disease":
        required_inputs = ["gene_name"]
    elif metadata.category == "Protein Structure":
        required_inputs = ["protein_id"]
    elif metadata.category == "Literature & Patents":
        required_inputs = ["query"]

    optional_params = {
        "use_cache": "bool - Whether to use cached results (default: true)",
        "timeout": "int - Request timeout in seconds",
    }

    output_schema = {
        "success": "bool - Whether execution succeeded",
        "data": "object - Adapter-specific output data",
        "error": "string - Error message if failed",
        "metadata": "object - Additional execution metadata"
    }

    return AdapterInfoResponse(
        name=metadata.name,
        display_name=metadata.display_name,
        adapter_type=metadata.adapter_type,
        category=metadata.category,
        description=metadata.description,
        version=metadata.version,
        enabled=metadata.enabled,
        required_inputs=required_inputs,
        optional_params=optional_params,
        output_schema=output_schema,
        benchmarks=metadata.benchmarks,
        config=metadata.config
    )


# ============================================================================
# API Endpoints
# ============================================================================


@router.get("/list", response_model=AdapterListResponse)
async def list_adapters():
    """
    List all registered adapters with metadata

    Returns comprehensive metadata for all 75+ adapters including:
    - Name and display name
    - Category and type
    - Status and version
    - Description

    Example:
        GET /api/adapters/list
    """
    try:
        # Get all registered adapters
        adapter_names = registry.list_adapters()

        adapters_metadata = []
        categories_count = {}

        for name in adapter_names:
            adapter = registry.get(name)
            if adapter:
                metadata = _get_adapter_metadata(adapter)
                adapters_metadata.append(metadata)

                # Count categories
                categories_count[metadata.category] = categories_count.get(metadata.category, 0) + 1

        logger.info(f"Listed {len(adapters_metadata)} adapters across {len(categories_count)} categories")

        return AdapterListResponse(
            adapters=adapters_metadata,
            total_count=len(adapters_metadata),
            categories=categories_count
        )

    except Exception as e:
        logger.error(f"Error listing adapters: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to list adapters: {str(e)}"
        )


@router.get("/{adapter_name}/info", response_model=AdapterInfoResponse)
async def get_adapter_info(adapter_name: str):
    """
    Get detailed information about a specific adapter

    Args:
        adapter_name: Name of the adapter (e.g., "pubchem", "vina_docking")

    Returns:
        Detailed adapter information including required inputs, parameters, and output schema

    Example:
        GET /api/adapters/pubchem/info
    """
    try:
        adapter = registry.get(adapter_name)

        if not adapter:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Adapter '{adapter_name}' not found. Available adapters: {', '.join(registry.list_adapters())}"
            )

        info = _get_detailed_adapter_info(adapter)
        logger.info(f"Retrieved info for adapter: {adapter_name}")

        return info

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error getting adapter info for {adapter_name}: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to get adapter info: {str(e)}"
        )


@router.post("/{adapter_name}/execute", response_model=ExecuteAdapterResponse)
async def execute_adapter(adapter_name: str, request: ExecuteAdapterRequest):
    """
    Execute a single adapter on input data

    Args:
        adapter_name: Name of the adapter to execute
        request: Execution request containing input_data and optional params

    Returns:
        Execution results including success status, data, and timing

    Example:
        POST /api/adapters/pubchem/execute
        {
            "input_data": {"smiles": "CCO"},
            "params": {},
            "use_cache": true
        }
    """
    start_time = time.time()

    try:
        # Get adapter
        adapter = registry.get(adapter_name)

        if not adapter:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Adapter '{adapter_name}' not found"
            )

        if not adapter.enabled:
            raise HTTPException(
                status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                detail=f"Adapter '{adapter_name}' is currently disabled"
            )

        # Extract input
        input_value = request.input_data.get('smiles') or request.input_data.get('input') or request.input_data

        # Extract protein_data if provided (for docking adapters)
        protein_data = request.input_data.get('protein_data')

        logger.info(f"Executing adapter {adapter_name} with input: {input_value}")
        if protein_data:
            logger.info(f"Protein data provided: {protein_data.get('type')} = {protein_data.get('value')}")

        # Merge protein_data into params if provided
        adapter_params = dict(request.params) if request.params else {}
        if protein_data:
            adapter_params['protein_data'] = protein_data

        # Execute adapter with cache support
        result: AdapterResult = await adapter(
            input_value,
            use_cache=request.use_cache,
            **adapter_params
        )

        execution_time = time.time() - start_time

        logger.info(f"Adapter {adapter_name} execution completed in {execution_time:.3f}s (cache_hit={result.cache_hit})")

        return ExecuteAdapterResponse(
            success=result.success,
            data=result.data,
            error=result.error,
            execution_time=execution_time,
            cache_hit=result.cache_hit,
            metadata=result.metadata
        )

    except HTTPException:
        raise
    except Exception as e:
        execution_time = time.time() - start_time
        logger.error(f"Error executing adapter {adapter_name}: {e}", exc_info=True)

        return ExecuteAdapterResponse(
            success=False,
            data=None,
            error=str(e),
            execution_time=execution_time,
            cache_hit=False,
            metadata={"error_type": type(e).__name__}
        )


@router.post("/batch", response_model=BatchExecuteResponse)
async def batch_execute_adapters(request: BatchExecuteRequest):
    """
    Execute multiple adapters on the same input data

    This endpoint allows you to run multiple adapters in parallel on the same input,
    which is useful for workflows that require data from multiple sources.

    Args:
        request: Batch execution request with input_data and list of adapter_names

    Returns:
        Results from all adapter executions with timing information

    Example:
        POST /api/adapters/batch
        {
            "input_data": {"smiles": "CCO"},
            "adapter_names": ["pubchem", "chembl", "rdkit_local"],
            "use_cache": true
        }
    """
    start_time = time.time()

    try:
        # Validate adapters exist
        missing_adapters = []
        for name in request.adapter_names:
            if not registry.get(name):
                missing_adapters.append(name)

        if missing_adapters:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Adapters not found: {', '.join(missing_adapters)}"
            )

        logger.info(f"Batch executing {len(request.adapter_names)} adapters")

        # Execute all adapters in parallel
        async def execute_single(adapter_name: str) -> BatchAdapterResult:
            adapter_start = time.time()
            try:
                adapter = registry.get(adapter_name)

                if not adapter or not adapter.enabled:
                    return BatchAdapterResult(
                        adapter=adapter_name,
                        success=False,
                        error="Adapter disabled or not found",
                        execution_time=0,
                        cache_hit=False
                    )

                # Extract input
                input_value = request.input_data.get('smiles') or request.input_data.get('input') or request.input_data

                result: AdapterResult = await adapter(
                    input_value,
                    use_cache=request.use_cache,
                    **request.params
                )

                return BatchAdapterResult(
                    adapter=adapter_name,
                    success=result.success,
                    data=result.data,
                    error=result.error,
                    execution_time=time.time() - adapter_start,
                    cache_hit=result.cache_hit
                )

            except Exception as e:
                logger.error(f"Error in batch execution for {adapter_name}: {e}")
                return BatchAdapterResult(
                    adapter=adapter_name,
                    success=False,
                    error=str(e),
                    execution_time=time.time() - adapter_start,
                    cache_hit=False
                )

        # Run all adapters concurrently
        results = await asyncio.gather(*[
            execute_single(name) for name in request.adapter_names
        ])

        total_time = time.time() - start_time
        successful = sum(1 for r in results if r.success)
        failed = len(results) - successful

        logger.info(f"Batch execution completed in {total_time:.3f}s: {successful} succeeded, {failed} failed")

        return BatchExecuteResponse(
            results=list(results),
            total_time=total_time,
            successful_count=successful,
            failed_count=failed
        )

    except HTTPException:
        raise
    except Exception as e:
        total_time = time.time() - start_time
        logger.error(f"Error in batch execution: {e}", exc_info=True)
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Batch execution failed: {str(e)}"
        )


@router.get("/health")
async def adapters_health():
    """
    Check health of adapter system

    Returns:
        Health status including adapter count and registry status
    """
    try:
        adapter_count = len(registry.list_adapters())

        return {
            "status": "healthy",
            "adapter_count": adapter_count,
            "registry_initialized": adapter_count > 0,
            "timestamp": time.time()
        }
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return {
            "status": "unhealthy",
            "error": str(e),
            "timestamp": time.time()
        }
