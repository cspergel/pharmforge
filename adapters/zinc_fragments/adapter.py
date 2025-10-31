"""
ZINC Fragments Adapter - Migrated to ChEMBL API for fragment library search
Provides fragment search, similarity search, substructure search, and drug properties
Updated: 2025-10 - Migrated from ZINC15/ZINC20 (500 errors) to ChEMBL REST API
Migration reason: ZINC15 API returning 500 errors, SmallWorld API (sw.docking.org) has 502 proxy errors
ChEMBL provides: similarity search, substructure search, full molecular properties, drug data
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ZINCFragmentsAdapter(AdapterProtocol):
    """
    Adapter for fragment library searches using ChEMBL API

    NOTE: This adapter was originally designed for ZINC15/ZINC20 but has been migrated to ChEMBL
    due to persistent 500 errors from ZINC15 API. The adapter maintains backward compatibility
    with the "zinc_fragments" name and interface.

    Capabilities:
    - Fragment search by SMILES (similarity and substructure)
    - Similarity search with configurable Tanimoto threshold
    - Substructure search
    - Comprehensive molecular properties (MW, LogP, HBA, HBD, etc.)
    - Drug/clinical data from ChEMBL (max_phase, indications, etc.)
    - Property filtering (MW, LogP, etc.)

    ChEMBL Advantages over ZINC:
    - Stable, well-maintained REST API
    - Rich drug and bioactivity data
    - No bot-detection issues
    - Comprehensive molecular properties
    - Active development and support
    """

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self):
        super().__init__(
            name="zinc_fragments",  # Maintain name for backward compatibility
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,  # ChEMBL is more permissive
                "timeout": 60,
                "default_similarity": 70,  # ChEMBL uses percentage (0-100)
                "max_results": 100,
                "fragment_mw_max": 300,  # Maximum MW for fragments
                "max_retries": 3,
                "backoff_base": 2
            }
        )
        self.version = "3.0.0"  # Major version bump - API migration

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string or search params

        Args:
            input_data: SMILES string or dict with search parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return "smiles" in input_data or "query" in input_data
        return False

    def _get_headers(self) -> Dict[str, str]:
        """
        Get HTTP headers for ChEMBL API requests

        Returns:
            Dictionary of HTTP headers
        """
        return {
            'User-Agent': 'PharmForge/3.0 (https://pharmforge.org)',
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        }

    async def _retry_request(self, func, *args, **kwargs):
        """
        Retry a request with exponential backoff

        Args:
            func: Async function to retry
            *args, **kwargs: Arguments to pass to func

        Returns:
            Result of func or None on failure
        """
        max_retries = self.config.get("max_retries", 3)
        backoff_base = self.config.get("backoff_base", 2)

        for attempt in range(max_retries):
            try:
                result = await func(*args, **kwargs)
                return result
            except (asyncio.TimeoutError, aiohttp.ClientError) as e:
                if attempt == max_retries - 1:
                    logger.error(f"ZINC: All {max_retries} retry attempts failed: {e}")
                    return None

                wait_time = backoff_base ** attempt  # Exponential backoff: 2s, 4s, 8s
                logger.warning(f"ZINC: Retry {attempt + 1}/{max_retries} after {wait_time}s delay (error: {type(e).__name__})")
                await asyncio.sleep(wait_time)

        return None

    async def _search_fragments_async(
        self,
        smiles: str,
        search_type: str = "similarity",
        similarity_threshold: int = 70,
        max_results: int = 100
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Search ChEMBL for fragments with retry logic

        Args:
            smiles: SMILES string to search
            search_type: "similarity" or "substructure"
            similarity_threshold: Tanimoto threshold for similarity (0-100 percentage)
            max_results: Maximum number of results

        Returns:
            List of fragment records or None on error
        """
        encoded_smiles = quote(smiles, safe='')
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))
        headers = self._get_headers()

        # Construct search URL based on type
        if search_type == "similarity":
            # ChEMBL similarity endpoint: /similarity/{smiles}/{similarity}.json
            url = f"{self.BASE_URL}/similarity/{encoded_smiles}/{similarity_threshold}.json"
        elif search_type == "substructure":
            # ChEMBL substructure endpoint: /substructure/{smiles}.json
            url = f"{self.BASE_URL}/substructure/{encoded_smiles}.json"
        else:
            logger.error(f"Invalid search type: {search_type}. ChEMBL supports 'similarity' or 'substructure'")
            return None

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()

                        # ChEMBL returns {"molecules": [...], "page_meta": {...}}
                        molecules = data.get("molecules", [])

                        # Limit results
                        molecules = molecules[:max_results]

                        logger.info(f"ChEMBL: Found {len(molecules)} fragments for {smiles} (search_type={search_type})")
                        return molecules
                    elif response.status == 404:
                        logger.info(f"ChEMBL: No fragments found for {smiles}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"ChEMBL API error {response.status}: {error_text[:200]}")
                        return None

        return await self._retry_request(make_request)

    async def _get_molecule_details_async(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information for a ChEMBL molecule with retry logic

        Args:
            chembl_id: ChEMBL ID (e.g., CHEMBL545)

        Returns:
            Dictionary with detailed molecule info or None
        """
        url = f"{self.BASE_URL}/molecule/{chembl_id}.json"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))
        headers = self._get_headers()

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()

                        # Extract additional details
                        details = {
                            "chembl_id": chembl_id,
                            "max_phase": data.get("max_phase"),
                            "first_approval": data.get("first_approval"),
                            "oral": data.get("oral"),
                            "parenteral": data.get("parenteral"),
                            "topical": data.get("topical"),
                            "black_box_warning": data.get("black_box_warning"),
                            "availability_type": data.get("availability_type"),
                            "therapeutic_flag": data.get("therapeutic_flag"),
                            "withdrawn_flag": data.get("withdrawn_flag"),
                            "atc_classifications": data.get("atc_classifications", [])
                        }

                        return details
                    else:
                        logger.warning(f"ChEMBL: Could not get details for {chembl_id} (status: {response.status})")
                        return None

        return await self._retry_request(make_request)

    def _filter_fragments(
        self,
        results: List[Dict[str, Any]],
        mw_max: Optional[float] = None,
        logp_max: Optional[float] = None,
        hbd_max: Optional[int] = None,
        hba_max: Optional[int] = None
    ) -> List[Dict[str, Any]]:
        """
        Filter fragments by property constraints

        Args:
            results: List of ChEMBL molecule results
            mw_max: Maximum molecular weight
            logp_max: Maximum LogP
            hbd_max: Maximum H-bond donors
            hba_max: Maximum H-bond acceptors

        Returns:
            Filtered list of fragments
        """
        filtered = []

        for result in results:
            # ChEMBL stores properties in nested 'molecule_properties' dict
            props = result.get("molecule_properties", {})

            # Extract properties
            mw = props.get("full_mwt") or props.get("mw_freebase")
            logp = props.get("alogp")
            hbd = props.get("hbd")
            hba = props.get("hba")

            # Apply filters
            if mw_max is not None and mw is not None and float(mw) > mw_max:
                continue
            if logp_max is not None and logp is not None and float(logp) > logp_max:
                continue
            if hbd_max is not None and hbd is not None and int(hbd) > hbd_max:
                continue
            if hba_max is not None and hba is not None and int(hba) > hba_max:
                continue

            filtered.append(result)

        return filtered

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute fragment search using ChEMBL API

        Args:
            input_data: SMILES string or dict with search parameters
            **kwargs: Additional parameters:
                - search_type: "similarity" or "substructure"
                - similarity_threshold: float (0-1 range, converted to 0-100 for ChEMBL)
                - max_results: int
                - get_details: bool (fetch additional molecule details)
                - mw_max: float
                - logp_max: float
                - hbd_max: int
                - hba_max: int

        Returns:
            AdapterResult containing fragment search results from ChEMBL
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: expected SMILES string or search dict"
            )

        # Parse input
        if isinstance(input_data, str):
            smiles = input_data
        else:
            smiles = input_data.get("smiles") or input_data.get("query")

        # Get search parameters
        search_type = kwargs.get("search_type", "similarity")

        # Handle similarity threshold - convert from 0-1 to 0-100 if needed
        similarity_threshold = kwargs.get(
            "similarity_threshold",
            self.config.get("default_similarity", 70)
        )

        # Convert 0-1 range to 0-100 percentage for ChEMBL
        if isinstance(similarity_threshold, float) and similarity_threshold <= 1.0:
            similarity_threshold = int(similarity_threshold * 100)
        else:
            similarity_threshold = int(similarity_threshold)

        max_results = kwargs.get("max_results", self.config.get("max_results", 100))
        get_details = kwargs.get("get_details", False)

        # Property filters
        mw_max = kwargs.get("mw_max", self.config.get("fragment_mw_max", 300))
        logp_max = kwargs.get("logp_max")
        hbd_max = kwargs.get("hbd_max")
        hba_max = kwargs.get("hba_max")

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        logger.debug(f"ChEMBL: Rate limiting - waiting {rate_delay}s before request")
        await asyncio.sleep(rate_delay)

        # Search ChEMBL
        results = await self._search_fragments_async(
            smiles,
            search_type=search_type,
            similarity_threshold=similarity_threshold,
            max_results=max_results
        )

        if results is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search ChEMBL fragment library",
                metadata={"source": "chembl", "smiles": smiles}
            )

        # Filter by properties
        filtered_results = self._filter_fragments(
            results,
            mw_max=mw_max,
            logp_max=logp_max,
            hbd_max=hbd_max,
            hba_max=hba_max
        )

        # Optionally get detailed info (with rate limiting between requests)
        if get_details and filtered_results:
            logger.info(f"ChEMBL: Fetching detailed info for {min(len(filtered_results), 20)} compounds")

            for i, result in enumerate(filtered_results[:20]):  # Limit to top 20
                chembl_id = result.get("molecule_chembl_id")
                if chembl_id:
                    # Add delay between individual detail requests
                    if i > 0:  # Don't delay before first request
                        await asyncio.sleep(self.config.get("rate_limit_delay", 0.5))

                    details = await self._get_molecule_details_async(chembl_id)
                    if details:
                        result["drug_details"] = details

        result_data = {
            "query_smiles": smiles,
            "search_type": search_type,
            "num_results": len(filtered_results),
            "num_total": len(results),
            "fragments": filtered_results,
            "filters_applied": {
                "mw_max": mw_max,
                "logp_max": logp_max,
                "hbd_max": hbd_max,
                "hba_max": hba_max
            }
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "chembl",
                "adapter_version": self.version,
                "search_type": search_type,
                "similarity_threshold": similarity_threshold if search_type == "similarity" else None,
                "migration_note": "Migrated from ZINC15 to ChEMBL due to API availability"
            }
        )
