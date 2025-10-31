"""
SureChEMBL Adapter - Patent chemistry search
Searches patent databases for compound mentions, similar structures, and publication data
Updated: 2025-10 - Updated to SureChEMBL 2.0 API with ChEMBL fallback
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class SureChEMBLAdapter(AdapterProtocol):
    """
    Adapter for SureChEMBL 2.0 patent chemistry database

    Capabilities:
    - Patent compound search by structure
    - Similarity search in patent literature
    - Patent ID lookup
    - Extract compounds from patents
    - Publication dates and patent families
    - ChEMBL API fallback for additional drug data

    Note: SureChEMBL 2.0 was released in May 2025 with improved pipeline and RDKit migration
    """

    # SureChEMBL 2.0 API endpoints
    BASE_URL = "https://www.surechembl.org/api/v1"
    # ChEMBL API as fallback (includes SureChEMBL data)
    CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self):
        super().__init__(
            name="surechembl",
            adapter_type="api",
            config={
                "rate_limit_delay": 1.0,  # Be respectful to SureChEMBL
                "timeout": 90,
                "default_similarity": 0.8,
                "max_results": 50,
                "max_retries": 3,
                "use_chembl_fallback": True  # Use ChEMBL API if SureChEMBL fails
            }
        )
        self.version = "2.0.0"  # Updated for SureChEMBL 2.0

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input (SMILES, patent ID, or search dict)

        Args:
            input_data: SMILES string, patent ID, or dict

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return "smiles" in input_data or "patent_id" in input_data
        return False

    def _get_headers(self) -> Dict[str, str]:
        """
        Get HTTP headers for API requests

        Returns:
            Dictionary of HTTP headers
        """
        return {
            'User-Agent': 'PharmForge/2.0 (Drug Discovery Platform; +https://pharmforge.org)',
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

        for attempt in range(max_retries):
            try:
                result = await func(*args, **kwargs)
                return result
            except (asyncio.TimeoutError, aiohttp.ClientError) as e:
                if attempt == max_retries - 1:
                    logger.error(f"SureChEMBL: All {max_retries} retry attempts failed: {e}")
                    return None

                wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
                logger.warning(f"SureChEMBL: Retry {attempt + 1}/{max_retries} after {wait_time}s delay")
                await asyncio.sleep(wait_time)

        return None

    async def _search_chembl_fallback(
        self,
        smiles: str,
        similarity_threshold: float = 0.8,
        max_results: int = 50
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Fallback to ChEMBL API for structure search (includes SureChEMBL data)

        Args:
            smiles: SMILES string
            similarity_threshold: Similarity threshold (70-100)
            max_results: Maximum results

        Returns:
            List of results from ChEMBL or None
        """
        url = f"{self.CHEMBL_BASE_URL}/similarity/{smiles}/{int(similarity_threshold * 100)}.json"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))
        headers = self._get_headers()

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        molecules = data.get("molecules", [])
                        logger.info(f"SureChEMBL/ChEMBL fallback: Found {len(molecules)} compounds")
                        return molecules[:max_results]
                    elif response.status == 404:
                        logger.info(f"SureChEMBL/ChEMBL fallback: No results found")
                        return []
                    else:
                        logger.warning(f"ChEMBL API error {response.status}")
                        return None

        return await self._retry_request(make_request)

    async def _search_by_structure_async(
        self,
        smiles: str,
        search_type: str = "similarity",
        similarity_threshold: float = 0.8,
        max_results: int = 50
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Search SureChEMBL by chemical structure with ChEMBL fallback

        Args:
            smiles: SMILES string
            search_type: "similarity", "substructure", or "exact"
            similarity_threshold: Tanimoto threshold (0-1)
            max_results: Maximum results

        Returns:
            List of patent compound records or None
        """
        encoded_smiles = quote(smiles)
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))
        headers = self._get_headers()

        # Construct search URL for SureChEMBL 2.0
        if search_type == "similarity":
            url = f"{self.BASE_URL}/search"
            params = {
                "q": smiles,
                "type": "similarity",
                "threshold": int(similarity_threshold * 100),
                "limit": max_results
            }
        elif search_type == "substructure":
            url = f"{self.BASE_URL}/search"
            params = {
                "q": smiles,
                "type": "substructure",
                "limit": max_results
            }
        elif search_type == "exact":
            url = f"{self.BASE_URL}/search"
            params = {
                "q": smiles,
                "type": "exact",
                "limit": max_results
            }
        else:
            logger.error(f"Invalid search type: {search_type}")
            return None

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()

                        # Parse results
                        if isinstance(data, dict):
                            results = data.get("results", data.get("compounds", []))
                        elif isinstance(data, list):
                            results = data
                        else:
                            results = []

                        logger.info(f"SureChEMBL: Found {len(results)} patent compounds for {smiles}")
                        return results

                    elif response.status == 404:
                        logger.info(f"SureChEMBL: No compounds found for {smiles}")
                        return []
                    elif response.status == 500:
                        logger.warning(f"SureChEMBL: Server error 500 - will try ChEMBL fallback")
                        raise aiohttp.ClientError("SureChEMBL server error")
                    else:
                        error_text = await response.text()
                        logger.error(f"SureChEMBL API error {response.status}: {error_text[:200]}")
                        return None

        # Try SureChEMBL first
        results = await self._retry_request(make_request)

        # If SureChEMBL fails and fallback is enabled, try ChEMBL
        if results is None and self.config.get("use_chembl_fallback", True):
            logger.info("SureChEMBL failed, trying ChEMBL API fallback...")
            results = await self._search_chembl_fallback(smiles, similarity_threshold, max_results)

            if results:
                # Add metadata to indicate this came from ChEMBL fallback
                for result in results:
                    result["source"] = "chembl_fallback"

        return results

    async def _get_patent_details_async(self, patent_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a patent

        Args:
            patent_id: Patent ID (e.g., US20190123456)

        Returns:
            Dictionary with patent details or None
        """
        url = f"{self.BASE_URL}/document/{patent_id}"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))
        headers = self._get_headers()

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()

                        patent_info = {
                            "patent_id": patent_id,
                            "title": data.get("title"),
                            "abstract": data.get("abstract"),
                            "publication_date": data.get("published") or data.get("publication_date"),
                            "applicants": data.get("applicants", []),
                            "inventors": data.get("inventors", []),
                            "patent_family": data.get("family") or data.get("family_id"),
                            "num_compounds": len(data.get("compounds", []))
                        }

                        return patent_info
                    elif response.status == 500:
                        logger.warning(f"SureChEMBL: Server error 500 for patent {patent_id}")
                        raise aiohttp.ClientError("SureChEMBL server error")
                    else:
                        logger.warning(f"SureChEMBL: Could not get patent details for {patent_id} (status: {response.status})")
                        return None

        return await self._retry_request(make_request)

    async def _extract_compounds_from_patent_async(
        self,
        patent_id: str,
        max_compounds: int = 100
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Extract all chemical compounds mentioned in a patent

        Args:
            patent_id: Patent ID
            max_compounds: Maximum compounds to extract

        Returns:
            List of compound structures or None
        """
        url = f"{self.BASE_URL}/document/{patent_id}/compounds"
        params = {"limit": max_compounds}
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))
        headers = self._get_headers()

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        compounds = data.get("compounds", []) if isinstance(data, dict) else data

                        logger.info(f"SureChEMBL: Extracted {len(compounds)} compounds from {patent_id}")
                        return compounds
                    elif response.status == 500:
                        logger.warning(f"SureChEMBL: Server error 500 extracting from {patent_id}")
                        raise aiohttp.ClientError("SureChEMBL server error")
                    else:
                        logger.warning(f"SureChEMBL: Could not extract compounds from {patent_id} (status: {response.status})")
                        return []

        return await self._retry_request(make_request)

    def _summarize_patent_results(
        self,
        results: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Summarize patent search results

        Args:
            results: List of patent compound records

        Returns:
            Summary statistics
        """
        if not results:
            return {
                "num_patents": 0,
                "num_compounds": 0,
                "patent_families": [],
                "date_range": None,
                "top_applicants": []
            }

        # Extract unique patents and patent families
        patents = set()
        families = set()
        applicants = {}
        dates = []

        for result in results:
            patent_id = result.get("patent_id")
            if patent_id:
                patents.add(patent_id)

            family_id = result.get("family_id")
            if family_id:
                families.add(family_id)

            applicant = result.get("applicant")
            if applicant:
                applicants[applicant] = applicants.get(applicant, 0) + 1

            pub_date = result.get("publication_date")
            if pub_date:
                dates.append(pub_date)

        # Get date range
        date_range = None
        if dates:
            dates.sort()
            date_range = {"earliest": dates[0], "latest": dates[-1]}

        # Top applicants
        top_applicants = sorted(applicants.items(), key=lambda x: x[1], reverse=True)[:10]

        return {
            "num_patents": len(patents),
            "num_compounds": len(results),
            "num_patent_families": len(families),
            "date_range": date_range,
            "top_applicants": [{"name": name, "count": count} for name, count in top_applicants]
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute SureChEMBL patent search

        Args:
            input_data: SMILES string, patent ID, or search dict
            **kwargs: Additional parameters:
                - search_type: "similarity", "substructure", or "exact"
                - similarity_threshold: float (0-1)
                - max_results: int
                - mode: "structure_search" or "patent_lookup" or "extract_compounds"
                - patent_id: str (for patent lookup/extraction)

        Returns:
            AdapterResult containing patent search results
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: expected SMILES, patent ID, or search dict"
            )

        # Parse input and mode
        mode = kwargs.get("mode", "structure_search")

        if isinstance(input_data, str):
            if mode == "patent_lookup" or mode == "extract_compounds":
                patent_id = input_data
                smiles = None
            else:
                smiles = input_data
                patent_id = None
        else:
            smiles = input_data.get("smiles")
            patent_id = input_data.get("patent_id")

        # Get search parameters
        search_type = kwargs.get("search_type", "similarity")
        similarity_threshold = kwargs.get(
            "similarity_threshold",
            self.config.get("default_similarity", 0.8)
        )
        max_results = kwargs.get("max_results", self.config.get("max_results", 50))

        # Rate limiting
        await asyncio.sleep(self.config.get("rate_limit_delay", 1.0))

        # Execute based on mode
        if mode == "structure_search" and smiles:
            # Search by chemical structure
            results = await self._search_by_structure_async(
                smiles,
                search_type=search_type,
                similarity_threshold=similarity_threshold,
                max_results=max_results
            )

            if results is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to search SureChEMBL",
                    metadata={"source": "surechembl", "smiles": smiles}
                )

            summary = self._summarize_patent_results(results)

            result_data = {
                "query_smiles": smiles,
                "search_type": search_type,
                "mode": mode,
                **summary,
                "patents": results
            }

        elif mode == "patent_lookup" and patent_id:
            # Lookup patent details
            patent_info = await self._get_patent_details_async(patent_id)

            if patent_info is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Failed to get patent details for {patent_id}",
                    metadata={"source": "surechembl", "patent_id": patent_id}
                )

            result_data = {
                "mode": mode,
                "patent_id": patent_id,
                **patent_info
            }

        elif mode == "extract_compounds" and patent_id:
            # Extract compounds from patent
            compounds = await self._extract_compounds_from_patent_async(
                patent_id,
                max_compounds=max_results
            )

            if compounds is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Failed to extract compounds from {patent_id}",
                    metadata={"source": "surechembl", "patent_id": patent_id}
                )

            result_data = {
                "mode": mode,
                "patent_id": patent_id,
                "num_compounds": len(compounds),
                "compounds": compounds
            }

        else:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Invalid mode or missing required input: mode={mode}"
            )

        # Check if results came from ChEMBL fallback
        used_fallback = False
        if isinstance(result_data.get("patents"), list):
            used_fallback = any(r.get("source") == "chembl_fallback" for r in result_data.get("patents", []))

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "surechembl_v2" if not used_fallback else "chembl_fallback",
                "adapter_version": self.version,
                "mode": mode,
                "search_type": search_type if mode == "structure_search" else None,
                "used_chembl_fallback": used_fallback
            }
        )
