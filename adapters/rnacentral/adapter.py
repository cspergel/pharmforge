"""
RNAcentral Adapter - Query RNA sequences and annotations from RNAcentral
Provides access to RNA data for RNA-targeting drug discovery
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class RNAcentralAdapter(AdapterProtocol):
    """
    Adapter for RNAcentral REST API
    Searches RNA sequences, annotations, and cross-references

    Useful for:
    - Finding target RNAs for therapeutic development
    - Getting RNA sequence and structure information
    - Retrieving RNA functional annotations
    - Accessing cross-references to other databases
    """

    BASE_URL = "https://rnacentral.org/api/v1"

    def __init__(self):
        super().__init__(
            name="rnacentral",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.3,  # Be respectful to RNAcentral API
                "timeout": 60,
                "max_results": 100
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid RNA ID, keyword, or query

        Args:
            input_data: RNA ID (URS format), keyword, or search query

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        return True

    async def _search_by_keyword_async(
        self,
        query: str,
        rna_type: Optional[str] = None,
        organism: Optional[str] = None
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Search RNAcentral by keyword

        Args:
            query: Search query (keyword, gene name, etc.)
            rna_type: Optional RNA type filter (e.g., "miRNA", "lncRNA")
            organism: Optional organism filter (e.g., "Homo sapiens")

        Returns:
            List of RNA entries or None on error
        """
        search_url = f"{self.BASE_URL}/rna"

        # Build search parameters
        params = {
            "search": query,
            "page_size": self.config.get("max_results", 100)
        }

        if rna_type:
            params["rna_type"] = rna_type

        if organism:
            params["organism"] = organism

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(search_url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("results", [])
                        logger.info(f"RNAcentral: Found {len(results)} results for query: {query}")
                        return results
                    elif response.status == 404:
                        logger.warning(f"RNAcentral: No results found for query: {query}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"RNAcentral API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"RNAcentral: Timeout searching for {query}")
            return None
        except Exception as e:
            logger.error(f"RNAcentral: Error searching for {query}: {e}")
            return None

    async def _get_rna_by_id_async(self, rna_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information for a specific RNA by RNAcentral ID

        Args:
            rna_id: RNAcentral URS identifier (e.g., "URS0000000001")

        Returns:
            RNA details or None on error
        """
        # Clean the RNA ID
        rna_id = rna_id.strip().upper()

        rna_url = f"{self.BASE_URL}/rna/{rna_id}"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(rna_url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        logger.info(f"RNAcentral: Retrieved RNA details for {rna_id}")
                        return data
                    elif response.status == 404:
                        logger.warning(f"RNAcentral: RNA ID not found: {rna_id}")
                        return None
                    else:
                        error_text = await response.text()
                        logger.error(f"RNAcentral API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"RNAcentral: Timeout fetching RNA {rna_id}")
            return None
        except Exception as e:
            logger.error(f"RNAcentral: Error fetching RNA {rna_id}: {e}")
            return None

    async def _get_xrefs_async(self, rna_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get cross-references for a specific RNA

        Args:
            rna_id: RNAcentral URS identifier

        Returns:
            List of cross-references or None on error
        """
        rna_id = rna_id.strip().upper()
        xref_url = f"{self.BASE_URL}/rna/{rna_id}/xrefs"

        params = {"page_size": self.config.get("max_results", 100)}
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(xref_url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        xrefs = data.get("results", [])
                        logger.info(f"RNAcentral: Found {len(xrefs)} cross-references for {rna_id}")
                        return xrefs
                    elif response.status == 404:
                        logger.warning(f"RNAcentral: No cross-references found for {rna_id}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"RNAcentral xref API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"RNAcentral: Timeout fetching cross-references for {rna_id}")
            return None
        except Exception as e:
            logger.error(f"RNAcentral: Error fetching cross-references for {rna_id}: {e}")
            return None

    def _is_rna_id(self, query: str) -> bool:
        """
        Check if query looks like an RNAcentral ID (URS format)

        Args:
            query: Input query string

        Returns:
            True if it looks like an RNA ID
        """
        query = query.strip().upper()
        # RNAcentral IDs start with "URS" followed by digits
        return query.startswith("URS") and len(query) >= 13

    def _summarize_rna_data(self, rna_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Summarize RNA data into standardized format

        Args:
            rna_data: Raw RNA data from RNAcentral

        Returns:
            Summarized RNA information
        """
        return {
            "rnacentral_id": rna_data.get("rnacentral_id"),
            "sequence": rna_data.get("sequence"),
            "length": rna_data.get("length"),
            "rna_type": rna_data.get("rna_type"),
            "description": rna_data.get("description"),
            "species": [
                {
                    "taxid": sp.get("taxid") if isinstance(sp, dict) else None,
                    "name": sp.get("name") if isinstance(sp, dict) else str(sp),
                    "common_name": sp.get("common_name") if isinstance(sp, dict) else None
                }
                for sp in (rna_data.get("xrefs", []) if isinstance(rna_data.get("xrefs"), list) else [])[:5]  # Top 5 species
            ] if "xrefs" in rna_data else [],
            "url": rna_data.get("url", f"https://rnacentral.org/rna/{rna_data.get('rnacentral_id')}")
        }

    def _summarize_search_results(self, results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Summarize search results into useful metrics

        Args:
            results: List of RNA search results

        Returns:
            Summarized search information
        """
        if not results:
            return {
                "num_results": 0,
                "rna_types": [],
                "organisms": [],
                "entries": []
            }

        # Extract RNA types and organisms
        rna_types = {}
        organisms = {}

        for result in results:
            # Count RNA types
            rna_type = result.get("rna_type", "unknown")
            rna_types[rna_type] = rna_types.get(rna_type, 0) + 1

            # Extract organisms from descriptions
            description = result.get("description", "")
            # Simple organism extraction (can be improved)
            if "Homo sapiens" in description:
                organisms["Homo sapiens"] = organisms.get("Homo sapiens", 0) + 1
            elif "Mus musculus" in description:
                organisms["Mus musculus"] = organisms.get("Mus musculus", 0) + 1

        # Format top results
        entries = [
            {
                "rnacentral_id": r.get("rnacentral_id"),
                "description": r.get("description"),
                "rna_type": r.get("rna_type"),
                "length": r.get("length"),
                "url": f"https://rnacentral.org/rna/{r.get('rnacentral_id')}"
            }
            for r in results[:10]  # Top 10 results
        ]

        return {
            "num_results": len(results),
            "rna_types": [{"type": k, "count": v} for k, v in sorted(rna_types.items(), key=lambda x: x[1], reverse=True)],
            "organisms": [{"name": k, "count": v} for k, v in sorted(organisms.items(), key=lambda x: x[1], reverse=True)],
            "entries": entries
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the RNAcentral query

        Args:
            input_data: RNA ID (URS format) or search query
            **kwargs: Additional parameters
                - rna_type: Filter by RNA type (e.g., "miRNA", "lncRNA")
                - organism: Filter by organism (e.g., "Homo sapiens")
                - get_xrefs: Whether to fetch cross-references (default: False)

        Returns:
            AdapterResult containing RNA data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be RNA ID or search query"
            )

        query = input_data.strip()

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.3)
        await asyncio.sleep(rate_delay)

        # Determine if this is an ID lookup or search
        is_id = self._is_rna_id(query)

        if is_id:
            # Direct ID lookup
            rna_data = await self._get_rna_by_id_async(query)

            if rna_data is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Failed to retrieve RNA with ID: {query}",
                    metadata={"source": "rnacentral", "query": query, "query_type": "id"}
                )

            # Optionally get cross-references
            xrefs = None
            if kwargs.get("get_xrefs", False):
                xrefs = await self._get_xrefs_async(query)

            # Summarize the data
            summary = self._summarize_rna_data(rna_data)

            if xrefs:
                summary["cross_references"] = xrefs[:20]  # Top 20 xrefs

            result_data = {
                **summary,
                "raw_data": rna_data  # Include full data for advanced use
            }

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "rnacentral",
                    "query": query,
                    "query_type": "id",
                    "adapter_version": self.version
                }
            )

        else:
            # Keyword search
            rna_type = kwargs.get("rna_type")
            organism = kwargs.get("organism")

            results = await self._search_by_keyword_async(query, rna_type, organism)

            if results is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Failed to search RNAcentral for: {query}",
                    metadata={
                        "source": "rnacentral",
                        "query": query,
                        "query_type": "search"
                    }
                )

            # Summarize search results
            summary = self._summarize_search_results(results)

            result_data = {
                **summary,
                "raw_results": results  # Include full results
            }

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "rnacentral",
                    "query": query,
                    "query_type": "search",
                    "rna_type_filter": rna_type,
                    "organism_filter": organism,
                    "adapter_version": self.version
                }
            )
