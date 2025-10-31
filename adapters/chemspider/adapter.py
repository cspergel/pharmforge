"""
ChemSpider Adapter - Search and retrieve chemical data from Royal Society of Chemistry database
Aggregates data from 100+ sources including PubChem, ChEMBL, DrugBank, etc.
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
import os

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ChemSpiderAdapter(AdapterProtocol):
    """
    Adapter for ChemSpider API (Royal Society of Chemistry)
    Provides chemical structure search, property lookup, and aggregated data from 100+ sources

    Features:
    - Search by SMILES, InChI, InChIKey, name, or formula
    - Retrieve compound properties and identifiers
    - Get synonyms and trade names
    - Access safety and regulatory information
    - Cross-database validation

    API Documentation: https://developer.rsc.org/compounds-v1/apis
    """

    BASE_URL = "https://api.rsc.org/compounds/v1"

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize ChemSpider adapter

        Args:
            api_key: ChemSpider API key (optional, will try env var CHEMSPIDER_API_KEY)
        """
        # Try to get API key from parameter or environment
        self.api_key = api_key or os.getenv("CHEMSPIDER_API_KEY")

        super().__init__(
            name="chemspider",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,  # 2 requests/second to be respectful
                "timeout": 60,
                "max_results": 10
            }
        )
        self.version = "1.0.0"

        if not self.api_key:
            logger.warning("ChemSpider API key not provided - some features may be limited")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data

        Args:
            input_data: Can be a SMILES string or dict with query parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            # String input (SMILES, InChI, etc.)
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            # Dict input with query_type and query
            return "query" in input_data and len(input_data["query"]) > 0
        return False

    def _get_headers(self) -> Dict[str, str]:
        """
        Get API request headers

        Returns:
            Dictionary of HTTP headers
        """
        headers = {
            "Content-Type": "application/json",
            "User-Agent": "PharmForge/1.0 (Drug Discovery Platform)"
        }
        if self.api_key:
            headers["apikey"] = self.api_key
        return headers

    async def _search_by_identifier_async(
        self,
        query: str,
        query_type: str = "smiles"
    ) -> Optional[List[int]]:
        """
        Search ChemSpider by various identifier types

        Args:
            query: Search term (SMILES, InChI, name, formula, etc.)
            query_type: Type of query (smiles, inchi, inchikey, name, formula)

        Returns:
            List of ChemSpider IDs or None on error
        """
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))
        headers = self._get_headers()

        # Map query types to API endpoints
        endpoint_map = {
            "smiles": f"{self.BASE_URL}/filter/smiles",
            "inchi": f"{self.BASE_URL}/filter/inchi",
            "inchikey": f"{self.BASE_URL}/filter/inchikey",
            "formula": f"{self.BASE_URL}/filter/formula",
            "name": f"{self.BASE_URL}/filter/name"
        }

        endpoint = endpoint_map.get(query_type.lower())
        if not endpoint:
            logger.error(f"Unknown query type: {query_type}")
            return None

        try:
            async with aiohttp.ClientSession() as session:
                # POST search query
                payload = {query_type.lower(): query}

                async with session.post(
                    endpoint,
                    json=payload,
                    headers=headers,
                    timeout=timeout
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        # ChemSpider returns a queryId for async results
                        query_id = data.get("queryId")

                        if not query_id:
                            logger.warning(f"ChemSpider: No queryId returned for {query}")
                            return []

                        # Poll for results (ChemSpider uses async pattern)
                        return await self._get_search_results_async(query_id)

                    elif response.status == 404:
                        logger.info(f"ChemSpider: No results for {query_type}: {query}")
                        return []
                    elif response.status == 401:
                        logger.error("ChemSpider: Authentication failed - check API key")
                        return None
                    else:
                        error_text = await response.text()
                        logger.error(f"ChemSpider search error {response.status}: {error_text}")
                        return None

        except asyncio.TimeoutError:
            logger.error(f"ChemSpider: Timeout searching for {query}")
            return None
        except Exception as e:
            logger.error(f"ChemSpider: Error searching for {query}: {e}")
            return None

    async def _get_search_results_async(self, query_id: str) -> Optional[List[int]]:
        """
        Poll for search results using queryId

        Args:
            query_id: ChemSpider query ID

        Returns:
            List of ChemSpider IDs
        """
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))
        headers = self._get_headers()
        results_url = f"{self.BASE_URL}/filter/{query_id}/results"

        max_attempts = 10
        attempt = 0

        try:
            async with aiohttp.ClientSession() as session:
                while attempt < max_attempts:
                    await asyncio.sleep(1)  # Wait between polls

                    async with session.get(results_url, headers=headers, timeout=timeout) as response:
                        if response.status == 200:
                            data = await response.json()
                            results = data.get("results", [])
                            logger.info(f"ChemSpider: Retrieved {len(results)} IDs")
                            return results[:self.config.get("max_results", 10)]
                        elif response.status == 404:
                            # Results not ready yet
                            attempt += 1
                            continue
                        else:
                            logger.error(f"ChemSpider results error: {response.status}")
                            return None

                logger.warning("ChemSpider: Results polling timed out")
                return None

        except Exception as e:
            logger.error(f"ChemSpider: Error polling results: {e}")
            return None

    async def _get_compound_details_async(self, chemspider_id: int) -> Optional[Dict[str, Any]]:
        """
        Get detailed compound information by ChemSpider ID

        Args:
            chemspider_id: ChemSpider compound ID

        Returns:
            Dictionary of compound details
        """
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))
        headers = self._get_headers()

        details_url = f"{self.BASE_URL}/records/{chemspider_id}/details"

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(details_url, headers=headers, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    else:
                        logger.warning(f"ChemSpider: Could not fetch details for ID {chemspider_id}")
                        return None

        except Exception as e:
            logger.error(f"ChemSpider: Error fetching details for ID {chemspider_id}: {e}")
            return None

    def _parse_compound_data(self, details: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse and standardize compound data from ChemSpider API response

        Args:
            details: Raw API response data

        Returns:
            Standardized compound data dictionary
        """
        # Extract key fields from ChemSpider response
        result = {
            "chemspider_id": str(details.get("id", "")),
            "common_name": details.get("commonName", ""),
            "formula": details.get("formula", ""),
            "molecular_weight": details.get("molecularWeight"),
            "smiles": details.get("smiles", ""),
            "inchi": details.get("inchi", ""),
            "inchikey": details.get("inchiKey", ""),
            "synonyms": [],
            "properties": {},
            "data_sources": []
        }

        # Extract synonyms
        if "names" in details:
            result["synonyms"] = details["names"][:20]  # Limit to first 20

        # Extract properties if available
        if "properties" in details:
            props = details["properties"]
            result["properties"] = {
                "alogp": props.get("aLogP"),
                "xlogp": props.get("xLogP"),
                "molecular_formula": props.get("molecularFormula"),
                "nominal_mass": props.get("nominalMass"),
                "average_mass": props.get("averageMass"),
                "monoisotopic_mass": props.get("monoisotopicMass")
            }

        # Extract data sources
        if "dataSources" in details:
            result["data_sources"] = details["dataSources"][:10]  # Top 10 sources

        return result

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute ChemSpider search and retrieval

        Args:
            input_data: Either a SMILES string or dict with:
                - query: Search term
                - query_type: Type of search (smiles, inchi, name, formula, inchikey)
                - include_properties: Whether to fetch detailed properties (default: True)
                - include_synonyms: Whether to fetch synonyms (default: True)
                - max_results: Maximum number of results (default: 10)
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing search results and compound details
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input - must be a string or dict with 'query' field"
            )

        # Parse input
        if isinstance(input_data, str):
            # Assume SMILES by default for string input
            query = input_data
            query_type = "smiles"
        else:
            query = input_data["query"]
            query_type = input_data.get("query_type", "smiles")

        include_properties = kwargs.get("include_properties", True)
        include_synonyms = kwargs.get("include_synonyms", True)
        max_results = kwargs.get("max_results", self.config.get("max_results", 10))

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Search for compounds
        chemspider_ids = await self._search_by_identifier_async(query, query_type)

        if chemspider_ids is None:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to search ChemSpider for {query_type}: {query}",
                metadata={"source": "chemspider", "query": query, "query_type": query_type}
            )

        if not chemspider_ids:
            # No results found - this is not an error, just empty results
            return AdapterResult(
                success=True,
                data={
                    "results": [],
                    "num_results": 0,
                    "warnings": [f"No ChemSpider entries found for {query_type}: {query}"]
                },
                metadata={
                    "source": "chemspider",
                    "query": query,
                    "query_type": query_type
                }
            )

        # Fetch details for each compound (up to max_results)
        results = []
        warnings = []

        for chemspider_id in chemspider_ids[:max_results]:
            await asyncio.sleep(0.2)  # Additional rate limiting for detail fetches

            details = await self._get_compound_details_async(chemspider_id)

            if details:
                parsed = self._parse_compound_data(details)
                results.append(parsed)
            else:
                warnings.append(f"Could not fetch details for ChemSpider ID: {chemspider_id}")

        logger.info(f"ChemSpider: Retrieved {len(results)} compounds for {query}")

        return AdapterResult(
            success=True,
            data={
                "results": results,
                "num_results": len(results),
                "total_found": len(chemspider_ids),
                "warnings": warnings
            },
            cache_hit=False,
            metadata={
                "source": "chemspider",
                "query": query,
                "query_type": query_type,
                "adapter_version": self.version,
                "chemspider_ids": chemspider_ids[:max_results]
            }
        )
