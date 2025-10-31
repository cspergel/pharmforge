"""
COCONUT Adapter - Natural Products Database
Searches and retrieves natural product structures and metadata from COCONUT
(COlleCtion of Open Natural prodUcTs)

Database: 400,000+ natural products
API: https://coconut.naturalproducts.net/api/
License: CC0 (public domain)
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class COCONUTAdapter(AdapterProtocol):
    """
    Adapter for COCONUT Natural Products Database

    Features:
    - Search by structure (SMILES, InChI)
    - Search by molecular properties (MW, LogP range)
    - Search by organism/taxonomy
    - Retrieve biological source information
    - Get molecular descriptors
    - Natural product annotations

    Use Cases:
    - Natural product library screening
    - Scaffold mining from nature
    - Biosource tracking
    - Taxonomic analysis
    - Natural product-inspired drug design
    """

    BASE_URL = "https://coconut.naturalproducts.net/api"

    # COCONUT API endpoints
    SEARCH_ENDPOINT = "/search"
    COMPOUND_ENDPOINT = "/compound"
    TAXONOMY_ENDPOINT = "/taxonomy"

    def __init__(self):
        super().__init__(
            name="coconut",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,  # Be respectful to public API
                "timeout": 60,
                "max_results": 50,
                "default_similarity_threshold": 0.8
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data structure

        Args:
            input_data: Dictionary containing query information or simple SMILES string

        Returns:
            True if valid, False otherwise
        """
        # Allow simple SMILES string
        if isinstance(input_data, str):
            return len(input_data) > 0

        # Or structured query dictionary
        if isinstance(input_data, dict):
            query_type = input_data.get("query_type")
            query = input_data.get("query")

            if not query_type or not query:
                return False

            valid_types = ["smiles", "name", "inchi", "inchikey", "properties", "organism", "id"]
            if query_type not in valid_types:
                return False

            return True

        return False

    def _parse_input(self, input_data: Any) -> Dict[str, Any]:
        """
        Parse and normalize input data

        Args:
            input_data: SMILES string or query dictionary

        Returns:
            Normalized query dictionary
        """
        # Simple SMILES string
        if isinstance(input_data, str):
            return {
                "query_type": "smiles",
                "query": input_data,
                "filters": {},
                "include_taxonomy": True,
                "include_activities": True,
                "max_results": self.config.get("max_results", 50)
            }

        # Already a dictionary - fill in defaults
        query_dict = input_data.copy()
        query_dict.setdefault("filters", {})
        query_dict.setdefault("include_taxonomy", True)
        query_dict.setdefault("include_activities", True)
        query_dict.setdefault("max_results", self.config.get("max_results", 50))

        return query_dict

    async def _search_by_smiles(self, smiles: str, max_results: int = 50) -> Optional[List[Dict[str, Any]]]:
        """
        Search COCONUT by SMILES similarity

        Args:
            smiles: SMILES string to search
            max_results: Maximum number of results to return

        Returns:
            List of matching compounds or None on error
        """
        # Note: COCONUT API structure search endpoint
        # This is a simplified implementation - actual API may vary
        url = f"{self.BASE_URL}/search/structure"

        params = {
            "smiles": smiles,
            "limit": max_results,
            "threshold": self.config.get("default_similarity_threshold", 0.8)
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("results", [])
                        logger.info(f"COCONUT: Found {len(results)} compounds for SMILES: {smiles}")
                        return results
                    elif response.status == 404:
                        logger.info(f"COCONUT: No compounds found for SMILES: {smiles}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"COCONUT API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"COCONUT: Timeout searching for {smiles}")
            return None
        except Exception as e:
            logger.error(f"COCONUT: Error searching for {smiles}: {e}")
            return None

    async def _search_by_name(self, name: str, max_results: int = 50) -> Optional[List[Dict[str, Any]]]:
        """
        Search COCONUT by compound name

        Args:
            name: Compound name to search
            max_results: Maximum number of results

        Returns:
            List of matching compounds or None on error
        """
        url = f"{self.BASE_URL}/search/name"

        params = {
            "query": name,
            "limit": max_results
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("results", [])
                        logger.info(f"COCONUT: Found {len(results)} compounds for name: {name}")
                        return results
                    else:
                        logger.warning(f"COCONUT: Name search failed with status {response.status}")
                        return []
        except Exception as e:
            logger.error(f"COCONUT: Error searching by name {name}: {e}")
            return None

    async def _search_by_organism(self, organism: str, max_results: int = 50) -> Optional[List[Dict[str, Any]]]:
        """
        Search COCONUT by organism/taxonomy

        Args:
            organism: Organism name or taxonomy filter
            max_results: Maximum number of results

        Returns:
            List of compounds from this organism or None on error
        """
        url = f"{self.BASE_URL}/search/organism"

        params = {
            "organism": organism,
            "limit": max_results
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("results", [])
                        logger.info(f"COCONUT: Found {len(results)} compounds for organism: {organism}")
                        return results
                    else:
                        logger.warning(f"COCONUT: Organism search failed with status {response.status}")
                        return []
        except Exception as e:
            logger.error(f"COCONUT: Error searching by organism {organism}: {e}")
            return None

    async def _get_compound_by_id(self, compound_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed compound information by COCONUT ID

        Args:
            compound_id: COCONUT compound ID (e.g., "CNP0123456")

        Returns:
            Compound details dictionary or None on error
        """
        url = f"{self.BASE_URL}/compound/{compound_id}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        logger.info(f"COCONUT: Retrieved compound {compound_id}")
                        return data
                    elif response.status == 404:
                        logger.warning(f"COCONUT: Compound {compound_id} not found")
                        return None
                    else:
                        logger.error(f"COCONUT: Failed to get compound {compound_id}, status {response.status}")
                        return None
        except Exception as e:
            logger.error(f"COCONUT: Error retrieving compound {compound_id}: {e}")
            return None

    def _apply_property_filters(self, compounds: List[Dict[str, Any]], filters: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Apply molecular property filters to compound list

        Args:
            compounds: List of compound dictionaries
            filters: Filter criteria (molecular_weight, logp, etc.)

        Returns:
            Filtered list of compounds
        """
        if not filters:
            return compounds

        filtered = []

        for compound in compounds:
            props = compound.get("properties", {})

            # Molecular weight filter
            mw_filter = filters.get("molecular_weight", {})
            if mw_filter:
                mw = props.get("molecular_weight")
                if mw:
                    if mw < mw_filter.get("min", 0) or mw > mw_filter.get("max", float('inf')):
                        continue

            # LogP filter
            logp_filter = filters.get("logp", {})
            if logp_filter:
                logp = props.get("logp")
                if logp:
                    if logp < logp_filter.get("min", -float('inf')) or logp > logp_filter.get("max", float('inf')):
                        continue

            # Organism filter
            organism_filter = filters.get("organism")
            if organism_filter:
                source = compound.get("natural_source", {})
                organism = source.get("organism", "")
                if organism_filter.lower() not in organism.lower():
                    continue

            filtered.append(compound)

        return filtered

    def _format_compound_result(self, raw_compound: Dict[str, Any], include_taxonomy: bool = True,
                               include_activities: bool = True) -> Dict[str, Any]:
        """
        Format raw COCONUT compound data into standardized output

        Args:
            raw_compound: Raw compound data from API
            include_taxonomy: Whether to include taxonomy information
            include_activities: Whether to include biological activities

        Returns:
            Formatted compound dictionary
        """
        result = {
            "coconut_id": raw_compound.get("id", raw_compound.get("coconut_id")),
            "name": raw_compound.get("name", raw_compound.get("iupac_name")),
            "smiles": raw_compound.get("smiles", raw_compound.get("canonical_smiles")),
            "inchi": raw_compound.get("inchi"),
            "inchikey": raw_compound.get("inchikey"),
            "molecular_formula": raw_compound.get("molecular_formula"),
            "molecular_weight": raw_compound.get("molecular_weight"),
            "properties": {
                "logp": raw_compound.get("xlogp", raw_compound.get("alogp")),
                "hba": raw_compound.get("hba", raw_compound.get("h_bond_acceptors")),
                "hbd": raw_compound.get("hbd", raw_compound.get("h_bond_donors")),
                "rotatable_bonds": raw_compound.get("rotatable_bonds"),
                "tpsa": raw_compound.get("tpsa")
            }
        }

        # Add natural source information
        if include_taxonomy and "natural_source" in raw_compound:
            result["natural_source"] = raw_compound["natural_source"]
        elif "organism" in raw_compound or "taxonomy" in raw_compound:
            result["natural_source"] = {
                "organism": raw_compound.get("organism"),
                "common_name": raw_compound.get("common_name"),
                "taxonomy": raw_compound.get("taxonomy")
            }

        # Add biological activities
        if include_activities and "biological_activities" in raw_compound:
            result["biological_activities"] = raw_compound["biological_activities"]
        elif "activities" in raw_compound:
            result["biological_activities"] = raw_compound["activities"]

        # Add references
        if "references" in raw_compound:
            result["references"] = raw_compound["references"]

        return result

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute COCONUT search

        Args:
            input_data: SMILES string or query dictionary
            **kwargs: Additional parameters (override config)

        Returns:
            AdapterResult containing natural product search results
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input. Expected SMILES string or query dictionary with 'query_type' and 'query'"
            )

        # Parse input into normalized query
        query = self._parse_input(input_data)
        query_type = query["query_type"]
        query_value = query["query"]
        filters = query.get("filters", {})
        max_results = kwargs.get("max_results", query.get("max_results", 50))
        include_taxonomy = query.get("include_taxonomy", True)
        include_activities = query.get("include_activities", True)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Execute appropriate search based on query type
        raw_results = None

        if query_type == "smiles":
            raw_results = await self._search_by_smiles(query_value, max_results)
        elif query_type == "name":
            raw_results = await self._search_by_name(query_value, max_results)
        elif query_type == "organism":
            raw_results = await self._search_by_organism(query_value, max_results)
        elif query_type == "id":
            single_result = await self._get_compound_by_id(query_value)
            raw_results = [single_result] if single_result else []
        elif query_type == "inchi" or query_type == "inchikey":
            # For InChI/InChIKey, treat similar to SMILES search
            # This would need proper endpoint implementation
            logger.warning(f"COCONUT: {query_type} search not fully implemented, falling back to name search")
            raw_results = await self._search_by_name(query_value, max_results)
        else:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported query type: {query_type}"
            )

        if raw_results is None:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to search COCONUT database for {query_type}: {query_value}",
                metadata={
                    "source": "coconut",
                    "query_type": query_type,
                    "query": query_value
                }
            )

        # Apply property filters if specified
        if filters:
            raw_results = self._apply_property_filters(raw_results, filters)

        # Format results
        formatted_results = [
            self._format_compound_result(compound, include_taxonomy, include_activities)
            for compound in raw_results
        ]

        # Build result summary
        result_data = {
            "results": formatted_results,
            "num_results": len(formatted_results),
            "query_type": query_type,
            "query": query_value,
            "warnings": []
        }

        # Add warnings if needed
        if len(formatted_results) == 0:
            result_data["warnings"].append(f"No natural products found for {query_type}: {query_value}")

        if len(raw_results) > max_results:
            result_data["warnings"].append(f"Results truncated to {max_results} compounds")

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "coconut",
                "adapter_version": self.version,
                "query_type": query_type,
                "num_results": len(formatted_results),
                "filters_applied": bool(filters)
            }
        )
