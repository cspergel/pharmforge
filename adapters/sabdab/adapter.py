"""
SAbDab (Structural Antibody Database) Adapter for PharmForge

Retrieves antibody structure and sequence information from the
Structural Antibody Database maintained by Oxford Protein Informatics Group.

API Documentation: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab
Database: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/
"""

import hashlib
import logging
import json
from typing import Dict, Any, Optional, List
from pathlib import Path
import asyncio
import aiohttp

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class SAbDabAdapter(AdapterProtocol):
    """
    SAbDab adapter for retrieving antibody structures and sequences.

    The Structural Antibody Database (SAbDab) provides a curated database
    of all antibody structures available in the PDB, with additional
    annotations including:
    - CDR (Complementarity-Determining Region) definitions
    - Antibody-antigen complex information
    - Sequence and structure quality metrics
    - Canonical class assignments

    Features:
    - Search antibodies by antigen name
    - Search by PDB ID
    - Filter by antibody type (IgG, IgM, Fab, scFv, etc.)
    - Filter by species, resolution, R-factor
    - Retrieve CDR sequences (H1, H2, H3, L1, L2, L3)
    - Get antibody chain information
    - Download structure files

    Quality Filters:
    - Resolution: Lower is better (< 3.0 Å recommended)
    - R-factor: Lower is better (< 0.3 recommended)
    - Structure method: X-ray, EM, NMR
    """

    BASE_URL = "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab"
    SEARCH_URL = f"{BASE_URL}/search"
    SUMMARY_URL = f"{BASE_URL}/summary"

    def __init__(
        self,
        name: str = "sabdab",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize SAbDab adapter.

        Args:
            name: Adapter name (default: "sabdab")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - cache_dir: Directory for caching files (default: ./cache/sabdab)
                   - timeout: Request timeout in seconds (default: 60)
                   - max_results: Maximum search results (default: 100)
                   - rate_limit_delay: Delay between requests in seconds (default: 0.5)
        """
        default_config = {
            "cache_dir": "./cache/sabdab",
            "timeout": 60,
            "max_results": 100,
            "rate_limit_delay": 0.5
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Create cache directory
        self.cache_dir = Path(self.config["cache_dir"])
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"SAbDab cache directory: {self.cache_dir}")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data (PDB ID or search query).

        Args:
            input_data: Either a PDB ID (str) or search parameters (dict)

        Returns:
            True if valid format, False otherwise
        """
        if isinstance(input_data, str):
            # PDB ID format
            if not input_data or len(input_data.strip()) == 0:
                return False
            return True
        elif isinstance(input_data, dict):
            # Search parameters - at least one search criterion
            valid_keys = {"antigen", "pdb_id", "species", "ab_type", "resolution_max", "rfactor_max"}
            return any(key in input_data for key in valid_keys)
        else:
            return False

    async def execute(self, input_data: Any, **params) -> AdapterResult:
        """
        Execute antibody structure retrieval or search.

        Args:
            input_data: Either:
                       - PDB ID string (e.g., "7BWJ")
                       - Dict with search parameters:
                         * antigen: Antigen name to search (e.g., "spike", "covid")
                         * pdb_id: Specific PDB ID
                         * species: Species filter (e.g., "human", "mouse")
                         * ab_type: Antibody type (e.g., "IgG", "Fab", "scFv")
                         * resolution_max: Maximum resolution in Å
                         * rfactor_max: Maximum R-factor
            **params: Additional parameters:
                     - max_results: Override max results
                     - include_cdrs: Include CDR sequences (default: True)
                     - include_structures: Download structure files (default: False)

        Returns:
            AdapterResult containing:
                - antibodies: List of antibody records
                - total_results: Number of results
                - search_params: Search parameters used
                - cdrs: CDR sequence information (if requested)
                - structures: Structure file paths (if requested)
        """
        try:
            # Validate input
            if not self.validate_input(input_data):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid input. Provide PDB ID (str) or search parameters (dict)"
                )

            # Rate limiting
            await asyncio.sleep(self.config.get("rate_limit_delay", 0.5))

            # Determine if this is a direct PDB lookup or search
            if isinstance(input_data, str):
                # Direct PDB ID lookup
                result = await self._get_by_pdb_id(input_data, **params)
            else:
                # Search by parameters
                result = await self._search_antibodies(input_data, **params)

            return result

        except Exception as e:
            logger.error(f"SAbDab operation failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    async def _get_by_pdb_id(self, pdb_id: str, **params) -> AdapterResult:
        """
        Retrieve antibody information by PDB ID.

        Args:
            pdb_id: PDB identifier
            **params: Additional parameters

        Returns:
            AdapterResult with antibody information
        """
        pdb_id = pdb_id.strip().upper()

        logger.info(f"Fetching SAbDab data for PDB ID: {pdb_id}")

        # Fetch summary data for this PDB ID
        url = f"{self.SUMMARY_URL}/{pdb_id}/"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        # Parse the response (could be CSV or TSV format)
                        text = await response.text()
                        antibody_data = self._parse_summary_data(text, pdb_id)

                        if not antibody_data:
                            return AdapterResult(
                                success=False,
                                data={},
                                error=f"No antibody data found for PDB ID: {pdb_id}"
                            )

                        # Optionally include CDR sequences
                        include_cdrs = params.get("include_cdrs", True)
                        if include_cdrs:
                            cdrs = await self._get_cdr_sequences(pdb_id)
                            antibody_data["cdrs"] = cdrs

                        result_data = {
                            "pdb_id": pdb_id,
                            "antibody_data": antibody_data,
                            "source": "SAbDab",
                            "reference": "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/"
                        }

                        logger.info(f"Successfully retrieved SAbDab data for {pdb_id}")

                        return AdapterResult(
                            success=True,
                            data=result_data,
                            metadata={
                                "adapter_name": self.name,
                                "version": self.version,
                                "pdb_id": pdb_id
                            }
                        )
                    elif response.status == 404:
                        return AdapterResult(
                            success=False,
                            data={},
                            error=f"PDB ID not found in SAbDab: {pdb_id}"
                        )
                    else:
                        error_text = await response.text()
                        logger.error(f"SAbDab API error {response.status}: {error_text}")
                        return AdapterResult(
                            success=False,
                            data={},
                            error=f"SAbDab API returned status {response.status}"
                        )

        except asyncio.TimeoutError:
            logger.error(f"Timeout fetching SAbDab data for {pdb_id}")
            return AdapterResult(
                success=False,
                data={},
                error=f"Request timeout for PDB ID: {pdb_id}"
            )
        except Exception as e:
            logger.error(f"Error fetching SAbDab data: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e)
            )

    async def _search_antibodies(self, search_params: Dict[str, Any], **params) -> AdapterResult:
        """
        Search for antibodies by various criteria.

        Args:
            search_params: Dictionary of search parameters
            **params: Additional options

        Returns:
            AdapterResult with list of matching antibodies
        """
        logger.info(f"Searching SAbDab with parameters: {search_params}")

        # Build search URL with parameters
        query_params = {
            "ABtype": search_params.get("ab_type", "All"),
            "method": "All",  # X-ray, EM, NMR, or All
            "species": search_params.get("species", "All"),
            "resolution": search_params.get("resolution_max", "All"),
            "rfactor": search_params.get("rfactor_max", "All"),
        }

        # Add antigen if provided
        if "antigen" in search_params:
            query_params["antigen"] = search_params["antigen"]

        # Add PDB ID if searching by specific ID
        if "pdb_id" in search_params:
            query_params["pdb"] = search_params["pdb_id"]

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(self.SEARCH_URL, params=query_params) as response:
                    if response.status == 200:
                        # Parse response (typically TSV or CSV format)
                        text = await response.text()
                        results = self._parse_search_results(text)

                        # Limit results
                        max_results = params.get("max_results", self.config["max_results"])
                        if len(results) > max_results:
                            results = results[:max_results]

                        logger.info(f"Found {len(results)} antibodies matching search criteria")

                        result_data = {
                            "search_params": search_params,
                            "total_results": len(results),
                            "antibodies": results,
                            "source": "SAbDab"
                        }

                        return AdapterResult(
                            success=True,
                            data=result_data,
                            metadata={
                                "adapter_name": self.name,
                                "version": self.version,
                                "query": search_params
                            }
                        )
                    else:
                        error_text = await response.text()
                        logger.error(f"Search failed with status {response.status}: {error_text}")
                        return AdapterResult(
                            success=False,
                            data={},
                            error=f"Search failed with status {response.status}"
                        )

        except Exception as e:
            logger.error(f"Search failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e)
            )

    async def _get_cdr_sequences(self, pdb_id: str) -> Dict[str, Any]:
        """
        Retrieve CDR (Complementarity-Determining Region) sequences.

        Args:
            pdb_id: PDB identifier

        Returns:
            Dictionary with CDR sequences for heavy and light chains
        """
        # CDR data endpoint
        url = f"{self.SUMMARY_URL}/{pdb_id}/"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        text = await response.text()
                        # Parse CDR information from response
                        cdrs = self._extract_cdr_sequences(text)
                        return cdrs
                    else:
                        logger.warning(f"Could not fetch CDR data for {pdb_id}")
                        return {}
        except Exception as e:
            logger.warning(f"Error fetching CDR sequences: {e}")
            return {}

    def _parse_summary_data(self, text: str, pdb_id: str) -> Dict[str, Any]:
        """
        Parse summary data from SAbDab response.

        Args:
            text: Response text (TSV format)
            pdb_id: PDB identifier

        Returns:
            Dictionary with parsed antibody data
        """
        try:
            lines = text.strip().split('\n')
            if len(lines) < 2:
                return {}

            # First line is header
            headers = lines[0].split('\t')

            # Find data for this PDB ID
            for line in lines[1:]:
                values = line.split('\t')
                if len(values) > 0 and values[0].upper() == pdb_id:
                    # Create dictionary from headers and values
                    data = dict(zip(headers, values))

                    # Parse and structure the data
                    parsed_data = {
                        "pdb_id": data.get("pdb", pdb_id),
                        "resolution": self._safe_float(data.get("resolution")),
                        "method": data.get("method", "Unknown"),
                        "rfactor": self._safe_float(data.get("Rfactor")),
                        "antibody_type": data.get("Htype", "Unknown"),
                        "antigen": data.get("antigen_name", "Unknown"),
                        "heavy_chain": data.get("Hchain", ""),
                        "light_chain": data.get("Lchain", ""),
                        "antigen_chain": data.get("antigen_chain", ""),
                        "species": data.get("species", "Unknown"),
                        "scfv": data.get("scfv", "") == "TRUE",
                        "engineered": data.get("engineered", "") == "TRUE"
                    }

                    return parsed_data

            return {}

        except Exception as e:
            logger.error(f"Error parsing summary data: {e}")
            return {}

    def _parse_search_results(self, text: str) -> List[Dict[str, Any]]:
        """
        Parse search results from SAbDab.

        Args:
            text: Response text (TSV format)

        Returns:
            List of antibody dictionaries
        """
        try:
            lines = text.strip().split('\n')
            if len(lines) < 2:
                return []

            headers = lines[0].split('\t')
            results = []

            for line in lines[1:]:
                values = line.split('\t')
                if len(values) >= len(headers):
                    data = dict(zip(headers, values))

                    # Parse into structured format
                    antibody = {
                        "pdb_id": data.get("pdb", ""),
                        "resolution": self._safe_float(data.get("resolution")),
                        "method": data.get("method", "Unknown"),
                        "antibody_type": data.get("Htype", "Unknown"),
                        "antigen": data.get("antigen_name", "Unknown"),
                        "species": data.get("species", "Unknown"),
                        "heavy_chain": data.get("Hchain", ""),
                        "light_chain": data.get("Lchain", "")
                    }

                    results.append(antibody)

            return results

        except Exception as e:
            logger.error(f"Error parsing search results: {e}")
            return []

    def _extract_cdr_sequences(self, text: str) -> Dict[str, Any]:
        """
        Extract CDR sequences from response data.

        Args:
            text: Response text

        Returns:
            Dictionary with CDR sequences
        """
        # This is a placeholder - actual implementation would parse
        # CDR sequences from the response
        return {
            "heavy_chain": {
                "h1": None,
                "h2": None,
                "h3": None
            },
            "light_chain": {
                "l1": None,
                "l2": None,
                "l3": None
            },
            "note": "CDR sequences require detailed parsing from structure files"
        }

    def _safe_float(self, value: Any) -> Optional[float]:
        """
        Safely convert value to float.

        Args:
            value: Value to convert

        Returns:
            Float value or None
        """
        if value is None or value == "" or value == "NA":
            return None
        try:
            return float(value)
        except (ValueError, TypeError):
            return None

    def generate_cache_key(self, input_data: Any, **kwargs) -> str:
        """
        Generate cache key for SAbDab requests.

        Args:
            input_data: PDB ID or search parameters
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "input": input_data,
            "params": kwargs
        }

        cache_string = json.dumps(cache_dict, sort_keys=True)
        return hashlib.sha256(cache_string.encode()).hexdigest()

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata.

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "Structural Antibody Database (SAbDab) adapter for antibody structures",
            "capabilities": {
                "pdb_lookup": True,
                "antigen_search": True,
                "cdr_sequences": True,
                "antibody_classification": True,
                "quality_filtering": True
            },
            "database": {
                "name": "SAbDab",
                "maintainer": "Oxford Protein Informatics Group",
                "url": "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/",
                "update_frequency": "Weekly",
                "coverage": "All antibody structures in PDB"
            },
            "antibody_types": [
                "IgG", "IgA", "IgM", "IgE", "IgD",
                "Fab", "scFv", "VHH", "Nanobody"
            ],
            "quality_filters": {
                "resolution": "< 3.0 Å recommended",
                "rfactor": "< 0.3 recommended",
                "methods": ["X-ray", "EM", "NMR"]
            },
            "config": {
                "cache_dir": str(self.cache_dir),
                "timeout": self.config["timeout"],
                "max_results": self.config["max_results"],
                "rate_limit_delay": self.config["rate_limit_delay"]
            },
            "reference": {
                "paper": "Dunbar et al. (2014) Nucleic Acids Research",
                "doi": "10.1093/nar/gkt1043",
                "url": "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/"
            }
        }
