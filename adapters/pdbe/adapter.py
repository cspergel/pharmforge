"""
PDBe (Protein Data Bank in Europe) Adapter for PharmForge

Retrieves protein structures and metadata from the European mirror of the
Protein Data Bank, providing alternative access and additional annotations.

API Documentation: https://www.ebi.ac.uk/pdbe/api/doc/
"""

import hashlib
import logging
import os
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import asyncio
import aiohttp
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class PDBEAdapter(AdapterProtocol):
    """
    PDBe adapter for retrieving protein structures from the European PDB.

    The Protein Data Bank in Europe (PDBe) is the European node of the
    worldwide Protein Data Bank. It provides the same structures as RCSB
    but may have additional annotations and faster access for European users.

    Features:
    - Download structures by PDB ID
    - Search structures by keywords, organism, compound
    - Extract ligand and binding information
    - Get functional annotations
    - Quality metrics (resolution, R-factors)
    - Multiple file formats (PDB, mmCIF)
    - European data center for faster access
    - Additional European annotations

    Quality Metrics:
    - Resolution: Lower is better (< 2.0 Å is high quality)
    - R-value/R-free: Measure of model quality (< 0.25 is good)
    - Experimental method: X-ray, NMR, Cryo-EM, etc.
    """

    # PDBe API endpoints
    API_BASE = "https://www.ebi.ac.uk/pdbe/api"
    FILES_BASE = "https://www.ebi.ac.uk/pdbe/entry-files/download"

    # Specific endpoints
    SUMMARY_ENDPOINT = f"{API_BASE}/pdb/entry/summary"
    MOLECULES_ENDPOINT = f"{API_BASE}/pdb/entry/molecules"
    LIGANDS_ENDPOINT = f"{API_BASE}/pdb/entry/ligand_monomers"
    BINDING_SITES_ENDPOINT = f"{API_BASE}/pdb/entry/binding_sites"
    ANNOTATIONS_ENDPOINT = f"{API_BASE}/graph-api/pdb/funpdbe_annotation"
    SEARCH_ENDPOINT = f"{API_BASE}/search/pdb/select"

    def __init__(
        self,
        name: str = "pdbe",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize PDBe adapter.

        Args:
            name: Adapter name (default: "pdbe")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - cache_dir: Directory for caching files (default: ./cache/pdbe)
                   - download_pdb: Download PDB format (default: True)
                   - download_cif: Download mmCIF format (default: False)
                   - include_ligands: Extract ligand info (default: True)
                   - include_binding_sites: Extract binding sites (default: True)
                   - include_annotations: Fetch functional annotations (default: False)
                   - timeout: Request timeout in seconds (default: 30)
                   - rate_limit_delay: Delay between requests in seconds (default: 0.5)
        """
        default_config = {
            "cache_dir": "./cache/pdbe",
            "download_pdb": True,
            "download_cif": False,
            "include_ligands": True,
            "include_binding_sites": True,
            "include_annotations": False,
            "timeout": 30,
            "rate_limit_delay": 0.5
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Create cache directory
        self.cache_dir = Path(self.config["cache_dir"])
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"PDBe cache directory: {self.cache_dir}")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data.

        Args:
            input_data: Either a PDB ID string or a dictionary with query parameters

        Returns:
            True if valid format, False otherwise
        """
        # Support both PDB ID strings and query dictionaries
        if isinstance(input_data, str):
            return self._validate_pdb_id(input_data)
        elif isinstance(input_data, dict):
            return self._validate_query_dict(input_data)
        else:
            return False

    def _validate_pdb_id(self, pdb_id: str) -> bool:
        """
        Validate PDB ID format.

        Args:
            pdb_id: PDB identifier (e.g., "1ABC")

        Returns:
            True if valid format, False otherwise
        """
        if not pdb_id or not isinstance(pdb_id, str):
            return False

        # PDB IDs are exactly 4 characters: 1 digit + 3 alphanumeric
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) != 4:
            return False

        # First character must be a digit
        if not pdb_id[0].isdigit():
            return False

        # Remaining characters must be alphanumeric
        if not pdb_id[1:].isalnum():
            return False

        return True

    def _validate_query_dict(self, query_dict: Dict) -> bool:
        """
        Validate query dictionary format.

        Args:
            query_dict: Dictionary with query parameters

        Returns:
            True if valid format, False otherwise
        """
        if "query_type" not in query_dict or "query" not in query_dict:
            return False

        valid_query_types = ["pdb_id", "keyword", "organism", "compound"]
        if query_dict["query_type"] not in valid_query_types:
            return False

        if not isinstance(query_dict["query"], str) or len(query_dict["query"]) == 0:
            return False

        return True

    async def execute(self, input_data: Any, **params) -> AdapterResult:
        """
        Retrieve PDB structure and metadata from PDBe.

        Args:
            input_data: Either PDB ID string (e.g., "1ABC", "6LU7") or query dictionary:
                       {
                           "query_type": "pdb_id"|"keyword"|"organism"|"compound",
                           "query": "search term",
                           "max_results": 10
                       }
            **params: Additional parameters:
                     - download_pdb: Override config
                     - download_cif: Override config
                     - include_ligands: Override config
                     - include_binding_sites: Override config
                     - include_annotations: Override config
                     - max_results: Maximum search results (default: 10)

        Returns:
            AdapterResult containing:
                - For PDB ID query:
                    - structure_data: PDB/CIF file content
                    - metadata: Structure metadata
                    - quality_metrics: Resolution, R-factors, etc.
                    - ligands: Ligand information
                    - binding_sites: Binding site data
                    - annotations: Functional annotations (if requested)
                    - file_paths: Paths to cached files
                - For search queries:
                    - entries: List of matching PDB entries with metadata
                    - num_results: Number of results found
        """
        try:
            # Validate input
            if not self.validate_input(input_data):
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"Invalid input format: {input_data}"
                )

            # Rate limiting
            await asyncio.sleep(self.config["rate_limit_delay"])

            # Handle different input types
            if isinstance(input_data, str):
                # Single PDB ID query
                return await self._execute_pdb_id(input_data, **params)
            elif isinstance(input_data, dict):
                # Search query
                return await self._execute_search(input_data, **params)
            else:
                return AdapterResult(
                    success=False,
                    data={},
                    error="Input must be PDB ID string or query dictionary"
                )

        except Exception as e:
            logger.error(f"PDBe adapter execution failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    async def _execute_pdb_id(self, pdb_id: str, **params) -> AdapterResult:
        """
        Execute PDB ID retrieval.

        Args:
            pdb_id: PDB identifier
            **params: Additional parameters

        Returns:
            AdapterResult with structure data
        """
        # Validate and normalize
        pdb_id = pdb_id.strip().upper()

        # Get parameters
        download_pdb = params.get("download_pdb", self.config["download_pdb"])
        download_cif = params.get("download_cif", self.config["download_cif"])
        include_ligands = params.get("include_ligands", self.config["include_ligands"])
        include_binding_sites = params.get("include_binding_sites", self.config["include_binding_sites"])
        include_annotations = params.get("include_annotations", self.config["include_annotations"])

        logger.info(f"Fetching PDBe structure: {pdb_id}")

        # Step 1: Get structure summary
        summary = await self._get_summary(pdb_id)

        if not summary:
            return AdapterResult(
                success=False,
                data={},
                error=f"PDB ID not found: {pdb_id}"
            )

        # Step 2: Download structure files
        structure_files = {}
        file_paths = {}

        if download_pdb:
            pdb_content, pdb_path = await self._download_structure(pdb_id, format="pdb")
            if pdb_content:
                structure_files["pdb"] = pdb_content
                file_paths["pdb"] = str(pdb_path)

        if download_cif:
            cif_content, cif_path = await self._download_structure(pdb_id, format="cif")
            if cif_content:
                structure_files["cif"] = cif_content
                file_paths["cif"] = str(cif_path)

        # Step 3: Extract quality metrics
        quality_metrics = self._extract_quality_metrics(summary)

        # Step 4: Get molecule information
        molecules = await self._get_molecules(pdb_id)

        # Step 5: Extract ligand information
        ligands = []
        if include_ligands:
            ligands = await self._get_ligands(pdb_id)

        # Step 6: Extract binding sites
        binding_sites = []
        if include_binding_sites:
            binding_sites = await self._get_binding_sites(pdb_id)

        # Step 7: Get functional annotations (optional)
        annotations = {}
        if include_annotations:
            annotations = await self._get_annotations(pdb_id)

        # Construct result
        result_data = {
            "pdb_id": pdb_id,
            "summary": summary,
            "structure_files": structure_files,
            "file_paths": file_paths,
            "quality_metrics": quality_metrics,
            "molecules": molecules,
            "ligands": ligands,
            "binding_sites": binding_sites,
            "annotations": annotations,
            "experimental_method": quality_metrics.get("experimental_method", "unknown"),
            "resolution": quality_metrics.get("resolution", None),
            "organism": self._extract_organism(summary),
            "num_chains": len(molecules) if isinstance(molecules, list) else (len(molecules.get("entities", [])) if molecules else 0),
            "num_residues": self._count_residues(summary),
            "reference": "Protein Data Bank in Europe (PDBe)"
        }

        logger.info(f"Successfully retrieved PDBe structure: {pdb_id}")
        logger.info(f"  Method: {quality_metrics.get('experimental_method', 'unknown')}")
        logger.info(f"  Resolution: {quality_metrics.get('resolution', 'N/A')} Å")
        logger.info(f"  Ligands: {len(ligands)}")
        logger.info(f"  Binding sites: {len(binding_sites)}")

        return AdapterResult(
            success=True,
            data=result_data,
            metadata={
                "adapter_name": self.name,
                "version": self.version,
                "source": "PDBe Europe"
            }
        )

    async def _execute_search(self, query_dict: Dict, **params) -> AdapterResult:
        """
        Execute search query.

        Args:
            query_dict: Query parameters
            **params: Additional parameters

        Returns:
            AdapterResult with search results
        """
        query_type = query_dict["query_type"]
        query = query_dict["query"]
        max_results = query_dict.get("max_results", params.get("max_results", 10))

        logger.info(f"Searching PDBe: {query_type} = {query}")

        # Build search URL based on query type
        if query_type == "pdb_id":
            # Direct PDB ID retrieval
            return await self._execute_pdb_id(query, **params)

        # For other query types, use the search API
        search_params = {
            "q": query,
            "rows": max_results,
            "wt": "json"
        }

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(self.SEARCH_ENDPOINT, params=search_params) as response:
                    if response.status == 200:
                        data = await response.json()
                        docs = data.get("response", {}).get("docs", [])

                        entries = []
                        for doc in docs[:max_results]:
                            entry = {
                                "pdb_id": doc.get("pdb_id", ""),
                                "title": doc.get("title", ""),
                                "organism": doc.get("organism_scientific_name", [""])[0] if doc.get("organism_scientific_name") else "",
                                "method": doc.get("experimental_method", [""])[0] if doc.get("experimental_method") else "",
                                "resolution": doc.get("resolution", None)
                            }
                            entries.append(entry)

                        logger.info(f"Found {len(entries)} results for: {query}")

                        return AdapterResult(
                            success=True,
                            data={
                                "query_type": query_type,
                                "query": query,
                                "entries": entries,
                                "num_results": len(entries),
                                "total_available": data.get("response", {}).get("numFound", 0)
                            },
                            metadata={
                                "adapter_name": self.name,
                                "version": self.version
                            }
                        )
                    else:
                        logger.error(f"Search failed: {response.status}")
                        return AdapterResult(
                            success=False,
                            data={},
                            error=f"Search API error: {response.status}"
                        )

        except Exception as e:
            logger.error(f"Search failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e)
            )

    async def _get_summary(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Get structure summary from PDBe API.

        Args:
            pdb_id: PDB identifier

        Returns:
            Summary dictionary or None if not found
        """
        # PDBe API uses lowercase PDB IDs
        pdb_id_lower = pdb_id.lower()
        url = f"{self.SUMMARY_ENDPOINT}/{pdb_id_lower}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        # PDBe returns data keyed by lowercase PDB ID
                        return data.get(pdb_id_lower, [{}])[0] if data.get(pdb_id_lower) else None
                    elif response.status == 404:
                        logger.warning(f"PDB ID not found: {pdb_id}")
                        return None
                    else:
                        logger.error(f"PDBe API error: {response.status}")
                        return None
        except Exception as e:
            logger.error(f"Failed to fetch summary: {e}")
            return None

    async def _get_molecules(self, pdb_id: str) -> Optional[Any]:
        """
        Get molecule information from PDBe API.

        Args:
            pdb_id: PDB identifier

        Returns:
            Molecule information (dict or list) or None
        """
        pdb_id_lower = pdb_id.lower()
        url = f"{self.MOLECULES_ENDPOINT}/{pdb_id_lower}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        # Return the data keyed by PDB ID (could be dict or list)
                        return data.get(pdb_id_lower, None)
                    else:
                        logger.warning(f"Failed to get molecules: {response.status}")
                        return None
        except Exception as e:
            logger.warning(f"Failed to get molecules: {e}")
            return None

    async def _get_ligands(self, pdb_id: str) -> List[Dict[str, Any]]:
        """
        Get ligand information from PDBe API.

        Args:
            pdb_id: PDB identifier

        Returns:
            List of ligand dictionaries
        """
        pdb_id_lower = pdb_id.lower()
        url = f"{self.LIGANDS_ENDPOINT}/{pdb_id_lower}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        ligand_data = data.get(pdb_id_lower, [])

                        ligands = []
                        for ligand in ligand_data:
                            ligand_info = {
                                "id": ligand.get("chem_comp_id", ""),
                                "name": ligand.get("chem_comp_name", ""),
                                "formula": ligand.get("formula", ""),
                                "type": ligand.get("type", ""),
                                "weight": ligand.get("formula_weight", None)
                            }
                            ligands.append(ligand_info)

                        return ligands
                    else:
                        logger.warning(f"Failed to get ligands: {response.status}")
                        return []
        except Exception as e:
            logger.warning(f"Failed to get ligands: {e}")
            return []

    async def _get_binding_sites(self, pdb_id: str) -> List[Dict[str, Any]]:
        """
        Get binding site information from PDBe API.

        Args:
            pdb_id: PDB identifier

        Returns:
            List of binding site dictionaries
        """
        pdb_id_lower = pdb_id.lower()
        url = f"{self.BINDING_SITES_ENDPOINT}/{pdb_id_lower}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        sites_data = data.get(pdb_id_lower, [])

                        binding_sites = []
                        for site in sites_data:
                            site_info = {
                                "site_id": site.get("site_id", ""),
                                "ligand_id": site.get("ligand_id", ""),
                                "residues": site.get("site_residues", []),
                                "num_residues": len(site.get("site_residues", []))
                            }
                            binding_sites.append(site_info)

                        return binding_sites
                    else:
                        logger.warning(f"Failed to get binding sites: {response.status}")
                        return []
        except Exception as e:
            logger.warning(f"Failed to get binding sites: {e}")
            return []

    async def _get_annotations(self, pdb_id: str) -> Dict[str, Any]:
        """
        Get functional annotations from PDBe API.

        Args:
            pdb_id: PDB identifier

        Returns:
            Annotations dictionary
        """
        pdb_id_lower = pdb_id.lower()
        url = f"{self.ANNOTATIONS_ENDPOINT}/{pdb_id_lower}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data.get(pdb_id_lower, {})
                    else:
                        logger.warning(f"Failed to get annotations: {response.status}")
                        return {}
        except Exception as e:
            logger.warning(f"Failed to get annotations: {e}")
            return {}

    async def _download_structure(
        self,
        pdb_id: str,
        format: str = "pdb"
    ) -> Tuple[Optional[str], Optional[Path]]:
        """
        Download structure file from PDBe.

        Args:
            pdb_id: PDB identifier
            format: File format ("pdb" or "cif")

        Returns:
            Tuple of (file content, cache file path)
        """
        # Construct filename and URL
        if format == "pdb":
            filename = f"pdb{pdb_id.lower()}.ent"
            url = f"{self.FILES_BASE}/pdb{pdb_id.lower()}.ent"
        else:
            filename = f"{pdb_id.lower()}.cif"
            url = f"{self.FILES_BASE}/{pdb_id.lower()}.cif"

        # Check cache
        cache_path = self.cache_dir / filename
        if cache_path.exists():
            logger.debug(f"Using cached file: {cache_path}")
            with open(cache_path, 'r') as f:
                return f.read(), cache_path

        # Download
        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        content = await response.text()

                        # Save to cache
                        with open(cache_path, 'w') as f:
                            f.write(content)

                        logger.info(f"Downloaded and cached: {filename}")
                        return content, cache_path
                    else:
                        logger.error(f"Download failed: {response.status}")
                        return None, None
        except Exception as e:
            logger.error(f"Failed to download structure: {e}")
            return None, None

    def _extract_quality_metrics(self, summary: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract quality metrics from summary data.

        Args:
            summary: Structure summary

        Returns:
            Dictionary of quality metrics
        """
        metrics = {}

        # Experimental method
        if "experimental_method" in summary:
            methods = summary["experimental_method"]
            if isinstance(methods, list) and len(methods) > 0:
                metrics["experimental_method"] = methods[0]
            else:
                metrics["experimental_method"] = str(methods)

        # Resolution
        if "resolution" in summary:
            try:
                metrics["resolution"] = round(float(summary["resolution"]), 2)
            except (ValueError, TypeError):
                metrics["resolution"] = None

        # R-factors
        if "r_factor" in summary:
            try:
                metrics["r_value"] = round(float(summary["r_factor"]), 4)
            except (ValueError, TypeError):
                pass

        if "r_free" in summary:
            try:
                metrics["r_free"] = round(float(summary["r_free"]), 4)
            except (ValueError, TypeError):
                pass

        # Dates
        if "deposition_date" in summary:
            metrics["deposit_date"] = summary["deposition_date"]

        if "release_date" in summary:
            metrics["release_date"] = summary["release_date"]

        # Title
        if "title" in summary:
            metrics["title"] = summary["title"]

        return metrics

    def _extract_organism(self, summary: Dict[str, Any]) -> str:
        """
        Extract organism information from summary.

        Args:
            summary: Structure summary

        Returns:
            Organism name
        """
        if "organism_scientific_name" in summary:
            organisms = summary["organism_scientific_name"]
            if isinstance(organisms, list) and len(organisms) > 0:
                return organisms[0]
            return str(organisms)
        return "unknown"

    def _count_residues(self, summary: Dict[str, Any]) -> int:
        """
        Count total number of residues in structure.

        Args:
            summary: Structure summary

        Returns:
            Total residue count
        """
        if "number_of_polymer_entities" in summary:
            try:
                return int(summary.get("number_of_entities", 0))
            except (ValueError, TypeError):
                pass
        return 0

    def generate_cache_key(self, input_data: Any, **kwargs) -> str:
        """
        Generate cache key for PDBe requests.

        Args:
            input_data: PDB ID or query dictionary
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        if isinstance(input_data, str):
            cache_dict = {
                "adapter": self.name,
                "version": self.version,
                "pdb_id": input_data.upper(),
                "download_pdb": kwargs.get("download_pdb", self.config["download_pdb"]),
                "download_cif": kwargs.get("download_cif", self.config["download_cif"]),
                "include_ligands": kwargs.get("include_ligands", self.config["include_ligands"]),
                "include_binding_sites": kwargs.get("include_binding_sites", self.config["include_binding_sites"]),
                "include_annotations": kwargs.get("include_annotations", self.config["include_annotations"])
            }
        else:
            cache_dict = {
                "adapter": self.name,
                "version": self.version,
                "query": input_data,
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
            "description": "Protein Data Bank in Europe (PDBe) - European mirror of PDB with additional annotations",
            "capabilities": {
                "structure_retrieval": True,
                "search": True,
                "ligand_info": True,
                "binding_sites": True,
                "quality_metrics": True,
                "functional_annotations": True,
                "multiple_formats": ["pdb", "cif"]
            },
            "database": {
                "name": "PDBe",
                "full_name": "Protein Data Bank in Europe",
                "url": "https://www.ebi.ac.uk/pdbe/",
                "location": "European Bioinformatics Institute (EMBL-EBI), UK",
                "structures": "200,000+ (mirrors RCSB PDB)",
                "experimental_methods": ["X-ray", "NMR", "Cryo-EM", "Electron diffraction"],
                "advantages": [
                    "European data center for faster access",
                    "Additional European annotations",
                    "FunPDBe functional annotations",
                    "Alternative to RCSB for redundancy"
                ]
            },
            "quality_thresholds": {
                "resolution_high": "< 2.0 Å (high quality)",
                "resolution_medium": "2.0-3.0 Å (medium quality)",
                "resolution_low": "> 3.0 Å (low quality)",
                "r_value_good": "< 0.25",
                "r_free_good": "< 0.28"
            },
            "config": {
                "cache_dir": str(self.cache_dir),
                "download_pdb": self.config["download_pdb"],
                "download_cif": self.config["download_cif"],
                "include_ligands": self.config["include_ligands"],
                "include_binding_sites": self.config["include_binding_sites"],
                "include_annotations": self.config["include_annotations"],
                "rate_limit_delay": self.config["rate_limit_delay"]
            },
            "reference": {
                "website": "https://www.ebi.ac.uk/pdbe/",
                "api_docs": "https://www.ebi.ac.uk/pdbe/api/doc/",
                "citation": "Protein Data Bank in Europe (PDBe)",
                "organization": "EMBL-EBI"
            }
        }
