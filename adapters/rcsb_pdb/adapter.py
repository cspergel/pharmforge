"""
RCSB PDB Adapter for PharmForge

Retrieves experimental protein structures from the Research Collaboratory
for Structural Bioinformatics Protein Data Bank.

API Documentation: https://data.rcsb.org/
Search API: https://search.rcsb.org/
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


class RCSBPDBAdapter(AdapterProtocol):
    """
    RCSB PDB adapter for retrieving experimental protein structures.

    The RCSB PDB is the primary repository for experimentally-determined
    3D structures of proteins, nucleic acids, and complex assemblies.

    Features:
    - Download structures by PDB ID
    - Search structures by protein name
    - Extract ligand information
    - Binding site analysis
    - Quality metrics (resolution, R-factors)
    - Multiple file formats (PDB, mmCIF, PDBML)

    Quality Metrics:
    - Resolution: Lower is better (< 2.0 Å is high quality)
    - R-value/R-free: Measure of model quality (< 0.25 is good)
    - Experimental method: X-ray, NMR, Cryo-EM, etc.
    """

    DATA_API_URL = "https://data.rcsb.org/rest/v1/core"
    FILES_URL = "https://files.rcsb.org/download"
    SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

    def __init__(
        self,
        name: str = "rcsb_pdb",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize RCSB PDB adapter.

        Args:
            name: Adapter name (default: "rcsb_pdb")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - cache_dir: Directory for caching files (default: ./cache/rcsb_pdb)
                   - download_pdb: Download PDB format (default: True)
                   - download_cif: Download mmCIF format (default: False)
                   - include_ligands: Extract ligand info (default: True)
                   - include_binding_sites: Extract binding sites (default: True)
                   - timeout: Request timeout in seconds (default: 30)
        """
        default_config = {
            "cache_dir": "./cache/rcsb_pdb",
            "download_pdb": True,
            "download_cif": False,
            "include_ligands": True,
            "include_binding_sites": True,
            "timeout": 30
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Create cache directory
        self.cache_dir = Path(self.config["cache_dir"])
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"RCSB PDB cache directory: {self.cache_dir}")

    def validate_input(self, pdb_id: str) -> bool:
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

    async def execute(self, pdb_id: str, **params) -> AdapterResult:
        """
        Retrieve PDB structure and metadata.

        Args:
            pdb_id: PDB identifier (e.g., "1ABC", "6LU7")
            **params: Additional parameters:
                     - download_pdb: Override config
                     - download_cif: Override config
                     - include_ligands: Override config
                     - include_binding_sites: Override config

        Returns:
            AdapterResult containing:
                - structure_data: PDB/CIF file content
                - metadata: Structure metadata
                - quality_metrics: Resolution, R-factors, etc.
                - ligands: Ligand information
                - binding_sites: Binding site data
                - file_paths: Paths to cached files
        """
        try:
            # Validate and normalize input
            if not self.validate_input(pdb_id):
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"Invalid PDB ID format: {pdb_id}"
                )

            pdb_id = pdb_id.strip().upper()

            # Get parameters
            download_pdb = params.get("download_pdb", self.config["download_pdb"])
            download_cif = params.get("download_cif", self.config["download_cif"])
            include_ligands = params.get("include_ligands", self.config["include_ligands"])
            include_binding_sites = params.get("include_binding_sites", self.config["include_binding_sites"])

            # Generate cache key
            cache_key = self.generate_cache_key(pdb_id, **params)

            logger.info(f"Fetching RCSB PDB structure: {pdb_id}")

            # Step 1: Get structure metadata
            metadata = await self._get_structure_metadata(pdb_id)

            if not metadata:
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
            quality_metrics = self._extract_quality_metrics(metadata)

            # Step 4: Extract ligand information
            ligands = []
            if include_ligands:
                ligands = await self._get_ligand_info(pdb_id)

            # Step 5: Extract binding sites
            binding_sites = []
            if include_binding_sites and structure_files.get("pdb"):
                binding_sites = self._extract_binding_sites(structure_files["pdb"])

            # Construct result
            result_data = {
                "pdb_id": pdb_id,
                "metadata": metadata,
                "structure_files": structure_files,
                "file_paths": file_paths,
                "quality_metrics": quality_metrics,
                "ligands": ligands,
                "binding_sites": binding_sites,
                "experimental_method": metadata.get("experimental_method", "unknown"),
                "resolution": quality_metrics.get("resolution", None),
                "organism": metadata.get("organism", "unknown"),
                "reference": "RCSB Protein Data Bank"
            }

            logger.info(f"Successfully retrieved PDB structure: {pdb_id}")
            logger.info(f"  Method: {metadata.get('experimental_method', 'unknown')}")
            logger.info(f"  Resolution: {quality_metrics.get('resolution', 'N/A')} Å")
            logger.info(f"  Ligands: {len(ligands)}")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "cache_key": cache_key,
                    "version": self.version
                }
            )

        except Exception as e:
            logger.error(f"RCSB PDB retrieval failed for {pdb_id}: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    async def search_structures(
        self,
        query: str,
        max_results: int = 10,
        **filters
    ) -> AdapterResult:
        """
        Search for PDB structures by protein name or keywords.

        Args:
            query: Search query (protein name, keyword, etc.)
            max_results: Maximum number of results to return
            **filters: Additional filters:
                      - resolution_min: Minimum resolution in Å
                      - resolution_max: Maximum resolution in Å
                      - method: Experimental method ("X-RAY", "NMR", "EM")
                      - has_ligand: Only structures with ligands

        Returns:
            AdapterResult with list of matching PDB IDs and metadata
        """
        try:
            # Construct search query
            search_query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "value": query
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "results_content_type": ["experimental"],
                    "sort": [{"sort_by": "score", "direction": "desc"}],
                    "scoring_strategy": "combined"
                }
            }

            # Add filters
            if filters:
                search_query["query"] = {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [search_query["query"]]
                }

                # Resolution filter
                if "resolution_min" in filters or "resolution_max" in filters:
                    res_filter = {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "range"
                        }
                    }
                    if "resolution_min" in filters:
                        res_filter["parameters"]["value"] = {
                            "from": filters["resolution_min"]
                        }
                    if "resolution_max" in filters:
                        if "value" not in res_filter["parameters"]:
                            res_filter["parameters"]["value"] = {}
                        res_filter["parameters"]["value"]["to"] = filters["resolution_max"]

                    search_query["query"]["nodes"].append(res_filter)

            # Execute search
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.post(
                    self.SEARCH_API_URL,
                    json=search_query,
                    params={"limit": max_results}
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("result_set", [])

                        # Extract PDB IDs
                        pdb_ids = [result["identifier"] for result in results]

                        logger.info(f"Search found {len(pdb_ids)} results for: {query}")

                        return AdapterResult(
                            success=True,
                            data={
                                "query": query,
                                "pdb_ids": pdb_ids,
                                "total_results": len(pdb_ids),
                                "results": results
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

    async def _get_structure_metadata(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Get structure metadata from RCSB API.

        Args:
            pdb_id: PDB identifier

        Returns:
            Metadata dictionary or None if not found
        """
        url = f"{self.DATA_API_URL}/entry/{pdb_id}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        return await response.json()
                    elif response.status == 404:
                        logger.warning(f"PDB ID not found: {pdb_id}")
                        return None
                    else:
                        logger.error(f"RCSB API error: {response.status}")
                        return None
        except Exception as e:
            logger.error(f"Failed to fetch metadata: {e}")
            return None

    async def _download_structure(
        self,
        pdb_id: str,
        format: str = "pdb"
    ) -> Tuple[Optional[str], Optional[Path]]:
        """
        Download structure file from RCSB.

        Args:
            pdb_id: PDB identifier
            format: File format ("pdb" or "cif")

        Returns:
            Tuple of (file content, cache file path)
        """
        # Construct filename and URL
        if format == "pdb":
            filename = f"{pdb_id}.pdb"
            url = f"{self.FILES_URL}/{pdb_id}.pdb"
        else:
            filename = f"{pdb_id}.cif"
            url = f"{self.FILES_URL}/{pdb_id}.cif"

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

    async def _get_ligand_info(self, pdb_id: str) -> List[Dict[str, Any]]:
        """
        Get ligand information for a structure.

        Args:
            pdb_id: PDB identifier

        Returns:
            List of ligand dictionaries
        """
        url = f"{self.DATA_API_URL}/nonpolymer_entity/{pdb_id}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()

                        # Extract ligand information
                        ligands = []
                        if isinstance(data, dict):
                            data = [data]

                        for entity in data:
                            ligand_info = {
                                "name": entity.get("rcsb_nonpolymer_entity", {}).get("pdbx_description", "Unknown"),
                                "formula": entity.get("rcsb_nonpolymer_entity", {}).get("formula_weight", None),
                                "type": entity.get("rcsb_nonpolymer_entity", {}).get("type", "Unknown")
                            }

                            # Try to get SMILES if available
                            chem_comp = entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
                            if "smiles" in chem_comp:
                                ligand_info["smiles"] = chem_comp["smiles"]

                            ligands.append(ligand_info)

                        return ligands
                    else:
                        return []
        except Exception as e:
            logger.warning(f"Failed to get ligand info: {e}")
            return []

    def _extract_binding_sites(self, pdb_content: str) -> List[Dict[str, Any]]:
        """
        Extract binding site information from PDB file.

        Args:
            pdb_content: PDB file content

        Returns:
            List of binding site dictionaries
        """
        try:
            binding_sites = []
            current_site = None

            for line in pdb_content.split('\n'):
                if line.startswith('SITE'):
                    # SITE records define binding sites
                    site_id = line[11:14].strip()

                    if current_site is None or current_site["id"] != site_id:
                        if current_site:
                            binding_sites.append(current_site)

                        current_site = {
                            "id": site_id,
                            "residues": []
                        }

                    # Parse residues in this site
                    # Residues are in positions 18-21, 23-26, 28-31, etc.
                    for i in range(18, 62, 11):
                        res_name = line[i:i+3].strip()
                        if res_name:
                            current_site["residues"].append(res_name)

            if current_site:
                binding_sites.append(current_site)

            return binding_sites

        except Exception as e:
            logger.warning(f"Failed to extract binding sites: {e}")
            return []

    def _extract_quality_metrics(self, metadata: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract quality metrics from metadata.

        Args:
            metadata: Structure metadata

        Returns:
            Dictionary of quality metrics
        """
        metrics = {}

        # Experimental method
        exptl = metadata.get("exptl", [])
        if exptl and isinstance(exptl, list) and len(exptl) > 0:
            metrics["experimental_method"] = exptl[0].get("method", "Unknown")

        # Resolution
        refine = metadata.get("refine", [])
        if refine and isinstance(refine, list) and len(refine) > 0:
            resolution = refine[0].get("ls_d_res_high")
            if resolution:
                metrics["resolution"] = round(float(resolution), 2)

            # R-factors
            r_value = refine[0].get("ls_R_factor_R_work")
            if r_value:
                metrics["r_value"] = round(float(r_value), 4)

            r_free = refine[0].get("ls_R_factor_R_free")
            if r_free:
                metrics["r_free"] = round(float(r_free), 4)

        # Release date
        rcsb_accession = metadata.get("rcsb_accession_info", {})
        if "deposit_date" in rcsb_accession:
            metrics["deposit_date"] = rcsb_accession["deposit_date"]
        if "initial_release_date" in rcsb_accession:
            metrics["release_date"] = rcsb_accession["initial_release_date"]

        # Structure title
        struct = metadata.get("struct", {})
        if "title" in struct:
            metrics["title"] = struct["title"]

        return metrics

    def generate_cache_key(self, pdb_id: str, **kwargs) -> str:
        """
        Generate cache key for PDB requests.

        Args:
            pdb_id: PDB identifier
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "pdb_id": pdb_id.upper(),
            "download_pdb": kwargs.get("download_pdb", self.config["download_pdb"]),
            "download_cif": kwargs.get("download_cif", self.config["download_cif"]),
            "include_ligands": kwargs.get("include_ligands", self.config["include_ligands"]),
            "include_binding_sites": kwargs.get("include_binding_sites", self.config["include_binding_sites"])
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
            "description": "RCSB Protein Data Bank experimental structures adapter",
            "capabilities": {
                "structure_retrieval": True,
                "search": True,
                "ligand_info": True,
                "binding_sites": True,
                "quality_metrics": True,
                "multiple_formats": ["pdb", "cif", "pdbml"]
            },
            "database": {
                "name": "RCSB PDB",
                "url": "https://www.rcsb.org/",
                "structures": "200,000+",
                "experimental_methods": ["X-ray", "NMR", "Cryo-EM", "Electron diffraction"]
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
                "include_binding_sites": self.config["include_binding_sites"]
            },
            "reference": {
                "website": "https://www.rcsb.org/",
                "data_api": "https://data.rcsb.org/",
                "search_api": "https://search.rcsb.org/"
            }
        }
