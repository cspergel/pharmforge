"""
PDB-REDO Adapter for PharmForge

Retrieves re-refined and optimized protein structures from PDB-REDO.
PDB-REDO improves crystallographic structure models through automated
re-refinement and validation.

Website: https://pdb-redo.eu/
Reference: Joosten et al., IUCrJ 2014, 1, 213-220
"""

import hashlib
import logging
from typing import Dict, Any, Optional, Tuple
from pathlib import Path
import aiohttp
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class PDBRedoAdapter(AdapterProtocol):
    """
    PDB-REDO adapter for retrieving re-refined protein structures.

    PDB-REDO automatically re-refines crystallographic structures from the PDB,
    often resulting in improved geometry, fit to density, and overall quality.

    Features:
    - Download re-refined structures by PDB ID
    - Quality metrics comparison (before/after)
    - Validation reports
    - Re-refinement statistics
    - Multiple file formats

    Quality Improvements:
    - Better R-factors
    - Improved geometry (bonds, angles)
    - Better fit to electron density
    - Fixed outliers and errors
    """

    BASE_URL = "https://pdb-redo.eu"
    API_URL = f"{BASE_URL}/api"

    def __init__(
        self,
        name: str = "pdb_redo",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize PDB-REDO adapter.

        Args:
            name: Adapter name (default: "pdb_redo")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - cache_dir: Directory for caching files (default: ./cache/pdb_redo)
                   - download_pdb: Download refined PDB (default: True)
                   - download_cif: Download refined mmCIF (default: False)
                   - download_mtz: Download structure factors (default: False)
                   - include_validation: Include validation report (default: True)
                   - timeout: Request timeout in seconds (default: 60)
        """
        default_config = {
            "cache_dir": "./cache/pdb_redo",
            "download_pdb": True,
            "download_cif": False,
            "download_mtz": False,
            "include_validation": True,
            "timeout": 60  # PDB-REDO can be slow
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Create cache directory
        self.cache_dir = Path(self.config["cache_dir"])
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"PDB-REDO cache directory: {self.cache_dir}")

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
        pdb_id = pdb_id.strip().lower()  # PDB-REDO uses lowercase
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
        Retrieve PDB-REDO re-refined structure.

        Args:
            pdb_id: PDB identifier (e.g., "1ABC", "6LU7")
            **params: Additional parameters:
                     - download_pdb: Override config
                     - download_cif: Override config
                     - download_mtz: Override config
                     - include_validation: Override config

        Returns:
            AdapterResult containing:
                - structure_data: Refined PDB/CIF file content
                - quality_metrics: Before/after refinement metrics
                - validation_report: Quality validation data
                - improvements: List of improvements made
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

            pdb_id = pdb_id.strip().lower()  # PDB-REDO uses lowercase

            # Get parameters
            download_pdb = params.get("download_pdb", self.config["download_pdb"])
            download_cif = params.get("download_cif", self.config["download_cif"])
            download_mtz = params.get("download_mtz", self.config["download_mtz"])
            include_validation = params.get("include_validation", self.config["include_validation"])

            # Generate cache key
            cache_key = self.generate_cache_key(pdb_id, **params)

            logger.info(f"Fetching PDB-REDO re-refined structure: {pdb_id}")

            # Step 1: Check if structure exists in PDB-REDO
            exists = await self._check_structure_exists(pdb_id)
            if not exists:
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"PDB-REDO entry not found for: {pdb_id}"
                )

            # Step 2: Get quality metrics
            quality_metrics = await self._get_quality_metrics(pdb_id)

            # Step 3: Download structure files
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

            if download_mtz:
                mtz_content, mtz_path = await self._download_structure(pdb_id, format="mtz")
                if mtz_content:
                    structure_files["mtz"] = mtz_content
                    file_paths["mtz"] = str(mtz_path)

            # Step 4: Get validation report
            validation_report = None
            if include_validation:
                validation_report = await self._get_validation_report(pdb_id)

            # Step 5: Calculate improvements
            improvements = self._calculate_improvements(quality_metrics)

            # Construct result
            result_data = {
                "pdb_id": pdb_id.upper(),
                "structure_files": structure_files,
                "file_paths": file_paths,
                "quality_metrics": quality_metrics,
                "validation_report": validation_report,
                "improvements": improvements,
                "refinement_version": quality_metrics.get("version", "unknown"),
                "reference": "Joosten et al., IUCrJ 2014, 1, 213-220"
            }

            logger.info(f"Successfully retrieved PDB-REDO structure: {pdb_id}")
            if quality_metrics:
                logger.info(f"  R-work improvement: {improvements.get('r_work_improvement', 'N/A')}")
                logger.info(f"  R-free improvement: {improvements.get('r_free_improvement', 'N/A')}")

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
            logger.error(f"PDB-REDO retrieval failed for {pdb_id}: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    async def _check_structure_exists(self, pdb_id: str) -> bool:
        """
        Check if PDB-REDO entry exists.

        Args:
            pdb_id: PDB identifier (lowercase)

        Returns:
            True if exists, False otherwise
        """
        # PDB-REDO organizes files in /db/ directory
        url = f"{self.BASE_URL}/db/{pdb_id}/{pdb_id}_final.pdb"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.head(url) as response:
                    return response.status == 200
        except Exception as e:
            logger.error(f"Failed to check structure existence: {e}")
            return False

    async def _download_structure(
        self,
        pdb_id: str,
        format: str = "pdb"
    ) -> Tuple[Optional[str], Optional[Path]]:
        """
        Download re-refined structure file.

        Args:
            pdb_id: PDB identifier (lowercase)
            format: File format ("pdb", "cif", or "mtz")

        Returns:
            Tuple of (file content, cache file path)
        """
        # Construct filename and URL
        if format == "pdb":
            filename = f"{pdb_id}_final.pdb"
        elif format == "cif":
            filename = f"{pdb_id}_final.cif"
        elif format == "mtz":
            filename = f"{pdb_id}_final.mtz"
        else:
            return None, None

        url = f"{self.BASE_URL}/db/{pdb_id}/{filename}"

        # Check cache
        cache_path = self.cache_dir / filename
        if cache_path.exists():
            logger.debug(f"Using cached file: {cache_path}")
            mode = 'rb' if format == 'mtz' else 'r'
            with open(cache_path, mode) as f:
                content = f.read()
                return content if isinstance(content, str) else content.decode('utf-8'), cache_path

        # Download
        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        if format == 'mtz':
                            content = await response.read()
                            with open(cache_path, 'wb') as f:
                                f.write(content)
                            return content.decode('utf-8', errors='ignore'), cache_path
                        else:
                            content = await response.text()
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

    async def _get_quality_metrics(self, pdb_id: str) -> Dict[str, Any]:
        """
        Get quality metrics from PDB-REDO.

        Args:
            pdb_id: PDB identifier (lowercase)

        Returns:
            Dictionary with before/after metrics
        """
        url = f"{self.BASE_URL}/db/{pdb_id}/{pdb_id}_final.json"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        return self._parse_quality_metrics(data)
                    else:
                        logger.warning(f"Quality metrics not available: {response.status}")
                        return {}
        except Exception as e:
            logger.warning(f"Failed to get quality metrics: {e}")
            return {}

    def _parse_quality_metrics(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse quality metrics from JSON data.

        Args:
            data: Raw JSON data from PDB-REDO

        Returns:
            Structured quality metrics
        """
        metrics = {
            "version": data.get("version", "unknown"),
        }

        # Extract before/after R-factors
        if "rfactor" in data:
            rfactor = data["rfactor"]
            metrics["original"] = {
                "r_work": rfactor.get("rwork_original"),
                "r_free": rfactor.get("rfree_original")
            }
            metrics["refined"] = {
                "r_work": rfactor.get("rwork_final"),
                "r_free": rfactor.get("rfree_final")
            }

        # Extract geometry metrics
        if "geometry" in data:
            geom = data["geometry"]
            metrics["geometry"] = {
                "bond_rmsd": geom.get("bond_rmsd"),
                "angle_rmsd": geom.get("angle_rmsd"),
                "ramachandran_outliers": geom.get("rama_outliers")
            }

        return metrics

    async def _get_validation_report(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Get validation report.

        Args:
            pdb_id: PDB identifier (lowercase)

        Returns:
            Validation report dictionary
        """
        url = f"{self.BASE_URL}/db/{pdb_id}/{pdb_id}_validation.json"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        return await response.json()
                    else:
                        logger.warning(f"Validation report not available: {response.status}")
                        return None
        except Exception as e:
            logger.warning(f"Failed to get validation report: {e}")
            return None

    def _calculate_improvements(self, quality_metrics: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate improvements from re-refinement.

        Args:
            quality_metrics: Quality metrics with before/after values

        Returns:
            Dictionary of improvements
        """
        improvements = {}

        if "original" in quality_metrics and "refined" in quality_metrics:
            original = quality_metrics["original"]
            refined = quality_metrics["refined"]

            # R-work improvement
            if original.get("r_work") and refined.get("r_work"):
                r_work_delta = original["r_work"] - refined["r_work"]
                improvements["r_work_improvement"] = round(r_work_delta, 4)
                improvements["r_work_improvement_pct"] = round(
                    (r_work_delta / original["r_work"]) * 100, 2
                )

            # R-free improvement
            if original.get("r_free") and refined.get("r_free"):
                r_free_delta = original["r_free"] - refined["r_free"]
                improvements["r_free_improvement"] = round(r_free_delta, 4)
                improvements["r_free_improvement_pct"] = round(
                    (r_free_delta / original["r_free"]) * 100, 2
                )

        return improvements

    def generate_cache_key(self, pdb_id: str, **kwargs) -> str:
        """
        Generate cache key for PDB-REDO requests.

        Args:
            pdb_id: PDB identifier
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "pdb_id": pdb_id.lower(),
            "download_pdb": kwargs.get("download_pdb", self.config["download_pdb"]),
            "download_cif": kwargs.get("download_cif", self.config["download_cif"]),
            "download_mtz": kwargs.get("download_mtz", self.config["download_mtz"]),
            "include_validation": kwargs.get("include_validation", self.config["include_validation"])
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
            "description": "PDB-REDO re-refined protein structures adapter",
            "capabilities": {
                "re_refined_structures": True,
                "quality_improvement": True,
                "validation_reports": True,
                "structure_factors": True,
                "multiple_formats": ["pdb", "cif", "mtz"]
            },
            "database": {
                "name": "PDB-REDO",
                "url": "https://pdb-redo.eu/",
                "coverage": "180,000+ re-refined structures",
                "update_frequency": "Weekly"
            },
            "improvements": {
                "r_factors": "Typically 2-5% improvement",
                "geometry": "Better bond lengths and angles",
                "density_fit": "Improved fit to electron density",
                "outliers": "Reduced Ramachandran and rotamer outliers"
            },
            "config": {
                "cache_dir": str(self.cache_dir),
                "download_pdb": self.config["download_pdb"],
                "download_cif": self.config["download_cif"],
                "download_mtz": self.config["download_mtz"],
                "include_validation": self.config["include_validation"]
            },
            "reference": {
                "paper": "Joosten et al., IUCrJ 2014, 1, 213-220",
                "website": "https://pdb-redo.eu/",
                "documentation": "https://pdb-redo.eu/about"
            }
        }
