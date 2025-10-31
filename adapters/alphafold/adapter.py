"""
AlphaFold DB Adapter for PharmForge

Retrieves predicted protein structures from the AlphaFold Protein Structure Database.
AlphaFold provides highly accurate protein structure predictions using deep learning.

API Documentation: https://alphafold.ebi.ac.uk/api-docs
Reference: Jumper et al., Nature 2021, 596, 583-589
"""

import hashlib
import logging
import os
from typing import Dict, Any, Optional, List
from pathlib import Path
import asyncio
import aiohttp
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class AlphaFoldAdapter(AdapterProtocol):
    """
    AlphaFold Database adapter for retrieving predicted protein structures.

    The AlphaFold DB provides access to protein structure predictions for
    hundreds of thousands of proteins from model organisms and humans.

    Features:
    - Download structures by UniProt accession
    - Access pLDDT confidence scores
    - Quality metrics (pAE, pTM)
    - Model version information
    - Local structure caching

    Quality Metrics:
    - pLDDT (predicted LDDT): Per-residue confidence (0-100, >90 is very high)
    - PAE (Predicted Aligned Error): Confidence in relative positions
    - pTM (predicted TM-score): Overall model confidence
    """

    BASE_URL = "https://alphafold.ebi.ac.uk/api"
    FILES_URL = "https://alphafold.ebi.ac.uk/files"

    def __init__(
        self,
        name: str = "alphafold",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize AlphaFold DB adapter.

        Args:
            name: Adapter name (default: "alphafold")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - cache_dir: Directory for caching structure files (default: ./cache/alphafold)
                   - download_pdb: Download PDB format (default: True)
                   - download_cif: Download mmCIF format (default: False)
                   - download_pae: Download PAE data (default: False)
                   - timeout: Request timeout in seconds (default: 30)
        """
        default_config = {
            "cache_dir": "./cache/alphafold",
            "download_pdb": True,
            "download_cif": False,
            "download_pae": False,
            "timeout": 30
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Create cache directory
        self.cache_dir = Path(self.config["cache_dir"])
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"AlphaFold cache directory: {self.cache_dir}")

    def validate_input(self, uniprot_id: str) -> bool:
        """
        Validate UniProt accession ID.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P12345")

        Returns:
            True if valid format, False otherwise
        """
        if not uniprot_id or not isinstance(uniprot_id, str):
            return False

        # Basic UniProt ID validation (6-10 alphanumeric characters)
        uniprot_id = uniprot_id.strip().upper()
        if len(uniprot_id) < 6 or len(uniprot_id) > 10:
            return False

        # Should start with letter or number
        if not uniprot_id[0].isalnum():
            return False

        return True

    async def execute(self, uniprot_id: str, **params) -> AdapterResult:
        """
        Retrieve AlphaFold structure prediction for a protein.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P12345", "Q9Y6K9")
            **params: Additional parameters:
                     - version: Specific AlphaFold version (default: latest)
                     - download_pdb: Override config (default: from config)
                     - download_cif: Override config (default: from config)
                     - download_pae: Override PAE download (default: from config)

        Returns:
            AdapterResult containing:
                - structure_data: PDB/CIF file content
                - confidence_scores: pLDDT, pTM scores
                - quality_metrics: Model quality assessment
                - file_paths: Paths to cached files
        """
        try:
            # Validate and normalize input
            if not self.validate_input(uniprot_id):
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"Invalid UniProt ID format: {uniprot_id}"
                )

            uniprot_id = uniprot_id.strip().upper()

            # Get parameters
            download_pdb = params.get("download_pdb", self.config["download_pdb"])
            download_cif = params.get("download_cif", self.config["download_cif"])
            download_pae = params.get("download_pae", self.config["download_pae"])
            version = params.get("version", None)

            # Generate cache key
            cache_key = self.generate_cache_key(uniprot_id, **params)

            logger.info(f"Fetching AlphaFold structure for UniProt ID: {uniprot_id}")

            # Step 1: Get structure metadata
            metadata = await self._get_structure_metadata(uniprot_id)

            if not metadata:
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"No AlphaFold prediction found for {uniprot_id}"
                )

            # Step 2: Download structure files
            structure_files = {}
            file_paths = {}

            # Get latest version from metadata if not specified
            if version is None:
                version = metadata.get("latestVersion")
                logger.debug(f"Using latest version from metadata: v{version}")

            if download_pdb:
                pdb_content, pdb_path = await self._download_structure(
                    uniprot_id, format="pdb", version=version
                )
                if pdb_content:
                    structure_files["pdb"] = pdb_content
                    file_paths["pdb"] = str(pdb_path)

            if download_cif:
                cif_content, cif_path = await self._download_structure(
                    uniprot_id, format="cif", version=version
                )
                if cif_content:
                    structure_files["cif"] = cif_content
                    file_paths["cif"] = str(cif_path)

            # Step 3: Download PAE data if requested
            pae_data = None
            if download_pae:
                pae_content, pae_path = await self._download_pae(uniprot_id, version=version)
                if pae_content:
                    pae_data = json.loads(pae_content)
                    file_paths["pae"] = str(pae_path)

            # Step 4: Extract quality metrics
            quality_metrics = self._extract_quality_metrics(metadata, pae_data)

            # Step 5: Parse pLDDT scores from PDB file
            plddt_scores = None
            if download_pdb and structure_files.get("pdb"):
                plddt_scores = self._extract_plddt_from_pdb(structure_files["pdb"])

            # Construct result
            result_data = {
                "uniprot_id": uniprot_id,
                "metadata": metadata,
                "structure_files": structure_files,
                "file_paths": file_paths,
                "quality_metrics": quality_metrics,
                "plddt_scores": plddt_scores,
                "model_version": metadata.get("modelVersion", "unknown"),
                "organism": metadata.get("organism", "unknown"),
                "reference": "Jumper et al., Nature 2021, 596, 583-589"
            }

            logger.info(f"Successfully retrieved AlphaFold structure for {uniprot_id}")
            logger.info(f"  Model version: {metadata.get('modelVersion', 'unknown')}")
            logger.info(f"  Average pLDDT: {quality_metrics.get('mean_plddt', 'N/A')}")

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
            logger.error(f"AlphaFold retrieval failed for {uniprot_id}: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    async def _get_structure_metadata(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get structure metadata from AlphaFold API.

        Args:
            uniprot_id: UniProt accession

        Returns:
            Metadata dictionary or None if not found
        """
        url = f"{self.BASE_URL}/prediction/{uniprot_id}"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data[0] if isinstance(data, list) and len(data) > 0 else data
                    elif response.status == 404:
                        logger.warning(f"No AlphaFold prediction found for {uniprot_id}")
                        return None
                    else:
                        logger.error(f"AlphaFold API error: {response.status}")
                        return None
        except Exception as e:
            logger.error(f"Failed to fetch metadata: {e}")
            return None

    async def _download_structure(
        self,
        uniprot_id: str,
        format: str = "pdb",
        version: Optional[int] = None
    ) -> tuple[Optional[str], Optional[Path]]:
        """
        Download structure file from AlphaFold.

        Args:
            uniprot_id: UniProt accession
            format: File format ("pdb" or "cif")
            version: Model version (default: tries v4, then falls back to v3, v2, v1)

        Returns:
            Tuple of (file content, cache file path)
        """
        # Try multiple versions if none specified (v6 is current as of Oct 2025, but some entries may use older versions)
        versions_to_try = [version] if version else [6, 5, 4, 3, 2, 1]

        for ver in versions_to_try:
            version_str = f"v{ver}"
            if format == "pdb":
                filename = f"AF-{uniprot_id}-F1-model_{version_str}.pdb"
            else:
                filename = f"AF-{uniprot_id}-F1-model_{version_str}.cif"

            # Check cache first
            cache_path = self.cache_dir / filename
            if cache_path.exists():
                logger.debug(f"Using cached file: {cache_path}")
                with open(cache_path, 'r') as f:
                    return f.read(), cache_path

            # Download from AlphaFold
            url = f"{self.FILES_URL}/{filename}"

            try:
                timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
                async with aiohttp.ClientSession(timeout=timeout) as session:
                    async with session.get(url) as response:
                        if response.status == 200:
                            content = await response.text()

                            # Save to cache
                            with open(cache_path, 'w') as f:
                                f.write(content)

                            logger.info(f"Downloaded and cached: {filename} (version {ver})")
                            return content, cache_path
                        elif response.status == 404:
                            logger.debug(f"Version {ver} not found for {uniprot_id}, trying next version...")
                            continue
                        else:
                            logger.warning(f"Download failed: {response.status} for {url}")
                            continue
            except Exception as e:
                logger.warning(f"Failed to download {filename}: {e}")
                continue

        # If we get here, all versions failed
        logger.error(f"Failed to download {format.upper()} for {uniprot_id} - tried versions: {versions_to_try}")
        return None, None

    async def _download_pae(
        self,
        uniprot_id: str,
        version: Optional[int] = None
    ) -> tuple[Optional[str], Optional[Path]]:
        """
        Download Predicted Aligned Error (PAE) data.

        Args:
            uniprot_id: UniProt accession
            version: Model version (default: tries v4, then falls back to v3, v2, v1)

        Returns:
            Tuple of (JSON content, cache file path)
        """
        # Try multiple versions if none specified
        versions_to_try = [version] if version else [6, 5, 4, 3, 2, 1]

        for ver in versions_to_try:
            version_str = f"v{ver}"
            filename = f"AF-{uniprot_id}-F1-predicted_aligned_error_{version_str}.json"

            # Check cache
            cache_path = self.cache_dir / filename
            if cache_path.exists():
                logger.debug(f"Using cached PAE file: {cache_path}")
                with open(cache_path, 'r') as f:
                    return f.read(), cache_path

            # Download
            url = f"{self.FILES_URL}/{filename}"

            try:
                timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
                async with aiohttp.ClientSession(timeout=timeout) as session:
                    async with session.get(url) as response:
                        if response.status == 200:
                            content = await response.text()

                            # Save to cache
                            with open(cache_path, 'w') as f:
                                f.write(content)

                            logger.info(f"Downloaded PAE data: {filename} (version {ver})")
                            return content, cache_path
                        elif response.status == 404:
                            logger.debug(f"PAE version {ver} not found for {uniprot_id}, trying next version...")
                            continue
                        else:
                            logger.warning(f"PAE download failed: {response.status} for {url}")
                            continue
            except Exception as e:
                logger.warning(f"Failed to download PAE {filename}: {e}")
                continue

        # If we get here, all versions failed
        logger.warning(f"PAE data not available for {uniprot_id} - tried versions: {versions_to_try}")
        return None, None

    def _extract_plddt_from_pdb(self, pdb_content: str) -> Dict[str, Any]:
        """
        Extract pLDDT scores from PDB file (B-factor column).

        Args:
            pdb_content: PDB file content

        Returns:
            Dictionary with pLDDT statistics
        """
        try:
            plddt_scores = []

            for line in pdb_content.split('\n'):
                if line.startswith('ATOM'):
                    # B-factor is in columns 61-66
                    try:
                        plddt = float(line[60:66].strip())
                        plddt_scores.append(plddt)
                    except ValueError:
                        continue

            if not plddt_scores:
                return {}

            # Calculate statistics
            import statistics

            return {
                "mean": round(statistics.mean(plddt_scores), 2),
                "median": round(statistics.median(plddt_scores), 2),
                "min": round(min(plddt_scores), 2),
                "max": round(max(plddt_scores), 2),
                "stdev": round(statistics.stdev(plddt_scores), 2) if len(plddt_scores) > 1 else 0.0,
                "num_residues": len(plddt_scores),
                "high_confidence_pct": round(
                    sum(1 for s in plddt_scores if s > 90) / len(plddt_scores) * 100, 2
                ),
                "medium_confidence_pct": round(
                    sum(1 for s in plddt_scores if 70 < s <= 90) / len(plddt_scores) * 100, 2
                ),
                "low_confidence_pct": round(
                    sum(1 for s in plddt_scores if s <= 70) / len(plddt_scores) * 100, 2
                )
            }
        except Exception as e:
            logger.warning(f"Failed to extract pLDDT scores: {e}")
            return {}

    def _extract_quality_metrics(
        self,
        metadata: Dict[str, Any],
        pae_data: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Extract and calculate quality metrics.

        Args:
            metadata: Structure metadata
            pae_data: Optional PAE data

        Returns:
            Dictionary of quality metrics
        """
        metrics = {
            "model_version": metadata.get("modelVersion", "unknown"),
            "created_date": metadata.get("modelCreatedDate", "unknown"),
        }

        # Extract global metric value (mean pLDDT)
        if "globalMetricValue" in metadata:
            metrics["mean_plddt"] = round(metadata["globalMetricValue"], 2)

        # Extract pLDDT confidence fractions
        if "fractionPlddtVeryHigh" in metadata:
            metrics["plddt_very_high_pct"] = round(metadata["fractionPlddtVeryHigh"] * 100, 2)
        if "fractionPlddtConfident" in metadata:
            metrics["plddt_confident_pct"] = round(metadata["fractionPlddtConfident"] * 100, 2)
        if "fractionPlddtLow" in metadata:
            metrics["plddt_low_pct"] = round(metadata["fractionPlddtLow"] * 100, 2)
        if "fractionPlddtVeryLow" in metadata:
            metrics["plddt_very_low_pct"] = round(metadata["fractionPlddtVeryLow"] * 100, 2)

        # Extract latest version info
        if "latestVersion" in metadata:
            metrics["latest_version"] = metadata["latestVersion"]

        # Extract pTM score if available
        if "ptm" in metadata:
            metrics["ptm_score"] = metadata["ptm"]

        # Add PAE metrics if available
        if pae_data:
            # Calculate mean PAE
            # Handle both list and dict formats (API may return list or dict)
            if isinstance(pae_data, list):
                pae_data = pae_data[0] if len(pae_data) > 0 else {}

            pae_values = pae_data.get("predicted_aligned_error", [])
            if pae_values and isinstance(pae_values[0], list):
                # Flatten nested list
                flat_pae = [val for row in pae_values for val in row]
                import statistics
                metrics["mean_pae"] = round(statistics.mean(flat_pae), 2)
                metrics["max_pae"] = round(max(flat_pae), 2)

        return metrics

    def generate_cache_key(self, uniprot_id: str, **kwargs) -> str:
        """
        Generate cache key for AlphaFold requests.

        Args:
            uniprot_id: UniProt accession
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "uniprot_id": uniprot_id.upper(),
            "download_pdb": kwargs.get("download_pdb", self.config["download_pdb"]),
            "download_cif": kwargs.get("download_cif", self.config["download_cif"]),
            "download_pae": kwargs.get("download_pae", self.config["download_pae"]),
            "model_version": kwargs.get("version", None)
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
            "description": "AlphaFold Protein Structure Database adapter",
            "capabilities": {
                "structure_prediction": True,
                "confidence_scores": True,
                "quality_metrics": True,
                "pae_data": True,
                "multiple_formats": ["pdb", "cif"]
            },
            "database": {
                "name": "AlphaFold DB",
                "url": "https://alphafold.ebi.ac.uk/",
                "coverage": "200M+ structures",
                "organisms": "Model organisms + Human proteome"
            },
            "quality_thresholds": {
                "plddt_high": "> 90 (very high confidence)",
                "plddt_medium": "70-90 (confident)",
                "plddt_low": "< 70 (low confidence)",
                "ptm_good": "> 0.5 (good model)"
            },
            "config": {
                "cache_dir": str(self.cache_dir),
                "download_pdb": self.config["download_pdb"],
                "download_cif": self.config["download_cif"],
                "download_pae": self.config["download_pae"]
            },
            "reference": {
                "paper": "Jumper et al., Nature 2021, 596, 583-589",
                "website": "https://alphafold.ebi.ac.uk/",
                "api_docs": "https://alphafold.ebi.ac.uk/api-docs"
            }
        }
