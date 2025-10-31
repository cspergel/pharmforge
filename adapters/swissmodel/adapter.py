"""
SWISS-MODEL Adapter for PharmForge

Retrieves homology models and quality assessments from the SWISS-MODEL Repository.
SWISS-MODEL provides automated protein structure homology modeling and quality
assessment tools.

Repository: https://swissmodel.expasy.org/repository
Reference: Waterhouse et al., Nucleic Acids Res. 2018, 46, D296-D303
"""

import hashlib
import logging
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import aiohttp
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class SwissModelAdapter(AdapterProtocol):
    """
    SWISS-MODEL Repository adapter for homology models.

    SWISS-MODEL provides automated comparative protein structure modeling
    and a repository of annotated 3D protein structure models.

    Features:
    - Download homology models by UniProt ID
    - Structure quality assessment (QMEAN, GMQE)
    - Template information
    - Sequence coverage analysis
    - Model confidence scores

    Quality Metrics:
    - QMEAN: Quality Model Energy Analysis (-4 to 0, higher is better)
    - GMQE: Global Model Quality Estimation (0-1, higher is better)
    - Sequence identity: Percentage identity to template
    """

    REPOSITORY_URL = "https://swissmodel.expasy.org/repository"
    API_URL = "https://swissmodel.expasy.org/repository/uniprot"

    def __init__(
        self,
        name: str = "swissmodel",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize SWISS-MODEL adapter.

        Args:
            name: Adapter name (default: "swissmodel")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - cache_dir: Directory for caching files (default: ./cache/swissmodel)
                   - download_pdb: Download PDB models (default: True)
                   - include_templates: Include template info (default: True)
                   - min_qmean: Minimum QMEAN score (default: -4.0)
                   - min_gmqe: Minimum GMQE score (default: 0.0)
                   - timeout: Request timeout in seconds (default: 30)
        """
        default_config = {
            "cache_dir": "./cache/swissmodel",
            "download_pdb": True,
            "include_templates": True,
            "min_qmean": -4.0,
            "min_gmqe": 0.0,
            "timeout": 30
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Create cache directory
        self.cache_dir = Path(self.config["cache_dir"])
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"SWISS-MODEL cache directory: {self.cache_dir}")

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
        Retrieve SWISS-MODEL homology models for a protein.

        Args:
            uniprot_id: UniProt accession ID (e.g., "P12345", "Q9Y6K9")
            **params: Additional parameters:
                     - download_pdb: Override config
                     - include_templates: Override config
                     - min_qmean: Minimum QMEAN score filter
                     - min_gmqe: Minimum GMQE score filter

        Returns:
            AdapterResult containing:
                - models: List of available homology models
                - best_model: Highest quality model
                - structure_data: PDB file content for best model
                - quality_metrics: QMEAN, GMQE scores
                - template_info: Template structure information
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
            include_templates = params.get("include_templates", self.config["include_templates"])
            min_qmean = params.get("min_qmean", self.config["min_qmean"])
            min_gmqe = params.get("min_gmqe", self.config["min_gmqe"])

            # Generate cache key
            cache_key = self.generate_cache_key(uniprot_id, **params)

            logger.info(f"Fetching SWISS-MODEL structures for UniProt ID: {uniprot_id}")

            # Step 1: Get available models
            models = await self._get_models(uniprot_id)

            if not models:
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"No SWISS-MODEL models found for {uniprot_id}"
                )

            # Step 2: Filter models by quality
            # Note: QMEAN is a nested dict with multiple scores; use qmean4_z_score as the main metric
            filtered_models = []
            for m in models:
                # Extract QMEAN score from nested structure
                qmean_value = self._extract_qmean_score(m.get("qmean"))
                gmqe_value = m.get("gmqe", 0)

                # Filter based on quality thresholds
                if qmean_value >= min_qmean and gmqe_value >= min_gmqe:
                    # Store extracted values for later use
                    m["_qmean_score"] = qmean_value
                    m["_gmqe_score"] = gmqe_value
                    filtered_models.append(m)

            if not filtered_models:
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"No models meet quality criteria (QMEAN>={min_qmean}, GMQE>={min_gmqe})"
                )

            # Step 3: Select best model (highest GMQE)
            best_model = max(filtered_models, key=lambda m: m.get("_gmqe_score", 0))

            logger.info(f"Found {len(filtered_models)} models, best GMQE: {best_model.get('_gmqe_score', 'N/A')}")

            # Step 4: Download structure files
            structure_files = {}
            file_paths = {}

            if download_pdb and "coordinates" in best_model:
                pdb_content, pdb_path = await self._download_structure(
                    best_model["coordinates"], uniprot_id
                )
                if pdb_content:
                    structure_files["pdb"] = pdb_content
                    file_paths["pdb"] = str(pdb_path)

            # Step 5: Get template information
            template_info = None
            if include_templates and "template" in best_model:
                template_info = best_model["template"]

            # Step 6: Extract quality metrics
            quality_metrics = self._extract_quality_metrics(best_model)

            # Construct result
            result_data = {
                "uniprot_id": uniprot_id,
                "models": filtered_models,
                "best_model": best_model,
                "structure_files": structure_files,
                "file_paths": file_paths,
                "quality_metrics": quality_metrics,
                "template_info": template_info,
                "model_count": len(filtered_models),
                "reference": "Waterhouse et al., Nucleic Acids Res. 2018, 46, D296-D303"
            }

            logger.info(f"Successfully retrieved SWISS-MODEL structures for {uniprot_id}")
            logger.info(f"  Best model QMEAN: {quality_metrics.get('qmean', 'N/A')}")
            logger.info(f"  Best model GMQE: {quality_metrics.get('gmqe', 'N/A')}")

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
            logger.error(f"SWISS-MODEL retrieval failed for {uniprot_id}: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    async def assess_sequence(self, sequence: str) -> AdapterResult:
        """
        Assess modeling feasibility for a protein sequence.

        Args:
            sequence: Protein sequence (FASTA format or raw sequence)

        Returns:
            AdapterResult with template search results and feasibility assessment
        """
        try:
            # Clean sequence (remove FASTA header if present)
            if sequence.startswith('>'):
                sequence = '\n'.join(sequence.split('\n')[1:])
            sequence = sequence.replace('\n', '').replace(' ', '').upper()

            logger.info(f"Assessing modeling feasibility for sequence ({len(sequence)} residues)")

            # For now, return a placeholder result
            # The actual SWISS-MODEL API for sequence submission requires
            # authentication and job submission, which is more complex

            result_data = {
                "sequence_length": len(sequence),
                "sequence": sequence[:50] + "..." if len(sequence) > 50 else sequence,
                "assessment": "To build models from sequence, use SWISS-MODEL web interface",
                "repository_search_recommended": True,
                "note": "This adapter focuses on retrieving pre-built models from the repository"
            }

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={"adapter_name": self.name}
            )

        except Exception as e:
            logger.error(f"Sequence assessment failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e)
            )

    async def _get_models(self, uniprot_id: str) -> List[Dict[str, Any]]:
        """
        Get available models from SWISS-MODEL Repository.

        Args:
            uniprot_id: UniProt accession

        Returns:
            List of model dictionaries
        """
        url = f"{self.API_URL}/{uniprot_id}.json"

        try:
            timeout = aiohttp.ClientTimeout(total=self.config["timeout"])
            async with aiohttp.ClientSession(timeout=timeout) as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        # Extract structures from response
                        structures = data.get("result", {}).get("structures", [])
                        return structures
                    elif response.status == 404:
                        logger.warning(f"No SWISS-MODEL models found for {uniprot_id}")
                        return []
                    else:
                        logger.error(f"SWISS-MODEL API error: {response.status}")
                        return []
        except Exception as e:
            logger.error(f"Failed to fetch models: {e}")
            return []

    async def _download_structure(
        self,
        coordinates_url: str,
        uniprot_id: str
    ) -> Tuple[Optional[str], Optional[Path]]:
        """
        Download structure file from SWISS-MODEL.

        Args:
            coordinates_url: URL to PDB coordinates
            uniprot_id: UniProt ID for naming

        Returns:
            Tuple of (file content, cache file path)
        """
        # Generate filename
        filename = f"{uniprot_id}_swissmodel.pdb"

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
                async with session.get(coordinates_url) as response:
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

    def _extract_qmean_score(self, qmean_data: Any) -> float:
        """
        Extract QMEAN score from potentially nested structure.

        QMEAN can be either:
        - A float/int (for older API versions or PDB structures)
        - A dict with multiple QMEAN scores (for SWISSMODEL homology models)

        Args:
            qmean_data: QMEAN data (float, int, dict, or None)

        Returns:
            Float QMEAN score (uses qmean4_z_score if dict, -10.0 if missing)
        """
        try:
            if qmean_data is None:
                return -10.0

            # If it's already a number, return it
            if isinstance(qmean_data, (int, float)):
                return float(qmean_data)

            # If it's a dict, extract the main QMEAN score
            if isinstance(qmean_data, dict):
                # Prefer qmean4_z_score as the primary metric (range: typically -4 to 0)
                if "qmean4_z_score" in qmean_data:
                    return float(qmean_data["qmean4_z_score"])
                # Fall back to qmean6_z_score if available
                elif "qmean6_z_score" in qmean_data:
                    return float(qmean_data["qmean6_z_score"])
                # Fall back to normalized scores (range: 0 to 1, convert to z-score approximation)
                elif "qmean4_norm_score" in qmean_data:
                    # Approximate z-score from normalized score (not perfect but reasonable)
                    norm_score = float(qmean_data["qmean4_norm_score"])
                    return (norm_score - 0.5) * 8 - 4  # Rough conversion
                else:
                    logger.warning(f"QMEAN dict missing expected scores: {qmean_data.keys()}")
                    return -10.0

            # Unknown type
            logger.warning(f"Unexpected QMEAN data type: {type(qmean_data)}")
            return -10.0

        except (ValueError, TypeError, KeyError) as e:
            logger.error(f"Error extracting QMEAN score: {e}")
            return -10.0

    def _extract_quality_metrics(self, model: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract quality metrics from model data.

        Args:
            model: Model dictionary

        Returns:
            Structured quality metrics
        """
        # Extract QMEAN score (handle nested structure)
        qmean_data = model.get("qmean")
        qmean_score = self._extract_qmean_score(qmean_data)

        metrics = {
            "qmean": qmean_score,
            "qmean_raw": qmean_data,  # Store raw data for reference
            "gmqe": model.get("gmqe"),
            "sequence_identity": model.get("identity"),
            "sequence_similarity": model.get("similarity"),
            "coverage": model.get("coverage"),
            "oligo_state": model.get("oligo-state", "monomer"),
            "method": model.get("method", "homology_modeling")
        }

        # Add range information
        if "from" in model and "to" in model:
            metrics["model_range"] = {
                "start": model["from"],
                "end": model["to"],
                "length": model["to"] - model["from"] + 1
            }

        # Interpret quality
        qmean = metrics.get("qmean")
        gmqe = metrics.get("gmqe")

        quality_assessment = []
        if qmean is not None and qmean > -10.0:  # -10.0 is our "missing" indicator
            if qmean >= -1.5:
                quality_assessment.append("Very high quality (QMEAN)")
            elif qmean >= -2.5:
                quality_assessment.append("High quality (QMEAN)")
            elif qmean >= -3.5:
                quality_assessment.append("Medium quality (QMEAN)")
            else:
                quality_assessment.append("Low quality (QMEAN)")

        if gmqe is not None and gmqe > 0:
            if gmqe >= 0.7:
                quality_assessment.append("High reliability (GMQE)")
            elif gmqe >= 0.4:
                quality_assessment.append("Medium reliability (GMQE)")
            else:
                quality_assessment.append("Low reliability (GMQE)")

        metrics["quality_assessment"] = quality_assessment

        return metrics

    def generate_cache_key(self, uniprot_id: str, **kwargs) -> str:
        """
        Generate cache key for SWISS-MODEL requests.

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
            "min_qmean": kwargs.get("min_qmean", self.config["min_qmean"]),
            "min_gmqe": kwargs.get("min_gmqe", self.config["min_gmqe"])
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
            "description": "SWISS-MODEL homology modeling repository adapter",
            "capabilities": {
                "homology_models": True,
                "quality_assessment": True,
                "template_info": True,
                "sequence_coverage": True,
                "quality_filtering": True
            },
            "database": {
                "name": "SWISS-MODEL Repository",
                "url": "https://swissmodel.expasy.org/repository",
                "coverage": "800,000+ models covering UniProtKB",
                "update_frequency": "Quarterly"
            },
            "quality_metrics": {
                "qmean": {
                    "range": "-4 to 0",
                    "interpretation": "Higher is better",
                    "high_quality": "> -1.5",
                    "medium_quality": "-2.5 to -1.5",
                    "low_quality": "< -2.5"
                },
                "gmqe": {
                    "range": "0 to 1",
                    "interpretation": "Higher is better",
                    "high_reliability": "> 0.7",
                    "medium_reliability": "0.4 to 0.7",
                    "low_reliability": "< 0.4"
                }
            },
            "config": {
                "cache_dir": str(self.cache_dir),
                "download_pdb": self.config["download_pdb"],
                "min_qmean": self.config["min_qmean"],
                "min_gmqe": self.config["min_gmqe"]
            },
            "reference": {
                "paper": "Waterhouse et al., Nucleic Acids Res. 2018, 46, D296-D303",
                "website": "https://swissmodel.expasy.org/",
                "documentation": "https://swissmodel.expasy.org/docs/repository_help"
            },
            "use_cases": [
                "Get homology models when experimental structures unavailable",
                "Quality assessment of protein structure models",
                "Template-based modeling feasibility analysis",
                "Comparative modeling for drug discovery"
            ]
        }
