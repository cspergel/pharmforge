"""
DiffDock Adapter for PharmForge

Performs diffusion-based molecular docking using DiffDock to predict binding poses
with high accuracy and confidence scores. DiffDock uses a diffusion generative model
for blind docking without requiring a binding site specification.

References:
    - DiffDock Paper: https://arxiv.org/abs/2210.01776
    - GitHub: https://github.com/gcorso/DiffDock
"""

import os
import tempfile
import asyncio
import logging
import json
import shutil
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional, List, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - install with: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DiffDockAdapter(AdapterProtocol):
    """
    Adapter for DiffDock diffusion-based molecular docking

    DiffDock uses a diffusion model to perform blind docking, predicting both
    binding poses and confidence scores without requiring binding site specification.

    Features:
        - Blind docking (no binding site required)
        - Multiple pose predictions with confidence scores
        - State-of-the-art accuracy on benchmarks
        - Fast inference compared to traditional docking

    Requirements:
        - DiffDock installation (Python environment or Docker)
        - RDKit for molecule handling
        - Protein PDB file
        - GPU recommended for faster inference

    Input:
        - Protein PDB file path
        - Ligand SMILES string

    Output:
        - Docked poses (SDF format)
        - Confidence scores for each pose
        - Binding affinity estimates
        - Position and rotation predictions
    """

    def __init__(
        self,
        name: str = "diffdock",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize DiffDock adapter

        Args:
            name: Adapter name (default: "diffdock")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - diffdock_path: Path to DiffDock installation
                   - model_checkpoint: Path to model weights (optional)
                   - inference_steps: Number of diffusion steps (default: 20)
                   - samples_per_complex: Number of poses to generate (default: 40)
                   - batch_size: Batch size for inference (default: 10)
                   - use_gpu: Whether to use GPU (default: True)
                   - python_env: Path to Python environment with DiffDock
        """
        super().__init__(name, adapter_type, config or {})
        self.version = "1.0.0"

        # DiffDock configuration
        self.diffdock_path = self.config.get('diffdock_path', None)
        self.model_checkpoint = self.config.get('model_checkpoint', None)
        self.inference_steps = self.config.get('inference_steps', 20)
        self.samples_per_complex = self.config.get('samples_per_complex', 40)
        self.batch_size = self.config.get('batch_size', 10)
        self.use_gpu = self.config.get('use_gpu', True)
        self.python_env = self.config.get('python_env', 'python')

        # Validate configuration
        if not self.diffdock_path:
            logger.warning(
                "No diffdock_path configured. Set this to the DiffDock installation directory."
            )

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string

        Args:
            input_data: SMILES string or dict with 'smiles' key

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE:
            logger.error("RDKit not available for validation")
            return False

        # Extract SMILES from dict or use directly
        if isinstance(input_data, dict):
            smiles = input_data.get('smiles')
        else:
            smiles = input_data

        if not isinstance(smiles, str) or len(smiles) == 0:
            return False

        # Validate with RDKit
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.debug(f"SMILES validation failed: {e}")
            return False

    def _prepare_input_csv(
        self,
        protein_path: str,
        smiles: str,
        csv_path: str,
        complex_name: str = "complex"
    ) -> bool:
        """
        Prepare CSV input file for DiffDock

        DiffDock expects a CSV file with columns: complex_name, protein_path, ligand

        Args:
            protein_path: Path to protein PDB file
            smiles: Ligand SMILES string
            csv_path: Output CSV file path
            complex_name: Name for this complex (default: "complex")

        Returns:
            True if successful, False otherwise
        """
        try:
            with open(csv_path, 'w') as f:
                f.write("complex_name,protein_path,ligand_description,protein_sequence\n")
                # Absolute path for protein
                abs_protein_path = os.path.abspath(protein_path)
                f.write(f"{complex_name},{abs_protein_path},{smiles},\n")
            return True
        except Exception as e:
            logger.error(f"Error creating CSV input: {e}")
            return False

    async def _run_diffdock(
        self,
        csv_path: str,
        output_dir: str
    ) -> Tuple[bool, str]:
        """
        Run DiffDock inference

        Args:
            csv_path: Path to input CSV file
            output_dir: Directory for output files

        Returns:
            (success, error_message) tuple
        """
        try:
            if not self.diffdock_path:
                return False, "DiffDock path not configured"

            # Build DiffDock command
            cmd = [
                self.python_env,
                os.path.join(self.diffdock_path, "inference.py"),
                "--protein_ligand_csv", csv_path,
                "--out_dir", output_dir,
                "--inference_steps", str(self.inference_steps),
                "--samples_per_complex", str(self.samples_per_complex),
                "--batch_size", str(self.batch_size),
            ]

            # Add model checkpoint if specified
            if self.model_checkpoint:
                cmd.extend(["--model_dir", self.model_checkpoint])

            # GPU configuration
            if not self.use_gpu:
                cmd.append("--no_cuda")

            logger.info(f"Running DiffDock: {' '.join(cmd)}")

            # Run DiffDock as subprocess
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=self.diffdock_path
            )

            stdout, stderr = await process.communicate()

            if process.returncode != 0:
                error = stderr.decode() if stderr else "Unknown error"
                logger.error(f"DiffDock failed: {error}")
                return False, error

            stdout_text = stdout.decode() if stdout else ""
            logger.debug(f"DiffDock output: {stdout_text}")

            return True, ""

        except FileNotFoundError:
            error = f"DiffDock not found at: {self.diffdock_path}"
            logger.error(error)
            return False, error
        except Exception as e:
            logger.error(f"Error running DiffDock: {e}")
            return False, str(e)

    def _parse_diffdock_results(
        self,
        output_dir: str,
        complex_name: str = "complex"
    ) -> List[Dict[str, Any]]:
        """
        Parse DiffDock output files

        DiffDock generates:
        - SDF files with docked poses (rank1.sdf, rank2.sdf, ...)
        - confidence scores in filenames or separate file

        Args:
            output_dir: DiffDock output directory
            complex_name: Name of the complex

        Returns:
            List of docking results with poses and confidence scores
        """
        results = []

        try:
            # Look for output directory for this complex
            complex_dir = os.path.join(output_dir, complex_name)
            if not os.path.exists(complex_dir):
                logger.warning(f"Complex directory not found: {complex_dir}")
                return results

            # Find all SDF files (ranked poses)
            sdf_files = sorted(
                [f for f in os.listdir(complex_dir) if f.endswith('.sdf')],
                key=lambda x: int(x.replace('rank', '').replace('.sdf', ''))
            )

            for sdf_file in sdf_files:
                sdf_path = os.path.join(complex_dir, sdf_file)

                # Extract rank from filename
                rank = int(sdf_file.replace('rank', '').replace('.sdf', ''))

                # Read SDF to get confidence (often stored in properties)
                try:
                    suppl = Chem.SDMolSupplier(sdf_path)
                    mol = next(suppl)
                    if mol is not None:
                        # Try to extract confidence from properties
                        confidence = mol.GetPropsAsDict().get('confidence', None)

                        # If not in properties, parse from filename or separate file
                        if confidence is None:
                            # DiffDock may store confidence in a separate file
                            confidence_file = os.path.join(complex_dir, 'confidence.txt')
                            if os.path.exists(confidence_file):
                                with open(confidence_file, 'r') as f:
                                    lines = f.readlines()
                                    if rank <= len(lines):
                                        confidence = float(lines[rank - 1].strip())

                        results.append({
                            "rank": rank,
                            "confidence": float(confidence) if confidence else None,
                            "sdf_path": sdf_path,
                            # Could extract more properties here
                        })
                except Exception as e:
                    logger.warning(f"Error reading SDF {sdf_file}: {e}")

        except Exception as e:
            logger.error(f"Error parsing DiffDock results: {e}")

        return results

    def _estimate_binding_affinity(self, confidence: float) -> Optional[float]:
        """
        Estimate binding affinity from DiffDock confidence score

        DiffDock confidence scores correlate with binding affinity but are not
        direct measurements. This provides a rough estimate.

        Args:
            confidence: DiffDock confidence score (higher = better)

        Returns:
            Estimated binding affinity in kcal/mol (more negative = better)
        """
        if confidence is None:
            return None

        # Empirical mapping from confidence to binding affinity
        # This is a rough approximation - DiffDock confidence is not calibrated
        # Higher confidence -> more negative binding affinity
        # Typical range: confidence 0-5, affinity -12 to -4 kcal/mol
        estimated_affinity = -4.0 - (confidence * 1.5)
        return estimated_affinity

    def _fetch_pdb_structure(self, pdb_id: str, output_path: str) -> bool:
        """
        Fetch PDB structure from RCSB PDB database

        Args:
            pdb_id: 4-character PDB ID (e.g., "1HSG")
            output_path: Path to save PDB file

        Returns:
            True if successful, False otherwise
        """
        try:
            pdb_id = pdb_id.upper().strip()
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

            logger.info(f"Fetching PDB structure {pdb_id} from RCSB...")
            urllib.request.urlretrieve(url, output_path)

            if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
                logger.error(f"Failed to fetch PDB {pdb_id}")
                return False

            logger.info(f"Downloaded PDB structure: {pdb_id}")
            return True

        except Exception as e:
            logger.error(f"Error fetching PDB {pdb_id}: {e}")
            return False

    def _fetch_alphafold_structure(self, uniprot_id: str, output_path: str) -> bool:
        """
        Fetch AlphaFold predicted structure from AlphaFold DB

        Args:
            uniprot_id: UniProt ID (e.g., "P04637")
            output_path: Path to save PDB file

        Returns:
            True if successful, False otherwise
        """
        try:
            uniprot_id = uniprot_id.upper().strip()

            # Try different AlphaFold DB versions (v4, v3, v2)
            versions = ['v4', 'v3', 'v2']

            for version in versions:
                url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{version}.pdb"

                try:
                    logger.info(f"Trying AlphaFold structure {uniprot_id} ({version})...")
                    urllib.request.urlretrieve(url, output_path)

                    if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                        logger.info(f"Downloaded AlphaFold structure: {uniprot_id} ({version})")
                        return True

                except urllib.error.HTTPError as e:
                    logger.debug(f"AlphaFold {version} not available: {e.code}")
                    continue

            logger.error(f"AlphaFold structure not found for {uniprot_id} (tried v2, v3, v4)")
            return False

        except Exception as e:
            logger.error(f"Error fetching AlphaFold structure for {uniprot_id}: {e}")
            return False

    def _calculate_protein_center(self, pdb_path: str) -> Optional[Tuple[float, float, float]]:
        """
        Calculate geometric center of protein from PDB file

        Args:
            pdb_path: Path to PDB file

        Returns:
            (x, y, z) coordinates of protein center, or None on failure
        """
        try:
            x_coords, y_coords, z_coords = [], [], []

            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            x_coords.append(x)
                            y_coords.append(y)
                            z_coords.append(z)
                        except (ValueError, IndexError):
                            continue

            if not x_coords:
                logger.error("No atom coordinates found in PDB file")
                return None

            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)
            center_z = sum(z_coords) / len(z_coords)

            logger.info(f"Calculated protein center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
            return (center_x, center_y, center_z)

        except Exception as e:
            logger.error(f"Error calculating protein center: {e}")
            return None

    def _process_protein_data(self, protein_data: Dict, temp_dir: str) -> Optional[str]:
        """
        Process protein_data dictionary to get protein PDB file

        Args:
            protein_data: Dictionary with 'type' and 'value' keys
            temp_dir: Temporary directory for file processing

        Returns:
            protein_pdb_path or None on failure
        """
        try:
            protein_type = protein_data.get('type')
            protein_value = protein_data.get('value')

            if not protein_type or not protein_value:
                logger.error("Invalid protein_data: missing 'type' or 'value'")
                return None

            pdb_path = os.path.join(temp_dir, "protein.pdb")

            # Fetch protein structure based on type
            if protein_type == 'pdb_id':
                if not self._fetch_pdb_structure(protein_value, pdb_path):
                    return None

            elif protein_type == 'uniprot_id':
                if not self._fetch_alphafold_structure(protein_value, pdb_path):
                    return None

            elif protein_type == 'pdb_file':
                # protein_value would contain the PDB file content
                with open(pdb_path, 'w') as f:
                    f.write(protein_value)
                logger.info("Saved uploaded PDB file")

            else:
                logger.error(f"Unknown protein_data type: {protein_type}")
                return None

            return pdb_path

        except Exception as e:
            logger.error(f"Error processing protein_data: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute DiffDock molecular docking

        Args:
            input_data: SMILES string or dict with 'smiles' and 'protein_path'
            **kwargs: Additional parameters:
                - protein_path: Path to protein PDB file (required if not in input_data)
                - complex_name: Name for this docking run (default: "complex")
                - num_poses: Override samples_per_complex

        Returns:
            AdapterResult containing docking results
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string or could not parse with RDKit"
            )

        # Extract SMILES
        if isinstance(input_data, dict):
            smiles = input_data.get('smiles')
        else:
            smiles = input_data

        complex_name = kwargs.get('complex_name', 'complex')

        # Override number of poses if specified
        if 'num_poses' in kwargs:
            original_samples = self.samples_per_complex
            self.samples_per_complex = kwargs['num_poses']

        # Create temporary directory for DiffDock files
        temp_dir = tempfile.mkdtemp(prefix="diffdock_")

        # Handle protein_data if provided
        protein_data = kwargs.get('protein_data')
        protein_path = None

        if protein_data:
            logger.info(f"Processing protein_data: {protein_data.get('type')} = {protein_data.get('value')}")
            protein_path = self._process_protein_data(protein_data, temp_dir)

            if not protein_path:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to process protein_data"
                )

        else:
            # Fallback to protein_path from input_data or kwargs
            if isinstance(input_data, dict):
                protein_path = input_data.get('protein_path')
            if not protein_path:
                protein_path = kwargs.get('protein_path')

            # Validate protein path
            if not protein_path:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No protein_path or protein_data provided. Pass in input_data dict, kwarg, or protein_data."
                )

            if not os.path.exists(protein_path):
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Protein file not found: {protein_path}"
                )

        try:
            # Step 1: Prepare CSV input
            csv_path = os.path.join(temp_dir, "input.csv")
            if not self._prepare_input_csv(protein_path, smiles, csv_path, complex_name):
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to prepare DiffDock input CSV"
                )

            # Step 2: Run DiffDock
            output_dir = os.path.join(temp_dir, "output")
            os.makedirs(output_dir, exist_ok=True)

            logger.info(f"Running DiffDock for: {smiles}")
            success, error_msg = await self._run_diffdock(csv_path, output_dir)

            if not success:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"DiffDock failed: {error_msg}"
                )

            # Step 3: Parse results
            poses = self._parse_diffdock_results(output_dir, complex_name)

            if not poses:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No docking poses found in DiffDock output"
                )

            # Get best pose (highest confidence)
            poses_sorted = sorted(poses, key=lambda x: x.get('confidence', 0), reverse=True)
            best_pose = poses_sorted[0]

            # Estimate binding affinity from confidence
            confidence = best_pose.get('confidence')
            estimated_affinity = self._estimate_binding_affinity(confidence)

            # Normalize confidence to 0-1 score (if available)
            binding_score = None
            if confidence is not None:
                # DiffDock confidence typically ranges 0-5, normalize to 0-1
                binding_score = min(confidence / 5.0, 1.0)

            # Prepare result data
            result_data = {
                "smiles": smiles,
                "protein": os.path.basename(protein_path),
                "confidence": confidence,
                "binding_score": binding_score,  # Normalized 0-1
                "estimated_affinity": estimated_affinity,  # kcal/mol (estimated)
                "best_pose": {
                    "rank": best_pose.get('rank'),
                    "confidence": confidence,
                    "sdf_path": best_pose.get('sdf_path')
                },
                "all_poses": poses,
                "num_poses": len(poses),
                "docking_method": "DiffDock (diffusion-based)"
            }

            logger.info(
                f"DiffDock complete: confidence={confidence:.3f}, "
                f"estimated affinity={estimated_affinity:.2f} kcal/mol"
            )

            # Restore original samples_per_complex if overridden
            if 'num_poses' in kwargs:
                self.samples_per_complex = original_samples

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "inference_steps": self.inference_steps,
                    "samples_per_complex": len(poses)
                }
            )

        except Exception as e:
            logger.error(f"DiffDock adapter error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=str(e)
            )

        finally:
            # Cleanup temporary directory
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                logger.warning(f"Failed to cleanup temp directory: {e}")

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "Diffusion-based molecular docking using DiffDock",
            "features": [
                "Blind docking (no binding site required)",
                "Multiple pose predictions",
                "Confidence scores",
                "State-of-the-art accuracy"
            ],
            "requirements": [
                "DiffDock installation",
                "RDKit (pip install rdkit)",
                "GPU recommended"
            ],
            "config": {
                "diffdock_path": self.diffdock_path,
                "inference_steps": self.inference_steps,
                "samples_per_complex": self.samples_per_complex,
                "use_gpu": self.use_gpu
            },
            "input": {
                "protein": "PDB file path",
                "ligand": "SMILES string"
            },
            "output": {
                "confidence": "DiffDock confidence score (higher = better)",
                "binding_score": "Normalized score 0-1 (higher = better)",
                "estimated_affinity": "Estimated kcal/mol (more negative = better)",
                "poses": "List of docked poses with SDF files"
            }
        }
