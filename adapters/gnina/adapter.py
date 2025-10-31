"""
GNINA Adapter for PharmForge

Performs molecular docking using GNINA (CNN-enhanced AutoDock Vina) to predict
binding affinity and poses with improved scoring using convolutional neural networks.

GNINA combines the speed of AutoDock Vina with CNN-based scoring for more
accurate binding affinity predictions.

References:
    - GNINA Paper: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2
    - GitHub: https://github.com/gnina/gnina
"""

import os
import tempfile
import asyncio
import logging
import shutil
import urllib.request
import subprocess
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
from backend.core.scoring_utils import vina_affinity_to01

logger = logging.getLogger(__name__)


class GNINAAdapter(AdapterProtocol):
    """
    Adapter for GNINA CNN-based molecular docking

    GNINA (formerly smina) is a fork of AutoDock Vina with CNN-based scoring
    for improved binding affinity predictions. It uses deep learning to score
    protein-ligand interactions while maintaining Vina's speed.

    Features:
        - CNN-based scoring for improved accuracy
        - Compatible with AutoDock Vina workflow
        - Both Vina scores and CNN scores
        - Fast docking with better predictions
        - Multiple scoring functions available

    Requirements:
        - GNINA binary (gnina or gnina.exe)
        - RDKit for molecule handling
        - Protein receptor in PDBQT format
        - GPU recommended for CNN scoring

    Input:
        - Protein PDBQT file
        - Ligand SMILES/SDF

    Output:
        - Docked poses
        - Vina scores (kcal/mol)
        - CNN affinity predictions
        - CNNscore and CNNaffinity values
    """

    GNINA_BINARY = "gnina"  # or "gnina.exe" on Windows

    def __init__(
        self,
        name: str = "gnina",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize GNINA adapter

        Args:
            name: Adapter name (default: "gnina")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - receptor_path: Path to receptor PDBQT file (required)
                   - center_x/y/z: Docking box center coordinates (required)
                   - size_x/y/z: Docking box dimensions (default: 25,25,25)
                   - exhaustiveness: Search exhaustiveness (default: 8)
                   - num_modes: Number of binding modes (default: 9)
                   - energy_range: Energy range for output (default: 3 kcal/mol)
                   - gnina_binary: Path to gnina executable (default: "gnina")
                   - cnn_scoring: Which CNN scoring to use (default: "default")
                   - minimize: Whether to minimize poses (default: True)
        """
        super().__init__(name, adapter_type, config or {})
        self.version = "1.0.0"

        # GNINA configuration
        self.receptor_path = self.config.get('receptor_path', None)
        self.center_x = self.config.get('center_x', None)
        self.center_y = self.config.get('center_y', None)
        self.center_z = self.config.get('center_z', None)
        self.size_x = self.config.get('size_x', 25)
        self.size_y = self.config.get('size_y', 25)
        self.size_z = self.config.get('size_z', 25)
        self.exhaustiveness = self.config.get('exhaustiveness', 8)
        self.num_modes = self.config.get('num_modes', 9)
        self.energy_range = self.config.get('energy_range', 3)
        self.gnina_binary = self.config.get('gnina_binary', self.GNINA_BINARY)
        self.cnn_scoring = self.config.get('cnn_scoring', 'default')
        self.minimize = self.config.get('minimize', True)

        # Validate configuration
        if not self.receptor_path:
            logger.warning("No receptor_path configured - adapter will require it at runtime")

        if any(x is None for x in [self.center_x, self.center_y, self.center_z]):
            logger.warning("Docking box center not configured - will use autobox")

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string

        Args:
            input_data: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE:
            logger.error("RDKit not available for validation")
            return False

        if not isinstance(input_data, str):
            return False

        if len(input_data) == 0:
            return False

        # Validate with RDKit
        try:
            mol = Chem.MolFromSmiles(input_data)
            return mol is not None
        except Exception as e:
            logger.debug(f"SMILES validation failed: {e}")
            return False

    def _smiles_to_3d(self, smiles: str) -> Optional[Any]:
        """
        Convert SMILES to 3D molecule with coordinates

        Args:
            smiles: SMILES string

        Returns:
            RDKit molecule with 3D coordinates, or None on failure
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Could not parse SMILES: {smiles}")
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                logger.warning("Standard embedding failed, trying with random coordinates")
                AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

            # Optimize geometry
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)

            return mol

        except Exception as e:
            logger.error(f"Error generating 3D structure: {e}")
            return None

    def _write_sdf(self, mol: Any, filepath: str) -> bool:
        """
        Write molecule to SDF file

        Args:
            mol: RDKit molecule
            filepath: Output SDF file path

        Returns:
            True if successful, False otherwise
        """
        try:
            writer = Chem.SDWriter(filepath)
            writer.write(mol)
            writer.close()
            return True
        except Exception as e:
            logger.error(f"Error writing SDF file: {e}")
            return False

    async def _run_gnina(
        self,
        receptor: str,
        ligand: str,
        output: str,
        log: str
    ) -> Tuple[bool, str]:
        """
        Run GNINA docking

        Args:
            receptor: Path to receptor PDBQT file
            ligand: Path to ligand SDF/PDBQT file
            output: Path for output SDF file
            log: Path for log file

        Returns:
            (success, error_message) tuple
        """
        try:
            # Build GNINA command
            cmd = [
                self.gnina_binary,
                "--receptor", receptor,
                "--ligand", ligand,
                "--out", output,
                "--log", log,
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.num_modes),
            ]

            # Add docking box
            if all(x is not None for x in [self.center_x, self.center_y, self.center_z]):
                cmd.extend([
                    "--center_x", str(self.center_x),
                    "--center_y", str(self.center_y),
                    "--center_z", str(self.center_z),
                    "--size_x", str(self.size_x),
                    "--size_y", str(self.size_y),
                    "--size_z", str(self.size_z),
                ])
            else:
                # Use autobox
                cmd.append("--autobox_ligand")
                cmd.append(ligand)

            # CNN scoring configuration
            if self.cnn_scoring != "none":
                cmd.extend(["--cnn_scoring", self.cnn_scoring])

            # Minimization
            if not self.minimize:
                cmd.append("--no_minimize")

            logger.info(f"Running GNINA: {' '.join(cmd)}")

            # Run GNINA as subprocess
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )

            stdout, stderr = await process.communicate()

            if process.returncode != 0:
                error = stderr.decode() if stderr else "Unknown error"
                logger.error(f"GNINA failed: {error}")
                return False, error

            return True, ""

        except FileNotFoundError:
            error = f"GNINA binary not found: {self.gnina_binary}"
            logger.error(error)
            return False, error
        except Exception as e:
            logger.error(f"Error running GNINA: {e}")
            return False, str(e)

    def _parse_gnina_log(self, log_path: str) -> List[Dict[str, Any]]:
        """
        Parse GNINA log file to extract scores

        GNINA outputs both Vina scores and CNN scores:
        - minimizedAffinity: Vina score after minimization
        - CNNscore: CNN-based pose scoring
        - CNNaffinity: CNN-based affinity prediction

        Args:
            log_path: Path to GNINA log file

        Returns:
            List of docking results with all scores
        """
        results = []

        try:
            with open(log_path, 'r') as f:
                lines = f.readlines()

            # Find results section
            # GNINA log format includes multiple score columns
            parsing = False
            for line in lines:
                # Header line contains score names
                if "mode" in line and "affinity" in line:
                    parsing = True
                    continue

                if parsing and line.strip():
                    # Parse result line
                    parts = line.split()
                    if len(parts) >= 4 and parts[0].isdigit():
                        try:
                            result = {
                                "mode": int(parts[0]),
                                "affinity": float(parts[1]),  # kcal/mol (Vina)
                            }

                            # GNINA may have additional columns for CNN scores
                            # Format varies, typically: mode | affinity | dist | CNNscore | CNNaffinity
                            if len(parts) >= 6:
                                result["dist"] = float(parts[2])
                                result["cnn_score"] = float(parts[3])
                                result["cnn_affinity"] = float(parts[4])

                            results.append(result)
                        except (ValueError, IndexError):
                            continue

        except Exception as e:
            logger.error(f"Error parsing GNINA log: {e}")

        return results

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

    def _pdb_to_pdbqt(self, pdb_path: str, pdbqt_path: str) -> bool:
        """
        Convert PDB receptor to PDBQT format using OpenBabel

        Args:
            pdb_path: Input PDB file path
            pdbqt_path: Output PDBQT file path

        Returns:
            True if successful, False otherwise
        """
        try:
            result = subprocess.run(
                ['obabel', pdb_path, '-O', pdbqt_path, '-xr'],
                capture_output=True,
                text=True,
                timeout=60
            )

            if result.returncode != 0:
                logger.error(f"PDB to PDBQT conversion failed: {result.stderr}")
                return False

            if not os.path.exists(pdbqt_path) or os.path.getsize(pdbqt_path) == 0:
                logger.error("PDBQT file was not created or is empty")
                return False

            logger.info("Converted receptor to PDBQT format")
            return True

        except subprocess.TimeoutExpired:
            logger.error("PDB to PDBQT conversion timed out")
            return False
        except FileNotFoundError:
            logger.error("OpenBabel (obabel) not found")
            return False
        except Exception as e:
            logger.error(f"Error converting PDB to PDBQT: {e}")
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

    def _process_protein_data(self, protein_data: Dict, temp_dir: str) -> Tuple[Optional[str], Optional[Tuple[float, float, float]]]:
        """
        Process protein_data dictionary to get receptor PDBQT and binding site center

        Args:
            protein_data: Dictionary with 'type' and 'value' keys
            temp_dir: Temporary directory for file processing

        Returns:
            (receptor_pdbqt_path, (center_x, center_y, center_z)) or (None, None) on failure
        """
        try:
            protein_type = protein_data.get('type')
            protein_value = protein_data.get('value')

            if not protein_type or not protein_value:
                logger.error("Invalid protein_data: missing 'type' or 'value'")
                return None, None

            pdb_path = os.path.join(temp_dir, "receptor.pdb")

            # Fetch protein structure based on type
            if protein_type == 'pdb_id':
                if not self._fetch_pdb_structure(protein_value, pdb_path):
                    return None, None

            elif protein_type == 'uniprot_id':
                if not self._fetch_alphafold_structure(protein_value, pdb_path):
                    return None, None

            elif protein_type == 'pdb_file':
                # protein_value would contain the PDB file content
                with open(pdb_path, 'w') as f:
                    f.write(protein_value)
                logger.info("Saved uploaded PDB file")

            else:
                logger.error(f"Unknown protein_data type: {protein_type}")
                return None, None

            # Calculate binding site center (geometric center of protein)
            center = self._calculate_protein_center(pdb_path)
            if not center:
                logger.warning("Could not calculate protein center, using default (0,0,0)")
                center = (0.0, 0.0, 0.0)

            # Convert PDB to PDBQT
            pdbqt_path = os.path.join(temp_dir, "receptor.pdbqt")
            if not self._pdb_to_pdbqt(pdb_path, pdbqt_path):
                return None, None

            return pdbqt_path, center

        except Exception as e:
            logger.error(f"Error processing protein_data: {e}")
            return None, None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute GNINA molecular docking

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters:
                - receptor_path: Override default receptor
                - center_x/y/z: Override docking box center
                - size_x/y/z: Override docking box size

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

        smiles = input_data

        # Create temporary directory for docking files
        temp_dir = tempfile.mkdtemp(prefix="gnina_")

        # Handle protein_data if provided
        protein_data = kwargs.get('protein_data')
        receptor_path = None
        center_coords = None

        if protein_data:
            logger.info(f"Processing protein_data: {protein_data.get('type')} = {protein_data.get('value')}")
            receptor_path, center_coords = self._process_protein_data(protein_data, temp_dir)

            if not receptor_path or not center_coords:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to process protein_data"
                )

            # Use calculated center
            self.center_x, self.center_y, self.center_z = center_coords

        else:
            # Fallback to configured receptor_path
            receptor_path = kwargs.get('receptor_path', self.receptor_path)
            if not receptor_path:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No receptor_path or protein_data provided. Configure in adapter or pass as parameter."
                )

            if not os.path.exists(receptor_path):
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Receptor file not found: {receptor_path}"
                )

        try:
            # Step 1: Generate 3D structure
            logger.info(f"Generating 3D structure for: {smiles}")
            mol_3d = self._smiles_to_3d(smiles)
            if mol_3d is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to generate 3D structure from SMILES"
                )

            # Step 2: Write SDF (GNINA can read SDF directly)
            sdf_path = os.path.join(temp_dir, "ligand.sdf")
            if not self._write_sdf(mol_3d, sdf_path):
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to write SDF file"
                )

            # Step 3: Run GNINA docking
            output_path = os.path.join(temp_dir, "output.sdf")
            log_path = os.path.join(temp_dir, "gnina.log")

            logger.info(f"Running GNINA docking for: {smiles}")
            success, error_msg = await self._run_gnina(
                receptor=receptor_path,
                ligand=sdf_path,
                output=output_path,
                log=log_path
            )

            if not success:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"GNINA docking failed: {error_msg}"
                )

            # Step 4: Parse results
            poses = self._parse_gnina_log(log_path)

            if not poses:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No docking poses found in GNINA output"
                )

            # Get best pose
            best_pose = poses[0]
            vina_affinity = best_pose["affinity"]  # kcal/mol

            # Use CNN affinity if available, otherwise use Vina
            cnn_affinity = best_pose.get("cnn_affinity")
            cnn_score = best_pose.get("cnn_score")

            # Primary binding affinity (prefer CNN if available)
            binding_affinity = cnn_affinity if cnn_affinity is not None else vina_affinity

            # Normalize score (0-1, higher=better)
            binding_score = vina_affinity_to01(binding_affinity)

            # Prepare result data
            result_data = {
                "smiles": smiles,
                "vina_affinity": vina_affinity,  # Vina score (kcal/mol)
                "cnn_affinity": cnn_affinity,  # CNN affinity (kcal/mol)
                "cnn_score": cnn_score,  # CNN pose score
                "binding_affinity": binding_affinity,  # Primary score
                "binding_score": binding_score,  # Normalized 0-1
                "best_pose": best_pose,
                "all_poses": poses,
                "num_poses": len(poses),
                "receptor": os.path.basename(receptor_path),
                "docking_box": {
                    "center": [self.center_x, self.center_y, self.center_z],
                    "size": [self.size_x, self.size_y, self.size_z]
                },
                "scoring_method": "GNINA (CNN + Vina)"
            }

            logger.info(
                f"GNINA docking complete: Vina={vina_affinity:.2f}, "
                f"CNN={cnn_affinity:.2f if cnn_affinity else 'N/A'} kcal/mol "
                f"(score: {binding_score:.3f})"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "exhaustiveness": self.exhaustiveness,
                    "num_modes": self.num_modes,
                    "cnn_scoring": self.cnn_scoring
                }
            )

        except Exception as e:
            logger.error(f"GNINA adapter error: {e}", exc_info=True)
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
            "description": "CNN-enhanced molecular docking using GNINA",
            "features": [
                "CNN-based scoring",
                "Vina-compatible docking",
                "Multiple scoring functions",
                "Fast and accurate"
            ],
            "requirements": [
                "GNINA binary",
                "RDKit (pip install rdkit)",
                "Receptor PDBQT file",
                "GPU recommended for CNN scoring"
            ],
            "config": {
                "receptor_path": self.receptor_path,
                "docking_box_center": [self.center_x, self.center_y, self.center_z],
                "docking_box_size": [self.size_x, self.size_y, self.size_z],
                "exhaustiveness": self.exhaustiveness,
                "num_modes": self.num_modes,
                "gnina_binary": self.gnina_binary,
                "cnn_scoring": self.cnn_scoring
            },
            "scoring": {
                "vina_score": "kcal/mol (Vina scoring function)",
                "cnn_affinity": "kcal/mol (CNN prediction)",
                "cnn_score": "CNN pose score",
                "binding_score": "Normalized 0-1 (higher is better)",
                "normalization_function": "vina_affinity_to01()"
            }
        }
