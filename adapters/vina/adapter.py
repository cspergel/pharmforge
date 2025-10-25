"""
Vina Molecular Docking Adapter for PharmForge

Performs molecular docking using AutoDock Vina to predict binding affinity
and poses of ligands to protein targets.
"""

import os
import tempfile
import asyncio
import logging
from pathlib import Path
from typing import Any, Dict, Optional, List, Tuple, TYPE_CHECKING
import shutil

if TYPE_CHECKING:
    from rdkit import Chem

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    Chem = None  # type: ignore
    logging.warning("RDKit not available - install with: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult
from backend.core.scoring_utils import vina_affinity_to01

logger = logging.getLogger(__name__)


class VinaAdapter(AdapterProtocol):
    """
    Adapter for AutoDock Vina molecular docking

    Performs protein-ligand docking to predict binding affinity and poses.

    Requirements:
        - AutoDock Vina binary (vina or vina.exe)
        - RDKit for molecule handling
        - Open Babel (optional, for format conversions)

    Workflow:
        1. Validate SMILES input
        2. Generate 3D conformer for ligand
        3. Prepare ligand in PDBQT format
        4. Run Vina docking against receptor
        5. Parse results and normalize scores

    Returns:
        - Binding affinity (kcal/mol, raw score)
        - Normalized score (0-1, higher=better)
        - Top docking poses with coordinates
    """

    DEFAULT_RECEPTOR = "default_receptor.pdbqt"  # Placeholder
    VINA_BINARY = "vina"  # or "vina.exe" on Windows

    def __init__(self, name: str = "vina_docking", adapter_type: str = "local", config: Optional[Dict[str, Any]] = None):
        """
        Initialize Vina adapter

        Args:
            name: Adapter name (default: "vina_docking")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - receptor_path: Path to receptor PDBQT file (required)
                   - center_x/y/z: Docking box center coordinates (required)
                   - size_x/y/z: Docking box dimensions (default: 25,25,25)
                   - exhaustiveness: Search exhaustiveness (default: 8)
                   - num_modes: Number of binding modes (default: 9)
                   - energy_range: Energy range for output (default: 3 kcal/mol)
                   - vina_binary: Path to vina executable (default: "vina")
        """
        super().__init__(name, adapter_type, config)
        self.version = "1.0.0"

        # Vina configuration
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
        self.vina_binary = self.config.get('vina_binary', self.VINA_BINARY)

        # Validate configuration
        self._validate_configuration()

        # Check dependencies
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

    def _validate_configuration(self):
        """Validate that required configuration is present"""
        if not self.receptor_path:
            logger.warning("No receptor_path configured - adapter will require it at runtime")

        if any(x is None for x in [self.center_x, self.center_y, self.center_z]):
            logger.warning("Docking box center not configured - will use receptor center")

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

        # Try to parse with RDKit
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

            # Optimize geometry with UFF
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)

            return mol

        except Exception as e:
            logger.error(f"Error generating 3D structure: {e}")
            return None

    def _write_pdb(self, mol: Any, filepath: str) -> bool:
        """
        Write molecule to PDB file

        Args:
            mol: RDKit molecule
            filepath: Output PDB file path

        Returns:
            True if successful, False otherwise
        """
        try:
            writer = Chem.PDBWriter(filepath)
            writer.write(mol)
            writer.close()
            return True
        except Exception as e:
            logger.error(f"Error writing PDB file: {e}")
            return False

    def _prepare_ligand_pdbqt(self, pdb_path: str, pdbqt_path: str) -> bool:
        """
        Convert PDB to PDBQT format (required by Vina)

        Uses simple charge assignment since we don't have Open Babel/MGLTools.
        In production, you'd use: obabel or prepare_ligand4.py

        Args:
            pdb_path: Input PDB file path
            pdbqt_path: Output PDBQT file path

        Returns:
            True if successful, False otherwise
        """
        try:
            # Simple PDBQT generation (assigns Gasteiger charges via RDKit)
            mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
            if mol is None:
                logger.error("Could not read PDB file")
                return False

            # Compute Gasteiger charges
            AllChem.ComputeGasteigerCharges(mol)

            # Write PDBQT (simplified format)
            with open(pdbqt_path, 'w') as f:
                conf = mol.GetConformer()
                for atom in mol.GetAtoms():
                    idx = atom.GetIdx()
                    pos = conf.GetAtomPosition(idx)
                    charge = atom.GetDoubleProp('_GasteigerCharge')

                    # Simple PDBQT ATOM line
                    atom_type = atom.GetSymbol()
                    f.write(f"ATOM  {idx+1:5d}  {atom_type:<3s} LIG     1    "
                           f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
                           f"  1.00  0.00    {charge:6.3f} {atom_type}\n")

                f.write("ENDMDL\n")

            return True

        except Exception as e:
            logger.error(f"Error preparing PDBQT: {e}")
            return False

    async def _run_vina(
        self,
        receptor: str,
        ligand: str,
        output: str,
        log: str
    ) -> Tuple[bool, str]:
        """
        Run AutoDock Vina docking

        Args:
            receptor: Path to receptor PDBQT file
            ligand: Path to ligand PDBQT file
            output: Path for output PDBQT file
            log: Path for log file

        Returns:
            (success, error_message) tuple
        """
        try:
            # Build Vina command
            cmd = [
                self.vina_binary,
                "--receptor", receptor,
                "--ligand", ligand,
                "--out", output,
                "--log", log,
                "--center_x", str(self.center_x),
                "--center_y", str(self.center_y),
                "--center_z", str(self.center_z),
                "--size_x", str(self.size_x),
                "--size_y", str(self.size_y),
                "--size_z", str(self.size_z),
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.num_modes),
                "--energy_range", str(self.energy_range)
            ]

            # Run Vina as subprocess
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )

            stdout, stderr = await process.communicate()

            if process.returncode != 0:
                error = stderr.decode() if stderr else "Unknown error"
                logger.error(f"Vina failed: {error}")
                return False, error

            return True, ""

        except FileNotFoundError:
            error = f"Vina binary not found: {self.vina_binary}"
            logger.error(error)
            return False, error

        except Exception as e:
            logger.error(f"Error running Vina: {e}")
            return False, str(e)

    def _parse_vina_log(self, log_path: str) -> List[Dict[str, Any]]:
        """
        Parse Vina log file to extract binding affinities

        Args:
            log_path: Path to Vina log file

        Returns:
            List of docking results (mode, affinity, rmsd_lb, rmsd_ub)
        """
        results = []

        try:
            with open(log_path, 'r') as f:
                lines = f.readlines()

            # Find results section
            parsing = False
            for line in lines:
                if "-----+------------+----------+----------" in line:
                    parsing = True
                    continue

                if parsing and line.strip():
                    # Parse result line
                    # Format: "   1        -7.5      0.000      0.000"
                    parts = line.split()
                    if len(parts) >= 4 and parts[0].isdigit():
                        try:
                            result = {
                                "mode": int(parts[0]),
                                "affinity": float(parts[1]),  # kcal/mol
                                "rmsd_lb": float(parts[2]),
                                "rmsd_ub": float(parts[3])
                            }
                            results.append(result)
                        except ValueError:
                            continue

        except Exception as e:
            logger.error(f"Error parsing Vina log: {e}")

        return results

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Vina molecular docking

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

        # Override config with kwargs
        receptor_path = kwargs.get('receptor_path', self.receptor_path)
        if not receptor_path:
            return AdapterResult(
                success=False,
                data=None,
                error="No receptor_path provided. Configure in adapter or pass as parameter."
            )

        if not os.path.exists(receptor_path):
            return AdapterResult(
                success=False,
                data=None,
                error=f"Receptor file not found: {receptor_path}"
            )

        # Create temporary directory for docking files
        temp_dir = tempfile.mkdtemp(prefix="vina_")

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

            # Step 2: Write PDB
            pdb_path = os.path.join(temp_dir, "ligand.pdb")
            if not self._write_pdb(mol_3d, pdb_path):
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to write PDB file"
                )

            # Step 3: Prepare PDBQT
            pdbqt_path = os.path.join(temp_dir, "ligand.pdbqt")
            if not self._prepare_ligand_pdbqt(pdb_path, pdbqt_path):
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to prepare PDBQT file"
                )

            # Step 4: Run Vina docking
            output_path = os.path.join(temp_dir, "output.pdbqt")
            log_path = os.path.join(temp_dir, "vina.log")

            logger.info(f"Running Vina docking for: {smiles}")
            success, error_msg = await self._run_vina(
                receptor=receptor_path,
                ligand=pdbqt_path,
                output=output_path,
                log=log_path
            )

            if not success:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Vina docking failed: {error_msg}"
                )

            # Step 5: Parse results
            poses = self._parse_vina_log(log_path)

            if not poses:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No docking poses found in Vina output"
                )

            # Get best pose
            best_pose = poses[0]
            binding_affinity = best_pose["affinity"]  # kcal/mol (more negative = better)

            # Normalize score using scoring_utils (0-1, higher=better)
            binding_score = vina_affinity_to01(binding_affinity)

            # Prepare result data
            result_data = {
                "smiles": smiles,
                "binding_affinity": binding_affinity,  # Raw kcal/mol
                "binding_score": binding_score,  # Normalized 0-1
                "best_pose": best_pose,
                "all_poses": poses,
                "num_poses": len(poses),
                "receptor": os.path.basename(receptor_path),
                "docking_box": {
                    "center": [self.center_x, self.center_y, self.center_z],
                    "size": [self.size_x, self.size_y, self.size_z]
                }
            }

            logger.info(
                f"✓ Docking complete: {binding_affinity:.2f} kcal/mol "
                f"(score: {binding_score:.3f})"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "vina_exhaustiveness": self.exhaustiveness,
                    "num_modes": self.num_modes
                }
            )

        except Exception as e:
            logger.error(f"Vina adapter error: {e}", exc_info=True)
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
            "description": "Molecular docking using AutoDock Vina",
            "requirements": [
                "AutoDock Vina binary",
                "RDKit (pip install rdkit)",
                "Receptor PDBQT file"
            ],
            "config": {
                "receptor_path": self.receptor_path,
                "docking_box_center": [self.center_x, self.center_y, self.center_z],
                "docking_box_size": [self.size_x, self.size_y, self.size_z],
                "exhaustiveness": self.exhaustiveness,
                "num_modes": self.num_modes,
                "vina_binary": self.vina_binary
            },
            "scoring": {
                "raw_unit": "kcal/mol",
                "raw_range": "typically -12 to -4",
                "normalized_range": "0 to 1 (higher is better)",
                "normalization_function": "vina_affinity_to01()"
            }
        }
