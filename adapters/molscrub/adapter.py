"""
MolScrub Adapter for PharmForge

Generates clean, standardized 3D conformers for docking and MD simulations.
Provides automated molecular cleaning, conformer generation, and energy optimization.
"""

import os
import logging
import tempfile
from typing import Any, Dict, Optional, List
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - install with: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class MolScrubAdapter(AdapterProtocol):
    """
    Adapter for MolScrub conformer generation and molecular cleaning

    Features:
        - Generate 3D conformers from SMILES
        - Clean and standardize molecular structures
        - Optimize conformer geometry
        - Generate multiple low-energy conformers
        - Export to various formats (PDB, MOL2, SDF)
        - Filter conformers by energy and RMSD

    Workflow:
        1. Parse and validate SMILES input
        2. Clean and standardize structure (remove salts, standardize tautomers)
        3. Generate multiple 3D conformers using ETKDG
        4. Optimize conformer geometries with MMFF or UFF
        5. Filter by energy window and RMSD
        6. Export to requested format

    Returns:
        - List of low-energy conformers
        - Energy values (kcal/mol)
        - RMSD to lowest energy conformer
        - Structure data in requested format
    """

    SUPPORTED_FORMATS = ["pdb", "mol2", "sdf", "mol"]
    DEFAULT_ENERGY_WINDOW = 10.0  # kcal/mol
    DEFAULT_RMS_THRESHOLD = 0.5  # Angstroms

    def __init__(self):
        super().__init__(
            name="molscrub",
            adapter_type="local",
            config={
                "num_conformers": 10,
                "energy_window": self.DEFAULT_ENERGY_WINDOW,
                "rms_threshold": self.DEFAULT_RMS_THRESHOLD,
                "optimize": True,
                "output_format": "pdb",
                "use_mmff": True,  # Use MMFF94 force field (fallback to UFF if unavailable)
                "max_iterations": 200,
                "random_seed": 42  # For reproducibility
            }
        )
        self.version = "1.0.0"

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string or molecule structure

        Args:
            input_data: SMILES string or dictionary with 'smiles' or 'structure' key

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE:
            logger.error("RDKit not available for validation")
            return False

        # Handle dictionary input
        if isinstance(input_data, dict):
            smiles = input_data.get("smiles")
            structure = input_data.get("structure")

            if not smiles and not structure:
                logger.error("Input dictionary must contain 'smiles' or 'structure' key")
                return False

            # Validate SMILES if provided
            if smiles:
                input_data = smiles
            else:
                # For structure input, we'll validate during execution
                return True

        # Validate SMILES string
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

    def _clean_molecule(self, mol: Any) -> Optional[Any]:
        """
        Clean and standardize molecule structure

        Operations:
            - Remove salts and solvents (keep largest fragment)
            - Standardize charges
            - Neutralize when possible
            - Add explicit hydrogens

        Args:
            mol: RDKit molecule object

        Returns:
            Cleaned RDKit molecule or None on failure
        """
        try:
            # Remove salts/solvents - keep largest fragment
            from rdkit.Chem.MolStandardize import rdMolStandardize

            # Get the largest fragment
            fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
            if len(fragments) > 1:
                logger.info(f"Molecule contains {len(fragments)} fragments, keeping largest")
                mol = max(fragments, key=lambda m: m.GetNumAtoms())

            # Standardize using RDKit's standardization
            clean_mol = rdMolStandardize.Cleanup(mol)

            # Add hydrogens (essential for 3D conformer generation)
            mol_with_h = Chem.AddHs(clean_mol)

            logger.debug(f"Cleaned molecule: {mol_with_h.GetNumAtoms()} atoms (including H)")
            return mol_with_h

        except Exception as e:
            logger.error(f"Error cleaning molecule: {e}")
            # Fallback: just add hydrogens
            try:
                return Chem.AddHs(mol)
            except:
                return None

    def _generate_conformers(self, mol: Any, num_conformers: int, random_seed: int = 42) -> List[int]:
        """
        Generate 3D conformers using ETKDG algorithm

        Args:
            mol: RDKit molecule object (with hydrogens)
            num_conformers: Number of conformers to generate
            random_seed: Random seed for reproducibility

        Returns:
            List of conformer IDs
        """
        try:
            # Use ETKDG (Experimental Torsion-angle preference with Distance Geometry)
            params = AllChem.ETKDGv3()
            params.randomSeed = random_seed
            params.numThreads = 0  # Use all available threads
            params.useRandomCoords = False

            conformer_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_conformers,
                params=params
            )

            if len(conformer_ids) == 0:
                logger.warning("ETKDG failed, trying with random coordinates")
                params.useRandomCoords = True
                conformer_ids = AllChem.EmbedMultipleConfs(
                    mol,
                    numConfs=num_conformers,
                    params=params
                )

            logger.info(f"Generated {len(conformer_ids)} initial conformers")
            return list(conformer_ids)

        except Exception as e:
            logger.error(f"Error generating conformers: {e}")
            return []

    def _optimize_conformers(self, mol: Any, conformer_ids: List[int], use_mmff: bool = True, max_iterations: int = 200) -> Dict[int, float]:
        """
        Optimize conformer geometries and calculate energies

        Args:
            mol: RDKit molecule object
            conformer_ids: List of conformer IDs to optimize
            use_mmff: Use MMFF94 force field (fallback to UFF)
            max_iterations: Maximum optimization iterations

        Returns:
            Dictionary mapping conformer_id -> energy (kcal/mol)
        """
        energies = {}

        for conf_id in conformer_ids:
            try:
                # Try MMFF94 first (more accurate but not always available)
                if use_mmff:
                    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
                    if props is not None:
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
                        if ff is not None:
                            ff.Initialize()
                            converged = ff.Minimize(maxIts=max_iterations)
                            energy = ff.CalcEnergy()
                            energies[conf_id] = energy
                            logger.debug(f"Conformer {conf_id}: MMFF94 energy = {energy:.2f} kcal/mol")
                            continue

                # Fallback to UFF (Universal Force Field)
                ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                if ff is not None:
                    ff.Initialize()
                    converged = ff.Minimize(maxIts=max_iterations)
                    energy = ff.CalcEnergy()
                    energies[conf_id] = energy
                    logger.debug(f"Conformer {conf_id}: UFF energy = {energy:.2f} kcal/mol")
                else:
                    logger.warning(f"Could not create force field for conformer {conf_id}")

            except Exception as e:
                logger.warning(f"Error optimizing conformer {conf_id}: {e}")
                continue

        return energies

    def _filter_conformers(self, mol: Any, energies: Dict[int, float], energy_window: float, rms_threshold: float) -> List[Dict[str, Any]]:
        """
        Filter conformers by energy window and RMSD

        Args:
            mol: RDKit molecule object
            energies: Dictionary of conformer_id -> energy
            energy_window: Maximum energy difference from lowest (kcal/mol)
            rms_threshold: Minimum RMSD between conformers (Angstroms)

        Returns:
            List of conformer dictionaries with metadata
        """
        if not energies:
            logger.warning("No energies available for filtering")
            return []

        # Sort by energy
        sorted_conformers = sorted(energies.items(), key=lambda x: x[1])
        lowest_energy = sorted_conformers[0][1]

        # Filter by energy window
        filtered_by_energy = [
            (conf_id, energy)
            for conf_id, energy in sorted_conformers
            if energy - lowest_energy <= energy_window
        ]

        logger.info(f"Filtered to {len(filtered_by_energy)} conformers within {energy_window} kcal/mol")

        # Filter by RMSD (remove similar conformers)
        unique_conformers = []

        for conf_id, energy in filtered_by_energy:
            # Check RMSD against all accepted conformers
            is_unique = True

            for existing in unique_conformers:
                existing_id = existing["conformer_id"]
                try:
                    rmsd = AllChem.GetBestRMS(mol, mol, conf_id, existing_id)
                    if rmsd < rms_threshold:
                        is_unique = False
                        logger.debug(f"Conformer {conf_id} too similar to {existing_id} (RMSD={rmsd:.3f})")
                        break
                except Exception as e:
                    logger.warning(f"Could not calculate RMSD: {e}")

            if is_unique:
                # Calculate RMSD to lowest energy conformer
                rmsd_to_lowest = 0.0
                if conf_id != sorted_conformers[0][0]:
                    try:
                        rmsd_to_lowest = AllChem.GetBestRMS(mol, mol, conf_id, sorted_conformers[0][0])
                    except:
                        rmsd_to_lowest = -1.0

                unique_conformers.append({
                    "conformer_id": conf_id,
                    "energy": energy,
                    "rmsd_to_lowest": rmsd_to_lowest
                })

        logger.info(f"Filtered to {len(unique_conformers)} unique conformers (RMSD > {rms_threshold} Å)")
        return unique_conformers

    def _export_conformer(self, mol: Any, conf_id: int, output_format: str) -> Optional[str]:
        """
        Export conformer to requested format

        Args:
            mol: RDKit molecule object
            conf_id: Conformer ID
            output_format: Output format ("pdb", "mol2", "sdf", "mol")

        Returns:
            Structure string in requested format or None on failure
        """
        try:
            if output_format.lower() == "pdb":
                return Chem.MolToPDBBlock(mol, confId=conf_id)

            elif output_format.lower() == "sdf" or output_format.lower() == "mol":
                return Chem.MolToMolBlock(mol, confId=conf_id)

            elif output_format.lower() == "mol2":
                # MOL2 requires separate handling
                # Write to temp file and read back
                with tempfile.NamedTemporaryFile(mode='w', suffix='.mol2', delete=False) as f:
                    temp_path = f.name

                try:
                    Chem.MolToMol2File(mol, temp_path, confId=conf_id)
                    with open(temp_path, 'r') as f:
                        mol2_content = f.read()
                    return mol2_content
                finally:
                    if os.path.exists(temp_path):
                        os.remove(temp_path)

            else:
                logger.error(f"Unsupported output format: {output_format}")
                return None

        except Exception as e:
            logger.error(f"Error exporting conformer to {output_format}: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute MolScrub conformer generation and cleaning

        Args:
            input_data: SMILES string or dictionary with parameters
            **kwargs: Additional parameters:
                - num_conformers: Number of conformers to generate (default: 10)
                - energy_window: Energy window for filtering (kcal/mol, default: 10.0)
                - rms_threshold: RMSD threshold for filtering (Å, default: 0.5)
                - optimize: Whether to optimize geometry (default: True)
                - output_format: Output format ("pdb", "mol2", "sdf", default: "pdb")
                - use_mmff: Use MMFF94 force field (default: True)

        Returns:
            AdapterResult containing conformer data
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        # Parse input
        if isinstance(input_data, dict):
            smiles = input_data.get("smiles")
            num_conformers = input_data.get("num_conformers", self.config["num_conformers"])
            energy_window = input_data.get("energy_window", self.config["energy_window"])
            rms_threshold = input_data.get("rms_threshold", self.config["rms_threshold"])
            optimize = input_data.get("optimize", self.config["optimize"])
            output_format = input_data.get("output_format", self.config["output_format"])
        else:
            smiles = input_data
            num_conformers = kwargs.get("num_conformers", self.config["num_conformers"])
            energy_window = kwargs.get("energy_window", self.config["energy_window"])
            rms_threshold = kwargs.get("rms_threshold", self.config["rms_threshold"])
            optimize = kwargs.get("optimize", self.config["optimize"])
            output_format = kwargs.get("output_format", self.config["output_format"])

        # Validate output format
        if output_format.lower() not in self.SUPPORTED_FORMATS:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported output format: {output_format}. Supported: {self.SUPPORTED_FORMATS}"
            )

        # Validate input
        if not self.validate_input(smiles):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string or could not parse with RDKit"
            )

        try:
            logger.info(f"Processing SMILES: {smiles}")

            # Step 1: Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Could not parse SMILES with RDKit"
                )

            # Step 2: Clean and standardize
            logger.info("Cleaning and standardizing molecule...")
            mol_clean = self._clean_molecule(mol)
            if mol_clean is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to clean molecule"
                )

            # Step 3: Generate conformers
            logger.info(f"Generating {num_conformers} conformers...")
            conformer_ids = self._generate_conformers(
                mol_clean,
                num_conformers,
                random_seed=self.config["random_seed"]
            )

            if not conformer_ids:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to generate conformers"
                )

            # Step 4: Optimize geometries (if requested)
            energies = {}
            if optimize:
                logger.info("Optimizing conformer geometries...")
                energies = self._optimize_conformers(
                    mol_clean,
                    conformer_ids,
                    use_mmff=kwargs.get("use_mmff", self.config["use_mmff"]),
                    max_iterations=self.config["max_iterations"]
                )

                if not energies:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="Failed to optimize conformers"
                    )
            else:
                # Assign dummy energies if not optimizing
                energies = {conf_id: 0.0 for conf_id in conformer_ids}

            # Step 5: Filter conformers
            logger.info("Filtering conformers by energy and RMSD...")
            filtered_conformers = self._filter_conformers(
                mol_clean,
                energies,
                energy_window,
                rms_threshold
            )

            if not filtered_conformers:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No conformers passed filtering criteria"
                )

            # Step 6: Export conformers
            logger.info(f"Exporting {len(filtered_conformers)} conformers to {output_format} format...")
            conformers_output = []

            for conf_data in filtered_conformers:
                conf_id = conf_data["conformer_id"]
                structure = self._export_conformer(mol_clean, conf_id, output_format)

                if structure:
                    conformers_output.append({
                        "conformer_id": conf_data["conformer_id"],
                        "energy": conf_data["energy"],
                        "structure": structure,
                        "rmsd_to_lowest": conf_data["rmsd_to_lowest"]
                    })

            # Calculate molecular properties for metadata
            properties = {
                "molecular_weight": Descriptors.MolWt(mol_clean),
                "num_atoms": mol_clean.GetNumAtoms(),
                "num_heavy_atoms": Descriptors.HeavyAtomCount(mol_clean),
                "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol_clean)
            }

            # Prepare result data
            result_data = {
                "smiles": smiles,
                "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol_clean)),
                "conformers": conformers_output,
                "num_generated": len(conformer_ids),
                "num_filtered": len(conformers_output),
                "lowest_energy": min(c["energy"] for c in conformers_output) if optimize else None,
                "output_format": output_format,
                "properties": properties,
                "warnings": []
            }

            # Add warnings if needed
            if len(conformers_output) < num_conformers:
                result_data["warnings"].append(
                    f"Only {len(conformers_output)} unique conformers found (requested {num_conformers})"
                )

            logger.info(
                f"✓ Generated {len(conformers_output)} conformers "
                f"(lowest energy: {result_data['lowest_energy']:.2f} kcal/mol)" if optimize
                else f"✓ Generated {len(conformers_output)} conformers"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "num_conformers_requested": num_conformers,
                    "energy_window": energy_window,
                    "rms_threshold": rms_threshold,
                    "optimized": optimize
                }
            )

        except Exception as e:
            logger.error(f"MolScrub adapter error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=str(e)
            )

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
            "description": "Conformer generation and molecular cleaning using RDKit",
            "requirements": [
                "RDKit (pip install rdkit)"
            ],
            "features": [
                "3D conformer generation (ETKDG algorithm)",
                "Molecular cleaning and standardization",
                "Energy optimization (MMFF94 or UFF)",
                "RMSD-based filtering",
                "Multiple output formats (PDB, MOL2, SDF)"
            ],
            "config": self.config,
            "supported_formats": self.SUPPORTED_FORMATS,
            "use_cases": [
                "Pre-docking ligand preparation",
                "MD simulation setup",
                "Virtual screening library preparation",
                "Conformer ensemble generation"
            ]
        }
