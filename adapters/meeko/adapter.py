"""
Meeko Ligand Preparation Adapter for PharmForge

Prepares ligands and receptors for AutoDock Vina/GPU using Meeko.
Handles PDBQT generation, macrocycles, flexible residues, and conformer generation.
"""

import logging
from typing import Any, Dict, Optional, List, TYPE_CHECKING
import tempfile
import os

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

try:
    from meeko import MoleculePreparation, PDBQTMolecule, RDKitMolCreate
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False
    logging.warning("Meeko not available - install with: pip install meeko")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class MeekoAdapter(AdapterProtocol):
    """
    Adapter for Meeko ligand preparation

    Prepares small molecules for docking with AutoDock Vina/GPU by generating
    PDBQT format files with proper atom types, charges, and torsion trees.

    Features:
        - SMILES to PDBQT conversion
        - Macrocycle support with ring flexibility
        - Multiple conformer generation
        - Flexible residue preparation
        - Hydrated docking support
        - Proper handling of atom types and charges

    Requirements:
        - Meeko (pip install meeko)
        - RDKit for molecule handling

    Workflow:
        1. Parse SMILES with RDKit
        2. Generate 3D conformers
        3. Prepare molecule with Meeko (detect macrocycles, assign charges)
        4. Generate PDBQT format
        5. Return PDBQT string and metadata

    Returns:
        - PDBQT string (ready for Vina)
        - Number of rotatable bonds
        - Macrocycle information
        - Conformer details
    """

    def __init__(
        self,
        name: str = "meeko",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize Meeko adapter

        Args:
            name: Adapter name (default: "meeko")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - rigid_macrocycles: Keep macrocycles rigid (default: False)
                   - keep_nonpolar_hydrogens: Keep nonpolar H atoms (default: False)
                   - merge_these_atom_types: Merge atom types list (default: None)
                   - hydrate: Add hydration sites (default: False)
                   - flexible_amides: Allow amide flexibility (default: False)
                   - rigidify_bonds_smarts: SMARTS patterns for rigid bonds (default: [])
                   - rigidify_bonds_indices: Bond indices to rigidify (default: [])
                   - double_bond_penalty: Penalty for double bond torsions (default: 50.0)
                   - min_ring_size: Minimum ring size for flexibility (default: 7)
                   - num_conformers: Number of conformers to generate (default: 1)
        """
        super().__init__(name, adapter_type, config)
        self.version = "1.0.0"

        # Meeko configuration parameters
        self.rigid_macrocycles = self.config.get('rigid_macrocycles', False)
        self.keep_nonpolar_hydrogens = self.config.get('keep_nonpolar_hydrogens', False)
        self.merge_these_atom_types = self.config.get('merge_these_atom_types', None)
        self.hydrate = self.config.get('hydrate', False)
        self.flexible_amides = self.config.get('flexible_amides', False)
        self.rigidify_bonds_smarts = self.config.get('rigidify_bonds_smarts', [])
        self.rigidify_bonds_indices = self.config.get('rigidify_bonds_indices', [])
        self.double_bond_penalty = self.config.get('double_bond_penalty', 50.0)
        self.min_ring_size = self.config.get('min_ring_size', 7)
        self.num_conformers = self.config.get('num_conformers', 1)

        # Check dependencies
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")
        if not MEEKO_AVAILABLE:
            logger.error("Meeko is not installed! Install with: pip install meeko")

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

            # Optimize geometry with UFF or MMFF
            try:
                # Try MMFF first (usually better)
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except:
                # Fall back to UFF
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)

            return mol

        except Exception as e:
            logger.error(f"Error generating 3D structure: {e}")
            return None

    def _generate_conformers(self, mol: Any, num_conformers: int) -> List[Any]:
        """
        Generate multiple conformers for a molecule

        Args:
            mol: RDKit molecule
            num_conformers: Number of conformers to generate

        Returns:
            List of RDKit molecules with different conformers
        """
        if num_conformers <= 1:
            return [mol]

        try:
            # Generate conformers
            conf_ids = AllChem.EmbedMultipleConfs(
                mol,
                numConfs=num_conformers,
                randomSeed=42,
                pruneRmsThresh=0.5  # Remove similar conformers
            )

            # Optimize each conformer
            results = []
            for conf_id in conf_ids:
                try:
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, maxIters=200)
                except:
                    AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=200)

            # Create separate molecule objects for each conformer
            conformers = []
            for conf_id in conf_ids:
                mol_copy = Chem.Mol(mol)
                mol_copy.RemoveAllConformers()
                mol_copy.AddConformer(mol.GetConformer(conf_id), assignId=True)
                conformers.append(mol_copy)

            logger.info(f"Generated {len(conformers)} conformers")
            return conformers

        except Exception as e:
            logger.warning(f"Conformer generation failed: {e}, using single conformer")
            return [mol]

    def _prepare_with_meeko(self, mol: Any) -> Optional[Dict[str, Any]]:
        """
        Prepare molecule with Meeko for docking

        Args:
            mol: RDKit molecule with 3D coordinates

        Returns:
            Dictionary with PDBQT string and metadata, or None on failure
        """
        try:
            # Create Meeko molecule preparation object
            preparator = MoleculePreparation(
                rigid_macrocycles=self.rigid_macrocycles,
                keep_nonpolar_hydrogens=self.keep_nonpolar_hydrogens,
                merge_these_atom_types=self.merge_these_atom_types,
                hydrate=self.hydrate,
                flexible_amides=self.flexible_amides,
                rigidify_bonds_smarts=self.rigidify_bonds_smarts,
                rigidify_bonds_indices=self.rigidify_bonds_indices,
                double_bond_penalty=self.double_bond_penalty,
                min_ring_size=self.min_ring_size
            )

            # Prepare the molecule
            preparator.prepare(mol)

            # Get the prepared molecule as PDBQT
            pdbqt_string = preparator.write_pdbqt_string()

            # Extract metadata
            metadata = {
                "num_torsions": len(preparator.setup.rotatable_bonds) if hasattr(preparator, 'setup') else 0,
                "macrocycles_detected": preparator.setup.macrocycle_list if hasattr(preparator, 'setup') else [],
                "num_macrocycles": len(preparator.setup.macrocycle_list) if hasattr(preparator, 'setup') and hasattr(preparator.setup, 'macrocycle_list') else 0,
                "num_atoms": mol.GetNumAtoms(),
                "num_heavy_atoms": mol.GetNumHeavyAtoms(),
                "rigid_macrocycles": self.rigid_macrocycles,
                "flexible_amides": self.flexible_amides,
                "hydrated": self.hydrate
            }

            return {
                "pdbqt_string": pdbqt_string,
                "metadata": metadata
            }

        except Exception as e:
            logger.error(f"Error preparing molecule with Meeko: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Meeko ligand preparation

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters:
                - num_conformers: Number of conformers to generate (overrides config)
                - rigid_macrocycles: Keep macrocycles rigid (overrides config)
                - hydrate: Add hydration sites (overrides config)
                - flexible_amides: Allow amide flexibility (overrides config)
                - output_format: 'string' (default) or 'file' (returns temp file path)

        Returns:
            AdapterResult containing PDBQT data and preparation metadata
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        if not MEEKO_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="Meeko is not installed. Install with: pip install meeko"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string or could not parse with RDKit"
            )

        smiles = input_data

        # Override config with kwargs if provided
        num_conformers = kwargs.get('num_conformers', self.num_conformers)
        self.rigid_macrocycles = kwargs.get('rigid_macrocycles', self.rigid_macrocycles)
        self.hydrate = kwargs.get('hydrate', self.hydrate)
        self.flexible_amides = kwargs.get('flexible_amides', self.flexible_amides)
        output_format = kwargs.get('output_format', 'string')

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

            # Step 2: Generate conformers if requested
            if num_conformers > 1:
                logger.info(f"Generating {num_conformers} conformers")
                conformers = self._generate_conformers(mol_3d, num_conformers)
            else:
                conformers = [mol_3d]

            # Step 3: Prepare each conformer with Meeko
            prepared_ligands = []
            for idx, conformer in enumerate(conformers):
                logger.info(f"Preparing conformer {idx + 1}/{len(conformers)} with Meeko")
                prepared = self._prepare_with_meeko(conformer)

                if prepared is None:
                    logger.warning(f"Failed to prepare conformer {idx + 1}, skipping")
                    continue

                prepared_ligands.append({
                    "conformer_id": idx,
                    "pdbqt_string": prepared["pdbqt_string"],
                    "metadata": prepared["metadata"]
                })

            if not prepared_ligands:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to prepare any conformers with Meeko"
                )

            # Step 4: Format output
            if output_format == 'file':
                # Write to temporary file
                temp_file = tempfile.NamedTemporaryFile(
                    mode='w',
                    suffix='.pdbqt',
                    delete=False,
                    prefix='meeko_ligand_'
                )
                temp_file.write(prepared_ligands[0]["pdbqt_string"])
                temp_file.close()

                result_data = {
                    "smiles": smiles,
                    "pdbqt_file": temp_file.name,
                    "num_conformers": len(prepared_ligands),
                    "conformers": prepared_ligands
                }
            else:
                # Return as string (default)
                result_data = {
                    "smiles": smiles,
                    "pdbqt_string": prepared_ligands[0]["pdbqt_string"],
                    "num_conformers": len(prepared_ligands),
                    "conformers": prepared_ligands
                }

            # Add summary metadata
            primary_metadata = prepared_ligands[0]["metadata"]

            logger.info(
                f"Ligand preparation complete: "
                f"{primary_metadata['num_torsions']} rotatable bonds, "
                f"{primary_metadata['num_macrocycles']} macrocycles"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "num_conformers": len(prepared_ligands),
                    "primary_conformer": primary_metadata,
                    "meeko_config": {
                        "rigid_macrocycles": self.rigid_macrocycles,
                        "keep_nonpolar_hydrogens": self.keep_nonpolar_hydrogens,
                        "hydrate": self.hydrate,
                        "flexible_amides": self.flexible_amides,
                        "min_ring_size": self.min_ring_size
                    }
                }
            )

        except Exception as e:
            logger.error(f"Meeko adapter error: {e}", exc_info=True)
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
            "description": "Ligand preparation for AutoDock Vina/GPU using Meeko",
            "requirements": [
                "Meeko (pip install meeko)",
                "RDKit (pip install rdkit)"
            ],
            "features": [
                "SMILES to PDBQT conversion",
                "Macrocycle flexibility support",
                "Multiple conformer generation",
                "Hydrated docking preparation",
                "Flexible amide bonds",
                "Custom bond rigidification"
            ],
            "config": {
                "rigid_macrocycles": self.rigid_macrocycles,
                "keep_nonpolar_hydrogens": self.keep_nonpolar_hydrogens,
                "hydrate": self.hydrate,
                "flexible_amides": self.flexible_amides,
                "min_ring_size": self.min_ring_size,
                "num_conformers": self.num_conformers,
                "double_bond_penalty": self.double_bond_penalty
            },
            "outputs": {
                "pdbqt_string": "PDBQT formatted ligand (string)",
                "pdbqt_file": "Optional temporary file path",
                "num_torsions": "Number of rotatable bonds",
                "macrocycles": "Detected macrocycle information",
                "conformers": "List of prepared conformer data"
            },
            "integration": {
                "vina_adapter": "Compatible with VinaAdapter for docking workflows",
                "format": "PDBQT format ready for AutoDock Vina/GPU"
            }
        }
