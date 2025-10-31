"""
ChemML Adapter - Machine learning workflows for chemical/materials data
Provides molecular representations and feature engineering for chemistry ML
"""
from typing import Any, Dict, Optional, List, Union
import logging
import numpy as np

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - required for ChemML")

try:
    from chemml.chem import Molecule
    from chemml.chem import CoulombMatrix as ChemMLCoulombMatrix
    from chemml.chem import BagofBonds
    from chemml.datasets import load_xyz_polarizability
    CHEMML_AVAILABLE = True
except (ImportError, RuntimeError) as e:
    CHEMML_AVAILABLE = False
    logging.warning(f"ChemML not available - {str(e)[:100]}... Install with: pip install chemml")
    # Define dummy classes to prevent NameError
    Molecule = None
    ChemMLCoulombMatrix = None
    BagofBonds = None

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ChemMLAdapter(AdapterProtocol):
    """
    Adapter for ChemML molecular representations and ML workflows
    Generates specialized chemical features for machine learning models

    Supports multiple representation methods:
    - coulomb_matrix: Coulomb matrix representation of molecular structure
    - bag_of_bonds: Bag of bonds representation
    - morgan_fingerprint: Morgan fingerprints via RDKit integration
    - rdkit_descriptors: Calculate RDKit descriptors with ChemML wrapper
    - feature_engineering: Automated feature engineering pipelines

    ChemML bridges chemistry and machine learning with specialized
    representations that capture molecular structure and properties.
    """

    SUPPORTED_OPERATIONS = {
        "coulomb_matrix": {
            "description": "Generate Coulomb matrix representation",
            "requires_3d": True,
            "default_params": {
                "max_n_atoms": 23,
                "sorting": "row_norm",
                "padding": 0.0
            }
        },
        "bag_of_bonds": {
            "description": "Generate bag of bonds representation",
            "requires_3d": True,
            "default_params": {
                "const": 1.0
            }
        },
        "morgan_fingerprint": {
            "description": "Generate Morgan fingerprints using RDKit",
            "requires_3d": False,
            "default_params": {
                "radius": 2,
                "n_bits": 2048
            }
        },
        "rdkit_descriptors": {
            "description": "Calculate RDKit descriptors",
            "requires_3d": False,
            "default_params": {}
        },
        "feature_standardization": {
            "description": "Standardize molecular features",
            "requires_3d": False,
            "default_params": {
                "with_mean": True,
                "with_std": True
            }
        }
    }

    def __init__(self):
        super().__init__(
            name="chemml",
            adapter_type="local",
            config={
                "timeout": 60,  # ChemML computations can take time
                "default_operation": "coulomb_matrix",
                "max_molecules": 1000,  # Limit for batch processing
                "default_max_atoms": 23  # Default for Coulomb matrix
            }
        )
        self.version = "1.0.0"

        # Check availability
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Required for ChemML. Install with: pip install rdkit")
        if not CHEMML_AVAILABLE:
            logger.error("ChemML is not installed! Install with: pip install chemml")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is valid for ChemML processing

        Args:
            input_data: Dictionary with 'smiles' and 'operation', or SMILES string/list

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE:
            return False

        # Handle dictionary input
        if isinstance(input_data, dict):
            if "smiles" not in input_data:
                return False
            smiles_data = input_data["smiles"]
        else:
            smiles_data = input_data

        # Validate SMILES string
        if isinstance(smiles_data, str):
            if len(smiles_data) == 0:
                return False
            try:
                mol = Chem.MolFromSmiles(smiles_data)
                return mol is not None
            except Exception:
                return False

        # Validate list of SMILES
        if isinstance(smiles_data, list):
            if len(smiles_data) == 0:
                return False
            for smiles in smiles_data:
                if not isinstance(smiles, str) or len(smiles) == 0:
                    return False
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        return False
                except Exception:
                    return False
            return True

        return False

    def _smiles_to_xyz(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to XYZ format with 3D coordinates

        Args:
            smiles: SMILES string

        Returns:
            XYZ format string or None if failed
        """
        if not RDKIT_AVAILABLE:
            return None

        try:
            from rdkit.Chem import AllChem

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                # Try with random coordinate generation
                AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)

            # Optimize geometry
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                # If MMFF fails, try UFF
                try:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                except Exception:
                    logger.warning(f"Could not optimize geometry for {smiles}")

            # Convert to XYZ format
            conf = mol.GetConformer()
            atoms = mol.GetAtoms()

            xyz_lines = [str(mol.GetNumAtoms())]
            xyz_lines.append(f"Generated from SMILES: {smiles}")

            for i, atom in enumerate(atoms):
                pos = conf.GetAtomPosition(i)
                symbol = atom.GetSymbol()
                xyz_lines.append(f"{symbol} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}")

            return "\n".join(xyz_lines)

        except Exception as e:
            logger.error(f"ChemML: Error converting SMILES to XYZ: {e}")
            return None

    def _generate_coulomb_matrix(
        self,
        smiles_list: List[str],
        max_n_atoms: int = 23,
        sorting: str = "row_norm",
        padding: float = 0.0
    ) -> Optional[Dict[str, Any]]:
        """
        Generate Coulomb matrix representations

        Args:
            smiles_list: List of SMILES strings
            max_n_atoms: Maximum number of atoms
            sorting: Sorting method ("row_norm", "unsorted", "random")
            padding: Padding value for smaller molecules

        Returns:
            Dictionary with representations and metadata
        """
        if not CHEMML_AVAILABLE or not RDKIT_AVAILABLE:
            raise Exception("ChemML and RDKit are required")

        try:
            representations = []
            failed_molecules = []

            for idx, smiles in enumerate(smiles_list):
                try:
                    # Convert SMILES to XYZ
                    xyz_str = self._smiles_to_xyz(smiles)
                    if xyz_str is None:
                        failed_molecules.append(idx)
                        representations.append(None)
                        continue

                    # Create ChemML Molecule object from XYZ string
                    mol = Molecule(xyz_str, input_type='xyz')

                    # Generate Coulomb matrix
                    cm = ChemMLCoulombMatrix(
                        max_n_atoms=max_n_atoms,
                        sorting=sorting,
                        padding=padding
                    )

                    cm_matrix = cm.represent(mol)

                    # Flatten the matrix for ML
                    cm_flat = cm_matrix.flatten()
                    representations.append(cm_flat)

                except Exception as e:
                    logger.warning(f"ChemML: Failed to generate Coulomb matrix for molecule {idx}: {e}")
                    failed_molecules.append(idx)
                    representations.append(None)

            # Filter out None values
            valid_representations = [r for r in representations if r is not None]

            if not valid_representations:
                return None

            # Convert to numpy array
            representations_array = np.array(valid_representations)

            return {
                "representations": representations_array.tolist(),
                "representation_type": "coulomb_matrix",
                "shape": list(representations_array.shape),
                "num_molecules": len(smiles_list),
                "num_successful": len(valid_representations),
                "failed_indices": failed_molecules,
                "feature_dimension": representations_array.shape[1] if len(representations_array.shape) > 1 else len(representations_array[0]),
                "parameters": {
                    "max_n_atoms": max_n_atoms,
                    "sorting": sorting,
                    "padding": padding
                }
            }

        except Exception as e:
            logger.error(f"ChemML: Error generating Coulomb matrices: {e}")
            return None

    def _generate_bag_of_bonds(
        self,
        smiles_list: List[str],
        const: float = 1.0
    ) -> Optional[Dict[str, Any]]:
        """
        Generate bag of bonds representations

        Args:
            smiles_list: List of SMILES strings
            const: Constant for BoB calculation

        Returns:
            Dictionary with representations and metadata
        """
        if not CHEMML_AVAILABLE or not RDKIT_AVAILABLE:
            raise Exception("ChemML and RDKit are required")

        try:
            representations = []
            failed_molecules = []

            for idx, smiles in enumerate(smiles_list):
                try:
                    # Convert SMILES to XYZ
                    xyz_str = self._smiles_to_xyz(smiles)
                    if xyz_str is None:
                        failed_molecules.append(idx)
                        representations.append(None)
                        continue

                    # Create ChemML Molecule object
                    mol = Molecule(xyz_str, input_type='xyz')

                    # Generate bag of bonds
                    bob = BagofBonds(const=const)
                    bob_vector = bob.represent(mol)

                    representations.append(bob_vector)

                except Exception as e:
                    logger.warning(f"ChemML: Failed to generate BoB for molecule {idx}: {e}")
                    failed_molecules.append(idx)
                    representations.append(None)

            # Filter out None values
            valid_representations = [r for r in representations if r is not None]

            if not valid_representations:
                return None

            # Pad representations to same length
            max_length = max(len(r) for r in valid_representations)
            padded_reps = []
            for r in valid_representations:
                if len(r) < max_length:
                    padded = np.pad(r, (0, max_length - len(r)), mode='constant', constant_values=0.0)
                    padded_reps.append(padded)
                else:
                    padded_reps.append(r)

            representations_array = np.array(padded_reps)

            return {
                "representations": representations_array.tolist(),
                "representation_type": "bag_of_bonds",
                "shape": list(representations_array.shape),
                "num_molecules": len(smiles_list),
                "num_successful": len(valid_representations),
                "failed_indices": failed_molecules,
                "feature_dimension": representations_array.shape[1] if len(representations_array.shape) > 1 else len(representations_array[0]),
                "parameters": {
                    "const": const
                }
            }

        except Exception as e:
            logger.error(f"ChemML: Error generating bag of bonds: {e}")
            return None

    def _generate_morgan_fingerprints(
        self,
        smiles_list: List[str],
        radius: int = 2,
        n_bits: int = 2048
    ) -> Optional[Dict[str, Any]]:
        """
        Generate Morgan fingerprints using RDKit

        Args:
            smiles_list: List of SMILES strings
            radius: Fingerprint radius
            n_bits: Number of bits in fingerprint

        Returns:
            Dictionary with representations and metadata
        """
        if not RDKIT_AVAILABLE:
            raise Exception("RDKit is required")

        try:
            from rdkit.Chem import AllChem

            representations = []
            failed_molecules = []

            for idx, smiles in enumerate(smiles_list):
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        failed_molecules.append(idx)
                        representations.append(None)
                        continue

                    # Generate Morgan fingerprint
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
                    fp_array = np.array(fp)
                    representations.append(fp_array)

                except Exception as e:
                    logger.warning(f"ChemML: Failed to generate Morgan FP for molecule {idx}: {e}")
                    failed_molecules.append(idx)
                    representations.append(None)

            # Filter out None values
            valid_representations = [r for r in representations if r is not None]

            if not valid_representations:
                return None

            representations_array = np.array(valid_representations)

            return {
                "representations": representations_array.tolist(),
                "representation_type": "morgan_fingerprint",
                "shape": list(representations_array.shape),
                "num_molecules": len(smiles_list),
                "num_successful": len(valid_representations),
                "failed_indices": failed_molecules,
                "feature_dimension": n_bits,
                "parameters": {
                    "radius": radius,
                    "n_bits": n_bits
                }
            }

        except Exception as e:
            logger.error(f"ChemML: Error generating Morgan fingerprints: {e}")
            return None

    def _calculate_rdkit_descriptors(
        self,
        smiles_list: List[str]
    ) -> Optional[Dict[str, Any]]:
        """
        Calculate RDKit descriptors for molecules

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Dictionary with descriptors and metadata
        """
        if not RDKIT_AVAILABLE:
            raise Exception("RDKit is required")

        try:
            from rdkit.Chem import Descriptors
            from rdkit.ML.Descriptors import MoleculeDescriptors

            # Get all available descriptors
            descriptor_names = [x[0] for x in Descriptors._descList]
            calc = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)

            descriptors = []
            failed_molecules = []

            for idx, smiles in enumerate(smiles_list):
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        failed_molecules.append(idx)
                        descriptors.append(None)
                        continue

                    # Calculate descriptors
                    desc_values = calc.CalcDescriptors(mol)
                    descriptors.append(desc_values)

                except Exception as e:
                    logger.warning(f"ChemML: Failed to calculate descriptors for molecule {idx}: {e}")
                    failed_molecules.append(idx)
                    descriptors.append(None)

            # Filter out None values
            valid_descriptors = [d for d in descriptors if d is not None]

            if not valid_descriptors:
                return None

            descriptors_array = np.array(valid_descriptors)

            return {
                "descriptors": descriptors_array.tolist(),
                "descriptor_names": descriptor_names,
                "representation_type": "rdkit_descriptors",
                "shape": list(descriptors_array.shape),
                "num_molecules": len(smiles_list),
                "num_successful": len(valid_descriptors),
                "failed_indices": failed_molecules,
                "num_descriptors": len(descriptor_names)
            }

        except Exception as e:
            logger.error(f"ChemML: Error calculating RDKit descriptors: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute ChemML operations

        Args:
            input_data: Dictionary with structure:
                {
                    "smiles": "CCO" or ["CCO", "CC(=O)O", ...],
                    "operation": "coulomb_matrix",
                    "parameters": {...}
                }
                Or simply a SMILES string/list (uses default operation)
            **kwargs: Additional parameters that override input_data["parameters"]

        Returns:
            AdapterResult containing molecular representations
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        if not CHEMML_AVAILABLE:
            # For operations that don't strictly need ChemML
            operation = None
            if isinstance(input_data, dict):
                operation = input_data.get("operation", self.config.get("default_operation"))

            if operation not in ["morgan_fingerprint", "rdkit_descriptors"]:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="ChemML is not installed. Install with: pip install chemml"
                )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input. Expected SMILES string(s) or dict with 'smiles' key"
            )

        # Parse input
        if isinstance(input_data, dict):
            smiles_data = input_data.get("smiles")
            operation = input_data.get("operation", self.config.get("default_operation"))
            parameters = input_data.get("parameters", {})
        else:
            smiles_data = input_data
            operation = kwargs.get("operation", self.config.get("default_operation"))
            parameters = {}

        # Override with kwargs
        parameters.update(kwargs)

        # Normalize to list
        if isinstance(smiles_data, str):
            smiles_list = [smiles_data]
            single_molecule = True
        else:
            smiles_list = smiles_data
            single_molecule = False

        # Check molecule limit
        max_molecules = self.config.get("max_molecules", 1000)
        if len(smiles_list) > max_molecules:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Too many molecules ({len(smiles_list)}). Maximum is {max_molecules}"
            )

        # Validate operation
        if operation not in self.SUPPORTED_OPERATIONS:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported operation: {operation}. "
                      f"Supported: {list(self.SUPPORTED_OPERATIONS.keys())}"
            )

        # Execute operation
        try:
            if operation == "coulomb_matrix":
                max_n_atoms = parameters.get("max_n_atoms", 23)
                sorting = parameters.get("sorting", "row_norm")
                padding = parameters.get("padding", 0.0)

                result_data = self._generate_coulomb_matrix(
                    smiles_list,
                    max_n_atoms=max_n_atoms,
                    sorting=sorting,
                    padding=padding
                )

            elif operation == "bag_of_bonds":
                const = parameters.get("const", 1.0)

                result_data = self._generate_bag_of_bonds(
                    smiles_list,
                    const=const
                )

            elif operation == "morgan_fingerprint":
                radius = parameters.get("radius", 2)
                n_bits = parameters.get("n_bits", 2048)

                result_data = self._generate_morgan_fingerprints(
                    smiles_list,
                    radius=radius,
                    n_bits=n_bits
                )

            elif operation == "rdkit_descriptors":
                result_data = self._calculate_rdkit_descriptors(smiles_list)

            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Operation {operation} not yet implemented"
                )

        except Exception as e:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to execute {operation}: {str(e)}",
                metadata={
                    "source": "chemml",
                    "operation": operation
                }
            )

        if result_data is None:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to generate {operation} representations",
                metadata={
                    "source": "chemml",
                    "operation": operation,
                    "num_molecules": len(smiles_list)
                }
            )

        # Add input SMILES to result
        result_data["smiles"] = smiles_list[0] if single_molecule else smiles_list

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "chemml",
                "adapter_version": self.version,
                "computation_type": "local",
                "operation": operation,
                "num_molecules": len(smiles_list),
                "single_molecule": single_molecule,
                "chemml_available": CHEMML_AVAILABLE,
                "rdkit_available": RDKIT_AVAILABLE
            }
        )
