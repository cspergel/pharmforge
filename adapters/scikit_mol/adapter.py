"""
scikit-mol Adapter - Generates molecular fingerprints and descriptors for scikit-learn ML
Bridges RDKit and scikit-learn for molecular machine learning workflows
"""
from typing import Any, Dict, Optional, List
import logging
import numpy as np

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - required for scikit-mol")

try:
    from scikit_mol.fingerprints import (
        MorganFingerprintTransformer,
        MACCSKeysFingerprintTransformer,
        RDKitFingerprintTransformer,
        AtomPairFingerprintTransformer,
        TopologicalTorsionFingerprintTransformer
    )
    from scikit_mol.descriptors import MolecularDescriptorTransformer
    SCIKIT_MOL_AVAILABLE = True
except ImportError:
    SCIKIT_MOL_AVAILABLE = False
    logging.warning("scikit-mol not available - install with: pip install scikit-mol")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ScikitMolAdapter(AdapterProtocol):
    """
    Adapter for scikit-mol molecular fingerprints and descriptors
    Generates feature vectors compatible with scikit-learn ML pipelines

    Supports multiple fingerprint types:
    - Morgan (circular) fingerprints
    - MACCS keys (166-bit structural keys)
    - RDKit fingerprints
    - Atom pair fingerprints
    - Topological torsion fingerprints
    - Molecular descriptors
    """

    # Mapping of fingerprint types to their default parameters
    FINGERPRINT_TYPES = {
        "morgan": {"radius": 2, "nBits": 2048},
        "maccs": {},
        "rdkit": {"fpSize": 2048},
        "atom_pair": {"nBits": 2048},
        "topological_torsion": {"nBits": 2048},
        "descriptors": {}
    }

    def __init__(self):
        super().__init__(
            name="scikit_mol",
            adapter_type="local",
            config={
                "timeout": 30,  # Local calculations should be reasonably fast
                "default_fingerprint": "morgan",  # Default to Morgan fingerprints
                "default_radius": 2,  # Morgan fingerprint radius
                "default_nbits": 2048  # Default fingerprint size
            }
        )
        self.version = "1.0.0"

        # Check availability
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Required for scikit-mol. Install with: pip install rdkit")
        if not SCIKIT_MOL_AVAILABLE:
            logger.error("scikit-mol is not installed! Install with: pip install scikit-mol")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string or list of SMILES

        Args:
            input_data: SMILES string or list of SMILES strings to validate

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE:
            return False

        # Handle single SMILES string
        if isinstance(input_data, str):
            if len(input_data) == 0:
                return False
            try:
                mol = Chem.MolFromSmiles(input_data)
                return mol is not None
            except Exception:
                return False

        # Handle list of SMILES strings
        if isinstance(input_data, list):
            if len(input_data) == 0:
                return False
            # Validate all SMILES in the list
            for smiles in input_data:
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

    def _get_fingerprint_transformer(
        self,
        fingerprint_type: str,
        **kwargs
    ):
        """
        Get the appropriate scikit-mol fingerprint transformer

        Args:
            fingerprint_type: Type of fingerprint ("morgan", "maccs", "rdkit", etc.)
            **kwargs: Fingerprint-specific parameters

        Returns:
            scikit-mol transformer instance

        Raises:
            ValueError: If fingerprint type is not supported
        """
        if not SCIKIT_MOL_AVAILABLE:
            raise Exception("scikit-mol is not installed")

        # Get default parameters for this fingerprint type
        default_params = self.FINGERPRINT_TYPES.get(fingerprint_type, {})

        # Merge with user-provided kwargs
        params = {**default_params, **kwargs}

        # Create appropriate transformer
        if fingerprint_type == "morgan":
            radius = params.get("radius", self.config.get("default_radius", 2))
            nBits = params.get("nBits", self.config.get("default_nbits", 2048))
            return MorganFingerprintTransformer(radius=radius, nBits=nBits)

        elif fingerprint_type == "maccs":
            return MACCSKeysFingerprintTransformer()

        elif fingerprint_type == "rdkit":
            fpSize = params.get("fpSize", self.config.get("default_nbits", 2048))
            return RDKitFingerprintTransformer(fpSize=fpSize)

        elif fingerprint_type == "atom_pair":
            nBits = params.get("nBits", self.config.get("default_nbits", 2048))
            return AtomPairFingerprintTransformer(nBits=nBits)

        elif fingerprint_type == "topological_torsion":
            nBits = params.get("nBits", self.config.get("default_nbits", 2048))
            return TopologicalTorsionFingerprintTransformer(nBits=nBits)

        elif fingerprint_type == "descriptors":
            return MolecularDescriptorTransformer()

        else:
            raise ValueError(
                f"Unsupported fingerprint type: {fingerprint_type}. "
                f"Supported types: {list(self.FINGERPRINT_TYPES.keys())}"
            )

    def _generate_fingerprints(
        self,
        smiles: List[str],
        fingerprint_type: str,
        **kwargs
    ) -> Optional[np.ndarray]:
        """
        Generate fingerprints for a list of SMILES strings

        Args:
            smiles: List of SMILES strings
            fingerprint_type: Type of fingerprint to generate
            **kwargs: Fingerprint-specific parameters

        Returns:
            NumPy array of fingerprints (n_molecules x n_features)
        """
        if not RDKIT_AVAILABLE or not SCIKIT_MOL_AVAILABLE:
            raise Exception("RDKit and scikit-mol are required")

        try:
            # Convert SMILES to RDKit molecules
            mols = [Chem.MolFromSmiles(s) for s in smiles]

            # Check for invalid molecules
            if None in mols:
                invalid_indices = [i for i, mol in enumerate(mols) if mol is None]
                logger.warning(
                    f"scikit-mol: Could not parse {len(invalid_indices)} SMILES strings: "
                    f"indices {invalid_indices[:5]}{'...' if len(invalid_indices) > 5 else ''}"
                )
                # Filter out None values
                mols = [mol for mol in mols if mol is not None]
                if not mols:
                    return None

            # Get transformer
            transformer = self._get_fingerprint_transformer(fingerprint_type, **kwargs)

            # Generate fingerprints
            fingerprints = transformer.transform(mols)

            return fingerprints

        except Exception as e:
            logger.error(f"scikit-mol: Error generating fingerprints: {e}")
            return None

    def _get_feature_info(
        self,
        fingerprint_type: str,
        fingerprints: np.ndarray
    ) -> Dict[str, Any]:
        """
        Get metadata about the generated features

        Args:
            fingerprint_type: Type of fingerprint generated
            fingerprints: Generated fingerprint array

        Returns:
            Dictionary with feature information
        """
        return {
            "fingerprint_type": fingerprint_type,
            "n_molecules": fingerprints.shape[0],
            "n_features": fingerprints.shape[1],
            "feature_shape": fingerprints.shape,
            "dtype": str(fingerprints.dtype),
            "is_binary": np.array_equal(fingerprints, fingerprints.astype(bool)),
            "sparsity": 1.0 - (np.count_nonzero(fingerprints) / fingerprints.size)
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute scikit-mol fingerprint/descriptor generation

        Args:
            input_data: SMILES string or list of SMILES strings
            **kwargs: Additional parameters:
                - fingerprint_type: Type of fingerprint ("morgan", "maccs", "rdkit",
                  "atom_pair", "topological_torsion", "descriptors")
                - radius: Morgan fingerprint radius (default: 2)
                - nBits: Fingerprint size (default: 2048)
                - return_array: If True, return NumPy array; if False, return list
                  (default: False for compatibility)

        Returns:
            AdapterResult containing feature vectors
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        if not SCIKIT_MOL_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="scikit-mol is not installed. Install with: pip install scikit-mol"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string(s) or could not parse with RDKit"
            )

        # Normalize input to list
        if isinstance(input_data, str):
            smiles_list = [input_data]
            single_molecule = True
        else:
            smiles_list = input_data
            single_molecule = False

        # Get fingerprint type from kwargs
        fingerprint_type = kwargs.get(
            "fingerprint_type",
            self.config.get("default_fingerprint", "morgan")
        )

        # Validate fingerprint type
        if fingerprint_type not in self.FINGERPRINT_TYPES:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported fingerprint type: {fingerprint_type}. "
                      f"Supported types: {list(self.FINGERPRINT_TYPES.keys())}"
            )

        # Extract fingerprint parameters from kwargs
        fp_params = {}
        if fingerprint_type == "morgan":
            if "radius" in kwargs:
                fp_params["radius"] = kwargs["radius"]
            if "nBits" in kwargs:
                fp_params["nBits"] = kwargs["nBits"]
        elif fingerprint_type in ["rdkit"]:
            if "fpSize" in kwargs:
                fp_params["fpSize"] = kwargs["fpSize"]
            elif "nBits" in kwargs:
                fp_params["fpSize"] = kwargs["nBits"]
        elif fingerprint_type in ["atom_pair", "topological_torsion"]:
            if "nBits" in kwargs:
                fp_params["nBits"] = kwargs["nBits"]

        # Generate fingerprints
        try:
            fingerprints = self._generate_fingerprints(
                smiles_list,
                fingerprint_type,
                **fp_params
            )
        except Exception as e:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to generate fingerprints: {str(e)}",
                metadata={
                    "source": "scikit_mol",
                    "fingerprint_type": fingerprint_type
                }
            )

        if fingerprints is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to generate fingerprints with scikit-mol",
                metadata={
                    "source": "scikit_mol",
                    "smiles": smiles_list[:5],  # First 5 for debugging
                    "fingerprint_type": fingerprint_type
                }
            )

        # Get feature metadata
        feature_info = self._get_feature_info(fingerprint_type, fingerprints)

        # Prepare result data
        return_array = kwargs.get("return_array", False)

        if return_array:
            # Return as NumPy array (useful for ML pipelines)
            result_data = {
                "fingerprints": fingerprints,
                "smiles": smiles_list,
                **feature_info
            }
        else:
            # Return as list (default, more JSON-friendly)
            if single_molecule:
                # For single molecule, return 1D array/list
                result_data = {
                    "fingerprint": fingerprints[0].tolist(),
                    "smiles": smiles_list[0],
                    **feature_info
                }
            else:
                # For multiple molecules, return 2D array/list
                result_data = {
                    "fingerprints": fingerprints.tolist(),
                    "smiles": smiles_list,
                    **feature_info
                }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "scikit_mol",
                "adapter_version": self.version,
                "computation_type": "local",
                "fingerprint_type": fingerprint_type,
                "n_molecules": len(smiles_list),
                "parameters": fp_params
            }
        )
