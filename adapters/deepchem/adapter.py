"""
DeepChem Adapter - Generates molecular features using DeepChem featurizers
Provides deep learning-ready molecular representations for chemistry/biology ML
"""
from typing import Any, Dict, Optional, List
import logging
import numpy as np

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - required for DeepChem")

try:
    import deepchem as dc
    from deepchem.feat import (
        CircularFingerprint,
        RDKitDescriptors,
        MorganFingerprint,
        MACCSKeysFingerprint,
        CoulombMatrix,
        CoulombMatrixEig,
        AtomicCoordinates,
        BPSymmetryFunctionInput,
        WeaveFeaturizer,
        ConvMolFeaturizer,
        MolGraphConvFeaturizer,
    )
    DEEPCHEM_AVAILABLE = True
except (ImportError, RuntimeError) as e:
    DEEPCHEM_AVAILABLE = False
    logging.warning(f"DeepChem not available - {str(e)[:100]}... This is expected if there are PyTorch/torchvision compatibility issues")
    # Define dummy classes to prevent NameError in class definition
    CircularFingerprint = None
    RDKitDescriptors = None
    MorganFingerprint = None
    MACCSKeysFingerprint = None
    CoulombMatrix = None
    CoulombMatrixEig = None
    AtomicCoordinates = None
    BPSymmetryFunctionInput = None
    WeaveFeaturizer = None
    ConvMolFeaturizer = None
    MolGraphConvFeaturizer = None
    dc = None

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DeepChemAdapter(AdapterProtocol):
    """
    Adapter for DeepChem molecular featurization
    Generates feature vectors suitable for deep learning models in chemistry/biology

    Supports multiple featurization methods:
    - CircularFingerprint: Circular fingerprints (similar to ECFP)
    - MorganFingerprint: Morgan fingerprints with customizable radius
    - MACCSKeys: 166-bit MACCS structural keys
    - RDKitDescriptors: 200+ physicochemical descriptors
    - WeaveFeaturizer: Graph convolution features for Weave models
    - ConvMolFeaturizer: Graph convolution features for GraphConv models
    - MolGraphConvFeaturizer: Molecular graph features for GCN models
    - CoulombMatrix: Coulomb matrix representation (requires 3D coordinates)
    - CoulombMatrixEig: Eigenvalues of Coulomb matrix
    """

    # Mapping of featurizer types to their classes and default parameters
    FEATURIZER_TYPES = {
        "circular": {
            "class": CircularFingerprint,
            "params": {"radius": 2, "size": 2048},
            "requires_3d": False
        },
        "morgan": {
            "class": MorganFingerprint,
            "params": {"radius": 2, "size": 2048},
            "requires_3d": False
        },
        "maccs": {
            "class": MACCSKeysFingerprint,
            "params": {},
            "requires_3d": False
        },
        "rdkit_descriptors": {
            "class": RDKitDescriptors,
            "params": {},
            "requires_3d": False
        },
        "weave": {
            "class": WeaveFeaturizer,
            "params": {},
            "requires_3d": False
        },
        "graph_conv": {
            "class": ConvMolFeaturizer,
            "params": {},
            "requires_3d": False
        },
        "mol_graph_conv": {
            "class": MolGraphConvFeaturizer,
            "params": {"use_edges": True},
            "requires_3d": False
        },
        "coulomb_matrix": {
            "class": CoulombMatrix,
            "params": {"max_atoms": 100},
            "requires_3d": True
        },
        "coulomb_matrix_eig": {
            "class": CoulombMatrixEig,
            "params": {"max_atoms": 100},
            "requires_3d": True
        },
    }

    def __init__(self):
        super().__init__(
            name="deepchem",
            adapter_type="local",
            config={
                "timeout": 60,  # Deep learning featurization can be slower
                "default_featurizer": "morgan",  # Default to Morgan fingerprints
                "default_radius": 2,  # Morgan/Circular fingerprint radius
                "default_size": 2048,  # Default fingerprint size
                "max_atoms": 100  # Maximum atoms for 3D representations
            }
        )
        self.version = "1.0.0"

        # Check availability
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Required for DeepChem. Install with: pip install rdkit")
        if not DEEPCHEM_AVAILABLE:
            logger.error("DeepChem is not installed! Install with: pip install deepchem")

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

    def _get_featurizer(
        self,
        featurizer_type: str,
        **kwargs
    ):
        """
        Get the appropriate DeepChem featurizer instance

        Args:
            featurizer_type: Type of featurizer ("morgan", "circular", "rdkit_descriptors", etc.)
            **kwargs: Featurizer-specific parameters

        Returns:
            DeepChem featurizer instance

        Raises:
            ValueError: If featurizer type is not supported
        """
        if not DEEPCHEM_AVAILABLE:
            raise Exception("DeepChem is not installed")

        if featurizer_type not in self.FEATURIZER_TYPES:
            raise ValueError(
                f"Unsupported featurizer type: {featurizer_type}. "
                f"Supported types: {list(self.FEATURIZER_TYPES.keys())}"
            )

        # Get featurizer configuration
        featurizer_config = self.FEATURIZER_TYPES[featurizer_type]
        featurizer_class = featurizer_config["class"]
        default_params = featurizer_config["params"].copy()

        # Merge with user-provided kwargs
        params = {**default_params, **kwargs}

        # Filter out None values
        params = {k: v for k, v in params.items() if v is not None}

        # Create featurizer instance
        try:
            featurizer = featurizer_class(**params)
            return featurizer
        except Exception as e:
            logger.error(f"DeepChem: Error creating featurizer {featurizer_type}: {e}")
            raise

    def _featurize_molecules(
        self,
        smiles: List[str],
        featurizer_type: str,
        **kwargs
    ) -> Optional[np.ndarray]:
        """
        Generate features for a list of SMILES strings using DeepChem

        Args:
            smiles: List of SMILES strings
            featurizer_type: Type of featurizer to use
            **kwargs: Featurizer-specific parameters

        Returns:
            NumPy array of features or list of feature objects (depends on featurizer)
        """
        if not RDKIT_AVAILABLE or not DEEPCHEM_AVAILABLE:
            raise Exception("RDKit and DeepChem are required")

        try:
            # Get featurizer
            featurizer = self._get_featurizer(featurizer_type, **kwargs)

            # Featurize molecules
            features = featurizer.featurize(smiles)

            # Handle failed featurizations (DeepChem returns None for failures)
            if features is None:
                logger.warning("DeepChem: Featurization returned None")
                return None

            # Count failures
            none_count = sum(1 for f in features if f is None)
            if none_count > 0:
                logger.warning(
                    f"DeepChem: {none_count}/{len(features)} molecules failed featurization"
                )

            return features

        except Exception as e:
            logger.error(f"DeepChem: Error during featurization: {e}")
            return None

    def _get_feature_info(
        self,
        featurizer_type: str,
        features: Any
    ) -> Dict[str, Any]:
        """
        Get metadata about the generated features

        Args:
            featurizer_type: Type of featurizer used
            features: Generated features

        Returns:
            Dictionary with feature information
        """
        info = {
            "featurizer_type": featurizer_type,
            "requires_3d": self.FEATURIZER_TYPES[featurizer_type]["requires_3d"]
        }

        # Handle different feature types
        if isinstance(features, np.ndarray):
            info.update({
                "n_molecules": features.shape[0] if features.ndim > 0 else 1,
                "feature_shape": features.shape if features.ndim > 1 else (1, features.shape[0]),
                "dtype": str(features.dtype),
                "is_array": True
            })

            # Additional stats for numeric arrays
            if features.ndim > 1:
                info["n_features"] = features.shape[1]
                # Check if binary
                if features.dtype in [np.int32, np.int64, np.float32, np.float64]:
                    unique_vals = np.unique(features)
                    info["is_binary"] = len(unique_vals) == 2 and set(unique_vals).issubset({0, 1})
                    if info["is_binary"]:
                        info["sparsity"] = 1.0 - (np.count_nonzero(features) / features.size)

        elif isinstance(features, list):
            info.update({
                "n_molecules": len(features),
                "is_array": False,
                "is_graph": featurizer_type in ["weave", "graph_conv", "mol_graph_conv"]
            })

            # For graph-based features, get info from first non-None feature
            if info["is_graph"]:
                valid_features = [f for f in features if f is not None]
                if valid_features:
                    first_feature = valid_features[0]
                    info["feature_type"] = type(first_feature).__name__

        return info

    def _convert_features_to_serializable(
        self,
        features: Any,
        featurizer_type: str,
        single_molecule: bool
    ) -> Any:
        """
        Convert DeepChem features to JSON-serializable format

        Args:
            features: Raw DeepChem features
            featurizer_type: Type of featurizer used
            single_molecule: Whether this is a single molecule

        Returns:
            JSON-serializable representation of features
        """
        # For array-based features (fingerprints, descriptors)
        if isinstance(features, np.ndarray):
            if single_molecule:
                return features[0].tolist() if features.ndim > 1 else features.tolist()
            else:
                return features.tolist()

        # For list-based features (graph representations)
        if isinstance(features, list):
            # Graph-based features are not easily serializable
            # Return a summary instead
            if featurizer_type in ["weave", "graph_conv", "mol_graph_conv"]:
                return {
                    "feature_type": "graph",
                    "note": "Graph features are not JSON-serializable. Use return_raw=True to get raw features.",
                    "n_features": len(features),
                    "valid_features": sum(1 for f in features if f is not None)
                }
            else:
                # Try to convert to list
                try:
                    if single_molecule:
                        return features[0].tolist() if hasattr(features[0], 'tolist') else features[0]
                    else:
                        return [f.tolist() if hasattr(f, 'tolist') else f for f in features]
                except Exception:
                    return str(features)

        return features

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute DeepChem featurization

        Args:
            input_data: SMILES string or list of SMILES strings
            **kwargs: Additional parameters:
                - featurizer_type: Type of featurizer ("morgan", "circular", "rdkit_descriptors",
                  "maccs", "weave", "graph_conv", "mol_graph_conv", "coulomb_matrix", etc.)
                - radius: Fingerprint radius for Morgan/Circular (default: 2)
                - size: Fingerprint size for Morgan/Circular (default: 2048)
                - max_atoms: Maximum atoms for 3D representations (default: 100)
                - return_raw: If True, return raw DeepChem features; if False, return
                  JSON-serializable format (default: False)

        Returns:
            AdapterResult containing molecular features
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        if not DEEPCHEM_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="DeepChem is not installed. Install with: pip install deepchem"
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

        # Get featurizer type from kwargs
        featurizer_type = kwargs.get(
            "featurizer_type",
            self.config.get("default_featurizer", "morgan")
        )

        # Validate featurizer type
        if featurizer_type not in self.FEATURIZER_TYPES:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported featurizer type: {featurizer_type}. "
                      f"Supported types: {list(self.FEATURIZER_TYPES.keys())}"
            )

        # Extract featurizer parameters from kwargs
        featurizer_params = {}

        # Parameters for fingerprint-based featurizers
        if featurizer_type in ["morgan", "circular"]:
            if "radius" in kwargs:
                featurizer_params["radius"] = kwargs["radius"]
            if "size" in kwargs:
                featurizer_params["size"] = kwargs["size"]

        # Parameters for 3D featurizers
        if featurizer_type in ["coulomb_matrix", "coulomb_matrix_eig"]:
            if "max_atoms" in kwargs:
                featurizer_params["max_atoms"] = kwargs["max_atoms"]

        # Parameters for graph featurizers
        if featurizer_type == "mol_graph_conv":
            if "use_edges" in kwargs:
                featurizer_params["use_edges"] = kwargs["use_edges"]

        # Generate features
        try:
            features = self._featurize_molecules(
                smiles_list,
                featurizer_type,
                **featurizer_params
            )
        except Exception as e:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to generate features: {str(e)}",
                metadata={
                    "source": "deepchem",
                    "featurizer_type": featurizer_type
                }
            )

        if features is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to generate features with DeepChem",
                metadata={
                    "source": "deepchem",
                    "smiles": smiles_list[:5],  # First 5 for debugging
                    "featurizer_type": featurizer_type
                }
            )

        # Get feature metadata
        feature_info = self._get_feature_info(featurizer_type, features)

        # Prepare result data
        return_raw = kwargs.get("return_raw", False)

        if return_raw:
            # Return raw DeepChem features (may not be JSON-serializable)
            result_data = {
                "features": features,
                "smiles": smiles_list if not single_molecule else smiles_list[0],
                **feature_info
            }
        else:
            # Return JSON-serializable format
            serializable_features = self._convert_features_to_serializable(
                features,
                featurizer_type,
                single_molecule
            )

            if single_molecule:
                result_data = {
                    "features": serializable_features,
                    "smiles": smiles_list[0],
                    **feature_info
                }
            else:
                result_data = {
                    "features": serializable_features,
                    "smiles": smiles_list,
                    **feature_info
                }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "deepchem",
                "adapter_version": self.version,
                "computation_type": "local",
                "featurizer_type": featurizer_type,
                "n_molecules": len(smiles_list),
                "parameters": featurizer_params,
                "deepchem_version": dc.__version__ if DEEPCHEM_AVAILABLE else None
            }
        )
