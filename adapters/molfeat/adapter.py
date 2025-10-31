"""
MolFeat Adapter - Unified molecular featurization library
Provides 100+ featurizers with pre-trained models and embeddings
"""
from typing import Any, Dict, Optional, List, Union
import logging
import numpy as np

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - required for MolFeat. Install with: pip install rdkit")

try:
    from molfeat.calc import FPCalculator
    from molfeat.trans import MoleculeTransformer
    from molfeat.trans.pretrained import PretrainedMolTransformer
    MOLFEAT_AVAILABLE = True
except ImportError:
    MOLFEAT_AVAILABLE = False
    logging.warning("MolFeat not available - install with: pip install molfeat")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class MolFeatAdapter(AdapterProtocol):
    """
    Adapter for MolFeat molecular featurization
    Provides unified access to 100+ molecular featurizers including:
    - 2D fingerprints (Morgan, MACCS, RDKit, etc.)
    - 3D descriptors (USR, USRCAT, etc.)
    - Learned representations (pre-trained models)
    - Physicochemical descriptors

    Supports both single molecules and batches with pre-trained model embeddings
    from MolFeat's extensive model zoo.
    """

    # Common featurizer types with their default configurations
    FEATURIZER_TYPES = {
        # 2D Fingerprints
        "morgan": {
            "calculator": "FPCalculator",
            "params": {"name": "morgan", "radius": 2, "nBits": 2048},
            "requires_3d": False,
            "description": "Morgan (circular) fingerprints"
        },
        "ecfp": {
            "calculator": "FPCalculator",
            "params": {"name": "ecfp", "radius": 2, "nBits": 2048},
            "requires_3d": False,
            "description": "Extended Connectivity Fingerprints"
        },
        "fcfp": {
            "calculator": "FPCalculator",
            "params": {"name": "fcfp", "radius": 2, "nBits": 2048},
            "requires_3d": False,
            "description": "Functional Connectivity Fingerprints"
        },
        "maccs": {
            "calculator": "FPCalculator",
            "params": {"name": "maccs"},
            "requires_3d": False,
            "description": "MACCS keys (166-bit structural keys)"
        },
        "rdkit": {
            "calculator": "FPCalculator",
            "params": {"name": "rdkit", "nBits": 2048},
            "requires_3d": False,
            "description": "RDKit fingerprints"
        },
        "atompair": {
            "calculator": "FPCalculator",
            "params": {"name": "atompair", "nBits": 2048},
            "requires_3d": False,
            "description": "Atom pair fingerprints"
        },
        "topological": {
            "calculator": "FPCalculator",
            "params": {"name": "topological", "nBits": 2048},
            "requires_3d": False,
            "description": "Topological torsion fingerprints"
        },
        "avalon": {
            "calculator": "FPCalculator",
            "params": {"name": "avalon", "nBits": 2048},
            "requires_3d": False,
            "description": "Avalon fingerprints"
        },
        "erg": {
            "calculator": "FPCalculator",
            "params": {"name": "erg"},
            "requires_3d": False,
            "description": "ErG fingerprints"
        },
        "estate": {
            "calculator": "FPCalculator",
            "params": {"name": "estate"},
            "requires_3d": False,
            "description": "EState fingerprints"
        },
        # Descriptors
        "desc2d": {
            "calculator": "FPCalculator",
            "params": {"name": "desc2d"},
            "requires_3d": False,
            "description": "2D physicochemical descriptors"
        },
        "desc3d": {
            "calculator": "FPCalculator",
            "params": {"name": "desc3d"},
            "requires_3d": True,
            "description": "3D physicochemical descriptors"
        },
        # 3D Shape-based
        "usr": {
            "calculator": "FPCalculator",
            "params": {"name": "usr"},
            "requires_3d": True,
            "description": "Ultrafast Shape Recognition"
        },
        "usrcat": {
            "calculator": "FPCalculator",
            "params": {"name": "usrcat"},
            "requires_3d": True,
            "description": "USR with CREDO Atom Types"
        },
        "electroshape": {
            "calculator": "FPCalculator",
            "params": {"name": "electroshape"},
            "requires_3d": True,
            "description": "Electroshape descriptors"
        },
    }

    # Pre-trained model names (examples from MolFeat model zoo)
    PRETRAINED_MODELS = {
        "ChemBERTa-77M-MLM": "ChemBERTa pre-trained on 77M molecules",
        "ChemBERTa-77M-MTR": "ChemBERTa multi-task regression",
        "gin_supervised_contextpred": "GIN supervised with context prediction",
        "gin_supervised_infomax": "GIN supervised with infomax",
        "gin_supervised_edgepred": "GIN supervised with edge prediction",
        "gin_supervised_masking": "GIN supervised with masking",
        "MolT5": "MolT5 transformer model",
        "ChemGPT-1.2B": "ChemGPT large language model",
    }

    def __init__(self):
        super().__init__(
            name="molfeat",
            adapter_type="local",
            config={
                "timeout": 120,  # Pre-trained models may need more time
                "default_featurizer": "morgan",
                "default_radius": 2,
                "default_nbits": 2048,
                "batch_size": 32,  # Batch size for featurization
            }
        )
        self.version = "1.0.0"

        # Check availability
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Required for MolFeat. Install with: pip install rdkit")
        if not MOLFEAT_AVAILABLE:
            logger.error("MolFeat is not installed! Install with: pip install molfeat")

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
        pretrained_model: Optional[str] = None,
        **kwargs
    ):
        """
        Get the appropriate MolFeat featurizer or transformer

        Args:
            featurizer_type: Type of featurizer or "pretrained" for pre-trained models
            pretrained_model: Name of pre-trained model (if featurizer_type="pretrained")
            **kwargs: Featurizer-specific parameters

        Returns:
            MolFeat featurizer or transformer instance

        Raises:
            ValueError: If featurizer type is not supported
        """
        if not MOLFEAT_AVAILABLE:
            raise Exception("MolFeat is not installed")

        # Handle pre-trained models
        if featurizer_type == "pretrained":
            if pretrained_model is None:
                raise ValueError("pretrained_model must be specified when featurizer_type='pretrained'")

            try:
                logger.info(f"Loading pre-trained model: {pretrained_model}")
                transformer = PretrainedMolTransformer(kind=pretrained_model, **kwargs)
                return transformer
            except Exception as e:
                logger.error(f"MolFeat: Error loading pre-trained model {pretrained_model}: {e}")
                raise

        # Handle standard featurizers
        if featurizer_type not in self.FEATURIZER_TYPES:
            raise ValueError(
                f"Unsupported featurizer type: {featurizer_type}. "
                f"Supported types: {list(self.FEATURIZER_TYPES.keys())} or 'pretrained'"
            )

        # Get featurizer configuration
        featurizer_config = self.FEATURIZER_TYPES[featurizer_type]
        calculator_type = featurizer_config["calculator"]
        default_params = featurizer_config["params"].copy()

        # Merge with user-provided kwargs
        params = {**default_params, **kwargs}

        # Filter out None values and non-FPCalculator params
        params = {k: v for k, v in params.items() if v is not None}

        # Create featurizer instance
        try:
            if calculator_type == "FPCalculator":
                # FPCalculator is for fingerprints and descriptors
                calculator = FPCalculator(**params)
                logger.info(f"Created FPCalculator with params: {params}")
                return calculator
            else:
                # Use MoleculeTransformer for other types
                transformer = MoleculeTransformer(featurizer=featurizer_type, **params)
                logger.info(f"Created MoleculeTransformer for {featurizer_type}")
                return transformer
        except Exception as e:
            logger.error(f"MolFeat: Error creating featurizer {featurizer_type}: {e}")
            raise

    def _featurize_molecules(
        self,
        smiles: List[str],
        featurizer_type: str,
        pretrained_model: Optional[str] = None,
        **kwargs
    ) -> Optional[np.ndarray]:
        """
        Generate features for a list of SMILES strings using MolFeat

        Args:
            smiles: List of SMILES strings
            featurizer_type: Type of featurizer to use
            pretrained_model: Name of pre-trained model (if using pretrained)
            **kwargs: Featurizer-specific parameters

        Returns:
            NumPy array of features or list of feature objects
        """
        if not RDKIT_AVAILABLE or not MOLFEAT_AVAILABLE:
            raise Exception("RDKit and MolFeat are required")

        try:
            # Get featurizer
            featurizer = self._get_featurizer(
                featurizer_type,
                pretrained_model=pretrained_model,
                **kwargs
            )

            # Featurize molecules
            logger.info(f"Featurizing {len(smiles)} molecules with {featurizer_type}")

            # MolFeat can handle both single and batch inputs
            if isinstance(featurizer, PretrainedMolTransformer):
                # Pre-trained models use transform
                features = featurizer(smiles)
            elif hasattr(featurizer, '__call__'):
                # FPCalculator and other calculators are callable
                features = featurizer(smiles)
            else:
                raise Exception(f"Featurizer {featurizer_type} is not callable")

            if features is None:
                logger.warning("MolFeat: Featurization returned None")
                return None

            # Convert to numpy array if not already
            if not isinstance(features, np.ndarray):
                try:
                    features = np.array(features)
                except Exception as e:
                    logger.warning(f"Could not convert features to numpy array: {e}")

            return features

        except Exception as e:
            logger.error(f"MolFeat: Error during featurization: {e}")
            return None

    def _get_feature_info(
        self,
        featurizer_type: str,
        features: Any,
        pretrained_model: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Get metadata about the generated features

        Args:
            featurizer_type: Type of featurizer used
            features: Generated features
            pretrained_model: Name of pre-trained model (if used)

        Returns:
            Dictionary with feature information
        """
        info = {
            "featurizer_type": featurizer_type,
            "requires_3d": self.FEATURIZER_TYPES.get(featurizer_type, {}).get("requires_3d", False),
        }

        if pretrained_model:
            info["pretrained_model"] = pretrained_model
            info["model_description"] = self.PRETRAINED_MODELS.get(
                pretrained_model,
                "Custom pre-trained model"
            )

        # Handle numpy arrays
        if isinstance(features, np.ndarray):
            info.update({
                "n_molecules": features.shape[0] if features.ndim > 0 else 1,
                "feature_shape": list(features.shape),
                "dtype": str(features.dtype),
                "is_array": True
            })

            # Additional stats for numeric arrays
            if features.ndim > 1:
                info["n_features"] = features.shape[1]

                # Check if binary
                if features.dtype in [np.int32, np.int64, np.float32, np.float64]:
                    unique_vals = np.unique(features)
                    info["is_binary"] = len(unique_vals) <= 2 and set(unique_vals).issubset({0, 1, 0.0, 1.0})

                    if info["is_binary"]:
                        info["sparsity"] = float(1.0 - (np.count_nonzero(features) / features.size))

        elif isinstance(features, list):
            info.update({
                "n_molecules": len(features),
                "is_array": False
            })

        return info

    def _convert_features_to_serializable(
        self,
        features: Any,
        single_molecule: bool
    ) -> Any:
        """
        Convert MolFeat features to JSON-serializable format

        Args:
            features: Raw MolFeat features
            single_molecule: Whether this is a single molecule

        Returns:
            JSON-serializable representation of features
        """
        # For numpy arrays
        if isinstance(features, np.ndarray):
            if single_molecule:
                return features[0].tolist() if features.ndim > 1 else features.tolist()
            else:
                return features.tolist()

        # For lists
        if isinstance(features, list):
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
        Execute MolFeat featurization

        Args:
            input_data: SMILES string or list of SMILES strings
            **kwargs: Additional parameters:
                - featurizer_type: Type of featurizer (default: "morgan")
                  Options: "morgan", "ecfp", "fcfp", "maccs", "rdkit", "atompair",
                  "topological", "avalon", "erg", "estate", "desc2d", "desc3d",
                  "usr", "usrcat", "electroshape", "pretrained"
                - pretrained_model: Name of pre-trained model (when featurizer_type="pretrained")
                  Examples: "ChemBERTa-77M-MLM", "gin_supervised_contextpred", "MolT5"
                - radius: Fingerprint radius for circular fingerprints (default: 2)
                - nBits: Fingerprint size in bits (default: 2048)
                - return_raw: If True, return raw MolFeat features (default: False)

        Returns:
            AdapterResult containing molecular features and metadata
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        if not MOLFEAT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="MolFeat is not installed. Install with: pip install molfeat"
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

        # Get pre-trained model name if using pretrained
        pretrained_model = kwargs.get("pretrained_model")

        # Validate featurizer type
        if featurizer_type not in self.FEATURIZER_TYPES and featurizer_type != "pretrained":
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported featurizer type: {featurizer_type}. "
                      f"Supported types: {list(self.FEATURIZER_TYPES.keys())} or 'pretrained'"
            )

        # Extract featurizer parameters from kwargs
        featurizer_params = {}

        # Parameters for fingerprint-based featurizers
        if featurizer_type in ["morgan", "ecfp", "fcfp", "rdkit", "atompair", "topological", "avalon"]:
            if "radius" in kwargs:
                featurizer_params["radius"] = kwargs["radius"]
            if "nBits" in kwargs:
                featurizer_params["nBits"] = kwargs["nBits"]

        # Generate features
        try:
            features = self._featurize_molecules(
                smiles_list,
                featurizer_type,
                pretrained_model=pretrained_model,
                **featurizer_params
            )
        except Exception as e:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to generate features: {str(e)}",
                metadata={
                    "source": "molfeat",
                    "featurizer_type": featurizer_type,
                    "pretrained_model": pretrained_model
                }
            )

        if features is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to generate features with MolFeat",
                metadata={
                    "source": "molfeat",
                    "smiles": smiles_list[:5],  # First 5 for debugging
                    "featurizer_type": featurizer_type,
                    "pretrained_model": pretrained_model
                }
            )

        # Get feature metadata
        feature_info = self._get_feature_info(
            featurizer_type,
            features,
            pretrained_model=pretrained_model
        )

        # Prepare result data
        return_raw = kwargs.get("return_raw", False)

        if return_raw:
            # Return raw MolFeat features (may not be JSON-serializable)
            result_data = {
                "features": features,
                "smiles": smiles_list if not single_molecule else smiles_list[0],
                **feature_info
            }
        else:
            # Return JSON-serializable format
            serializable_features = self._convert_features_to_serializable(
                features,
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
                "source": "molfeat",
                "adapter_version": self.version,
                "computation_type": "local",
                "featurizer_type": featurizer_type,
                "pretrained_model": pretrained_model,
                "n_molecules": len(smiles_list),
                "parameters": featurizer_params
            }
        )
