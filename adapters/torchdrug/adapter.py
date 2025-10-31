"""
TorchDrug Adapter - GNN-based molecular modeling and feature extraction
Provides molecular embeddings and features using pre-trained GNN models from TorchDrug
"""
from typing import Any, Dict, Optional, List
import logging
import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    logging.warning("PyTorch not available - required for TorchDrug")

try:
    import torchdrug
    from torchdrug import data, models, tasks
    from torchdrug.core import Registry
    TORCHDRUG_AVAILABLE = True
except ImportError:
    TORCHDRUG_AVAILABLE = False
    logging.warning("TorchDrug not available - install with: pip install torchdrug")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class TorchDrugAdapter(AdapterProtocol):
    """
    Adapter for TorchDrug molecular feature extraction and GNN-based embeddings

    Supports multiple GNN architectures:
    - GCN (Graph Convolutional Network)
    - GAT (Graph Attention Network)
    - GIN (Graph Isomorphism Network)
    - SchNet (Continuous-filter convolutional neural network)
    - RGCN (Relational Graph Convolutional Network)
    - GraphSAGE (Graph Sample and Aggregate)
    - NFP (Neural Fingerprint)
    - MPNN (Message Passing Neural Network)

    Features:
    - Molecular graph generation from SMILES
    - GNN-based feature extraction
    - Pre-trained model loading
    - Multiple architecture support
    - Embedding generation for molecular property prediction
    """

    # Supported GNN architectures and their default configurations
    SUPPORTED_MODELS = {
        "gcn": {
            "class_name": "GCN",
            "default_params": {
                "input_dim": 69,  # Default atom feature dimension
                "hidden_dims": [256, 256, 256],
                "edge_input_dim": 4,  # Default bond feature dimension
                "batch_norm": True,
                "activation": "relu"
            },
            "description": "Graph Convolutional Network - learns node representations via neighborhood aggregation"
        },
        "gat": {
            "class_name": "GAT",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "edge_input_dim": 4,
                "num_head": 4,
                "batch_norm": True,
                "activation": "relu"
            },
            "description": "Graph Attention Network - uses attention mechanisms for neighbor aggregation"
        },
        "gin": {
            "class_name": "GIN",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "edge_input_dim": 4,
                "batch_norm": True,
                "activation": "relu",
                "eps": 0
            },
            "description": "Graph Isomorphism Network - powerful GNN for graph classification"
        },
        "schnet": {
            "class_name": "SchNet",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "edge_input_dim": 4,
                "cutoff": 5.0,
                "num_gaussian": 100,
                "batch_norm": True
            },
            "description": "SchNet - continuous-filter convolutional neural network for molecular modeling"
        },
        "rgcn": {
            "class_name": "RelationalGraphConv",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "num_relation": 4,
                "batch_norm": True,
                "activation": "relu"
            },
            "description": "Relational GCN - handles multiple edge types in molecular graphs"
        },
        "graphsage": {
            "class_name": "GraphSAGE",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "batch_norm": True,
                "activation": "relu"
            },
            "description": "GraphSAGE - inductive learning via sampling and aggregating"
        },
        "nfp": {
            "class_name": "NeuralFingerprint",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "edge_input_dim": 4,
                "batch_norm": True,
                "activation": "relu"
            },
            "description": "Neural Fingerprint - learnable molecular fingerprints"
        },
        "mpnn": {
            "class_name": "MPNN",
            "default_params": {
                "input_dim": 69,
                "hidden_dims": [256, 256, 256],
                "edge_input_dim": 4,
                "batch_norm": True,
                "activation": "relu"
            },
            "description": "Message Passing Neural Network - general framework for graph neural networks"
        }
    }

    def __init__(self):
        super().__init__(
            name="torchdrug",
            adapter_type="local",
            config={
                "timeout": 60,  # GNN inference can take time
                "default_model": "gin",  # Default to GIN (powerful and widely used)
                "default_hidden_dims": [256, 256, 256],
                "embedding_dim": 256,  # Output embedding dimension
                "batch_size": 32,
                "device": "cuda" if TORCH_AVAILABLE and torch.cuda.is_available() else "cpu"
            }
        )
        self.version = "1.0.0"
        self._model = None
        self._model_type = None

        # Check availability
        if not TORCH_AVAILABLE:
            logger.error("PyTorch is not installed! Required for TorchDrug. Install with: pip install torch")
        if not TORCHDRUG_AVAILABLE:
            logger.error("TorchDrug is not installed! Install with: pip install torchdrug")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string or list of SMILES

        Args:
            input_data: SMILES string or list of SMILES strings to validate

        Returns:
            True if valid, False otherwise
        """
        if not TORCHDRUG_AVAILABLE:
            return False

        # Handle single SMILES string
        if isinstance(input_data, str):
            if len(input_data) == 0:
                return False
            # Basic SMILES validation
            return self._validate_smiles(input_data)

        # Handle list of SMILES strings
        if isinstance(input_data, list):
            if len(input_data) == 0:
                return False
            # Validate all SMILES in the list
            return all(isinstance(s, str) and len(s) > 0 and self._validate_smiles(s) for s in input_data)

        return False

    def _validate_smiles(self, smiles: str) -> bool:
        """
        Validate a single SMILES string using TorchDrug

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        try:
            # Try to create a molecule from SMILES
            mol = data.Molecule.from_smiles(smiles, atom_feature="default", bond_feature="default")
            return mol is not None
        except Exception:
            return False

    def _get_model(self, model_type: str, **model_params) -> Optional[Any]:
        """
        Get or create a TorchDrug model instance

        Args:
            model_type: Type of GNN model ("gcn", "gat", "gin", etc.)
            **model_params: Model-specific parameters

        Returns:
            TorchDrug model instance

        Raises:
            ValueError: If model type is not supported
        """
        if not TORCHDRUG_AVAILABLE:
            raise Exception("TorchDrug is not installed")

        if model_type not in self.SUPPORTED_MODELS:
            raise ValueError(
                f"Unsupported model type: {model_type}. "
                f"Supported types: {list(self.SUPPORTED_MODELS.keys())}"
            )

        # Check if we can reuse the existing model
        if self._model is not None and self._model_type == model_type:
            return self._model

        # Get model configuration
        model_config = self.SUPPORTED_MODELS[model_type]
        class_name = model_config["class_name"]
        default_params = model_config["default_params"].copy()

        # Merge with user-provided parameters
        params = {**default_params, **model_params}

        # Filter out None values
        params = {k: v for k, v in params.items() if v is not None}

        try:
            # Get model class from TorchDrug registry
            if hasattr(models, class_name):
                model_class = getattr(models, class_name)
            else:
                raise ValueError(f"Model class {class_name} not found in TorchDrug")

            # Create model instance
            model = model_class(**params)

            # Move to appropriate device
            device = self.config.get("device", "cpu")
            model = model.to(device)

            # Set to evaluation mode
            model.eval()

            # Cache the model
            self._model = model
            self._model_type = model_type

            logger.info(f"Created TorchDrug {model_type.upper()} model with params: {params}")

            return model

        except Exception as e:
            logger.error(f"TorchDrug: Error creating model {model_type}: {e}")
            raise

    def load_pretrained_model(self, model_path: str, model_type: str) -> bool:
        """
        Load a pre-trained TorchDrug model from file

        Args:
            model_path: Path to the pre-trained model checkpoint (.pth file)
            model_type: Type of model architecture

        Returns:
            True if model loaded successfully, False otherwise
        """
        if not TORCH_AVAILABLE or not TORCHDRUG_AVAILABLE:
            logger.error("PyTorch and TorchDrug are required")
            return False

        try:
            # Create model architecture
            model = self._get_model(model_type)

            # Load checkpoint
            checkpoint = torch.load(model_path, map_location=self.config.get("device", "cpu"))

            # Load state dict
            if isinstance(checkpoint, dict) and "model_state_dict" in checkpoint:
                model.load_state_dict(checkpoint["model_state_dict"])
            else:
                model.load_state_dict(checkpoint)

            # Set to evaluation mode
            model.eval()

            # Update cached model
            self._model = model
            self._model_type = model_type

            logger.info(f"TorchDrug: Loaded pre-trained {model_type} model from {model_path}")
            return True

        except Exception as e:
            logger.error(f"TorchDrug: Error loading model from {model_path}: {e}")
            return False

    def _create_molecular_graph(self, smiles: str) -> Optional[Any]:
        """
        Create a TorchDrug molecular graph from SMILES

        Args:
            smiles: SMILES string

        Returns:
            TorchDrug Molecule object or None if failed
        """
        if not TORCHDRUG_AVAILABLE:
            raise Exception("TorchDrug is not installed")

        try:
            # Create molecule with default features
            mol = data.Molecule.from_smiles(
                smiles,
                atom_feature="default",
                bond_feature="default",
                mol_feature=None
            )

            if mol is None or mol.num_node == 0:
                logger.warning(f"TorchDrug: Could not create valid molecule from SMILES: {smiles}")
                return None

            return mol

        except Exception as e:
            logger.error(f"TorchDrug: Error creating molecular graph for {smiles}: {e}")
            return None

    def _extract_graph_features(self, mol: Any) -> Dict[str, Any]:
        """
        Extract basic graph features from a TorchDrug molecule

        Args:
            mol: TorchDrug Molecule object

        Returns:
            Dictionary of graph features
        """
        features = {
            "num_atoms": int(mol.num_node),
            "num_bonds": int(mol.num_edge),
            "num_residues": int(mol.num_residue) if hasattr(mol, "num_residue") else 0,
            "atom_feature_dim": int(mol.atom_feature.shape[-1]) if hasattr(mol, "atom_feature") else 0,
            "bond_feature_dim": int(mol.bond_feature.shape[-1]) if hasattr(mol, "bond_feature") else 0,
        }

        # Add molecular weight if available
        if hasattr(mol, "molecular_weight"):
            features["molecular_weight"] = float(mol.molecular_weight)

        # Add atom types
        if hasattr(mol, "atom_type"):
            atom_types = mol.atom_type.cpu().numpy() if torch.is_tensor(mol.atom_type) else mol.atom_type
            features["unique_atom_types"] = int(len(np.unique(atom_types)))

        # Add bond types
        if hasattr(mol, "bond_type"):
            bond_types = mol.bond_type.cpu().numpy() if torch.is_tensor(mol.bond_type) else mol.bond_type
            features["unique_bond_types"] = int(len(np.unique(bond_types)))

        return features

    def _generate_embeddings(
        self,
        molecules: List[Any],
        model: Any,
        batch_size: int = 32
    ) -> Optional[np.ndarray]:
        """
        Generate molecular embeddings using a GNN model

        Args:
            molecules: List of TorchDrug Molecule objects
            model: TorchDrug GNN model
            batch_size: Batch size for inference

        Returns:
            NumPy array of embeddings [n_molecules, embedding_dim]
        """
        if not TORCH_AVAILABLE or not TORCHDRUG_AVAILABLE:
            raise Exception("PyTorch and TorchDrug are required")

        try:
            device = self.config.get("device", "cpu")
            embeddings_list = []

            # Process in batches
            for i in range(0, len(molecules), batch_size):
                batch = molecules[i:i + batch_size]

                # Create batch
                batch_data = data.Molecule.pack(batch)
                batch_data = batch_data.to(device)

                # Generate embeddings
                with torch.no_grad():
                    # Forward pass through model
                    output = model(batch_data, batch_data.node_feature.float())

                    # Extract graph-level representation
                    # Most TorchDrug models return node features, so we pool them
                    if hasattr(output, "shape") and len(output.shape) > 1:
                        # Pool node features to graph features (mean pooling)
                        graph_feature = output["graph_feature"] if isinstance(output, dict) else output

                        # If we have node-level features, aggregate them
                        if len(graph_feature.shape) == 2:
                            # Already graph-level features
                            batch_embeddings = graph_feature
                        else:
                            # Need to pool node features
                            # Use scatter_mean to pool per graph
                            batch_embeddings = []
                            for mol in batch:
                                mol_mask = (batch_data.node2graph == batch.index(mol))
                                mol_features = graph_feature[mol_mask]
                                mol_embedding = mol_features.mean(dim=0, keepdim=True)
                                batch_embeddings.append(mol_embedding)
                            batch_embeddings = torch.cat(batch_embeddings, dim=0)
                    else:
                        batch_embeddings = output

                    # Convert to numpy
                    batch_embeddings = batch_embeddings.cpu().numpy()
                    embeddings_list.append(batch_embeddings)

            # Concatenate all batches
            embeddings = np.vstack(embeddings_list)

            return embeddings

        except Exception as e:
            logger.error(f"TorchDrug: Error generating embeddings: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute TorchDrug molecular feature extraction and embedding generation

        Args:
            input_data: SMILES string or list of SMILES strings
            **kwargs: Additional parameters:
                - model_type: GNN architecture to use ("gcn", "gat", "gin", etc.)
                - generate_embeddings: Whether to generate GNN embeddings (default: True)
                - pretrained_model_path: Path to pre-trained model checkpoint
                - hidden_dims: List of hidden layer dimensions
                - batch_size: Batch size for embedding generation
                - return_graph_objects: Return raw TorchDrug molecule objects (not serializable)
                - include_graph_features: Include basic graph statistics (default: True)

        Returns:
            AdapterResult containing molecular embeddings and/or features
        """
        # Check dependencies
        if not TORCH_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="PyTorch is not installed. Install with: pip install torch"
            )

        if not TORCHDRUG_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="TorchDrug is not installed. Install with: pip install torchdrug"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string(s) or could not parse with TorchDrug"
            )

        # Normalize input to list
        if isinstance(input_data, str):
            smiles_list = [input_data]
            single_molecule = True
        else:
            smiles_list = input_data
            single_molecule = False

        # Extract parameters
        model_type = kwargs.get("model_type", self.config.get("default_model", "gin"))
        generate_embeddings = kwargs.get("generate_embeddings", True)
        pretrained_model_path = kwargs.get("pretrained_model_path")
        include_graph_features = kwargs.get("include_graph_features", True)
        return_graph_objects = kwargs.get("return_graph_objects", False)
        batch_size = kwargs.get("batch_size", self.config.get("batch_size", 32))

        # Validate model type
        if model_type not in self.SUPPORTED_MODELS:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported model type: {model_type}. "
                      f"Supported types: {list(self.SUPPORTED_MODELS.keys())}"
            )

        # Create molecular graphs
        molecules = []
        failed_smiles = []

        for smiles in smiles_list:
            mol = self._create_molecular_graph(smiles)
            if mol is not None:
                molecules.append(mol)
            else:
                failed_smiles.append(smiles)

        if len(molecules) == 0:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to create molecular graphs for all provided SMILES",
                metadata={
                    "source": "torchdrug",
                    "failed_smiles": failed_smiles
                }
            )

        # Prepare result data
        result_data = {
            "smiles": smiles_list if not single_molecule else smiles_list[0],
            "n_molecules": len(molecules),
            "n_failed": len(failed_smiles)
        }

        if len(failed_smiles) > 0:
            result_data["failed_smiles"] = failed_smiles

        # Extract basic graph features if requested
        if include_graph_features:
            graph_features = [self._extract_graph_features(mol) for mol in molecules]

            if single_molecule:
                result_data["graph_features"] = graph_features[0]
            else:
                result_data["graph_features"] = graph_features

        # Generate embeddings if requested
        if generate_embeddings:
            try:
                # Get model parameters from kwargs
                model_params = {}
                if "hidden_dims" in kwargs:
                    model_params["hidden_dims"] = kwargs["hidden_dims"]

                # Load pre-trained model if provided
                if pretrained_model_path:
                    success = self.load_pretrained_model(pretrained_model_path, model_type)
                    if not success:
                        return AdapterResult(
                            success=False,
                            data=None,
                            error=f"Failed to load pre-trained model from {pretrained_model_path}"
                        )
                    model = self._model
                else:
                    # Create new model
                    model = self._get_model(model_type, **model_params)

                # Generate embeddings
                embeddings = self._generate_embeddings(molecules, model, batch_size)

                if embeddings is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="Failed to generate embeddings with TorchDrug",
                        metadata={
                            "source": "torchdrug",
                            "model_type": model_type,
                            "smiles": smiles_list[:5]  # First 5 for debugging
                        }
                    )

                # Add embeddings to result
                if single_molecule:
                    result_data["embeddings"] = embeddings[0].tolist()
                    result_data["embedding_dim"] = int(embeddings.shape[1])
                else:
                    result_data["embeddings"] = embeddings.tolist()
                    result_data["embedding_shape"] = list(embeddings.shape)
                    result_data["embedding_dim"] = int(embeddings.shape[1])

                # Add embedding statistics
                result_data["embedding_stats"] = {
                    "mean": float(np.mean(embeddings)),
                    "std": float(np.std(embeddings)),
                    "min": float(np.min(embeddings)),
                    "max": float(np.max(embeddings))
                }

            except Exception as e:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Error during embedding generation: {str(e)}",
                    metadata={
                        "source": "torchdrug",
                        "model_type": model_type
                    }
                )

        # Optionally include raw graph objects
        if return_graph_objects:
            result_data["molecule_objects"] = molecules if not single_molecule else molecules[0]

        # Build metadata
        metadata = {
            "source": "torchdrug",
            "adapter_version": self.version,
            "computation_type": "local",
            "model_type": model_type,
            "model_description": self.SUPPORTED_MODELS[model_type]["description"],
            "n_molecules": len(molecules),
            "n_failed": len(failed_smiles),
            "embeddings_generated": generate_embeddings,
            "device": self.config.get("device", "cpu"),
            "torchdrug_version": torchdrug.__version__ if TORCHDRUG_AVAILABLE else None,
            "torch_version": torch.__version__ if TORCH_AVAILABLE else None
        }

        if pretrained_model_path:
            metadata["pretrained_model_path"] = pretrained_model_path

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata=metadata
        )
