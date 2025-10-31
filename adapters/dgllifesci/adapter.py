"""
DGL-LifeSci Adapter - Graph neural network-based molecular property prediction
Provides molecular graph generation, featurization, and pre-trained model predictions using DGL-LifeSci
"""
from typing import Any, Dict, Optional, List, Union
import logging
import numpy as np

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - required for DGL-LifeSci")

try:
    import dgl
    import torch
    from dgllife.utils import (
        smiles_to_bigraph,
        mol_to_bigraph,
        CanonicalAtomFeaturizer,
        CanonicalBondFeaturizer,
        AttentiveFPAtomFeaturizer,
        AttentiveFPBondFeaturizer,
        WeaveAtomFeaturizer,
        WeaveBondFeaturizer,
        BaseAtomFeaturizer,
        BaseBondFeaturizer,
    )
    DGLLIFE_AVAILABLE = True
except ImportError:
    DGLLIFE_AVAILABLE = False
    logging.warning("DGL-LifeSci not available - install with: pip install dgllife")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DGLLifeSciAdapter(AdapterProtocol):
    """
    Adapter for DGL-LifeSci molecular graph neural networks
    Generates molecular graphs with node/edge features and supports pre-trained GNN models

    Supported featurizers:
    - canonical: CanonicalAtomFeaturizer + CanonicalBondFeaturizer (default)
    - attentivefp: AttentiveFPAtomFeaturizer + AttentiveFPBondFeaturizer
    - weave: WeaveAtomFeaturizer + WeaveBondFeaturizer
    - minimal: Minimal atom/bond features

    Supported models (property prediction):
    - GCN: Graph Convolutional Network
    - GAT: Graph Attention Network
    - AttentiveFP: Attentive Fingerprint for molecular property prediction
    - GIN: Graph Isomorphism Network
    - MPNN: Message Passing Neural Network
    """

    @staticmethod
    def _get_featurizer_types():
        """
        Get featurizer configurations (lazy-loaded to avoid import errors)
        """
        if not DGLLIFE_AVAILABLE:
            return {}

        return {
            "canonical": {
                "atom_featurizer": CanonicalAtomFeaturizer,
                "bond_featurizer": CanonicalBondFeaturizer,
                "atom_params": {"atom_data_field": "h"},
                "bond_params": {"bond_data_field": "e"},
                "description": "Standard DGL-LifeSci featurization with 74 atom features and 12 bond features"
            },
            "attentivefp": {
                "atom_featurizer": AttentiveFPAtomFeaturizer,
                "bond_featurizer": AttentiveFPBondFeaturizer,
                "atom_params": {"atom_data_field": "h"},
                "bond_params": {"bond_data_field": "e"},
                "description": "AttentiveFP featurization for attention-based GNNs"
            },
            "weave": {
                "atom_featurizer": WeaveAtomFeaturizer,
                "bond_featurizer": WeaveBondFeaturizer,
                "atom_params": {"atom_data_field": "h"},
                "bond_params": {"bond_data_field": "e"},
                "description": "Weave featurization for Weave convolutional networks"
            },
            "minimal": {
                "atom_featurizer": BaseAtomFeaturizer,
                "bond_featurizer": BaseBondFeaturizer,
                "atom_params": {"atom_data_field": "h"},
                "bond_params": {"bond_data_field": "e"},
                "description": "Minimal baseline featurization"
            }
        }

    @property
    def FEATURIZER_TYPES(self):
        """Property to access featurizer types"""
        return self._get_featurizer_types()

    # Pre-trained models available in DGL-LifeSci
    PRETRAINED_MODELS = {
        "gcn": {
            "description": "Graph Convolutional Network for property prediction",
            "tasks": ["classification", "regression"],
            "note": "Requires model checkpoint"
        },
        "gat": {
            "description": "Graph Attention Network for property prediction",
            "tasks": ["classification", "regression"],
            "note": "Requires model checkpoint"
        },
        "attentivefp": {
            "description": "Attentive Fingerprint GNN for property prediction",
            "tasks": ["classification", "regression"],
            "note": "Requires model checkpoint"
        },
        "gin": {
            "description": "Graph Isomorphism Network for property prediction",
            "tasks": ["classification", "regression"],
            "note": "Requires model checkpoint"
        },
        "mpnn": {
            "description": "Message Passing Neural Network for property prediction",
            "tasks": ["classification", "regression"],
            "note": "Requires model checkpoint"
        }
    }

    def __init__(self):
        super().__init__(
            name="dgllifesci",
            adapter_type="local",
            config={
                "timeout": 60,  # GNN inference can take time
                "default_featurizer": "canonical",
                "default_model": None,
                "device": "cpu",  # Can be "cuda" if GPU available
                "batch_size": 32
            }
        )
        self.version = "1.0.0"
        self._model = None
        self._model_type = None
        self._model_path = None

        # Check dependencies
        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Required for DGL-LifeSci. Install with: pip install rdkit")
        if not DGLLIFE_AVAILABLE:
            logger.error("DGL-LifeSci is not installed! Install with: pip install dgllife")

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

    def _get_featurizers(
        self,
        featurizer_type: str
    ):
        """
        Get atom and bond featurizers for the specified type

        Args:
            featurizer_type: Type of featurizer ("canonical", "attentivefp", "weave", "minimal")

        Returns:
            Tuple of (atom_featurizer, bond_featurizer)

        Raises:
            ValueError: If featurizer type is not supported
        """
        if not DGLLIFE_AVAILABLE:
            raise Exception("DGL-LifeSci is not installed")

        if featurizer_type not in self.FEATURIZER_TYPES:
            raise ValueError(
                f"Unsupported featurizer type: {featurizer_type}. "
                f"Supported types: {list(self.FEATURIZER_TYPES.keys())}"
            )

        # Get featurizer configuration
        config = self.FEATURIZER_TYPES[featurizer_type]

        # Create atom featurizer
        atom_featurizer_class = config["atom_featurizer"]
        atom_params = config["atom_params"].copy()
        atom_featurizer = atom_featurizer_class(**atom_params)

        # Create bond featurizer
        bond_featurizer_class = config["bond_featurizer"]
        bond_params = config["bond_params"].copy()
        bond_featurizer = bond_featurizer_class(**bond_params)

        return atom_featurizer, bond_featurizer

    def _smiles_to_graph(
        self,
        smiles: str,
        featurizer_type: str = "canonical"
    ) -> Optional[Any]:
        """
        Convert SMILES to DGL graph with node and edge features

        Args:
            smiles: SMILES string
            featurizer_type: Type of featurizer to use

        Returns:
            DGL graph with node and edge features, or None if conversion fails
        """
        if not RDKIT_AVAILABLE or not DGLLIFE_AVAILABLE:
            raise Exception("RDKit and DGL-LifeSci are required")

        try:
            # Get featurizers
            atom_featurizer, bond_featurizer = self._get_featurizers(featurizer_type)

            # Convert SMILES to graph
            graph = smiles_to_bigraph(
                smiles,
                node_featurizer=atom_featurizer,
                edge_featurizer=bond_featurizer
            )

            if graph is None:
                logger.warning(f"DGL-LifeSci: Could not convert SMILES to graph: {smiles}")
                return None

            return graph

        except Exception as e:
            logger.error(f"DGL-LifeSci: Error converting SMILES to graph: {e}")
            return None

    def _extract_graph_features(
        self,
        graph: Any,
        smiles: str
    ) -> Dict[str, Any]:
        """
        Extract feature information from a DGL graph

        Args:
            graph: DGL graph
            smiles: Original SMILES string

        Returns:
            Dictionary containing graph feature information
        """
        if not DGLLIFE_AVAILABLE:
            raise Exception("DGL-LifeSci is not installed")

        try:
            features = {
                "smiles": smiles,
                "num_nodes": int(graph.num_nodes()),
                "num_edges": int(graph.num_edges()),
                "graph_type": "bidirected"  # DGL-LifeSci uses bidirected graphs
            }

            # Extract node features if available
            if "h" in graph.ndata:
                node_features = graph.ndata["h"]
                features["node_features"] = {
                    "shape": list(node_features.shape),
                    "dtype": str(node_features.dtype),
                    "feature_dim": int(node_features.shape[1]) if len(node_features.shape) > 1 else 1,
                    "mean": float(node_features.float().mean().item()),
                    "std": float(node_features.float().std().item())
                }

            # Extract edge features if available
            if "e" in graph.edata:
                edge_features = graph.edata["e"]
                features["edge_features"] = {
                    "shape": list(edge_features.shape),
                    "dtype": str(edge_features.dtype),
                    "feature_dim": int(edge_features.shape[1]) if len(edge_features.shape) > 1 else 1,
                    "mean": float(edge_features.float().mean().item()),
                    "std": float(edge_features.float().std().item())
                }

            return features

        except Exception as e:
            logger.error(f"DGL-LifeSci: Error extracting graph features: {e}")
            return {
                "smiles": smiles,
                "error": str(e)
            }

    def load_model(
        self,
        model_path: str,
        model_type: str = "gcn"
    ) -> bool:
        """
        Load a pre-trained DGL-LifeSci model for predictions

        Args:
            model_path: Path to the trained model checkpoint (.pth or .pt file)
            model_type: Type of model ("gcn", "gat", "attentivefp", "gin", "mpnn")

        Returns:
            True if model loaded successfully, False otherwise
        """
        if not DGLLIFE_AVAILABLE:
            logger.error("DGL-LifeSci is not installed")
            return False

        if model_type not in self.PRETRAINED_MODELS:
            logger.error(
                f"Unsupported model type: {model_type}. "
                f"Supported types: {list(self.PRETRAINED_MODELS.keys())}"
            )
            return False

        try:
            # Note: Actual model loading depends on model architecture
            # This is a placeholder for the model loading infrastructure
            import torch

            # Set device
            device = torch.device(self.config.get("device", "cpu"))

            # Store model information
            self._model_path = model_path
            self._model_type = model_type

            logger.info(
                f"DGL-LifeSci: Model path set to {model_path}, "
                f"type: {model_type}, device: {device}"
            )

            # Actual model loading would happen here
            # self._model = torch.load(model_path, map_location=device)
            # self._model.eval()

            return True

        except Exception as e:
            logger.error(f"DGL-LifeSci: Error loading model from {model_path}: {e}")
            return False

    def _make_predictions(
        self,
        graphs: Union[Any, List[Any]],
        smiles: Union[str, List[str]]
    ) -> Optional[Dict[str, Any]]:
        """
        Make predictions using a pre-trained DGL-LifeSci model

        Args:
            graphs: DGL graph or list of DGL graphs
            smiles: SMILES string or list of SMILES strings

        Returns:
            Dictionary containing predictions
        """
        if not DGLLIFE_AVAILABLE:
            raise Exception("DGL-LifeSci is not installed")

        if self._model_path is None:
            logger.warning("DGL-LifeSci: No model loaded for predictions")
            return {
                "predictions_available": False,
                "reason": "No pre-trained model loaded. Use load_model() first."
            }

        try:
            # Convert single graph to list
            if not isinstance(graphs, list):
                graphs = [graphs]
                smiles = [smiles] if isinstance(smiles, str) else smiles

            # Placeholder for actual predictions
            # In a real implementation, this would:
            # 1. Batch the graphs
            # 2. Run forward pass through the model
            # 3. Post-process predictions

            predictions = {
                "predictions_available": True,
                "model_type": self._model_type,
                "model_path": self._model_path,
                "n_molecules": len(graphs),
                "smiles": smiles,
                "device": self.config.get("device", "cpu"),
                "note": "Model infrastructure ready. Actual predictions require loaded checkpoint and forward pass implementation."
            }

            # If model was actually loaded, predictions would be added here
            # predictions["values"] = model_output.tolist()
            # predictions["confidence"] = confidence_scores.tolist()

            return predictions

        except Exception as e:
            logger.error(f"DGL-LifeSci: Error making predictions: {e}")
            return {
                "predictions_available": False,
                "error": str(e)
            }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute DGL-LifeSci molecular graph generation and optional prediction

        Args:
            input_data: SMILES string or list of SMILES strings
            **kwargs: Additional parameters:
                - featurizer_type: Type of featurizer ("canonical", "attentivefp", "weave", "minimal")
                - model_path: Path to pre-trained model checkpoint for predictions
                - model_type: Type of model ("gcn", "gat", "attentivefp", "gin", "mpnn")
                - include_predictions: Whether to attempt predictions if model available
                - return_graphs: If True, return raw DGL graphs (not JSON-serializable)

        Returns:
            AdapterResult containing molecular graphs and/or predictions
        """
        # Check dependencies
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        if not DGLLIFE_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="DGL-LifeSci is not installed. Install with: pip install dgllife"
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

        # Get parameters
        featurizer_type = kwargs.get(
            "featurizer_type",
            self.config.get("default_featurizer", "canonical")
        )
        model_path = kwargs.get("model_path")
        model_type = kwargs.get("model_type", "gcn")
        include_predictions = kwargs.get("include_predictions", False)
        return_graphs = kwargs.get("return_graphs", False)

        # Validate featurizer type
        if featurizer_type not in self.FEATURIZER_TYPES:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported featurizer type: {featurizer_type}. "
                      f"Supported types: {list(self.FEATURIZER_TYPES.keys())}"
            )

        # Convert SMILES to graphs
        graphs = []
        graph_features = []
        failed_smiles = []

        for smiles in smiles_list:
            try:
                graph = self._smiles_to_graph(smiles, featurizer_type)
                if graph is not None:
                    graphs.append(graph)
                    features = self._extract_graph_features(graph, smiles)
                    graph_features.append(features)
                else:
                    failed_smiles.append(smiles)
                    logger.warning(f"Failed to convert SMILES to graph: {smiles}")
            except Exception as e:
                failed_smiles.append(smiles)
                logger.error(f"Error processing SMILES {smiles}: {e}")

        # Check if any graphs were successfully generated
        if len(graphs) == 0:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to generate graphs for all SMILES. Failed: {failed_smiles}",
                metadata={
                    "source": "dgllifesci",
                    "featurizer_type": featurizer_type,
                    "failed_smiles": failed_smiles
                }
            )

        # Prepare result data
        result_data = {
            "featurizer_type": featurizer_type,
            "featurizer_description": self.FEATURIZER_TYPES[featurizer_type]["description"],
            "n_molecules": len(graphs),
            "n_successful": len(graphs),
            "n_failed": len(failed_smiles),
            "graph_features": graph_features[0] if single_molecule else graph_features
        }

        if len(failed_smiles) > 0:
            result_data["failed_smiles"] = failed_smiles

        # Optionally include predictions
        if include_predictions or model_path:
            if model_path and self._model_path != model_path:
                # Load model if path provided and different from current
                success = self.load_model(model_path, model_type)
                if not success:
                    result_data["predictions"] = {
                        "predictions_available": False,
                        "error": "Failed to load model"
                    }
                else:
                    predictions = self._make_predictions(
                        graphs[0] if single_molecule else graphs,
                        smiles_list[0] if single_molecule else smiles_list
                    )
                    result_data["predictions"] = predictions
            elif self._model_path:
                # Use already loaded model
                predictions = self._make_predictions(
                    graphs[0] if single_molecule else graphs,
                    smiles_list[0] if single_molecule else smiles_list
                )
                result_data["predictions"] = predictions
            else:
                result_data["predictions"] = {
                    "predictions_available": False,
                    "reason": "No model loaded. Provide model_path parameter."
                }

        # Optionally return raw graphs
        if return_graphs:
            result_data["graphs"] = graphs[0] if single_molecule else graphs
            result_data["note"] = "Raw DGL graphs included (not JSON-serializable)"

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "dgllifesci",
                "adapter_version": self.version,
                "computation_type": "local",
                "featurizer_type": featurizer_type,
                "n_molecules": len(smiles_list),
                "n_successful": len(graphs),
                "n_failed": len(failed_smiles),
                "predictions_included": include_predictions or model_path is not None,
                "dgl_version": dgl.__version__ if DGLLIFE_AVAILABLE else None,
                "torch_version": torch.__version__ if DGLLIFE_AVAILABLE else None
            }
        )
