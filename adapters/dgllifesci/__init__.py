"""
DGL-LifeSci Adapter for PharmForge
Provides graph neural network-based molecular property prediction using DGL-LifeSci

Features:
- SMILES to DGL graph conversion with multiple featurizers
- Molecular graph featurization (canonical, AttentiveFP, Weave, minimal)
- Support for pre-trained GNN models (GCN, GAT, AttentiveFP, GIN, MPNN)
- Property prediction capabilities

Installation:
    pip install dgllife

Example Usage:
    from adapters.dgllifesci.adapter import DGLLifeSciAdapter

    # Initialize adapter
    adapter = DGLLifeSciAdapter()

    # Generate molecular graph with canonical featurization
    result = await adapter.execute(
        "CCO",
        featurizer_type="canonical"
    )

    # Generate graph with AttentiveFP featurization
    result = await adapter.execute(
        "CCO",
        featurizer_type="attentivefp"
    )

    # Generate graphs for multiple molecules
    result = await adapter.execute(
        ["CCO", "c1ccccc1", "CC(=O)O"],
        featurizer_type="canonical"
    )

    # Load pre-trained model and make predictions
    adapter.load_model("path/to/model.pt", model_type="gcn")
    result = await adapter.execute(
        "CCO",
        featurizer_type="canonical",
        include_predictions=True
    )

Supported Featurizers:
- canonical: CanonicalAtomFeaturizer + CanonicalBondFeaturizer (74 atom features, 12 bond features)
- attentivefp: AttentiveFPAtomFeaturizer + AttentiveFPBondFeaturizer (optimized for AttentiveFP models)
- weave: WeaveAtomFeaturizer + WeaveBondFeaturizer (for Weave convolutional networks)
- minimal: BaseAtomFeaturizer + BaseBondFeaturizer (minimal baseline features)

Supported Models:
- GCN: Graph Convolutional Network
- GAT: Graph Attention Network
- AttentiveFP: Attentive Fingerprint GNN
- GIN: Graph Isomorphism Network
- MPNN: Message Passing Neural Network

Note: Pre-trained models require model checkpoints to be loaded separately.
"""

from .adapter import DGLLifeSciAdapter

__all__ = ["DGLLifeSciAdapter"]
