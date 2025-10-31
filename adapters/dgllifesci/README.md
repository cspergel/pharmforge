# DGL-LifeSci Adapter

Graph neural network-based molecular property prediction using DGL-LifeSci.

## Overview

The DGL-LifeSci adapter provides:
- **Molecular Graph Generation**: Convert SMILES to DGL graphs with rich node and edge features
- **Multiple Featurizers**: Support for canonical, AttentiveFP, Weave, and minimal featurizations
- **Pre-trained Models**: Integration with GNN models (GCN, GAT, AttentiveFP, GIN, MPNN)
- **Property Prediction**: Molecular property prediction using trained models
- **Batch Processing**: Efficient processing of multiple molecules

## Installation

```bash
pip install dgllife
```

**Dependencies:**
- `dgllife`: DGL-LifeSci library
- `dgl`: Deep Graph Library
- `torch`: PyTorch
- `rdkit`: RDKit for molecular processing

## Usage

### Basic Graph Generation

```python
from adapters.dgllifesci.adapter import DGLLifeSciAdapter

# Initialize adapter
adapter = DGLLifeSciAdapter()

# Generate molecular graph with canonical featurization
result = await adapter.execute("CCO", featurizer_type="canonical")

if result.success:
    data = result.data
    print(f"Nodes: {data['graph_features']['num_nodes']}")
    print(f"Edges: {data['graph_features']['num_edges']}")
    print(f"Node feature dim: {data['graph_features']['node_features']['feature_dim']}")
    print(f"Edge feature dim: {data['graph_features']['edge_features']['feature_dim']}")
```

### Multiple Featurizers

```python
# Canonical featurization (default)
result = await adapter.execute("CCO", featurizer_type="canonical")

# AttentiveFP featurization (optimized for attention-based GNNs)
result = await adapter.execute("CCO", featurizer_type="attentivefp")

# Weave featurization (for Weave convolutional networks)
result = await adapter.execute("CCO", featurizer_type="weave")

# Minimal featurization (baseline features)
result = await adapter.execute("CCO", featurizer_type="minimal")
```

### Batch Processing

```python
# Process multiple molecules at once
smiles_list = [
    "CCO",                    # Ethanol
    "c1ccccc1",              # Benzene
    "CC(=O)O",               # Acetic acid
    "CC(C)CC1=CC=C(C=C1)C(C)C"  # Ibuprofen
]

result = await adapter.execute(smiles_list, featurizer_type="canonical")

if result.success:
    data = result.data
    print(f"Successfully processed: {data['n_successful']} molecules")
    print(f"Failed: {data['n_failed']} molecules")

    # Access individual graph features
    for features in data['graph_features']:
        print(f"SMILES: {features['smiles']}")
        print(f"Nodes: {features['num_nodes']}, Edges: {features['num_edges']}")
```

### Pre-trained Model Predictions

```python
# Load a pre-trained model
adapter.load_model("path/to/gcn_model.pt", model_type="gcn")

# Make predictions
result = await adapter.execute(
    "CCO",
    featurizer_type="canonical",
    include_predictions=True
)

if result.success:
    predictions = result.data['predictions']
    print(f"Predictions available: {predictions['predictions_available']}")
    print(f"Model type: {predictions['model_type']}")
```

### Return Raw DGL Graphs

```python
# Get raw DGL graph objects (not JSON-serializable)
result = await adapter.execute(
    "CCO",
    featurizer_type="canonical",
    return_graphs=True
)

if result.success:
    graph = result.data['graphs']  # DGL graph object

    # Access DGL graph methods
    print(f"Number of nodes: {graph.num_nodes()}")
    print(f"Number of edges: {graph.num_edges()}")

    # Access node features
    node_features = graph.ndata['h']
    print(f"Node features shape: {node_features.shape}")

    # Access edge features
    edge_features = graph.edata['e']
    print(f"Edge features shape: {edge_features.shape}")
```

## Featurizer Types

### Canonical Featurizer
- **Node features**: 74 dimensions
- **Edge features**: 12 dimensions
- **Description**: Standard DGL-LifeSci featurization with comprehensive atom and bond features
- **Use case**: General-purpose molecular property prediction

```python
result = await adapter.execute("CCO", featurizer_type="canonical")
```

### AttentiveFP Featurizer
- **Node features**: Optimized for attention mechanisms
- **Edge features**: Bond features for attention-based GNNs
- **Description**: Specialized featurization for AttentiveFP models
- **Use case**: Attention-based graph neural networks

```python
result = await adapter.execute("CCO", featurizer_type="attentivefp")
```

### Weave Featurizer
- **Node features**: Designed for Weave architecture
- **Edge features**: Bond features for Weave convolutions
- **Description**: Featurization for Weave convolutional networks
- **Use case**: Weave-specific models

```python
result = await adapter.execute("CCO", featurizer_type="weave")
```

### Minimal Featurizer
- **Node features**: Minimal baseline features
- **Edge features**: Minimal baseline bond features
- **Description**: Lightweight featurization with basic atom and bond information
- **Use case**: Baseline models or when minimal features are sufficient

```python
result = await adapter.execute("CCO", featurizer_type="minimal")
```

## Supported Models

The adapter supports the following pre-trained GNN architectures:

| Model | Description | Tasks |
|-------|-------------|-------|
| **GCN** | Graph Convolutional Network | Classification, Regression |
| **GAT** | Graph Attention Network | Classification, Regression |
| **AttentiveFP** | Attentive Fingerprint GNN | Classification, Regression |
| **GIN** | Graph Isomorphism Network | Classification, Regression |
| **MPNN** | Message Passing Neural Network | Classification, Regression |

### Loading Pre-trained Models

```python
# Load GCN model
adapter.load_model("models/gcn_solubility.pt", model_type="gcn")

# Load GAT model
adapter.load_model("models/gat_toxicity.pt", model_type="gat")

# Load AttentiveFP model
adapter.load_model("models/attentivefp_bioactivity.pt", model_type="attentivefp")

# Make predictions
result = await adapter.execute(
    "CCO",
    featurizer_type="canonical",
    include_predictions=True
)
```

## Response Format

### Successful Graph Generation

```python
{
    "success": True,
    "data": {
        "featurizer_type": "canonical",
        "featurizer_description": "Standard DGL-LifeSci featurization...",
        "n_molecules": 1,
        "n_successful": 1,
        "n_failed": 0,
        "graph_features": {
            "smiles": "CCO",
            "num_nodes": 9,
            "num_edges": 16,
            "graph_type": "bidirected",
            "node_features": {
                "shape": [9, 74],
                "dtype": "torch.float32",
                "feature_dim": 74,
                "mean": 0.234,
                "std": 0.412
            },
            "edge_features": {
                "shape": [16, 12],
                "dtype": "torch.float32",
                "feature_dim": 12,
                "mean": 0.156,
                "std": 0.298
            }
        }
    },
    "metadata": {
        "source": "dgllifesci",
        "adapter_version": "1.0.0",
        "computation_type": "local",
        "featurizer_type": "canonical",
        "n_molecules": 1,
        "n_successful": 1,
        "n_failed": 0,
        "predictions_included": False,
        "dgl_version": "1.1.0",
        "torch_version": "2.0.0"
    }
}
```

### With Predictions

```python
{
    "success": True,
    "data": {
        "featurizer_type": "canonical",
        "n_molecules": 1,
        "graph_features": { ... },
        "predictions": {
            "predictions_available": True,
            "model_type": "gcn",
            "model_path": "models/gcn_model.pt",
            "n_molecules": 1,
            "smiles": ["CCO"],
            "device": "cpu",
            "note": "Model infrastructure ready..."
        }
    }
}
```

### Error Response

```python
{
    "success": False,
    "data": None,
    "error": "DGL-LifeSci is not installed. Install with: pip install dgllife",
    "metadata": {
        "source": "dgllifesci"
    }
}
```

## Configuration

The adapter can be configured through its config dictionary:

```python
adapter = DGLLifeSciAdapter()

# Update configuration
adapter.config.update({
    "timeout": 120,                    # Increase timeout for large molecules
    "default_featurizer": "attentivefp",  # Change default featurizer
    "device": "cuda",                  # Use GPU if available
    "batch_size": 64                   # Batch size for predictions
})
```

## Graph Structure

DGL-LifeSci generates **bidirected graphs** where:
- Each chemical bond is represented by two directed edges (one in each direction)
- Node features are stored in `graph.ndata['h']`
- Edge features are stored in `graph.edata['e']`

### Node Features (Canonical)
- Atom type (one-hot encoded)
- Degree
- Formal charge
- Hybridization
- Aromaticity
- Number of hydrogen atoms
- Chirality
- Mass

### Edge Features (Canonical)
- Bond type (single, double, triple, aromatic)
- Conjugation
- Ring membership
- Stereochemistry

## Advanced Usage

### Custom Featurizer Parameters

```python
# Note: This adapter uses predefined featurizers
# For custom featurization, use the featurizer classes directly from dgllife.utils
```

### GPU Acceleration

```python
# Enable GPU if available
adapter.config["device"] = "cuda"

# Load model on GPU
adapter.load_model("model.pt", model_type="gcn")

# Make predictions (automatically uses GPU)
result = await adapter.execute("CCO", include_predictions=True)
```

### Error Handling

```python
result = await adapter.execute(smiles_list, featurizer_type="canonical")

if result.success:
    data = result.data

    # Check for partial failures
    if data['n_failed'] > 0:
        print(f"Warning: {data['n_failed']} molecules failed")
        print(f"Failed SMILES: {data.get('failed_smiles', [])}")

    # Process successful graphs
    for features in data['graph_features']:
        print(f"Processing: {features['smiles']}")
else:
    print(f"Error: {result.error}")
```

## Comparison with Other Adapters

| Feature | DGL-LifeSci | Chemprop | DeepChem |
|---------|-------------|----------|----------|
| Graph Generation | ✓ Bidirected DGL | ✓ Directed | ✓ Various formats |
| Featurizers | 4 types | 1 (built-in) | 8+ types |
| Pre-trained Models | GCN, GAT, AttentiveFP, etc. | MPNN | Various |
| GPU Support | ✓ | ✓ | ✓ |
| Batch Processing | ✓ | ✓ | ✓ |

## Common Use Cases

### 1. Molecular Property Prediction

```python
# Load pre-trained solubility model
adapter.load_model("models/gcn_solubility.pt", model_type="gcn")

# Predict solubility for compounds
compounds = ["CCO", "c1ccccc1", "CC(=O)O"]
result = await adapter.execute(
    compounds,
    featurizer_type="canonical",
    include_predictions=True
)
```

### 2. Graph Feature Extraction

```python
# Extract graph features for downstream ML
result = await adapter.execute(
    smiles_list,
    featurizer_type="canonical",
    return_graphs=False  # Get JSON-serializable features
)

# Use features in other models
features = result.data['graph_features']
```

### 3. Model Development

```python
# Generate graphs for model training
result = await adapter.execute(
    training_smiles,
    featurizer_type="canonical",
    return_graphs=True  # Get raw DGL graphs
)

graphs = result.data['graphs']

# Use graphs for training custom GNN models
# (training code not shown)
```

## Performance Considerations

- **Featurization**: Canonical featurization is generally fast (<100ms per molecule)
- **GPU**: Use `device="cuda"` for faster model inference on large batches
- **Batch Size**: Process multiple molecules together for better efficiency
- **Memory**: DGL graphs can be memory-intensive for large molecules

## Troubleshooting

### Import Errors

```python
# If you see "DGL-LifeSci is not installed"
pip install dgllife

# If you see "RDKit not available"
pip install rdkit
```

### Graph Conversion Failures

```python
# Some SMILES may fail to convert to graphs
result = await adapter.execute(smiles_list, featurizer_type="canonical")

# Check for failures
if result.data['n_failed'] > 0:
    print(f"Failed SMILES: {result.data['failed_smiles']}")
```

### Model Loading Issues

```python
# Ensure model file exists and is compatible
import os
if not os.path.exists("model.pt"):
    print("Model file not found")

# Check model type matches architecture
adapter.load_model("model.pt", model_type="gcn")  # Must match trained architecture
```

## References

- [DGL-LifeSci Documentation](https://lifesci.dgl.ai/)
- [DGL Documentation](https://docs.dgl.ai/)
- [Paper: MoleculeNet (Benchmarking dataset)](https://arxiv.org/abs/1703.00564)
- [Paper: Pushing the Boundaries of Molecular Representation](https://arxiv.org/abs/2011.03230)

## License

This adapter follows the PharmForge project license. DGL-LifeSci is licensed under Apache-2.0.
