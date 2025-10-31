# DGL-LifeSci Adapter - Quick Start Guide

## Installation

```bash
pip install dgllife
```

## Basic Usage

```python
from adapters.dgllifesci.adapter import DGLLifeSciAdapter

# Initialize
adapter = DGLLifeSciAdapter()

# Generate molecular graph
result = await adapter.execute("CCO", featurizer_type="canonical")

if result.success:
    print(f"Nodes: {result.data['graph_features']['num_nodes']}")
    print(f"Edges: {result.data['graph_features']['num_edges']}")
```

## Key Features

### 1. Multiple Featurizers

```python
# Canonical (default) - 74 atom features, 12 bond features
result = await adapter.execute("CCO", featurizer_type="canonical")

# AttentiveFP - optimized for attention-based GNNs
result = await adapter.execute("CCO", featurizer_type="attentivefp")

# Weave - for Weave convolutional networks
result = await adapter.execute("CCO", featurizer_type="weave")

# Minimal - baseline features
result = await adapter.execute("CCO", featurizer_type="minimal")
```

### 2. Batch Processing

```python
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]
result = await adapter.execute(smiles_list, featurizer_type="canonical")

print(f"Processed: {result.data['n_successful']} molecules")
```

### 3. Pre-trained Models

```python
# Load model
adapter.load_model("path/to/gcn_model.pt", model_type="gcn")

# Make predictions
result = await adapter.execute(
    "CCO",
    featurizer_type="canonical",
    include_predictions=True
)
```

### 4. Raw DGL Graphs

```python
result = await adapter.execute(
    "CCO",
    featurizer_type="canonical",
    return_graphs=True
)

graph = result.data['graphs']
node_features = graph.ndata['h']  # Access node features
edge_features = graph.edata['e']  # Access edge features
```

## Supported Featurizers

| Featurizer | Node Features | Edge Features | Use Case |
|------------|---------------|---------------|----------|
| canonical | 74 | 12 | General-purpose |
| attentivefp | Variable | Variable | Attention-based GNNs |
| weave | Variable | Variable | Weave networks |
| minimal | Minimal | Minimal | Baseline |

## Supported Models

- **GCN**: Graph Convolutional Network
- **GAT**: Graph Attention Network
- **AttentiveFP**: Attentive Fingerprint GNN
- **GIN**: Graph Isomorphism Network
- **MPNN**: Message Passing Neural Network

## Common Patterns

### Property Prediction

```python
# Load solubility model
adapter.load_model("models/gcn_solubility.pt", model_type="gcn")

# Predict for multiple compounds
result = await adapter.execute(
    ["CCO", "c1ccccc1"],
    featurizer_type="canonical",
    include_predictions=True
)
```

### Feature Extraction

```python
# Extract graph features for ML
result = await adapter.execute(
    smiles_list,
    featurizer_type="canonical"
)

features = result.data['graph_features']
# Use features in downstream models
```

### Error Handling

```python
result = await adapter.execute(smiles_list)

if result.success:
    if result.data['n_failed'] > 0:
        print(f"Failed SMILES: {result.data['failed_smiles']}")
else:
    print(f"Error: {result.error}")
```

## Response Format

```python
{
    "success": True,
    "data": {
        "featurizer_type": "canonical",
        "n_molecules": 1,
        "n_successful": 1,
        "n_failed": 0,
        "graph_features": {
            "smiles": "CCO",
            "num_nodes": 9,
            "num_edges": 16,
            "node_features": {
                "shape": [9, 74],
                "feature_dim": 74,
                "mean": 0.234,
                "std": 0.412
            },
            "edge_features": {
                "shape": [16, 12],
                "feature_dim": 12,
                "mean": 0.156,
                "std": 0.298
            }
        }
    },
    "metadata": {
        "source": "dgllifesci",
        "featurizer_type": "canonical",
        "n_molecules": 1,
        "predictions_included": False
    }
}
```

## Configuration

```python
adapter.config.update({
    "timeout": 120,
    "default_featurizer": "attentivefp",
    "device": "cuda",  # Use GPU
    "batch_size": 64
})
```

## Tips

1. **Start with canonical**: Use the canonical featurizer for most tasks
2. **Batch processing**: Process multiple molecules together for efficiency
3. **GPU acceleration**: Set `device="cuda"` for faster inference
4. **Error handling**: Always check for failed molecules in batch processing
5. **Model compatibility**: Ensure featurizer matches your model's training featurizer

## Examples

See `example_usage.py` for comprehensive examples including:
- Basic graph generation
- Multiple featurizers comparison
- Batch processing
- Raw graph access
- Model predictions
- Error handling

## Troubleshooting

**"DGL-LifeSci is not installed"**
```bash
pip install dgllife
```

**"RDKit not available"**
```bash
pip install rdkit
```

**Some molecules fail to convert**
```python
# Check failed SMILES
if result.data['n_failed'] > 0:
    print(result.data['failed_smiles'])
```

## References

- [DGL-LifeSci Documentation](https://lifesci.dgl.ai/)
- [Full README](README.md)
- [Example Usage](example_usage.py)
