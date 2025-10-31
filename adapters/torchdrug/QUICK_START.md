# TorchDrug Adapter - Quick Start Guide

## Installation

```bash
pip install torch torchdrug
```

## 30-Second Start

```python
from adapters.torchdrug import TorchDrugAdapter

# Initialize
adapter = TorchDrugAdapter()

# Generate embeddings
result = await adapter.execute("CCO")  # Ethanol

if result.success:
    embeddings = result.data["embeddings"]
    print(f"Embedding dimension: {result.data['embedding_dim']}")
```

## Common Use Cases

### 1. Single Molecule Embedding

```python
result = await adapter.execute("c1ccccc1")  # Benzene
embeddings = result.data["embeddings"]
```

### 2. Batch Processing

```python
molecules = ["CCO", "CC(=O)O", "c1ccccc1"]
result = await adapter.execute(molecules)
embeddings = result.data["embeddings"]  # Shape: [3, 256]
```

### 3. Different GNN Models

```python
# Graph Isomorphism Network (default)
result = await adapter.execute("CCO", model_type="gin")

# Graph Attention Network
result = await adapter.execute("CCO", model_type="gat")

# Graph Convolutional Network
result = await adapter.execute("CCO", model_type="gcn")
```

### 4. Graph Features Only

```python
result = await adapter.execute(
    "CCO",
    generate_embeddings=False,
    include_graph_features=True
)

features = result.data["graph_features"]
# Contains: num_atoms, num_bonds, atom_feature_dim, etc.
```

### 5. Pre-trained Model

```python
result = await adapter.execute(
    "CCO",
    model_type="gin",
    pretrained_model_path="path/to/model.pth"
)
```

## Available Models

- **gin** - Graph Isomorphism Network (default, recommended)
- **gat** - Graph Attention Network
- **gcn** - Graph Convolutional Network
- **schnet** - SchNet (3D modeling)
- **rgcn** - Relational GCN
- **graphsage** - GraphSAGE
- **nfp** - Neural Fingerprint
- **mpnn** - Message Passing NN

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `model_type` | "gin" | GNN architecture |
| `generate_embeddings` | True | Generate embeddings |
| `include_graph_features` | True | Include graph stats |
| `batch_size` | 32 | Batch size |
| `hidden_dims` | [256, 256, 256] | Hidden layers |

## Output Format

```python
{
    "smiles": "CCO",
    "embeddings": [...],  # List of floats
    "embedding_dim": 256,
    "embedding_stats": {
        "mean": 0.123,
        "std": 0.456,
        "min": -1.234,
        "max": 2.345
    },
    "graph_features": {
        "num_atoms": 9,
        "num_bonds": 8,
        "atom_feature_dim": 69,
        "bond_feature_dim": 4
    }
}
```

## Error Handling

```python
result = await adapter.execute("INVALID")

if not result.success:
    print(f"Error: {result.error}")
```

## Integration with PharmForge

```python
# With caching
result = await adapter("CCO", use_cache=True)

# Register with PharmForge
from backend.core.adapters.protocol import registry
registry.register(TorchDrugAdapter())
```

## GPU Support

Automatically uses CUDA if available:

```python
adapter = TorchDrugAdapter()
# Check device in metadata
result = await adapter.execute("CCO")
print(result.metadata["device"])  # "cuda" or "cpu"
```

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Check [example_usage.py](example_usage.py) for more examples
- See [adapter.py](adapter.py) for implementation details

## Troubleshooting

**ImportError: TorchDrug not installed**
```bash
pip install torchdrug
```

**CUDA out of memory**
```python
# Reduce batch size
result = await adapter.execute(molecules, batch_size=8)
```

**Invalid SMILES**
```python
# Check result.success before using data
if result.success:
    # Use result.data
else:
    # Handle error
    print(result.error)
```
