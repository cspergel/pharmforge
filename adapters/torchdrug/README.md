# TorchDrug Adapter

GNN-based molecular modeling and feature extraction adapter for PharmForge using TorchDrug library.

## Overview

The TorchDrug adapter provides molecular embeddings and features using Graph Neural Network (GNN) models from the TorchDrug library. It supports multiple GNN architectures and can load pre-trained models for molecular property prediction tasks.

## Installation

```bash
# Install PyTorch (required)
pip install torch

# Install TorchDrug
pip install torchdrug
```

## Features

- **Molecular graph generation from SMILES**
- **Multiple GNN architectures** (GCN, GAT, GIN, SchNet, RGCN, GraphSAGE, NFP, MPNN)
- **Pre-trained model loading**
- **Embedding generation** for molecular property prediction
- **Graph feature extraction** (atom counts, bond counts, feature dimensions)
- **Batch processing** for efficient computation
- **GPU support** (automatic CUDA detection)

## Supported GNN Models

| Model | Description | Use Case |
|-------|-------------|----------|
| **GCN** | Graph Convolutional Network | General molecular property prediction |
| **GAT** | Graph Attention Network | Tasks requiring attention mechanisms |
| **GIN** | Graph Isomorphism Network | Powerful for graph classification (default) |
| **SchNet** | Continuous-filter CNN | 3D molecular modeling |
| **RGCN** | Relational GCN | Multiple edge types in molecular graphs |
| **GraphSAGE** | Graph Sample and Aggregate | Inductive learning |
| **NFP** | Neural Fingerprint | Learnable molecular fingerprints |
| **MPNN** | Message Passing Neural Network | General GNN framework |

## Usage

### Basic Usage

```python
from adapters.torchdrug import TorchDrugAdapter

# Initialize adapter
adapter = TorchDrugAdapter()

# Generate embeddings for a single molecule
result = await adapter.execute("CCO")  # Ethanol

if result.success:
    embeddings = result.data["embeddings"]
    graph_features = result.data["graph_features"]
    print(f"Embedding dimension: {result.data['embedding_dim']}")
    print(f"Number of atoms: {graph_features['num_atoms']}")
```

### Batch Processing

```python
# Process multiple molecules
smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]  # Ethanol, Acetic acid, Benzene

result = await adapter.execute(smiles_list)

if result.success:
    embeddings = result.data["embeddings"]  # Shape: [3, embedding_dim]
    print(f"Generated embeddings for {result.data['n_molecules']} molecules")
```

### Using Different GNN Models

```python
# Use Graph Attention Network (GAT)
result = await adapter.execute(
    "CCO",
    model_type="gat"
)

# Use SchNet for 3D modeling
result = await adapter.execute(
    "CCO",
    model_type="schnet"
)

# Use GCN with custom hidden dimensions
result = await adapter.execute(
    "CCO",
    model_type="gcn",
    hidden_dims=[128, 128, 128]
)
```

### Loading Pre-trained Models

```python
# Load a pre-trained model checkpoint
result = await adapter.execute(
    "CCO",
    model_type="gin",
    pretrained_model_path="/path/to/model_checkpoint.pth"
)
```

### Graph Features Only (No Embeddings)

```python
# Extract only graph features without generating embeddings
result = await adapter.execute(
    "CCO",
    generate_embeddings=False,
    include_graph_features=True
)

if result.success:
    features = result.data["graph_features"]
    print(f"Atoms: {features['num_atoms']}")
    print(f"Bonds: {features['num_bonds']}")
    print(f"Atom feature dim: {features['atom_feature_dim']}")
    print(f"Bond feature dim: {features['bond_feature_dim']}")
```

### Advanced Configuration

```python
# Full configuration example
result = await adapter.execute(
    smiles_list,
    model_type="gin",
    generate_embeddings=True,
    include_graph_features=True,
    batch_size=64,
    hidden_dims=[256, 256, 256],
    pretrained_model_path=None,
    return_graph_objects=False  # Set True to get raw TorchDrug molecules
)
```

## Return Format

### Successful Result

```python
AdapterResult(
    success=True,
    data={
        "smiles": "CCO",  # or list for batch
        "n_molecules": 1,
        "n_failed": 0,
        "embeddings": [...],  # List of float values
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
            "bond_feature_dim": 4,
            "unique_atom_types": 2,
            "unique_bond_types": 1
        }
    },
    metadata={
        "source": "torchdrug",
        "adapter_version": "1.0.0",
        "computation_type": "local",
        "model_type": "gin",
        "model_description": "Graph Isomorphism Network...",
        "n_molecules": 1,
        "n_failed": 0,
        "embeddings_generated": True,
        "device": "cuda",  # or "cpu"
        "torchdrug_version": "0.2.0",
        "torch_version": "2.0.0"
    }
)
```

### Error Result

```python
AdapterResult(
    success=False,
    data=None,
    error="TorchDrug is not installed. Install with: pip install torchdrug"
)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_type` | str | "gin" | GNN architecture to use |
| `generate_embeddings` | bool | True | Whether to generate GNN embeddings |
| `include_graph_features` | bool | True | Include basic graph statistics |
| `pretrained_model_path` | str | None | Path to pre-trained model checkpoint |
| `hidden_dims` | List[int] | [256, 256, 256] | Hidden layer dimensions |
| `batch_size` | int | 32 | Batch size for embedding generation |
| `return_graph_objects` | bool | False | Return raw TorchDrug molecule objects |

## Integration with PharmForge

### Registration

```python
from backend.core.adapters.protocol import registry
from adapters.torchdrug import TorchDrugAdapter

# Register the adapter
adapter = TorchDrugAdapter()
registry.register(adapter)
```

### Using with Caching

```python
# Adapter automatically uses PharmForge's caching system
result = await adapter(
    "CCO",
    use_cache=True,  # Enable caching
    model_type="gin"
)

if result.cache_hit:
    print("Result retrieved from cache!")
```

## Performance Considerations

1. **GPU Acceleration**: Automatically uses CUDA if available
2. **Batch Processing**: Use lists of SMILES for better throughput
3. **Model Caching**: Models are cached after first creation
4. **Embedding Caching**: Results are cached via PharmForge's cache system

## Example Use Cases

### Molecular Similarity Search

```python
# Generate embeddings for a database of molecules
database_smiles = ["CCO", "CC(=O)O", "c1ccccc1", ...]
result = await adapter.execute(database_smiles, model_type="gin")

embeddings = np.array(result.data["embeddings"])

# Generate embedding for query molecule
query_result = await adapter.execute("CC(C)O", model_type="gin")
query_embedding = np.array(query_result.data["embeddings"])

# Calculate cosine similarity
from sklearn.metrics.pairwise import cosine_similarity
similarities = cosine_similarity([query_embedding], embeddings)[0]
```

### Property Prediction Pipeline

```python
# Generate features for downstream ML models
molecules = ["CCO", "CC(=O)O", "c1ccccc1"]

# Extract embeddings
result = await adapter.execute(
    molecules,
    model_type="gin",
    generate_embeddings=True
)

embeddings = np.array(result.data["embeddings"])

# Use embeddings for property prediction
from sklearn.ensemble import RandomForestRegressor
model = RandomForestRegressor()
# ... train model with embeddings as features
```

### Transfer Learning

```python
# Use pre-trained model for specific task
result = await adapter.execute(
    molecules,
    model_type="gin",
    pretrained_model_path="pretrained_models/molecular_properties.pth"
)

# Fine-tune or use as feature extractor
embeddings = result.data["embeddings"]
```

## Error Handling

```python
result = await adapter.execute("INVALID_SMILES")

if not result.success:
    print(f"Error: {result.error}")
    if result.metadata and "failed_smiles" in result.metadata:
        print(f"Failed SMILES: {result.metadata['failed_smiles']}")
```

## Dependencies

- **PyTorch** >= 1.9.0
- **TorchDrug** >= 0.2.0
- **NumPy** >= 1.19.0

## Limitations

1. Requires PyTorch and TorchDrug installation
2. GPU recommended for large-scale processing
3. Some models may require 3D coordinates (e.g., SchNet)
4. Pre-trained models must match the specified architecture

## References

- [TorchDrug Documentation](https://torchdrug.ai/)
- [TorchDrug Paper](https://arxiv.org/abs/2202.08320)
- [Graph Neural Networks for Drug Discovery](https://arxiv.org/abs/1904.01561)

## Version History

- **1.0.0** (2025-10-30): Initial release
  - Multiple GNN architecture support
  - Pre-trained model loading
  - Batch processing
  - GPU acceleration
  - Graph feature extraction
