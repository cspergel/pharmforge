# DeepChem Adapter

## Overview
The DeepChem adapter provides molecular featurization using the DeepChem library, a popular deep learning library for chemistry and biology. It generates feature vectors suitable for machine learning models in drug discovery and molecular property prediction.

## Installation
```bash
pip install deepchem
```

Note: DeepChem requires RDKit, which will be installed automatically as a dependency.

## Supported Featurizers

### Fingerprint-Based Featurizers
- **morgan**: Morgan fingerprints (ECFP-like) with customizable radius and size
- **circular**: Circular fingerprints (similar to Morgan)
- **maccs**: 166-bit MACCS structural keys

### Descriptor-Based Featurizers
- **rdkit_descriptors**: 200+ physicochemical descriptors from RDKit

### Graph-Based Featurizers (for GNN models)
- **weave**: Graph convolution features for Weave models
- **graph_conv**: Graph convolution features for GraphConv models
- **mol_graph_conv**: Molecular graph features for GCN models

### 3D Structure-Based Featurizers (requires 3D coordinates)
- **coulomb_matrix**: Coulomb matrix representation
- **coulomb_matrix_eig**: Eigenvalues of Coulomb matrix

## Usage Examples

### Basic Usage (Morgan Fingerprints)
```python
from adapters.deepchem.adapter import DeepChemAdapter

adapter = DeepChemAdapter()

# Single molecule
result = await adapter.execute("CCO")  # Ethanol
if result.success:
    features = result.data["features"]  # 2048-bit Morgan fingerprint
    print(f"Generated {len(features)} features")

# Multiple molecules
result = await adapter.execute(["CCO", "CC(=O)O", "c1ccccc1"])
if result.success:
    features = result.data["features"]  # List of fingerprints
    print(f"Generated features for {len(features)} molecules")
```

### Customizing Featurizer Type
```python
# Use MACCS keys instead
result = await adapter.execute(
    "CCO",
    featurizer_type="maccs"
)

# Use RDKit descriptors
result = await adapter.execute(
    "CCO",
    featurizer_type="rdkit_descriptors"
)

# Use graph convolution features
result = await adapter.execute(
    "CCO",
    featurizer_type="mol_graph_conv"
)
```

### Customizing Fingerprint Parameters
```python
# Custom Morgan fingerprint parameters
result = await adapter.execute(
    "CCO",
    featurizer_type="morgan",
    radius=3,      # Larger radius
    size=4096      # More bits
)
```

### Working with Graph Features (Raw Format)
```python
# Get raw DeepChem graph features (not JSON-serializable)
result = await adapter.execute(
    "CCO",
    featurizer_type="mol_graph_conv",
    return_raw=True
)

if result.success:
    raw_features = result.data["features"]
    # Use these features directly with DeepChem models
```

## Result Format

The adapter returns an `AdapterResult` with the following structure:

```python
{
    "success": True,
    "data": {
        "features": [...],        # Generated features (array or list)
        "smiles": "CCO",          # Input SMILES
        "featurizer_type": "morgan",
        "n_molecules": 1,
        "feature_shape": (1, 2048),
        "dtype": "int64",
        "is_binary": True,
        "sparsity": 0.95          # For binary fingerprints
    },
    "metadata": {
        "source": "deepchem",
        "adapter_version": "1.0.0",
        "computation_type": "local",
        "featurizer_type": "morgan",
        "n_molecules": 1,
        "parameters": {"radius": 2, "size": 2048},
        "deepchem_version": "2.8.0"
    }
}
```

## Integration with Machine Learning

### Using with scikit-learn
```python
# Generate features for training data
smiles_train = ["CCO", "CC(=O)O", "c1ccccc1", ...]
result = await adapter.execute(
    smiles_train,
    featurizer_type="morgan"
)

X_train = np.array(result.data["features"])
y_train = [...]  # Your target values

# Train a model
from sklearn.ensemble import RandomForestClassifier
model = RandomForestClassifier()
model.fit(X_train, y_train)
```

### Using with DeepChem Models
```python
# Generate graph features for GNN
result = await adapter.execute(
    smiles_train,
    featurizer_type="mol_graph_conv",
    return_raw=True
)

# Use with DeepChem GraphConvModel
import deepchem as dc
model = dc.models.GraphConvModel(...)
```

## Performance Notes

- **Local computation**: All featurization is done locally (no API calls)
- **Speed**: Fingerprints are fast; graph features are slower
- **Memory**: Graph features can be memory-intensive for large molecules
- **Caching**: Results are automatically cached for repeated queries

## Error Handling

The adapter handles errors gracefully:

```python
# Missing dependency
result = await adapter.execute("CCO")
if not result.success:
    print(result.error)  # "DeepChem is not installed. Install with: pip install deepchem"

# Invalid SMILES
result = await adapter.execute("invalid_smiles")
if not result.success:
    print(result.error)  # "Invalid SMILES string(s) or could not parse with RDKit"

# Unsupported featurizer
result = await adapter.execute("CCO", featurizer_type="unknown")
if not result.success:
    print(result.error)  # "Unsupported featurizer type: unknown. Supported types: [...]"
```

## Adapter Metadata

- **Name**: `deepchem`
- **Type**: `local`
- **Version**: `1.0.0`
- **Dependencies**: `deepchem`, `rdkit`, `numpy`

## Comparison with Other Adapters

### vs. RDKit Adapter
- **DeepChem**: Provides more featurization options (graphs, 3D features)
- **RDKit**: Focused on molecular properties and simple fingerprints

### vs. scikit-mol Adapter
- **DeepChem**: Better for deep learning workflows, graph neural networks
- **scikit-mol**: Better for traditional ML with scikit-learn integration

### vs. Chemprop Adapter
- **DeepChem**: More general-purpose, multiple featurization options
- **Chemprop**: Specialized for message-passing neural networks

## References

- [DeepChem Documentation](https://deepchem.readthedocs.io/)
- [DeepChem Featurizers](https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html)
- [DeepChem GitHub](https://github.com/deepchem/deepchem)
