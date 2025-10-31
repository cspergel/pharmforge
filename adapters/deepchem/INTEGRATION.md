# DeepChem Adapter Integration Guide

## Quick Start

### 1. Register the Adapter

The DeepChem adapter should be automatically registered when PharmForge loads adapters. To manually register:

```python
from adapters.deepchem.adapter import DeepChemAdapter
from backend.core.adapters.protocol import registry

# Create and register adapter instance
deepchem_adapter = DeepChemAdapter()
registry.register(deepchem_adapter)
```

### 2. Use the Adapter

```python
# Get adapter from registry
adapter = registry.get("deepchem")

# Execute with default settings (Morgan fingerprints)
result = await adapter("CCO")  # Uses __call__ method with caching

# Or use execute directly
result = await adapter.execute("CCO", featurizer_type="rdkit_descriptors")
```

## Configuration

The adapter can be configured when initializing:

```python
adapter = DeepChemAdapter()
adapter.config.update({
    "timeout": 120,  # Increase timeout for large molecules
    "default_featurizer": "mol_graph_conv",  # Change default
    "default_radius": 3,
    "default_size": 4096
})
```

## Caching

The adapter supports automatic caching through the PharmForge caching system:

```python
# First call - computes features and caches
result1 = await adapter("CCO", featurizer_type="morgan")
print(result1.cache_hit)  # False

# Second call - retrieves from cache
result2 = await adapter("CCO", featurizer_type="morgan")
print(result2.cache_hit)  # True

# Disable caching for specific call
result3 = await adapter("CCO", use_cache=False)
```

## Validation

The adapter validates input before execution:

```python
# Check if input is valid
is_valid = adapter.validate_input("CCO")  # True
is_valid = adapter.validate_input("invalid")  # False
is_valid = adapter.validate_input("")  # False
is_valid = adapter.validate_input(123)  # False
```

## API Integration

### REST API Endpoint

If PharmForge has an adapter API endpoint, you can use:

```bash
# POST request to execute adapter
curl -X POST http://localhost:8000/api/adapters/deepchem/execute \
  -H "Content-Type: application/json" \
  -d '{
    "input_data": "CCO",
    "featurizer_type": "morgan",
    "radius": 2,
    "size": 2048
  }'
```

### Response Format

```json
{
  "success": true,
  "data": {
    "features": [0, 1, 0, 1, ...],
    "smiles": "CCO",
    "featurizer_type": "morgan",
    "n_molecules": 1,
    "feature_shape": [1, 2048],
    "dtype": "int64",
    "is_binary": true,
    "sparsity": 0.95
  },
  "cache_hit": false,
  "metadata": {
    "source": "deepchem",
    "adapter_version": "1.0.0",
    "computation_type": "local",
    "featurizer_type": "morgan",
    "n_molecules": 1,
    "parameters": {
      "radius": 2,
      "size": 2048
    },
    "deepchem_version": "2.8.0",
    "cache_key": "a1b2c3d4..."
  }
}
```

## Batch Processing

Process multiple molecules efficiently:

```python
# Batch processing
smiles_list = [
    "CCO",           # Ethanol
    "CC(=O)O",       # Acetic acid
    "c1ccccc1",      # Benzene
    "CC(C)CC1=CC=C(C=C1)C(C)C"  # Ibuprofen
]

result = await adapter.execute(
    smiles_list,
    featurizer_type="morgan",
    radius=2,
    size=2048
)

if result.success:
    features = result.data["features"]  # List of feature vectors
    print(f"Processed {len(features)} molecules")
```

## Error Handling

The adapter returns structured errors:

```python
result = await adapter.execute("invalid_smiles")

if not result.success:
    print(f"Error: {result.error}")
    print(f"Metadata: {result.metadata}")
```

Common errors:
- `"DeepChem is not installed. Install with: pip install deepchem"`
- `"RDKit is not installed. Install with: pip install rdkit"`
- `"Invalid SMILES string(s) or could not parse with RDKit"`
- `"Unsupported featurizer type: xyz. Supported types: [...]"`
- `"Failed to generate features: ..."`

## Workflow Integration

### With Other Adapters

Combine DeepChem features with other adapters:

```python
from backend.core.adapters.protocol import registry

# Get multiple adapters
deepchem = registry.get("deepchem")
rdkit = registry.get("rdkit_local")
pubchem = registry.get("pubchem")

# Combine data
smiles = "CCO"

# Get features from DeepChem
features_result = await deepchem(smiles, featurizer_type="morgan")

# Get properties from RDKit
props_result = await rdkit(smiles)

# Get data from PubChem
pubchem_result = await pubchem(smiles)

# Combine results
combined = {
    "features": features_result.data["features"],
    "properties": props_result.data,
    "pubchem_data": pubchem_result.data
}
```

### Machine Learning Pipeline

```python
# Step 1: Feature generation
training_smiles = [...]  # Your training data
feature_result = await deepchem.execute(
    training_smiles,
    featurizer_type="morgan",
    radius=3,
    size=4096
)

X_train = np.array(feature_result.data["features"])

# Step 2: Train model
from sklearn.ensemble import RandomForestClassifier
y_train = [...]  # Your labels

model = RandomForestClassifier(n_estimators=100)
model.fit(X_train, y_train)

# Step 3: Predict on new molecules
test_smiles = [...]
feature_result = await deepchem.execute(
    test_smiles,
    featurizer_type="morgan",
    radius=3,
    size=4096
)

X_test = np.array(feature_result.data["features"])
predictions = model.predict(X_test)
```

## Metadata Access

Get adapter information:

```python
# Get adapter metadata
metadata = adapter.get_metadata()
print(metadata)
# {
#     "name": "deepchem",
#     "type": "local",
#     "version": "1.0.0",
#     "enabled": True,
#     "config": {...}
# }

# Check if adapter is enabled
if adapter.enabled:
    result = await adapter(smiles)
```

## Advanced Usage

### Custom Featurizer Parameters

```python
# Graph convolution with custom settings
result = await adapter.execute(
    "CCO",
    featurizer_type="mol_graph_conv",
    use_edges=True,
    return_raw=True  # Get raw DeepChem features
)

# 3D features (requires 3D coordinates)
result = await adapter.execute(
    "CCO",
    featurizer_type="coulomb_matrix",
    max_atoms=50
)
```

### Cache Key Generation

```python
# Generate cache key for specific parameters
cache_key = adapter.generate_cache_key(
    "CCO",
    featurizer_type="morgan",
    radius=2,
    size=2048
)
print(cache_key)  # SHA256 hash
```

## Testing

Example test cases:

```python
import pytest

@pytest.mark.asyncio
async def test_deepchem_adapter():
    adapter = DeepChemAdapter()

    # Test basic execution
    result = await adapter.execute("CCO")
    assert result.success
    assert "features" in result.data

    # Test invalid input
    result = await adapter.execute("invalid")
    assert not result.success
    assert "Invalid SMILES" in result.error

    # Test different featurizers
    for ftype in ["morgan", "maccs", "rdkit_descriptors"]:
        result = await adapter.execute("CCO", featurizer_type=ftype)
        assert result.success

@pytest.mark.asyncio
async def test_deepchem_batch():
    adapter = DeepChemAdapter()

    smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
    result = await adapter.execute(smiles_list)

    assert result.success
    assert len(result.data["features"]) == 3
```

## Performance Tips

1. **Use appropriate featurizer**:
   - Fingerprints for similarity/classification tasks
   - Graph features for GNN models
   - Descriptors for interpretability

2. **Batch processing**:
   - Process molecules in batches for better performance
   - Use list input instead of individual calls

3. **Caching**:
   - Enable caching for repeated queries
   - Use same parameters for cache hits

4. **Feature size**:
   - Larger fingerprints = more information but slower
   - Start with 2048 bits, increase if needed

## Troubleshooting

### Import Errors
```python
# Check if dependencies are available
from adapters.deepchem.adapter import DEEPCHEM_AVAILABLE, RDKIT_AVAILABLE

if not DEEPCHEM_AVAILABLE:
    print("Install: pip install deepchem")
if not RDKIT_AVAILABLE:
    print("Install: pip install rdkit")
```

### Memory Issues
```python
# For large datasets, process in chunks
def process_in_chunks(smiles_list, chunk_size=100):
    results = []
    for i in range(0, len(smiles_list), chunk_size):
        chunk = smiles_list[i:i+chunk_size]
        result = await adapter.execute(chunk)
        if result.success:
            results.extend(result.data["features"])
    return results
```

### Graph Features Serialization
```python
# Graph features are not JSON-serializable
# Use return_raw=True and handle separately
result = await adapter.execute(
    "CCO",
    featurizer_type="mol_graph_conv",
    return_raw=True
)

# Don't try to JSON serialize graph features
# Use them directly with DeepChem models instead
```

## References

- [PharmForge Adapter Protocol](../../backend/core/adapters/protocol.py)
- [DeepChem Documentation](https://deepchem.readthedocs.io/)
- [Example Adapters](../rdkit_local/adapter.py)
