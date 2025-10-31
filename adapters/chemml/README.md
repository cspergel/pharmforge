# ChemML Adapter

Machine learning workflows for chemical and materials data with specialized molecular representations.

## Overview

The ChemML adapter provides access to ChemML's powerful molecular representation methods that bridge chemistry and machine learning. It generates specialized features optimized for ML models in chemistry and materials science.

## Installation

```bash
pip install chemml
pip install rdkit  # Required dependency
```

## Features

### Molecular Representations

1. **Coulomb Matrix**
   - Representation of molecular structure based on nuclear charges and distances
   - Captures 3D molecular geometry
   - Invariant to translation and rotation

2. **Bag of Bonds**
   - Groups atoms by element type and bond distances
   - Captures local chemical environment
   - Efficient for property prediction

3. **Morgan Fingerprints**
   - Circular fingerprints via RDKit integration
   - Configurable radius and bit length
   - Fast and widely used

4. **RDKit Descriptors**
   - 200+ physicochemical descriptors
   - Comprehensive molecular properties
   - No 3D coordinates required

## Usage

### Basic Usage

```python
from adapters.chemml import ChemMLAdapter

adapter = ChemMLAdapter()

# Single molecule
input_data = {
    "smiles": "CCO",
    "operation": "coulomb_matrix",
    "parameters": {
        "max_n_atoms": 23,
        "sorting": "row_norm"
    }
}

result = await adapter(input_data)

if result.success:
    print(f"Shape: {result.data['shape']}")
    print(f"Representations: {result.data['representations']}")
```

### Batch Processing

```python
# Multiple molecules
input_data = {
    "smiles": ["CCO", "CC(=O)O", "c1ccccc1"],
    "operation": "morgan_fingerprint",
    "parameters": {
        "radius": 2,
        "n_bits": 2048
    }
}

result = await adapter(input_data)

if result.success:
    print(f"Generated {result.data['num_successful']} representations")
    print(f"Feature dimension: {result.data['feature_dimension']}")
```

## Supported Operations

### 1. Coulomb Matrix

Generates Coulomb matrix representation of molecular structure.

**Requires**: 3D coordinates (generated automatically from SMILES)

```python
{
    "smiles": "CCO",
    "operation": "coulomb_matrix",
    "parameters": {
        "max_n_atoms": 23,        # Maximum atoms in matrix
        "sorting": "row_norm",     # "row_norm", "unsorted", or "random"
        "padding": 0.0            # Padding value for smaller molecules
    }
}
```

**Output**:
```python
{
    "representations": [[...], [...]],  # Flattened matrices
    "representation_type": "coulomb_matrix",
    "shape": [10, 529],                 # 10 molecules, 23x23 flattened
    "num_molecules": 10,
    "num_successful": 10,
    "failed_indices": [],
    "feature_dimension": 529,
    "parameters": {...}
}
```

### 2. Bag of Bonds

Generates bag of bonds representation.

**Requires**: 3D coordinates (generated automatically from SMILES)

```python
{
    "smiles": ["CCO", "CC(=O)O"],
    "operation": "bag_of_bonds",
    "parameters": {
        "const": 1.0  # Constant for BoB calculation
    }
}
```

### 3. Morgan Fingerprints

Generates Morgan fingerprints using RDKit.

**Requires**: Only SMILES (no 3D coordinates needed)

```python
{
    "smiles": ["CCO", "CC(=O)O"],
    "operation": "morgan_fingerprint",
    "parameters": {
        "radius": 2,      # Fingerprint radius
        "n_bits": 2048    # Number of bits
    }
}
```

### 4. RDKit Descriptors

Calculates 200+ RDKit molecular descriptors.

**Requires**: Only SMILES (no 3D coordinates needed)

```python
{
    "smiles": ["CCO", "CC(=O)O"],
    "operation": "rdkit_descriptors"
}
```

**Output**:
```python
{
    "descriptors": [[...], [...]],
    "descriptor_names": ["MolWt", "LogP", ...],
    "representation_type": "rdkit_descriptors",
    "shape": [2, 208],
    "num_descriptors": 208
}
```

## Input Format

### Dictionary Format (Recommended)

```python
{
    "smiles": str or List[str],           # Required
    "operation": str,                      # Optional, default: "coulomb_matrix"
    "parameters": Dict[str, Any]           # Optional, operation-specific params
}
```

### Simple Format

```python
# Single SMILES string
"CCO"

# List of SMILES
["CCO", "CC(=O)O", "c1ccccc1"]
```

When using simple format, specify operation via kwargs:
```python
result = await adapter("CCO", operation="morgan_fingerprint", radius=3)
```

## Output Structure

All operations return:

```python
{
    "representations": List[List[float]],  # or "descriptors" for RDKit descriptors
    "representation_type": str,
    "shape": List[int],
    "num_molecules": int,
    "num_successful": int,
    "failed_indices": List[int],           # Indices of failed molecules
    "feature_dimension": int,
    "parameters": Dict[str, Any],
    "smiles": str or List[str]
}
```

## Error Handling

The adapter handles various error conditions:

1. **Missing Dependencies**
   ```python
   if not result.success:
       print(result.error)  # "RDKit is not installed..."
   ```

2. **Invalid SMILES**
   - Invalid molecules are skipped
   - `failed_indices` contains their positions
   - Only successful molecules are returned

3. **3D Generation Failures**
   - For 3D operations, some molecules may fail coordinate generation
   - These are logged and skipped
   - Check `failed_indices` for details

## Performance Considerations

1. **3D Operations are Slower**
   - Coulomb matrix and bag of bonds require 3D coordinate generation
   - ~0.1-1 second per molecule depending on size
   - Consider using Morgan fingerprints for faster processing

2. **Batch Processing**
   - Process multiple molecules at once for efficiency
   - Maximum 1000 molecules per batch (configurable)

3. **Memory Usage**
   - Coulomb matrices can be large for big molecules
   - Adjust `max_n_atoms` based on your dataset
   - Typical values: 23 for small molecules, 50+ for larger molecules

## Integration with ML Pipelines

### Scikit-learn Integration

```python
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import numpy as np

# Generate features
input_data = {
    "smiles": smiles_list,
    "operation": "morgan_fingerprint"
}

result = await adapter(input_data)

if result.success:
    X = np.array(result.data['representations'])
    y = property_values  # Your target property

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

    model = RandomForestRegressor()
    model.fit(X_train, y_train)
    score = model.score(X_test, y_test)
```

### Feature Engineering Pipeline

```python
# Generate multiple representations
operations = ["morgan_fingerprint", "rdkit_descriptors", "coulomb_matrix"]

all_features = []
for operation in operations:
    result = await adapter({
        "smiles": smiles_list,
        "operation": operation
    })
    if result.success:
        all_features.append(result.data['representations'])

# Concatenate features
X_combined = np.hstack(all_features)
```

## Comparison with Other Adapters

| Adapter | Focus | 3D Required | Speed | Best For |
|---------|-------|-------------|-------|----------|
| ChemML | Specialized ML representations | Some ops | Medium | Property prediction, ML pipelines |
| DeepChem | Deep learning features | No | Fast | Neural networks, graph models |
| Mordred | Comprehensive descriptors | No | Slow | Descriptor-based QSAR |
| RDKit | General cheminformatics | No | Fast | Quick calculations, screening |

## Troubleshooting

### ChemML Not Available

```python
from adapters.chemml import CHEMML_AVAILABLE

if not CHEMML_AVAILABLE:
    print("Install ChemML: pip install chemml")
```

### 3D Generation Failures

If molecules consistently fail 3D generation:
- Check SMILES validity
- Try simpler molecules first
- Use 2D operations (Morgan fingerprints, RDKit descriptors)

### Memory Issues

For large datasets:
- Process in smaller batches
- Use Morgan fingerprints instead of Coulomb matrices
- Reduce `max_n_atoms` for Coulomb matrices

## References

1. ChemML Documentation: https://chemml.io/
2. Coulomb Matrix: Rupp et al., PRL 108, 058301 (2012)
3. Bag of Bonds: Hansen et al., J. Phys. Chem. Lett. 6, 2326 (2015)

## Version History

- 1.0.0 - Initial release with Coulomb matrix, bag of bonds, Morgan fingerprints, and RDKit descriptors

## License

Part of PharmForge adapter collection.
