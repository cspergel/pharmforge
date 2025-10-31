# Oloren ChemEngine Adapter

State-of-the-art molecular property prediction using modern machine learning models with uncertainty quantification.

## Overview

The Oloren ChemEngine adapter provides access to advanced graph neural network models for molecular property prediction. It features pre-trained models for common ADMET properties, uncertainty quantification for reliability assessment, and support for batch predictions.

## Features

- **Pre-trained Models**: Ready-to-use models for 10+ common ADMET properties
- **Uncertainty Quantification**: Confidence estimates for every prediction
- **Multiple Architectures**: Support for GCN, AttentiveFP, and MPNN models
- **Batch Processing**: Efficient prediction for multiple molecules
- **Transfer Learning**: Fine-tune models on custom datasets
- **Modern ML**: State-of-the-art graph neural networks

## Installation

```bash
pip install olorenchemengine
```

## Supported Properties

### Physicochemical
- **solubility**: Aqueous solubility (log mol/L)
- **logp**: Lipophilicity (octanol-water partition coefficient)

### Absorption
- **permeability**: Membrane permeability (log cm/s)
- **caco2**: Caco-2 cell permeability (log cm/s)
- **bioavailability**: Oral bioavailability (probability)

### Distribution
- **bbb_permeability**: Blood-brain barrier permeability (log BB ratio)

### Metabolism
- **clearance**: Hepatic clearance (mL/min/kg)

### Excretion
- **half_life**: Plasma half-life (hours)

### Toxicity
- **herg**: hERG channel inhibition (pIC50)
- **ames**: Ames mutagenicity (probability)

## Usage

### Basic Usage

```python
from adapters.olorenchemengine.adapter import OlorenChemEngineAdapter

# Initialize adapter
adapter = OlorenChemEngineAdapter()

# Single molecule prediction
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin

print(result.data["predictions"][0]["properties"])
# {
#   "solubility": {"value": -2.5, "unit": "log(mol/L)", "uncertainty": 0.3},
#   "logp": {"value": 1.2, "unit": "unitless", "uncertainty": 0.2},
#   ...
# }
```

### Batch Predictions

```python
# Multiple molecules
smiles_list = [
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    "CCO",                     # Ethanol
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
]

result = await adapter.execute(smiles_list)

for pred in result.data["predictions"]:
    print(f"SMILES: {pred['smiles']}")
    print(f"Solubility: {pred['properties']['solubility']}")
```

### Specific Properties

```python
# Predict only specific properties
result = await adapter.execute(
    "CC(=O)Oc1ccccc1C(=O)O",
    properties=["solubility", "logp", "permeability"]
)
```

### Custom Model Architecture

```python
# Use specific model architecture
adapter = OlorenChemEngineAdapter(
    config={
        "model": "attentivefp",  # or "gcn", "mpnn"
        "include_uncertainty": True,
        "batch_size": 64
    }
)

result = await adapter.execute("CCO")
```

### Without Uncertainty

```python
# Faster predictions without uncertainty quantification
result = await adapter.execute(
    "CC(=O)Oc1ccccc1C(=O)O",
    include_uncertainty=False
)
```

## Integration with PharmForge Pipeline

```python
from backend.core.pipeline import Pipeline

# Create pipeline with Oloren predictions
pipeline = Pipeline()
pipeline.add_step("property_prediction", "olorenchemengine", {
    "properties": ["solubility", "permeability", "herg"],
    "include_uncertainty": True
})

# Execute pipeline
results = await pipeline.execute("CC(=O)Oc1ccccc1C(=O)O")
```

## Output Format

### Successful Prediction

```json
{
  "predictions": [
    {
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "properties": {
        "solubility": {
          "value": -2.5,
          "unit": "log(mol/L)",
          "uncertainty": 0.3,
          "confidence": 0.7
        },
        "logp": {
          "value": 1.2,
          "unit": "unitless",
          "uncertainty": 0.2,
          "confidence": 0.8
        }
      }
    }
  ],
  "num_molecules": 1,
  "properties_predicted": ["solubility", "logp"],
  "model_used": "default",
  "uncertainty_included": true,
  "metadata": {
    "model_version": "Oloren ChemEngine v1.0",
    "architecture": "default"
  }
}
```

## Property Categories

Access properties by category:

```python
adapter = OlorenChemEngineAdapter()
categories = adapter.get_property_categories()

# {
#   "physicochemical": ["solubility", "logp"],
#   "absorption": ["permeability", "caco2", "bioavailability"],
#   "distribution": ["bbb_permeability"],
#   "metabolism": ["clearance"],
#   "excretion": ["half_life"],
#   "toxicity": ["herg", "ames"]
# }
```

## Comparison with Other Adapters

| Feature | Oloren ChemEngine | ADMET-ai | Chemprop |
|---------|-------------------|----------|----------|
| Pre-trained models | ✅ Yes | ✅ Yes | ⚠️ Need training |
| Uncertainty quantification | ✅ Yes | ❌ No | ⚠️ Limited |
| Model architectures | ✅ Multiple | ⚠️ Single | ✅ Multiple |
| Transfer learning | ✅ Easy | ❌ No | ✅ Yes |
| Property coverage | ⚠️ 10+ | ✅ 49 | ✅ Custom |
| Batch efficiency | ✅ Optimized | ✅ Good | ✅ Good |

## Use Cases

### 1. Quick ADMET Screening
```python
# Screen library for drug-like properties
compounds = get_compound_library()
result = await adapter.execute(
    compounds,
    properties=["solubility", "permeability", "bioavailability"]
)
```

### 2. Risk Assessment with Uncertainty
```python
# Identify high-confidence predictions
result = await adapter.execute(smiles, include_uncertainty=True)
for pred in result.data["predictions"]:
    for prop, value in pred["properties"].items():
        if value["confidence"] > 0.8:
            print(f"High confidence: {prop} = {value['value']}")
```

### 3. Ensemble with ADMET-ai
```python
# Combine predictions from multiple models
from adapters.admet_ai.adapter import ADMETaiAdapter

oce_result = await oce_adapter.execute(smiles)
admet_result = await admet_adapter.execute(smiles)

# Average predictions where properties overlap
```

### 4. Early-Stage Filtering
```python
# Fast filtering before expensive docking
compounds = ["CCO", "CC(=O)O", ...]
result = await adapter.execute(
    compounds,
    properties=["herg", "ames"],  # Toxicity only
    include_uncertainty=False  # Faster
)

safe_compounds = [
    pred["smiles"] for pred in result.data["predictions"]
    if pred["properties"]["herg"]["value"] < 5.5
]
```

## Advanced Configuration

### Custom Property Defaults

```python
adapter = OlorenChemEngineAdapter(
    config={
        "properties": ["solubility", "logp", "permeability"],
        "model": "attentivefp",
        "include_uncertainty": True,
        "batch_size": 128
    }
)
```

### Model Selection by Property

Different models may perform better for different properties:

- **GCN**: Fast, good for physicochemical properties
- **AttentiveFP**: Best accuracy for ADMET properties
- **MPNN**: Balanced performance, good generalization

## Error Handling

```python
result = await adapter.execute("invalid_smiles")
if not result.success:
    print(f"Error: {result.error}")
    # Output: "Invalid SMILES string(s)"

# Check for unsupported properties
result = await adapter.execute(
    "CCO",
    properties=["unknown_property"]
)
# Returns error with list of supported properties
```

## Performance Considerations

### Batch Size Optimization
- Small molecules (<30 atoms): batch_size=64
- Large molecules (>30 atoms): batch_size=32
- Memory constrained: batch_size=16

### Caching
Predictions are automatically cached based on:
- SMILES (canonicalized)
- Properties requested
- Model architecture
- Uncertainty setting

## Limitations

1. **Property Coverage**: Limited to 10+ pre-trained properties (vs ADMET-ai's 49)
2. **Model Size**: Larger models require more memory
3. **Speed**: Uncertainty quantification adds ~20% overhead
4. **Dependencies**: Requires PyTorch and graph neural network libraries

## Future Enhancements

- [ ] Custom model training interface
- [ ] Active learning for targeted predictions
- [ ] Multi-task learning across properties
- [ ] Integration with molecular optimization
- [ ] Explainability features (attention visualization)

## References

- Oloren ChemEngine GitHub: https://github.com/Oloren-AI/olorenchemengine
- Paper: [Link to publication when available]
- Documentation: https://docs.oloren.ai/

## Support

For issues specific to this adapter:
- Open an issue in the PharmForge repository

For Oloren ChemEngine issues:
- Visit: https://github.com/Oloren-AI/olorenchemengine/issues

## License

This adapter: MIT License (PharmForge)
Oloren ChemEngine: Check upstream repository for license terms
