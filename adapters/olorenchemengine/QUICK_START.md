# Oloren ChemEngine Adapter - Quick Start Guide

Get up and running with the Oloren ChemEngine adapter in 5 minutes.

## Installation

```bash
# Install Oloren ChemEngine
pip install olorenchemengine

# Verify installation
python -c "import olorenchemengine; print('Oloren ChemEngine installed successfully')"
```

## Basic Usage

### Single Molecule Prediction

```python
import asyncio
from adapters.olorenchemengine.adapter import OlorenChemEngineAdapter

async def predict_properties():
    adapter = OlorenChemEngineAdapter()

    # Predict properties for aspirin
    result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

    if result.success:
        prediction = result.data["predictions"][0]
        for prop, values in prediction["properties"].items():
            print(f"{prop}: {values['value']:.3f} {values['unit']}")

asyncio.run(predict_properties())
```

### Batch Predictions

```python
async def batch_predict():
    adapter = OlorenChemEngineAdapter()

    compounds = [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CCO"  # Ethanol
    ]

    result = await adapter.execute(compounds)

    print(f"Predicted {result.data['num_molecules']} molecules")
    for pred in result.data["predictions"]:
        print(f"\nSMILES: {pred['smiles']}")
        print(f"Solubility: {pred['properties']['solubility']['value']:.3f}")

asyncio.run(batch_predict())
```

### Specific Properties Only

```python
async def specific_properties():
    adapter = OlorenChemEngineAdapter()

    result = await adapter.execute(
        "CC(=O)Oc1ccccc1C(=O)O",
        properties=["solubility", "logp", "permeability"]
    )

    # Only these 3 properties will be predicted
    print(result.data["predictions"][0]["properties"].keys())

asyncio.run(specific_properties())
```

## Common Use Cases

### 1. ADMET Screening

```python
# Quick ADMET profile
adapter = OlorenChemEngineAdapter(
    config={
        "properties": ["solubility", "permeability", "bioavailability", "herg"]
    }
)

result = await adapter.execute(your_smiles)
```

### 2. Toxicity Assessment

```python
# Check toxicity flags
result = await adapter.execute(
    your_smiles,
    properties=["herg", "ames"]
)

prediction = result.data["predictions"][0]
herg = prediction["properties"]["herg"]["value"]
ames = prediction["properties"]["ames"]["value"]

if herg > 6.0 or ames > 0.5:
    print("WARNING: Potential toxicity issues")
```

### 3. High-Confidence Filtering

```python
# Only trust high-confidence predictions
result = await adapter.execute(your_smiles, include_uncertainty=True)

for prop, values in result.data["predictions"][0]["properties"].items():
    if values["confidence"] > 0.8:
        print(f"High confidence: {prop} = {values['value']}")
```

## Configuration Options

### Custom Model Architecture

```python
# Use AttentiveFP for better ADMET accuracy
adapter = OlorenChemEngineAdapter(
    config={"model": "attentivefp"}
)
```

### Disable Uncertainty (Faster)

```python
# 20% faster without uncertainty quantification
result = await adapter.execute(
    smiles,
    include_uncertainty=False
)
```

### Batch Size Tuning

```python
# Optimize for your hardware
adapter = OlorenChemEngineAdapter(
    config={"batch_size": 64}  # Default: 32
)
```

## Property Reference

### Quick Property Lookup

```python
adapter = OlorenChemEngineAdapter()

# See all available properties
print(adapter.SUPPORTED_PROPERTIES.keys())

# Get properties by category
categories = adapter.get_property_categories()
print(categories["toxicity"])  # ['herg', 'ames']

# Get property details
info = adapter.get_property_info("solubility")
print(info["unit"])  # "log(mol/L)"
```

### Available Properties by Category

**Physicochemical**
- solubility, logp

**Absorption**
- permeability, caco2, bioavailability

**Distribution**
- bbb_permeability

**Metabolism**
- clearance

**Excretion**
- half_life

**Toxicity**
- herg, ames

## Integration with PharmForge

### Add to Pipeline

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()
pipeline.add_step("property_prediction", "olorenchemengine", {
    "properties": ["solubility", "permeability", "herg"],
    "include_uncertainty": True
})

results = await pipeline.execute("your_smiles_here")
```

### Combine with Other Adapters

```python
# Ensemble predictions
from adapters.admet_ai.adapter import ADMETaiAdapter

oce = OlorenChemEngineAdapter()
admet_ai = ADMETaiAdapter()

oce_result = await oce.execute(smiles)
admet_result = await admet_ai.execute(smiles)

# Compare or average predictions
```

## Error Handling

```python
result = await adapter.execute(smiles)

if not result.success:
    print(f"Error: {result.error}")
    # Handle error appropriately
else:
    # Process predictions
    pass
```

## Performance Tips

### 1. Use Batch Processing
```python
# ✅ Good: Single batch call
result = await adapter.execute([smiles1, smiles2, smiles3])

# ❌ Bad: Multiple single calls
for smiles in [smiles1, smiles2, smiles3]:
    result = await adapter.execute(smiles)
```

### 2. Leverage Caching
```python
# Predictions are automatically cached
# Repeated calls with same SMILES are instant
result1 = await adapter.execute(smiles)  # Computes
result2 = await adapter.execute(smiles)  # Cache hit
```

### 3. Choose Right Model
```python
# Fast, good for physicochemical properties
adapter = OlorenChemEngineAdapter(config={"model": "gcn"})

# Best accuracy for ADMET
adapter = OlorenChemEngineAdapter(config={"model": "attentivefp"})

# Balanced
adapter = OlorenChemEngineAdapter(config={"model": "default"})
```

## Troubleshooting

### Issue: "Oloren ChemEngine not installed"
```bash
pip install olorenchemengine
```

### Issue: "Invalid SMILES string"
```python
# Validate SMILES first
if adapter.validate_input(smiles):
    result = await adapter.execute(smiles)
else:
    print("Invalid SMILES")
```

### Issue: "Unsupported properties"
```python
# Check supported properties
print(list(adapter.SUPPORTED_PROPERTIES.keys()))
```

### Issue: Out of memory with large batches
```python
# Reduce batch size
adapter = OlorenChemEngineAdapter(config={"batch_size": 16})
```

## Next Steps

- **Full Documentation**: See [README.md](README.md) for complete feature list
- **Examples**: Run [example_usage.py](example_usage.py) for comprehensive demos
- **Tests**: Check [test_olorenchemengine_adapter.py](../../backend/tests/test_olorenchemengine_adapter.py)

## Quick Reference Card

```python
from adapters.olorenchemengine.adapter import OlorenChemEngineAdapter

# Initialize
adapter = OlorenChemEngineAdapter()

# Single prediction
result = await adapter.execute("CCO")

# Batch prediction
result = await adapter.execute(["CCO", "CC(=O)O"])

# Specific properties
result = await adapter.execute("CCO", properties=["solubility", "logp"])

# Without uncertainty (faster)
result = await adapter.execute("CCO", include_uncertainty=False)

# Check result
if result.success:
    predictions = result.data["predictions"]
else:
    print(result.error)
```

## Support

- PharmForge Issues: https://github.com/yourusername/pharmforge/issues
- Oloren Issues: https://github.com/Oloren-AI/olorenchemengine/issues
- Documentation: See README.md in this directory
