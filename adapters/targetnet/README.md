# TargetNet Adapter

Deep learning-based protein target prediction for PharmForge.

## Overview

TargetNet predicts protein targets for small molecules using deep learning on molecular fingerprints. This adapter enables target identification, off-target analysis, and polypharmacology profiling.

## Features

- **Target Prediction**: Predict protein targets for molecules
- **Off-Target Analysis**: Identify potential off-target interactions
- **Polypharmacology**: Analyze multi-target profiles
- **Confidence Scoring**: Probabilistic target predictions
- **Target Database**: Integration with ChEMBL/UniProt

## Installation

```bash
# Install dependencies
pip install torch rdkit numpy

# Optional: Install target prediction models
# pip install targetnet  # If available
```

## Usage

### Basic Target Prediction

```python
from adapters.targetnet import TargetNetAdapter

# Initialize adapter
adapter = TargetNetAdapter(config={
    'confidence_threshold': 0.5,
    'max_targets': 10
})

# Predict targets
smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen

result = await adapter.execute(
    input_data=smiles,
    max_targets=10
)

primary_targets = result.data['primary_targets']
off_targets = result.data['off_targets']
```

### Off-Target Analysis

```python
# Include off-target predictions
result = await adapter.execute(
    input_data=smiles,
    include_off_targets=True,
    confidence_threshold=0.4  # Lower threshold for off-targets
)

for target in result.data['off_targets']:
    print(f"{target['target_name']}: {target['confidence']:.2f}")
```

### Polypharmacology Profiling

```python
# Analyze polypharmacology
result = await adapter.execute(
    input_data=smiles,
    max_targets=20
)

poly = result.data['polypharmacology']
print(f"Polypharmacological: {poly['is_polypharmacological']}")
print(f"Target diversity: {poly['target_diversity']:.2f}")
print(f"Selectivity: {poly['selectivity_score']:.2f}")
```

### Batch Prediction

```python
# Predict for multiple molecules
smiles_list = [
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
]

results = await adapter.predict_batch(
    smiles_list,
    max_targets=5
)
```

## Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_path` | str | None | Path to pretrained model |
| `confidence_threshold` | float | 0.5 | Minimum prediction confidence |
| `max_targets` | int | 10 | Maximum targets to return |
| `use_local` | bool | True | Use local model vs API |

## Parameters

### execute() Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_data` | str | - | SMILES string or fingerprint |
| `confidence_threshold` | float | 0.5 | Minimum confidence |
| `max_targets` | int | 10 | Maximum targets |
| `target_types` | list | None | Filter by type (['Kinase', 'GPCR']) |
| `include_off_targets` | bool | True | Include off-target predictions |

## Output Format

```python
{
    'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
    'primary_targets': [
        {
            'target_id': 'CHEMBL1862',
            'target_name': 'Cyclooxygenase-2',
            'uniprot_id': 'P35354',
            'target_type': 'Enzyme',
            'confidence': 0.89,
            'disease_areas': ['Inflammation', 'Pain']
        },
        ...
    ],
    'off_targets': [
        {
            'target_id': 'CHEMBL217',
            'target_name': 'Acetylcholinesterase',
            'confidence': 0.58,
            ...
        },
        ...
    ],
    'polypharmacology': {
        'is_polypharmacological': True,
        'target_count': 5,
        'target_diversity': 0.75,
        'selectivity_score': 0.6
    },
    'statistics': {
        'total_targets': 8,
        'primary_targets': 3,
        'off_targets': 5,
        'mean_confidence': 0.68
    }
}
```

## Target Database

The adapter includes a curated target database with:

- **ChEMBL IDs**: Standard target identifiers
- **UniProt IDs**: Protein database cross-references
- **Target Types**: Kinase, GPCR, Enzyme, Ion Channel, etc.
- **Disease Areas**: Associated therapeutic areas

### Included Targets

- Dopamine D2 receptor
- Acetylcholinesterase
- Cyclooxygenase-2 (COX-2)
- VEGFR2
- EGFR
- Beta-2 adrenergic receptor
- Monoamine oxidase A
- CYP3A4
- Serotonin 2a receptor
- hERG

## Confidence Thresholds

Recommended thresholds:

- **High confidence**: > 0.7 (primary targets)
- **Medium confidence**: 0.5-0.7 (likely targets)
- **Low confidence**: 0.3-0.5 (off-targets)
- **Very low**: < 0.3 (exclude)

## Target Types

Supported target types:

- **Kinase**: Protein kinases
- **GPCR**: G-protein coupled receptors
- **Enzyme**: Enzymes (non-kinase)
- **Ion Channel**: Ion channels
- **Transporter**: Membrane transporters
- **Nuclear Receptor**: Nuclear receptors

## Polypharmacology Analysis

Metrics:

- **is_polypharmacological**: Multiple high-confidence targets
- **target_count**: Number of predicted targets
- **target_diversity**: Target type diversity (0-1)
- **selectivity_score**: Inverse of promiscuity (0-1)

## Use Cases

### 1. Drug Repurposing

```python
# Find new uses for existing drugs
known_drug = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

result = await adapter.execute(
    input_data=known_drug,
    max_targets=20,
    confidence_threshold=0.4
)

# Look for unexpected high-confidence targets
```

### 2. Safety Profiling

```python
# Identify potential off-target liabilities
candidate = "NEW_DRUG_SMILES"

result = await adapter.execute(
    input_data=candidate,
    include_off_targets=True
)

# Check for problematic targets (e.g., hERG)
herg_hits = [
    t for t in result.data['all_predictions']
    if 'hERG' in t['target_name']
]
```

### 3. Multi-Target Drug Design

```python
# Design polypharmacological drugs
result = await adapter.execute(
    input_data=molecule,
    max_targets=15
)

if result.data['polypharmacology']['is_polypharmacological']:
    print("Suitable for multi-target approach")
```

## Integration

### With De Novo Generation

```python
from adapters.denovo import DeNovoAdapter
from adapters.targetnet import TargetNetAdapter

# Generate molecules
denovo = DeNovoAdapter()
gen_result = await denovo.execute(input_data=None, num_molecules=10)

# Predict targets for generated molecules
target_adapter = TargetNetAdapter()

for mol in gen_result.data['molecules']:
    target_result = await target_adapter.execute(mol['smiles'])
    # Filter molecules with desired targets
```

### With ADMET Prediction

```python
from adapters.admet_ai import ADMETaiAdapter
from adapters.targetnet import TargetNetAdapter

# Predict both targets and ADMET
admet = ADMETaiAdapter()
target = TargetNetAdapter()

admet_result = await admet.execute(smiles)
target_result = await target.execute(smiles)

# Combine for comprehensive profiling
```

## Performance

- **Prediction time**: 0.1-0.5s per molecule
- **Accuracy**: Depends on model and target
- **Coverage**: ~100+ protein targets
- **Batch processing**: Efficient for multiple molecules

## Advantages

1. **Fast predictions**: No experimental assays needed
2. **Comprehensive**: Multiple targets simultaneously
3. **Early identification**: Detect issues early in discovery
4. **Drug repurposing**: Find new uses for compounds

## Limitations

1. **Model dependent**: Accuracy varies by target
2. **Bias**: Limited by training data
3. **Novel targets**: May miss new/orphan targets
4. **Confidence calibration**: Probabilities may not be perfectly calibrated

## References

- Target prediction methods: ChEMBL, BindingDB
- Molecular fingerprints: Morgan/ECFP fingerprints
- Deep learning: Neural network architectures for target prediction

## Troubleshooting

### Low confidence predictions

- Lower confidence_threshold
- Increase max_targets
- Check molecule validity
- Verify molecular fingerprint

### Missing targets

- Use lower threshold for off-targets
- Check target_types filter
- Verify target in database
- Consider batch prediction

### Performance issues

- Reduce max_targets
- Use batch processing
- Enable caching
- Consider API mode

## Future Enhancements

- Expanded target database
- Structure-based prediction
- Binding site prediction
- Target-specific models
- Activity prediction

## License

Adapter code: MIT License
Target databases: See respective licenses (ChEMBL, UniProt)
