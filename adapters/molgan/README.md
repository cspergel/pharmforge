# MolGAN Adapter

Generative Adversarial Network for molecular generation in PharmForge.

## Overview

MolGAN (Molecular Generative Adversarial Network) uses adversarial training to generate novel drug-like molecules. This adapter provides access to GAN-based molecular generation with property conditioning and drug-likeness optimization.

## Features

- **GAN-based Generation**: Generate molecules using adversarial training
- **Property Conditioning**: Guide generation with target properties
- **Molecular Graphs**: Direct molecular graph generation
- **Drug-likeness**: Built-in drug-likeness optimization
- **Validity Filtering**: Automatic validation of generated structures

## Installation

### Option 1: Local MolGAN Installation

```bash
# Install MolGAN dependencies
pip install torch
pip install molgan  # If available

# Or install from source
git clone https://github.com/nicola-decao/MolGAN
cd MolGAN
pip install -e .
```

### Option 2: Fallback Mode

The adapter includes RDKit-based fallback generation when MolGAN is not available.

## Usage

### Basic Generation

```python
from adapters.molgan import MolGANAdapter

# Initialize adapter
adapter = MolGANAdapter(config={
    'temperature': 0.8,
    'max_atoms': 38
})

# Generate molecules
result = await adapter.execute(
    input_data=None,
    num_molecules=100,
    validate=True
)

molecules = result.data['molecules']
statistics = result.data['statistics']
```

### Property-Guided Generation

```python
# Target CNS drug-like properties
target_properties = {
    'molecular_weight': 350.0,
    'logp': 2.5,
    'tpsa': 60.0,  # CNS penetration
    'qed': 0.75
}

result = await adapter.execute(
    input_data=None,
    num_molecules=100,
    target_properties=target_properties,
    temperature=0.8
)
```

### Conditional Generation

```python
# Generate variations of seed molecule
seed_smiles = "CC(=O)Oc1ccccc1C(=O)O"

result = await adapter.execute(
    input_data=seed_smiles,
    num_molecules=50,
    temperature=0.7
)
```

## Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model_path` | str | None | Path to pretrained model |
| `use_local` | bool | True | Use local model vs API |
| `temperature` | float | 0.8 | Sampling temperature |
| `max_atoms` | int | 38 | Maximum atoms per molecule |

## Parameters

### execute() Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_data` | dict/str/None | - | Generation parameters or seed molecule |
| `num_molecules` | int | 100 | Number to generate |
| `target_properties` | dict | {} | Target property values |
| `seed_molecules` | list | [] | Seed SMILES for conditioning |
| `temperature` | float | 0.8 | Sampling temperature (overrides config) |
| `validate` | bool | True | Validate generated molecules |

## Output Format

```python
{
    'molecules': [
        {
            'smiles': 'CCO...',
            'valid': True,
            'properties': {
                'molecular_weight': 345.2,
                'qed': 0.68,
                ...
            },
            'validity_score': 0.85,
            'source': 'MolGAN'
        },
        ...
    ],
    'statistics': {
        'total_generated': 100,
        'valid_molecules': 87,
        'validity_ratio': 0.87,
        'uniqueness_ratio': 0.93,
        'mean_qed': 0.65
    },
    'parameters': {...},
    'model': 'MolGAN'
}
```

## Temperature Parameter

The temperature parameter controls generation diversity:

- **Low (0.3-0.5)**: More conservative, higher validity
- **Medium (0.6-0.8)**: Balanced diversity and validity
- **High (0.9-1.2)**: More diverse, may reduce validity

## Validation

The adapter performs multiple validation checks:

1. **Structural validity**: Valid SMILES and RDKit molecule
2. **Atom count**: Within max_atoms limit
3. **Valence**: Proper valence for all atoms
4. **Sanitization**: Passes RDKit sanitization

## Performance Metrics

Typical performance:

- **Validity**: 70-90%
- **Uniqueness**: 85-95%
- **Novelty**: 90-95%
- **Generation time**: 0.1-0.5s per molecule

## Examples

### Example 1: Basic Drug-like Generation

```python
adapter = MolGANAdapter()

result = await adapter.execute(
    input_data=None,
    num_molecules=100,
    target_properties={'qed': 0.7}
)

# Filter for high QED
druglike = [
    m for m in result.data['molecules']
    if m['properties']['qed'] > 0.6
]
```

### Example 2: CNS Drug Generation

```python
# Target CNS properties
cns_properties = {
    'molecular_weight': 350.0,
    'logp': 2.5,
    'tpsa': 60.0,
    'hbd': 1.0
}

result = await adapter.execute(
    input_data=None,
    num_molecules=50,
    target_properties=cns_properties
)
```

### Example 3: Fragment Elaboration

```python
# Start with fragment
fragment = "c1ccccc1"  # Benzene

result = await adapter.execute(
    input_data=fragment,
    num_molecules=30,
    temperature=0.9  # Higher diversity
)
```

## Integration

Works with:

- **ADMET-AI**: Predict ADMET properties
- **TargetNet**: Predict protein targets
- **RDKit**: Property calculation
- **AiZynthfinder**: Synthetic accessibility

## Molecular Graph Representation

MolGAN generates molecules as molecular graphs:

- **Nodes**: Atoms (C, N, O, S, F, Cl, Br)
- **Edges**: Bonds (single, double, triple)
- **Max size**: Configurable via max_atoms

## Advantages

1. **Direct graph generation**: No SMILES grammar constraints
2. **Property conditioning**: Guide generation with targets
3. **Adversarial training**: Learn realistic molecular distributions
4. **Fast sampling**: Quick generation from latent space

## Limitations

1. **Validity**: Lower than grammar-based methods
2. **Large molecules**: Limited by max_atoms parameter
3. **Training required**: Pre-trained models needed
4. **Memory**: Can be memory-intensive

## References

- MolGAN paper: [arXiv:1805.11973](https://arxiv.org/abs/1805.11973)
- GitHub: [MolGAN Implementation](https://github.com/nicola-decao/MolGAN)

## Troubleshooting

### Low validity rates

- Decrease temperature
- Increase validation threshold
- Use pre-trained models
- Check max_atoms parameter

### Poor property targeting

- Adjust target_properties values
- Generate more molecules
- Use property optimization post-processing

### Memory issues

- Reduce num_molecules
- Lower max_atoms
- Use batch generation
- Clear GPU memory between runs

## License

Adapter code: MIT License
MolGAN: See original repository for licensing
