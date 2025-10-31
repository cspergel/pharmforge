# REINVENT Adapter

Reinforcement learning-based molecular design adapter for PharmForge.

## Overview

REINVENT (Reinforcement Learning for Structural Evolution) is a state-of-the-art approach for de novo drug design using reinforcement learning. This adapter provides access to REINVENT's capabilities for generating novel molecules, optimizing existing structures, and decorating scaffolds.

## Features

- **De Novo Generation**: Generate novel molecules with desired properties
- **Property Optimization**: Optimize molecules for multiple objectives
- **Scaffold Decoration**: Add functional groups to core scaffolds
- **Multi-Objective Design**: Balance multiple property constraints
- **Diversity Control**: Ensure structural diversity in generated sets

## Installation

### Option 1: Local REINVENT Installation

```bash
# Install REINVENT (requires specific dependencies)
pip install reinvent

# Note: REINVENT may have specific version requirements
```

### Option 2: Fallback Mode

The adapter includes a fallback mode using RDKit for basic generation when REINVENT is not available.

## Usage

### Basic Generation

```python
from adapters.reinvent import REINVENTAdapter

# Initialize adapter
adapter = REINVENTAdapter(config={
    'mode': 'generate',
    'use_local': True
})

# Define target properties
target_properties = {
    'molecular_weight': 400.0,
    'logp': 3.0,
    'qed': 0.7
}

# Generate molecules
result = await adapter.execute(
    input_data={},
    num_molecules=100,
    target_properties=target_properties,
    optimization_steps=100
)

molecules = result.data['molecules']
```

### Molecule Optimization

```python
# Optimize existing molecule
adapter = REINVENTAdapter(config={'mode': 'optimize'})

seed_smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

result = await adapter.execute(
    input_data=seed_smiles,
    num_molecules=50,
    target_properties={'qed': 0.8}
)
```

### Scaffold Decoration

```python
# Decorate a scaffold
adapter = REINVENTAdapter(config={'mode': 'decorate'})

scaffold = "c1ccccc1"  # Benzene

result = await adapter.execute(
    input_data={},
    scaffold=scaffold,
    num_molecules=50,
    target_properties={'qed': 0.7}
)
```

## Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | str | 'generate' | Operation mode: 'generate', 'optimize', 'decorate' |
| `model_path` | str | None | Path to pretrained REINVENT model |
| `use_local` | bool | True | Use local installation vs API |
| `api_url` | str | None | URL for REINVENT API service |

## Parameters

### execute() Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_data` | dict/str | - | Target properties dict or seed SMILES |
| `num_molecules` | int | 100 | Number of molecules to generate |
| `target_properties` | dict | {} | Target property values |
| `constraints` | dict | {} | Molecular constraints |
| `optimization_steps` | int | 100 | Number of RL optimization steps |
| `diversity_filter` | bool | True | Apply diversity filtering |

## Output Format

```python
{
    'molecules': [
        {
            'smiles': 'CCO...',
            'score': 0.85,
            'properties': {
                'molecular_weight': 398.5,
                'logp': 2.95,
                'qed': 0.72,
                ...
            }
        },
        ...
    ],
    'statistics': {
        'total_molecules': 100,
        'valid_ratio': 0.95,
        'mean_score': 0.78,
        ...
    },
    'target_properties': {...},
    'model': 'REINVENT'
}
```

## Property Targeting

Supported target properties:

- `molecular_weight` (MW): Target molecular weight
- `logp`: Target partition coefficient
- `qed`: Target quantitative estimate of drug-likeness
- `tpsa`: Target topological polar surface area

Properties are scored using Gaussian functions centered at target values.

## Constraints

Supported constraints:

- `mw_range`: (min, max) molecular weight range
- `logp_range`: (min, max) LogP range
- `max_atoms`: Maximum number of atoms

## Performance

- Generation time: ~0.5-2 seconds per molecule (local)
- Validity rate: Typically 85-95%
- Diversity: Configurable via diversity_filter parameter

## Examples

See `example_workflows.py` for complete examples:

1. Basic de novo generation
2. Property-guided optimization
3. Scaffold hopping
4. Multi-objective optimization
5. Integration with other adapters

## References

- REINVENT 4.0: [GitHub](https://github.com/MolecularAI/REINVENT4)
- Original paper: Molecular De Novo Design through Deep Reinforcement Learning

## Troubleshooting

### REINVENT not installed

If REINVENT is not installed, the adapter will automatically use RDKit-based fallback generation.

### Low validity scores

- Adjust optimization_steps (increase for better results)
- Relax constraints
- Use diversity_filter=False for more exploration

### Memory issues

- Reduce num_molecules
- Use batch processing
- Consider API mode instead of local

## Integration

Works seamlessly with:

- **ADMET-AI**: Property prediction for generated molecules
- **TargetNet**: Target prediction for generated molecules
- **AiZynthFinder**: Retrosynthetic analysis
- **Vina**: Molecular docking

## License

Adapter code: MIT License
REINVENT software: See REINVENT repository for licensing
