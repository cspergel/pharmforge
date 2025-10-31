# Molecular Docking and Binding Prediction Adapters

Complete implementation of molecular docking and binding affinity prediction adapters for PharmForge.

## Overview

This collection provides four specialized adapters for molecular docking and binding prediction workflows:

| Adapter | Type | Purpose | Speed | Accuracy |
|---------|------|---------|-------|----------|
| **DiffDock** | Local | Blind docking with diffusion models | Medium | High |
| **GNINA** | Local | CNN-enhanced docking | Fast | Good |
| **BindingDB** | API | Experimental binding data | Fast | N/A (Reference) |
| **SwissTarget** | API | Target prediction | Medium | Good |

## Quick Start

```python
import asyncio
from adapters.gnina.adapter import GNINAAdapter

async def dock_molecule():
    # Initialize adapter
    gnina = GNINAAdapter(config={
        'receptor_path': 'protein.pdbqt',
        'center_x': 10.0, 'center_y': 20.0, 'center_z': 30.0
    })

    # Run docking
    result = await gnina.execute("CC(=O)Oc1ccccc1C(=O)O")

    if result.success:
        print(f"Binding affinity: {result.data['binding_affinity']:.2f} kcal/mol")
        print(f"Confidence score: {result.data['binding_score']:.3f}")

asyncio.run(dock_molecule())
```

## Directory Structure

```
adapters/
├── diffdock/
│   ├── __init__.py
│   └── adapter.py          # DiffDock diffusion-based docking
├── gnina/
│   ├── __init__.py
│   └── adapter.py          # GNINA CNN-enhanced docking
├── bindingdb/
│   ├── __init__.py
│   └── adapter.py          # BindingDB experimental data
├── swisstarget/
│   ├── __init__.py
│   └── adapter.py          # SwissTargetPrediction
└── DOCKING_README.md       # This file
```

## Adapter Details

### 1. DiffDock Adapter

**Diffusion-based blind molecular docking**

- **No binding site required**: Automatically finds binding poses
- **State-of-the-art accuracy**: ~37% success rate (RMSD < 2Å)
- **Confidence scoring**: Provides pose confidence estimates
- **GPU accelerated**: Requires GPU for practical use

**Installation:**
```bash
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
conda env create -f environment.yml
```

**Usage:**
```python
from adapters.diffdock.adapter import DiffDockAdapter

diffdock = DiffDockAdapter(config={
    'diffdock_path': '/path/to/DiffDock',
    'inference_steps': 20,
    'samples_per_complex': 40,
    'use_gpu': True
})

result = await diffdock.execute({
    'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'protein_path': 'protein.pdb'
})

print(f"Confidence: {result.data['confidence']:.3f}")
print(f"Estimated affinity: {result.data['estimated_affinity']:.2f} kcal/mol")
```

**Key Features:**
- Blind docking (no site specification)
- Multiple pose predictions
- Confidence scores
- High accuracy on benchmarks

### 2. GNINA Adapter

**CNN-enhanced AutoDock Vina**

- **Fast docking**: Similar speed to Vina (~10-30s per ligand)
- **Improved scoring**: CNN-based affinity prediction
- **Dual scores**: Both Vina and CNN scores
- **GPU optional**: Works on CPU, faster with GPU

**Installation:**
```bash
# Download binary
wget https://github.com/gnina/gnina/releases/download/v1.0.3/gnina
chmod +x gnina
sudo mv gnina /usr/local/bin/
```

**Usage:**
```python
from adapters.gnina.adapter import GNINAAdapter

gnina = GNINAAdapter(config={
    'receptor_path': 'receptor.pdbqt',
    'center_x': 10.0,
    'center_y': 20.0,
    'center_z': 30.0,
    'size_x': 25,
    'size_y': 25,
    'size_z': 25,
    'exhaustiveness': 8,
    'cnn_scoring': 'default'
})

result = await gnina.execute("CC(C)Cc1ccc(cc1)C(C)C(=O)O")

print(f"Vina score: {result.data['vina_affinity']:.2f} kcal/mol")
print(f"CNN score: {result.data['cnn_affinity']:.2f} kcal/mol")
```

**Key Features:**
- CNN-based scoring
- Vina-compatible workflow
- Fast execution
- Multiple scoring functions

### 3. BindingDB Adapter

**Experimental binding affinity database**

- **2.5M+ measurements**: Ki, Kd, IC50, EC50 values
- **Literature references**: PubMed citations
- **Assay conditions**: pH, temperature, organism
- **Validation data**: Compare with computational predictions

**Usage:**
```python
from adapters.bindingdb.adapter import BindingDBAdapter

bindingdb = BindingDBAdapter(config={
    'max_results': 100,
    'timeout': 30
})

result = await bindingdb.execute({
    'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
    'target': 'COX-2'
})

for entry in result.data['entries']:
    print(f"{entry['type']}: {entry['value']} {entry['unit']}")
    print(f"  Target: {entry['target']}")
    print(f"  Reference: {entry['reference']}")
```

**Key Features:**
- Experimental binding data
- Multiple measurement types
- Assay metadata
- Literature references

### 4. SwissTargetPrediction Adapter

**Target prediction for small molecules**

- **Target discovery**: Find potential protein targets
- **Probability scores**: Confidence estimates
- **Multiple organisms**: Human, mouse, rat, etc.
- **Target classification**: Enzyme, GPCR, kinase, etc.

**Usage:**
```python
from adapters.swisstarget.adapter import SwissTargetAdapter

swisstarget = SwissTargetAdapter(config={
    'organism': 'Homo sapiens',
    'min_probability': 0.5,
    'max_targets': 15
})

result = await swisstarget.execute("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

for target in result.data['targets']:
    print(f"{target['target_name']}: {target['probability']:.2f}")
    print(f"  Class: {target['target_class']}")
    print(f"  Known actives: {target['known_actives']}")
```

**Key Features:**
- Target prediction
- Probability scoring
- Multi-organism support
- Target classification

## Workflow Examples

### Complete Drug Discovery Workflow

```python
async def drug_discovery_workflow(smiles: str):
    """
    Complete workflow: target prediction → experimental data → docking
    """
    # Step 1: Predict potential targets
    swiss = SwissTargetAdapter()
    targets = await swiss.execute(smiles)

    top_target = targets.data['targets'][0]
    print(f"Top target: {top_target['target_name']}")

    # Step 2: Query experimental binding data
    bindingdb = BindingDBAdapter()
    exp_data = await bindingdb.execute({
        'smiles': smiles,
        'target': top_target['target_name']
    })

    if exp_data.data['count'] > 0:
        ki_values = [e for e in exp_data.data['entries'] if e['type'] == 'Ki']
        if ki_values:
            print(f"Experimental Ki: {ki_values[0]['value']} {ki_values[0]['unit']}")

    # Step 3: Computational docking (requires receptor)
    # gnina = GNINAAdapter(config={'receptor_path': '...'})
    # docking = await gnina.execute(smiles)

    # Step 4: Compare experimental vs computational
    # ...
```

### Batch Docking

```python
async def batch_dock(compounds: dict, receptor_path: str):
    """Dock multiple compounds efficiently"""
    gnina = GNINAAdapter(config={'receptor_path': receptor_path})

    results = []
    for name, smiles in compounds.items():
        result = await gnina.execute(smiles)
        if result.success:
            results.append({
                'name': name,
                'smiles': smiles,
                'affinity': result.data['binding_affinity'],
                'score': result.data['binding_score']
            })

    # Rank by affinity
    ranked = sorted(results, key=lambda x: x['affinity'])
    return ranked
```

### Virtual Screening

```python
async def virtual_screening(
    library: List[str],
    receptor_path: str,
    threshold: float = -7.0
):
    """Screen compound library and filter by binding affinity"""
    gnina = GNINAAdapter(config={'receptor_path': receptor_path})

    hits = []
    for smiles in library:
        result = await gnina.execute(smiles)

        if result.success:
            affinity = result.data['binding_affinity']
            if affinity <= threshold:  # More negative = better
                hits.append({
                    'smiles': smiles,
                    'affinity': affinity,
                    'score': result.data['binding_score']
                })

    return sorted(hits, key=lambda x: x['affinity'])
```

## Integration with Other Adapters

### With Structure Preparation (OpenMM)

```python
from adapters.openmm.adapter import OpenMMAdapter

# 1. Prepare and minimize protein
openmm = OpenMMAdapter()
prepared = await openmm.execute('protein.pdb', task='minimize')

# 2. Dock ligand
gnina = GNINAAdapter(config={'receptor_path': prepared.data['output_path']})
result = await gnina.execute(smiles)
```

### With Retrosynthesis Planning

```python
from adapters.aizynthfinder.adapter import AiZynthFinderAdapter

# 1. Dock and evaluate
docking = await gnina.execute(smiles)

# 2. Plan synthesis for promising hits
if docking.data['binding_score'] > 0.7:
    aizf = AiZynthFinderAdapter()
    synthesis = await aizf.execute(smiles)
    print(f"Synthesis steps: {synthesis.data['num_steps']}")
```

## Performance Guidelines

### Speed Comparison

| Adapter | Time/Ligand | Best For |
|---------|-------------|----------|
| GNINA | 10-30s | Large screens |
| DiffDock | 30-60s (GPU) | High accuracy |
| BindingDB | 1-5s | Data lookup |
| SwissTarget | 30-60s | Target discovery |

### Resource Requirements

| Adapter | CPU | GPU | Memory |
|---------|-----|-----|--------|
| GNINA | ✓ | Optional | 2-4 GB |
| DiffDock | ✗ | Required | 4-8 GB |
| BindingDB | ✓ | ✗ | < 1 GB |
| SwissTarget | ✓ | ✗ | < 1 GB |

### Accuracy Benchmarks

**Pose Prediction (RMSD < 2Å on PDBbind)**
- DiffDock: ~37% success rate
- GNINA: ~25% success rate
- AutoDock Vina: ~20% success rate

**Affinity Prediction (Correlation with experimental)**
- DiffDock: R² ~0.82
- GNINA: R² ~0.65
- AutoDock Vina: R² ~0.52

## Testing

Run the test suite:

```bash
# All docking adapter tests
pytest backend/tests/test_docking_adapters.py -v

# Specific adapter
pytest backend/tests/test_docking_adapters.py::TestGNINAAdapter -v

# Integration tests
pytest backend/tests/test_docking_adapters.py::TestIntegration -v
```

Run benchmarks:

```bash
# Performance benchmarks
python benchmarks/docking_benchmarks.py

# With custom receptor
python benchmarks/docking_benchmarks.py --receptor /path/to/receptor.pdbqt
```

## Troubleshooting

### Common Issues

**1. DiffDock not found**
```python
# Ensure DiffDock path is correct
config = {'diffdock_path': '/absolute/path/to/DiffDock'}
```

**2. GNINA binary not found**
```bash
# Add to PATH
export PATH="/path/to/gnina:$PATH"

# Or specify full path
config = {'gnina_binary': '/usr/local/bin/gnina'}
```

**3. Receptor preparation**
```bash
# Convert PDB to PDBQT
prepare_receptor4.py -r protein.pdb -o protein.pdbqt
```

**4. BindingDB timeout**
```python
# Increase timeout for slow connections
config = {'timeout': 60}
```

## Best Practices

### 1. Receptor Preparation

- Remove water molecules
- Add hydrogens
- Assign partial charges
- Minimize structure

### 2. Binding Site Definition

- Use crystal structure ligand center
- Ensure box size covers binding site
- Allow ~5Å margin around ligand

### 3. Validation

- Compare with experimental data
- Test on known actives/inactives
- Cross-validate with multiple methods

### 4. Virtual Screening Strategy

1. Fast filtering (GNINA)
2. Accurate re-docking (DiffDock)
3. Experimental validation (BindingDB)
4. Synthesis planning (AiZynthFinder)

## References

### Papers

1. **DiffDock**: Corso, G., et al. (2023). "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking." ICLR 2023.
   - https://arxiv.org/abs/2210.01776

2. **GNINA**: McNutt, A.T., et al. (2021). "GNINA 1.0: molecular docking with deep learning." Journal of Cheminformatics, 13(1), 43.
   - https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00522-2

3. **BindingDB**: Gilson, M.K., et al. (2016). "BindingDB in 2015: A public database for medicinal chemistry, computational chemistry and systems pharmacology." Nucleic Acids Research, 44(D1), D1045-D1053.

4. **SwissTargetPrediction**: Daina, A., et al. (2019). "SwissTargetPrediction: updated data and new features for efficient prediction of protein targets of small molecules." Nucleic Acids Research, 47(W1), W357-W364.

### Documentation Links

- DiffDock GitHub: https://github.com/gcorso/DiffDock
- GNINA GitHub: https://github.com/gnina/gnina
- BindingDB: https://www.bindingdb.org/
- SwissTargetPrediction: http://www.swisstargetprediction.ch/

## Contributing

See main PharmForge documentation for contribution guidelines.

## License

See LICENSE file in repository root.

## Support

- Issues: https://github.com/your-repo/PharmForge/issues
- Documentation: https://pharmforge.readthedocs.io/
- Email: support@pharmforge.org
