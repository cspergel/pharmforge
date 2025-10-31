# xTB Adapter - Quick Start Guide

## Installation

```bash
# Install via conda (REQUIRED - not available on pip)
conda install -c conda-forge xtb-python rdkit

# Or use mamba (faster)
mamba install -c conda-forge xtb-python rdkit
```

## Basic Usage

```python
from adapters.xtb.adapter import xTBAdapter

# Create adapter
xtb = xTBAdapter()

# Run calculation
result = await xtb.execute("CCO")  # Ethanol

# Get results
if result.success:
    print(f"Energy: {result.data['energy_kcal_mol']:.2f} kcal/mol")
    print(f"HOMO-LUMO gap: {result.data['homo_lumo_gap_eV']:.2f} eV")
```

## Common Use Cases

### 1. Quick Property Calculation
```python
result = await xtb.execute("c1ccccc1")  # Benzene
# Returns: energy, HOMO/LUMO, dipole moment
```

### 2. Geometry Optimization
```python
result = await xtb.execute("CCO", optimize=True)
coords = result.data['optimized_coordinates']
```

### 3. Solvent Effects
```python
# Water solvation
result = await xtb.execute("CCO", solvent="water")
```

### 4. Charged Molecules
```python
# Acetate anion
result = await xtb.execute("CC(=O)[O-]", charge=-1)
```

### 5. Fast Screening (Use GFN1-xTB)
```python
xtb = xTBAdapter(config={"method": "GFN1-xTB"})
for smiles in molecule_list:
    result = await xtb.execute(smiles)
```

## Output Properties

| Property | Description | Units |
|----------|-------------|-------|
| `energy_kcal_mol` | Total energy | kcal/mol |
| `homo_energy_eV` | HOMO energy | eV |
| `lumo_energy_eV` | LUMO energy | eV |
| `homo_lumo_gap_eV` | HOMO-LUMO gap | eV |
| `dipole_moment_debye` | Dipole moment | Debye |
| `converged` | Optimization status | bool |

## GFN Method Selection

| Method | Speed | Accuracy | Best For |
|--------|-------|----------|----------|
| GFN2-xTB | Moderate | Highest | Property predictions |
| GFN1-xTB | Fast | Good | Large molecules |
| GFN0-xTB | Fastest | Lower | High-throughput screening |

## Configuration Options

```python
config = {
    "method": "GFN2-xTB",      # GFN method
    "optimize": True,           # Optimize geometry
    "charge": 0,                # Molecular charge
    "uhf": 0,                   # Unpaired electrons
    "solvent": "water",         # GBSA solvent
    "accuracy": 1.0             # Numerical accuracy
}

xtb = xTBAdapter(config=config)
```

## Performance Tips

1. **Use GFN1-xTB for large molecules** (>50 atoms)
2. **Disable optimization** for faster calculations (`optimize=False`)
3. **Caching is automatic** - identical calculations are instant
4. **Batch processing** - use GFN1-xTB or GFN0-xTB for speed

## Common Solvents

- `"water"` - Aqueous solution
- `"methanol"` - Methanol
- `"acetonitrile"` - Acetonitrile
- `"dmso"` - DMSO
- `"chloroform"` - Chloroform
- `None` - Gas phase (default)

## Error Handling

```python
result = await xtb.execute("invalid_smiles")
if not result.success:
    print(f"Error: {result.error}")
```

## Typical Runtime

| Molecule Size | GFN2-xTB | GFN1-xTB | GFN0-xTB |
|--------------|----------|----------|----------|
| 10-20 atoms | 2-3 sec | 1-2 sec | <1 sec |
| 20-50 atoms | 3-7 sec | 2-5 sec | 1-2 sec |
| 50-100 atoms | 7-15 sec | 5-10 sec | 2-5 sec |

## Example: Drug Molecule Screening

```python
# Screen drug candidates for electronic properties
drugs = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}

xtb = xTBAdapter(config={"method": "GFN1-xTB"})

for name, smiles in drugs.items():
    result = await xtb.execute(smiles, solvent="water")
    if result.success:
        gap = result.data['homo_lumo_gap_eV']
        print(f"{name}: HOMO-LUMO gap = {gap:.2f} eV")
```

## Troubleshooting

**xtb-python not found:**
```bash
conda install -c conda-forge xtb-python
```

**Slow calculations:**
- Use GFN1-xTB instead of GFN2-xTB
- Set `optimize=False` for single-point calculations

**Convergence issues:**
```python
xtb = xTBAdapter(config={"accuracy": 2.0})
```

## See Also

- [Full Documentation](README.md)
- [Example Usage](example_usage.py)
- [xTB Official Docs](https://xtb-docs.readthedocs.io/)
