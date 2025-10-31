# xTB Adapter for PharmForge

Fast semiempirical quantum mechanical calculations for molecular properties using the xTB (extended Tight-Binding) method.

## Overview

The xTB adapter provides rapid quantum chemistry calculations optimized for drug-like organic molecules. It uses the GFN (Geometry, Frequency, Noncovalent) family of semiempirical methods developed by the Grimme group.

**Key Features:**
- ⚡ Fast calculations (1-10 seconds per molecule)
- 🔬 Geometry optimization
- ⚛️ Electronic structure properties (HOMO/LUMO, orbital energies)
- 💧 Implicit solvent models (GBSA)
- 🎯 Accurate for drug-like molecules
- 📊 Molecular properties (dipole moment, polarizability)

## Installation

xTB requires installation via conda:

```bash
# Install xtb-python and dependencies
conda install -c conda-forge xtb-python rdkit

# Or using mamba (faster)
mamba install -c conda-forge xtb-python rdkit
```

**Note:** xtb-python is NOT available via pip. You must use conda/mamba.

## Usage

### Basic Usage

```python
from adapters.xtb.adapter import xTBAdapter

# Create adapter instance
xtb = xTBAdapter()

# Run calculation
result = await xtb.execute("CCO")  # Ethanol

# Access results
if result.success:
    data = result.data
    print(f"Energy: {data['energy_kcal_mol']:.2f} kcal/mol")
    print(f"HOMO-LUMO gap: {data['homo_lumo_gap_eV']:.2f} eV")
    print(f"Dipole moment: {data['dipole_moment_debye']:.2f} Debye")
```

### Advanced Configuration

```python
# Custom configuration
config = {
    "method": "GFN2-xTB",      # GFN0-xTB, GFN1-xTB, or GFN2-xTB
    "optimize": True,           # Optimize geometry
    "charge": 0,                # Molecular charge
    "uhf": 0,                   # Unpaired electrons
    "solvent": "water",         # GBSA solvent model
    "accuracy": 1.0             # Numerical accuracy
}

xtb = xTBAdapter(config=config)

# Run with custom parameters
result = await xtb.execute(
    "CCO",
    method="GFN1-xTB",
    optimize=True,
    solvent="water"
)
```

## GFN Methods

The adapter supports three GFN methods:

### GFN2-xTB (Default)
- **Most accurate** for drug-like molecules
- **Best for:** ADMET predictions, property calculations
- **Speed:** Moderate (2-10 seconds)
- **Use when:** Accuracy is important

### GFN1-xTB
- **Good balance** of speed and accuracy
- **Best for:** Large molecules (>100 atoms)
- **Speed:** Fast (1-5 seconds)
- **Use when:** Processing many molecules

### GFN0-xTB
- **Fastest** method
- **Best for:** High-throughput screening
- **Speed:** Very fast (<1 second)
- **Use when:** Speed is critical

## Output Properties

The adapter returns the following properties:

### Energy Properties
- `energy_hartree`: Total electronic energy (Hartree)
- `energy_kcal_mol`: Total electronic energy (kcal/mol)

### Electronic Structure
- `homo_energy_eV`: HOMO orbital energy (eV)
- `lumo_energy_eV`: LUMO orbital energy (eV)
- `homo_lumo_gap_eV`: HOMO-LUMO energy gap (eV)

### Molecular Properties
- `dipole_moment_debye`: Molecular dipole moment (Debye)

### Structural Information
- `n_atoms`: Number of atoms
- `optimized_coordinates`: 3D coordinates after optimization (if optimize=True)
- `converged`: Whether optimization converged

## Solvent Models

xTB supports implicit solvation via the GBSA (Generalized Born + Solvent Accessible Surface Area) model:

```python
# Common solvents
result = await xtb.execute("CCO", solvent="water")
result = await xtb.execute("CCO", solvent="methanol")
result = await xtb.execute("CCO", solvent="acetonitrile")
result = await xtb.execute("CCO", solvent="dmso")
result = await xtb.execute("CCO", solvent="chloroform")
```

Available solvents: water, methanol, acetonitrile, dmso, chloroform, and more.

## Use Cases

### 1. HOMO-LUMO Gap Analysis
```python
# Screen compounds for electronic properties
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]

for smiles in smiles_list:
    result = await xtb.execute(smiles)
    if result.success:
        gap = result.data['homo_lumo_gap_eV']
        print(f"{smiles}: HOMO-LUMO gap = {gap:.2f} eV")
```

### 2. Geometry Optimization
```python
# Optimize geometry and get final structure
result = await xtb.execute("CCO", optimize=True)

if result.success:
    coords = result.data['optimized_coordinates']
    energy = result.data['energy_kcal_mol']
    print(f"Optimized energy: {energy:.2f} kcal/mol")
```

### 3. Solvent Effects
```python
# Compare gas phase vs. aqueous phase
result_gas = await xtb.execute("CCO", solvent=None)
result_water = await xtb.execute("CCO", solvent="water")

if result_gas.success and result_water.success:
    delta_G_solv = (result_water.data['energy_kcal_mol'] -
                    result_gas.data['energy_kcal_mol'])
    print(f"Solvation free energy: {delta_G_solv:.2f} kcal/mol")
```

### 4. Charged Species
```python
# Calculate properties of ions
result = await xtb.execute("CC(=O)[O-]", charge=-1)  # Acetate anion
result = await xtb.execute("CC[NH3+]", charge=+1)    # Ethylammonium cation
```

## Performance

Typical calculation times on a modern CPU:

| Molecule Size | GFN0-xTB | GFN1-xTB | GFN2-xTB |
|--------------|----------|----------|----------|
| Small (10-20 atoms) | <1 sec | 1-2 sec | 2-3 sec |
| Medium (20-50 atoms) | 1-2 sec | 2-5 sec | 3-7 sec |
| Large (50-100 atoms) | 2-5 sec | 5-10 sec | 7-15 sec |

**Note:** Times include geometry optimization. Single-point calculations are ~5x faster.

## Caching

The adapter automatically caches results based on:
- SMILES (canonicalized)
- GFN method
- Optimization settings
- Charge and spin state
- Solvent model

```python
# First call: calculates from scratch
result1 = await xtb.execute("CCO")  # ~2 seconds

# Second call: retrieved from cache
result2 = await xtb.execute("CCO")  # <0.01 seconds
```

## Error Handling

```python
result = await xtb.execute("invalid_smiles")

if not result.success:
    print(f"Error: {result.error}")
else:
    # Process results
    data = result.data
```

Common errors:
- Invalid SMILES string
- xtb-python not installed
- Convergence failure (rare)

## Integration with PharmForge

The xTB adapter integrates seamlessly with PharmForge pipelines:

```python
from backend.core.pipeline import Pipeline

# Create pipeline with xTB calculations
pipeline = Pipeline()
pipeline.add_step("xtb", {
    "method": "GFN2-xTB",
    "optimize": True,
    "solvent": "water"
})

# Run pipeline
results = await pipeline.execute("CCO")
```

## Scientific Background

**Reference:**
Grimme, S., Bannwarth, C., & Shushkov, P. (2017).
*A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1–86).*
Journal of Chemical Theory and Computation, 13(5), 1989-2009.

**DOI:** 10.1021/acs.jctc.7b00118

## Limitations

- **Organic molecules only:** Best for C, H, N, O, S, P, halogens
- **Size limit:** Practical limit ~500 atoms
- **Not for:** Transition metals, highly exotic chemistry
- **Accuracy:** Good for trends, not benchmark-quality energies

## Troubleshooting

### xtb-python not installed
```bash
conda install -c conda-forge xtb-python
```

### Convergence failures
Try using a looser accuracy setting:
```python
xtb = xTBAdapter(config={"accuracy": 2.0})
```

### Memory issues with large molecules
Use GFN1-xTB or GFN0-xTB instead of GFN2-xTB:
```python
result = await xtb.execute("large_smiles", method="GFN1-xTB")
```

## See Also

- [xTB Documentation](https://xtb-docs.readthedocs.io/)
- [xtb-python GitHub](https://github.com/grimme-lab/xtb-python)
- [GFN Methods Paper](https://doi.org/10.1021/acs.jctc.7b00118)

## License

The xTB adapter follows PharmForge's MIT license. The underlying xtb program is licensed under LGPL-3.0.
