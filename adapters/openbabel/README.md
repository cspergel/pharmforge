# OpenBabel Adapter

The OpenBabel adapter provides molecular format conversion and ligand preparation capabilities for PharmForge.

## Installation

Install OpenBabel using one of the following methods:

```bash
# Using pip
pip install openbabel-wheel

# Using conda
conda install -c conda-forge openbabel
```

## Features

- **Format Conversion**: Convert between 20+ molecular formats (SMILES, InChI, MOL, MOL2, SDF, PDB, etc.)
- **3D Coordinate Generation**: Generate 3D structures from 1D/2D representations
- **Hydrogen Addition**: Add explicit hydrogens to molecular structures
- **Geometry Optimization**: Optimize molecular geometry using force fields (MMFF94, UFF, GAFF, Ghemical)
- **Property Calculation**: Compute basic molecular properties

## Usage

### Basic Format Conversion

```python
from adapters.openbabel import OpenBabelAdapter

adapter = OpenBabelAdapter()

# Convert SMILES to MOL2 format
result = await adapter.execute(
    "CCO",  # SMILES for ethanol
    output_format="mol2"
)

if result.success:
    mol2_structure = result.data["structure"]
    properties = result.data["properties"]
```

### 3D Structure Generation

```python
# Generate 3D coordinates with hydrogens
result = await adapter.execute(
    "CCO",
    output_format="pdb",
    gen_3d=True,
    add_hydrogens=True
)
```

### Geometry Optimization

```python
# Optimize geometry with MMFF94 force field
result = await adapter.execute(
    "CCO",
    output_format="mol2",
    gen_3d=True,
    add_hydrogens=True,
    optimize=True,
    force_field="mmff94"
)
```

### Convert from Non-SMILES Format

```python
# Input with explicit format specification
result = await adapter.execute(
    {
        "structure": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "format": "inchi"
    },
    output_format="sdf"
)
```

## Parameters

### Input Data

- **String**: SMILES or InChI string (auto-detected)
- **Dictionary**: `{"structure": "...", "format": "smi|inchi|mol|..."}`

### Keyword Arguments

- `output_format` (str): Target format (default: "mol2")
  - Supported: smi, smiles, inchi, inchikey, mol, mol2, sdf, pdb, xyz, cml, can
- `gen_3d` (bool): Generate 3D coordinates (default: True)
- `add_hydrogens` (bool): Add hydrogen atoms (default: True)
- `optimize` (bool): Optimize geometry (default: False)
- `force_field` (str): Force field for optimization (default: "mmff94")
  - Options: mmff94, uff, gaff, ghemical

## Output Format

```python
AdapterResult(
    success=True,
    data={
        "structure": "...",  # Converted structure string
        "format": "mol2",
        "properties": {
            "molecular_formula": "C2H6O",
            "molecular_weight": 46.069,
            "num_atoms": 9,
            "num_heavy_atoms": 3,
            "num_bonds": 8,
            "num_rotors": 1,
            "exact_mass": 46.0419,
            "has_3d": True,
            "energy": -12.345  # Only if optimized
        },
        "canonical_smiles": "CCO",  # Added automatically
        "inchi": "InChI=1S/C2H6O/...",  # Added automatically
        "inchikey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    },
    metadata={
        "source": "openbabel",
        "input_format": "smi",
        "output_format": "mol2",
        "operations": {
            "gen_3d": True,
            "add_hydrogens": True,
            "optimize": False
        }
    }
)
```

## Common Use Cases

### Ligand Preparation for Docking

```python
# Prepare ligand for molecular docking
result = await adapter.execute(
    "CCO",
    output_format="pdbqt",
    gen_3d=True,
    add_hydrogens=True,
    optimize=True,
    force_field="mmff94"
)
```

### Batch Format Conversion

```python
smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]

for smiles in smiles_list:
    result = await adapter.execute(
        smiles,
        output_format="sdf",
        gen_3d=True
    )
    if result.success:
        print(f"Converted: {smiles}")
```

### InChI to SMILES Conversion

```python
result = await adapter.execute(
    {
        "structure": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "format": "inchi"
    },
    output_format="can"
)

canonical_smiles = result.data["structure"]
```

## Error Handling

```python
result = await adapter.execute("invalid_smiles")

if not result.success:
    print(f"Error: {result.error}")
    # Error: Invalid molecular structure or format
```

## Integration with Other Adapters

OpenBabel adapter works well with other PharmForge adapters:

```python
# 1. Generate 3D structure with OpenBabel
ob_result = await openbabel_adapter.execute(
    "CCO",
    output_format="mol2",
    gen_3d=True,
    optimize=True
)

# 2. Use for molecular docking with Vina
vina_result = await vina_adapter.execute(
    ligand=ob_result.data["structure"],
    receptor="protein.pdb"
)

# 3. Calculate properties with RDKit
rdkit_result = await rdkit_adapter.execute(
    ob_result.data["canonical_smiles"]
)
```

## Supported Formats

### Input Formats
- SMILES (smi, smiles)
- InChI (inchi)
- InChIKey (inchikey)
- MDL MOL (mol)
- Mol2 (mol2)
- SDF (sdf)
- PDB (pdb)
- XYZ (xyz)
- CML (cml)
- Canonical SMILES (can)

### Output Formats
Same as input formats plus additional specialized formats depending on OpenBabel installation.

## Notes

- 3D coordinate generation and optimization can take longer for larger molecules
- MMFF94 is recommended for small molecules; UFF is more general but less accurate
- Energy values are in kcal/mol
- Geometry optimization uses up to 500 steps by default (configurable)

## References

- OpenBabel Documentation: http://openbabel.org/
- Force Fields: http://openbabel.org/wiki/ForceFields
