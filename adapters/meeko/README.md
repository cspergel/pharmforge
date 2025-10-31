# Meeko Ligand Preparation Adapter

A PharmForge adapter for preparing ligands for molecular docking using [Meeko](https://github.com/forlilab/Meeko).

## Overview

The Meeko adapter converts SMILES strings into PDBQT format files ready for docking with AutoDock Vina or Vina-GPU. It handles complex molecular features including macrocycles, flexible residues, and multiple conformers.

## Features

- **SMILES to PDBQT Conversion**: Automatic conversion from SMILES to docking-ready PDBQT format
- **Macrocycle Support**: Intelligent handling of macrocyclic compounds with ring flexibility
- **Multiple Conformers**: Generate and prepare multiple conformational states
- **Flexible Amides**: Optional amide bond flexibility for improved docking
- **Hydrated Docking**: Support for explicit hydration sites
- **Custom Rigidification**: Control over which bonds remain rigid using SMARTS patterns
- **Proper Atom Types**: Correct AutoDock atom type assignment and charge calculation

## Installation

```bash
pip install meeko rdkit
```

## Usage

### Basic Usage

```python
from adapters.meeko import MeekoAdapter

# Initialize adapter
meeko_adapter = MeekoAdapter()

# Prepare a ligand from SMILES
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
result = await meeko_adapter.execute(smiles)

if result.success:
    pdbqt_string = result.data["pdbqt_string"]
    num_torsions = result.data["conformers"][0]["metadata"]["num_torsions"]
    print(f"Prepared ligand with {num_torsions} rotatable bonds")
    print(pdbqt_string)
```

### Advanced Configuration

```python
# Configure adapter with custom settings
config = {
    "rigid_macrocycles": False,        # Allow macrocycle flexibility
    "keep_nonpolar_hydrogens": False,  # Remove nonpolar hydrogens
    "hydrate": False,                  # Don't add hydration sites
    "flexible_amides": False,          # Keep amides rigid
    "min_ring_size": 7,                # Minimum ring size for flexibility
    "num_conformers": 1,               # Number of conformers to generate
    "double_bond_penalty": 50.0        # Penalty for double bond torsions
}

meeko_adapter = MeekoAdapter(config=config)
```

### Multiple Conformers

```python
# Generate multiple conformers for ensemble docking
result = await meeko_adapter.execute(
    smiles,
    num_conformers=10  # Generate 10 different conformations
)

if result.success:
    for conformer in result.data["conformers"]:
        conformer_id = conformer["conformer_id"]
        pdbqt = conformer["pdbqt_string"]
        torsions = conformer["metadata"]["num_torsions"]
        print(f"Conformer {conformer_id}: {torsions} torsions")
```

### Macrocycle Handling

```python
# Prepare a macrocyclic ligand with flexible rings
macrocycle_smiles = "C1CCCCCCCCCCC1"  # Simple macrocycle

result = await meeko_adapter.execute(
    macrocycle_smiles,
    rigid_macrocycles=False  # Allow ring flexibility
)

if result.success:
    metadata = result.data["conformers"][0]["metadata"]
    print(f"Macrocycles detected: {metadata['num_macrocycles']}")
    print(f"Macrocycle info: {metadata['macrocycles_detected']}")
```

### Integration with Vina Adapter

```python
from adapters.meeko import MeekoAdapter
from adapters.vina import VinaAdapter

# Step 1: Prepare ligand with Meeko
meeko = MeekoAdapter()
ligand_result = await meeko.execute(smiles)

if ligand_result.success:
    pdbqt_string = ligand_result.data["pdbqt_string"]

    # Save to temporary file for Vina
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as f:
        f.write(pdbqt_string)
        ligand_path = f.name

    # Step 2: Dock with Vina
    vina = VinaAdapter(config={
        "receptor_path": "receptor.pdbqt",
        "center_x": 10.0,
        "center_y": 15.0,
        "center_z": 20.0
    })

    docking_result = await vina.execute(smiles, ligand_pdbqt=ligand_path)
```

### File Output

```python
# Save to file instead of returning string
result = await meeko_adapter.execute(
    smiles,
    output_format='file'  # Returns file path instead of string
)

if result.success:
    pdbqt_file = result.data["pdbqt_file"]
    print(f"PDBQT saved to: {pdbqt_file}")
```

## Configuration Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `rigid_macrocycles` | bool | `False` | Keep macrocycles rigid (no ring flexibility) |
| `keep_nonpolar_hydrogens` | bool | `False` | Retain nonpolar hydrogen atoms |
| `merge_these_atom_types` | list | `None` | List of atom types to merge |
| `hydrate` | bool | `False` | Add explicit hydration sites |
| `flexible_amides` | bool | `False` | Allow amide bond rotation |
| `rigidify_bonds_smarts` | list | `[]` | SMARTS patterns for bonds to rigidify |
| `rigidify_bonds_indices` | list | `[]` | Specific bond indices to rigidify |
| `double_bond_penalty` | float | `50.0` | Energy penalty for double bond torsions |
| `min_ring_size` | int | `7` | Minimum ring size for flexibility |
| `num_conformers` | int | `1` | Number of conformers to generate |

## Output Format

### Success Response

```python
{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "pdbqt_string": "REMARK ... (full PDBQT content)",
    "num_conformers": 1,
    "conformers": [
        {
            "conformer_id": 0,
            "pdbqt_string": "REMARK ...",
            "metadata": {
                "num_torsions": 4,
                "macrocycles_detected": [],
                "num_macrocycles": 0,
                "num_atoms": 21,
                "num_heavy_atoms": 13,
                "rigid_macrocycles": False,
                "flexible_amides": False,
                "hydrated": False
            }
        }
    ]
}
```

### Metadata

```python
result.metadata = {
    "adapter_name": "meeko",
    "adapter_version": "1.0.0",
    "num_conformers": 1,
    "primary_conformer": {
        "num_torsions": 4,
        "macrocycles_detected": [],
        "num_macrocycles": 0,
        "num_atoms": 21,
        "num_heavy_atoms": 13
    },
    "meeko_config": {
        "rigid_macrocycles": False,
        "keep_nonpolar_hydrogens": False,
        "hydrate": False,
        "flexible_amides": False,
        "min_ring_size": 7
    }
}
```

## Error Handling

The adapter gracefully handles errors and missing dependencies:

```python
# Missing dependencies
if not result.success:
    if "RDKit is not installed" in result.error:
        print("Install RDKit: pip install rdkit")
    elif "Meeko is not installed" in result.error:
        print("Install Meeko: pip install meeko")
    else:
        print(f"Error: {result.error}")
```

## Integration Examples

### Complete Docking Workflow

```python
from adapters.meeko import MeekoAdapter
from adapters.vina import VinaAdapter
import tempfile

async def dock_molecule(smiles, receptor_pdbqt, center, box_size):
    """Complete docking workflow with Meeko preparation"""

    # Step 1: Prepare ligand
    meeko = MeekoAdapter(config={"num_conformers": 10})
    prep_result = await meeko.execute(smiles, num_conformers=10)

    if not prep_result.success:
        return {"error": prep_result.error}

    # Step 2: Dock each conformer
    vina = VinaAdapter(config={
        "receptor_path": receptor_pdbqt,
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": box_size[0],
        "size_y": box_size[1],
        "size_z": box_size[2]
    })

    best_score = None
    best_pose = None

    for conformer in prep_result.data["conformers"]:
        # Save conformer to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as f:
            f.write(conformer["pdbqt_string"])
            ligand_path = f.name

        # Dock conformer
        dock_result = await vina.execute(smiles, ligand_pdbqt=ligand_path)

        if dock_result.success:
            affinity = dock_result.data["binding_affinity"]
            if best_score is None or affinity < best_score:
                best_score = affinity
                best_pose = dock_result.data

    return {
        "smiles": smiles,
        "best_affinity": best_score,
        "best_pose": best_pose,
        "num_conformers_tested": len(prep_result.data["conformers"])
    }

# Usage
result = await dock_molecule(
    smiles="CC(=O)Oc1ccccc1C(=O)O",
    receptor_pdbqt="protein.pdbqt",
    center=[10.0, 15.0, 20.0],
    box_size=[25, 25, 25]
)
```

### Batch Ligand Preparation

```python
async def prepare_ligand_library(smiles_list):
    """Prepare multiple ligands for docking"""

    meeko = MeekoAdapter()
    prepared_ligands = []

    for smiles in smiles_list:
        result = await meeko.execute(smiles)

        if result.success:
            prepared_ligands.append({
                "smiles": smiles,
                "pdbqt": result.data["pdbqt_string"],
                "torsions": result.data["conformers"][0]["metadata"]["num_torsions"]
            })
        else:
            print(f"Failed to prepare {smiles}: {result.error}")

    return prepared_ligands

# Usage
library = await prepare_ligand_library([
    "CC(=O)Oc1ccccc1C(=O)O",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
])
```

## Troubleshooting

### Common Issues

1. **ImportError: No module named 'meeko'**
   ```bash
   pip install meeko
   ```

2. **ImportError: No module named 'rdkit'**
   ```bash
   pip install rdkit
   ```

3. **Invalid SMILES string**
   - Ensure SMILES string is valid
   - Check for special characters or syntax errors
   - Use RDKit to validate: `Chem.MolFromSmiles(smiles)`

4. **Conformer generation fails**
   - Some molecules may be difficult to embed in 3D
   - Try reducing `num_conformers`
   - Check for unusual molecular structures

5. **Macrocycle flexibility issues**
   - Set `rigid_macrocycles=False` for flexible macrocycles
   - Adjust `min_ring_size` parameter
   - Check Meeko documentation for macrocycle handling

## References

- [Meeko GitHub Repository](https://github.com/forlilab/Meeko)
- [Meeko Documentation](https://meeko.readthedocs.io/)
- [AutoDock Vina](http://vina.scripps.edu/)
- [RDKit](https://www.rdkit.org/)

## License

This adapter follows the PharmForge licensing terms. Meeko is licensed under the Apache License 2.0.

## Contributing

Contributions are welcome! Please ensure:
- Code follows PharmForge adapter protocol
- All dependencies are properly handled
- Error cases are gracefully managed
- Documentation is updated
