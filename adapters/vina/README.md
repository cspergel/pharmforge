# Vina Adapter - Molecular Docking

AutoDock Vina adapter for PharmForge that performs protein-ligand docking to predict binding affinity and poses.

## Overview

The Vina adapter integrates AutoDock Vina molecular docking into PharmForge pipelines. It takes SMILES strings as input, generates 3D conformers, prepares them in PDBQT format, and performs docking against a specified protein receptor.

## Features

- **SMILES to 3D conversion** using RDKit
- **Automated ligand preparation** with Gasteiger charge assignment
- **Flexible docking box configuration** (center and size)
- **Score normalization** (0-1 scale, higher = better binding)
- **Multiple pose generation** with RMSD clustering
- **Caching support** for reproducibility
- **Async execution** for pipeline integration

## Requirements

### Software

1. **AutoDock Vina** (v1.2.0+)
   - Download: https://github.com/ccsb-scripps/AutoDock-Vina/releases
   - Install and add to PATH

2. **Python packages:**
   ```bash
   pip install rdkit
   ```

### Data

- **Receptor file** in PDBQT format
- **Docking box coordinates** (center x,y,z and size x,y,z)

## Installation

### 1. Install AutoDock Vina

**Linux:**
```bash
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64
sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina
```

**macOS:**
```bash
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_macos_x86_64
chmod +x vina_1.2.5_macos_x86_64
sudo mv vina_1.2.5_macos_x86_64 /usr/local/bin/vina
```

**Windows:**
- Download `vina_1.2.5_windows_x86_64.exe`
- Rename to `vina.exe`
- Add directory to PATH

### 2. Verify Installation

```bash
vina --version
```

## Configuration

The adapter requires configuration for the receptor and docking box:

```python
config = {
    # Required
    'receptor_path': '/path/to/receptor.pdbqt',
    'center_x': 10.0,  # Docking box center X coordinate
    'center_y': 20.0,  # Docking box center Y coordinate
    'center_z': 15.0,  # Docking box center Z coordinate

    # Optional (with defaults)
    'size_x': 25,              # Box size X (Angstroms)
    'size_y': 25,              # Box size Y (Angstroms)
    'size_z': 25,              # Box size Z (Angstroms)
    'exhaustiveness': 8,       # Search exhaustiveness (1-32)
    'num_modes': 9,            # Number of binding modes
    'energy_range': 3,         # Energy range for output (kcal/mol)
    'vina_binary': 'vina'      # Path to vina executable
}
```

### Finding Docking Box Coordinates

Use PyMOL, Chimera, or VMD to visualize your receptor and identify the binding site:

1. **Load receptor** in molecular viewer
2. **Identify binding site** (e.g., known ligand position, catalytic residues)
3. **Measure coordinates** of binding site center
4. **Set box size** to cover binding site (typically 20-30 Å)

Example using PyMOL:
```python
# In PyMOL console
get_extent binding_site
# Returns: [[x_min, y_min, z_min], [x_max, y_max, z_max]]
# Calculate center: (min + max) / 2
```

## Usage

### Basic Usage

```python
import asyncio
from adapters.vina.adapter import VinaAdapter

async def dock_molecule():
    # Configure adapter
    adapter = VinaAdapter(config={
        'receptor_path': '/path/to/receptor.pdbqt',
        'center_x': 10.0,
        'center_y': 20.0,
        'center_z': 15.0
    })

    # Dock a molecule
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    result = await adapter.execute(smiles)

    if result.success:
        print(f"Binding affinity: {result.data['binding_affinity']:.2f} kcal/mol")
        print(f"Normalized score: {result.data['binding_score']:.3f}")
        print(f"Number of poses: {result.data['num_poses']}")
    else:
        print(f"Docking failed: {result.error}")

asyncio.run(dock_molecule())
```

### Batch Screening

```python
async def screen_library():
    adapter = VinaAdapter(config={...})

    compounds = [
        ("CHEMBL1", "CCO"),
        ("CHEMBL2", "CC(=O)Oc1ccccc1C(=O)O"),
        # ... more compounds
    ]

    results = []
    for compound_id, smiles in compounds:
        result = await adapter.execute(smiles)
        if result.success:
            results.append({
                'id': compound_id,
                'smiles': smiles,
                'affinity': result.data['binding_affinity'],
                'score': result.data['binding_score']
            })

    # Rank by binding affinity
    ranked = sorted(results, key=lambda x: x['affinity'])

    for i, compound in enumerate(ranked[:10], 1):
        print(f"{i}. {compound['id']}: {compound['affinity']:.2f} kcal/mol")
```

### Custom Docking Box per Molecule

```python
async def dock_with_custom_box():
    adapter = VinaAdapter(config={'receptor_path': '...'})

    # Override docking box for specific molecule
    result = await adapter.execute(
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        center_x=15.0,  # Custom center
        center_y=22.0,
        center_z=18.0,
        size_x=20,      # Custom size
        size_y=20,
        size_z=20
    )
```

## Output Format

The adapter returns an `AdapterResult` with the following data structure:

```python
{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "binding_affinity": -7.5,        # Raw score (kcal/mol)
    "binding_score": 0.85,           # Normalized 0-1 (higher=better)
    "best_pose": {
        "mode": 1,
        "affinity": -7.5,
        "rmsd_lb": 0.0,
        "rmsd_ub": 0.0
    },
    "all_poses": [
        {"mode": 1, "affinity": -7.5, "rmsd_lb": 0.0, "rmsd_ub": 0.0},
        {"mode": 2, "affinity": -7.2, "rmsd_lb": 1.2, "rmsd_ub": 2.5},
        # ... more poses
    ],
    "num_poses": 9,
    "receptor": "receptor.pdbqt",
    "docking_box": {
        "center": [10.0, 20.0, 15.0],
        "size": [25, 25, 25]
    }
}
```

## Score Normalization

The adapter uses `vina_affinity_to01()` from `backend.core.scoring_utils` to normalize binding affinities:

- **Input range:** -12 to -4 kcal/mol (typical Vina range)
- **Output range:** 0 to 1 (higher = better binding)
- **Mapping:**
  - -12 kcal/mol → 1.0 (excellent binding)
  - -8 kcal/mol → 0.5 (moderate binding)
  - -4 kcal/mol → 0.0 (weak binding)

This normalization allows Vina scores to be compared with other objectives in multi-objective optimization.

## Preparing Receptor Files

### Option 1: Using MGLTools (Recommended)

```bash
# Download MGLTools
wget https://ccsb.scripps.edu/mgltools/downloads/

# Install and prepare receptor
pythonsh prepare_receptor4.py -r protein.pdb -o receptor.pdbqt
```

### Option 2: Using Open Babel

```bash
# Convert PDB to PDBQT
obabel protein.pdb -O receptor.pdbqt -xr
```

### Option 3: Manual Preparation

1. Remove water molecules and heteroatoms
2. Add hydrogens (pH 7.4)
3. Assign Gasteiger charges
4. Convert to PDBQT format

## Performance Considerations

### Exhaustiveness

- **Low (1-4):** Fast but less accurate (~2-5 min/ligand)
- **Medium (8):** Balanced (default, ~5-10 min/ligand)
- **High (16-32):** Slow but thorough (~15-30 min/ligand)

### Box Size

- **Small (15-20 Å):** Fast, for known binding sites
- **Medium (25-30 Å):** Balanced (default)
- **Large (40+ Å):** Slow, for blind docking

### Parallelization

For large virtual screens, run multiple adapters in parallel:

```python
import asyncio

async def screen_parallel(compounds, adapter, max_concurrent=4):
    semaphore = asyncio.Semaphore(max_concurrent)

    async def dock_with_semaphore(smiles):
        async with semaphore:
            return await adapter.execute(smiles)

    tasks = [dock_with_semaphore(smi) for smi in compounds]
    results = await asyncio.gather(*tasks)
    return results
```

## Testing

Run the test suite:

```bash
cd backend
pytest tests/test_vina_adapter.py -v
```

Run example usage:

```bash
python adapters/vina/example_usage.py
```

## Troubleshooting

### "Vina binary not found"

- Ensure Vina is installed and in PATH
- Or specify full path: `'vina_binary': '/usr/local/bin/vina'`

### "Could not parse SMILES"

- Check SMILES syntax is valid
- Use canonical SMILES from RDKit or PubChem

### "Receptor file not found"

- Verify receptor_path is correct
- Use absolute paths
- Ensure file is in PDBQT format

### "Docking failed" errors

- Check receptor and ligand files are valid PDBQT
- Verify docking box encompasses binding site
- Increase box size if ligand is large
- Check Vina version compatibility

### Slow performance

- Reduce exhaustiveness (e.g., 4 instead of 8)
- Reduce num_modes (e.g., 3 instead of 9)
- Decrease box size if possible
- Use parallel execution for batch screening

## References

- **AutoDock Vina:** https://vina.scripps.edu/
- **Vina GitHub:** https://github.com/ccsb-scripps/AutoDock-Vina
- **Original Paper:** Trott & Olson (2010) J. Comput. Chem. 31:455-461
- **Documentation:** https://autodock-vina.readthedocs.io/

## License

This adapter follows PharmForge's MIT license. AutoDock Vina is licensed under Apache 2.0.

## Support

For issues with:
- **Adapter code:** Open issue on PharmForge GitHub
- **Vina software:** See AutoDock Vina documentation
- **Receptor preparation:** Consult MGLTools or Open Babel docs
