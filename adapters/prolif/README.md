# ProLIF Adapter for PharmForge

## Overview

The ProLIF (Protein-Ligand Interaction Fingerprints) adapter provides comprehensive protein-ligand interaction analysis for PharmForge. It identifies and quantifies various types of non-covalent interactions between proteins and ligands, making it ideal for analyzing docking poses or molecular dynamics trajectories.

## Features

- **Comprehensive Interaction Detection**: H-bonds, hydrophobic contacts, pi-stacking, salt bridges, halogen bonds, and more
- **Multi-Frame Support**: Analyze single structures or entire MD trajectories
- **Residue-Level Analysis**: Identify which specific residues interact with the ligand
- **Frequency Calculations**: Compute interaction occurrence frequencies across trajectory frames
- **Flexible Input**: Support for separate files, complex files, or MDAnalysis Universe objects

## Installation

```bash
# Required
pip install prolif
pip install MDAnalysis

# Optional (for enhanced features)
pip install pandas
```

## Supported Interaction Types

- **HBDonor/HBAcceptor**: Hydrogen bonds
- **Hydrophobic**: Hydrophobic contacts
- **PiStacking**: Pi-pi stacking interactions
- **PiCation/CationPi**: Pi-cation interactions
- **Anionic/Cationic**: Salt bridges
- **XBDonor/XBAcceptor**: Halogen bonds
- **MetalDonor/MetalAcceptor**: Metal coordination

## Input Formats

### 1. Separate Protein and Ligand Files

```python
input_data = {
    'protein_file': '/path/to/protein.pdb',
    'ligand_file': '/path/to/ligand.mol2'
}
```

### 2. Complex File with Ligand Selection

```python
input_data = {
    'complex_file': '/path/to/complex.pdb',
    'ligand_selection': 'resname LIG'  # MDAnalysis selection syntax
}
```

### 3. MDAnalysis Universe Object

```python
import MDAnalysis as mda

universe = mda.Universe('topology.pdb', 'trajectory.dcd')
input_data = {
    'universe': universe,
    'ligand_selection': 'resname LIG'
}
```

## Usage Examples

### Basic Analysis (Single Structure)

```python
from adapters.prolif.adapter import ProLIFAdapter

# Initialize adapter
adapter = ProLIFAdapter()

# Analyze a docking pose
input_data = {
    'complex_file': '/path/to/docked_complex.pdb',
    'ligand_selection': 'resname LIG'
}

result = await adapter.execute(input_data)

if result.success:
    print(f"Total interactions: {result.data['total_interactions']}")
    print(f"Interaction counts: {result.data['interaction_counts']}")
    print(f"Unique residue pairs: {result.data['n_unique_interactions']}")
else:
    print(f"Error: {result.error}")
```

### Custom Configuration

```python
adapter = ProLIFAdapter(config={
    'protein_selection': 'protein and not resname HOH',  # Exclude water
    'ligand_selection': 'resname LIG',
    'compute_frequency': True,
    'interactions': ['HBDonor', 'HBAcceptor', 'Hydrophobic', 'PiStacking']  # Specific interactions
})
```

### Trajectory Analysis (Multiple Frames)

```python
import MDAnalysis as mda

# Load trajectory
universe = mda.Universe('protein.pdb', 'trajectory.dcd')

input_data = {
    'universe': universe,
    'ligand_selection': 'resname LIG'
}

result = await adapter.execute(input_data)

if result.success:
    print(f"Analyzed {result.data['n_frames']} frames")
    print(f"Average interactions per frame: {result.data['avg_interactions_per_frame']:.2f}")
    print(f"Interaction frequencies: {result.data['interaction_frequencies']}")
```

### Runtime Parameter Overrides

```python
result = await adapter.execute(
    input_data,
    protein_selection='protein and resid 1-100',  # Override protein selection
    ligand_selection='resname LIG or resname MOL',  # Override ligand selection
    interactions=['HBDonor', 'HBAcceptor']  # Override interaction types
)
```

## Output Format

The adapter returns an `AdapterResult` with the following data structure:

```python
{
    "n_frames": 100,  # Number of frames analyzed
    "interaction_counts": {
        "HBDonor": 45,
        "HBAcceptor": 38,
        "Hydrophobic": 62,
        "PiStacking": 12
    },
    "interaction_frequencies": {  # Average per frame
        "HBDonor": 0.45,
        "HBAcceptor": 0.38,
        "Hydrophobic": 0.62,
        "PiStacking": 0.12
    },
    "residue_pairs": [
        {"residue": "ASP123", "interaction_type": "HBDonor"},
        {"residue": "PHE45", "interaction_type": "PiStacking"},
        ...
    ],
    "n_unique_interactions": 15,  # Unique residue-interaction pairs
    "total_interactions": 157,  # Total across all frames
    "avg_interactions_per_frame": 1.57,
    "most_common_interactions": [
        {"type": "Hydrophobic", "count": 62},
        {"type": "HBDonor", "count": 45},
        ...
    ],
    "input_type": "complex_file",
    "protein_selection": "protein",
    "ligand_selection": "resname LIG"
}
```

## Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protein_selection` | str | `"protein"` | MDAnalysis selection for protein |
| `ligand_selection` | str | `"resname LIG"` | MDAnalysis selection for ligand |
| `interactions` | list | `None` | List of interaction types (None = all) |
| `compute_frequency` | bool | `True` | Calculate interaction frequencies |
| `distance_cutoff` | float | `None` | Distance cutoff in Angstroms (None = ProLIF defaults) |

## Integration with PharmForge Workflow

### 1. Post-Docking Analysis

```python
# After Vina docking
vina_adapter = VinaAdapter(...)
docking_result = await vina_adapter.execute(smiles)

# Analyze interactions in best pose
prolif_adapter = ProLIFAdapter()
interaction_result = await prolif_adapter.execute({
    'complex_file': docking_result.data['output_file'],
    'ligand_selection': 'resname LIG'
})
```

### 2. MD Trajectory Analysis

```python
# After MD simulation with OpenMM
mdanalysis_adapter = MDAnalysisAdapter(...)
md_result = await mdanalysis_adapter.execute({
    'topology': 'protein.pdb',
    'trajectory': 'trajectory.dcd'
})

# Analyze interactions across trajectory
prolif_adapter = ProLIFAdapter()
interaction_result = await prolif_adapter.execute({
    'complex_file': 'protein.pdb',
    'ligand_selection': 'resname LIG'
})
```

## Error Handling

The adapter includes comprehensive error handling:

```python
result = await adapter.execute(input_data)

if not result.success:
    if "not installed" in result.error:
        print("Missing dependencies - install ProLIF and MDAnalysis")
    elif "does not exist" in result.error:
        print("File not found - check your file paths")
    else:
        print(f"Analysis failed: {result.error}")
```

## Performance Considerations

- **Single Frame**: Fast (seconds for typical protein-ligand complex)
- **Trajectory**: Linear scaling with number of frames
- **Large Proteins**: Consider using `protein_selection` to focus on binding site

### Optimization Tips

1. **Focus on Binding Site**: Use specific protein selection
   ```python
   config = {'protein_selection': 'protein and around 10 resname LIG'}
   ```

2. **Limit Interaction Types**: Only compute needed interactions
   ```python
   config = {'interactions': ['HBDonor', 'HBAcceptor', 'Hydrophobic']}
   ```

3. **Frame Subsampling**: For long trajectories, analyze every Nth frame
   ```python
   # In MDAnalysis
   universe = mda.Universe('top.pdb', 'traj.dcd')
   universe = universe.trajectory[::10]  # Every 10th frame
   ```

## Troubleshooting

### Common Issues

1. **ProLIF not installed**
   ```
   Error: ProLIF is not installed
   Solution: pip install prolif
   ```

2. **MDAnalysis not installed**
   ```
   Error: MDAnalysis is not installed
   Solution: pip install MDAnalysis
   ```

3. **No atoms selected**
   ```
   Error: Ligand selection returned no atoms
   Solution: Verify ligand_selection matches actual residue name
   ```

4. **Invalid file format**
   ```
   Error: Failed to load structures
   Solution: Ensure files are in supported format (PDB, MOL2, etc.)
   ```

## References

- ProLIF Documentation: https://prolif.readthedocs.io/
- MDAnalysis Documentation: https://docs.mdanalysis.org/
- ProLIF Paper: Bouysset & Fiorucci (2021) J. Cheminform. 13:72

## Adapter Protocol Compliance

This adapter follows the PharmForge AdapterProtocol:

- ✅ Inherits from `AdapterProtocol`
- ✅ Implements `validate_input()`
- ✅ Implements `execute()` (async)
- ✅ Returns `AdapterResult`
- ✅ Handles `ImportError` gracefully
- ✅ Implements `generate_cache_key()`
- ✅ Implements `get_metadata()`
- ✅ Provides comprehensive logging
- ✅ Thread-safe execution via asyncio

## License

This adapter is part of PharmForge and follows the same license terms.
