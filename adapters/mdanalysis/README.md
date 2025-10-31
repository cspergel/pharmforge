# MDAnalysis Adapter for PharmForge

This adapter provides molecular dynamics trajectory analysis using the [MDAnalysis](https://www.mdanalysis.org/) toolkit.

## Installation

```bash
pip install MDAnalysis
```

## Features

- **RMSD Calculation**: Root Mean Square Deviation over trajectory
- **RMSF Calculation**: Root Mean Square Fluctuation (per-atom or per-residue)
- **Radius of Gyration**: Track molecular compactness over time
- **Hydrogen Bond Analysis**: Identify and analyze H-bonds (optional)
- **Multiple Format Support**: DCD, XTC, TRR, NetCDF, H5MD, and more
- **Trajectory Alignment**: Automatic alignment before analysis

## Usage

### Basic Usage

```python
from adapters.mdanalysis import MDAnalysisAdapter

# Initialize adapter
adapter = MDAnalysisAdapter()

# Prepare input data
input_data = {
    "topology": "/path/to/topology.pdb",
    "trajectory": "/path/to/trajectory.dcd"
}

# Run analysis
result = await adapter.execute(input_data)

if result.success:
    print("RMSD:", result.data["rmsd"]["rmsd_mean"])
    print("RMSF:", result.data["rmsf"]["rmsf_mean"])
    print("Radius of gyration:", result.data["radius_of_gyration"]["rg_mean"])
```

### Custom Configuration

```python
# Configure which analyses to run
config = {
    "selection": "protein and name CA",  # C-alpha atoms only
    "compute_rmsd": True,
    "compute_rmsf": True,
    "compute_rg": True,
    "compute_hbonds": False,  # Hydrogen bonds (more expensive)
    "align_trajectory": True,
    "reference_frame": 0,  # First frame as reference
    "step": 10  # Analyze every 10th frame for speed
}

adapter = MDAnalysisAdapter(config=config)
```

### Single Structure Analysis

For analyzing a single structure without a trajectory:

```python
input_data = {
    "topology": "/path/to/structure.pdb"
    # No trajectory needed
}

result = await adapter.execute(input_data)
# Only radius of gyration will be computed (single-frame metric)
```

### Runtime Configuration Override

```python
# Override configuration at runtime
result = await adapter.execute(
    input_data,
    selection="backbone",
    compute_rmsd=True,
    compute_hbonds=True
)
```

## Input Format

The adapter expects a dictionary with the following keys:

- **topology** (required): Path to topology file
  - Supported formats: PDB, PSF, GRO, TPR, PRMTOP, TOP

- **trajectory** (optional): Path to trajectory file
  - Supported formats: DCD, XTC, TRR, NetCDF, H5MD, LAMMPS, CRD
  - If omitted, only single-frame analysis is performed

## Output Format

```python
{
    "system": {
        "n_atoms": 1000,
        "n_residues": 125,
        "n_frames": 500,
        "total_time": 50000.0,  # ps
        "dt": 100.0  # ps
    },
    "rmsd": {
        "n_frames": 500,
        "rmsd_mean": 2.5,  # Angstroms
        "rmsd_std": 0.8,
        "rmsd_min": 1.2,
        "rmsd_max": 4.3,
        "rmsd_values": [1.2, 1.5, ...],
        "time_values": [0.0, 100.0, ...],
        "selection": "protein and name CA",
        "n_atoms": 125
    },
    "rmsf": {
        "n_atoms": 125,
        "rmsf_mean": 1.5,  # Angstroms
        "rmsf_std": 0.5,
        "rmsf_min": 0.5,
        "rmsf_max": 3.2,
        "rmsf_values": [0.5, 0.8, ...],
        "selection": "protein and name CA"
    },
    "radius_of_gyration": {
        "n_frames": 500,
        "rg_mean": 15.2,  # Angstroms
        "rg_std": 0.3,
        "rg_min": 14.8,
        "rg_max": 15.8,
        "rg_values": [15.0, 15.1, ...],
        "time_values": [0.0, 100.0, ...],
        "selection": "protein and name CA",
        "n_atoms": 125
    },
    "hydrogen_bonds": {
        "n_hbonds_total": 1500,
        "n_unique_hbonds": 25,
        "avg_hbonds_per_frame": 3.0,
        "distance_cutoff": 3.0,  # Angstroms
        "angle_cutoff": 150.0  # Degrees
    }
}
```

## Atom Selection Syntax

MDAnalysis uses a powerful selection language. Common examples:

- `"protein"` - All protein atoms
- `"protein and name CA"` - C-alpha atoms
- `"backbone"` - Backbone atoms (N, CA, C, O)
- `"resid 1:50"` - Residues 1 through 50
- `"not protein"` - Everything except protein
- `"protein or nucleic"` - Protein or nucleic acids
- `"around 5.0 protein"` - Atoms within 5 Angstroms of protein

See [MDAnalysis Selection Documentation](https://docs.mdanalysis.org/stable/documentation_pages/selections.html) for more.

## Performance Tips

1. **Use step parameter**: Analyze every Nth frame to reduce computation time
   ```python
   config = {"step": 10}  # Analyze every 10th frame
   ```

2. **Select fewer atoms**: Use specific selections to reduce calculation time
   ```python
   config = {"selection": "protein and name CA"}  # C-alpha only
   ```

3. **Disable expensive analyses**: Turn off hydrogen bond analysis for faster results
   ```python
   config = {"compute_hbonds": False}
   ```

4. **Use aligned trajectories**: Pre-align trajectories externally if running multiple analyses
   ```python
   config = {"align_trajectory": False}  # Skip alignment
   ```

## Error Handling

The adapter gracefully handles errors:

```python
result = await adapter.execute(input_data)

if not result.success:
    print(f"Error: {result.error}")
    print(f"Help: {result.metadata.get('installation_help')}")
```

Common error scenarios:
- MDAnalysis not installed
- Invalid file paths
- Unsupported file formats
- No atoms matching selection
- Corrupted trajectory files

## References

- **MDAnalysis Paper 1**: Michaud-Agrawal et al., J. Comput. Chem. 2011, 32, 2319-2327
- **MDAnalysis Paper 2**: Gowers et al., Proc. of the 15th Python in Science Conf. 2016, 98-105
- **Website**: https://www.mdanalysis.org/
- **Documentation**: https://docs.mdanalysis.org/

## Computational Cost

| Analysis Type | Complexity | Typical Time |
|--------------|------------|--------------|
| RMSD | O(n_frames × n_atoms) | Fast |
| RMSF | O(n_frames × n_atoms) | Fast |
| Radius of Gyration | O(n_frames × n_atoms) | Fast |
| Hydrogen Bonds | O(n_frames × n_atoms²) | Moderate |

For a typical protein (1000 atoms, 1000 frames):
- RMSD/RMSF/Rg: 1-5 seconds
- Hydrogen bonds: 10-60 seconds

## Examples

### Example 1: Protein Stability Analysis

```python
# Analyze protein backbone stability
config = {
    "selection": "backbone",
    "compute_rmsd": True,
    "compute_rmsf": True,
    "compute_rg": True,
    "align_trajectory": True
}

adapter = MDAnalysisAdapter(config=config)
result = await adapter.execute({
    "topology": "protein.pdb",
    "trajectory": "md_simulation.dcd"
})

# Check if structure is stable (low RMSD fluctuation)
if result.success:
    rmsd_std = result.data["rmsd"]["rmsd_std"]
    if rmsd_std < 1.0:
        print("Structure is stable")
    else:
        print("Structure shows significant fluctuation")
```

### Example 2: Ligand Binding Site Analysis

```python
# Analyze binding site residues
config = {
    "selection": "resid 10:30 and name CA",  # Binding site residues
    "compute_rmsf": True  # Check flexibility
}

result = await adapter.execute(input_data, **config)

# Identify flexible residues
if result.success:
    rmsf_values = result.data["rmsf"]["rmsf_values"]
    flexible_residues = [i for i, v in enumerate(rmsf_values) if v > 2.0]
    print(f"Flexible residues: {flexible_residues}")
```

### Example 3: Quick Analysis (Every 10th Frame)

```python
# Fast analysis for large trajectories
adapter = MDAnalysisAdapter(config={"step": 10})

result = await adapter.execute({
    "topology": "large_system.gro",
    "trajectory": "long_trajectory.xtc"
})
```
