# gmx_MMPBSA Adapter for PharmForge

Calculate binding free energies from molecular dynamics trajectories using Molecular Mechanics Poisson-Boltzmann Surface Area (MM/PBSA) and Molecular Mechanics Generalized Born Surface Area (MM/GBSA) methods.

## Overview

The gmx_MMPBSA adapter integrates the [gmx_MMPBSA](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA) tool into PharmForge, enabling accurate binding free energy calculations from GROMACS molecular dynamics simulations. This adapter is essential for validating docking predictions and ranking protein-ligand complexes based on their thermodynamic stability.

## Features

- **MM/PBSA Calculations**: Poisson-Boltzmann solvation model for accurate polar solvation
- **MM/GBSA Calculations**: Generalized Born model for faster approximate calculations
- **Energy Decomposition**: Detailed breakdown of binding energy components:
  - Van der Waals energy
  - Electrostatic energy
  - Polar solvation energy
  - Non-polar solvation (SASA)
- **Per-Residue Decomposition**: Optional analysis of per-residue contributions
- **Multiple Trajectory Formats**: Supports GROMACS (.xtc, .trr) trajectories
- **Flexible Topology Support**: Compatible with AMBER, GROMACS, CHARMM, NAMD topologies

## Installation

### Step 1: Install AmberTools

gmx_MMPBSA requires AmberTools for the underlying MMPBSA.py calculations:

```bash
# Using conda (recommended)
conda install -c conda-forge ambertools

# Or using mamba (faster)
mamba install -c conda-forge ambertools
```

### Step 2: Install gmx_MMPBSA

```bash
pip install gmx_MMPBSA
```

### Step 3: Install GROMACS (Optional but Recommended)

For trajectory preprocessing and index file generation:

```bash
# Using conda
conda install -c conda-forge gromacs

# Or follow GROMACS installation instructions:
# https://manual.gromacs.org/documentation/current/install-guide/index.html
```

### Verification

Verify the installation:

```bash
gmx_MMPBSA --version
```

Expected output:
```
gmx_MMPBSA v1.5.x
```

## Usage

### Basic Usage

```python
from adapters.gmx_mmpbsa import GmxMMPBSAAdapter

# Initialize adapter
adapter = GmxMMPBSAAdapter(
    config={
        "method": "pb",  # Use MM/PBSA
        "saltcon": 0.15,  # 150 mM salt concentration
        "temperature": 298.15  # 298.15 K
    }
)

# Prepare input data
input_data = {
    "topology_file": "path/to/complex.prmtop",  # Topology file
    "trajectory_file": "path/to/trajectory.xtc",  # MD trajectory
    "index_file": "path/to/index.ndx",  # GROMACS index file
    "receptor_group": "Protein",  # Receptor group name
    "ligand_group": "Ligand"  # Ligand group name
}

# Run calculation
result = await adapter.execute(input_data)

if result.success:
    data = result.data
    print(f"Binding Free Energy: {data['binding_free_energy']:.2f} kcal/mol")
    print(f"ΔH: {data['delta_h']:.2f} kcal/mol")
    print(f"Energy Components:")
    print(f"  Van der Waals: {data['components']['van_der_waals']:.2f}")
    print(f"  Electrostatic: {data['components']['electrostatic']:.2f}")
    print(f"  Polar Solvation: {data['components']['polar_solvation']:.2f}")
    print(f"  SASA: {data['components']['sasa']:.2f}")
else:
    print(f"Error: {result.error}")
```

### Advanced: Per-Residue Decomposition

```python
# Enable per-residue decomposition
adapter = GmxMMPBSAAdapter(
    config={
        "method": "pb",
        "decomp": True,  # Enable decomposition
        "startframe": 10,  # Skip equilibration frames
        "endframe": -1,  # Use all remaining frames
        "interval": 5  # Analyze every 5th frame
    }
)

result = await adapter.execute(input_data)

if result.success:
    print(f"Binding Energy: {result.data['binding_free_energy']:.2f} kcal/mol")

    # Per-residue contributions
    if result.data['per_residue_decomposition']:
        print("\nKey residues contributing to binding:")
        for residue in result.data['per_residue_decomposition'][:10]:
            print(f"  {residue['name']}: {residue['energy']:.2f} kcal/mol")
```

### Using MM/GBSA (Faster)

```python
# Use Generalized Born solvation model
adapter = GmxMMPBSAAdapter(
    config={
        "method": "gb",  # Use MM/GBSA
        "igb": 5,  # GB model (1-8, 5 is recommended)
        "saltcon": 0.15
    }
)

result = await adapter.execute(input_data)
```

## Input Files

### Required Files

1. **Topology File** (`topology_file`):
   - Formats: `.prmtop` (AMBER), `.gro` (GROMACS), `.psf` (CHARMM)
   - Contains atomic coordinates and force field parameters
   - Must include the complete protein-ligand complex

2. **Trajectory File** (`trajectory_file`):
   - Formats: `.xtc`, `.trr` (GROMACS), `.nc` (NetCDF)
   - Contains molecular dynamics trajectory frames
   - Should be equilibrated and production MD frames

3. **Index File** (`index_file`):
   - Format: `.ndx` (GROMACS index file)
   - Defines atom groups for receptor and ligand
   - Generated using `gmx make_ndx` or during simulation setup

### Creating an Index File

If you don't have an index file, create one using GROMACS:

```bash
# Create index file with protein and ligand groups
gmx make_ndx -f complex.gro -o index.ndx

# In the interactive prompt:
# > 1 | 13  (combine Protein and Ligand into Complex)
# > name 20 Complex
# > q
```

Example index file structure:
```
[ Protein ]
   1    2    3    4  ...

[ Ligand ]
 1234 1235 1236 1237 ...

[ Complex ]
   1    2    3    4  ... 1234 1235 1236 1237 ...
```

## Configuration Options

### Method Selection

- **`method`**: `"pb"` (MM/PBSA) or `"gb"` (MM/GBSA)
  - MM/PBSA: More accurate but slower
  - MM/GBSA: Faster but less accurate
  - Default: `"pb"`

### Frame Selection

- **`startframe`**: First frame to analyze (default: `1`)
- **`endframe`**: Last frame to analyze (default: `-1`, all frames)
- **`interval`**: Analyze every Nth frame (default: `1`)

Example: Analyze frames 100-200, every 5th frame:
```python
config = {
    "startframe": 100,
    "endframe": 200,
    "interval": 5
}
```

### Solvation Parameters

- **`saltcon`**: Salt concentration in Molarity (default: `0.15`)
- **`temperature`**: Temperature in Kelvin (default: `298.15`)

### GB-Specific Options

- **`igb`**: GB model selection (1-8, default: `5`)
  - `1`: Hawkins, Cramer, Truhlar model
  - `2`: Modified GB model (mbondi radii)
  - `5`: Modified GB model (mbondi2 radii) - **Recommended**
  - `7`: GBneck model
  - `8`: GBneck2 model

### Decomposition

- **`decomp`**: Enable per-residue decomposition (default: `False`)
  - Provides energy contributions from individual residues
  - Significantly increases computation time
  - Useful for identifying hot spots

## Output Format

### Successful Result

```python
{
    "binding_free_energy": -25.3,  # kcal/mol (lower is better)
    "delta_h": -30.5,  # Enthalpy (kcal/mol)
    "delta_s": -5.2,  # Entropy contribution (-T*ΔS, kcal/mol)
    "components": {
        "van_der_waals": -45.2,  # kcal/mol
        "electrostatic": -12.3,  # kcal/mol
        "polar_solvation": 28.5,  # kcal/mol
        "sasa": -4.3  # kcal/mol (non-polar solvation)
    },
    "per_residue_decomposition": [],  # If decomp=True
    "warnings": [],
    "method": "PB",  # or "GB"
    "frames_analyzed": {
        "start": 1,
        "end": -1,
        "interval": 1
    },
    "receptor_group": "Protein",
    "ligand_group": "Ligand"
}
```

### Understanding the Output

- **`binding_free_energy`**: Total binding free energy (ΔG)
  - More negative = stronger binding
  - Typical range: -5 to -15 kcal/mol for good binders
  - Strong binders: < -10 kcal/mol

- **`delta_h`**: Enthalpy change
  - Favorable interactions: negative values

- **`delta_s`**: Entropy contribution (-T*ΔS)
  - Entropy loss upon binding: positive values
  - Entropy gain upon binding: negative values

- **Energy Components**:
  - `van_der_waals`: Shape complementarity and hydrophobic interactions
  - `electrostatic`: Charge-charge interactions
  - `polar_solvation`: Desolvation penalty (usually positive)
  - `sasa`: Non-polar solvation (surface area dependent)

## Workflow Integration

### Typical PharmForge Pipeline

```python
from adapters.vina import VinaAdapter
from adapters.gmx_mmpbsa import GmxMMPBSAAdapter

# Step 1: Initial docking with Vina
vina = VinaAdapter(config={"receptor_path": "receptor.pdbqt"})
docking_result = await vina.execute(smiles)

# Step 2: Convert docked pose to MD input
# (user prepares MD system, runs simulation)

# Step 3: Calculate binding free energy with MM/PBSA
mmpbsa = GmxMMPBSAAdapter()
energy_result = await mmpbsa.execute({
    "topology_file": "md_system.prmtop",
    "trajectory_file": "production.xtc",
    "index_file": "index.ndx",
    "receptor_group": "Protein",
    "ligand_group": "Ligand"
})

print(f"Vina Score: {docking_result.data['binding_affinity']:.2f} kcal/mol")
print(f"MM/PBSA ΔG: {energy_result.data['binding_free_energy']:.2f} kcal/mol")
```

## Troubleshooting

### Error: "gmx_MMPBSA not installed"

**Solution**: Install using conda and pip:
```bash
conda install -c conda-forge ambertools
pip install gmx_MMPBSA
```

### Error: "File not found"

**Solution**: Ensure all file paths are absolute or relative to the working directory:
```python
import os
input_data = {
    "topology_file": os.path.abspath("complex.prmtop"),
    "trajectory_file": os.path.abspath("traj.xtc"),
    "index_file": os.path.abspath("index.ndx"),
    # ...
}
```

### Error: "Group not found in index file"

**Solution**: Verify group names using:
```bash
gmx make_ndx -f complex.gro -n index.ndx
# Lists all available groups
```

Ensure `receptor_group` and `ligand_group` match exactly (case-sensitive).

### Warning: "Binding free energy not found in output"

**Solution**: Check gmx_MMPBSA output for errors:
- Insufficient frames (need at least ~50-100 frames)
- Incompatible topology/trajectory formats
- Missing force field parameters

### Calculation Taking Too Long

**Solutions**:
1. Increase frame interval:
   ```python
   config = {"interval": 10}  # Analyze every 10th frame
   ```

2. Reduce trajectory length:
   ```python
   config = {"startframe": 100, "endframe": 200}
   ```

3. Use MM/GBSA instead of MM/PBSA:
   ```python
   config = {"method": "gb"}
   ```

## Performance Considerations

### Computational Cost

- **MM/PBSA**: 1-30 minutes for 100 frames
  - ~1-5 seconds per frame
  - Scales linearly with number of frames

- **MM/GBSA**: 10x faster than MM/PBSA
  - ~0.1-0.5 seconds per frame
  - Recommended for initial screening

### Memory Usage

- Typical: 1-4 GB RAM
- Increases with system size and decomposition

### Optimization Tips

1. **Use fewer frames**: Analyze every 5-10th frame
2. **Skip equilibration**: Set `startframe` to skip equilibration
3. **Use GB for screening**: Switch to PB for final validation
4. **Disable decomposition**: Only enable when needed

## Scientific Background

### MM/PBSA Theory

The binding free energy is calculated as:

```
ΔG_bind = ΔH - TΔS
        = ΔE_MM + ΔG_solv - TΔS
```

Where:
- `ΔE_MM` = Van der Waals + Electrostatic (molecular mechanics)
- `ΔG_solv` = Polar solvation (PB/GB) + Non-polar solvation (SASA)
- `TΔS` = Entropy contribution (often approximated or omitted)

### When to Use MM/PBSA vs MM/GBSA

**Use MM/PBSA when**:
- High accuracy is required
- Final validation of lead candidates
- Publishing computational predictions
- System has complex electrostatics

**Use MM/GBSA when**:
- Screening many compounds
- Quick initial validation
- Limited computational resources
- Relative ranking is sufficient

### Accuracy and Limitations

**Strengths**:
- More accurate than docking scores
- Captures dynamic effects from MD
- Provides energy component breakdown

**Limitations**:
- Entropy calculations are approximate
- Sensitive to force field choice
- Requires equilibrated MD trajectory
- May not capture specific solvation effects

**Typical Accuracy**:
- Correlation with experiment: R² ~ 0.5-0.8
- MAE: 1-3 kcal/mol
- Better for relative rankings than absolute values

## References

1. **gmx_MMPBSA Paper**:
   Valdes-Tresanco et al., "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS"
   *J. Chem. Theory Comput.* 2021, 17(10), 6281-6291
   DOI: [10.1021/acs.jctc.1c00645](https://doi.org/10.1021/acs.jctc.1c00645)

2. **MM/PBSA Review**:
   Genheden & Ryde, "The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities"
   *Expert Opin. Drug Discov.* 2015, 10(5), 449-461

3. **MMPBSA.py Original**:
   Miller et al., "MMPBSA.py: An Efficient Program for End-State Free Energy Calculations"
   *J. Chem. Theory Comput.* 2012, 8(9), 3314-3321

## Links

- **GitHub**: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
- **Documentation**: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/
- **Tutorial**: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/examples/
- **GROMACS**: https://www.gromacs.org/
- **AmberTools**: https://ambermd.org/AmberTools.php

## License

This adapter is part of PharmForge and follows the MIT license. The underlying gmx_MMPBSA tool is GPL-3.0 licensed.

## Support

For issues with this adapter, please open an issue on the PharmForge GitHub repository.

For gmx_MMPBSA-specific questions, consult the [official documentation](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/) or GitHub issues.
