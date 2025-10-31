# gmx_MMPBSA Adapter - Quick Start Guide

## 5-Minute Setup

### 1. Install Dependencies

```bash
# Install AmberTools (required)
conda install -c conda-forge ambertools

# Install gmx_MMPBSA
pip install gmx_MMPBSA

# Verify installation
gmx_MMPBSA --version
```

### 2. Basic Usage

```python
from adapters.gmx_mmpbsa import GmxMMPBSAAdapter

# Initialize adapter
adapter = GmxMMPBSAAdapter()

# Prepare input
input_data = {
    "topology_file": "complex.prmtop",      # Your topology file
    "trajectory_file": "production.xtc",    # Your MD trajectory
    "index_file": "index.ndx",              # GROMACS index file
    "receptor_group": "Protein",            # From index file
    "ligand_group": "Ligand"                # From index file
}

# Run calculation
result = await adapter.execute(input_data)

# Get results
if result.success:
    print(f"ΔG = {result.data['binding_free_energy']:.2f} kcal/mol")
else:
    print(f"Error: {result.error}")
```

### 3. Fast Screening Mode

```python
# Use MM/GBSA for quick screening
adapter = GmxMMPBSAAdapter(
    config={
        "method": "gb",      # Fast Generalized Born
        "interval": 10,      # Every 10th frame
    }
)
```

## Common Use Cases

### Case 1: Validate Docking Results
```python
# After Vina docking, validate with MM/PBSA
mmpbsa = GmxMMPBSAAdapter(config={"method": "pb"})
result = await mmpbsa.execute(input_data)
```

### Case 2: Screen Multiple Compounds
```python
# Use fast GB method
adapter = GmxMMPBSAAdapter(config={"method": "gb", "interval": 10})

for compound in compounds:
    result = await adapter.execute(compound["input"])
    print(f"{compound['name']}: {result.data['binding_free_energy']:.2f}")
```

### Case 3: Identify Hotspots
```python
# Enable per-residue decomposition
adapter = GmxMMPBSAAdapter(config={"decomp": True})
result = await adapter.execute(input_data)

# Top contributing residues in result.data['per_residue_decomposition']
```

## Troubleshooting

### "gmx_MMPBSA not installed"
```bash
conda install -c conda-forge ambertools
pip install gmx_MMPBSA
```

### "File not found"
```python
import os
input_data = {
    "topology_file": os.path.abspath("complex.prmtop"),  # Use absolute paths
    # ...
}
```

### "Group not found in index file"
```bash
# Check available groups
gmx make_ndx -f complex.gro -n index.ndx
# Verify group names match exactly (case-sensitive)
```

## Next Steps

- Read full documentation: `adapters/gmx_mmpbsa/README.md`
- See examples: `adapters/gmx_mmpbsa/example_usage.py`
- Integration guide: `adapters/gmx_mmpbsa/integration_example.py`

## Quick Reference

### MM/PBSA (Accurate)
```python
{"method": "pb", "saltcon": 0.15}
```

### MM/GBSA (Fast)
```python
{"method": "gb", "igb": 5, "interval": 10}
```

### Per-Residue Decomposition
```python
{"method": "gb", "decomp": True}
```

### Frame Selection
```python
{"startframe": 10, "endframe": 100, "interval": 5}
```

---

**Ready to use!** 🚀

For detailed documentation, see `README.md` in this directory.
