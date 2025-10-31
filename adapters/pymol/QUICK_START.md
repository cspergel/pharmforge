# PyMOL Adapter - Quick Start Guide

Get started with publication-quality molecular visualization in 5 minutes.

## Installation

```bash
# Install PyMOL open-source
pip install pymol-open-source

# Install RDKit (required for SMILES support)
pip install rdkit
```

## Basic Usage

### 1. Simple Visualization

```python
from adapters.pymol import PyMOLAdapter

adapter = PyMOLAdapter()

# Visualize aspirin
result = await adapter.execute(
    {"structure": "CC(=O)Oc1ccccc1C(=O)O"}
)

print(f"Image saved: {result.data['image_path']}")
```

### 2. Custom Styling

```python
result = await adapter.execute(
    {"structure": "CC(=O)Oc1ccccc1C(=O)O"},
    style="spheres",
    color_scheme="rainbow",
    background="white",
    ray_trace=True
)
```

### 3. High Resolution

```python
result = await adapter.execute(
    {"structure": "CC(=O)Oc1ccccc1C(=O)O"},
    width=1920,
    height=1080,
    ray_trace=True,
    output_path="aspirin_hd.png"
)
```

## Available Options

### Styles
- `sticks` - Ball and stick (default)
- `spheres` - Space-filling
- `cartoon` - Cartoon (proteins)
- `surface` - Molecular surface
- `lines` - Wireframe
- `ribbon` - Ribbon (proteins)
- `mesh` - Mesh surface
- `dots` - Dot surface

### Color Schemes
- `by_element` - Atomic colors (default)
- `by_residue` - Rainbow by residue
- `rainbow` - Rainbow by count
- `spectrum` - By B-factor
- `white`, `gray`, `cyan` - Single colors

### Backgrounds
- `white` (default)
- `black`
- `gray`
- `transparent` (if supported)

## Common Use Cases

### Drug Molecule Gallery

```python
drugs = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

for name, smiles in drugs.items():
    result = await adapter.execute(
        {"structure": smiles},
        output_path=f"{name}.png"
    )
```

### Protein Visualization

```python
# Load PDB file
result = await adapter.execute(
    {"structure": "protein.pdb"},
    style="cartoon",
    color_scheme="rainbow"
)
```

### Publication Figure

```python
# High-quality render for papers
result = await adapter.execute(
    {"structure": "CC(=O)Oc1ccccc1C(=O)O"},
    style="sticks",
    color_scheme="by_element",
    background="white",
    width=1200,
    height=900,
    ray_trace=True,
    output_path="figure_1a.png"
)
```

## Output

All renders return:

```python
{
    "image_path": "/path/to/image.png",
    "width": 800,
    "height": 600,
    "format": "png",
    "file_size": 245678,
    "render_time": 2.3,
    "num_atoms": 21,
    "style": "sticks",
    "color_scheme": "by_element",
    "ray_traced": True
}
```

## Tips

### Fast Previews
```python
# Disable ray tracing for speed
result = await adapter.execute(
    {"structure": smiles},
    ray_trace=False
)
```

### Batch Processing
```python
structures = [
    {"structure": smiles1},
    {"structure": smiles2},
    {"structure": smiles3}
]

results = adapter.render_batch(structures, ray_trace=True)
```

### Error Handling
```python
result = await adapter.execute({"structure": smiles})

if result.success:
    print(f"Success: {result.data['image_path']}")
else:
    print(f"Error: {result.error}")
```

## Performance

| Resolution | Ray Trace | Time |
|-----------|-----------|------|
| 800x600   | No        | ~1s  |
| 800x600   | Yes       | ~3s  |
| 1920x1080 | Yes       | ~8s  |

## Need Help?

- **Full docs:** See `README.md`
- **Examples:** Run `example_usage.py`
- **Tests:** `pytest backend/tests/test_pymol_adapter.py`

## Next Steps

1. Try the examples in `example_usage.py`
2. Read the full documentation in `README.md`
3. Integrate into your PharmForge pipeline
4. Create publication-quality figures!
