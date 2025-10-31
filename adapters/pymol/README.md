# PyMOL Adapter

Publication-quality molecular visualization adapter for PharmForge using PyMOL.

## Overview

The PyMOL adapter generates high-quality molecular images with customizable styles, colors, and rendering options. It supports SMILES strings (via RDKit conversion), PDB files, and MOL files, and is designed for headless/server deployment without requiring a GUI.

## Features

- **Multiple Input Formats**: SMILES, PDB, MOL, SDF
- **Visualization Styles**: sticks, spheres, cartoon, surface, lines, ribbon, mesh, dots
- **Color Schemes**: by_element, by_residue, rainbow, spectrum, white, gray, cyan
- **Ray Tracing**: High-quality publication-ready renders
- **Batch Processing**: Render multiple molecules efficiently
- **Headless Operation**: Works on servers without GUI
- **Automatic Cleanup**: Manages temporary files automatically

## Installation

### PyMOL Open Source

```bash
# Install PyMOL open-source version
pip install pymol-open-source

# Also requires RDKit for SMILES support
pip install rdkit
```

### System Dependencies (Linux)

```bash
# Ubuntu/Debian
sudo apt-get install pymol

# CentOS/RHEL
sudo yum install pymol
```

## Usage

### Basic Usage

```python
from adapters.pymol import PyMOLAdapter

# Initialize adapter
adapter = PyMOLAdapter()

# Render from SMILES
result = await adapter.execute({
    "structure": "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
}, style="sticks", ray_trace=True)

if result.success:
    print(f"Image saved to: {result.data['image_path']}")
    print(f"Render time: {result.data['render_time']}s")
```

### From PDB File

```python
# Load from PDB file
result = await adapter.execute({
    "structure": "/path/to/protein.pdb"
}, style="cartoon", color_scheme="rainbow")
```

### From PDB String

```python
# Use PDB content directly
pdb_content = """
ATOM      1  N   MET A   1      20.154  29.699   5.276  1.00 49.05           N
ATOM      2  CA  MET A   1      21.379  30.420   5.768  1.00 49.05           C
...
"""

result = await adapter.execute({
    "structure": pdb_content
}, format="pdb", style="cartoon")
```

### Advanced Options

```python
result = await adapter.execute(
    {"structure": "CC(=O)Oc1ccccc1C(=O)O"},
    style="sticks",              # Visualization style
    color_scheme="by_element",   # Color scheme
    background="white",          # Background color
    width=1200,                  # Image width
    height=900,                  # Image height
    ray_trace=True,             # High-quality ray tracing
    output_path="/path/to/output.png"  # Custom output path
)
```

## Visualization Styles

### Available Styles

1. **sticks** - Ball-and-stick representation (default)
2. **spheres** - Space-filling spheres (CPK model)
3. **cartoon** - Cartoon representation (proteins)
4. **surface** - Molecular surface
5. **lines** - Thin lines (wireframe)
6. **ribbon** - Ribbon representation (proteins)
7. **mesh** - Mesh surface
8. **dots** - Dot surface

### Style Examples

```python
# Protein visualization
await adapter.execute(
    {"structure": "protein.pdb"},
    style="cartoon",
    color_scheme="rainbow"
)

# Small molecule visualization
await adapter.execute(
    {"structure": "CCO"},
    style="sticks",
    color_scheme="by_element"
)

# Surface representation
await adapter.execute(
    {"structure": "protein.pdb"},
    style="surface",
    color_scheme="white"
)
```

## Color Schemes

1. **by_element** - Color by atomic element (default)
2. **by_residue** - Color by residue (rainbow gradient)
3. **rainbow** - Rainbow gradient by atom count
4. **spectrum** - Color by B-factor
5. **white** - Single white color
6. **gray** - Single gray color
7. **cyan** - Single cyan color

## Batch Processing

Render multiple molecules efficiently:

```python
structures = [
    {"structure": "CC(=O)Oc1ccccc1C(=O)O"},  # Aspirin
    {"structure": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"},  # Ibuprofen
    {"structure": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"}  # Caffeine
]

results = adapter.render_batch(
    structures,
    style="sticks",
    ray_trace=True,
    width=800,
    height=600
)

for i, result in enumerate(results):
    if result.success:
        print(f"Molecule {i+1}: {result.data['image_path']}")
```

## Output Format

The adapter returns an `AdapterResult` with the following data structure:

```python
{
    "image_path": "/path/to/output.png",
    "width": 800,
    "height": 600,
    "format": "png",
    "file_size": 245678,  # bytes
    "render_time": 2.3,   # seconds
    "num_atoms": 21,
    "style": "sticks",
    "color_scheme": "by_element",
    "background": "white",
    "ray_traced": True,
    "warnings": []
}
```

## Configuration

Configure the adapter during initialization:

```python
adapter = PyMOLAdapter()

# Update configuration
adapter.config.update({
    "default_style": "spheres",
    "default_width": 1200,
    "default_height": 900,
    "default_ray_trace": True,
    "temp_dir": "/custom/temp/path",
    "cleanup_temp_files": True
})
```

## Integration with PharmForge Pipeline

Use within a PharmForge workflow:

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()

# Add visualization step
pipeline.add_step(
    name="visualize",
    adapter="pymol",
    params={
        "style": "sticks",
        "color_scheme": "by_element",
        "ray_trace": True,
        "width": 1200,
        "height": 900
    }
)

# Execute pipeline
result = await pipeline.execute("CC(=O)Oc1ccccc1C(=O)O")
```

## Performance Considerations

### Ray Tracing

Ray tracing produces high-quality images but is computationally expensive:

```python
# Fast preview (no ray tracing)
result = await adapter.execute(
    {"structure": smiles},
    ray_trace=False
)

# High-quality render (with ray tracing)
result = await adapter.execute(
    {"structure": smiles},
    ray_trace=True
)
```

### Batch Optimization

For batch processing, consider:

1. Disable ray tracing for previews
2. Reduce image dimensions for thumbnails
3. Use caching for repeated structures

```python
# Generate thumbnails quickly
for structure in structures:
    await adapter.execute(
        {"structure": structure},
        width=400,
        height=300,
        ray_trace=False
    )
```

## Troubleshooting

### PyMOL Not Available

```python
# Check if PyMOL is installed
from adapters.pymol import PyMOLAdapter

adapter = PyMOLAdapter()
if not adapter.enabled:
    print("PyMOL is not installed!")
```

### RDKit Not Available

For SMILES support, RDKit must be installed:

```bash
pip install rdkit
```

### Headless Mode Issues

Ensure PyMOL is running in headless mode:

```python
import pymol
pymol.finish_launching(['pymol', '-c'])  # -c for command-line mode
```

### File Permission Errors

Ensure write permissions to temp directory:

```python
adapter.config["temp_dir"] = "/path/with/write/access"
```

## Examples

### Example 1: Drug Molecule Gallery

Generate a gallery of drug molecules:

```python
drugs = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Paracetamol": "CC(=O)Nc1ccc(cc1)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

for name, smiles in drugs.items():
    result = await adapter.execute(
        {"structure": smiles},
        style="sticks",
        ray_trace=True,
        output_path=f"gallery/{name}.png"
    )
    print(f"{name}: {result.success}")
```

### Example 2: Protein-Ligand Complex

Visualize a protein-ligand complex:

```python
# Load complex PDB
result = await adapter.execute(
    {"structure": "complex.pdb"},
    style="cartoon",
    color_scheme="by_residue",
    ray_trace=True,
    width=1920,
    height=1080
)
```

### Example 3: Style Comparison

Generate images with different styles:

```python
styles = ["sticks", "spheres", "lines", "surface"]
smiles = "CC(=O)Oc1ccccc1C(=O)O"

for style in styles:
    result = await adapter.execute(
        {"structure": smiles},
        style=style,
        output_path=f"comparison_{style}.png"
    )
```

## API Reference

### PyMOLAdapter Methods

#### `__init__()`

Initialize the adapter with default configuration.

#### `async execute(input_data, **kwargs)`

Execute visualization and generate image.

**Parameters:**
- `input_data` (dict): Dictionary with 'structure' key
- `format` (str): Format hint ('smiles', 'pdb', 'mol')
- `style` (str): Visualization style
- `color_scheme` (str): Color scheme
- `background` (str): Background color
- `width` (int): Image width
- `height` (int): Image height
- `ray_trace` (bool): Enable ray tracing
- `output_path` (str): Custom output path

**Returns:**
- `AdapterResult`: Result object with image path and metadata

#### `render_batch(structures, **common_kwargs)`

Render multiple structures in batch.

**Parameters:**
- `structures` (list): List of structure dictionaries
- `**common_kwargs`: Common parameters for all structures

**Returns:**
- `list`: List of AdapterResult objects

## License

Part of PharmForge - MIT License

## Citation

If you use this adapter in your research, please cite:

```bibtex
@software{pharmforge_pymol_adapter,
  title = {PyMOL Adapter for PharmForge},
  author = {PharmForge Team},
  year = {2025},
  url = {https://github.com/your-repo/pharmforge}
}
```

## See Also

- [PyMOL Official Documentation](https://pymol.org/)
- [PyMOL Open Source](https://github.com/schrodinger/pymol-open-source)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [PharmForge Documentation](../../docs/)
