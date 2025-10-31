# py3Dmol Adapter - Quick Start Guide

Get started with 3D molecular visualization in 5 minutes.

## Installation

```bash
pip install py3Dmol rdkit
```

## Basic Usage

```python
from adapters.py3dmol import Py3DmolAdapter
import asyncio

async def visualize():
    # Initialize adapter
    adapter = Py3DmolAdapter()

    # Visualize aspirin
    result = await adapter.execute(
        {"structure": "CC(=O)Oc1ccccc1C(=O)O"},
        style="stick",
        width=500,
        height=400
    )

    # Save visualization
    if result.success:
        with open("molecule.html", "w") as f:
            f.write(result.data["html"])
        print("Saved to molecule.html")

# Run
asyncio.run(visualize())
```

## Input Formats

### SMILES String
```python
result = await adapter.execute({"structure": "CCO"})
```

### PDB File
```python
result = await adapter.execute({"structure": "/path/to/protein.pdb"})
```

### PDB String
```python
pdb = "ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00  0.00           C"
result = await adapter.execute({"structure": pdb})
```

## Visual Styles

```python
# Stick representation (default)
result = await adapter.execute({"structure": smiles}, style="stick")

# Sphere representation
result = await adapter.execute({"structure": smiles}, style="sphere")

# Cartoon representation (proteins)
result = await adapter.execute({"structure": pdb_file}, style="cartoon")

# Surface representation
result = await adapter.execute({"structure": smiles}, style="surface")
```

## Customize Dimensions

```python
result = await adapter.execute(
    {"structure": smiles},
    width=800,
    height=600,
    background_color="#f0f0f0"
)
```

## Output Structure

```python
{
    "html": "<div>...</div>",        # Embeddable HTML
    "javascript": "...",              # JavaScript code
    "viewer_id": "viewer_abc123",    # Unique ID
    "structure_type": "ligand",      # protein/ligand/complex
    "structure_format": "pdb",       # Format used
    "style": "stick",                # Applied style
    "width": 500,
    "height": 400
}
```

## Common Use Cases

### 1. Visualize Drug Molecule

```python
drugs = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

for name, smiles in drugs.items():
    result = await adapter.execute({"structure": smiles})
    with open(f"{name}.html", "w") as f:
        f.write(result.data["html"])
```

### 2. Visualize Protein

```python
result = await adapter.execute(
    {"structure": "protein.pdb"},
    style="cartoon",
    width=800,
    height=600
)
```

### 3. Visualize Protein-Ligand Complex

```python
# Adapter automatically detects protein + ligand
result = await adapter.execute(
    {"structure": "complex.pdb"},
    width=700,
    height=700
)
# Protein shown as cartoon, ligand as stick
```

## Frontend Integration

### Simple HTML

```html
<!DOCTYPE html>
<html>
<head>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
</head>
<body>
    <!-- Paste generated HTML here -->
</body>
</html>
```

### With Python Flask

```python
from flask import Flask, render_template

app = Flask(__name__)

@app.route('/molecule/<smiles>')
async def show_molecule(smiles):
    adapter = Py3DmolAdapter()
    result = await adapter.execute({"structure": smiles})
    return render_template('molecule.html',
                         visualization=result.data['html'])
```

### With React

```jsx
function MoleculeViewer({ htmlContent }) {
    return (
        <div dangerouslySetInnerHTML={{ __html: htmlContent }} />
    );
}
```

## Error Handling

```python
result = await adapter.execute({"structure": smiles})

if result.success:
    print("Success:", result.data["structure_type"])
else:
    print("Error:", result.error)
```

## Tips

1. **SMILES Support**: Requires RDKit installation
2. **Large Proteins**: May take a few seconds to render
3. **Multiple Structures**: Combine PDB strings with newlines
4. **Browser Compatibility**: Works in all modern browsers
5. **Interactive**: Viewer supports rotation, zoom, and clicking

## Examples

See [example_usage.py](example_usage.py) for more examples:

```bash
python example_usage.py
```

## Next Steps

- Read [README.md](README.md) for detailed documentation
- Check [INSTALLATION.md](INSTALLATION.md) for setup help
- Run [test_adapter.py](test_adapter.py) to verify installation

## Common Issues

### "py3Dmol not installed"
```bash
pip install py3Dmol
```

### "RDKit not available"
```bash
pip install rdkit
```

### "Viewer not rendering"
- Check jQuery is loaded
- Check browser console
- Verify CDN access

## Configuration

```python
adapter.config["default_style"] = "sphere"
adapter.config["default_width"] = 800
adapter.config["default_height"] = 600
```

## Adapter Info

```python
print(adapter.name)           # "py3dmol"
print(adapter.adapter_type)   # "local"
print(adapter.version)        # "1.0.0"
```

## Quick Reference

| Feature | Parameter | Values |
|---------|-----------|--------|
| Input | structure | SMILES, PDB file, PDB string |
| Style | style | stick, sphere, cartoon, surface, line, cross |
| Size | width, height | Pixels (default: 400x400) |
| Background | background_color | CSS color (default: "white") |
| Labels | show_labels | Boolean (default: False) |

## Performance

- **Small molecules**: < 1 second
- **Proteins**: 1-5 seconds
- **Large complexes**: 5-10 seconds
- **Rendering**: Browser-side (fast)

## Support

For help:
1. Check [README.md](README.md)
2. Review [example_usage.py](example_usage.py)
3. Run tests: `pytest test_adapter.py -v`
4. Check py3Dmol docs: https://3dmol.csb.pitt.edu/
