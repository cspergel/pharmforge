# py3Dmol Adapter

Interactive 3D molecular visualization adapter for PharmForge using py3Dmol.

## Overview

The py3Dmol adapter generates embeddable HTML/JavaScript visualizations of molecular structures. It supports multiple input formats (PDB, mol2, SDF, SMILES) and various visual styles for proteins, ligands, and protein-ligand complexes.

## Installation

```bash
pip install py3Dmol
```

For SMILES support (recommended):
```bash
pip install rdkit
```

## Features

- **Multiple Input Formats**: PDB files, PDB strings, mol2, SDF, SMILES
- **Visual Styles**: stick, sphere, cartoon, surface, line, cross
- **Structure Types**: Automatic detection of proteins, ligands, and complexes
- **Interactive Visualization**: Clickable atoms, rotation, zoom
- **Embeddable Output**: HTML/JavaScript ready for frontend integration
- **Smart Styling**: Automatic style selection based on structure type

## Usage

### Basic Usage

```python
from adapters.py3dmol import Py3DmolAdapter

# Initialize adapter
adapter = Py3DmolAdapter()

# Visualize from SMILES
result = await adapter.execute({
    "structure": "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
}, style="stick", width=600, height=400)

# Visualize from PDB file
result = await adapter.execute({
    "structure": "/path/to/protein.pdb"
}, style="cartoon")

# Visualize from PDB string
pdb_string = """
ATOM      1  CA  ALA A   1      -8.778   1.234   0.000  1.00  0.00           C
...
"""
result = await adapter.execute({
    "structure": pdb_string
}, style="stick")
```

### Input Format

```python
input_data = {
    "structure": str  # PDB string, file path, or SMILES
}

kwargs = {
    "style": str,              # Visual style (default: "stick")
    "width": int,              # Viewer width in pixels (default: 400)
    "height": int,             # Viewer height in pixels (default: 400)
    "background_color": str,   # Background color (default: "white")
    "show_labels": bool,       # Show atom labels (default: False)
    "format": str              # Format hint: 'pdb', 'mol2', 'sdf', 'smiles'
}
```

### Output Format

```python
{
    "success": True,
    "data": {
        "html": str,              # Complete HTML with embedded viewer
        "javascript": str,        # JavaScript code for viewer
        "viewer_id": str,         # Unique viewer ID
        "structure_type": str,    # "protein", "ligand", or "complex"
        "structure_format": str,  # "pdb", "mol2", "sdf"
        "style": str,             # Applied style
        "width": int,             # Viewer width
        "height": int             # Viewer height
    },
    "metadata": {
        "source": "py3dmol",
        "adapter_version": "1.0.0",
        "computation_type": "local",
        "structure_type": str,
        "structure_format": str,
        "style": str
    }
}
```

## Visual Styles

### Available Styles

1. **stick** (default): Ball-and-stick representation
2. **sphere**: Space-filling spheres
3. **cartoon**: Cartoon representation (ideal for proteins)
4. **surface**: Molecular surface
5. **line**: Line representation
6. **cross**: Cross representation

### Automatic Style Selection

The adapter intelligently applies styles based on structure type:

- **Proteins**: Cartoon representation with spectrum coloring
- **Ligands**: User-specified style (default: stick)
- **Complexes**: Cartoon for protein + stick for ligand

## Structure Type Detection

The adapter automatically detects:

- **Protein**: Contains standard amino acid residues
- **Ligand**: Small molecule (from SMILES or HETATM records)
- **Complex**: Protein with bound ligand(s)

## Examples

### Example 1: Drug Molecule (SMILES)

```python
result = await adapter.execute({
    "structure": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
}, style="sphere", width=500, height=500, background_color="#f0f0f0")

# Use the HTML in your frontend
html_output = result.data["html"]
```

### Example 2: Protein Structure

```python
result = await adapter.execute({
    "structure": "/path/to/1ubq.pdb"  # Ubiquitin
}, style="cartoon", width=800, height=600)

# Structure type will be automatically detected as "protein"
print(result.data["structure_type"])  # "protein"
```

### Example 3: Protein-Ligand Complex

```python
result = await adapter.execute({
    "structure": "/path/to/complex.pdb"
}, width=700, height=700)

# Automatically applies cartoon to protein, stick to ligand
print(result.data["structure_type"])  # "complex"
```

### Example 4: Multiple Structures

For multiple structures in one view, combine PDB strings:

```python
pdb_combined = protein_pdb + "\n" + ligand_pdb

result = await adapter.execute({
    "structure": pdb_combined
}, style="stick")
```

## Frontend Integration

### Embedding in HTML

```html
<!DOCTYPE html>
<html>
<head>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
</head>
<body>
    <!-- Insert the generated HTML here -->
    <div id="viewer_abc123" style="width: 400px; height: 400px;"></div>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <script>
        // Insert the generated JavaScript here
    </script>
</body>
</html>
```

### React Integration

```jsx
import React, { useEffect } from 'react';

function MoleculeViewer({ htmlContent, javascript }) {
    useEffect(() => {
        // Execute the visualization JavaScript
        eval(javascript);
    }, [javascript]);

    return (
        <div dangerouslySetInnerHTML={{ __html: htmlContent }} />
    );
}
```

## Error Handling

The adapter gracefully handles errors:

```python
result = await adapter.execute({
    "structure": "invalid_smiles"
})

if not result.success:
    print(result.error)  # Error message
    print(result.metadata)  # Additional context
```

## Dependencies

- **py3Dmol**: Required for visualization
- **RDKit**: Optional, required for SMILES support
- **jQuery**: Required in frontend (loaded via CDN)

## Supported Formats

| Format | Extension | Description | SMILES Support |
|--------|-----------|-------------|----------------|
| PDB | .pdb | Protein Data Bank format | With RDKit |
| mol2 | .mol2 | Tripos Mol2 format | With RDKit |
| SDF | .sdf, .mol | Structure Data File | With RDKit |
| SMILES | - | Simplified molecular input | With RDKit |

## Performance Considerations

- Local computation (no API calls)
- Fast for small molecules (< 1 second)
- Moderate for proteins (1-5 seconds)
- SMILES conversion requires 3D coordinate generation
- Visualization rendering happens in browser

## Limitations

- SMILES support requires RDKit installation
- Very large structures (>10,000 atoms) may be slow in browser
- Surface rendering can be computationally intensive
- Generated coordinates from SMILES are not necessarily biologically relevant

## Configuration

```python
adapter.config = {
    "timeout": 30,
    "default_style": "stick",
    "default_width": 400,
    "default_height": 400,
    "supported_styles": ["stick", "sphere", "cartoon", "surface", "line", "cross"],
    "supported_formats": ["pdb", "mol2", "sdf", "smiles"]
}
```

## Troubleshooting

### py3Dmol not installed
```
Error: py3Dmol is not installed. Install with: pip install py3Dmol
```
**Solution**: `pip install py3Dmol`

### SMILES conversion fails
```
Error: RDKit required for SMILES conversion
```
**Solution**: `pip install rdkit`

### Invalid structure
```
Error: Failed to prepare structure data
```
**Solution**: Verify file path exists or SMILES string is valid

### Viewer not rendering
- Ensure jQuery is loaded before py3Dmol
- Check browser console for JavaScript errors
- Verify 3Dmol.js CDN is accessible

## Advanced Features

### Custom Colors

Modify the generated JavaScript to customize colors:

```javascript
viewer.setStyle({stick: {colorscheme: 'greenCarbon', radius: 0.2}});
```

### Adding Labels

```python
result = await adapter.execute({
    "structure": smiles
}, show_labels=True)
```

### Surface Representation

```python
result = await adapter.execute({
    "structure": protein_pdb
}, style="surface", background_color="#000000")
```

## Integration with Other Adapters

### With RDKit Adapter

```python
# Get molecular properties
rdkit_result = await rdkit_adapter.execute(smiles)

# Visualize the molecule
viz_result = await py3dmol_adapter.execute({
    "structure": smiles
}, style="stick")

# Combine results
combined = {
    "properties": rdkit_result.data,
    "visualization": viz_result.data
}
```

### With Docking Adapters

```python
# Run docking
docking_result = await vina_adapter.execute(...)

# Visualize docked complex
viz_result = await py3dmol_adapter.execute({
    "structure": docking_result.data["complex_pdb"]
}, style="cartoon")
```

## API Reference

### Py3DmolAdapter

#### Methods

- `validate_input(input_data)`: Validate input structure
- `execute(input_data, **kwargs)`: Generate visualization
- `generate_cache_key(input_data, **kwargs)`: Generate cache key

#### Properties

- `name`: "py3dmol"
- `adapter_type`: "local"
- `version`: "1.0.0"

## License

Part of PharmForge project.

## References

- [py3Dmol Documentation](https://3dmol.csb.pitt.edu/)
- [3Dmol.js GitHub](https://github.com/3dmol/3Dmol.js)
- [RDKit Documentation](https://www.rdkit.org/docs/)
