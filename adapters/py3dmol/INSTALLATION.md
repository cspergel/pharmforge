# py3Dmol Adapter Installation Guide

This guide covers installation and setup of the py3Dmol adapter for PharmForge.

## Prerequisites

- Python 3.8 or higher
- pip package manager
- PharmForge backend core (adapters protocol)

## Installation

### 1. Install py3Dmol

```bash
pip install py3Dmol
```

### 2. Install RDKit (Optional but Recommended)

RDKit is required for SMILES support and 3D coordinate generation:

```bash
pip install rdkit
```

**Note**: If RDKit is not installed, the adapter will still work but SMILES input will not be supported.

### 3. Verify Installation

```python
import py3Dmol
print(f"py3Dmol version: {py3Dmol.__version__}")

# Optional: verify RDKit
try:
    from rdkit import Chem
    print("RDKit is installed")
except ImportError:
    print("RDKit is not installed (SMILES support disabled)")
```

## Quick Start

### Basic Usage

```python
from adapters.py3dmol import Py3DmolAdapter

# Initialize adapter
adapter = Py3DmolAdapter()

# Visualize molecule
result = await adapter.execute(
    {"structure": "CCO"},  # Ethanol SMILES
    style="stick",
    width=500,
    height=400
)

if result.success:
    print("HTML:", result.data["html"])
```

### Frontend Integration

The adapter generates HTML/JavaScript that requires:

1. **jQuery**: Include before py3Dmol
2. **3Dmol.js**: Loaded via CDN

```html
<!DOCTYPE html>
<html>
<head>
    <!-- jQuery (required) -->
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
</head>
<body>
    <!-- Insert generated HTML here -->
    <div id="viewer_xyz" style="width: 400px; height: 400px;"></div>

    <!-- 3Dmol.js (included in generated HTML) -->
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <script>
        // Insert generated JavaScript here
    </script>
</body>
</html>
```

## Platform-Specific Notes

### Windows

```bash
pip install py3Dmol
pip install rdkit
```

### macOS

```bash
pip install py3Dmol
pip install rdkit
```

### Linux

```bash
pip install py3Dmol
pip install rdkit
```

## Docker Installation

Add to your Dockerfile:

```dockerfile
# Install py3Dmol and RDKit
RUN pip install py3Dmol rdkit

# Copy adapter
COPY adapters/py3dmol /app/adapters/py3dmol
```

## Conda Installation

If using Conda:

```bash
conda install -c conda-forge rdkit
pip install py3Dmol
```

## Troubleshooting

### Issue: py3Dmol not found

```
ImportError: No module named 'py3Dmol'
```

**Solution**:
```bash
pip install py3Dmol
```

### Issue: RDKit not found (SMILES support)

```
Warning: RDKit not available - some features require it: pip install rdkit
```

**Solution**:
```bash
pip install rdkit
```

If pip installation fails, try conda:
```bash
conda install -c conda-forge rdkit
```

### Issue: Visualization not rendering in browser

**Solution**:
1. Ensure jQuery is loaded before 3Dmol.js
2. Check browser console for JavaScript errors
3. Verify 3Dmol.js CDN is accessible: https://3Dmol.csb.pitt.edu/build/3Dmol-min.js

### Issue: 3D coordinates not generating from SMILES

**Solution**:
- Ensure RDKit is installed
- Try a simpler SMILES string to test
- Check that molecule is valid

## Dependencies

### Required

- **py3Dmol**: 3D molecular visualization library
  - Version: >= 2.0.0
  - Used for generating interactive visualizations

### Optional

- **RDKit**: Cheminformatics toolkit
  - Version: >= 2022.03.1
  - Required for SMILES support
  - Provides 3D coordinate generation

### Frontend Dependencies

- **jQuery**: JavaScript library (loaded via CDN)
  - Version: >= 3.6.0
  - Required for 3Dmol.js

- **3Dmol.js**: JavaScript 3D visualization (loaded via CDN)
  - Version: Latest from CDN
  - Automatically included in generated HTML

## Version Compatibility

| Component | Minimum Version | Tested Version |
|-----------|----------------|----------------|
| Python | 3.8 | 3.11 |
| py3Dmol | 2.0.0 | 2.1.0 |
| RDKit | 2022.03.1 | 2023.09.1 |
| jQuery | 3.6.0 | 3.7.1 |
| 3Dmol.js | Latest | Latest |

## Testing Installation

Run the test suite to verify installation:

```bash
cd adapters/py3dmol
pytest test_adapter.py -v
```

Or run example usage:

```bash
python example_usage.py
```

## Environment Setup

### Development Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install py3Dmol rdkit pytest

# Verify installation
python -c "import py3Dmol; print('Success!')"
```

### Production Environment

```bash
# Install with production dependencies
pip install --no-cache-dir py3Dmol rdkit

# Set environment variables if needed
export PHARMFORGE_ADAPTERS_PATH=/path/to/adapters
```

## Configuration

The adapter can be configured via the config dictionary:

```python
adapter = Py3DmolAdapter()

# Modify configuration
adapter.config["timeout"] = 60
adapter.config["default_style"] = "sphere"
adapter.config["default_width"] = 800
adapter.config["default_height"] = 600
```

## Adapter Registration

Register the adapter with PharmForge:

```python
from backend.core.adapters.protocol import registry
from adapters.py3dmol import Py3DmolAdapter

# Register adapter
adapter = Py3DmolAdapter()
registry.register(adapter)

# Verify registration
print(registry.list_adapters())  # Should include 'py3dmol'
```

## Next Steps

After installation:

1. Read the [README.md](README.md) for detailed usage
2. Run [example_usage.py](example_usage.py) for examples
3. Run [test_adapter.py](test_adapter.py) to verify functionality
4. Integrate with your PharmForge workflow

## Support

For issues or questions:

1. Check the troubleshooting section above
2. Review the test file for examples
3. Check py3Dmol documentation: https://3dmol.csb.pitt.edu/
4. Check RDKit documentation: https://www.rdkit.org/docs/

## License

Part of PharmForge project.
