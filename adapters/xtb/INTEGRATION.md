# xTB Adapter Integration Guide

## Overview

The xTB adapter has been successfully integrated into PharmForge and is ready for use in quantum chemistry calculations.

## Registration Status

✅ **Registered in adapter registry** (`backend/core/adapter_registry.py`)
- Adapter name: `xtb`
- Display name: "xTB (Quantum Chemistry)"
- Category: ADMET & Properties
- Type: `local`

## Installation Requirements

The xTB adapter requires the `xtb-python` package, which is **only available via conda**:

```bash
# Install xtb-python (REQUIRED)
conda install -c conda-forge xtb-python

# Also ensure RDKit is installed (for SMILES handling)
conda install -c conda-forge rdkit
```

**Note:** These packages are NOT available via pip. You must use conda or mamba.

## API Usage

### Via Adapter Registry

```python
from backend.core.adapter_registry import registry

# Get the adapter
xtb = registry.get("xtb")

# Run calculation
result = await xtb.execute("CCO")

if result.success:
    print(f"Energy: {result.data['energy_kcal_mol']:.2f} kcal/mol")
```

### Direct Import

```python
from adapters.xtb.adapter import xTBAdapter

xtb = xTBAdapter()
result = await xtb.execute("CCO")
```

## Pipeline Integration

The xTB adapter can be integrated into PharmForge pipelines:

```python
from backend.core.pipeline import Pipeline

# Create pipeline
pipeline = Pipeline()

# Add quantum chemistry step
pipeline.add_step("xtb", {
    "adapter": "xtb",
    "method": "GFN2-xTB",
    "optimize": True,
    "solvent": "water"
})

# Execute pipeline
results = await pipeline.execute("CCO")
```

## REST API Endpoints

Once registered, the adapter is accessible via PharmForge's REST API:

### Execute xTB Calculation

```http
POST /api/adapters/xtb/execute
Content-Type: application/json

{
  "smiles": "CCO",
  "method": "GFN2-xTB",
  "optimize": true,
  "solvent": "water"
}
```

Response:
```json
{
  "success": true,
  "data": {
    "smiles": "CCO",
    "method": "GFN2-xTB",
    "energy_kcal_mol": -12345.67,
    "homo_energy_eV": -8.5,
    "lumo_energy_eV": 2.3,
    "homo_lumo_gap_eV": 10.8,
    "dipole_moment_debye": 1.69,
    "converged": true
  },
  "cache_hit": false,
  "metadata": {
    "adapter_name": "xtb",
    "cache_key": "abc123...",
    "version": "1.0.0"
  }
}
```

### Get Adapter Metadata

```http
GET /api/adapters/xtb/metadata
```

Response:
```json
{
  "name": "xtb",
  "type": "local",
  "version": "1.0.0",
  "enabled": true,
  "description": "Fast semiempirical quantum chemistry calculations using xTB",
  "methods": ["GFN0-xTB", "GFN1-xTB", "GFN2-xTB"],
  "capabilities": [
    "Geometry optimization",
    "Energy calculations",
    "HOMO/LUMO analysis",
    "Molecular properties",
    "Solvent models (GBSA)"
  ]
}
```

## Use Cases in Drug Discovery

### 1. Electronic Property Screening

Screen compound libraries for HOMO-LUMO gaps to identify reactive vs. stable compounds:

```python
compounds = ["CCO", "c1ccccc1", "CC(=O)O"]

for smiles in compounds:
    result = await xtb.execute(smiles)
    if result.success:
        gap = result.data['homo_lumo_gap_eV']
        if gap > 5.0:
            print(f"{smiles}: Stable (gap={gap:.2f} eV)")
```

### 2. Geometry Optimization

Optimize molecular geometries for docking or other calculations:

```python
# Optimize structure
result = await xtb.execute("CCO", optimize=True)

# Get optimized coordinates for docking
coords = result.data['optimized_coordinates']
```

### 3. Solvation Free Energy

Calculate solvation effects for ADME predictions:

```python
# Gas phase
result_gas = await xtb.execute("CCO", solvent=None)

# Aqueous phase
result_water = await xtb.execute("CCO", solvent="water")

# Solvation free energy
delta_G_solv = (result_water.data['energy_kcal_mol'] -
                result_gas.data['energy_kcal_mol'])
```

### 4. Charged Species

Calculate properties of ionized forms (important for pH-dependent properties):

```python
# Neutral form
result_neutral = await xtb.execute("CC(=O)O", charge=0)

# Anionic form (deprotonated)
result_anion = await xtb.execute("CC(=O)[O-]", charge=-1)
```

## Performance Characteristics

| Molecule Size | Typical Runtime (GFN2-xTB) |
|--------------|----------------------------|
| 10-20 atoms  | 2-3 seconds |
| 20-50 atoms  | 3-7 seconds |
| 50-100 atoms | 7-15 seconds |

**Optimization:** Use GFN1-xTB for large molecules or high-throughput screening.

## Caching Behavior

The adapter implements automatic caching:
- Cache key based on: SMILES (canonicalized), method, charge, solvent
- First calculation: Full computation (~5 seconds)
- Subsequent identical calculations: Retrieved from cache (<0.01 seconds)
- Cache is persistent across API calls

## Error Handling

Common error scenarios:

1. **xtb-python not installed:**
   ```json
   {
     "success": false,
     "error": "xtb-python not installed. Install with: conda install -c conda-forge xtb-python"
   }
   ```

2. **Invalid SMILES:**
   ```json
   {
     "success": false,
     "error": "Invalid SMILES string"
   }
   ```

3. **Convergence failure:**
   ```json
   {
     "success": true,
     "data": {
       "converged": false,
       ...
     }
   }
   ```

## Testing

Run the test suite:

```bash
cd adapters/xtb
python test_adapter.py
```

This will test:
- Adapter initialization
- Input validation
- Cache key generation
- Metadata
- Execution (if xtb-python is available)

## Integration Checklist

- [x] Adapter implements `AdapterProtocol`
- [x] Registered in `backend/core/adapter_registry.py`
- [x] Input validation implemented
- [x] Cache key generation implemented
- [x] Error handling implemented
- [x] Documentation complete (README, QUICK_START)
- [x] Example usage provided
- [x] Test suite created
- [x] Integration guide written

## Dependencies

```txt
# Required
xtb-python>=0.1.0  # conda-forge only
rdkit>=2023.0.0    # conda-forge recommended

# Optional (for testing)
pytest>=7.0.0
pytest-asyncio>=0.20.0
```

## Known Limitations

1. **Conda-only installation:** xtb-python is not available on PyPI
2. **Organic molecules:** Best for C, H, N, O, S, P, halogens
3. **Size limit:** Practical limit ~500 atoms
4. **Transition metals:** Limited support for metal complexes

## Future Enhancements

Potential additions:
- [ ] Frequency calculations (IR spectra)
- [ ] Thermochemistry (G, H, S)
- [ ] Excited state calculations
- [ ] Metal complexes support (if needed)
- [ ] Batch optimization for multiple conformers

## Support

For issues or questions:
1. Check [README.md](README.md) for usage examples
2. Check [QUICK_START.md](QUICK_START.md) for common patterns
3. Run test suite to verify installation
4. Check PharmForge documentation

## Citation

If using xTB in publications:

```
Grimme, S., Bannwarth, C., & Shushkov, P. (2017).
A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures,
Vibrational Frequencies, and Noncovalent Interactions of Large Molecular
Systems Parametrized for All spd-Block Elements (Z = 1–86).
Journal of Chemical Theory and Computation, 13(5), 1989-2009.
DOI: 10.1021/acs.jctc.7b00118
```

## License

The xTB adapter code follows PharmForge's MIT license. The underlying xtb program is licensed under LGPL-3.0.
