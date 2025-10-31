# AutoDock Vina - Quick Start Guide

## 🎯 Status: READY FOR USE

All components are installed, configured, and tested. The Vina docking adapter is operational and integrated into the PharmForge platform.

---

## ✅ What's Been Done

1. **Installed Software**:
   - AutoDock Vina 1.2.7
   - OpenBabel 3.1.1
   - Python bindings

2. **Configured Receptor**:
   - 1HSG HIV-1 Protease receptor prepared and validated
   - Docking box optimized for binding site

3. **Fixed Adapter Code**:
   - Updated for Vina 1.2.7 compatibility
   - Integrated OpenBabel for proper PDBQT conversion
   - Validated score normalization

4. **Tested Thoroughly**:
   - 5/5 validation checks passed
   - HIV protease inhibitor shows expected strong binding (-9.43 kcal/mol)
   - Controls show weak binding as expected

---

## 🚀 Quick Usage

### Basic Example

```python
from adapters.vina.adapter import VinaAdapter

# Configure adapter
config = {
    'receptor_path': '/app/receptors/1HSG_receptor.pdbqt',
    'center_x': 13.073,
    'center_y': 22.467,
    'center_z': 5.557,
    'size_x': 25,
    'size_y': 25,
    'size_z': 25
}

adapter = VinaAdapter(config=config)

# Dock a molecule
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin

if result.success:
    print(f"Binding Affinity: {result.data['binding_affinity']:.2f} kcal/mol")
    print(f"Normalized Score: {result.data['binding_score']:.3f}")
```

### Using Configuration File

```python
import json
from adapters.vina.adapter import VinaAdapter

# Load pre-configured receptor
with open('config/vina_receptor_1hsg.json') as f:
    receptor_config = json.load(f)

adapter = VinaAdapter(config={
    'receptor_path': receptor_config['receptor_path'],
    **receptor_config['docking_box'],
    **receptor_config['vina_parameters']
})
```

---

## 📊 Test Results Summary

| Molecule | Binding Affinity | Normalized Score | Interpretation |
|----------|------------------|------------------|----------------|
| Indinavir (HIV inhibitor) | **-9.43 kcal/mol** | **0.679** | Strong binding ✓ |
| Aspirin | -5.28 kcal/mol | 0.159 | Weak binding |
| Benzene | -4.09 kcal/mol | 0.011 | Very weak |
| Ethanol | -2.34 kcal/mol | 0.000 | No specific binding |

**Key Finding**: HIV protease inhibitor binds 2.3× stronger than benzene control, confirming biological relevance.

---

## 📁 Important Files

### Configuration
- **Receptor config**: `config/vina_receptor_1hsg.json`
- **1HSG receptor**: `/app/receptors/1HSG_receptor.pdbqt` (in Docker)

### Code
- **Adapter**: `adapters/vina/adapter.py`
- **Tests**: `backend/tests/test_vina_adapter.py`
- **Integration test**: `test_vina_real_docking.py`

### Documentation
- **Full documentation**: `VINA_SETUP_DOCUMENTATION.md`
- **This quick start**: `VINA_QUICK_START.md`

---

## 🧪 Running Tests

### Integration Test (Real Docking)

```bash
# Inside Docker container
docker-compose exec backend bash -c "cd /app && python3 test_vina_real_docking.py"
```

### Unit Tests (Mocked)

```bash
docker-compose exec backend pytest backend/tests/test_vina_adapter.py -v
```

---

## 🔧 Receptor Configuration

### Current Receptor: 1HSG (HIV-1 Protease)

```json
{
  "receptor_path": "/app/receptors/1HSG_receptor.pdbqt",
  "center_x": 13.073,
  "center_y": 22.467,
  "center_z": 5.557,
  "size_x": 25,
  "size_y": 25,
  "size_z": 25
}
```

**Best for**: HIV protease inhibitors, drug-like small molecules

---

## 📈 Score Interpretation

### Raw Binding Affinity (kcal/mol)

- **< -10**: Excellent binding (drug-like)
- **-8 to -10**: Good binding
- **-5 to -8**: Moderate binding
- **> -5**: Weak/non-specific binding

### Normalized Score (0-1)

- **0.8-1.0**: Excellent
- **0.6-0.8**: Good
- **0.4-0.6**: Moderate
- **0.2-0.4**: Weak
- **0.0-0.2**: Very weak

**Formula**: Maps -12 kcal/mol → 1.0, -4 kcal/mol → 0.0

---

## ⚙️ Advanced Configuration

### Adjust Search Exhaustiveness

```python
config = {
    'receptor_path': '/app/receptors/1HSG_receptor.pdbqt',
    'exhaustiveness': 16,  # More thorough (slower)
    'num_modes': 20,       # More binding poses
    # ... other config
}
```

### Adjust Docking Box

```python
config = {
    'center_x': 15.0,  # Move box
    'center_y': 20.0,
    'center_z': 10.0,
    'size_x': 30,      # Larger box
    'size_y': 30,
    'size_z': 30,
}
```

---

## 🐛 Troubleshooting

### Problem: "Vina binary not found"

**Solution**:
```bash
docker-compose exec backend apt-get update
docker-compose exec backend apt-get install -y autodock-vina
```

### Problem: "OpenBabel not found"

**Solution**:
```bash
docker-compose exec backend apt-get install -y openbabel python3-openbabel
```

### Problem: Unrealistic binding scores

**Check**:
1. Is the receptor properly prepared?
2. Is the docking box centered on the binding site?
3. Is the SMILES valid?

---

## 📚 Learn More

- **Full Documentation**: See `VINA_SETUP_DOCUMENTATION.md`
- **AutoDock Vina**: https://vina.scripps.edu/
- **1HSG Structure**: https://www.rcsb.org/structure/1HSG

---

## ✨ Next Steps

1. **Test with your molecules**: Replace SMILES in example code
2. **Try different receptors**: Download more PDB structures
3. **Integrate with pipeline**: Use in multi-adapter workflows
4. **Optimize parameters**: Adjust exhaustiveness for speed/accuracy trade-off

---

**Last Updated**: October 25, 2025
**Status**: ✓ Production Ready
**Version**: Vina 1.2.7, OpenBabel 3.1.1
