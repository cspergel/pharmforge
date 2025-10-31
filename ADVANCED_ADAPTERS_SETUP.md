# Advanced ML Tool Adapters - Setup and Testing Report

## Overview

This document provides the status, requirements, and setup instructions for three advanced ML tool adapters in PharmForge:

1. **DiffDock** - Diffusion-based molecular docking
2. **SwissTargetPrediction** - Target prediction from molecular structure
3. **OpenMM** - Molecular dynamics simulations

---

## 1. DiffDock Adapter

### Status: ⚠️ REQUIRES COMPLEX SETUP

**Location:** `adapters/diffdock/adapter.py`

### Description
DiffDock is a state-of-the-art diffusion-based molecular docking tool that performs blind docking without requiring binding site specification. It uses a diffusion generative model to predict binding poses with high accuracy.

### Requirements

#### Essential Dependencies
- ✅ **Python 3.8+**
- ❌ **RDKit** - Not currently installed
  - Install: `conda install -c conda-forge rdkit`
- ❌ **DiffDock Installation** - Requires manual setup
  - Clone: `git clone https://github.com/gcorso/DiffDock.git`
  - Install dependencies: `cd DiffDock && pip install -r requirements.txt`
- ❌ **Model Weights** - ~2-3GB download
  - Download from: https://github.com/gcorso/DiffDock/releases

#### Recommended
- **GPU with CUDA** - 5-50x faster than CPU
- **8GB+ RAM**
- **10GB+ disk space** (for models and dependencies)

### Configuration

The adapter requires configuration in the adapter initialization:

```python
config = {
    'diffdock_path': '/path/to/DiffDock',  # Path to DiffDock installation
    'model_checkpoint': '/path/to/weights',  # Optional: custom model weights
    'inference_steps': 20,  # Number of diffusion steps
    'samples_per_complex': 40,  # Number of poses to generate
    'batch_size': 10,
    'use_gpu': True,
    'python_env': 'python'  # Python executable with DiffDock env
}
```

### Setup Steps

1. **Install RDKit:**
   ```bash
   conda install -c conda-forge rdkit
   ```

2. **Clone and Install DiffDock:**
   ```bash
   git clone https://github.com/gcorso/DiffDock.git
   cd DiffDock
   pip install -r requirements.txt
   ```

3. **Download Model Weights:**
   - Visit: https://github.com/gcorso/DiffDock/releases
   - Download the pre-trained model weights
   - Extract to a known location

4. **Configure Adapter:**
   - Update adapter config with paths to DiffDock installation and model weights
   - Test with a simple protein-ligand pair

### Estimated Setup Time
**30-60 minutes** (depending on download speeds and GPU setup)

### Production Recommendations
- **Use Docker** - DiffDock provides Docker images for easier deployment
- **Dedicated GPU Instance** - AWS, GCP, or Azure with CUDA support
- **API Wrapper** - Consider wrapping DiffDock in a microservice for scalability

### Testing
```python
from adapters.diffdock.adapter import DiffDockAdapter

adapter = DiffDockAdapter(config={
    'diffdock_path': '/path/to/DiffDock',
    'use_gpu': True
})

result = await adapter.execute(
    input_data={
        'smiles': 'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'protein_path': '/path/to/protein.pdb'
    }
)
```

---

## 2. SwissTargetPrediction Adapter

### Status: ⚠️ REQUIRES WEB SCRAPING IMPLEMENTATION

**Location:** `adapters/swisstarget/adapter.py`

### Description
SwissTargetPrediction predicts protein targets for small molecules using 2D and 3D similarity to known bioactive compounds. Useful for reverse pharmacology and off-target prediction.

### Issues Found and Fixes Applied

#### Issue 1: HTTP URL Not Accessible
- **Problem:** Adapter used `http://www.swisstargetprediction.ch` which times out
- **Fix:** Updated to `https://www.swisstargetprediction.ch` ✅
- **Status:** Website is now accessible

#### Issue 2: No Public REST API
- **Problem:** SwissTargetPrediction doesn't provide a public REST API
- **Current Implementation:** Adapter returns example data structure
- **Required:** Web scraping implementation to parse HTML results
- **Status:** Adapter structure is valid, but needs full implementation

### Requirements

#### Essential Dependencies
- ✅ **aiohttp** - Already installed (v3.10.5)
- ✅ **Python 3.8+**

#### For Full Implementation (Optional)
- **BeautifulSoup4** or **lxml** - For HTML parsing
  - Install: `pip install beautifulsoup4 lxml`
- **Selenium** (if needed for JavaScript-heavy pages)
  - Install: `pip install selenium`

### Current Functionality

The adapter currently:
- ✅ Validates SMILES input
- ✅ Has correct HTTPS URL
- ✅ Returns example target prediction data
- ⚠️ Does NOT actually query the web service

### Example Data Structure

```python
{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "organism": "Homo sapiens",
    "targets": [
        {
            "target_name": "Cyclooxygenase-2 (COX-2)",
            "uniprot_id": "P35354",
            "target_class": "Enzyme",
            "probability": 0.95,
            "known_actives": 1523,
            "organism": "Homo sapiens"
        }
        # ... more targets
    ],
    "statistics": {
        "total_targets": 5,
        "max_probability": 0.95,
        "mean_probability": 0.75
    }
}
```

### Implementation Options

#### Option 1: Web Scraping (Recommended)
Implement HTML parsing to extract results from SwissTargetPrediction:

```python
async def _predict_targets(self, smiles: str, organism: str):
    # 1. POST SMILES to web form
    # 2. Parse result page or job ID
    # 3. Poll for completion
    # 4. Extract and parse target predictions
    # 5. Return structured data
```

Tools needed:
- `beautifulsoup4` for HTML parsing
- `aiohttp` for async HTTP requests

#### Option 2: Use Example Data (Current)
Continue using example data for demonstration purposes. Good for testing workflows without external dependencies.

#### Option 3: Alternative Services
Consider using similar services with REST APIs:
- **ChEMBL** - Has REST API for target queries
- **PubChem** - BioAssay data via REST API
- **TargetHunter** - Some target prediction tools with APIs

### Setup Steps

For full implementation:

1. **Install parsing libraries:**
   ```bash
   pip install beautifulsoup4 lxml
   ```

2. **Implement HTML parsing:**
   - Study SwissTargetPrediction's web interface
   - Implement form submission
   - Parse result tables
   - Extract target data

3. **Test with known compounds:**
   - Use well-characterized drugs
   - Verify predicted targets match literature

### Estimated Setup Time
- **Current (example data):** 0 minutes - works out of the box ✅
- **Full web scraping:** 4-8 hours of development

### Testing
```python
from adapters.swisstarget.adapter import SwissTargetAdapter

adapter = SwissTargetAdapter(config={
    'organism': 'Homo sapiens',
    'min_probability': 0.5,
    'max_targets': 15
})

result = await adapter.execute('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
print(result.data['targets'])
```

---

## 3. OpenMM Adapter

### Status: ❌ REQUIRES DEPENDENCIES

**Location:** `adapters/openmm/adapter.py`

### Description
OpenMM is a high-performance molecular dynamics simulation toolkit. This adapter performs energy minimization, MD simulations, and molecular property calculations for small molecules.

### Requirements

#### Essential Dependencies
- ❌ **RDKit** - Not currently installed
  - Install: `conda install -c conda-forge rdkit`
- ❌ **OpenMM** - Not currently installed
  - Install: `conda install -c conda-forge openmm` (recommended)
  - Or: `pip install openmm` (CPU only)
- ✅ **NumPy** - Already installed (v1.26.4)

#### Optional Dependencies
- **pdbfixer** - For PDB structure cleanup
  - Install: `conda install -c conda-forge pdbfixer`
- **openmmforcefields** - For GAFF/AMBER force fields
  - Install: `pip install openmmforcefields`
- **openmm-ml** - For machine learning potentials
  - Install: `pip install openmm-ml`

### GPU Support

OpenMM supports GPU acceleration via:
- **CUDA** - For NVIDIA GPUs (5-50x speedup)
- **OpenCL** - For AMD and other GPUs

GPU support is automatic if CUDA/OpenCL drivers are installed.

### Configuration

```python
config = {
    'force_field': 'gaff',  # Options: gaff, amber14, amber99
    'minimize_steps': 1000,
    'minimize_tolerance': 10.0,  # kJ/mol
    'run_md': False,  # Enable MD simulation
    'md_steps': 10000,
    'md_temperature': 300.0,  # Kelvin
    'md_timestep': 2.0,  # femtoseconds
    'platform': None,  # Auto-select (or 'CUDA', 'OpenCL', 'CPU')
}
```

### Setup Steps

1. **Install via Conda (Recommended):**
   ```bash
   conda install -c conda-forge rdkit openmm pdbfixer
   ```

   This installs:
   - RDKit for SMILES to 3D conversion
   - OpenMM with GPU support
   - PDBFixer for structure cleanup

2. **Or Install via Pip (CPU only):**
   ```bash
   pip install rdkit openmm pdbfixer
   ```

3. **Test Installation:**
   ```python
   import openmm
   print(f"OpenMM version: {openmm.version.version}")

   # Check available platforms
   for i in range(openmm.Platform.getNumPlatforms()):
       platform = openmm.Platform.getPlatform(i)
       print(f"Platform {i}: {platform.getName()}")
   ```

4. **Optional - Install Advanced Force Fields:**
   ```bash
   pip install openmmforcefields
   ```

### Estimated Setup Time
**10-20 minutes** (conda installation + testing)

### Features

The OpenMM adapter provides:

1. **SMILES to 3D Structure Conversion**
   - Uses RDKit's ETKDG algorithm
   - UFF pre-optimization

2. **Energy Minimization**
   - Configurable force fields
   - Returns initial/final energies
   - Convergence detection

3. **Molecular Dynamics (Optional)**
   - Langevin integrator
   - Temperature control
   - Trajectory recording
   - RMSD and radius of gyration calculation

4. **Stability Assessment**
   - Combines energy and structural metrics
   - Returns feasibility score (0-1)
   - Classification: high/medium/low stability

### Testing

```python
from adapters.openmm.adapter import OpenMMAdapter

# Basic energy minimization
adapter = OpenMMAdapter(config={
    'force_field': 'gaff',
    'minimize_steps': 1000,
    'run_md': False
})

result = await adapter.execute('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin

print(f"Final energy: {result.data['minimization']['final_energy']} kJ/mol")
print(f"Stability score: {result.data['stability_score']}")
print(f"Feasibility: {result.data['feasibility']}")

# With molecular dynamics
adapter_md = OpenMMAdapter(config={
    'run_md': True,
    'md_steps': 10000,
    'md_temperature': 300.0
})

result_md = await adapter_md.execute('CC(=O)Oc1ccccc1C(=O)O')
print(f"RMSD: {result_md.data['molecular_dynamics']['rmsd']} Å")
print(f"Radius of gyration: {result_md.data['molecular_dynamics']['radius_of_gyration']} Å")
```

### Platform Selection

OpenMM automatically selects the fastest available platform:
1. **CUDA** - Fastest, requires NVIDIA GPU
2. **OpenCL** - Fast, works with AMD/Intel GPUs
3. **CPU** - Slowest, but always available
4. **Reference** - Very slow, for validation only

To check your platform:
```python
import openmm

system = openmm.System()
integrator = openmm.LangevinIntegrator(300, 1.0, 2.0)
context = openmm.Context(system, integrator)
print(f"Using platform: {context.getPlatform().getName()}")
```

---

## Summary

### Adapter Status Matrix

| Adapter | Status | Dependencies OK | API Available | Ready to Use |
|---------|--------|----------------|---------------|--------------|
| **DiffDock** | ⚠️ Needs Setup | ❌ | N/A (Local) | ❌ |
| **SwissTargetPrediction** | ⚠️ Limited | ✅ | ⚠️ Partial | ⚠️ Example Data |
| **OpenMM** | ❌ Missing Deps | ❌ | N/A (Local) | ❌ |

### Immediate Actions Required

#### For Quick Testing (Recommended)
1. **SwissTargetPrediction**: ✅ Works with example data - no action needed
2. **OpenMM**: Install dependencies (10-20 min)
   ```bash
   conda install -c conda-forge rdkit openmm
   ```
3. **DiffDock**: Skip for now or use Docker deployment

#### For Production Use
1. **SwissTargetPrediction**: Implement web scraping (4-8 hours dev time)
2. **OpenMM**: Install with GPU support for performance
3. **DiffDock**: Set up dedicated GPU instance with Docker

### Fixes Applied

1. ✅ **SwissTargetPrediction URL**: Changed from HTTP to HTTPS
   - File: `adapters/swisstarget/adapter.py`
   - Line 53: `BASE_URL = "https://www.swisstargetprediction.ch"`

### Testing Results

Test script created: `test_advanced_adapters.py`

Results saved to: `adapter_test_results.json`

Run tests with:
```bash
python test_advanced_adapters.py
```

---

## Installation Quick Reference

### Minimum Install (for OpenMM + SwissTarget)
```bash
# Via Conda (recommended)
conda install -c conda-forge rdkit openmm

# Via Pip (CPU only)
pip install rdkit openmm
```

### Full Install (all adapters)
```bash
# Core dependencies
conda install -c conda-forge rdkit openmm pdbfixer

# Web scraping (for SwissTarget)
pip install beautifulsoup4 lxml

# DiffDock (requires manual setup)
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
pip install -r requirements.txt
# Download model weights from releases
```

### Docker Alternative (DiffDock)
```bash
docker pull ghcr.io/gcorso/diffdock:latest
# Configure adapter to use Docker container
```

---

## Support and References

### DiffDock
- Paper: https://arxiv.org/abs/2210.01776
- GitHub: https://github.com/gcorso/DiffDock
- Documentation: See GitHub README

### SwissTargetPrediction
- Website: https://www.swisstargetprediction.ch
- Paper: https://academic.oup.com/nar/article/47/W1/W357/5494775
- No official API documentation

### OpenMM
- Website: http://openmm.org/
- Documentation: http://docs.openmm.org/
- Paper: Eastman et al., PLoS Comput. Biol. 2017, 13(7), e1005659
- Conda: https://anaconda.org/conda-forge/openmm

---

*Last Updated: 2025-10-25*
*Test Environment: Windows 10, Python 3.12.7, Anaconda*
