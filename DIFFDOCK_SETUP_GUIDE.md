# DiffDock Setup Guide for PharmForge

**Date:** October 26, 2025
**Purpose:** Complete setup guide for DiffDock GPU-accelerated molecular docking

---

## Overview

**DiffDock** is a state-of-the-art diffusion model for molecular docking that significantly outperforms traditional methods like AutoDock Vina on complex binding scenarios.

**Key Features:**
- Deep learning-based pose prediction
- Handles conformational flexibility
- GPU-accelerated for fast inference
- High accuracy on challenging targets

**Requirements:**
- NVIDIA GPU with CUDA support
- ~5GB disk space for models
- Docker with GPU passthrough

---

## Current Status

✅ **Docker GPU Configuration:** Completed
- `docker-compose.yml` updated with GPU support
- `NVIDIA_VISIBLE_DEVICES=all` environment variable added
- GPU deployment resources configured

⚠️ **GPU Availability:** Currently not accessible in container
- Host has RTX 5080 GPU
- Docker needs GPU passthrough enabled

⏳ **DiffDock Installation:** Pending
- Package needs to be installed
- Pre-trained models need to be downloaded (~5GB)

---

## Prerequisites

### 1. Windows Docker Desktop GPU Setup

**For Windows with Docker Desktop:**

1. **Install NVIDIA Container Toolkit (WSL2 required)**

```powershell
# In PowerShell (Administrator)
# Ensure WSL2 is installed
wsl --install

# Update WSL2
wsl --update

# Install NVIDIA drivers for Windows
# Download from: https://www.nvidia.com/Download/index.aspx
```

2. **Enable GPU in Docker Desktop**

- Open Docker Desktop
- Go to Settings → Resources → WSL Integration
- Enable integration for your WSL distribution
- Restart Docker Desktop

3. **Verify GPU Access in WSL2**

```bash
# In WSL2 terminal
nvidia-smi
```

Expected output: GPU information displayed

4. **Test GPU in Docker**

```bash
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
```

---

### 2. Verify Current GPU Configuration

```bash
# Check if GPU is visible in Docker
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"
docker-compose exec -T backend python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}'); print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"None\"}')"
```

**Current Result:** `CUDA available: False`

**Solution:** Follow Prerequisites step 1 above to enable GPU passthrough

---

## Installation Steps

### Step 1: Install DiffDock Package

DiffDock requires several dependencies including PyTorch with CUDA support (already installed).

```bash
# Enter backend container
docker-compose exec -T backend bash

# Install DiffDock and dependencies
pip install e3nn==0.5.1
pip install spyrmsd==0.6.0
pip install prody==2.4.0
pip install biopython==1.81
pip install fair-esm==2.0.0

# Install DiffDock from GitHub
cd /tmp
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
pip install -e .
```

**Alternative (if available on PyPI):**
```bash
docker-compose exec -T backend pip install diffdock
```

---

### Step 2: Download Pre-trained Models

DiffDock uses pre-trained weights (~5GB total).

```bash
# Create model directory
docker-compose exec -T backend mkdir -p /app/diffdock_models

# Download models using Python script
docker-compose exec -T backend python << 'EOF'
import os
import requests
from pathlib import Path
from tqdm import tqdm

def download_file(url, dest_path):
    """Download file with progress bar"""
    print(f"Downloading {dest_path.name}...")
    response = requests.get(url, stream=True)
    response.raise_for_status()

    total_size = int(response.headers.get('content-length', 0))
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    with open(dest_path, 'wb') as f:
        if total_size == 0:
            f.write(response.content)
        else:
            with tqdm(total=total_size, unit='MB', unit_scale=True, unit_divisor=1024*1024) as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
                    pbar.update(len(chunk))

    print(f"✓ Downloaded: {dest_path.name} ({dest_path.stat().st_size / 1024 / 1024:.1f} MB)")

# Model URLs (update with actual URLs from DiffDock repository)
models = {
    'score_model': 'https://github.com/gcorso/DiffDock/releases/download/v1.1/score_model.pt',
    'confidence_model': 'https://github.com/gcorso/DiffDock/releases/download/v1.1/confidence_model.pt',
}

model_dir = Path('/app/diffdock_models')

for name, url in models.items():
    dest = model_dir / f"{name}.pt"
    if dest.exists():
        print(f"EXISTS: {dest.name}")
    else:
        download_file(url, dest)

print("\n✓ DiffDock models downloaded successfully")
print(f"Location: {model_dir}")
EOF
```

**Note:** Model URLs may change. Check the official DiffDock repository for current download links:
https://github.com/gcorso/DiffDock#pre-trained-models

---

### Step 3: Configure DiffDock Adapter

Update adapter configuration to use downloaded models.

**File:** `adapters/diffdock/config.yaml`

```yaml
# DiffDock Adapter Configuration
name: diffdock
version: 1.0.0
type: docking

# Model paths
models:
  score_model: /app/diffdock_models/score_model.pt
  confidence_model: /app/diffdock_models/confidence_model.pt

# GPU configuration
gpu:
  enabled: true
  device: cuda:0  # Use first GPU

# Docking parameters
parameters:
  inference_steps: 20
  samples_per_complex: 40
  batch_size: 10
  confidence_cutoff: 0.0

# Cache settings
cache:
  enabled: true
  ttl: 86400  # 24 hours
```

---

### Step 4: Update Environment Variables

Add DiffDock configuration to `.env`:

```env
# DiffDock Configuration
DIFFDOCK_MODEL_DIR=/app/diffdock_models
DIFFDOCK_USE_GPU=true
DIFFDOCK_DEVICE=cuda:0
DIFFDOCK_BATCH_SIZE=10
```

---

### Step 5: Restart Docker Services

```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"
docker-compose restart backend celery-worker
```

---

## Verification

### Test 1: Verify GPU Access

```bash
docker-compose exec -T backend python -c "
import torch
print(f'PyTorch version: {torch.__version__}')
print(f'CUDA available: {torch.cuda.is_available()}')
if torch.cuda.is_available():
    print(f'CUDA version: {torch.version.cuda}')
    print(f'GPU count: {torch.cuda.device_count()}')
    print(f'GPU name: {torch.cuda.get_device_name(0)}')
    print(f'GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB')
"
```

**Expected Output:**
```
PyTorch version: 2.5.0+cu124
CUDA available: True
CUDA version: 12.4
GPU count: 1
GPU name: NVIDIA GeForce RTX 5080
GPU memory: 16.0 GB
```

---

### Test 2: Test DiffDock Import

```bash
docker-compose exec -T backend python -c "
try:
    import diffdock
    print('✓ DiffDock imported successfully')
    print(f'  Version: {diffdock.__version__ if hasattr(diffdock, \"__version__\") else \"unknown\"}')
except ImportError as e:
    print(f'✗ DiffDock import failed: {e}')
"
```

---

### Test 3: Test DiffDock Adapter

```bash
docker-compose exec -T backend python -c "
from adapters.diffdock.adapter import DiffDockAdapter

adapter = DiffDockAdapter()
print(f'✓ DiffDock adapter initialized')
print(f'  Name: {adapter.name}')
print(f'  Version: {adapter.version}')
print(f'  GPU enabled: {adapter.config.get(\"gpu\", {}).get(\"enabled\", False)}')
"
```

---

### Test 4: Run Sample Docking

```bash
docker-compose exec -T backend python -c "
from adapters.diffdock.adapter import DiffDockAdapter

# Simple docking test
adapter = DiffDockAdapter()

# Example: Dock aspirin (CCO) to a protein
result = adapter.execute(
    protein_pdb='1abc',  # Example PDB ID
    ligand_smiles='CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
    use_gpu=True
)

print(f'✓ Docking completed')
print(f'  Poses generated: {len(result.get(\"poses\", []))}')
print(f'  Best confidence: {result.get(\"best_confidence\", 0):.3f}')
"
```

---

## Performance Benchmarks

### Expected Performance (RTX 5080)

| Metric | Value |
|--------|-------|
| **Single complex docking** | ~10-30 seconds |
| **40 poses generation** | ~2-5 minutes |
| **Batch (10 complexes)** | ~10-20 minutes |
| **GPU Memory Usage** | ~4-8 GB |

### Comparison to Vina

| Method | Speed | Accuracy (Top-1) | Accuracy (Top-5) |
|--------|-------|------------------|------------------|
| **AutoDock Vina** | Fast (1-5 min) | ~30-40% | ~50-60% |
| **DiffDock (CPU)** | Slow (30+ min) | ~50-60% | ~70-80% |
| **DiffDock (GPU)** | Fast (2-5 min) | ~50-60% | ~70-80% |

**Recommendation:** Use DiffDock for challenging targets, Vina for high-throughput screening

---

## Troubleshooting

### GPU Not Detected in Container

**Problem:** `torch.cuda.is_available()` returns `False`

**Solutions:**

1. **Enable GPU in Docker Desktop** (Windows)
   - Settings → Resources → WSL Integration
   - Enable for your distribution
   - Restart Docker Desktop

2. **Verify nvidia-docker2 installed** (Linux)
   ```bash
   docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
   ```

3. **Check docker-compose.yml GPU configuration**
   ```yaml
   deploy:
     resources:
       reservations:
         devices:
           - driver: nvidia
             count: all
             capabilities: [gpu]
   ```

4. **Rebuild containers**
   ```bash
   docker-compose down
   docker-compose build --no-cache backend
   docker-compose up -d
   ```

---

### Out of Memory Errors

**Problem:** `CUDA out of memory` error during docking

**Solutions:**

1. **Reduce batch size**
   ```env
   DIFFDOCK_BATCH_SIZE=5  # Default: 10
   ```

2. **Reduce samples per complex**
   ```python
   result = adapter.execute(
       protein_pdb='1abc',
       ligand_smiles='CCO',
       samples_per_complex=20  # Default: 40
   )
   ```

3. **Clear GPU cache**
   ```python
   import torch
   torch.cuda.empty_cache()
   ```

---

### Model Download Fails

**Problem:** Cannot download pre-trained models

**Solutions:**

1. **Manual download**
   - Visit: https://github.com/gcorso/DiffDock/releases
   - Download model files manually
   - Copy to `/app/diffdock_models/` in container

2. **Use wget inside container**
   ```bash
   docker-compose exec backend bash
   cd /app/diffdock_models
   wget https://github.com/gcorso/DiffDock/releases/download/v1.1/score_model.pt
   wget https://github.com/gcorso/DiffDock/releases/download/v1.1/confidence_model.pt
   ```

---

### DiffDock Import Error

**Problem:** `ImportError: No module named 'diffdock'`

**Solution:**

1. **Reinstall DiffDock**
   ```bash
   docker-compose exec backend pip uninstall diffdock -y
   docker-compose exec backend bash -c "cd /tmp && git clone https://github.com/gcorso/DiffDock.git && cd DiffDock && pip install -e ."
   ```

2. **Check dependencies**
   ```bash
   docker-compose exec backend pip list | grep -E "torch|e3nn|spyrmsd|prody|biopython"
   ```

---

## Alternative: CPU-Only Mode

If GPU setup is problematic, DiffDock can run on CPU (slower).

**Update `.env`:**
```env
DIFFDOCK_USE_GPU=false
DIFFDOCK_DEVICE=cpu
```

**Performance Impact:**
- 10-20x slower than GPU
- 30-60 minutes per complex
- Suitable for testing only

**Recommendation:** Use AutoDock Vina for CPU-based docking

---

## Integration with PharmForge

### Using DiffDock in Pipelines

```python
# Example pipeline: Protein-ligand docking with DiffDock
from backend.core.pipeline import Pipeline

pipeline = Pipeline([
    {
        "adapter": "alphafold",
        "params": {"protein_sequence": "MKTAYIAKQ..."}
    },
    {
        "adapter": "diffdock",
        "params": {
            "ligand_smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "use_gpu": True,
            "samples_per_complex": 40
        }
    }
])

result = await pipeline.execute()
```

---

### Comparing DiffDock vs Vina

```python
# Run both docking methods for comparison
from adapters.diffdock.adapter import DiffDockAdapter
from adapters.vina.adapter import VinaAdapter

# DiffDock (GPU-accelerated)
diffdock = DiffDockAdapter()
dd_result = diffdock.execute(
    protein_pdb='1abc',
    ligand_smiles='CCO'
)

# AutoDock Vina (CPU)
vina = VinaAdapter()
vina_result = vina.execute(
    receptor_pdb='1abc.pdb',
    ligand_smiles='CCO',
    center=[10.0, 10.0, 10.0],
    size=[20.0, 20.0, 20.0]
)

# Compare results
print(f"DiffDock best confidence: {dd_result['best_confidence']:.3f}")
print(f"Vina best affinity: {vina_result['affinity']:.3f} kcal/mol")
```

---

## Storage Requirements

| Component | Size | Location |
|-----------|------|----------|
| **DiffDock package** | ~500 MB | Container Python env |
| **Score model** | ~3 GB | `/app/diffdock_models/` |
| **Confidence model** | ~2 GB | `/app/diffdock_models/` |
| **ESM embeddings cache** | Variable | `/app/cache/diffdock/` |
| **Total** | ~5.5 GB+ | - |

---

## Next Steps

### Immediate (Required for GPU acceleration)

1. ✅ Docker GPU configuration complete
2. ⏳ Enable GPU passthrough in Docker Desktop (Windows)
3. ⏳ Install DiffDock package
4. ⏳ Download pre-trained models (~5GB)
5. ⏳ Test GPU accessibility
6. ⏳ Run sample docking

### Optional (Performance optimization)

1. Configure model caching
2. Optimize batch sizes for RTX 5080
3. Benchmark against Vina on test set
4. Implement ensemble docking (DiffDock + Vina)

---

## References

- **DiffDock Paper:** https://arxiv.org/abs/2210.01776
- **GitHub Repository:** https://github.com/gcorso/DiffDock
- **Pre-trained Models:** https://github.com/gcorso/DiffDock/releases
- **NVIDIA Docker:** https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html

---

**Document Version:** 1.0
**Last Updated:** October 26, 2025
**GPU Configuration Status:** ✅ Configured, ⏳ Testing required
