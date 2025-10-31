# DiffDock Adapter Setup Instructions

## Overview

DiffDock is a state-of-the-art diffusion-based molecular docking tool that predicts protein-ligand binding poses without requiring binding site specification (blind docking).

**Status**: NEEDS_SETUP
- GPU: Available on host (RTX 5080), but NOT configured in Docker
- RDKit: ✅ Installed in Docker
- PyTorch: ✅ Installed (CPU-only, CUDA compiled but not accessible)
- DiffDock: ❌ Not installed
- ESM Models: ❌ Not downloaded (~2.5GB)

## Requirements

### Hardware
- **GPU**: NVIDIA GPU with 8GB+ VRAM (highly recommended)
  - RTX 3060 or better
  - Host GPU detected: RTX 5080 (16GB VRAM) ✅
- **RAM**: 16GB+ recommended
- **Storage**: ~5GB for models and dependencies

### Software
- Python 3.9-3.11
- CUDA 11.7+ (for GPU acceleration)
- Docker with NVIDIA Container Runtime (for containerized deployment)

## Installation Steps

### Option 1: Local Installation (Recommended for Development)

#### Step 1: Clone DiffDock Repository

```bash
cd /path/to/your/workspace
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
```

#### Step 2: Create Conda Environment

```bash
# Create environment from DiffDock's environment.yml
conda env create -f environment.yml
conda activate diffdock
```

Or create manually:

```bash
conda create -n diffdock python=3.9
conda activate diffdock

# Install PyTorch with CUDA support
conda install pytorch==1.11.0 pytorch-cuda=11.7 -c pytorch -c nvidia

# Install PyTorch Geometric
conda install pyg -c pyg

# Install other dependencies
pip install -r requirements.txt
```

#### Step 3: Download Model Weights

```bash
# Download pre-trained models (~1.5GB)
# Option A: From the DiffDock releases
wget https://github.com/gcorso/DiffDock/releases/download/v1.0/score_model.pt
wget https://github.com/gcorso/DiffDock/releases/download/v1.0/confidence_model.pt

# Or use the provided download script
python download_models.py
```

Model files should be placed in:
- `workdir/score_model/best_model.pt`
- `workdir/confidence_model/best_model.pt`

#### Step 4: Install ESM (Protein Language Model)

```bash
# Install ESM from Facebook Research
pip install fair-esm

# Or clone and install from source
git clone https://github.com/facebookresearch/esm.git
cd esm
pip install -e .
```

ESM models will be automatically downloaded (~2.5GB) on first use.

#### Step 5: Test Installation

```bash
# Run inference on example
python inference.py \
    --protein_ligand_csv data/example.csv \
    --out_dir results/example \
    --inference_steps 20 \
    --samples_per_complex 40
```

#### Step 6: Configure PharmForge Adapter

Update your PharmForge configuration:

```python
# In your config file or environment
DIFFDOCK_CONFIG = {
    'diffdock_path': '/path/to/DiffDock',
    'model_checkpoint': None,  # Use default models in workdir/
    'inference_steps': 20,
    'samples_per_complex': 40,
    'use_gpu': True,
    'python_env': '/path/to/conda/envs/diffdock/bin/python'
}
```

### Option 2: Docker Installation with GPU Support

#### Step 1: Install NVIDIA Container Toolkit

```bash
# For Ubuntu/Debian
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
    sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

For Windows with WSL2:
1. Ensure you have WSL2 with GPU support
2. Install Docker Desktop for Windows with WSL2 backend
3. Enable GPU support in Docker Desktop settings

#### Step 2: Update docker-compose.yml

Add GPU support to the backend service:

```yaml
backend:
  build:
    context: .
    dockerfile: Dockerfile.backend
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            count: 1
            capabilities: [gpu]
  environment:
    - NVIDIA_VISIBLE_DEVICES=all
```

#### Step 3: Create DiffDock Dockerfile

Create `adapters/diffdock/Dockerfile`:

```dockerfile
FROM nvidia/cuda:11.7.1-cudnn8-runtime-ubuntu22.04

# Install Python and dependencies
RUN apt-get update && apt-get install -y \
    python3.9 python3-pip git wget \
    && rm -rf /var/lib/apt/lists/*

# Clone DiffDock
WORKDIR /app
RUN git clone https://github.com/gcorso/DiffDock.git

# Install Python dependencies
WORKDIR /app/DiffDock
RUN pip3 install torch==1.11.0+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
RUN pip3 install torch-geometric pyg-lib torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-1.11.0+cu117.html
RUN pip3 install -r requirements.txt

# Download models
RUN python3 download_models.py

# Install ESM
RUN pip3 install fair-esm

WORKDIR /app
```

#### Step 4: Build and Run

```bash
docker build -t pharmforge-diffdock -f adapters/diffdock/Dockerfile .
docker run --gpus all pharmforge-diffdock nvidia-smi
```

## Model Files

### Score Model (~1.2GB)
- Location: `workdir/score_model/best_model.pt`
- Purpose: Predicts binding poses through diffusion process
- Parameters: ~20M

### Confidence Model (~300MB)
- Location: `workdir/confidence_model/best_model.pt`
- Purpose: Scores predicted poses for ranking
- Parameters: ~5M

### ESM Protein Model (~2.5GB)
- Model: esm2_t33_650M_UR50D
- Auto-downloaded to: `~/.cache/torch/hub/checkpoints/`
- Purpose: Protein sequence embeddings

**Total Storage**: ~5GB

## Performance Estimates

### With GPU (RTX 3060 or better)
- Single docking run (40 poses): **2-5 minutes**
- Memory usage: **4-6GB VRAM**
- Throughput: ~10-20 complexes/hour

### CPU-only (Not Recommended)
- Single docking run (40 poses): **30-60 minutes**
- Memory usage: **8-12GB RAM**
- Throughput: ~1-2 complexes/hour

### Inference Parameters
- `inference_steps`: 20 (default, can reduce to 10 for faster but less accurate)
- `samples_per_complex`: 40 (number of poses, can reduce to 10 for testing)
- `batch_size`: 10 (GPU batch size)

## Testing the Adapter

### Run Unit Tests

```bash
# From PharmForge root
cd claude-code-agents-wizard-v2
python -m pytest backend/tests/test_diffdock_adapter.py -v
```

Test results:
- ✅ **13 passed**: Basic functionality works
- ⚠️ **10 skipped**: RDKit not installed in test environment
- ❌ **1 failed**: SDF parsing (needs RDKit)

### Run with Sample Data

```python
import asyncio
from adapters.diffdock.adapter import DiffDockAdapter

async def test_docking():
    config = {
        'diffdock_path': '/path/to/DiffDock',
        'use_gpu': True,
        'inference_steps': 20,
        'samples_per_complex': 10  # Fewer for testing
    }

    adapter = DiffDockAdapter(config=config)

    # Example: Aspirin docked to protein
    result = await adapter.execute({
        'smiles': 'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'protein_path': '/path/to/protein.pdb'
    })

    if result.success:
        print(f"Confidence: {result.data['confidence']:.3f}")
        print(f"Estimated Affinity: {result.data['estimated_affinity']:.2f} kcal/mol")
        print(f"Number of Poses: {result.data['num_poses']}")
    else:
        print(f"Error: {result.error}")

asyncio.run(test_docking())
```

## Current Docker Container Status

### Host System
- GPU: NVIDIA RTX 5080 (16GB VRAM) ✅
- CUDA: 12.9 ✅
- Driver: 576.88 ✅

### Docker Container (pharmforge-backend)
- PyTorch: 2.5.0+cu124 (CUDA compiled) ✅
- CUDA Available: ❌ False
- RDKit: 2023.09.6 ✅
- GPU Access: ❌ Not configured
- nvidia-smi: ❌ Not available in container

**Issue**: Docker container needs GPU runtime configuration

## Setup Time Estimates

### Quick Setup (CPU-only for testing)
- **Time**: 30-45 minutes
- **Steps**: Clone repo, install dependencies, download models
- **Usable for**: Testing, development, small-scale docking

### Full GPU Setup (Recommended)
- **Time**: 1-2 hours
- **Steps**:
  1. Clone repo (5 min)
  2. Setup conda environment (10 min)
  3. Download models (10-20 min depending on connection)
  4. Configure Docker GPU support (20-30 min)
  5. Testing and validation (15-30 min)
- **Usable for**: Production, large-scale docking, real-time applications

### Docker GPU Integration
- **Time**: 2-3 hours
- **Steps**: Install NVIDIA Container Toolkit, rebuild images, configure compose
- **Complexity**: Medium (requires Docker and NVIDIA driver knowledge)

## Common Issues

### 1. CUDA Out of Memory
**Solution**: Reduce `batch_size` or `samples_per_complex`

### 2. ESM Model Download Fails
**Solution**: Manually download from Hugging Face:
```bash
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt
```

### 3. PyTorch Geometric Version Conflicts
**Solution**: Install compatible versions:
```bash
pip install torch-geometric -f https://data.pyg.org/whl/torch-${TORCH}+${CUDA}.html
```

### 4. Docker GPU Not Accessible
**Solution**:
- Verify: `nvidia-docker version`
- Check runtime: `docker run --gpus all nvidia/cuda:11.7.1-base nvidia-smi`
- Update docker-compose.yml with GPU configuration

## References

- **DiffDock Paper**: [arXiv:2210.01776](https://arxiv.org/abs/2210.01776)
- **GitHub Repository**: https://github.com/gcorso/DiffDock
- **ESM (Protein LM)**: https://github.com/facebookresearch/esm
- **PyTorch Geometric**: https://pytorch-geometric.readthedocs.io/

## Next Steps

1. **Install DiffDock** following Option 1 or Option 2 above
2. **Download model weights** (~1.5GB)
3. **Install ESM** for protein embeddings
4. **Configure GPU access** in Docker (if using containers)
5. **Update adapter config** with correct paths
6. **Run tests** to verify installation
7. **Test with sample protein-ligand pair**

Estimated total time: **2-3 hours** for full GPU setup
