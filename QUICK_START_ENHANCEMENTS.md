# PharmForge Quick-Start Enhancements Guide

**Date:** October 26, 2025
**Purpose:** Fast-track guide for DiffDock and REINVENT full installation

---

## Current Status Summary

### ✅ Configured and Ready
- **OpenAI API Key:** Added to `.env` (LLM retrosynthesis ready)
- **Google Custom Search Engine ID:** `50db04cab2a8b4961` added to `.env`
- **Docker GPU Configuration:** Completed (needs GPU enabled in Docker Desktop)
- **Patent Search APIs:** Deferred until launch (Lens.org costs money)
- **23/33 adapters:** Production-ready, no additional setup needed

### ⏳ Next: Full DiffDock and REINVENT Installation

This guide walks you through installing DiffDock and REINVENT for maximum performance.

---

## Option A: DiffDock Full Installation (GPU-Accelerated Docking)

### Prerequisites

**1. Enable GPU in Docker Desktop** (Windows) - **REQUIRED**

```bash
# Steps:
1. Open Docker Desktop
2. Go to Settings → Resources → WSL Integration
3. Enable "Enable integration with my default WSL distro"
4. Enable for your specific WSL distribution
5. Click "Apply & Restart"

# Verify GPU access:
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi
```

**Expected Output:** Should show your RTX 5080 GPU details

---

### Installation Steps

**Step 1: Install DiffDock Dependencies** (~5 minutes)

```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"

# Install required packages
docker-compose exec -T backend pip install e3nn==0.5.1
docker-compose exec -T backend pip install spyrmsd==0.6.0
docker-compose exec -T backend pip install prody==2.4.0
docker-compose exec -T backend pip install biopython==1.81
docker-compose exec -T backend pip install fair-esm==2.0.0

# Install DiffDock from GitHub
docker-compose exec backend bash -c "cd /tmp && git clone https://github.com/gcorso/DiffDock.git && cd DiffDock && pip install -e ."
```

---

**Step 2: Download Pre-trained Models** (~15 minutes, 5GB download)

```bash
# Create model directory
docker-compose exec backend mkdir -p /app/diffdock_models

# Download models (check GitHub for latest URLs)
docker-compose exec backend bash -c "
cd /app/diffdock_models
wget https://github.com/gcorso/DiffDock/releases/download/v1.1/score_model.pt
wget https://github.com/gcorso/DiffDock/releases/download/v1.1/confidence_model.pt
"
```

**Note:** Check https://github.com/gcorso/DiffDock/releases for current model download links

---

**Step 3: Verify Installation**

```bash
# Test GPU access
docker-compose exec -T backend python -c "
import torch
print(f'CUDA available: {torch.cuda.is_available()}')
print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else \"No GPU\"}')
"

# Test DiffDock import
docker-compose exec -T backend python -c "
try:
    import diffdock
    print('✅ DiffDock installed successfully')
except ImportError:
    print('❌ DiffDock not found')
"
```

**Expected:**
- CUDA available: True
- GPU: NVIDIA GeForce RTX 5080
- DiffDock: Installed

---

### Expected Performance

| Task | RTX 5080 (GPU) | CPU Only |
|------|----------------|----------|
| Single docking (40 poses) | 2-5 minutes | 30-60 minutes |
| Batch (10 complexes) | 10-20 minutes | 5-10 hours |
| GPU Memory Usage | 4-8 GB | N/A |

**Recommendation:** GPU essential for practical use. CPU mode only for testing.

---

## Option B: REINVENT Full Installation (Advanced Molecular Generation)

### Current Status
- **Fallback mode:** ✅ Working (uses RDKit for basic generation)
- **Full REINVENT:** Not installed (requires GitHub source)

### When to Install Full REINVENT
- **Use fallback mode if:** You need basic molecular generation with property filtering
- **Install full REINVENT if:** You need RL-based optimization, scaffold hopping, advanced generation

---

### Installation Steps (~30 minutes)

**Step 1: Install REINVENT from GitHub**

```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"

# Clone and install
docker-compose exec backend bash -c "
cd /tmp
git clone https://github.com/MolecularAI/REINVENT4.git
cd REINVENT4
pip install -e .
"
```

---

**Step 2: Download Pre-trained Models** (Optional, ~2GB)

```bash
# Create model directory
docker-compose exec backend mkdir -p /app/reinvent_models

# Download ChEMBL prior (most common use case)
docker-compose exec backend bash -c "
cd /app/reinvent_models
# Check REINVENT4 GitHub for current model URLs
# wget <model_url>
"
```

**Note:** REINVENT can work without pre-trained models (random initialization). Models are optional but improve quality.

---

**Step 3: Verify Installation**

```bash
# Test REINVENT import
docker-compose exec -T backend python -c "
try:
    import reinvent
    print('✅ REINVENT installed successfully')
    print(f'   Version: {reinvent.__version__ if hasattr(reinvent, \"__version__\") else \"Unknown\"}')
except ImportError:
    print('❌ REINVENT not installed')
    print('   Fallback mode still available')
"
```

---

### REINVENT Usage Comparison

| Mode | Setup Required | Capabilities | Quality |
|------|----------------|--------------|---------|
| **Fallback (current)** | None | Random generation, property filtering | Basic |
| **Full REINVENT** | GitHub install | RL optimization, scaffold hopping | Advanced |
| **With models** | + Download models | Transfer learning, faster convergence | Best |

**Recommendation:** Keep using fallback mode unless you need RL-based optimization. The quality improvement requires significant computational time.

---

## Quick Installation (Both Tools)

### One-Command Setup

```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"

# Install DiffDock dependencies
docker-compose exec -T backend bash << 'EOF'
pip install e3nn==0.5.1 spyrmsd==0.6.0 prody==2.4.0 biopython==1.81 fair-esm==2.0.0
cd /tmp && git clone https://github.com/gcorso/DiffDock.git && cd DiffDock && pip install -e .
mkdir -p /app/diffdock_models
echo "✅ DiffDock dependencies installed"
EOF

# Install REINVENT (optional)
docker-compose exec -T backend bash << 'EOF'
cd /tmp && git clone https://github.com/MolecularAI/REINVENT4.git && cd REINVENT4 && pip install -e .
mkdir -p /app/reinvent_models
echo "✅ REINVENT installed"
EOF
```

**Time Required:** 10-15 minutes
**Models Download:** Additional 15-30 minutes (5GB for DiffDock)

---

## Troubleshooting

### DiffDock: GPU Not Detected

**Problem:** `torch.cuda.is_available()` returns `False`

**Solutions:**

1. **Enable Docker Desktop GPU**
   - Settings → Resources → WSL Integration
   - Restart Docker Desktop

2. **Check Windows NVIDIA Drivers**
   ```bash
   # In Windows PowerShell
   nvidia-smi
   ```
   Should show RTX 5080

3. **Rebuild Containers**
   ```bash
   docker-compose down
   docker-compose build --no-cache backend
   docker-compose up -d
   ```

---

### REINVENT: Import Error

**Problem:** `ModuleNotFoundError: No module named 'reinvent'`

**Solution:**
```bash
# Reinstall
docker-compose exec backend bash
cd /tmp/REINVENT4
pip install -e .
```

**Alternative:** Use fallback mode (already working)

---

### Out of Disk Space

**Problem:** Insufficient space for models

**Storage Requirements:**
- DiffDock: ~5.5 GB (package + models)
- REINVENT: ~2.2 GB (package + models)
- Total: ~8 GB

**Solution:**
- Free up disk space
- Use external storage for model directory
- Skip optional pre-trained models

---

## Verification Checklist

After installation, verify everything works:

```bash
# 1. GPU Access
docker-compose exec -T backend python -c "import torch; print(f'GPU: {torch.cuda.is_available()}')"
# Expected: True

# 2. DiffDock
docker-compose exec -T backend python -c "import diffdock; print('DiffDock: OK')"
# Expected: DiffDock: OK

# 3. REINVENT
docker-compose exec -T backend python -c "
try:
    import reinvent
    print('REINVENT: Full')
except ImportError:
    print('REINVENT: Fallback mode')
"
# Expected: REINVENT: Full (or Fallback mode is OK too)

# 4. Environment Variables
docker-compose exec -T backend python -c "
import os
print(f'OpenAI API: {\"✅\" if os.getenv(\"OPENAI_API_KEY\") else \"❌\"}')
print(f'Google CSE: {\"✅\" if os.getenv(\"GOOGLE_CSE_ID\") else \"❌\"}')
"
# Expected: Both ✅
```

---

## What's Been Done

### ✅ Completed This Session

1. **Environment Configuration**
   - Updated `.env` with OpenAI API key
   - Added Google Custom Search Engine ID: `50db04cab2a8b4961`
   - Structured .env file with all API key placeholders
   - Marked Lens.org as deferred (cost consideration)

2. **Docker GPU Configuration**
   - Updated `docker-compose.yml` with GPU support
   - Added NVIDIA environment variables
   - Configured GPU for backend and celery-worker

3. **Documentation Created**
   - `API_KEYS_SETUP_GUIDE.md` (500+ lines)
   - `DIFFDOCK_SETUP_GUIDE.md` (600+ lines)
   - `BACKEND_ENHANCEMENTS_SUMMARY.md` (comprehensive overview)
   - This quick-start guide

4. **Package Installations**
   - OpenMM 8.3.1 installed
   - OpenAI package (requires rebuild to persist)

---

## Next Steps (Your Choice)

### Option 1: Install Enhancements Now (30-60 minutes)

Follow the installation steps above for:
1. DiffDock (requires GPU enabled first)
2. REINVENT (optional, fallback mode works)

### Option 2: Test Current Setup (~5 minutes)

```bash
# Quick test of current functionality
docker-compose exec -T backend python -c "
# Test fallback mode adapters
from adapters.reinvent.adapter import REINVENTAdapter
from adapters.molgan.adapter import MolGANAdapter

reinvent = REINVENTAdapter()
print(f'✅ REINVENT (fallback): {reinvent.name}')

molgan = MolGANAdapter()
print(f'✅ MolGAN (fallback): {molgan.name}')
"
```

### Option 3: Continue to Frontend Integration

Proceed with building the Streamlit/React frontend per `FRONTEND_INTEGRATION_ROADMAP.md`

---

## API Key Status

| Service | Status | Purpose | Cost |
|---------|--------|---------|------|
| **OpenAI** | ✅ Configured | LLM retrosynthesis | Pay-as-you-go (~$0.03/1K tokens) |
| **Google CSE** | ✅ Configured | Patent search | Free (with limits) |
| **Google Patents API** | ⚠️ Key needed | Patent API access | Free |
| **Lens.org** | ⏭️ Deferred | Patent search | $99/month (free tier 100 req/day) |
| **Anthropic** | ❌ Not configured | Alternative to OpenAI | Pay-as-you-go |

**Note:** Add Google Patents API key to `.env` when available. Line 48 in `.env` file.

---

## Recommended Priority

1. **High Priority** (Do first):
   - Enable GPU in Docker Desktop
   - Install DiffDock (essential for production docking)

2. **Medium Priority** (If needed):
   - Install REINVENT full version (fallback works for now)
   - Add Google Patents API key to `.env`

3. **Low Priority** (Defer):
   - Lens.org API key (costs money, defer until launch)
   - Anthropic API key (OpenAI already configured)

---

## Summary

**What You Have Now:**
- ✅ 23/33 adapters production-ready
- ✅ OpenAI API for LLM retrosynthesis
- ✅ Google CSE for patent search
- ✅ Docker GPU configuration complete
- ✅ All ML adapters working (fallback modes)
- ✅ Comprehensive documentation

**What You Need:**
- ⏳ Enable GPU in Docker Desktop (15 min)
- ⏳ Install DiffDock (30 min + 5GB download)
- ⏳ Optionally install REINVENT (30 min)

**Total Time to Full Setup:** 1-2 hours

---

**Guide Created:** October 26, 2025
**Your Setup:** RTX 5080, Windows 10, Docker Desktop
**Status:** Ready to proceed with enhancements
**Next:** Enable GPU → Install DiffDock → (Optional) Install REINVENT
