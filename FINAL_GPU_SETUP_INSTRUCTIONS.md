# Final GPU Setup Instructions for RTX 5080

## Current Status

**Second rebuild in progress** - Installing PyTorch 2.7.0 dev (nightly) with sm_120 support for your RTX 5080.

### Why Two Builds?

1. **First build**: Installed PyTorch nightly, then `admet-ai` reinstalled stable PyTorch over it
2. **Second build** (current): Install requirements first, THEN force-reinstall PyTorch nightly last

This ensures PyTorch nightly stays installed.

---

## When Build Completes

### Step 1: Restart Services
```bash
cd claude-code-agents-wizard-v2
docker-compose up -d
```

### Step 2: Verify GPU Support
```bash
docker exec pharmforge-backend python -c "
import torch
print('PyTorch:', torch.__version__)
print('CUDA available:', torch.cuda.is_available())
print('GPU:', torch.cuda.get_device_name(0))
props = torch.cuda.get_device_properties(0)
print(f'Compute capability: {props.major}.{props.minor}')
print('Supported archs:', torch.cuda.get_arch_list())
if 'sm_120' in torch.cuda.get_arch_list():
    print('✅ SUCCESS! sm_120 support confirmed!')
"
```

**Expected output:**
```
PyTorch: 2.7.0.dev20250226+cu124  ← Should be 2.7.x dev
CUDA available: True
GPU: NVIDIA GeForce RTX 5080
Compute capability: 12.0
Supported archs: [..., 'sm_120', ...]  ← KEY: Must include sm_120!
✅ SUCCESS! sm_120 support confirmed!
```

### Step 3: Test ADMET AI

1. Open http://localhost:8501
2. Go to **"🧪 Compound Testing"**
3. Enter SMILES: `CCO`
4. Select **"Admet Ai"** checkbox
5. Click **"Run Adapters"**

**Expected result:**
- ✅ **Success** status
- ADMET predictions displayed (toxicity, solubility, etc.)
- Execution time: 10-20 seconds first run (model loading), then 1-3 seconds cached

---

## If Build Fails or Takes Too Long

### Alternative: Quick CPU-Only Fix

If you want to test the system NOW without GPU, you can use CPU mode temporarily:

```bash
cd claude-code-agents-wizard-v2
docker exec pharmforge-backend pip install torch --index-url https://download.pytorch.org/whl/cpu
docker-compose restart backend
```

**Trade-offs:**
- ✅ Works immediately
- ✅ ADMET AI will run (on CPU)
- ❌ Much slower (20-30 seconds per compound vs 1-3 seconds)
- ❌ Can't use large batches

---

## Checking Build Progress

### Is build still running?
```bash
docker ps -a | grep building
```

### View live build output:
```bash
cd claude-code-agents-wizard-v2
docker-compose logs --follow backend
```

### Estimated time remaining:
- Most dependencies already cached: ✅
- PyTorch nightly download: ~2-3 minutes
- PyTorch nightly install: ~1-2 minutes
- **Total**: ~3-5 minutes from when it started

---

## What We've Accomplished So Far

### ✅ Frontend Fixed:
- Compound testing works perfectly
- 39 adapters accessible
- Professional UI with status indicators
- Marketplace displays correctly
- All pages functional

### ✅ GPU Detected:
- RTX 5080 visible to Docker
- CUDA 12.9 driver working
- nvidia-smi shows GPU in container

### 🔄 GPU Support (In Progress):
- PyTorch 2.7.0 dev being installed
- Will enable sm_120 (Blackwell) support
- ADMET AI and ML adapters will work with GPU acceleration

---

## Summary

**What works NOW:**
- ✅ 30+ non-GPU adapters (databases, literature, proteins, etc.)
- ✅ Compound testing interface
- ✅ Full frontend functionality

**What will work AFTER build:**
- 🔄 ADMET AI (GPU accelerated)
- 🔄 All ML/deep learning adapters
- 🔄 Fast batch processing with your RTX 5080

**Your hardware:**
- 💪 RTX 5080 with 16GB VRAM
- 💪 CUDA 12.9 driver
- 💪 Blackwell architecture (cutting edge!)

---

## Next Steps

1. **Wait** for build to complete (~3-5 min from start)
2. **Restart**: `docker-compose up -d`
3. **Verify**: Run the GPU verification command above
4. **Test**: Try ADMET AI in compound testing

I'll monitor the build and let you know when it completes!

---

## If You're Impatient

You can test the non-GPU adapters RIGHT NOW while build runs:

1. Go to http://localhost:8501
2. Click "🧪 Compound Testing"
3. Enter: `CCO`
4. Select these (all CPU-based, work perfectly):
   - ✅ Rdkit Local
   - ✅ PubChem
   - ✅ ChEMBL
5. Click "Run Adapters"

These will work instantly and show you the system is functional!
