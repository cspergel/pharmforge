# Vina Adapter - Pipeline Integration Guide

Quick reference for integrating the Vina adapter into PharmForge pipelines.

## Quick Start

### 1. Import and Register

The adapter is automatically registered in `backend/core/adapter_registry.py`:

```python
from adapters.vina.adapter import VinaAdapter

# In register_all_adapters():
registry.register(VinaAdapter())
```

### 2. Access via Registry

```python
from backend.core.adapters.protocol import registry

# Get the adapter
vina = registry.get("vina_docking")

# Execute docking
result = await vina.execute("CCO", receptor_path="/path/to/receptor.pdbqt")
```

## Pipeline Integration Patterns

### Pattern 1: Property → Docking Pipeline

Filter molecules by properties, then dock the best candidates:

```python
from backend.core.adapters.protocol import registry

async def property_docking_pipeline(smiles_list, receptor_path):
    rdkit = registry.get("rdkit_local")
    admet = registry.get("admet_ai")
    vina = registry.get("vina_docking")

    results = []

    for smiles in smiles_list:
        # Step 1: Calculate properties
        props = await rdkit.execute(smiles)
        if not props.success:
            continue

        # Step 2: Filter by Lipinski's Rule of 5
        if props.data['lipinski_violations'] > 1:
            continue  # Skip compounds with >1 violation

        # Step 3: ADMET prediction
        admet_result = await admet.execute(smiles)
        if not admet_result.success:
            continue

        # Step 4: Filter by ADMET
        if admet_result.data['properties'].get('hERG', 1.0) > 0.5:
            continue  # Skip if likely hERG liability

        # Step 5: Dock filtered compounds
        docking = await vina.execute(
            smiles,
            receptor_path=receptor_path,
            center_x=10.0,
            center_y=20.0,
            center_z=15.0
        )

        if docking.success:
            results.append({
                'smiles': smiles,
                'properties': props.data,
                'admet': admet_result.data,
                'docking': docking.data
            })

    # Rank by binding affinity
    ranked = sorted(results, key=lambda x: x['docking']['binding_affinity'])
    return ranked
```

### Pattern 2: Multi-Target Docking

Dock against multiple receptors and aggregate scores:

```python
async def multi_target_docking(smiles, receptors):
    vina = registry.get("vina_docking")

    docking_results = []

    for receptor_info in receptors:
        result = await vina.execute(
            smiles,
            receptor_path=receptor_info['path'],
            center_x=receptor_info['center_x'],
            center_y=receptor_info['center_y'],
            center_z=receptor_info['center_z']
        )

        if result.success:
            docking_results.append({
                'target': receptor_info['name'],
                'affinity': result.data['binding_affinity'],
                'score': result.data['binding_score']
            })

    # Calculate consensus score
    if docking_results:
        avg_score = sum(r['score'] for r in docking_results) / len(docking_results)
        best_affinity = min(r['affinity'] for r in docking_results)

        return {
            'smiles': smiles,
            'multi_target_results': docking_results,
            'consensus_score': avg_score,
            'best_affinity': best_affinity
        }

    return None
```

### Pattern 3: Docking-Guided Optimization

Use docking scores to guide molecular optimization:

```python
async def docking_guided_optimization(initial_smiles, receptor_path, generations=5):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    vina = registry.get("vina_docking")

    current_smiles = initial_smiles
    history = []

    for gen in range(generations):
        # Dock current molecule
        result = await vina.execute(current_smiles, receptor_path=receptor_path)

        if not result.success:
            break

        history.append({
            'generation': gen,
            'smiles': current_smiles,
            'affinity': result.data['binding_affinity'],
            'score': result.data['binding_score']
        })

        # Generate variants (simplified - use real optimization)
        mol = Chem.MolFromSmiles(current_smiles)
        # ... mutation/crossover logic ...

        # Select best variant for next generation
        # current_smiles = best_variant

    return history
```

## Multi-Objective Ranking with Docking

Combine docking scores with other objectives:

```python
from backend.core.scoring_utils import vina_affinity_to01

async def multi_objective_ranking(smiles_list, receptor_path):
    rdkit = registry.get("rdkit_local")
    admet = registry.get("admet_ai")
    vina = registry.get("vina_docking")

    compounds = []

    for smiles in smiles_list:
        # Get all predictions
        props = await rdkit.execute(smiles)
        admet_pred = await admet.execute(smiles)
        docking = await vina.execute(smiles, receptor_path=receptor_path)

        if all([props.success, admet_pred.success, docking.success]):
            # Normalize all objectives to 0-1 (higher=better)
            objectives = {
                'binding_score': docking.data['binding_score'],  # Already normalized
                'druglikeness': 1.0 - (props.data['lipinski_violations'] / 4.0),
                'solubility': admet_pred.data['properties'].get('Solubility_AqSolDB', 0.5),
                'safety': 1.0 - admet_pred.data['properties'].get('hERG', 0.5)
            }

            # Calculate weighted score
            weights = {
                'binding_score': 0.4,    # 40% weight on binding
                'druglikeness': 0.2,     # 20% weight on drug-likeness
                'solubility': 0.2,       # 20% weight on solubility
                'safety': 0.2            # 20% weight on safety
            }

            overall_score = sum(
                objectives[key] * weights[key]
                for key in objectives
            )

            compounds.append({
                'smiles': smiles,
                'objectives': objectives,
                'overall_score': overall_score,
                'binding_affinity': docking.data['binding_affinity']  # For reference
            })

    # Rank by overall score
    ranked = sorted(compounds, key=lambda x: x['overall_score'], reverse=True)
    return ranked
```

## Configuration Management

### Environment-Based Configuration

```python
import os

def get_vina_config(target_name: str):
    """Load Vina configuration from environment or config file"""
    return {
        'receptor_path': os.getenv(f'VINA_RECEPTOR_{target_name.upper()}'),
        'center_x': float(os.getenv(f'VINA_CENTER_X_{target_name.upper()}', 0)),
        'center_y': float(os.getenv(f'VINA_CENTER_Y_{target_name.upper()}', 0)),
        'center_z': float(os.getenv(f'VINA_CENTER_Z_{target_name.upper()}', 0)),
        'size_x': int(os.getenv('VINA_BOX_SIZE_X', 25)),
        'size_y': int(os.getenv('VINA_BOX_SIZE_Y', 25)),
        'size_z': int(os.getenv('VINA_BOX_SIZE_Z', 25)),
        'exhaustiveness': int(os.getenv('VINA_EXHAUSTIVENESS', 8))
    }

# Usage
config = get_vina_config('CDK2')
vina = VinaAdapter(config=config)
```

### YAML Configuration

```yaml
# config/docking_targets.yaml
targets:
  CDK2:
    receptor_path: /data/receptors/cdk2_prepared.pdbqt
    center: [10.5, 20.3, 15.7]
    size: [25, 25, 25]
    exhaustiveness: 8

  EGFR:
    receptor_path: /data/receptors/egfr_prepared.pdbqt
    center: [12.1, 18.9, 22.4]
    size: [30, 30, 30]
    exhaustiveness: 16
```

```python
import yaml

def load_docking_config(config_path: str, target: str):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    target_config = config['targets'][target]
    return {
        'receptor_path': target_config['receptor_path'],
        'center_x': target_config['center'][0],
        'center_y': target_config['center'][1],
        'center_z': target_config['center'][2],
        'size_x': target_config['size'][0],
        'size_y': target_config['size'][1],
        'size_z': target_config['size'][2],
        'exhaustiveness': target_config['exhaustiveness']
    }
```

## Caching for Performance

The Vina adapter supports caching via the protocol's `generate_cache_key()` method:

```python
# Cache key is automatically generated based on:
# - SMILES (canonicalized)
# - Receptor path
# - Docking box parameters
# - Vina settings

cache_key = vina.generate_cache_key(
    smiles,
    receptor_path=receptor_path,
    center_x=10.0,
    center_y=20.0,
    center_z=15.0
)

# Check cache before running
cached_result = redis_client.get(cache_key)
if cached_result:
    return AdapterResult(**json.loads(cached_result))

# Otherwise run docking
result = await vina.execute(smiles, receptor_path=receptor_path)

# Cache the result
redis_client.setex(
    cache_key,
    3600,  # 1 hour TTL
    json.dumps(result.__dict__)
)
```

## Error Handling Best Practices

```python
async def robust_docking(smiles, receptor_path, max_retries=3):
    vina = registry.get("vina_docking")

    for attempt in range(max_retries):
        try:
            result = await vina.execute(
                smiles,
                receptor_path=receptor_path,
                center_x=10.0,
                center_y=20.0,
                center_z=15.0
            )

            if result.success:
                return result
            else:
                logger.warning(f"Docking failed (attempt {attempt+1}): {result.error}")

        except Exception as e:
            logger.error(f"Docking error (attempt {attempt+1}): {e}")

        # Wait before retry
        await asyncio.sleep(2 ** attempt)  # Exponential backoff

    return AdapterResult(
        success=False,
        data=None,
        error=f"Docking failed after {max_retries} attempts"
    )
```

## FastAPI Endpoint Integration

```python
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

router = APIRouter()

class DockingRequest(BaseModel):
    smiles: str
    receptor_path: str
    center_x: float
    center_y: float
    center_z: float

@router.post("/dock")
async def dock_molecule(request: DockingRequest):
    """Dock a molecule against a receptor"""
    vina = registry.get("vina_docking")

    result = await vina.execute(
        request.smiles,
        receptor_path=request.receptor_path,
        center_x=request.center_x,
        center_y=request.center_y,
        center_z=request.center_z
    )

    if not result.success:
        raise HTTPException(status_code=500, detail=result.error)

    return {
        "smiles": request.smiles,
        "binding_affinity": result.data['binding_affinity'],
        "binding_score": result.data['binding_score'],
        "num_poses": result.data['num_poses'],
        "poses": result.data['all_poses']
    }
```

## Performance Tips

### 1. Batch Processing with Concurrency Limits

```python
import asyncio

async def batch_dock_with_limit(compounds, receptor_path, max_concurrent=4):
    vina = registry.get("vina_docking")
    semaphore = asyncio.Semaphore(max_concurrent)

    async def dock_one(smiles):
        async with semaphore:
            return await vina.execute(smiles, receptor_path=receptor_path)

    tasks = [dock_one(smi) for smi in compounds]
    results = await asyncio.gather(*tasks, return_exceptions=True)

    return [r for r in results if isinstance(r, AdapterResult) and r.success]
```

### 2. Fast Pre-filtering

```python
async def filter_then_dock(smiles_list, receptor_path):
    rdkit = registry.get("rdkit_local")
    vina = registry.get("vina_docking")

    # Fast filter (< 1ms per molecule)
    filtered = []
    for smiles in smiles_list:
        props = await rdkit.execute(smiles)
        if props.success and props.data['molecular_weight'] <= 500:
            filtered.append(smiles)

    # Slow docking only on filtered set (minutes per molecule)
    docking_results = []
    for smiles in filtered:
        result = await vina.execute(smiles, receptor_path=receptor_path)
        if result.success:
            docking_results.append(result)

    return docking_results
```

### 3. Early Termination

```python
async def dock_until_hit(smiles_list, receptor_path, affinity_threshold=-8.0):
    vina = registry.get("vina_docking")

    for smiles in smiles_list:
        result = await vina.execute(smiles, receptor_path=receptor_path)

        if result.success:
            affinity = result.data['binding_affinity']
            if affinity <= affinity_threshold:
                return result  # Found a hit, stop searching

    return None  # No hits found
```

## Testing Integration

```python
import pytest

@pytest.mark.asyncio
async def test_vina_in_pipeline():
    """Test Vina adapter in a full pipeline"""
    from backend.core.adapters.protocol import registry

    # Setup
    test_receptor = "/path/to/test/receptor.pdbqt"
    test_smiles = "CCO"

    # Get adapters
    rdkit = registry.get("rdkit_local")
    vina = registry.get("vina_docking")

    # Run pipeline
    props = await rdkit.execute(test_smiles)
    assert props.success

    docking = await vina.execute(
        test_smiles,
        receptor_path=test_receptor,
        center_x=0.0,
        center_y=0.0,
        center_z=0.0
    )

    # Verify
    if docking.success:  # May fail if Vina not installed
        assert 'binding_affinity' in docking.data
        assert 'binding_score' in docking.data
        assert 0 <= docking.data['binding_score'] <= 1
```

## Common Pitfalls

### ❌ Don't: Use relative paths

```python
# BAD
vina.execute(smiles, receptor_path="./receptor.pdbqt")
```

### ✅ Do: Use absolute paths

```python
# GOOD
import os
receptor_path = os.path.abspath("receptor.pdbqt")
vina.execute(smiles, receptor_path=receptor_path)
```

### ❌ Don't: Forget to check success

```python
# BAD
result = await vina.execute(smiles, receptor_path=path)
affinity = result.data['binding_affinity']  # May crash if failed
```

### ✅ Do: Always check success flag

```python
# GOOD
result = await vina.execute(smiles, receptor_path=path)
if result.success:
    affinity = result.data['binding_affinity']
else:
    logger.error(f"Docking failed: {result.error}")
```

### ❌ Don't: Mix up score directions

```python
# BAD - binding_affinity is negative (lower=better)
ranked = sorted(results, key=lambda x: x['binding_affinity'], reverse=True)
```

### ✅ Do: Use normalized score or sort correctly

```python
# GOOD - Use normalized score (higher=better)
ranked = sorted(results, key=lambda x: x['binding_score'], reverse=True)

# OR sort affinity correctly (lower=better)
ranked = sorted(results, key=lambda x: x['binding_affinity'])
```

## Next Steps

- See `adapters/vina/README.md` for detailed documentation
- See `adapters/vina/example_usage.py` for complete examples
- See `backend/tests/test_vina_adapter.py` for test cases
- See PharmForge pipeline documentation for advanced workflows
