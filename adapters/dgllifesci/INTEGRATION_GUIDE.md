# DGL-LifeSci Adapter Integration Guide

## Quick Integration

### 1. Import and Register

```python
from adapters.dgllifesci import DGLLifeSciAdapter
from backend.core.adapters.protocol import registry

# Create and register adapter
adapter = DGLLifeSciAdapter()
registry.register(adapter)
```

### 2. Use via Registry

```python
# Get adapter from registry
adapter = registry.get("dgllifesci")

# Execute
result = await adapter("CCO", featurizer_type="canonical")
```

### 3. Direct Usage

```python
from adapters.dgllifesci import DGLLifeSciAdapter

# Create adapter instance
adapter = DGLLifeSciAdapter()

# Use directly
result = await adapter.execute("CCO", featurizer_type="canonical")
```

## Workflow Integration Examples

### Example 1: Property Prediction Pipeline

```python
async def predict_molecular_properties(smiles_list):
    """
    Predict properties for a list of molecules using GNN
    """
    from adapters.dgllifesci import DGLLifeSciAdapter

    # Initialize adapter
    adapter = DGLLifeSciAdapter()

    # Load pre-trained model
    adapter.load_model("models/gcn_solubility.pt", model_type="gcn")

    # Process molecules
    result = await adapter.execute(
        smiles_list,
        featurizer_type="canonical",
        include_predictions=True
    )

    if result.success:
        return result.data['predictions']
    else:
        raise Exception(f"Prediction failed: {result.error}")
```

### Example 2: Feature Extraction for ML

```python
async def extract_graph_features(smiles_list):
    """
    Extract graph features for downstream machine learning
    """
    from adapters.dgllifesci import DGLLifeSciAdapter

    adapter = DGLLifeSciAdapter()

    # Generate graphs with features
    result = await adapter.execute(
        smiles_list,
        featurizer_type="canonical",
        return_graphs=False  # Get JSON-serializable features
    )

    if result.success:
        # Extract features for ML model
        features = []
        for graph_feat in result.data['graph_features']:
            feature_vector = {
                'num_nodes': graph_feat['num_nodes'],
                'num_edges': graph_feat['num_edges'],
                'node_feature_mean': graph_feat['node_features']['mean'],
                'node_feature_std': graph_feat['node_features']['std'],
                'edge_feature_mean': graph_feat['edge_features']['mean'],
                'edge_feature_std': graph_feat['edge_features']['std'],
            }
            features.append(feature_vector)

        return features
    else:
        raise Exception(f"Feature extraction failed: {result.error}")
```

### Example 3: Multi-Adapter Workflow

```python
async def comprehensive_molecular_analysis(smiles):
    """
    Combine multiple adapters for comprehensive analysis
    """
    from adapters.dgllifesci import DGLLifeSciAdapter
    from adapters.rdkit_local.adapter import RDKitAdapter
    from adapters.deepchem.adapter import DeepChemAdapter

    # Initialize adapters
    dgl_adapter = DGLLifeSciAdapter()
    rdkit_adapter = RDKitAdapter()
    deepchem_adapter = DeepChemAdapter()

    # Gather results from all adapters
    results = {}

    # RDKit properties
    rdkit_result = await rdkit_adapter.execute(smiles)
    if rdkit_result.success:
        results['rdkit_properties'] = rdkit_result.data

    # DGL-LifeSci graph features
    dgl_result = await dgl_adapter.execute(smiles, featurizer_type="canonical")
    if dgl_result.success:
        results['dgl_graph_features'] = dgl_result.data['graph_features']

    # DeepChem features
    deepchem_result = await deepchem_adapter.execute(smiles, featurizer_type="morgan")
    if deepchem_result.success:
        results['deepchem_fingerprint'] = deepchem_result.data['features']

    return results
```

### Example 4: Batch Screening Pipeline

```python
async def screen_compound_library(smiles_library, model_path):
    """
    Screen a large library of compounds using GNN model
    """
    from adapters.dgllifesci import DGLLifeSciAdapter
    import asyncio

    adapter = DGLLifeSciAdapter()
    adapter.load_model(model_path, model_type="gcn")

    # Process in batches
    batch_size = 100
    results = []

    for i in range(0, len(smiles_library), batch_size):
        batch = smiles_library[i:i+batch_size]

        result = await adapter.execute(
            batch,
            featurizer_type="canonical",
            include_predictions=True
        )

        if result.success:
            results.append(result.data)
        else:
            print(f"Batch {i//batch_size} failed: {result.error}")

        # Rate limiting
        await asyncio.sleep(0.1)

    return results
```

### Example 5: Model Comparison

```python
async def compare_featurizers(smiles):
    """
    Compare different featurizers for the same molecule
    """
    from adapters.dgllifesci import DGLLifeSciAdapter

    adapter = DGLLifeSciAdapter()
    featurizers = ["canonical", "attentivefp", "weave", "minimal"]

    results = {}

    for featurizer in featurizers:
        result = await adapter.execute(smiles, featurizer_type=featurizer)

        if result.success:
            results[featurizer] = {
                'num_nodes': result.data['graph_features']['num_nodes'],
                'num_edges': result.data['graph_features']['num_edges'],
                'node_feature_dim': result.data['graph_features']['node_features']['feature_dim'],
                'edge_feature_dim': result.data['graph_features']['edge_features']['feature_dim'],
            }

    return results
```

## API Integration

### REST API Endpoint Example

```python
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from adapters.dgllifesci import DGLLifeSciAdapter

app = FastAPI()
adapter = DGLLifeSciAdapter()

class MoleculeRequest(BaseModel):
    smiles: str
    featurizer_type: str = "canonical"
    include_predictions: bool = False

@app.post("/api/dgllifesci/featurize")
async def featurize_molecule(request: MoleculeRequest):
    """Generate molecular graph features"""
    result = await adapter.execute(
        request.smiles,
        featurizer_type=request.featurizer_type,
        include_predictions=request.include_predictions
    )

    if result.success:
        return {
            "status": "success",
            "data": result.data,
            "metadata": result.metadata
        }
    else:
        raise HTTPException(status_code=400, detail=result.error)

@app.post("/api/dgllifesci/batch-featurize")
async def batch_featurize(smiles_list: list[str], featurizer_type: str = "canonical"):
    """Batch featurize molecules"""
    result = await adapter.execute(smiles_list, featurizer_type=featurizer_type)

    if result.success:
        return {
            "status": "success",
            "n_successful": result.data['n_successful'],
            "n_failed": result.data['n_failed'],
            "data": result.data
        }
    else:
        raise HTTPException(status_code=400, detail=result.error)
```

### GraphQL Integration

```python
import strawberry
from adapters.dgllifesci import DGLLifeSciAdapter

@strawberry.type
class GraphFeatures:
    smiles: str
    num_nodes: int
    num_edges: int
    node_feature_dim: int
    edge_feature_dim: int

@strawberry.type
class Query:
    @strawberry.field
    async def molecular_graph(self, smiles: str, featurizer_type: str = "canonical") -> GraphFeatures:
        adapter = DGLLifeSciAdapter()
        result = await adapter.execute(smiles, featurizer_type=featurizer_type)

        if result.success:
            gf = result.data['graph_features']
            return GraphFeatures(
                smiles=gf['smiles'],
                num_nodes=gf['num_nodes'],
                num_edges=gf['num_edges'],
                node_feature_dim=gf['node_features']['feature_dim'],
                edge_feature_dim=gf['edge_features']['feature_dim']
            )
        else:
            raise Exception(result.error)

schema = strawberry.Schema(query=Query)
```

## Configuration Management

### Environment-based Configuration

```python
import os
from adapters.dgllifesci import DGLLifeSciAdapter

# Configure adapter based on environment
adapter = DGLLifeSciAdapter()

adapter.config.update({
    'device': 'cuda' if os.getenv('USE_GPU', 'false').lower() == 'true' else 'cpu',
    'timeout': int(os.getenv('ADAPTER_TIMEOUT', '60')),
    'default_featurizer': os.getenv('DEFAULT_FEATURIZER', 'canonical'),
    'batch_size': int(os.getenv('BATCH_SIZE', '32'))
})
```

### Dynamic Model Loading

```python
from pathlib import Path
from adapters.dgllifesci import DGLLifeSciAdapter

class DGLLifeSciModelManager:
    def __init__(self, model_dir: str):
        self.model_dir = Path(model_dir)
        self.adapter = DGLLifeSciAdapter()
        self.loaded_models = {}

    def load_model(self, task: str, model_type: str = "gcn"):
        """Load model for specific task"""
        model_path = self.model_dir / f"{task}_{model_type}.pt"

        if model_path.exists():
            success = self.adapter.load_model(str(model_path), model_type=model_type)
            if success:
                self.loaded_models[task] = {
                    'path': str(model_path),
                    'type': model_type
                }
                return True

        return False

    async def predict(self, smiles, task: str):
        """Make predictions for a task"""
        if task not in self.loaded_models:
            if not self.load_model(task):
                raise ValueError(f"Model for task '{task}' not found")

        result = await self.adapter.execute(
            smiles,
            featurizer_type="canonical",
            include_predictions=True
        )

        return result
```

## Error Handling Best Practices

### Robust Error Handling

```python
from adapters.dgllifesci import DGLLifeSciAdapter
import logging

logger = logging.getLogger(__name__)

async def safe_featurize(smiles, featurizer_type="canonical"):
    """Safely featurize molecules with comprehensive error handling"""
    try:
        adapter = DGLLifeSciAdapter()
        result = await adapter.execute(smiles, featurizer_type=featurizer_type)

        if result.success:
            # Check for partial failures in batch processing
            if isinstance(smiles, list) and result.data.get('n_failed', 0) > 0:
                logger.warning(
                    f"Partial batch failure: {result.data['n_failed']} molecules failed"
                )
                logger.warning(f"Failed SMILES: {result.data.get('failed_smiles', [])}")

            return result.data
        else:
            logger.error(f"Featurization failed: {result.error}")
            return None

    except ImportError as e:
        logger.error(f"Required package not installed: {e}")
        logger.info("Install with: pip install dgllife")
        return None

    except Exception as e:
        logger.error(f"Unexpected error during featurization: {e}", exc_info=True)
        return None
```

## Performance Optimization

### Caching Strategy

```python
from adapters.dgllifesci import DGLLifeSciAdapter
from functools import lru_cache

class CachedDGLAdapter:
    def __init__(self):
        self.adapter = DGLLifeSciAdapter()

    @lru_cache(maxsize=1000)
    async def get_features(self, smiles: str, featurizer_type: str):
        """Cache molecular features"""
        result = await self.adapter.execute(smiles, featurizer_type=featurizer_type)
        if result.success:
            return result.data
        return None
```

### Parallel Processing

```python
import asyncio
from adapters.dgllifesci import DGLLifeSciAdapter

async def parallel_featurize(smiles_lists, featurizer_type="canonical"):
    """Process multiple batches in parallel"""
    adapter = DGLLifeSciAdapter()

    tasks = [
        adapter.execute(batch, featurizer_type=featurizer_type)
        for batch in smiles_lists
    ]

    results = await asyncio.gather(*tasks, return_exceptions=True)

    successful = [r for r in results if isinstance(r, AdapterResult) and r.success]
    failed = [r for r in results if not isinstance(r, AdapterResult) or not r.success]

    return successful, failed
```

## Testing Integration

### Unit Test Example

```python
import pytest
from adapters.dgllifesci import DGLLifeSciAdapter

@pytest.mark.asyncio
async def test_basic_featurization():
    """Test basic molecular graph generation"""
    adapter = DGLLifeSciAdapter()
    result = await adapter.execute("CCO", featurizer_type="canonical")

    assert result.success
    assert result.data['n_successful'] == 1
    assert result.data['graph_features']['num_nodes'] > 0

@pytest.mark.asyncio
async def test_batch_processing():
    """Test batch processing"""
    adapter = DGLLifeSciAdapter()
    smiles = ["CCO", "c1ccccc1", "CC(=O)O"]

    result = await adapter.execute(smiles, featurizer_type="canonical")

    assert result.success
    assert result.data['n_molecules'] == 3
    assert len(result.data['graph_features']) == 3

@pytest.mark.asyncio
async def test_invalid_smiles():
    """Test handling of invalid SMILES"""
    adapter = DGLLifeSciAdapter()
    result = await adapter.execute("INVALID", featurizer_type="canonical")

    # Should fail validation or processing
    assert not result.success or result.data['n_failed'] > 0
```

### Integration Test Example

```python
import pytest
from adapters.dgllifesci import DGLLifeSciAdapter
from backend.core.adapters.protocol import registry

@pytest.mark.asyncio
async def test_registry_integration():
    """Test integration with adapter registry"""
    # Register adapter
    adapter = DGLLifeSciAdapter()
    registry.register(adapter)

    # Retrieve from registry
    retrieved = registry.get("dgllifesci")
    assert retrieved is not None
    assert retrieved.name == "dgllifesci"

    # Use via registry
    result = await retrieved("CCO", featurizer_type="canonical")
    assert result.success
```

## Monitoring and Logging

### Structured Logging

```python
import logging
import json
from adapters.dgllifesci import DGLLifeSciAdapter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def monitored_featurize(smiles, featurizer_type="canonical"):
    """Featurize with structured logging"""
    adapter = DGLLifeSciAdapter()

    logger.info(json.dumps({
        'event': 'featurization_start',
        'n_molecules': len(smiles) if isinstance(smiles, list) else 1,
        'featurizer_type': featurizer_type
    }))

    result = await adapter.execute(smiles, featurizer_type=featurizer_type)

    logger.info(json.dumps({
        'event': 'featurization_complete',
        'success': result.success,
        'n_successful': result.data.get('n_successful') if result.success else 0,
        'n_failed': result.data.get('n_failed') if result.success else 0,
        'cache_hit': result.cache_hit
    }))

    return result
```

## Production Checklist

- [ ] Install required dependencies (`pip install dgllife`)
- [ ] Configure adapter settings (device, timeout, batch_size)
- [ ] Load pre-trained models if needed
- [ ] Implement error handling and logging
- [ ] Set up monitoring and metrics
- [ ] Configure caching strategy
- [ ] Test with production data
- [ ] Document API endpoints
- [ ] Set up rate limiting if needed
- [ ] Configure resource limits (memory, GPU)

## Support and Resources

- **Documentation**: See `README.md` for full documentation
- **Examples**: See `example_usage.py` for code examples
- **Quick Start**: See `QUICK_START.md` for quick reference
- **DGL-LifeSci Docs**: https://lifesci.dgl.ai/
- **PharmForge Docs**: See main project documentation
