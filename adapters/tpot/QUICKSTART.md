# TPOT Adapter - Quick Start Guide

## Installation

```bash
pip install tpot scikit-learn numpy
```

## Basic Usage (30 seconds)

```python
from adapters.tpot import TPOTAdapter
import numpy as np

# Your data
X_train = np.array([[...], [...]])  # Features
y_train = np.array([0, 1, 0, 1])    # Labels

# Create and run
adapter = TPOTAdapter()
result = await adapter.execute({
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification",
    "generations": 5,
    "population_size": 20
})

# Get results
if result.success:
    print(f"Score: {result.data['best_score']}")
    print(f"Pipeline: {result.data['best_pipeline']}")
```

## Common Use Cases

### 1. Quick Classification Model
```python
result = await adapter.execute({
    "X_train": X,
    "y_train": y,
    "task": "classification"
})
```

### 2. Regression Model
```python
result = await adapter.execute({
    "X_train": X,
    "y_train": y,
    "task": "regression",
    "scoring": "r2"
})
```

### 3. Fast Mode (TPOT Light)
```python
result = await adapter.execute({
    "X_train": X,
    "y_train": y,
    "task": "classification",
    "config_dict": "TPOT light",
    "generations": 3,
    "population_size": 10
})
```

### 4. Time-Limited Optimization
```python
result = await adapter.execute({
    "X_train": X,
    "y_train": y,
    "task": "classification",
    "max_time_mins": 10  # Stop after 10 minutes
})
```

### 5. With Test Set
```python
result = await adapter.execute({
    "X_train": X_train,
    "y_train": y_train,
    "X_test": X_test,
    "y_test": y_test,
    "task": "classification"
})
print(f"Test score: {result.data['test_score']}")
```

## With Molecular Descriptors

```python
from adapters.mordred import MordredAdapter
from adapters.tpot import TPOTAdapter

# Generate descriptors
mordred = MordredAdapter()
desc_result = await mordred.execute({
    "smiles": ["CCO", "c1ccccc1", ...],
    "descriptors": ["all"]
})

# Optimize model
tpot = TPOTAdapter()
result = await tpot.execute({
    "X_train": desc_result.data["descriptors"],
    "y_train": [1, 0, 1, ...],
    "task": "classification"
})
```

## Export Pipeline Code

```python
result = await adapter.execute({...})
if result.success:
    # Get executable Python code
    code = result.data['pipeline_code']

    # Save to file
    with open('optimized_pipeline.py', 'w') as f:
        f.write(code)
```

## Scoring Metrics

### Classification
- `accuracy` - Overall accuracy
- `roc_auc` - ROC AUC (default for binary)
- `f1` - F1 score
- `f1_weighted` - Weighted F1
- `balanced_accuracy` - Balanced accuracy

### Regression
- `r2` - R² score
- `neg_mean_squared_error` - Negative MSE (default)
- `neg_mean_absolute_error` - Negative MAE
- `neg_root_mean_squared_error` - Negative RMSE

## Performance Tips

### Fast (Minutes)
```python
{
    "config_dict": "TPOT light",
    "generations": 3,
    "population_size": 10,
    "cv": 3
}
```

### Balanced (10-30 minutes)
```python
{
    "generations": 5,
    "population_size": 20,
    "cv": 5,
    "max_time_mins": 30
}
```

### Thorough (Hours)
```python
{
    "generations": 20,
    "population_size": 100,
    "cv": 10,
    "max_time_mins": 120
}
```

## Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `X_train` | array | required | Feature matrix |
| `y_train` | array | required | Target values |
| `task` | str | "classification" | "classification" or "regression" |
| `generations` | int | 5 | Number of generations |
| `population_size` | int | 20 | Population size |
| `scoring` | str | auto | Scoring metric |
| `cv` | int | 5 | Cross-validation folds |
| `max_time_mins` | int | None | Maximum runtime |
| `random_state` | int | None | Random seed |
| `n_jobs` | int | -1 | Parallel jobs |
| `verbosity` | int | 2 | 0=silent, 3=debug |

## Output Fields

```python
result.data = {
    "best_pipeline": str,           # Pipeline description
    "best_score": float,            # CV score
    "test_score": float,            # Test score (if provided)
    "pipeline_code": str,           # Exportable code
    "fitted_pipeline": object,      # Fitted sklearn pipeline
    "pipeline_steps": List[Dict],   # Step details
    "optimization_history": [...],  # Generation stats
    "metadata": {...}               # Optimization metadata
}
```

## Troubleshooting

### "TPOT is not installed"
```bash
pip install tpot
```

### "Out of memory"
- Reduce `population_size` to 10-20
- Use `config_dict="TPOT light"`
- Sample your dataset

### "Taking too long"
- Set `max_time_mins=10`
- Reduce `generations` to 3-5
- Use `config_dict="TPOT light"`

### "Poor performance"
- Increase `generations` to 10-20
- Increase `population_size` to 50-100
- Try different `scoring` metric
- Check data quality and preprocessing

## Next Steps

1. **Read full documentation**: `README.md`
2. **See integrations**: `INTEGRATION.md`
3. **Run examples**: `python example_usage.py`
4. **Combine with Optuna**: For meta-optimization

## Support

- Documentation: `README.md`
- Examples: `example_usage.py`
- Integration guide: `INTEGRATION.md`
- TPOT docs: http://epistasislab.github.io/tpot/
