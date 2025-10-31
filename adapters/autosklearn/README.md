# Auto-sklearn Adapter

Automated Machine Learning adapter using Auto-sklearn with Bayesian optimization and meta-learning.

## Overview

Auto-sklearn is an automated machine learning toolkit that:
- Automatically selects ML algorithms
- Optimizes hyperparameters using Bayesian optimization (SMAC3)
- Performs feature preprocessing
- Builds weighted ensembles of best models
- Leverages meta-learning from past runs

## Installation

```bash
pip install auto-sklearn
```

**Note:** Auto-sklearn has best support on Linux and macOS. Windows support is limited.

### System Requirements

- Python 3.7+
- Linux/macOS (recommended)
- Build tools: `gcc`, `g++`, `swig`
- At least 3GB RAM

### Linux Installation

```bash
sudo apt-get install build-essential swig
pip install auto-sklearn
```

### macOS Installation

```bash
brew install swig
pip install auto-sklearn
```

## Features

### Automated Pipeline Optimization

Auto-sklearn automatically:
1. **Algorithm Selection**: Tests classification/regression algorithms
2. **Hyperparameter Tuning**: Uses Bayesian optimization (SMAC3)
3. **Feature Preprocessing**: Applies scaling, PCA, feature selection
4. **Ensemble Building**: Creates weighted ensembles of top models
5. **Meta-Learning**: Warm-starts from similar past problems

### Key Differences from TPOT

| Feature | Auto-sklearn | TPOT |
|---------|-------------|------|
| **Optimization** | Bayesian (SMAC3) | Genetic Programming |
| **Meta-Learning** | Yes (warm-start) | No |
| **Ensembling** | Weighted ensemble | Single pipeline |
| **Speed** | Faster convergence | Broader search |
| **Output** | Multiple models | Single pipeline |

## Usage

### Basic Classification

```python
from adapters.autosklearn import AutoSklearnAdapter

adapter = AutoSklearnAdapter()

input_data = {
    "task": "classification",
    "X_train": [[5.1, 3.5, 1.4, 0.2], [4.9, 3.0, 1.4, 0.2], ...],
    "y_train": [0, 0, 1, 1, 2, 2, ...],
    "time_limit": 300,  # 5 minutes
    "ensemble_size": 50,
    "metric": "accuracy"
}

result = await adapter.execute(input_data)

if result.success:
    data = result.data
    print(f"Training Score: {data['performance']['train_score']:.4f}")
    print(f"Best Model: {data['best_model']}")
    print(f"Ensemble Size: {data['ensemble']['size']}")
```

### Regression

```python
input_data = {
    "task": "regression",
    "X_train": [[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], ...],
    "y_train": [1.5, 2.5, 3.5, ...],
    "time_limit": 600,  # 10 minutes
    "metric": "r2"
}

result = await adapter.execute(input_data)
```

### With Test Set Evaluation

```python
input_data = {
    "task": "classification",
    "X_train": train_features,
    "y_train": train_labels,
    "X_test": test_features,
    "y_test": test_labels,
    "time_limit": 300,
    "ensemble_size": 50
}

result = await adapter.execute(input_data)

if result.success:
    print(f"Train Score: {result.data['performance']['train_score']:.4f}")
    print(f"Test Score: {result.data['performance']['test_score']:.4f}")
```

### Custom Resampling Strategy

```python
input_data = {
    "task": "classification",
    "X_train": features,
    "y_train": labels,
    "time_limit": 300,
    "resampling_strategy": "cv",  # Use cross-validation
    "resampling_strategy_arguments": {"folds": 5}
}

result = await adapter.execute(input_data)
```

### Include Specific Preprocessors

```python
input_data = {
    "task": "classification",
    "X_train": features,
    "y_train": labels,
    "time_limit": 300,
    "include_preprocessors": ["PCA", "StandardScaler"]
}

result = await adapter.execute(input_data)
```

### Export Trained Model

```python
input_data = {
    "task": "classification",
    "X_train": features,
    "y_train": labels,
    "time_limit": 300,
    "export_model_path": "/path/to/model.pkl"
}

result = await adapter.execute(input_data)

if result.success:
    print(f"Model saved to: {result.data['model_path']}")
```

## Input Parameters

### Required Parameters

- **`X_train`** (list/array): Feature matrix for training
- **`y_train`** (list/array): Target values for training
- **`task`** (str): `"classification"` or `"regression"` (default: `"classification"`)

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `time_limit` | int | 300 | Total time budget in seconds |
| `per_run_time_limit` | int | 30 | Time limit per model in seconds |
| `ensemble_size` | int | 50 | Number of models in ensemble |
| `metric` | str | auto | Metric to optimize (see below) |
| `memory_limit` | int | 3072 | Memory limit in MB |
| `n_jobs` | int | 1 | Number of parallel jobs (-1 = all cores) |
| `seed` | int | None | Random seed for reproducibility |
| `resampling_strategy` | str | "holdout-iterative-fit" | CV strategy |
| `resampling_strategy_arguments` | dict | None | CV configuration |
| `include_preprocessors` | list | None | Preprocessors to include |
| `X_test` | array | None | Test features |
| `y_test` | array | None | Test labels |
| `export_model_path` | str | None | Path to save model |

### Supported Metrics

**Classification:**
- `accuracy` (default)
- `roc_auc`
- `f1`
- `precision`
- `recall`
- `log_loss`

**Regression:**
- `r2` (default)
- `mse`
- `mae`
- `rmse`

### Resampling Strategies

- `holdout` - Single train/validation split
- `holdout-iterative-fit` - Iterative training (default, fastest)
- `cv` - Cross-validation
- `partial-cv` - Partial cross-validation

## Output Format

```python
{
    "model_id": "autosklearn_12345",
    "task": "classification",
    "best_model": {
        "weight": 0.45,
        "model_id": "...",
        "algorithm": "RandomForest"
    },
    "ensemble": {
        "size": 50,
        "models": [
            {"weight": 0.45, "algorithm": "RandomForest"},
            {"weight": 0.30, "algorithm": "GradientBoosting"},
            {"weight": 0.25, "algorithm": "ExtraTrees"}
        ]
    },
    "predictions": {
        "train": [0, 1, 0, 1, ...],
        "test": [0, 1, 1, ...]  # If X_test provided
    },
    "performance": {
        "train_score": 0.98,
        "test_score": 0.95,  # If X_test provided
        "metric": "accuracy"
    },
    "leaderboard": [
        {"rank": 1, "score": 0.98, "model_id": 42},
        {"rank": 2, "score": 0.97, "model_id": 15},
        ...
    ],
    "optimization_history": {
        "total_models_evaluated": 127,
        "meta_learning_used": true
    },
    "feature_importance": [0.25, 0.18, 0.15, ...],  # If available
    "configuration": {
        "time_limit": 300,
        "per_run_time_limit": 30,
        "ensemble_size": 50,
        "resampling_strategy": "holdout-iterative-fit"
    },
    "warnings": []
}
```

## Drug Discovery Use Cases

### 1. QSAR Model Building

Build predictive models for molecular properties:

```python
# Features: molecular descriptors
# Target: bioactivity (IC50, pIC50, etc.)

input_data = {
    "task": "regression",
    "X_train": molecular_descriptors,
    "y_train": pic50_values,
    "time_limit": 600,
    "metric": "r2",
    "resampling_strategy": "cv",
    "resampling_strategy_arguments": {"folds": 5}
}

result = await adapter.execute(input_data)
```

### 2. Activity Classification

Classify compounds as active/inactive:

```python
input_data = {
    "task": "classification",
    "X_train": fingerprints,
    "y_train": activity_labels,  # [0, 1, 0, 1, ...]
    "time_limit": 300,
    "metric": "roc_auc",
    "ensemble_size": 100
}

result = await adapter.execute(input_data)
```

### 3. Multi-Class Target Prediction

Predict which target a compound binds to:

```python
input_data = {
    "task": "classification",
    "X_train": compound_features,
    "y_train": target_ids,  # [0, 1, 2, 3, ...]
    "time_limit": 900,
    "metric": "f1"
}

result = await adapter.execute(input_data)
```

### 4. Ensemble Learning for Robustness

Build robust models with multiple algorithms:

```python
input_data = {
    "task": "regression",
    "X_train": descriptors,
    "y_train": solubility_values,
    "time_limit": 600,
    "ensemble_size": 100,  # Large ensemble for stability
    "per_run_time_limit": 60
}

result = await adapter.execute(input_data)

# Ensemble reduces overfitting and improves generalization
print(f"Ensemble size: {result.data['ensemble']['size']}")
print(f"Models: {[m['algorithm'] for m in result.data['ensemble']['models']]}")
```

## Advanced Usage

### Parallel Training

Use all CPU cores for faster optimization:

```python
input_data = {
    "task": "classification",
    "X_train": features,
    "y_train": labels,
    "time_limit": 600,
    "n_jobs": -1,  # Use all cores
    "per_run_time_limit": 30
}
```

### Memory-Constrained Environments

Limit memory usage:

```python
input_data = {
    "task": "classification",
    "X_train": features,
    "y_train": labels,
    "time_limit": 300,
    "memory_limit": 2048  # 2GB limit
}
```

### Reproducible Results

Set random seed for reproducibility:

```python
input_data = {
    "task": "classification",
    "X_train": features,
    "y_train": labels,
    "time_limit": 300,
    "seed": 42
}
```

## Performance Tips

### 1. Time Budget Allocation

- **Quick exploration**: 60-300 seconds
- **Standard optimization**: 300-900 seconds (5-15 minutes)
- **Thorough search**: 1800-3600 seconds (30-60 minutes)
- **Production models**: 7200+ seconds (2+ hours)

### 2. Per-Run Time Limit

- **Fast iterations**: 10-30 seconds
- **Standard**: 30-60 seconds
- **Complex models**: 60-300 seconds

### 3. Ensemble Size

- **Quick models**: 10-20 models
- **Standard**: 50 models
- **Robust ensembles**: 100+ models

### 4. Cross-Validation vs Holdout

- **`holdout-iterative-fit`**: Fastest, good for large datasets
- **`holdout`**: Fast, simple split
- **`cv`**: More robust, slower
- **`partial-cv`**: Compromise between speed and robustness

## Comparison with TPOT

### When to Use Auto-sklearn

- You want faster convergence
- Meta-learning is available (similar problems solved before)
- You prefer weighted ensembles over single pipelines
- You have limited time budget

### When to Use TPOT

- You want to explore a broader search space
- You prefer genetic programming approach
- You want a single, interpretable pipeline
- You need Windows compatibility

### Using Both (Recommended)

Run both adapters and ensemble their predictions:

```python
# Run Auto-sklearn
autosklearn_result = await autosklearn_adapter.execute(input_data)

# Run TPOT
tpot_result = await tpot_adapter.execute(input_data)

# Ensemble predictions
final_predictions = (
    0.6 * autosklearn_predictions +
    0.4 * tpot_predictions
)
```

## Troubleshooting

### ImportError: auto-sklearn not found

```bash
# Install build dependencies first (Linux)
sudo apt-get install build-essential swig

# Install auto-sklearn
pip install auto-sklearn
```

### Memory Limit Errors

Reduce `memory_limit` parameter:

```python
input_data = {
    ...
    "memory_limit": 2048  # Reduce from 3072
}
```

### Time Limit Too Short

Increase `time_limit` and reduce `per_run_time_limit`:

```python
input_data = {
    ...
    "time_limit": 900,
    "per_run_time_limit": 20
}
```

### No Ensemble Built

Increase `time_limit` to allow more models to be evaluated:

```python
input_data = {
    ...
    "time_limit": 600,  # More time for more models
    "ensemble_size": 25  # Reduce ensemble size requirement
}
```

### Windows Compatibility Issues

Auto-sklearn has limited Windows support. Consider:
- Using Windows Subsystem for Linux (WSL)
- Using Docker container with Linux
- Using TPOT adapter instead (better Windows support)

## Integration with PharmForge

### Register Adapter

Add to `backend/core/adapter_registry.py`:

```python
from adapters.autosklearn import AutoSklearnAdapter

def register_all_adapters():
    registry = get_registry()

    # ... other adapters ...

    # Auto-sklearn adapter
    registry.register(AutoSklearnAdapter())

    logger.info("All adapters registered")
```

### Use in Pipeline

```python
from backend.core.adapter_registry import get_registry

registry = get_registry()
autosklearn = registry.get("autosklearn")

# Prepare molecular features
from rdkit import Chem
from rdkit.Chem import Descriptors

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
    ]

X_train = [compute_descriptors(s) for s in smiles_list]
y_train = activity_values

result = await autosklearn(
    {
        "task": "regression",
        "X_train": X_train,
        "y_train": y_train,
        "time_limit": 300
    }
)
```

## References

- **Auto-sklearn Paper**: [Efficient and Robust Automated Machine Learning (NeurIPS 2015)](https://papers.nips.cc/paper/2015/hash/11d0e6287202fced83f79975ec59a3a6-Abstract.html)
- **Documentation**: https://automl.github.io/auto-sklearn/
- **GitHub**: https://github.com/automl/auto-sklearn
- **SMAC3**: https://github.com/automl/SMAC3

## License

Auto-sklearn is licensed under the BSD 3-Clause License (commercial friendly).

## Version

Adapter Version: 1.0.0
Auto-sklearn Version: 0.15+
Last Updated: 2025-10-30
