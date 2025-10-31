# TPOT Adapter

## Overview

The TPOT (Tree-based Pipeline Optimization Tool) adapter provides automated machine learning pipeline optimization using genetic programming. TPOT automatically selects the best machine learning pipeline for your data by exploring thousands of possible configurations.

## Features

- **Automated Pipeline Optimization**: Automatically discovers the best ML pipeline for your data
- **Classification & Regression**: Supports both classification and regression tasks
- **Feature Engineering**: Automatically includes feature preprocessing and selection
- **Model Selection**: Tests multiple ML algorithms (Random Forests, XGBoost, SVM, etc.)
- **Hyperparameter Tuning**: Optimizes hyperparameters using genetic programming
- **Pipeline Export**: Export optimized pipeline as executable Python code
- **Multi-objective Optimization**: Balances accuracy with model complexity
- **Cross-validation**: Built-in cross-validation for robust evaluation

## Installation

```bash
pip install tpot scikit-learn numpy
```

Optional dependencies for enhanced functionality:
```bash
pip install xgboost  # For XGBoost models
pip install dask distributed  # For parallel processing
```

## Input Format

```python
{
    "X_train": [[...], [...]],  # Feature matrix (2D array)
    "y_train": [...],           # Target values (1D array)
    "task": "classification",   # or "regression"
    "generations": 5,           # Number of generations (default: 5)
    "population_size": 20,      # Population size (default: 20)
    "scoring": "roc_auc",       # Scoring metric (optional)
    "cv": 5,                    # Cross-validation folds (default: 5)
    "random_state": 42,         # Random seed (optional)
    "verbosity": 2,             # Verbosity level 0-3 (default: 2)
    "max_time_mins": 60,        # Maximum runtime in minutes (optional)
    "max_eval_time_mins": 5,    # Max time per pipeline eval (optional)
    "n_jobs": -1,               # Number of parallel jobs (default: -1)
    "config_dict": None,        # TPOT configuration preset (optional)
    "X_test": [[...], [...]],   # Test features (optional)
    "y_test": [...]             # Test targets (optional)
}
```

## Output Structure

```python
{
    "best_pipeline": "Pipeline(steps=[...])",  # String representation
    "best_score": 0.95,                        # Best CV score
    "test_score": 0.93,                        # Test score (if provided)
    "pipeline_code": "# Python code...",       # Exportable code
    "fitted_pipeline": pipeline_object,        # Fitted sklearn pipeline
    "pipeline_steps": [                        # List of pipeline steps
        {
            "name": "preprocessor",
            "type": "StandardScaler",
            "parameters": {...}
        },
        {
            "name": "classifier",
            "type": "RandomForestClassifier",
            "parameters": {...}
        }
    ],
    "optimization_history": [                  # Generation-by-generation stats
        {
            "generation": 0,
            "pipelines_evaluated": 20,
            "best_score": 0.85,
            "mean_score": 0.78,
            "std_score": 0.05
        }
    ],
    "metadata": {
        "task": "classification",
        "generations_evaluated": 5,
        "population_size": 20,
        "pipelines_tested": 100,
        "best_generation": 3,
        "best_generation_score": 0.95,
        "scoring": "roc_auc",
        "cv_folds": 5,
        "n_samples": 1000,
        "n_features": 20,
        "n_test_samples": 200
    }
}
```

## Usage Examples

### Basic Classification

```python
from adapters.tpot import TPOTAdapter
import numpy as np

# Create adapter
adapter = TPOTAdapter()

# Prepare data
X_train = np.array([[...], [...]])
y_train = np.array([0, 1, 0, 1, ...])

# Configure TPOT
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification",
    "generations": 5,
    "population_size": 20,
    "scoring": "roc_auc"
}

# Run optimization
result = await adapter.execute(input_data)

if result.success:
    print(f"Best CV score: {result.data['best_score']:.4f}")
    print(f"Pipeline: {result.data['best_pipeline']}")

    # Export optimized code
    print(result.data['pipeline_code'])
```

### Regression Task

```python
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "task": "regression",
    "generations": 10,
    "population_size": 50,
    "scoring": "neg_mean_squared_error",
    "cv": 10
}

result = await adapter.execute(input_data)
```

### With Test Set Evaluation

```python
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "X_test": X_test,
    "y_test": y_test,
    "task": "classification",
    "generations": 5,
    "population_size": 20
}

result = await adapter.execute(input_data)

if result.success:
    print(f"Training score: {result.data['best_score']:.4f}")
    print(f"Test score: {result.data['test_score']:.4f}")
```

### Using Configuration Presets

TPOT supports several configuration presets:

```python
# TPOT light - Faster, simpler models
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification",
    "config_dict": "TPOT light"
}

# TPOT MDR - For feature construction
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification",
    "config_dict": "TPOT MDR"
}

# TPOT sparse - For sparse data
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification",
    "config_dict": "TPOT sparse"
}
```

### Time-constrained Optimization

```python
input_data = {
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification",
    "max_time_mins": 30,           # Total runtime limit
    "max_eval_time_mins": 2,       # Per-pipeline time limit
    "generations": None,           # Let time limit control
    "population_size": 20
}
```

## Scoring Metrics

### Classification Metrics
- `accuracy` - Accuracy score
- `roc_auc` - ROC AUC score (default for binary classification)
- `roc_auc_ovr` - ROC AUC for multiclass
- `f1` - F1 score
- `f1_weighted` - Weighted F1 score
- `precision` - Precision score
- `recall` - Recall score
- `balanced_accuracy` - Balanced accuracy

### Regression Metrics
- `neg_mean_squared_error` - Negative MSE (default for regression)
- `neg_mean_absolute_error` - Negative MAE
- `r2` - R² score
- `neg_root_mean_squared_error` - Negative RMSE
- `neg_median_absolute_error` - Negative median AE

## Integration with Other Adapters

### TPOT + Optuna (Nested Optimization)

```python
# Use Optuna to optimize TPOT hyperparameters
from adapters.optuna import OptunaAdapter
from adapters.tpot import TPOTAdapter

optuna_adapter = OptunaAdapter()
tpot_adapter = TPOTAdapter()

# Define search space for TPOT parameters
search_space = {
    "generations": {"type": "int", "low": 3, "high": 10},
    "population_size": {"type": "int", "low": 10, "high": 50},
}

# Objective function using TPOT
async def objective(params, trial):
    input_data = {
        "X_train": X_train,
        "y_train": y_train,
        "task": "classification",
        "generations": params["generations"],
        "population_size": params["population_size"],
        "verbosity": 0
    }
    result = await tpot_adapter.execute(input_data)
    return result.data["best_score"] if result.success else 0.0

# Optimize
optuna_input = {
    "search_space": search_space,
    "objective_function": objective,
    "n_trials": 20,
    "direction": "maximize"
}

result = await optuna_adapter.execute(
    optuna_input,
    custom_objective=objective
)
```

### TPOT + Feature Engineering Adapters

```python
# Use Mordred or scikit-mol for feature generation, then TPOT
from adapters.mordred import MordredAdapter
from adapters.tpot import TPOTAdapter

# Generate molecular descriptors
mordred = MordredAdapter()
features_result = await mordred.execute({
    "smiles": ["CCO", "c1ccccc1", ...],
    "descriptors": ["all"]
})

X = features_result.data["descriptors"]
y = [0, 1, 0, ...]  # Your labels

# Optimize ML pipeline with TPOT
tpot = TPOTAdapter()
result = await tpot.execute({
    "X_train": X,
    "y_train": y,
    "task": "classification"
})
```

## Best Practices

### 1. Start Small
Begin with small generations and population sizes to get quick results:
```python
{
    "generations": 3,
    "population_size": 10,
    "max_eval_time_mins": 1
}
```

### 2. Use Time Limits
Instead of fixed generations, use time limits for better control:
```python
{
    "max_time_mins": 60,
    "max_eval_time_mins": 5
}
```

### 3. Choose Appropriate Scoring
Select metrics that match your problem:
- Imbalanced classification: `f1_weighted`, `balanced_accuracy`
- Binary classification: `roc_auc`
- Regression: `r2`, `neg_mean_squared_error`

### 4. Use Configuration Presets
For faster results, use `TPOT light`:
```python
{"config_dict": "TPOT light"}
```

### 5. Set Random State
For reproducible results:
```python
{"random_state": 42}
```

### 6. Monitor Progress
Use verbosity to track optimization:
```python
{"verbosity": 2}  # 0=silent, 1=minimal, 2=standard, 3=debug
```

## Advantages Over Manual ML

1. **Automation**: No need to manually try different algorithms and parameters
2. **Feature Engineering**: Automatically includes preprocessing steps
3. **Exploration**: Tests thousands of configurations
4. **Reproducibility**: Export exact pipeline as code
5. **Best Practices**: Includes cross-validation automatically
6. **Model Diversity**: Tests wide range of algorithms

## Limitations

1. **Computational Cost**: Can be time-consuming for large datasets
2. **Memory Usage**: Stores population of pipelines in memory
3. **Black Box**: Genetic programming is less interpretable
4. **Overfitting Risk**: May overfit on small datasets despite CV
5. **No Deep Learning**: Focuses on traditional ML algorithms

## Performance Tips

1. **Use `n_jobs=-1`**: Enables parallel processing
2. **Reduce CV folds**: Use `cv=3` for faster evaluation
3. **Limit search space**: Use configuration presets
4. **Sample large datasets**: Use subset for optimization
5. **Set time limits**: Prevent excessive runtime

## Troubleshooting

### TPOT not found
```bash
pip install tpot
```

### Out of memory
- Reduce `population_size`
- Reduce dataset size
- Use `config_dict="TPOT light"`

### Taking too long
- Reduce `generations`
- Reduce `population_size`
- Set `max_time_mins`
- Use `config_dict="TPOT light"`

### Poor performance
- Increase `generations`
- Increase `population_size`
- Try different `scoring` metric
- Check data quality

## References

- [TPOT Documentation](http://epistasislab.github.io/tpot/)
- [TPOT Paper](https://doi.org/10.1007/978-3-319-31204-0_9)
- [GitHub Repository](https://github.com/EpistasisLab/tpot)

## Version History

- **1.0.0** (2024): Initial implementation with classification and regression support
