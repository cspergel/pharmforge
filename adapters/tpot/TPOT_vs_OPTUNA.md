# TPOT vs Optuna: When to Use Which?

## Overview

PharmForge includes two powerful optimization adapters: **TPOT** and **Optuna**. While both optimize machine learning workflows, they serve different purposes and use different approaches.

## Quick Comparison

| Aspect | TPOT | Optuna |
|--------|------|--------|
| **Algorithm** | Genetic Programming | Bayesian Optimization (TPE) |
| **Optimizes** | Complete ML pipelines | Hyperparameters only |
| **Output** | Trained model + code | Best parameters |
| **Model Selection** | Automatic | Manual |
| **Feature Engineering** | Included | Manual |
| **Preprocessing** | Automatic | Manual |
| **Search Space** | Predefined (sklearn) | Fully customizable |
| **Interpretability** | Low | Medium |
| **Speed** | Slower | Faster |
| **Best For** | End-to-end automation | Fine-tuning |
| **Expertise Required** | None | Some ML knowledge |

## When to Use TPOT

### ✅ Use TPOT When:

1. **You don't know which algorithm to use**
   - TPOT will test dozens of algorithms automatically
   - Good for exploratory analysis

2. **You want complete automation**
   - TPOT handles everything: preprocessing, feature selection, model selection
   - Minimal ML expertise required

3. **You need a complete pipeline**
   - TPOT outputs a ready-to-deploy sklearn pipeline
   - Includes all preprocessing steps

4. **You want exportable code**
   - TPOT exports the optimized pipeline as Python code
   - Great for reproducibility and deployment

5. **You're starting a new project**
   - TPOT quickly establishes a baseline
   - Explores the solution space broadly

6. **Feature engineering is needed**
   - TPOT automatically includes preprocessing
   - Tests different feature transformations

### Example TPOT Use Cases:

```python
# Use Case 1: New QSAR model, don't know best algorithm
tpot = TPOTAdapter()
result = await tpot.execute({
    "X_train": molecular_descriptors,
    "y_train": activity_labels,
    "task": "classification"
})
# TPOT tries Random Forests, SVM, XGBoost, etc. automatically

# Use Case 2: Need complete pipeline with preprocessing
result = await tpot.execute({
    "X_train": raw_features,  # Unnormalized, mixed scales
    "y_train": labels,
    "task": "classification"
})
# TPOT includes StandardScaler, PCA, feature selection, etc.

# Use Case 3: Quick baseline for new dataset
result = await tpot.execute({
    "X_train": X,
    "y_train": y,
    "config_dict": "TPOT light",  # Fast mode
    "generations": 3
})
```

## When to Use Optuna

### ✅ Use Optuna When:

1. **You know which algorithm to use**
   - You want Random Forest, just need best hyperparameters
   - Already decided on model architecture

2. **You need fine-grained control**
   - Custom search spaces
   - Domain-specific constraints
   - Non-standard optimizations

3. **You want faster optimization**
   - Bayesian optimization is more efficient than genetic programming
   - Fewer evaluations needed

4. **You have custom objectives**
   - Multi-objective optimization
   - Complex custom metrics
   - Business constraints

5. **You're refining an existing model**
   - Already have a working pipeline
   - Want to squeeze out extra performance

6. **You need advanced features**
   - Pruning of unpromising trials
   - Parallel optimization
   - Distributed optimization

### Example Optuna Use Cases:

```python
# Use Case 1: Optimize known Random Forest model
optuna = OptunaAdapter()
result = await optuna.execute({
    "search_space": {
        "n_estimators": {"type": "int", "low": 50, "high": 500},
        "max_depth": {"type": "int", "low": 3, "high": 20},
        "min_samples_split": {"type": "int", "low": 2, "high": 20}
    },
    "objective_function": rf_objective,
    "n_trials": 100
})

# Use Case 2: Optimize custom neural network
result = await optuna.execute({
    "search_space": {
        "learning_rate": {"type": "float", "low": 1e-5, "high": 1e-2, "log": True},
        "batch_size": {"type": "categorical", "choices": [32, 64, 128, 256]},
        "dropout": {"type": "float", "low": 0.1, "high": 0.5}
    },
    "objective_function": nn_objective,
    "n_trials": 50
})

# Use Case 3: Multi-objective optimization (accuracy + speed)
async def multi_objective(params, trial):
    model = train_model(params)
    accuracy = evaluate_accuracy(model)
    inference_time = measure_speed(model)
    return accuracy, inference_time  # Returns tuple
```

## When to Use Both Together

### 🔥 Use TPOT + Optuna for Maximum Performance:

The most powerful approach is to use them sequentially or nested:

### Strategy 1: Sequential Optimization

```python
# Step 1: Use TPOT to find best pipeline architecture
tpot = TPOTAdapter()
tpot_result = await tpot.execute({
    "X_train": X,
    "y_train": y,
    "task": "classification"
})

# Step 2: Use Optuna to fine-tune TPOT's hyperparameters
# Extract the best model type from TPOT
best_model = tpot_result.data['fitted_pipeline']

# Now optimize that specific model with Optuna
optuna = OptunaAdapter()
# ... optimize the TPOT-discovered pipeline
```

### Strategy 2: Meta-Optimization (Optuna optimizes TPOT)

```python
# Use Optuna to find best TPOT configuration
optuna = OptunaAdapter()
tpot = TPOTAdapter()

async def objective(params, trial):
    # Run TPOT with different configurations
    result = await tpot.execute({
        "X_train": X,
        "y_train": y,
        "generations": params["generations"],
        "population_size": params["population_size"],
        "verbosity": 0
    })
    return result.data["best_score"]

result = await optuna.execute({
    "search_space": {
        "generations": {"type": "int", "low": 3, "high": 10},
        "population_size": {"type": "int", "low": 10, "high": 50}
    },
    "objective_function": objective,
    "n_trials": 20
}, custom_objective=objective)

# Get best TPOT configuration
best_tpot_config = result.data["best_params"]
```

### Strategy 3: Parallel Exploration

```python
# Run TPOT and Optuna in parallel, compare results
async def optimize_both():
    # TPOT: Broad exploration
    tpot_result = await tpot.execute({
        "X_train": X,
        "y_train": y,
        "task": "classification",
        "generations": 5
    })

    # Optuna: Optimize specific model (e.g., XGBoost)
    async def xgb_objective(params, trial):
        # ... XGBoost with params
        return score

    optuna_result = await optuna.execute({
        "search_space": xgb_search_space,
        "objective_function": xgb_objective,
        "n_trials": 100
    })

    # Compare and choose best
    if tpot_result.data["best_score"] > optuna_result.data["best_value"]:
        return tpot_result
    return optuna_result
```

## Decision Tree

```
Do you know which algorithm to use?
│
├─ NO → Use TPOT
│       └─ TPOT will explore and find the best algorithm
│
└─ YES → Is it a standard sklearn algorithm?
         │
         ├─ NO (custom model/neural network) → Use Optuna
         │                                       └─ Full control over search space
         │
         └─ YES → Do you need preprocessing/feature engineering?
                  │
                  ├─ YES → Use TPOT
                  │        └─ Automatic pipeline construction
                  │
                  └─ NO → Use Optuna
                           └─ Faster hyperparameter optimization
```

## Practical Examples

### Example 1: New Drug Discovery Project

**Problem**: Building QSAR model for new target, no prior knowledge

**Solution**: Use TPOT
```python
# TPOT explores all possibilities
tpot = TPOTAdapter()
result = await tpot.execute({
    "X_train": mordred_descriptors,
    "y_train": activity_labels,
    "task": "classification",
    "generations": 10,
    "population_size": 50
})

# Get complete pipeline
print(result.data['pipeline_code'])  # Ready to deploy!
```

**Why TPOT**: No idea which algorithm works best, need exploration

---

### Example 2: Optimizing Existing Random Forest Model

**Problem**: Have working RF model, want better hyperparameters

**Solution**: Use Optuna
```python
# Already know Random Forest works, just optimize it
optuna = OptunaAdapter()
result = await optuna.execute({
    "search_space": {
        "n_estimators": {"type": "int", "low": 100, "high": 1000},
        "max_depth": {"type": "int", "low": 5, "high": 30},
        "min_samples_leaf": {"type": "int", "low": 1, "high": 10}
    },
    "objective_function": rf_objective,
    "n_trials": 100,
    "sampler": "TPE"
})
```

**Why Optuna**: Known algorithm, just need hyperparameter tuning

---

### Example 3: Virtual Screening Pipeline

**Problem**: Need best possible model for screening millions of compounds

**Solution**: Use TPOT then Optuna
```python
# Phase 1: TPOT finds best pipeline (exploratory)
tpot_result = await tpot.execute({
    "X_train": X_train,
    "y_train": y_train,
    "task": "classification"
})

# Phase 2: Optuna fine-tunes the TPOT pipeline
# Extract model from TPOT and create Optuna objective
# ... detailed optimization

# Phase 3: Deploy the best of both
final_model = best_of(tpot_result, optuna_result)
```

**Why Both**: Maximum performance needed for critical application

---

### Example 4: Multi-Target Optimization

**Problem**: Optimize for multiple related targets

**Solution**: TPOT for each target, then ensemble
```python
# Use TPOT to find best pipeline for each target
targets = ["target_A", "target_B", "target_C"]
models = {}

for target in targets:
    tpot = TPOTAdapter()
    result = await tpot.execute({
        "X_train": X,
        "y_train": labels[target],
        "task": "classification"
    })
    models[target] = result.data['fitted_pipeline']
```

**Why TPOT**: Each target may need different algorithm/preprocessing

## Performance Comparison

### Speed

| Dataset Size | TPOT Time | Optuna Time | Winner |
|--------------|-----------|-------------|--------|
| Small (100s) | 5-10 min | 1-2 min | Optuna |
| Medium (1000s) | 30-60 min | 5-10 min | Optuna |
| Large (10000s+) | Hours | 30-60 min | Optuna |

**Note**: TPOT is slower because it explores more (algorithms + hyperparameters)

### Performance

| Scenario | TPOT Score | Optuna Score | Winner |
|----------|------------|--------------|--------|
| Unknown problem | 0.85 | 0.82 | TPOT |
| Known algorithm | 0.87 | 0.89 | Optuna |
| With preprocessing | 0.88 | 0.84 | TPOT |
| Fine-tuning | 0.86 | 0.91 | Optuna |

**Note**: Results vary by dataset and problem

## Cost Comparison

### Computational Cost

- **TPOT**: Higher (tests many algorithms × hyperparameters)
- **Optuna**: Lower (optimizes single algorithm)

### Time Cost

- **TPOT**: Hours for thorough optimization
- **Optuna**: Minutes to hours

### Expertise Cost

- **TPOT**: Low (little ML knowledge needed)
- **Optuna**: Medium (need to define search space and objective)

## Summary

### Choose TPOT if:
- ✅ Starting new project
- ✅ Don't know best algorithm
- ✅ Need complete pipeline
- ✅ Want exportable code
- ✅ Minimal ML expertise
- ✅ Feature engineering needed

### Choose Optuna if:
- ✅ Known algorithm
- ✅ Need fine control
- ✅ Want speed
- ✅ Custom objectives
- ✅ Refining existing model
- ✅ Advanced features needed

### Use Both if:
- ✅ Need maximum performance
- ✅ Have time for thorough optimization
- ✅ Critical application
- ✅ Want to compare approaches

## Getting Started

### TPOT Quick Start
```python
from adapters.tpot import TPOTAdapter
adapter = TPOTAdapter()
result = await adapter.execute({
    "X_train": X,
    "y_train": y,
    "task": "classification"
})
```

### Optuna Quick Start
```python
from adapters.optuna import OptunaAdapter
adapter = OptunaAdapter()
result = await adapter.execute({
    "search_space": {...},
    "objective_function": objective,
    "n_trials": 100
})
```

### Combined Approach
See `INTEGRATION.md` for detailed examples of using TPOT and Optuna together.

## Resources

- **TPOT Documentation**: `adapters/tpot/README.md`
- **Optuna Documentation**: `adapters/optuna/README.md`
- **Integration Guide**: `adapters/tpot/INTEGRATION.md`
- **TPOT Examples**: `adapters/tpot/example_usage.py`
- **Optuna Examples**: `adapters/optuna/example_usage.py`
