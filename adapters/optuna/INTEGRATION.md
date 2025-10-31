# Optuna Adapter Integration Guide

## Quick Start

### 1. Installation

```bash
pip install optuna
```

### 2. Basic Usage

```python
from adapters.optuna import OptunaAdapter
import asyncio

async def main():
    # Initialize adapter
    adapter = OptunaAdapter()

    # Define objective function
    def objective(params, trial):
        # Your optimization logic here
        x = params["x"]
        return (x - 2) ** 2  # Minimize

    # Configure optimization
    config = {
        "search_space": {
            "x": {"type": "float", "low": -10.0, "high": 10.0}
        },
        "objective_function": "custom",
        "n_trials": 20,
        "direction": "minimize"
    }

    # Run optimization
    result = await adapter.execute(config, custom_objective=objective)

    if result.success:
        print(f"Best value: {result.data['best_params']['x']}")

asyncio.run(main())
```

## Integration with PharmForge Adapters

### Optimizing Molecular Docking Parameters

```python
from adapters.vina import VinaAdapter
from adapters.optuna import OptunaAdapter

async def optimize_docking():
    vina = VinaAdapter()
    optuna = OptunaAdapter()

    async def docking_objective(params, trial):
        """Optimize Vina docking parameters"""
        result = await vina.execute(
            receptor_pdbqt="protein.pdbqt",
            ligand_pdbqt="ligand.pdbqt",
            exhaustiveness=params["exhaustiveness"],
            num_modes=params["num_modes"],
            energy_range=params["energy_range"]
        )

        if result.success:
            return result.data["best_affinity"]
        return float('inf')

    config = {
        "search_space": {
            "exhaustiveness": {"type": "int", "low": 8, "high": 32},
            "num_modes": {"type": "int", "low": 5, "high": 20},
            "energy_range": {"type": "int", "low": 3, "high": 10}
        },
        "n_trials": 30,
        "direction": "minimize"
    }

    result = await optuna.execute(config, custom_objective=docking_objective)
    return result.data['best_params']
```

### Optimizing ADMET Prediction Models

```python
from adapters.admet_ai import ADMETAdapter
from adapters.optuna import OptunaAdapter
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

async def optimize_admet_model():
    optuna = OptunaAdapter()

    def model_objective(params, trial):
        """Optimize random forest for ADMET prediction"""
        # Create model with suggested parameters
        model = RandomForestClassifier(
            n_estimators=params["n_estimators"],
            max_depth=params["max_depth"],
            min_samples_split=params["min_samples_split"],
            random_state=42
        )

        # Evaluate with cross-validation
        scores = cross_val_score(model, X_train, y_train, cv=5)
        return -scores.mean()  # Negative because we minimize

    config = {
        "search_space": {
            "n_estimators": {"type": "int", "low": 50, "high": 500},
            "max_depth": {"type": "int", "low": 3, "high": 20},
            "min_samples_split": {"type": "int", "low": 2, "high": 10}
        },
        "n_trials": 50,
        "direction": "minimize"
    }

    result = await optuna.execute(config, custom_objective=model_objective)
    return result.data['best_params']
```

### Optimizing ChemProp Model Hyperparameters

```python
from adapters.chemprop import ChemPropAdapter
from adapters.optuna import OptunaAdapter

async def optimize_chemprop():
    chemprop = ChemPropAdapter()
    optuna = OptunaAdapter()

    async def chemprop_objective(params, trial):
        """Optimize ChemProp training parameters"""
        result = await chemprop.execute(
            smiles_list=training_smiles,
            targets=training_targets,
            config={
                "learning_rate": params["learning_rate"],
                "hidden_size": params["hidden_size"],
                "depth": params["depth"],
                "dropout": params["dropout"],
                "epochs": 50  # Fixed for faster optimization
            }
        )

        if result.success:
            return result.data["validation_rmse"]
        return float('inf')

    config = {
        "search_space": {
            "learning_rate": {
                "type": "float",
                "low": 1e-5,
                "high": 1e-2,
                "log": True
            },
            "hidden_size": {
                "type": "int",
                "low": 128,
                "high": 1024
            },
            "depth": {
                "type": "int",
                "low": 2,
                "high": 6
            },
            "dropout": {
                "type": "float",
                "low": 0.0,
                "high": 0.5
            }
        },
        "n_trials": 100,
        "direction": "minimize",
        "sampler": "TPE"
    }

    result = await optuna.execute(config, custom_objective=chemprop_objective)
    return result.data
```

## Advanced Usage Patterns

### Multi-Stage Optimization

Optimize in multiple stages for complex problems:

```python
async def multi_stage_optimization():
    optuna = OptunaAdapter()

    # Stage 1: Coarse search
    coarse_config = {
        "search_space": {
            "param1": {"type": "float", "low": 0.0, "high": 100.0},
            "param2": {"type": "int", "low": 1, "high": 100}
        },
        "n_trials": 20,
        "sampler": "Random"
    }

    coarse_result = await optuna.execute(coarse_config, custom_objective=objective)

    # Stage 2: Fine-tune around best parameters
    best_param1 = coarse_result.data['best_params']['param1']
    fine_config = {
        "search_space": {
            "param1": {
                "type": "float",
                "low": max(0.0, best_param1 - 10),
                "high": min(100.0, best_param1 + 10)
            },
            "param2": {"type": "int", "low": 1, "high": 100}
        },
        "n_trials": 50,
        "sampler": "TPE"
    }

    fine_result = await optuna.execute(fine_config, custom_objective=objective)
    return fine_result
```

### Ensemble Optimization

Optimize parameters for multiple models simultaneously:

```python
async def optimize_ensemble():
    optuna = OptunaAdapter()

    def ensemble_objective(params, trial):
        """Optimize weights and parameters for ensemble"""
        # Train multiple models
        models = []
        weights = []

        for i in range(3):
            weight = params[f"weight_{i}"]
            model_param = params[f"model_{i}_param"]

            model = train_model(model_param)
            models.append(model)
            weights.append(weight)

        # Evaluate ensemble
        predictions = ensemble_predict(models, weights, validation_data)
        return calculate_loss(predictions, validation_labels)

    config = {
        "search_space": {
            "weight_0": {"type": "float", "low": 0.0, "high": 1.0},
            "weight_1": {"type": "float", "low": 0.0, "high": 1.0},
            "weight_2": {"type": "float", "low": 0.0, "high": 1.0},
            "model_0_param": {"type": "float", "low": 0.01, "high": 10.0},
            "model_1_param": {"type": "float", "low": 0.01, "high": 10.0},
            "model_2_param": {"type": "float", "low": 0.01, "high": 10.0}
        },
        "n_trials": 100,
        "direction": "minimize"
    }

    result = await optuna.execute(config, custom_objective=ensemble_objective)
    return result.data
```

## Best Practices

### 1. Start Small, Scale Up

```python
# First, test with few trials
test_config = {**config, "n_trials": 10}
test_result = await adapter.execute(test_config, custom_objective=objective)

# If successful, run full optimization
full_config = {**config, "n_trials": 100}
full_result = await adapter.execute(full_config, custom_objective=objective)
```

### 2. Use Appropriate Samplers

- **TPE**: Best for most cases (< 100 dimensions)
- **Random**: Good baseline, parallel-friendly
- **Grid**: Exhaustive search for small spaces

### 3. Log Search Space

Use logarithmic scale for parameters spanning orders of magnitude:

```python
"learning_rate": {
    "type": "float",
    "low": 1e-5,
    "high": 1e-1,
    "log": True  # Important for learning rates
}
```

### 4. Handle Failures Gracefully

```python
def robust_objective(params, trial):
    try:
        result = expensive_computation(params)
        return result
    except Exception as e:
        logger.warning(f"Trial {trial.number} failed: {e}")
        return float('inf')  # Return worst possible value
```

### 5. Save Results

```python
result = await optuna.execute(config, custom_objective=objective)

if result.success:
    # Save best parameters
    import json
    with open('best_params.json', 'w') as f:
        json.dump(result.data['best_params'], f, indent=2)

    # Save trial history for visualization
    with open('trial_history.json', 'w') as f:
        json.dump(result.data['trial_history'], f, indent=2)
```

## Troubleshooting

### ImportError: optuna not installed

```bash
pip install optuna
```

### Optimization taking too long

```python
# Reduce number of trials
config["n_trials"] = 20  # Instead of 100

# Or use timeout
config["timeout"] = 300  # 5 minutes
```

### Poor optimization results

1. Check search space bounds
2. Try different samplers
3. Increase number of trials
4. Verify objective function correctness

## Performance Tips

1. **Parallelize when possible**: Optuna supports parallel optimization (future enhancement)
2. **Use pruning**: Enable pruning for expensive objectives
3. **Warm start**: Use previous study results as starting point
4. **Cache objective evaluations**: Avoid redundant computations

## Support

For issues or questions:
- Check the README.md for detailed documentation
- Review example_usage.py for code examples
- Consult Optuna documentation: https://optuna.readthedocs.io/
