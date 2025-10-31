# Optuna Adapter

Hyperparameter optimization adapter using the Optuna framework for PharmForge.

## Overview

The Optuna adapter provides automated hyperparameter tuning capabilities for optimizing adapter parameters, model configurations, and other numerical parameters. It supports various sampling strategies and can be used to optimize any objective function.

## Installation

```bash
pip install optuna
```

## Features

- **Multiple Samplers**: TPE (Tree-structured Parzen Estimator), Random, Grid Search
- **Flexible Search Spaces**: Support for float, int, and categorical parameters
- **Trial History**: Complete optimization history with all trials
- **Pruning Support**: Early stopping for unpromising trials
- **Minimize/Maximize**: Support for both minimization and maximization objectives

## Usage

### Basic Example

```python
from adapters.optuna import OptunaAdapter

# Initialize adapter
adapter = OptunaAdapter()

# Define optimization configuration
config = {
    "search_space": {
        "learning_rate": {
            "type": "float",
            "low": 1e-5,
            "high": 1e-1,
            "log": True
        },
        "batch_size": {
            "type": "int",
            "low": 16,
            "high": 128
        },
        "optimizer": {
            "type": "categorical",
            "choices": ["adam", "sgd", "rmsprop"]
        }
    },
    "objective_function": "custom",
    "n_trials": 50,
    "sampler": "TPE",
    "direction": "minimize"
}

# Define custom objective function
def my_objective(params, trial):
    """
    Custom objective function to optimize

    Args:
        params: Dictionary of suggested parameters
        trial: Optuna trial object (for reporting intermediate values)

    Returns:
        float: Objective value to minimize/maximize
    """
    # Your optimization logic here
    # For example, train a model and return validation loss
    lr = params["learning_rate"]
    batch_size = params["batch_size"]
    optimizer = params["optimizer"]

    # Simulate training
    score = train_and_evaluate(lr, batch_size, optimizer)

    return score

# Run optimization
result = await adapter.execute(config, custom_objective=my_objective)

if result.success:
    print(f"Best parameters: {result.data['best_params']}")
    print(f"Best value: {result.data['best_value']}")
    print(f"Number of trials: {result.data['n_trials']}")
else:
    print(f"Error: {result.error}")
```

### Search Space Configuration

#### Float Parameters

```python
"parameter_name": {
    "type": "float",
    "low": 0.001,
    "high": 1.0,
    "log": True  # Use log scale (optional)
}
```

#### Integer Parameters

```python
"parameter_name": {
    "type": "int",
    "low": 1,
    "high": 100,
    "log": False  # Use log scale (optional)
}
```

#### Categorical Parameters

```python
"parameter_name": {
    "type": "categorical",
    "choices": ["option1", "option2", "option3"]
}
```

### Sampler Options

#### TPE Sampler (Default)

Tree-structured Parzen Estimator - efficient Bayesian optimization

```python
config = {
    "sampler": "TPE",
    "sampler_kwargs": {
        "n_startup_trials": 10,  # Random trials before TPE
        "n_ei_candidates": 24,   # Number of candidates
        "seed": 42               # Random seed
    }
}
```

#### Random Sampler

Random search over the parameter space

```python
config = {
    "sampler": "Random",
    "sampler_kwargs": {
        "seed": 42
    }
}
```

#### Grid Sampler

Exhaustive grid search

```python
config = {
    "sampler": "Grid",
    "sampler_kwargs": {
        "search_space": {
            "learning_rate": [0.001, 0.01, 0.1],
            "batch_size": [16, 32, 64]
        }
    }
}
```

### Result Structure

```python
{
    "best_params": {
        "learning_rate": 0.01,
        "batch_size": 32,
        "optimizer": "adam"
    },
    "best_value": 0.1234,
    "best_trial": 42,
    "n_trials": 50,
    "trial_history": [
        {
            "number": 0,
            "value": 0.5,
            "params": {...},
            "state": "COMPLETE",
            "datetime_start": "2025-01-01T00:00:00",
            "datetime_complete": "2025-01-01T00:01:00",
            "duration": 60.0
        },
        ...
    ],
    "sampler": "TPE",
    "direction": "minimize"
}
```

## Use Cases

### 1. Optimizing Docking Parameters

```python
config = {
    "search_space": {
        "exhaustiveness": {"type": "int", "low": 8, "high": 32},
        "num_modes": {"type": "int", "low": 1, "high": 20},
        "energy_range": {"type": "int", "low": 3, "high": 10}
    },
    "n_trials": 30,
    "direction": "minimize"  # Minimize binding energy
}
```

### 2. Optimizing ML Model Hyperparameters

```python
config = {
    "search_space": {
        "learning_rate": {"type": "float", "low": 1e-5, "high": 1e-2, "log": True},
        "hidden_dim": {"type": "int", "low": 64, "high": 512},
        "dropout": {"type": "float", "low": 0.0, "high": 0.5},
        "activation": {"type": "categorical", "choices": ["relu", "tanh", "elu"]}
    },
    "n_trials": 100,
    "direction": "maximize"  # Maximize accuracy
}
```

### 3. Optimizing ADMET Prediction Parameters

```python
config = {
    "search_space": {
        "n_estimators": {"type": "int", "low": 50, "high": 500},
        "max_depth": {"type": "int", "low": 3, "high": 20},
        "min_samples_split": {"type": "int", "low": 2, "high": 10}
    },
    "n_trials": 50,
    "direction": "maximize"  # Maximize prediction accuracy
}
```

## Advanced Features

### Pruning

Enable pruning to stop unpromising trials early:

```python
result = await adapter.execute(config, pruner=True)

# In your objective function, report intermediate values:
def objective(params, trial):
    for epoch in range(n_epochs):
        score = train_one_epoch(params, epoch)

        # Report intermediate value
        trial.report(score, epoch)

        # Prune if necessary
        if trial.should_prune():
            raise optuna.TrialPruned()

    return final_score
```

### Multi-Objective Optimization

For optimizing multiple objectives simultaneously, you can use weighted combinations:

```python
def multi_objective(params, trial):
    accuracy = train_model(params)
    inference_time = measure_inference_time(params)

    # Weight accuracy and speed
    return -accuracy + 0.1 * inference_time  # Minimize
```

## Integration with Other Adapters

The Optuna adapter can be used to optimize parameters for other PharmForge adapters:

```python
from adapters.vina import VinaAdapter
from adapters.optuna import OptunaAdapter

vina_adapter = VinaAdapter()
optuna_adapter = OptunaAdapter()

def optimize_docking(params, trial):
    """Optimize Vina docking parameters"""
    result = await vina_adapter.execute(
        receptor_pdbqt="protein.pdbqt",
        ligand_pdbqt="ligand.pdbqt",
        exhaustiveness=params["exhaustiveness"],
        num_modes=params["num_modes"]
    )

    if result.success:
        # Return best binding affinity
        return result.data["best_affinity"]
    else:
        return float('inf')

config = {
    "search_space": {
        "exhaustiveness": {"type": "int", "low": 8, "high": 32},
        "num_modes": {"type": "int", "low": 5, "high": 20}
    },
    "n_trials": 20
}

result = await optuna_adapter.execute(config, custom_objective=optimize_docking)
```

## Error Handling

The adapter handles various error conditions:

- **Missing Optuna**: Returns error if Optuna is not installed
- **Invalid Configuration**: Validates search space and parameters
- **Objective Failures**: Catches exceptions in objective function
- **Empty Search Space**: Validates non-empty search space

## Performance Considerations

- **Number of Trials**: Start with 50-100 trials for complex problems
- **Sampler Choice**: TPE is generally most efficient for < 100 dimensions
- **Pruning**: Use pruning for expensive objective functions
- **Parallel Execution**: Optuna supports parallel optimization (not implemented in this adapter version)

## API Reference

### OptunaAdapter

**Methods:**
- `validate_input(input_data)`: Validate optimization configuration
- `execute(input_data, **kwargs)`: Run optimization study

**Configuration:**
- `timeout`: Maximum execution time (default: 3600 seconds)
- `default_n_trials`: Default number of trials (default: 100)
- `default_sampler`: Default sampler type (default: "TPE")

## References

- [Optuna Documentation](https://optuna.readthedocs.io/)
- [Hyperparameter Optimization Best Practices](https://optuna.readthedocs.io/en/stable/tutorial/index.html)
- [TPE Algorithm Paper](https://papers.nips.cc/paper/2011/hash/86e8f7ab32cfd12577bc2619bc635690-Abstract.html)
