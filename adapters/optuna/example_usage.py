"""
Example usage of the Optuna Adapter for PharmForge

This script demonstrates how to use the Optuna adapter for hyperparameter optimization
"""
import asyncio
import logging
from adapter import OptunaAdapter

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


async def example_simple_optimization():
    """
    Simple example: Optimize a mathematical function
    Find parameters that minimize: (x - 2)^2 + (y - 3)^2
    """
    print("\n" + "="*60)
    print("Example 1: Simple Mathematical Optimization")
    print("="*60)

    adapter = OptunaAdapter()

    # Define a simple objective function
    def simple_objective(params, trial):
        """Minimize: (x - 2)^2 + (y - 3)^2"""
        x = params["x"]
        y = params["y"]
        result = (x - 2)**2 + (y - 3)**2
        logger.info(f"Trial {trial.number}: x={x:.4f}, y={y:.4f}, value={result:.4f}")
        return result

    # Configuration
    config = {
        "search_space": {
            "x": {
                "type": "float",
                "low": -10.0,
                "high": 10.0
            },
            "y": {
                "type": "float",
                "low": -10.0,
                "high": 10.0
            }
        },
        "objective_function": "custom",
        "n_trials": 20,
        "sampler": "TPE",
        "direction": "minimize",
        "study_name": "simple_math_optimization"
    }

    # Run optimization
    result = await adapter.execute(config, custom_objective=simple_objective)

    if result.success:
        print("\n✓ Optimization completed successfully!")
        print(f"  Best parameters: x={result.data['best_params']['x']:.4f}, "
              f"y={result.data['best_params']['y']:.4f}")
        print(f"  Best value: {result.data['best_value']:.6f}")
        print(f"  Best trial: {result.data['best_trial']}")
        print(f"  Expected: x≈2.0, y≈3.0, value≈0.0")
    else:
        print(f"\n✗ Optimization failed: {result.error}")


async def example_ml_hyperparameters():
    """
    Example: Optimize machine learning hyperparameters
    Simulates optimizing a neural network configuration
    """
    print("\n" + "="*60)
    print("Example 2: Machine Learning Hyperparameter Optimization")
    print("="*60)

    adapter = OptunaAdapter()

    def ml_objective(params, trial):
        """
        Simulate ML model training with different hyperparameters
        Returns a simulated validation loss
        """
        lr = params["learning_rate"]
        batch_size = params["batch_size"]
        hidden_dim = params["hidden_dim"]
        dropout = params["dropout"]

        # Simulate training (in reality, you would train an actual model)
        # This is a dummy formula to demonstrate the concept
        simulated_loss = (
            0.5 / (lr * 100) +  # Lower LR → higher loss
            abs(batch_size - 64) / 100 +  # Optimal batch size around 64
            abs(hidden_dim - 256) / 1000 +  # Optimal hidden dim around 256
            dropout * 0.5  # Penalty for high dropout
        )

        logger.info(f"Trial {trial.number}: lr={lr:.5f}, batch={batch_size}, "
                   f"hidden={hidden_dim}, dropout={dropout:.3f}, loss={simulated_loss:.4f}")

        return simulated_loss

    config = {
        "search_space": {
            "learning_rate": {
                "type": "float",
                "low": 1e-5,
                "high": 1e-2,
                "log": True
            },
            "batch_size": {
                "type": "int",
                "low": 16,
                "high": 128
            },
            "hidden_dim": {
                "type": "int",
                "low": 64,
                "high": 512
            },
            "dropout": {
                "type": "float",
                "low": 0.0,
                "high": 0.5
            }
        },
        "objective_function": "custom",
        "n_trials": 30,
        "sampler": "TPE",
        "direction": "minimize",
        "study_name": "ml_hyperparameter_optimization"
    }

    result = await adapter.execute(config, custom_objective=ml_objective)

    if result.success:
        print("\n✓ Optimization completed successfully!")
        print(f"  Best parameters:")
        for param, value in result.data['best_params'].items():
            print(f"    {param}: {value}")
        print(f"  Best validation loss: {result.data['best_value']:.6f}")
        print(f"  Total trials: {result.data['n_trials']}")
    else:
        print(f"\n✗ Optimization failed: {result.error}")


async def example_categorical_optimization():
    """
    Example: Optimize with categorical parameters
    Choose best algorithm and its parameters
    """
    print("\n" + "="*60)
    print("Example 3: Optimization with Categorical Parameters")
    print("="*60)

    adapter = OptunaAdapter()

    def algorithm_objective(params, trial):
        """
        Simulate comparing different algorithms
        """
        algorithm = params["algorithm"]
        param_value = params["algorithm_param"]

        # Simulate different algorithm performances
        if algorithm == "algorithm_a":
            score = (param_value - 5)**2 + 1.0
        elif algorithm == "algorithm_b":
            score = abs(param_value - 10) + 2.0
        else:  # algorithm_c
            score = (param_value - 7.5)**2 + 0.5

        logger.info(f"Trial {trial.number}: {algorithm} with param={param_value:.2f}, score={score:.4f}")
        return score

    config = {
        "search_space": {
            "algorithm": {
                "type": "categorical",
                "choices": ["algorithm_a", "algorithm_b", "algorithm_c"]
            },
            "algorithm_param": {
                "type": "float",
                "low": 0.0,
                "high": 15.0
            }
        },
        "objective_function": "custom",
        "n_trials": 25,
        "sampler": "TPE",
        "direction": "minimize",
        "study_name": "algorithm_selection"
    }

    result = await adapter.execute(config, custom_objective=algorithm_objective)

    if result.success:
        print("\n✓ Optimization completed successfully!")
        print(f"  Best algorithm: {result.data['best_params']['algorithm']}")
        print(f"  Best parameter: {result.data['best_params']['algorithm_param']:.4f}")
        print(f"  Best score: {result.data['best_value']:.6f}")
    else:
        print(f"\n✗ Optimization failed: {result.error}")


async def example_with_history():
    """
    Example: Access trial history for visualization
    """
    print("\n" + "="*60)
    print("Example 4: Accessing Trial History")
    print("="*60)

    adapter = OptunaAdapter()

    def sphere_function(params, trial):
        """N-dimensional sphere function"""
        return sum(params[f"x{i}"]**2 for i in range(3))

    config = {
        "search_space": {
            "x0": {"type": "float", "low": -5.0, "high": 5.0},
            "x1": {"type": "float", "low": -5.0, "high": 5.0},
            "x2": {"type": "float", "low": -5.0, "high": 5.0},
        },
        "objective_function": "custom",
        "n_trials": 15,
        "sampler": "Random",
        "direction": "minimize",
        "study_name": "sphere_optimization"
    }

    result = await adapter.execute(config, custom_objective=sphere_function)

    if result.success:
        print("\n✓ Optimization completed!")
        print(f"  Best value: {result.data['best_value']:.6f}")

        # Show trial history
        print("\n  Trial History (first 5 and best):")
        history = result.data['trial_history']
        for i, trial in enumerate(history[:5]):
            print(f"    Trial {trial['number']}: value={trial['value']:.6f}, "
                  f"duration={trial['duration']:.2f}s")

        best_trial_num = result.data['best_trial']
        best_trial = history[best_trial_num]
        print(f"\n  Best Trial {best_trial_num}:")
        print(f"    Value: {best_trial['value']:.6f}")
        print(f"    Parameters: {best_trial['params']}")
    else:
        print(f"\n✗ Optimization failed: {result.error}")


async def main():
    """Run all examples"""
    print("\n" + "="*60)
    print("Optuna Adapter Examples for PharmForge")
    print("="*60)

    try:
        # Run examples
        await example_simple_optimization()
        await example_ml_hyperparameters()
        await example_categorical_optimization()
        await example_with_history()

        print("\n" + "="*60)
        print("All examples completed!")
        print("="*60)

    except ImportError as e:
        print(f"\n✗ Error: {e}")
        print("Please install Optuna: pip install optuna")
    except Exception as e:
        logger.error(f"Unexpected error: {e}", exc_info=True)


if __name__ == "__main__":
    asyncio.run(main())
