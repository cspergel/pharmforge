"""
TPOT Adapter - Example Usage
Demonstrates various use cases for the TPOT adapter
"""
import asyncio
import numpy as np
from sklearn.datasets import load_breast_cancer, load_diabetes
from sklearn.model_selection import train_test_split

from adapters.tpot import TPOTAdapter


async def example_basic_classification():
    """
    Example 1: Basic binary classification
    """
    print("\n" + "="*60)
    print("Example 1: Basic Binary Classification")
    print("="*60)

    # Load breast cancer dataset
    data = load_breast_cancer()
    X_train, X_test, y_train, y_test = train_test_split(
        data.data, data.target, test_size=0.2, random_state=42
    )

    # Create adapter
    adapter = TPOTAdapter()

    # Configure TPOT for quick demo
    input_data = {
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "task": "classification",
        "generations": 3,  # Small for demo
        "population_size": 10,
        "scoring": "roc_auc",
        "cv": 3,
        "random_state": 42,
        "verbosity": 2
    }

    print(f"\nOptimizing classification pipeline...")
    print(f"Training samples: {len(X_train)}")
    print(f"Features: {X_train.shape[1]}")

    result = await adapter.execute(input_data)

    if result.success:
        print(f"\nSuccess! Best pipeline found:")
        print(f"  Training CV score: {result.data['best_score']:.4f}")
        print(f"  Test score: {result.data['test_score']:.4f}")
        print(f"\nPipeline steps:")
        for step in result.data['pipeline_steps']:
            print(f"  - {step['name']}: {step['type']}")

        print(f"\nOptimization summary:")
        print(f"  Generations evaluated: {result.data['metadata']['generations_evaluated']}")
        if result.data['optimization_history']:
            print(f"  Total pipelines tested: {result.data['metadata']['pipelines_tested']}")
            print(f"  Best generation: {result.data['metadata']['best_generation']}")

        print(f"\nExported pipeline code (first 500 chars):")
        print(result.data['pipeline_code'][:500] + "...")

    else:
        print(f"Failed: {result.error}")


async def example_regression():
    """
    Example 2: Regression task
    """
    print("\n" + "="*60)
    print("Example 2: Regression Task")
    print("="*60)

    # Load diabetes dataset
    data = load_diabetes()
    X_train, X_test, y_train, y_test = train_test_split(
        data.data, data.target, test_size=0.2, random_state=42
    )

    adapter = TPOTAdapter()

    input_data = {
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "task": "regression",
        "generations": 3,
        "population_size": 10,
        "scoring": "r2",
        "cv": 3,
        "random_state": 42,
        "verbosity": 1
    }

    print(f"\nOptimizing regression pipeline...")
    result = await adapter.execute(input_data)

    if result.success:
        print(f"\nSuccess! Best regression pipeline found:")
        print(f"  Training CV R² score: {result.data['best_score']:.4f}")
        print(f"  Test R² score: {result.data['test_score']:.4f}")
        print(f"\nBest pipeline: {result.data['best_pipeline'][:100]}...")
    else:
        print(f"Failed: {result.error}")


async def example_time_constrained():
    """
    Example 3: Time-constrained optimization
    """
    print("\n" + "="*60)
    print("Example 3: Time-Constrained Optimization")
    print("="*60)

    # Load data
    data = load_breast_cancer()
    X_train, X_test, y_train, y_test = train_test_split(
        data.data, data.target, test_size=0.2, random_state=42
    )

    adapter = TPOTAdapter()

    # Use time limits instead of fixed generations
    input_data = {
        "X_train": X_train,
        "y_train": y_train,
        "task": "classification",
        "max_time_mins": 5,  # Total runtime limit: 5 minutes
        "max_eval_time_mins": 1,  # Per-pipeline limit: 1 minute
        "population_size": 10,
        "random_state": 42,
        "verbosity": 1
    }

    print(f"\nOptimizing with 5-minute time limit...")
    result = await adapter.execute(input_data)

    if result.success:
        print(f"\nCompleted within time limit!")
        print(f"  Best score: {result.data['best_score']:.4f}")
        print(f"  Generations completed: {result.data['metadata']['generations_evaluated']}")
    else:
        print(f"Failed: {result.error}")


async def example_configuration_presets():
    """
    Example 4: Using TPOT configuration presets
    """
    print("\n" + "="*60)
    print("Example 4: Configuration Presets")
    print("="*60)

    # Load data
    data = load_breast_cancer()
    X_train, y_train = data.data[:400], data.target[:400]

    adapter = TPOTAdapter()

    # Test different presets
    presets = [None, "TPOT light"]

    for preset in presets:
        preset_name = preset if preset else "Default (full)"
        print(f"\nTesting {preset_name} configuration...")

        input_data = {
            "X_train": X_train,
            "y_train": y_train,
            "task": "classification",
            "generations": 2,
            "population_size": 5,
            "cv": 3,
            "random_state": 42,
            "verbosity": 0
        }

        if preset:
            input_data["config_dict"] = preset

        result = await adapter.execute(input_data)

        if result.success:
            print(f"  Score: {result.data['best_score']:.4f}")
            print(f"  Pipeline: {result.data['best_pipeline'][:80]}...")
        else:
            print(f"  Failed: {result.error}")


async def example_custom_scoring():
    """
    Example 5: Custom scoring metrics
    """
    print("\n" + "="*60)
    print("Example 5: Different Scoring Metrics")
    print("="*60)

    # Load data
    data = load_breast_cancer()
    X_train, y_train = data.data[:300], data.target[:300]

    adapter = TPOTAdapter()

    # Test different scoring metrics
    scoring_metrics = ["accuracy", "f1", "roc_auc", "balanced_accuracy"]

    for scoring in scoring_metrics:
        print(f"\nOptimizing for {scoring}...")

        input_data = {
            "X_train": X_train,
            "y_train": y_train,
            "task": "classification",
            "generations": 2,
            "population_size": 5,
            "scoring": scoring,
            "cv": 3,
            "random_state": 42,
            "verbosity": 0
        }

        result = await adapter.execute(input_data)

        if result.success:
            print(f"  Best {scoring}: {result.data['best_score']:.4f}")
        else:
            print(f"  Failed: {result.error}")


async def example_synthetic_data():
    """
    Example 6: Synthetic molecular data scenario
    """
    print("\n" + "="*60)
    print("Example 6: Synthetic Molecular Data")
    print("="*60)

    # Simulate molecular descriptor data
    # In real use, this would come from Mordred, RDKit, etc.
    np.random.seed(42)

    n_samples = 200
    n_features = 50  # Simulating molecular descriptors

    # Generate synthetic features
    X = np.random.randn(n_samples, n_features)

    # Generate synthetic activity labels (binary classification)
    # Simulate that activity depends on some features
    true_coefs = np.random.randn(n_features)
    y_continuous = X @ true_coefs + np.random.randn(n_samples) * 0.5
    y = (y_continuous > np.median(y_continuous)).astype(int)

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    print(f"\nSimulated molecular dataset:")
    print(f"  Compounds: {n_samples}")
    print(f"  Descriptors: {n_features}")
    print(f"  Active compounds: {np.sum(y)} ({100*np.mean(y):.1f}%)")

    adapter = TPOTAdapter()

    input_data = {
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "task": "classification",
        "generations": 5,
        "population_size": 15,
        "scoring": "roc_auc",
        "cv": 5,
        "random_state": 42,
        "verbosity": 1
    }

    print(f"\nOptimizing QSAR model...")
    result = await adapter.execute(input_data)

    if result.success:
        print(f"\nQSAR Model Results:")
        print(f"  Training ROC-AUC: {result.data['best_score']:.4f}")
        print(f"  Test ROC-AUC: {result.data['test_score']:.4f}")
        print(f"\nOptimized Pipeline:")
        for step in result.data['pipeline_steps']:
            print(f"  {step['name']}: {step['type']}")

        print(f"\nOptimization History:")
        for gen_stats in result.data['optimization_history'][:3]:
            print(f"  Gen {gen_stats['generation']}: "
                  f"best={gen_stats['best_score']:.4f}, "
                  f"mean={gen_stats['mean_score']:.4f}")

        # Show exportable code
        print(f"\nExportable Python Code:")
        print("-" * 60)
        print(result.data['pipeline_code'][:800])
        print("...")

    else:
        print(f"Failed: {result.error}")


async def example_error_handling():
    """
    Example 7: Error handling and validation
    """
    print("\n" + "="*60)
    print("Example 7: Error Handling")
    print("="*60)

    adapter = TPOTAdapter()

    # Test 1: Missing required fields
    print("\nTest 1: Missing required fields")
    result = await adapter.execute({})
    print(f"  Expected failure: {not result.success}")
    if not result.success:
        print(f"  Error: {result.error}")

    # Test 2: Invalid task type
    print("\nTest 2: Invalid task type")
    result = await adapter.execute({
        "X_train": [[1, 2], [3, 4]],
        "y_train": [0, 1],
        "task": "invalid_task"
    })
    print(f"  Expected failure: {not result.success}")

    # Test 3: Mismatched dimensions
    print("\nTest 3: Mismatched X and y dimensions")
    result = await adapter.execute({
        "X_train": [[1, 2], [3, 4]],
        "y_train": [0, 1, 2],  # Wrong length
        "task": "classification"
    })
    print(f"  Expected failure: {not result.success}")
    if not result.success:
        print(f"  Error: {result.error}")


async def main():
    """
    Run all examples
    """
    print("\n" + "="*60)
    print("TPOT Adapter - Example Usage")
    print("="*60)

    try:
        # Run examples
        await example_basic_classification()
        await example_regression()
        await example_time_constrained()
        await example_configuration_presets()
        await example_custom_scoring()
        await example_synthetic_data()
        await example_error_handling()

        print("\n" + "="*60)
        print("All examples completed!")
        print("="*60)

    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Check if TPOT is available
    from adapters.tpot import TPOT_AVAILABLE

    if not TPOT_AVAILABLE:
        print("ERROR: TPOT is not installed!")
        print("Install with: pip install tpot scikit-learn")
    else:
        print("TPOT is available. Running examples...")
        asyncio.run(main())
