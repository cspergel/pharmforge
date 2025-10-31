"""
Auto-sklearn Adapter - Example Usage
Demonstrates various use cases for automated ML in drug discovery
"""

import asyncio
import numpy as np
from adapters.autosklearn import AutoSklearnAdapter


async def example_classification():
    """
    Example: Binary classification for active/inactive compounds
    """
    print("\n" + "="*70)
    print("Example 1: Binary Classification (Active/Inactive)")
    print("="*70)

    adapter = AutoSklearnAdapter()

    # Simulate molecular fingerprint data (binary features)
    # In practice, use RDKit to generate fingerprints
    np.random.seed(42)
    n_samples = 200
    n_features = 50

    # Generate synthetic data
    X_train = np.random.randint(0, 2, size=(n_samples, n_features)).tolist()
    y_train = np.random.randint(0, 2, size=n_samples).tolist()  # Active (1) or Inactive (0)

    # Split for test set
    split_idx = int(0.8 * n_samples)
    X_test = X_train[split_idx:]
    y_test = y_train[split_idx:]
    X_train = X_train[:split_idx]
    y_train = y_train[:split_idx]

    input_data = {
        "task": "classification",
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "time_limit": 120,  # 2 minutes
        "ensemble_size": 25,
        "metric": "roc_auc",
        "seed": 42
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\nResults:")
        print(f"  Training Score: {data['performance']['train_score']:.4f}")
        print(f"  Test Score: {data['performance']['test_score']:.4f}")
        print(f"  Metric: {data['performance']['metric']}")
        print(f"  Ensemble Size: {data['ensemble']['size']}")
        print(f"  Models Evaluated: {data['optimization_history']['total_models_evaluated']}")
        print(f"  Meta-Learning Used: {data['optimization_history']['meta_learning_used']}")

        if data['best_model']:
            print(f"\nBest Model:")
            print(f"  Algorithm: {data['best_model']['algorithm']}")
            print(f"  Weight in Ensemble: {data['best_model']['weight']:.4f}")

        if data['warnings']:
            print(f"\nWarnings:")
            for warning in data['warnings']:
                print(f"  - {warning}")
    else:
        print(f"Error: {result.error}")


async def example_regression():
    """
    Example: Regression for predicting IC50 values
    """
    print("\n" + "="*70)
    print("Example 2: Regression (IC50 Prediction)")
    print("="*70)

    adapter = AutoSklearnAdapter()

    # Simulate molecular descriptor data
    # In practice, compute descriptors using RDKit
    np.random.seed(42)
    n_samples = 150
    n_descriptors = 20

    # Generate synthetic data
    X_train = np.random.randn(n_samples, n_descriptors).tolist()
    y_train = (np.random.randn(n_samples) * 2 + 6).tolist()  # pIC50 values ~4-8

    # Split for test set
    split_idx = int(0.8 * n_samples)
    X_test = X_train[split_idx:]
    y_test = y_train[split_idx:]
    X_train = X_train[:split_idx]
    y_train = y_train[:split_idx]

    input_data = {
        "task": "regression",
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "time_limit": 120,  # 2 minutes
        "ensemble_size": 50,
        "metric": "r2",
        "resampling_strategy": "cv",
        "resampling_strategy_arguments": {"folds": 3},
        "seed": 42
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\nResults:")
        print(f"  Training R²: {data['performance']['train_score']:.4f}")
        print(f"  Test R²: {data['performance']['test_score']:.4f}")
        print(f"  Ensemble Size: {data['ensemble']['size']}")

        # Show top 5 models in ensemble
        print(f"\nTop 5 Models in Ensemble:")
        for i, model in enumerate(data['ensemble']['models'][:5], 1):
            print(f"  {i}. {model['algorithm']} (weight: {model['weight']:.4f})")

        # Show leaderboard
        if data['leaderboard']:
            print(f"\nLeaderboard (Top 5):")
            for entry in data['leaderboard'][:5]:
                print(f"  Rank {entry['rank']}: Score {entry['score']:.4f}")
    else:
        print(f"Error: {result.error}")


async def example_multiclass():
    """
    Example: Multi-class classification for target prediction
    """
    print("\n" + "="*70)
    print("Example 3: Multi-Class Classification (Target Prediction)")
    print("="*70)

    adapter = AutoSklearnAdapter()

    # Simulate data for predicting which target a compound binds to
    np.random.seed(42)
    n_samples = 300
    n_features = 40
    n_classes = 5  # 5 different targets

    X_train = np.random.randn(n_samples, n_features).tolist()
    y_train = np.random.randint(0, n_classes, size=n_samples).tolist()

    input_data = {
        "task": "classification",
        "X_train": X_train,
        "y_train": y_train,
        "time_limit": 150,  # 2.5 minutes
        "ensemble_size": 50,
        "metric": "accuracy",
        "seed": 42
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\nResults:")
        print(f"  Training Accuracy: {data['performance']['train_score']:.4f}")
        print(f"  Number of Classes: {n_classes}")
        print(f"  Ensemble Size: {data['ensemble']['size']}")
        print(f"  Models Evaluated: {data['optimization_history']['total_models_evaluated']}")
    else:
        print(f"Error: {result.error}")


async def example_with_preprocessing():
    """
    Example: Using specific preprocessors
    """
    print("\n" + "="*70)
    print("Example 4: Custom Preprocessing Pipeline")
    print("="*70)

    adapter = AutoSklearnAdapter()

    np.random.seed(42)
    n_samples = 200
    n_features = 30

    # Generate data with different scales (common in molecular descriptors)
    X_train = np.random.randn(n_samples, n_features) * np.random.randint(1, 100, n_features)
    y_train = np.random.randint(0, 2, size=n_samples)

    input_data = {
        "task": "classification",
        "X_train": X_train.tolist(),
        "y_train": y_train.tolist(),
        "time_limit": 120,
        "ensemble_size": 25,
        "include_preprocessors": ["StandardScaler", "PCA"],  # Force these preprocessors
        "metric": "f1",
        "seed": 42
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\nResults:")
        print(f"  Training F1: {data['performance']['train_score']:.4f}")
        print(f"  Preprocessors Requested: StandardScaler, PCA")
        print(f"  Ensemble Size: {data['ensemble']['size']}")
    else:
        print(f"Error: {result.error}")


async def example_model_export():
    """
    Example: Export trained model for later use
    """
    print("\n" + "="*70)
    print("Example 5: Model Export and Persistence")
    print("="*70)

    adapter = AutoSklearnAdapter()

    np.random.seed(42)
    n_samples = 150
    n_features = 25

    X_train = np.random.randn(n_samples, n_features).tolist()
    y_train = np.random.randint(0, 2, size=n_samples).tolist()

    import tempfile
    import os

    # Create temporary file for model export
    temp_model_path = os.path.join(tempfile.gettempdir(), "autosklearn_model.pkl")

    input_data = {
        "task": "classification",
        "X_train": X_train,
        "y_train": y_train,
        "time_limit": 90,
        "ensemble_size": 20,
        "export_model_path": temp_model_path,
        "seed": 42
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\nResults:")
        print(f"  Training Score: {data['performance']['train_score']:.4f}")

        if "model_path" in data:
            print(f"  Model Exported: {data['model_path']}")
            model_size = os.path.getsize(data['model_path']) / 1024  # KB
            print(f"  Model Size: {model_size:.2f} KB")

            # Load the model
            import pickle
            with open(data['model_path'], 'rb') as f:
                loaded_model = pickle.load(f)
            print(f"  Model Successfully Loaded: {type(loaded_model).__name__}")

            # Clean up
            os.remove(data['model_path'])
            print(f"  Temporary file cleaned up")
    else:
        print(f"Error: {result.error}")


async def example_quick_vs_thorough():
    """
    Example: Comparing quick vs thorough optimization
    """
    print("\n" + "="*70)
    print("Example 6: Quick vs Thorough Optimization")
    print("="*70)

    adapter = AutoSklearnAdapter()

    np.random.seed(42)
    n_samples = 200
    n_features = 30

    X_train = np.random.randn(n_samples, n_features).tolist()
    y_train = np.random.randn(n_samples).tolist()

    # Quick optimization (1 minute)
    print("\nQuick Optimization (60 seconds):")
    quick_input = {
        "task": "regression",
        "X_train": X_train,
        "y_train": y_train,
        "time_limit": 60,
        "per_run_time_limit": 10,
        "ensemble_size": 10,
        "seed": 42
    }

    quick_result = await adapter.execute(quick_input)

    if quick_result.success:
        data = quick_result.data
        print(f"  Training R²: {data['performance']['train_score']:.4f}")
        print(f"  Models Evaluated: {data['optimization_history']['total_models_evaluated']}")
        print(f"  Ensemble Size: {data['ensemble']['size']}")

    # Thorough optimization (3 minutes)
    print("\nThorough Optimization (180 seconds):")
    thorough_input = {
        "task": "regression",
        "X_train": X_train,
        "y_train": y_train,
        "time_limit": 180,
        "per_run_time_limit": 30,
        "ensemble_size": 50,
        "seed": 42
    }

    thorough_result = await adapter.execute(thorough_input)

    if thorough_result.success:
        data = thorough_result.data
        print(f"  Training R²: {data['performance']['train_score']:.4f}")
        print(f"  Models Evaluated: {data['optimization_history']['total_models_evaluated']}")
        print(f"  Ensemble Size: {data['ensemble']['size']}")

        # Compare
        improvement = (
            thorough_result.data['performance']['train_score'] -
            quick_result.data['performance']['train_score']
        )
        print(f"\nImprovement: {improvement:.4f} (thorough vs quick)")


async def example_rdkit_integration():
    """
    Example: Integration with RDKit for molecular descriptors
    """
    print("\n" + "="*70)
    print("Example 7: RDKit Integration (Molecular Descriptors)")
    print("="*70)

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        print("RDKit not available - skipping example")
        return

    adapter = AutoSklearnAdapter()

    # Sample SMILES strings
    smiles_list = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "C1=CC=C(C=C1)C2=CC=CC=C2",  # Biphenyl
        "CC(C)(C)NCC(C1=CC(=C(C=C1)O)CO)O",  # Salbutamol
    ] * 30  # Repeat for more samples

    # Compute molecular descriptors
    def compute_descriptors(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return [
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.NumAromaticRings(mol),
        ]

    X_train = [compute_descriptors(s) for s in smiles_list]
    X_train = [x for x in X_train if x is not None]  # Remove failed conversions

    # Simulate activity values
    np.random.seed(42)
    y_train = np.random.randn(len(X_train)) * 2 + 6  # pIC50 values

    input_data = {
        "task": "regression",
        "X_train": X_train,
        "y_train": y_train.tolist(),
        "time_limit": 120,
        "ensemble_size": 25,
        "metric": "r2",
        "seed": 42
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\nResults:")
        print(f"  Molecules: {len(X_train)}")
        print(f"  Descriptors: {len(X_train[0])}")
        print(f"  Training R²: {data['performance']['train_score']:.4f}")
        print(f"  Ensemble Size: {data['ensemble']['size']}")

        # Show which descriptors are most important (if available)
        if "feature_importance" in data:
            descriptor_names = ["MW", "LogP", "HBD", "HBA", "TPSA", "RotBonds", "ArRings"]
            importances = data['feature_importance']
            print(f"\nFeature Importance:")
            for name, importance in sorted(
                zip(descriptor_names, importances),
                key=lambda x: x[1],
                reverse=True
            ):
                print(f"  {name}: {importance:.4f}")
    else:
        print(f"Error: {result.error}")


async def main():
    """
    Run all examples
    """
    print("\n" + "="*70)
    print("Auto-sklearn Adapter - Comprehensive Examples")
    print("="*70)

    try:
        await example_classification()
        await example_regression()
        await example_multiclass()
        await example_with_preprocessing()
        await example_model_export()
        await example_quick_vs_thorough()
        await example_rdkit_integration()

        print("\n" + "="*70)
        print("All examples completed successfully!")
        print("="*70 + "\n")

    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
