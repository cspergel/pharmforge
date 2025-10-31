# ChemML Adapter Workflows

Complete workflows and examples for using the ChemML adapter in drug discovery and materials science.

## Table of Contents

1. [Property Prediction Pipeline](#property-prediction-pipeline)
2. [Multi-Representation Learning](#multi-representation-learning)
3. [Chemical Space Analysis](#chemical-space-analysis)
4. [QSAR Model Development](#qsar-model-development)
5. [Large-Scale Screening](#large-scale-screening)
6. [Integration with Other Adapters](#integration-with-other-adapters)

---

## Property Prediction Pipeline

Complete pipeline for predicting molecular properties using ChemML representations.

### Workflow: Solubility Prediction

```python
import asyncio
from adapters.chemml import ChemMLAdapter
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
import numpy as np

async def solubility_prediction_workflow():
    """
    Predict aqueous solubility using Coulomb matrices
    """
    adapter = ChemMLAdapter()

    # Training data (SMILES and solubility values)
    training_data = {
        "smiles": [
            "CCO",           # Ethanol
            "CC(=O)O",       # Acetic acid
            "c1ccccc1",      # Benzene
            "CC(C)O",        # Isopropanol
            "CCC(=O)O",      # Propionic acid
            # ... more molecules
        ],
        "solubility": [-0.77, 1.16, -1.48, -0.16, 0.88]  # log(mol/L)
    }

    # Step 1: Generate Coulomb matrices
    print("Generating Coulomb matrix representations...")
    result = await adapter({
        "smiles": training_data["smiles"],
        "operation": "coulomb_matrix",
        "parameters": {
            "max_n_atoms": 30,
            "sorting": "row_norm"
        }
    })

    if not result.success:
        print(f"Error: {result.error}")
        return

    # Step 2: Prepare features and targets
    X = np.array(result.data['representations'])
    y = np.array(training_data['solubility'])

    # Handle failed molecules
    if result.data['failed_indices']:
        print(f"Warning: {len(result.data['failed_indices'])} molecules failed")
        # Remove corresponding y values
        y = np.delete(y, result.data['failed_indices'])

    print(f"Feature shape: {X.shape}")
    print(f"Target shape: {y.shape}")

    # Step 3: Train model with cross-validation
    print("\nTraining Random Forest model...")
    model = RandomForestRegressor(
        n_estimators=100,
        max_depth=10,
        random_state=42
    )

    scores = cross_val_score(model, X, y, cv=5, scoring='r2')
    print(f"Cross-validation R² scores: {scores}")
    print(f"Mean R²: {scores.mean():.3f} (+/- {scores.std() * 2:.3f})")

    # Step 4: Train final model
    model.fit(X, y)

    # Step 5: Predict on new molecules
    new_molecules = ["CCCO", "CC(C)(C)O"]  # Propanol, tert-butanol

    print("\nPredicting solubility for new molecules...")
    new_result = await adapter({
        "smiles": new_molecules,
        "operation": "coulomb_matrix",
        "parameters": {
            "max_n_atoms": 30,
            "sorting": "row_norm"
        }
    })

    if new_result.success:
        X_new = np.array(new_result.data['representations'])
        predictions = model.predict(X_new)

        for smiles, pred in zip(new_molecules, predictions):
            print(f"{smiles}: {pred:.2f} log(mol/L)")

    return model

# Run workflow
asyncio.run(solubility_prediction_workflow())
```

---

## Multi-Representation Learning

Combine multiple molecular representations for improved predictions.

### Workflow: Ensemble Feature Engineering

```python
import asyncio
from adapters.chemml import ChemMLAdapter
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
import numpy as np

async def multi_representation_workflow():
    """
    Use multiple representations for robust predictions
    """
    adapter = ChemMLAdapter()

    smiles_list = [
        "CCO", "CC(=O)O", "c1ccccc1", "CC(C)O", "CCC(=O)O",
        "CCCO", "CC(C)CO", "c1ccc(O)cc1", "CCN", "CC(=O)N"
    ]

    target_property = [2.1, 1.5, 3.2, 2.3, 1.8, 2.0, 2.5, 3.0, 1.2, 0.9]

    # Step 1: Generate multiple representations
    representations = {}

    # Morgan Fingerprints (fast, 2D)
    print("Generating Morgan fingerprints...")
    morgan_result = await adapter({
        "smiles": smiles_list,
        "operation": "morgan_fingerprint",
        "parameters": {"radius": 2, "n_bits": 1024}
    })
    if morgan_result.success:
        representations['morgan'] = np.array(morgan_result.data['representations'])

    # RDKit Descriptors (comprehensive)
    print("Calculating RDKit descriptors...")
    rdkit_result = await adapter({
        "smiles": smiles_list,
        "operation": "rdkit_descriptors"
    })
    if rdkit_result.success:
        representations['rdkit'] = np.array(rdkit_result.data['descriptors'])

    # Coulomb Matrix (3D structure)
    print("Generating Coulomb matrices...")
    coulomb_result = await adapter({
        "smiles": smiles_list,
        "operation": "coulomb_matrix",
        "parameters": {"max_n_atoms": 25}
    })
    if coulomb_result.success:
        representations['coulomb'] = np.array(coulomb_result.data['representations'])

    # Step 2: Standardize each representation
    scaled_representations = {}
    scalers = {}

    for name, features in representations.items():
        scaler = StandardScaler()
        scaled = scaler.fit_transform(features)
        scaled_representations[name] = scaled
        scalers[name] = scaler
        print(f"{name}: {scaled.shape}")

    # Step 3: Train individual models
    print("\nTraining individual models...")
    models = {}

    for name, X in scaled_representations.items():
        model = GradientBoostingRegressor(n_estimators=50, random_state=42)
        model.fit(X, target_property)
        score = model.score(X, target_property)
        models[name] = model
        print(f"{name} R²: {score:.3f}")

    # Step 4: Ensemble predictions
    print("\nEnsemble predictions:")

    # Predict on training set
    predictions = {}
    for name, model in models.items():
        X = scaled_representations[name]
        predictions[name] = model.predict(X)

    # Average predictions (simple ensemble)
    ensemble_pred = np.mean(list(predictions.values()), axis=0)

    # Weighted ensemble (by individual performance)
    weights = {}
    for name, model in models.items():
        weights[name] = model.score(scaled_representations[name], target_property)

    total_weight = sum(weights.values())
    weighted_pred = np.zeros_like(ensemble_pred)

    for name, weight in weights.items():
        weighted_pred += predictions[name] * (weight / total_weight)

    # Calculate ensemble R²
    from sklearn.metrics import r2_score
    simple_r2 = r2_score(target_property, ensemble_pred)
    weighted_r2 = r2_score(target_property, weighted_pred)

    print(f"Simple ensemble R²: {simple_r2:.3f}")
    print(f"Weighted ensemble R²: {weighted_r2:.3f}")

    return models, scalers

# Run workflow
asyncio.run(multi_representation_workflow())
```

---

## Chemical Space Analysis

Analyze and visualize chemical space using ChemML representations.

### Workflow: Chemical Space Mapping

```python
import asyncio
from adapters.chemml import ChemMLAdapter
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np

async def chemical_space_analysis():
    """
    Map chemical space using dimensionality reduction
    """
    adapter = ChemMLAdapter()

    # Diverse molecule set
    molecules = {
        "alcohols": ["CCO", "CCCO", "CC(C)O", "CC(C)(C)O"],
        "acids": ["CC(=O)O", "CCC(=O)O", "CCCC(=O)O"],
        "aromatics": ["c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1"],
        "amines": ["CCN", "CCCN", "CC(C)N"]
    }

    all_smiles = []
    labels = []

    for category, smiles_list in molecules.items():
        all_smiles.extend(smiles_list)
        labels.extend([category] * len(smiles_list))

    # Generate Morgan fingerprints
    print("Generating molecular fingerprints...")
    result = await adapter({
        "smiles": all_smiles,
        "operation": "morgan_fingerprint",
        "parameters": {"radius": 2, "n_bits": 2048}
    })

    if not result.success:
        print(f"Error: {result.error}")
        return

    X = np.array(result.data['representations'])
    print(f"Feature matrix shape: {X.shape}")

    # PCA for initial dimensionality reduction
    print("\nPerforming PCA...")
    pca = PCA(n_components=50)
    X_pca = pca.fit_transform(X)
    print(f"Explained variance (50 components): {pca.explained_variance_ratio_.sum():.3f}")

    # t-SNE for visualization
    print("Performing t-SNE...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=5)
    X_tsne = tsne.fit_transform(X_pca)

    # Plot chemical space
    plt.figure(figsize=(10, 8))

    colors = {'alcohols': 'blue', 'acids': 'red', 'aromatics': 'green', 'amines': 'orange'}

    for category in molecules.keys():
        mask = [label == category for label in labels]
        plt.scatter(
            X_tsne[mask, 0],
            X_tsne[mask, 1],
            c=colors[category],
            label=category,
            s=100,
            alpha=0.6
        )

    plt.xlabel('t-SNE Component 1')
    plt.ylabel('t-SNE Component 2')
    plt.title('Chemical Space Visualization')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('chemical_space.png', dpi=300)
    print("\nPlot saved as 'chemical_space.png'")

    # Calculate similarity matrix
    from sklearn.metrics.pairwise import cosine_similarity

    similarity_matrix = cosine_similarity(X)

    print("\nMost similar molecule pairs:")
    for i in range(len(all_smiles)):
        for j in range(i + 1, len(all_smiles)):
            if similarity_matrix[i, j] > 0.7:  # Threshold
                print(f"{all_smiles[i]} <-> {all_smiles[j]}: {similarity_matrix[i, j]:.3f}")

# Run workflow
asyncio.run(chemical_space_analysis())
```

---

## QSAR Model Development

Complete QSAR (Quantitative Structure-Activity Relationship) model development.

### Workflow: Bioactivity Prediction

```python
import asyncio
from adapters.chemml import ChemMLAdapter
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import numpy as np
import pandas as pd

async def qsar_development_workflow():
    """
    Develop QSAR model for bioactivity prediction
    """
    adapter = ChemMLAdapter()

    # Example dataset: IC50 values (nM)
    dataset = pd.DataFrame({
        'smiles': [
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
            "CC(=O)Oc1ccccc1C(=O)O",        # Aspirin
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine
            "CC(C)NCC(COc1ccccc1)O",        # Propranolol
            # Add more molecules...
        ],
        'ic50': [10.2, 25.3, 150.5, 5.8],  # nM
        'activity_class': ['active', 'active', 'inactive', 'active']
    })

    # Convert IC50 to pIC50 (log scale)
    dataset['pic50'] = -np.log10(dataset['ic50'] * 1e-9)

    print(f"Dataset size: {len(dataset)} molecules")
    print(f"pIC50 range: {dataset['pic50'].min():.2f} - {dataset['pic50'].max():.2f}")

    # Step 1: Generate features
    print("\nGenerating molecular representations...")

    # Try multiple representations
    feature_sets = {}

    # Coulomb matrices
    coulomb_result = await adapter({
        "smiles": dataset['smiles'].tolist(),
        "operation": "coulomb_matrix",
        "parameters": {"max_n_atoms": 40, "sorting": "row_norm"}
    })

    if coulomb_result.success:
        feature_sets['coulomb'] = np.array(coulomb_result.data['representations'])
        print(f"Coulomb matrices: {feature_sets['coulomb'].shape}")

    # Morgan fingerprints
    morgan_result = await adapter({
        "smiles": dataset['smiles'].tolist(),
        "operation": "morgan_fingerprint",
        "parameters": {"radius": 3, "n_bits": 2048}
    })

    if morgan_result.success:
        feature_sets['morgan'] = np.array(morgan_result.data['representations'])
        print(f"Morgan fingerprints: {feature_sets['morgan'].shape}")

    # RDKit descriptors
    rdkit_result = await adapter({
        "smiles": dataset['smiles'].tolist(),
        "operation": "rdkit_descriptors"
    })

    if rdkit_result.success:
        feature_sets['rdkit'] = np.array(rdkit_result.data['descriptors'])
        print(f"RDKit descriptors: {feature_sets['rdkit'].shape}")

    # Step 2: Model development for each feature set
    y = dataset['pic50'].values

    results = {}

    for feature_name, X in feature_sets.items():
        print(f"\n{'='*60}")
        print(f"Training model with {feature_name} features")
        print(f"{'='*60}")

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )

        # Hyperparameter tuning
        param_grid = {
            'n_estimators': [50, 100, 200],
            'max_depth': [5, 10, 15],
            'min_samples_split': [2, 5, 10]
        }

        rf = RandomForestRegressor(random_state=42)

        grid_search = GridSearchCV(
            rf, param_grid, cv=3, scoring='neg_mean_squared_error', n_jobs=-1
        )

        grid_search.fit(X_train, y_train)

        # Best model
        best_model = grid_search.best_estimator_

        # Predictions
        y_train_pred = best_model.predict(X_train)
        y_test_pred = best_model.predict(X_test)

        # Metrics
        train_r2 = r2_score(y_train, y_train_pred)
        test_r2 = r2_score(y_test, y_test_pred)
        train_rmse = np.sqrt(mean_squared_error(y_train, y_train_pred))
        test_rmse = np.sqrt(mean_squared_error(y_test, y_test_pred))
        test_mae = mean_absolute_error(y_test, y_test_pred)

        results[feature_name] = {
            'model': best_model,
            'train_r2': train_r2,
            'test_r2': test_r2,
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'test_mae': test_mae,
            'best_params': grid_search.best_params_
        }

        print(f"\nBest parameters: {grid_search.best_params_}")
        print(f"Training R²: {train_r2:.3f}")
        print(f"Test R²: {test_r2:.3f}")
        print(f"Test RMSE: {test_rmse:.3f}")
        print(f"Test MAE: {test_mae:.3f}")

    # Step 3: Compare models
    print(f"\n{'='*60}")
    print("Model Comparison")
    print(f"{'='*60}")

    comparison_df = pd.DataFrame([
        {
            'Feature Set': name,
            'Test R²': res['test_r2'],
            'Test RMSE': res['test_rmse'],
            'Test MAE': res['test_mae']
        }
        for name, res in results.items()
    ])

    print(comparison_df.to_string(index=False))

    # Select best model
    best_feature_set = max(results.items(), key=lambda x: x[1]['test_r2'])
    print(f"\nBest feature set: {best_feature_set[0]}")

    return results

# Run workflow
asyncio.run(qsar_development_workflow())
```

---

## Large-Scale Screening

Efficient processing of large molecular libraries.

### Workflow: Virtual Screening

```python
import asyncio
from adapters.chemml import ChemMLAdapter
import numpy as np
from typing import List
import pickle

async def virtual_screening_workflow(
    library_smiles: List[str],
    trained_model,
    scaler,
    threshold: float = 7.0
):
    """
    Screen large molecular library for active compounds

    Args:
        library_smiles: Large list of SMILES strings
        trained_model: Pre-trained ML model
        scaler: Feature scaler
        threshold: Activity threshold (e.g., pIC50 > 7.0)
    """
    adapter = ChemMLAdapter()

    batch_size = 100  # Process in batches
    hits = []

    print(f"Screening {len(library_smiles)} molecules...")
    print(f"Batch size: {batch_size}")

    for i in range(0, len(library_smiles), batch_size):
        batch = library_smiles[i:i + batch_size]
        print(f"\nProcessing batch {i//batch_size + 1}/{(len(library_smiles)-1)//batch_size + 1}")

        # Generate features
        result = await adapter({
            "smiles": batch,
            "operation": "morgan_fingerprint",
            "parameters": {"radius": 2, "n_bits": 2048}
        })

        if not result.success:
            print(f"Warning: Batch failed - {result.error}")
            continue

        X = np.array(result.data['representations'])

        # Scale features
        X_scaled = scaler.transform(X)

        # Predict
        predictions = trained_model.predict(X_scaled)

        # Filter hits
        for j, (smiles, pred) in enumerate(zip(batch, predictions)):
            if pred >= threshold:
                hits.append({
                    'smiles': smiles,
                    'predicted_pic50': pred,
                    'predicted_ic50_nm': 10**(-pred) * 1e9
                })

        print(f"Found {len(hits)} hits so far...")

    # Sort by predicted activity
    hits.sort(key=lambda x: x['predicted_pic50'], reverse=True)

    print(f"\n{'='*60}")
    print(f"Screening complete!")
    print(f"Total hits (pIC50 >= {threshold}): {len(hits)}")
    print(f"Hit rate: {len(hits)/len(library_smiles)*100:.2f}%")
    print(f"{'='*60}")

    # Top 10 hits
    print("\nTop 10 predicted actives:")
    for i, hit in enumerate(hits[:10], 1):
        print(f"{i}. {hit['smiles']}")
        print(f"   Predicted pIC50: {hit['predicted_pic50']:.2f}")
        print(f"   Predicted IC50: {hit['predicted_ic50_nm']:.2f} nM")

    return hits

# Example usage
async def run_screening():
    # Load pre-trained model
    with open('qsar_model.pkl', 'rb') as f:
        model_data = pickle.load(f)

    model = model_data['model']
    scaler = model_data['scaler']

    # Load screening library (could be from ZINC, ChEMBL, etc.)
    # For demo, generate random SMILES
    library = [
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "CCOc1ccc(cc1)C(=O)c2ccccc2",
        # ... thousands more
    ]

    hits = await virtual_screening_workflow(
        library,
        model,
        scaler,
        threshold=7.0
    )

    return hits

# Run
# asyncio.run(run_screening())
```

---

## Integration with Other Adapters

Combine ChemML with other PharmForge adapters for complete workflows.

### Workflow: Discovery Pipeline

```python
import asyncio
from adapters.chemml import ChemMLAdapter
from adapters.rdkit_local import RDKitAdapter
from adapters.vina import VinaAdapter

async def integrated_discovery_pipeline():
    """
    Complete pipeline: Representation -> Prediction -> Validation
    """
    chemml = ChemMLAdapter()
    rdkit = RDKitAdapter()
    vina = VinaAdapter()

    candidate_molecules = [
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "CCOc1ccc(cc1)C(=O)c2ccccc2",
        "c1ccc2c(c1)ccc3c2cccc3"
    ]

    results = []

    for smiles in candidate_molecules:
        print(f"\nProcessing: {smiles}")

        # Step 1: Generate ChemML features
        chemml_result = await chemml({
            "smiles": smiles,
            "operation": "coulomb_matrix"
        })

        # Step 2: Calculate basic properties
        rdkit_result = await rdkit(smiles)

        # Step 3: If promising, do docking
        if rdkit_result.success:
            props = rdkit_result.data

            # Filter by Lipinski's rules
            if props['lipinski_violations'] <= 1:
                print("  Passes Lipinski's rules - proceeding to docking")

                # Docking would go here
                # docking_result = await vina.dock(...)

                results.append({
                    'smiles': smiles,
                    'properties': props,
                    'chemml_features': chemml_result.data if chemml_result.success else None,
                    'passed_filters': True
                })
            else:
                print(f"  Failed Lipinski's rules ({props['lipinski_violations']} violations)")
                results.append({
                    'smiles': smiles,
                    'properties': props,
                    'passed_filters': False
                })

    return results

# Run
asyncio.run(integrated_discovery_pipeline())
```

---

## Best Practices

### 1. Feature Selection

- Start with Morgan fingerprints for quick prototyping
- Use Coulomb matrices when 3D effects are important
- Combine multiple representations for robustness

### 2. Performance Optimization

- Process molecules in batches (50-100 at a time)
- Cache results for repeated analyses
- Use simpler representations for large-scale screening

### 3. Model Validation

- Always use train/test split or cross-validation
- Check for overfitting (compare train vs test performance)
- Validate on external test sets when possible

### 4. Error Handling

- Always check `result.success` before using data
- Monitor `failed_indices` for problematic molecules
- Have fallback representations for 3D generation failures

---

## Next Steps

1. Explore combination with other adapters (DeepChem, Mordred)
2. Implement custom ML pipelines
3. Scale to larger datasets
4. Integrate with experimental data

For more information, see the main README.md file.
