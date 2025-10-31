# TPOT Integration Guide

## Overview

This guide demonstrates how to integrate TPOT with other PharmForge adapters to create powerful automated machine learning workflows.

## Integration Patterns

### 1. Feature Generation + TPOT

Use molecular feature generators with TPOT for automated drug discovery pipelines.

#### With Mordred Descriptors

```python
from adapters.mordred import MordredAdapter
from adapters.tpot import TPOTAdapter

async def optimize_qsar_model(smiles_list, activity_labels):
    """
    Generate molecular descriptors and optimize QSAR model
    """
    # Step 1: Generate molecular descriptors
    mordred = MordredAdapter()
    descriptor_result = await mordred.execute({
        "smiles": smiles_list,
        "descriptors": ["all"],
        "standardize": True
    })

    if not descriptor_result.success:
        return None

    X = descriptor_result.data["descriptors"]

    # Step 2: Optimize ML pipeline with TPOT
    tpot = TPOTAdapter()
    tpot_result = await tpot.execute({
        "X_train": X,
        "y_train": activity_labels,
        "task": "classification",
        "generations": 5,
        "population_size": 20,
        "scoring": "roc_auc",
        "cv": 5
    })

    return tpot_result
```

#### With Scikit-mol Fingerprints

```python
from adapters.scikit_mol import ScikitMolAdapter
from adapters.tpot import TPOTAdapter

async def optimize_with_fingerprints(smiles_list, labels):
    """
    Use molecular fingerprints with TPOT
    """
    # Generate Morgan fingerprints
    scikitmol = ScikitMolAdapter()
    fp_result = await scikitmol.execute({
        "smiles": smiles_list,
        "operation": "fingerprints",
        "fp_type": "morgan",
        "radius": 2,
        "n_bits": 2048
    })

    X = fp_result.data["fingerprints"]

    # Optimize with TPOT
    tpot = TPOTAdapter()
    result = await tpot.execute({
        "X_train": X,
        "y_train": labels,
        "task": "classification",
        "generations": 10,
        "population_size": 50
    })

    return result
```

#### With Molfeat Feature Engineering

```python
from adapters.molfeat import MolfeatAdapter
from adapters.tpot import TPOTAdapter

async def optimize_with_molfeat(smiles_list, labels):
    """
    Use Molfeat's advanced features with TPOT
    """
    # Generate features using Molfeat
    molfeat = MolfeatAdapter()
    feat_result = await molfeat.execute({
        "smiles": smiles_list,
        "featurizer": "ecfp",
        "params": {"radius": 3, "length": 2048}
    })

    X = feat_result.data["features"]

    # Optimize with TPOT
    tpot = TPOTAdapter()
    result = await tpot.execute({
        "X_train": X,
        "y_train": labels,
        "task": "regression",
        "generations": 5,
        "scoring": "r2"
    })

    return result
```

### 2. TPOT + Optuna (Meta-Optimization)

Use Optuna to optimize TPOT's hyperparameters for maximum performance.

```python
from adapters.optuna import OptunaAdapter
from adapters.tpot import TPOTAdapter

async def meta_optimize(X_train, y_train, X_val, y_val):
    """
    Use Optuna to find best TPOT configuration
    """
    tpot = TPOTAdapter()

    # Define search space for TPOT parameters
    search_space = {
        "generations": {
            "type": "int",
            "low": 3,
            "high": 10
        },
        "population_size": {
            "type": "int",
            "low": 10,
            "high": 50
        },
        "cv": {
            "type": "int",
            "low": 3,
            "high": 10
        },
        "max_eval_time_mins": {
            "type": "int",
            "low": 1,
            "high": 10
        }
    }

    # Objective function
    async def objective(params, trial):
        result = await tpot.execute({
            "X_train": X_train,
            "y_train": y_train,
            "X_test": X_val,
            "y_test": y_val,
            "task": "classification",
            "generations": params["generations"],
            "population_size": params["population_size"],
            "cv": params["cv"],
            "max_eval_time_mins": params["max_eval_time_mins"],
            "verbosity": 0  # Quiet mode
        })

        if result.success:
            return result.data["test_score"]
        return 0.0

    # Run meta-optimization
    optuna = OptunaAdapter()
    optuna_result = await optuna.execute(
        {
            "search_space": search_space,
            "objective_function": objective,
            "n_trials": 20,
            "direction": "maximize"
        },
        custom_objective=objective
    )

    # Get best TPOT configuration
    best_params = optuna_result.data["best_params"]

    # Run TPOT with best parameters on full dataset
    final_result = await tpot.execute({
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_val,
        "y_test": y_val,
        "task": "classification",
        **best_params
    })

    return final_result, optuna_result
```

### 3. Database Query + Feature Generation + TPOT

Complete end-to-end pipeline from database to optimized model.

```python
from adapters.chembl import ChEMBLAdapter
from adapters.mordred import MordredAdapter
from adapters.tpot import TPOTAdapter

async def build_target_model(target_id, activity_threshold=6.5):
    """
    Build optimized QSAR model for ChEMBL target
    """
    # Step 1: Query ChEMBL for bioactivity data
    chembl = ChEMBLAdapter()
    activity_result = await chembl.execute({
        "operation": "get_target_activities",
        "target_id": target_id,
        "min_activities": 100
    })

    if not activity_result.success:
        return None

    activities = activity_result.data["activities"]

    # Extract SMILES and create labels
    smiles_list = [a["molecule_smiles"] for a in activities]
    pchembl_values = [a.get("pchembl_value", 0) for a in activities]
    labels = [1 if v >= activity_threshold else 0 for v in pchembl_values]

    # Step 2: Generate molecular descriptors
    mordred = MordredAdapter()
    descriptor_result = await mordred.execute({
        "smiles": smiles_list,
        "descriptors": ["all"],
        "standardize": True
    })

    X = descriptor_result.data["descriptors"]

    # Step 3: Split data
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X, labels, test_size=0.2, random_state=42
    )

    # Step 4: Optimize with TPOT
    tpot = TPOTAdapter()
    model_result = await tpot.execute({
        "X_train": X_train,
        "y_train": y_train,
        "X_test": X_test,
        "y_test": y_test,
        "task": "classification",
        "generations": 10,
        "population_size": 50,
        "scoring": "roc_auc",
        "cv": 5
    })

    return {
        "target_id": target_id,
        "n_compounds": len(smiles_list),
        "model_result": model_result,
        "activity_data": activities
    }
```

### 4. Virtual Screening with TPOT Models

Use TPOT models for virtual screening of compound libraries.

```python
from adapters.tpot import TPOTAdapter
from adapters.mordred import MordredAdapter
from adapters.zinc import ZINCAdapter

async def virtual_screening(training_smiles, training_labels, n_candidates=1000):
    """
    Train model and screen ZINC compounds
    """
    # Step 1: Generate descriptors for training set
    mordred = MordredAdapter()
    train_desc = await mordred.execute({
        "smiles": training_smiles,
        "descriptors": ["all"]
    })

    X_train = train_desc.data["descriptors"]

    # Step 2: Train TPOT model
    tpot = TPOTAdapter()
    model_result = await tpot.execute({
        "X_train": X_train,
        "y_train": training_labels,
        "task": "classification",
        "generations": 10,
        "population_size": 50
    })

    if not model_result.success:
        return None

    fitted_pipeline = model_result.data["fitted_pipeline"]

    # Step 3: Get ZINC compounds
    zinc = ZINCAdapter()
    zinc_result = await zinc.execute({
        "operation": "search",
        "query_type": "similarity",
        "query_smiles": training_smiles[0],  # Use first active compound
        "limit": n_candidates
    })

    screening_smiles = [c["smiles"] for c in zinc_result.data["compounds"]]

    # Step 4: Generate descriptors for screening set
    screen_desc = await mordred.execute({
        "smiles": screening_smiles,
        "descriptors": ["all"]
    })

    X_screen = screen_desc.data["descriptors"]

    # Step 5: Predict activities
    predictions = fitted_pipeline.predict_proba(X_screen)[:, 1]

    # Step 6: Rank compounds
    results = []
    for i, (smiles, score) in enumerate(zip(screening_smiles, predictions)):
        results.append({
            "smiles": smiles,
            "predicted_score": float(score),
            "rank": i + 1
        })

    # Sort by predicted score
    results.sort(key=lambda x: x["predicted_score"], reverse=True)

    return {
        "model": model_result,
        "candidates": results,
        "n_screened": len(results)
    }
```

### 5. Multi-Task Learning Setup

Optimize models for multiple related tasks.

```python
from adapters.tpot import TPOTAdapter
from adapters.mordred import MordredAdapter

async def multitask_optimization(smiles_list, task_labels_dict):
    """
    Optimize separate models for multiple tasks
    task_labels_dict: {"task1": [labels], "task2": [labels], ...}
    """
    # Generate features once
    mordred = MordredAdapter()
    desc_result = await mordred.execute({
        "smiles": smiles_list,
        "descriptors": ["all"]
    })

    X = desc_result.data["descriptors"]

    # Optimize model for each task
    results = {}
    tpot = TPOTAdapter()

    for task_name, labels in task_labels_dict.items():
        print(f"Optimizing model for {task_name}...")

        result = await tpot.execute({
            "X_train": X,
            "y_train": labels,
            "task": "classification",
            "generations": 5,
            "population_size": 20,
            "scoring": "roc_auc"
        })

        results[task_name] = result

    return results
```

### 6. TPOT + Visualization

Combine TPOT with visualization adapters for model interpretation.

```python
from adapters.tpot import TPOTAdapter
from adapters.chemplot import ChemplotAdapter

async def optimize_and_visualize(smiles_list, labels):
    """
    Optimize model and visualize chemical space
    """
    # Optimize with TPOT
    tpot = TPOTAdapter()
    result = await tpot.execute({
        "X_train": X_train,
        "y_train": labels,
        "task": "classification"
    })

    # Visualize chemical space
    chemplot = ChemplotAdapter()
    viz_result = await chemplot.execute({
        "smiles": smiles_list,
        "target": labels,
        "viz_type": "pca",
        "title": "TPOT Optimized Model - Chemical Space"
    })

    return {
        "model": result,
        "visualization": viz_result
    }
```

### 7. ADMET Prediction Pipeline

Combine ADMET prediction with TPOT.

```python
from adapters.admet_ai import ADMETaiAdapter
from adapters.tpot import TPOTAdapter

async def admet_guided_optimization(smiles_list, activity_labels):
    """
    Include ADMET properties as features for model optimization
    """
    # Step 1: Predict ADMET properties
    admet = ADMETaiAdapter()
    admet_result = await admet.execute({
        "smiles": smiles_list,
        "properties": ["Caco2", "Solubility", "CYP3A4_Substrate"]
    })

    admet_features = admet_result.data["predictions"]

    # Step 2: Get molecular descriptors
    mordred = MordredAdapter()
    desc_result = await mordred.execute({
        "smiles": smiles_list,
        "descriptors": ["basic"]
    })

    # Step 3: Combine features
    import numpy as np
    X_descriptors = np.array(desc_result.data["descriptors"])
    X_admet = np.array(admet_features)
    X_combined = np.hstack([X_descriptors, X_admet])

    # Step 4: Optimize with TPOT
    tpot = TPOTAdapter()
    result = await tpot.execute({
        "X_train": X_combined,
        "y_train": activity_labels,
        "task": "classification",
        "generations": 10,
        "population_size": 50
    })

    return result
```

## Best Integration Practices

### 1. Cache Intermediate Results

```python
# Generate features once and reuse
descriptor_result = await mordred.execute({...})
X = descriptor_result.data["descriptors"]

# Try different TPOT configurations with same features
for config in configurations:
    result = await tpot.execute({
        "X_train": X,
        "y_train": y,
        **config
    })
```

### 2. Pipeline Error Handling

```python
async def robust_pipeline(smiles_list, labels):
    """
    Pipeline with comprehensive error handling
    """
    try:
        # Feature generation
        mordred = MordredAdapter()
        desc_result = await mordred.execute({
            "smiles": smiles_list,
            "descriptors": ["all"]
        })

        if not desc_result.success:
            print(f"Feature generation failed: {desc_result.error}")
            return None

        # Model optimization
        tpot = TPOTAdapter()
        model_result = await tpot.execute({
            "X_train": desc_result.data["descriptors"],
            "y_train": labels,
            "task": "classification"
        })

        if not model_result.success:
            print(f"Model optimization failed: {model_result.error}")
            return None

        return model_result

    except Exception as e:
        print(f"Pipeline failed: {e}")
        return None
```

### 3. Resource Management

```python
# Use time limits to prevent excessive computation
config = {
    "max_time_mins": 30,
    "max_eval_time_mins": 2,
    "n_jobs": -1,  # Use all cores
    "verbosity": 1  # Reduce logging
}
```

### 4. Progressive Optimization

```python
async def progressive_optimization(X, y):
    """
    Start fast, then refine
    """
    tpot = TPOTAdapter()

    # Quick exploration
    quick_result = await tpot.execute({
        "X_train": X,
        "y_train": y,
        "task": "classification",
        "generations": 3,
        "population_size": 10,
        "config_dict": "TPOT light"
    })

    # Detailed optimization
    detailed_result = await tpot.execute({
        "X_train": X,
        "y_train": y,
        "task": "classification",
        "generations": 20,
        "population_size": 100,
        "max_time_mins": 120
    })

    return quick_result, detailed_result
```

## Example Workflows

### Complete Drug Discovery Pipeline

```python
async def full_pipeline(target_name, activity_threshold=6.5):
    """
    End-to-end drug discovery pipeline
    """
    # 1. Query bioactivity data
    chembl = ChEMBLAdapter()
    # ... query ChEMBL

    # 2. Generate molecular features
    mordred = MordredAdapter()
    # ... generate descriptors

    # 3. Optimize ML model
    tpot = TPOTAdapter()
    # ... optimize pipeline

    # 4. Virtual screening
    zinc = ZINCAdapter()
    # ... screen compounds

    # 5. ADMET filtering
    admet = ADMETaiAdapter()
    # ... filter by ADMET

    return final_candidates
```

This comprehensive integration guide shows how TPOT can be combined with other PharmForge adapters to create sophisticated automated machine learning workflows for drug discovery.
