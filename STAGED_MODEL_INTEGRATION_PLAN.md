# Staged Model Integration Plan
## PharmForge ADMET System

**Philosophy**: Validate before integrating. Quality over quantity.

---

## Stage 1: ADMET-AI Foundation (Week 3 - Current) ✅

### Goal
Establish a **working, reliable baseline** with ADMET-AI

### Deliverables
- ✅ ADMET-AI installed and tested (49 properties)
- 🔄 `ADMETaiAdapter` implementation (clean, simple)
- 🔄 Unit tests for ADMET-AI adapter
- 🔄 Integration with adapter registry

### Why ADMET-AI First?
- **Production-ready**: Well-tested, actively maintained
- **Comprehensive**: 49 ADMET properties out of the box
- **Fast**: Optimized Chemprop models
- **Validated**: Published benchmark results
- **No training needed**: Pre-trained models

### Success Criteria
- ✅ Can predict 49 ADMET properties for any valid SMILES
- ✅ Predictions complete in <5 seconds per molecule
- ✅ Adapter follows `AdapterProtocol`
- ✅ Full test coverage

---

## Stage 2: Custom Model Evaluation (Week 3.5-4)

### Goal
**Rigorously evaluate** your custom models before integration

### Process

#### Step 1: Model Inventory
```bash
# Scan for model files
find adapters/custom -name "*.pkl" -o -name "*.pt"
```

Expected inventory:
- **Absorption**: Caco-2 (GCN, RF, XGBoost), HIA (RF), solubility (RF)
- **Distribution**: BBB (2x RF), P-gp (GNN, RF), protein binding (2x RF)
- **Metabolism**: CYP models
- **Excretion**: Hepatic clearance (2x models)
- **Toxicity**: Various toxicity endpoints

#### Step 2: Document Each Model
For each model file, create a metadata file:

**Example**: `adapters/custom/caco2_gcn_metadata.json`
```json
{
    "model_file": "caco2_gcn.pt",
    "property": "caco2_permeability",
    "model_type": "Graph Convolutional Network",
    "training_dataset": "TDC Caco-2 Wang",
    "training_date": "2024-XX-XX",
    "framework": "PyTorch",
    "input_format": "SMILES",
    "output_format": "log_permeability",
    "normalization": "0-1, higher=better",
    "performance": {
        "test_r2": null,
        "test_mae": null,
        "test_rmse": null,
        "dataset_size": null
    },
    "notes": "Needs evaluation"
}
```

#### Step 3: Create Evaluation Framework
```python
# adapters/custom/evaluate_models.py

class ModelEvaluator:
    """Evaluates custom models against ADMET-AI and benchmarks"""

    def evaluate_model(self, model_path, test_smiles, reference_values):
        """
        Args:
            model_path: Path to .pkl or .pt file
            test_smiles: List of SMILES for testing
            reference_values: Ground truth or ADMET-AI predictions

        Returns:
            metrics: {
                'r2': correlation with reference,
                'mae': mean absolute error,
                'agreement': % predictions within 20% of reference
            }
        """
        # Load model
        model = self.load_model(model_path)

        # Predict on test set
        predictions = model.predict(test_smiles)

        # Compare to reference
        metrics = self.calculate_metrics(predictions, reference_values)

        return metrics
```

#### Step 4: Run Evaluation
```bash
python adapters/custom/evaluate_models.py \
    --test-set data/test_molecules.csv \
    --reference admet_ai \
    --output evaluation_results.json
```

**Evaluation criteria**:
- **R² > 0.7**: Good correlation with reference
- **Agreement > 80%**: Most predictions within reasonable range
- **No errors**: Model loads and runs successfully

#### Step 5: Tier Models by Quality
- **Tier 1** (R² > 0.85, agreement > 90%): Use for ensemble
- **Tier 2** (R² > 0.7, agreement > 80%): Use as backup
- **Tier 3** (R² < 0.7 or errors): Retrain or discard

---

## Stage 3: Selective Integration (Week 4)

### Goal
Integrate **only validated, high-quality** custom models

### Approach: Tier-Based Integration

#### Tier 1 Models → Ensemble with ADMET-AI
```python
class HybridADMETAdapter:
    def predict(self, smiles):
        # Get ADMET-AI prediction
        ai_pred = self.admet_ai.predict(smiles)

        # Get custom model prediction (Tier 1 only)
        if 'caco2' in self.tier1_models:
            custom_pred = self.tier1_models['caco2'].predict(smiles)

            # Average the two
            final_pred = (ai_pred['Caco2_Wang'] + custom_pred) / 2
            confidence = 'high'  # Both models agree
        else:
            final_pred = ai_pred['Caco2_Wang']
            confidence = 'medium'  # ADMET-AI only

        return {'caco2': final_pred, 'confidence': confidence}
```

#### Tier 2 Models → Fallback Only
```python
# Only use if ADMET-AI fails or for specialized properties
if admet_ai_failed:
    return tier2_model.predict(smiles)
```

#### Tier 3 Models → Document and Archive
- Document why model didn't meet criteria
- Archive for potential retraining
- Don't integrate into production system

### Integration Checklist
- [ ] Model loads successfully
- [ ] Predictions match evaluation performance
- [ ] Error handling implemented
- [ ] Unit tests written
- [ ] Performance acceptable (<1s per prediction)
- [ ] Documented in adapter metadata

---

## Stage 4: Ensemble Optimization (Week 5+)

### Goal
Optimize ensemble methods for **maximum accuracy**

### Ensemble Strategies

#### 1. Simple Average (Default)
```python
ensemble_pred = (admet_ai_pred + custom_pred) / 2
```
**Pros**: Simple, works well if models have similar accuracy
**Cons**: Doesn't weight by model quality

#### 2. Weighted Average (Better)
```python
# Weight by validation R²
weight_ai = 0.9  # ADMET-AI R² = 0.9
weight_custom = 0.85  # Custom model R² = 0.85

ensemble_pred = (
    (admet_ai_pred * weight_ai) +
    (custom_pred * weight_custom)
) / (weight_ai + weight_custom)
```
**Pros**: Accounts for model quality
**Cons**: Requires validation metrics

#### 3. Confidence-Based Selection (Advanced)
```python
if abs(admet_ai_pred - custom_pred) < 0.1:
    # Models agree - high confidence
    return average(admet_ai_pred, custom_pred), confidence='high'
elif abs(admet_ai_pred - custom_pred) < 0.3:
    # Models disagree slightly - medium confidence
    return admet_ai_pred, confidence='medium'  # Trust ADMET-AI
else:
    # Models disagree significantly - low confidence
    return admet_ai_pred, confidence='low', flag='review_needed'
```

---

## Timeline

| Week | Stage | Focus | Deliverable |
|------|-------|-------|-------------|
| **Week 3** | Stage 1 | ADMET-AI adapter | Working baseline adapter |
| **Week 3.5** | Stage 2 | Model inventory | Complete model catalog |
| **Week 4** | Stage 2 | Model evaluation | Quality metrics for all models |
| **Week 4.5** | Stage 3 | Selective integration | Tier 1 models integrated |
| **Week 5** | Stage 4 | Ensemble optimization | Optimized hybrid adapter |

---

## Decision Framework

### Should I integrate this custom model?

```
┌─────────────────────────────────────┐
│ Does ADMET-AI cover this property?  │
└─────────────────┬───────────────────┘
                  │
         ┌────────┴────────┐
         │                 │
        YES               NO
         │                 │
         ▼                 ▼
    ┌─────────┐      ┌──────────┐
    │ R² > 0.85? │      │ Evaluate │
    └────┬────┘      │ & integrate│
         │           └──────────┘
    ┌────┴────┐
    │         │
   YES       NO
    │         │
    ▼         ▼
┌─────────┐ ┌──────────┐
│ Ensemble│ │ Archive  │
│  it!    │ │ for now  │
└─────────┘ └──────────┘
```

---

## Current Status

### ✅ Completed
- Docker optimization
- ADMET-AI installation and testing
- Coverage matrix documentation
- Staged integration plan

### 🔄 In Progress
- ADMET-AI adapter implementation

### ⏳ Next Steps
1. Complete ADMET-AI adapter
2. Move custom models to `adapters/custom/`
3. Create model inventory
4. Build evaluation framework
5. Evaluate all custom models
6. Integrate Tier 1 models only

---

## Philosophy

> **"Make it work, make it right, make it fast"**
>
> We're at the "make it work" stage. ADMET-AI gives us a working baseline.
> Next, we'll "make it right" by validating custom models.
> Finally, we'll "make it fast" by optimizing ensembles.

**Don't rush to ensemble.** Validate first. Quality over quantity.
