# TDC ADMET Adapter Investigation Report

**Date:** 2025-10-24
**Investigator:** Claude Code (Coder Agent)
**Status:** CRITICAL ISSUE IDENTIFIED

---

## Executive Summary

**FINDING:** The TDC (Therapeutics Data Commons) library does NOT provide pre-trained prediction models. It only provides datasets for training your own models.

**IMPACT:** The current TDC ADMET adapter (`adapters/tdc_admet/adapter.py`) will ALWAYS fail because it calls `model.get_pred()` which does not exist in the TDC library.

**RECOMMENDATION:** Remove the TDC ADMET adapter and rely on ADMET-AI adapter which is already working.

---

## Detailed Investigation

### What is TDC?

TDC (Therapeutics Data Commons) is a **data repository and benchmark platform** for machine learning in drug discovery. It provides:

- Curated datasets for training ML models
- Benchmark splits for reproducible evaluation
- Data loaders and preprocessing utilities
- **NOT** pre-trained prediction models

### Code Analysis

The current TDC adapter (`adapters/tdc_admet/adapter.py`) attempts to use TDC for predictions:

```python
# Line 174 in adapter.py
predictions = model.get_pred([smiles])  # THIS METHOD DOES NOT EXIST!
```

### Actual TDC API

When you load a TDC dataset:

```python
from tdc.single_pred import ADME
data = ADME(name='Caco2_Wang')
```

You get a dataset object with these methods:

- `get_data()` - Returns pandas DataFrame with training data
- `get_split()` - Returns train/val/test splits
- `get_label_meaning()` - Explains what the labels mean
- `convert_to_log()`, `binarize()` - Data preprocessing
- **NO** `get_pred()` or `predict()` methods!

### Experimental Verification

Tested live in Docker container:

```bash
$ docker-compose exec backend python -c "from tdc.single_pred import ADME; ..."
```

**Results:**
- Available methods: 16 methods for data manipulation
- Prediction methods: **NONE**
- Conclusion: TDC is dataset-only

### Data Structure Example

TDC provides datasets in this format:

```
   Drug_ID                                Drug         Y
0  (-)-epicatechin                        ...  -6.220000
1  (2E,4Z,8Z)-N-isobutyldodeca-2,4...    ...  -3.860000
2  codeine                                ...  -4.090000
```

Where:
- `Drug_ID` = compound name
- `Drug` = SMILES string
- `Y` = experimental measurement (e.g., Caco-2 permeability)

This is **training data**, not a prediction model.

---

## Why Was This Missed?

The TDC adapter was likely written based on:

1. **Assumption:** TDC name includes "Commons" - sounds like it provides ready-to-use tools
2. **Similar APIs:** Libraries like RDKit and ChemBL provide direct predictions
3. **Incomplete documentation review:** TDC docs focus on datasets, easy to misinterpret

However, TDC's actual purpose is:

> "Democratizing AI-driven drug discovery by providing benchmarks and datasets"

Not providing prediction services.

---

## Current State Assessment

### Working ADMET Solutions

We already have **ADMET-AI** adapter working (`adapters/admet_ai/adapter.py`):

- Provides 41 pre-trained ADMET prediction models
- Uses deep learning (Chemprop framework)
- Covers absorption, distribution, metabolism, excretion, toxicity
- **Status:** WORKING and tested

### TDC Adapter Status

- **Lines of code:** 363 lines
- **Actual functionality:** ZERO (will always fail)
- **Models claimed:** 11 ADME + 4 Tox models
- **Models that work:** 0

---

## Options Going Forward

### Option 1: REMOVE TDC Adapter (RECOMMENDED)

**Pros:**
- Immediate solution
- We already have ADMET-AI working
- Reduces maintenance burden
- Removes broken code from codebase

**Cons:**
- None - TDC doesn't work anyway

**Action Items:**
1. Delete `adapters/tdc_admet/` directory
2. Remove from adapter registry
3. Update documentation to remove TDC references
4. Keep ADMET-AI as primary ADMET solution

### Option 2: Train Models on TDC Datasets

**Pros:**
- Could leverage TDC's curated datasets
- Would be aligned with TDC's actual purpose

**Cons:**
- Requires training 15+ ML models
- Need labeled training data for each property
- Model deployment complexity
- Time-consuming (weeks of work)
- ADMET-AI already covers this use case

**Estimated Effort:** 2-4 weeks for 2 developers

### Option 3: Keep as Dataset Provider Only

**Pros:**
- Could use TDC for validation benchmarking
- Provides test data for ADMET-AI validation

**Cons:**
- Not an "adapter" in the PharmForge sense
- Confusing to users expecting predictions
- Should be a separate utility, not an adapter

---

## Recommendation

**REMOVE THE TDC ADMET ADAPTER**

### Reasoning

1. **Already have working solution:** ADMET-AI provides 41 ADMET models
2. **TDC can't predict:** No pre-trained models available
3. **Clear scope:** PharmForge adapters should provide predictions, not just data
4. **Reduce confusion:** Users expect adapters to work, not fail silently

### Implementation Plan

1. **Immediate:**
   - Remove `adapters/tdc_admet/` directory
   - Remove TDC from adapter registry
   - Update tests to remove TDC references

2. **Documentation:**
   - Note in docs that ADMET predictions use ADMET-AI
   - Explain TDC's actual purpose (datasets, not predictions)
   - Reference TDC for users who want to train custom models

3. **Future (Optional):**
   - Create benchmarking utilities using TDC datasets
   - Validate ADMET-AI predictions against TDC test sets
   - Publish validation results

---

## Comparison: ADMET-AI vs TDC

| Feature | ADMET-AI | TDC |
|---------|----------|-----|
| **Pre-trained models** | 41 models | 0 models |
| **Prediction capability** | YES | NO |
| **Absorption** | Caco-2, HIA, F20%, F30% | Datasets only |
| **Distribution** | BBB, VDss, fu | Datasets only |
| **Metabolism** | CL, CYP2C9, CYP2D6, CYP3A4 | Datasets only |
| **Excretion** | Half-life, CL-Hepa | Datasets only |
| **Toxicity** | hERG, AMES, DILI, LD50 | Datasets only |
| **Status** | WORKING | BROKEN |
| **Maintenance** | External package | Would need our models |

---

## References

- TDC GitHub: https://github.com/mims-harvard/TDC
- TDC Documentation: https://tdcommons.ai/
- ADMET-AI: https://github.com/swansonk14/admet_ai
- PharmForge ADMET-AI adapter: `adapters/admet_ai/adapter.py`

---

## Conclusion

The TDC ADMET adapter was written with incorrect assumptions about TDC's capabilities. TDC is a valuable resource for **dataset access** and **benchmarking**, but it does NOT provide prediction models.

**Action:** Remove the TDC ADMET adapter and rely on the working ADMET-AI adapter.

**Timeline:** This can be done immediately - it's a deletion, not a migration.

**Risk:** None - TDC adapter doesn't work anyway.

---

**Report Status:** COMPLETE
**Next Steps:** Await human decision on removal vs training custom models
