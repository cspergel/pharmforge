# PharmForge Custom Filters & Utilities

**Version:** 1.0.0
**Status:** Production-ready (requires setup for SAScore)

## Overview

These 4 tools complement the existing ADMET-AI and RDKit adapters by providing specialized filtering and scoring that isn't covered by the baseline stack.

---

## Tools Included

### 1. PAINS (Pan-Assay Interference Compounds)
**File:** `filters/toxicity/pains.py`

**Purpose:** Detects compounds that produce false positives in biochemical assays.

**Why Keep:** ADMET-AI doesn't filter for assay interference - this is an industry-standard safety check.

**Usage:**
```python
from adapters.custom.filters.toxicity import PAINSEvaluator

evaluator = PAINSEvaluator()
result = evaluator.evaluate(mol)  # mol = RDKit Mol object
# result = {"score": 1.0, "pass": True, "tags": ["pains_pass"], ...}
```

**Output:**
- `score`: 1.0 if clean, penalty for alerts (0.2 per alert)
- `pass`: True if no alerts
- `tags`: ["pains_pass"] or ["pains_fail", "pains_alert_1", ...]
- `raw.alert_descriptions`: Names of matched PAINS filters

**Reference:** Baell & Holloway, J. Med. Chem. 2010

---

### 2. Brenk Alerts (Toxicophore Filter)
**File:** `filters/toxicity/brenk.py`

**Purpose:** Flags reactive, toxic, or metabolically unstable substructures.

**Why Keep:** ADMET-AI predicts toxicity outcomes, but doesn't flag specific reactive groups. Brenk provides structural alerts.

**Usage:**
```python
from adapters.custom.filters.toxicity import BrenkAlertEvaluator

evaluator = BrenkAlertEvaluator()
result = evaluator.evaluate(mol)
# result = {"score": 1.0, "pass": True, "tags": ["brenk_low_risk"], ...}
```

**Output:**
- `score`: 1.0 if clean, decreases with severity
- `pass`: True if risk_level == "low"
- `tags`: ["brenk_low_risk"], ["brenk_medium_risk"], or ["brenk_high_risk"]
- `raw.alerts`: List of alerts with categories (reactivity, metabolic, toxicity, aggregation)
- `raw.categories`: Grouped by category with max severity

**Alert Categories:**
- **Reactivity:** Michael acceptors, aldehydes, quinones
- **Metabolic:** Nitro groups, sulfoxides
- **Toxicity:** Halogens, epoxides, aziridines
- **Aggregation:** Crown ethers, polycycles

**Reference:** Brenk et al., ChemMedChem 2008

---

### 3. Veber Rules (Oral Bioavailability)
**File:** `filters/druglikeness/veber.py`

**Purpose:** Predicts oral bioavailability using TPSA and rotatable bonds.

**Why Keep:** More predictive than Lipinski alone. RDKit has basic Lipinski, but Veber adds TPSA criterion which is critical for CNS drugs.

**Usage:**
```python
from adapters.custom.filters.druglikeness import VeberEvaluator

evaluator = VeberEvaluator()
result = evaluator.evaluate(mol)
# result = {"score": 0.85, "pass": True, "tags": ["PassVeber"], ...}
```

**Output:**
- `score`: Normalized (0-1) based on TPSA and rotatable bonds
- `pass`: True if all 3 rules satisfied
- `tags`: ["PassVeber"], ["PartialVeber"], or ["FailVeber"]
- `summary.oral_bioavailability`: "high", "moderate", or "low"
- `raw.rules_check`: Boolean for each rule

**Rules:**
1. Rotatable bonds ≤ 10
2. TPSA ≤ 140 Ų
3. Total H-bonds (donors + acceptors) ≤ 12

**Reference:** Veber et al., J. Med. Chem. 2002

---

### 4. SAScore (Synthetic Accessibility Score)
**File:** `synthesis/sascore.py`

**Purpose:** Estimates synthetic difficulty from 1 (easy) to 10 (very hard).

**Why Keep:** ADMET-AI doesn't predict synthesis difficulty. Critical for prioritizing synthetically accessible hits.

**⚠️ REQUIRES SETUP - See below**

**Usage:**
```python
from adapters.custom.synthesis import SynthesisAccessibilityEvaluator

evaluator = SynthesisAccessibilityEvaluator()
result = evaluator.evaluate([mol])  # Takes list of UnifiedMolecule objects
# result = [{"score": 2.5, "tags": ["EasySynthesis"], ...}]
```

**Output:**
- `score`: 1-10 (1=easy, 10=very hard)
- `tags`: ["EasySynthesis"] (< 3), ["ModerateSynthesis"] (3-7), or ["HardToSynthesize"] (> 7)
- `summary.synthesis_difficulty`: "low", "moderate", or "high"
- `raw.sascore`: Exact score
- `raw.fallback_fragments`: Fragments without known scores
- `raw.total_fragments`: Total fragments analyzed

**Algorithm:**
Combines 3 components:
1. **Fragment score:** Frequency of molecular fragments in known compounds
2. **Complexity penalties:** Atom count, chiral centers, spiro/bridgehead atoms, macrocycles
3. **Symmetry bonus:** Bit diversity vs atom count

**Reference:** Ertl & Schuffenhauer, J. Cheminform. 2009

---

## Setup Instructions

### Quick Setup (All Tools Except SAScore)

No setup needed! PAINS, Brenk, and Veber use RDKit's built-in filter catalogs.

```bash
# Test PAINS
docker-compose exec backend python -c "
from adapters.custom.filters.toxicity import PAINSEvaluator
from rdkit import Chem
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
print(PAINSEvaluator().evaluate(mol))
"
```

---

### SAScore Setup (Required Data File)

**Status:** ❌ Missing `fpscores.pkl.gz`

**Impact:** SAScore will crash on initialization without this file.

**Solution:**

1. **Download the fragment frequency database:**

```bash
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SAscore/fpscores.pkl.gz

# Or download from backup:
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SA_Score/fpscores.pkl.gz
```

2. **Place in resources directory:**

```bash
mkdir -p adapters/custom/resources
mv fpscores.pkl.gz adapters/custom/resources/
```

3. **Verify file is correct:**

```bash
# Should be ~1 MB compressed
ls -lh adapters/custom/resources/fpscores.pkl.gz
# Output: -rw-r--r-- 1 user user 1.0M fpscores.pkl.gz
```

4. **Test SAScore:**

```bash
docker-compose exec backend python -c "
from adapters.custom.synthesis import SynthesisAccessibilityEvaluator
from pharmforge_core.core.molecule import UnifiedMolecule

evaluator = SynthesisAccessibilityEvaluator()
mol = UnifiedMolecule('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
result = evaluator.evaluate([mol])
print(f'SAScore: {result[0][\"score\"]} (1=easy, 10=hard)')
"
```

**Expected Output:**
```
[SAScore] Loaded 85479 fragment scores.
SAScore: 2.3 (1=easy, 10=hard)
```

---

## Integration with PharmForge Pipeline

### Option 1: Standalone Usage

Use filters individually:

```python
from adapters.custom.filters.toxicity import PAINSEvaluator, BrenkAlertEvaluator
from adapters.custom.filters.druglikeness import VeberEvaluator

pains = PAINSEvaluator()
brenk = BrenkAlertEvaluator()
veber = VeberEvaluator()

# Screen molecule
pains_result = pains.evaluate(mol)
if not pains_result["pass"]:
    print("FAIL: PAINS interference detected")

brenk_result = brenk.evaluate(mol)
if brenk_result["raw"]["risk_assessment"]["level"] == "high":
    print("FAIL: High-risk toxicophore")

veber_result = veber.evaluate(mol)
if veber_result["summary"]["oral_bioavailability"] == "low":
    print("WARN: Poor oral bioavailability")
```

---

### Option 2: Add to Adapter Registry

Register as adapters:

```python
# backend/core/adapter_registry.py

from adapters.custom.filters.toxicity import PAINSEvaluator, BrenkAlertEvaluator
from adapters.custom.filters.druglikeness import VeberEvaluator
from adapters.custom.synthesis import SynthesisAccessibilityEvaluator

def register_custom_filters():
    registry.register(PAINSEvaluator())
    registry.register(BrenkAlertEvaluator())
    registry.register(VeberEvaluator())
    registry.register(SynthesisAccessibilityEvaluator())
```

---

## Recommended Filtering Pipeline

**For high-throughput screening:**

```python
def filter_compound(mol):
    """Return True if compound passes all filters"""

    # 1. PAINS filter (critical - eliminates false positives)
    if not pains.evaluate(mol)["pass"]:
        return False, "PAINS alert"

    # 2. Brenk filter (safety)
    brenk_result = brenk.evaluate(mol)
    if brenk_result["raw"]["risk_assessment"]["level"] == "high":
        return False, "High toxicity risk"

    # 3. Veber filter (oral bioavailability)
    veber_result = veber.evaluate(mol)
    if veber_result["summary"]["veber_violations"] > 1:
        return False, "Poor oral bioavailability"

    # 4. SAScore filter (synthetic accessibility)
    sascore_result = sascore.evaluate([mol])[0]
    if sascore_result["score"] > 6:
        return False, "Difficult to synthesize"

    return True, "PASS"
```

---

## File Structure

```
adapters/custom/
├── README.md                           # This file
├── __init__.py
├── filters/
│   ├── __init__.py
│   ├── toxicity/
│   │   ├── __init__.py
│   │   ├── pains.py                   # PAINS filter
│   │   └── brenk.py                   # Brenk alerts
│   └── druglikeness/
│       ├── __init__.py
│       └── veber.py                   # Veber rules
├── synthesis/
│   ├── __init__.py
│   └── sascore.py                     # Synthetic accessibility
└── resources/
    └── fpscores.pkl.gz                # ⚠️ DOWNLOAD REQUIRED
```

---

## Comparison to Existing Adapters

| Feature | RDKit | ADMET-AI | Custom Tools |
|---------|-------|----------|--------------|
| **Basic descriptors** | ✅ MW, LogP, TPSA | ✅ | ❌ Redundant |
| **Lipinski check** | ✅ Basic | ❌ | ❌ Redundant |
| **PAINS filtering** | ❌ | ❌ | ✅ **Unique** |
| **Toxicophore alerts** | ❌ | ⚠️ Outcomes only | ✅ **Unique** |
| **Veber bioavailability** | ⚠️ Partial | ❌ | ✅ **Unique** |
| **Synthesis difficulty** | ❌ | ❌ | ✅ **Unique** |
| **ADMET predictions** | ❌ | ✅ 49 properties | ❌ |

**Verdict:** These 4 tools fill specific gaps not covered by RDKit or ADMET-AI.

---

## Performance Notes

- **PAINS:** ~0.01s per molecule (RDKit substructure matching)
- **Brenk:** ~0.01s per molecule (RDKit filter catalog)
- **Veber:** ~0.001s per molecule (simple descriptor calculations)
- **SAScore:** ~0.05s per molecule (fragment fingerprint lookup)

**Total overhead:** ~0.07s per molecule for all 4 filters

---

## Troubleshooting

### SAScore fails with "FileNotFoundError: fpscores.pkl.gz not found"

**Fix:** Download the file (see Setup Instructions above)

```bash
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SAscore/fpscores.pkl.gz
mv fpscores.pkl.gz adapters/custom/resources/
```

---

### Import errors (cannot import PAINSEvaluator)

**Fix:** Ensure __init__.py files are present in all directories

```bash
# Check structure
find adapters/custom -name "__init__.py"

# Should show:
# adapters/custom/__init__.py
# adapters/custom/filters/__init__.py
# adapters/custom/filters/toxicity/__init__.py
# adapters/custom/filters/druglikeness/__init__.py
# adapters/custom/synthesis/__init__.py
```

---

### Deprecated project dependencies (pharmforge_core)

**Issue:** Old evaluators reference `pharmforge_core.evaluators.base`

**Fix:** These have been cleaned up. The 4 keeper tools are standalone and only depend on RDKit.

---

## Future Enhancements

**Not included (but could add later):**
- Novelty scoring (Tanimoto distance to ChEMBL)
- Similarity searches
- Additional drug-likeness rules (Ghose, Muegge, CNS MPO)
- 3D shape descriptors

**For now:** Focus on core 4 - they provide the most value.

---

## References

1. **PAINS:** Baell JB, Holloway GA. "New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays" *J. Med. Chem.* 2010, 53, 7, 2719–2740

2. **Brenk:** Brenk R, et al. "Lessons Learnt from Assembling Screening Libraries for Drug Discovery for Neglected Diseases" *ChemMedChem* 2008, 3, 435-444

3. **Veber:** Veber DF, et al. "Molecular Properties That Influence the Oral Bioavailability of Drug Candidates" *J. Med. Chem.* 2002, 45, 12, 2615–2623

4. **SAScore:** Ertl P, Schuffenhauer A. "Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions" *J. Cheminform.* 2009, 1:8

---

## License

These tools are derived from published algorithms and RDKit implementations (BSD license).

---

## Support

For issues or questions, check:
1. This README
2. PharmForge main documentation
3. RDKit documentation for filter catalogs
