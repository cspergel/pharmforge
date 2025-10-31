# Custom Tools: Keep or Toss?

## Current PharmForge Coverage (What You Already Have)

### Adapters in Production
1. **PubChem** - Molecular properties (API)
2. **RDKit** - Local descriptor calculations
3. **ChEMBL** - Bioactivity data
4. **ADMET-AI** - 49 ADMET properties (ML-based)

### What They Provide

**RDKit Adapter (Current):**
```python
{
  "molecular_weight": 180.16,
  "logp": 1.23,
  "tpsa": 63.6,
  "h_bond_donors": 1,
  "h_bond_acceptors": 4,
  "rotatable_bonds": 3,
  "num_atoms": 21,
  "num_heavy_atoms": 13,
  "num_aromatic_rings": 1,
  "num_rings": 1,
  "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "lipinski_violations": 0  # Basic Lipinski check
}
```

**ADMET-AI Adapter (Current):**
- **49 properties** covering:
  - Solubility, Caco-2, HIA, P-gp, Bioavailability
  - BBB penetration, Protein binding
  - CYP inhibition/substrates (6 isoforms)
  - Clearance, Half-life
  - Toxicity (AMES, hERG, DILI, LD50, etc.)
  - Physicochemical (MW, LogP, TPSA, QED)

---

## Custom Tools Analysis: Redundant vs Unique

### ❌ **REDUNDANT** (Already Covered)

| Tool | Provides | Already Have | Verdict |
|------|----------|--------------|---------|
| `molecular_weight.py` | MW with thresholds | RDKit + PubChem | ❌ TOSS |
| `logp.py` | LogP scoring | RDKit + PubChem + ADMET-AI | ❌ TOSS |
| `clogp.py` | Crippen LogP | RDKit (uses Crippen) | ❌ TOSS |
| `tpsa.py` | TPSA scoring | RDKit + PubChem + ADMET-AI | ❌ TOSS |
| `hbond.py` | H-bond count | RDKit + PubChem | ❌ TOSS |
| `rotatable_bonds.py` | Flexibility | RDKit + PubChem | ❌ TOSS |
| `ring_count_evaluator.py` | Ring counting | RDKit (num_rings) | ❌ TOSS |
| `basic_descriptors.py` | 15 basic props | RDKit + PubChem | ❌ TOSS |
| `core_descriptors.py` | Aggregated descriptors | RDKit (subset) | ❌ TOSS |

**Total Redundant:** 9/33 tools

---

### ✅ **UNIQUE** (Not in Current Stack)

#### High-Value Additions

| Tool | What It Adds | Why Useful | Keep? |
|------|-------------|-----------|-------|
| **PAINS** | Pan-assay interference detection | Filters false positives in screens | ✅ KEEP |
| **Brenk** | Toxicophore alerts (reactivity, metabolic) | Safety red flags | ✅ KEEP |
| **Veber** | Oral bioavailability rules (TPSA, rot bonds) | Better than basic Lipinski | ✅ KEEP |
| **SAScore** | Synthetic accessibility (1-10 scale) | Retrosynthesis prioritization | ✅ KEEP |

**Rationale:**
- **PAINS/Brenk:** ADMET-AI doesn't flag interference compounds or reactive groups
- **Veber:** More predictive than Lipinski for oral drugs (adds TPSA + rotatable bonds criteria)
- **SAScore:** ADMET-AI doesn't predict synthesis difficulty

#### Medium-Value Additions

| Tool | What It Adds | Why Useful | Keep? |
|------|-------------|-----------|-------|
| **Ghose** | Drug-likeness (broader than Lipinski) | Alternative filter | ⚠️ MAYBE |
| **Muegge** | Drug-likeness (Bayer rules) | Alternative filter | ⚠️ MAYBE |
| **Egan** | CNS penetration filter | If CNS drugs | ⚠️ MAYBE |
| **Lead-likeness** | Lead optimization rules (MW<350) | If early discovery | ⚠️ MAYBE |
| **RO3** | Fragment screening rules | If fragment-based | ⚠️ MAYBE |
| **CNS MPO** | 6-parameter CNS score | If CNS drugs | ⚠️ MAYBE |
| **Fsp3** | Aliphatic fraction (3D shape) | Novelty metric | ⚠️ MAYBE |
| **Aromatic proportion** | Aromaticity ratio | Flat molecule filter | ⚠️ MAYBE |
| **Macrocycle filter** | Flags large rings | Specialized chemistry | ⚠️ MAYBE |
| **Advanced 3D descriptors** | Shape (asphericity, PMI, etc.) | 3D structure-based work | ⚠️ MAYBE |
| **BDDCS** | Biopharmaceutics classification | Formulation insights | ⚠️ MAYBE |

**Rationale:**
- **Use case dependent:** If you're doing CNS drugs → keep CNS MPO, Egan
- **Discovery stage:** If lead optimization → keep Lead-likeness, Fsp3
- **Default:** Can skip these for general-purpose screening

#### Low-Value / Broken

| Tool | Status | Why Skip | Keep? |
|------|--------|----------|-------|
| **novelty.py** | Empty stub | Not implemented | ❌ TOSS |
| **similarity.py** | Empty stub | Not implemented | ❌ TOSS |
| **synthetic_pathway_complexity** | Basic topology | SAScore is better | ❌ TOSS |
| **utility_utils.py** | Empty | Not implemented | ❌ TOSS |
| **rule_utils.py** | Empty | Not implemented | ❌ TOSS |

**Total Unique:** 20/33 tools (4 high-value, 11 medium-value, 5 low-value)

---

## Decision Matrix

### Scenario 1: General-Purpose Drug Discovery (Default)

**Keep Only:**
- ✅ PAINS (interference detection)
- ✅ Brenk (toxicophore alerts)
- ✅ Veber (oral bioavailability)
- ✅ SAScore (synthetic accessibility)

**Total:** 4 tools
**Effort:** ~3 hours (need fpscores.pkl.gz for SAScore)
**Value:** High - covers safety + feasibility gaps in ADMET-AI

---

### Scenario 2: CNS Drug Discovery

**Keep:**
- ✅ PAINS, Brenk, Veber, SAScore (core 4)
- ✅ CNS MPO (6-parameter CNS penetration score)
- ✅ Egan (CNS bioavailability filter)

**Total:** 6 tools
**Effort:** ~4 hours
**Value:** Medium-High for CNS-focused projects

---

### Scenario 3: Lead Optimization Focus

**Keep:**
- ✅ PAINS, Brenk, Veber, SAScore (core 4)
- ✅ Lead-likeness (MW<350, simpler rules)
- ✅ Fsp3 (3D character for binding)

**Total:** 6 tools
**Effort:** ~4 hours
**Value:** Medium for early-stage optimization

---

### Scenario 4: Comprehensive Filtering (Max Coverage)

**Keep:**
- ✅ Core 4 (PAINS, Brenk, Veber, SAScore)
- ✅ Drug-likeness variants (Ghose, Muegge, Egan, Lead-likeness, RO3, BDDCS)
- ✅ Structural metrics (CNS MPO, Fsp3, Aromatic proportion, Macrocycle)
- ✅ 3D descriptors (Advanced descriptors)

**Total:** 16 tools
**Effort:** ~10 hours
**Value:** Medium - adds filtering diversity but some redundancy

---

## Recommended Path

### ✅ **Conservative Recommendation (Week 1)**

**Integrate ONLY the 4 high-value tools:**

1. **PAINS** - Prevents assay interference compounds
2. **Brenk** - Flags reactive/toxic substructures
3. **Veber** - Superior oral bioavailability filter vs basic Lipinski
4. **SAScore** - Synthetic feasibility scoring (requires fpscores.pkl.gz)

**Why These 4?**
- Fill gaps ADMET-AI doesn't cover (interference, reactivity, synthesis)
- Proven utility in industry (PAINS/Brenk are standard filters)
- Low integration effort (3 hours + data file)
- No redundancy with current stack

**Archive the rest** unless you have specific use case (CNS, fragments, etc.)

---

### ⚠️ **If You Need More (Week 2+)**

**Add as needed:**
- CNS projects → CNS MPO + Egan
- Lead optimization → Lead-likeness + Fsp3
- Fragment screening → RO3
- Diversity library → Ghose + Muegge + BDDCS

---

## Implementation Plan

### Option A: Minimal (Recommended)

```bash
# 1. Get SAScore data file
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SAscore/fpscores.pkl.gz
mkdir -p adapters/custom/resources/
mv fpscores.pkl.gz adapters/custom/resources/

# 2. Keep only 4 tools
cd adapters/custom/rule_based/
mv pains.py brenk.py veber.py ../filters_keep/
cd ../utility/
mv sascore.py ../filters_keep/

# 3. Archive the rest
mv adapters/custom/* adapters/custom_DEPRECATED/

# 4. Integrate into pipeline
# (register 4 evaluators in adapter registry)
```

**Effort:** 3-4 hours
**Files to keep:** 4
**Files to archive:** 29

---

### Option B: Toss Everything

If you want to **simplify and focus on ADMET-AI**:

```bash
# Archive entire custom folder
mv adapters/custom adapters/custom_DEPRECATED

# Rationale:
# - ADMET-AI covers 90% of ADMET needs
# - RDKit covers basic descriptors
# - PubChem/ChEMBL cover properties/bioactivity
# - Focus on core pipeline (docking, retrosynthesis, ranking)
```

**Effort:** 0 hours
**Benefit:** Cleaner codebase, less maintenance

---

## Cost-Benefit Analysis

### Keeping 4 Core Tools

**Benefits:**
- PAINS: Prevent ~5-10% false positives in HTS
- Brenk: Flag reactive groups missed by ADMET-AI
- Veber: Better oral bioavailability prediction
- SAScore: Prioritize synthetically accessible hits

**Costs:**
- 3-4 hours integration
- 1 MB data file (fpscores.pkl.gz)
- 4 more modules to maintain
- Testing overhead

**ROI:** Medium-High - worth it if doing HTS or lead prioritization

---

### Tossing Everything

**Benefits:**
- Zero integration time
- Cleaner codebase
- Focus on core PharmForge features

**Costs:**
- Miss PAINS/Brenk filtering (industry standard)
- No synthetic accessibility scoring
- Slightly less rigorous oral bioavailability check

**ROI:** High if prioritizing speed-to-MVP

---

## My Honest Recommendation

### 🎯 **For MVP: TOSS EVERYTHING**

**Why?**
1. Your current stack (PubChem + RDKit + ChEMBL + ADMET-AI) is **already comprehensive**
2. ADMET-AI covers toxicity (AMES, hERG, DILI, LD50) - overlaps with Brenk
3. RDKit has basic Lipinski - close enough to Veber for MVP
4. You can add PAINS/SAScore later if users request it
5. **Focus on core differentiators**: NL orchestration, docking, retrosynthesis, ranking

### 🔧 **Post-MVP: Add Core 4 if Needed**

If beta users say:
- "I'm getting false positives in my assays" → Add PAINS
- "Need to filter reactive compounds" → Add Brenk
- "Want synthesis difficulty scoring" → Add SAScore
- "Better oral drug filtering" → Add Veber

**Then integrate the 4 core tools** (3-hour task).

---

## Action Items

### Option A: Keep Core 4 (Conservative)
```bash
# 1. Download SAScore data
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SAscore/fpscores.pkl.gz

# 2. Keep: pains.py, brenk.py, veber.py, sascore.py
# 3. Archive: remaining 29 files

# 4. Test with aspirin:
python -c "
from adapters.custom.rule_based.pains import PAINSEvaluator
from adapters.custom.rule_based.brenk import BrenkAlertEvaluator
print(PAINSEvaluator().evaluate('CC(=O)Oc1ccccc1C(=O)O'))
"
```

### Option B: Toss Everything (Recommended for MVP)
```bash
# Archive entire folder
mv adapters/custom adapters/custom_DEPRECATED_20251024

# Document decision
echo "Archived: Redundant with RDKit + ADMET-AI. Revisit post-MVP if users need PAINS/Brenk/SAScore." > adapters/CUSTOM_TOOLS_ARCHIVED.txt
```

---

## Bottom Line

**Current stack covers 90% of needs.** These custom tools add **marginal value** for MVP.

**My vote:** Archive everything now, focus on core pipeline. Add PAINS/Brenk/SAScore later only if users explicitly request toxicophore filtering or synthesis scoring.

**If you disagree:** Keep the core 4 (PAINS, Brenk, Veber, SAScore) and archive the rest.

**Either way:** Don't integrate all 33 tools - that's complexity without clear ROI.
