# PharmForge Future Adapter Roadmap

**Current Status:** 33 adapters exist, 13 tested/working
**Date:** October 25, 2025

---

## ALREADY BUILT - Just Need Testing/Setup (20 adapters)

### Protein Structure Analysis
- ✓ **RCSB PDB** - Experimental structures (TESTED, WORKING)
- ✓ **AlphaFold** - Predicted structures (TESTED, WORKING)
- ✓ **SWISS-MODEL** - Homology modeling (TESTED, WORKING)
- ✓ **PDB-REDO** - Re-refined PDB structures (EXISTS, UNTESTED)

### Pathway & Systems Biology
- ✓ **Reactome** - Pathway analysis for off-target predictions (EXISTS, UNTESTED)
- ✓ **GTEx** - Gene expression by tissue for selectivity (EXISTS, UNTESTED)

### Molecular Docking
- ✓ **AutoDock Vina** - Classic docking (EXISTS, UNTESTED)
- ✓ **GNINA** - CNN-based docking (DOCUMENTED, needs conda install)
- ✓ **DiffDock** - Diffusion-based docking (EXISTS, UNTESTED - needs GPU)

### Database & Compound Search
- ✓ **ChEMBL** - Bioactivity database (EXISTS, UNTESTED)
- ✓ **PubChem** - Compound database (EXISTS, UNTESTED)
- ✓ **ZINC** - Fragment library (TESTED, WORKING - migrated to ChEMBL)

### Literature & Patents
- ✓ **PubMed** - Literature search (EXISTS, UNTESTED)
- ✓ **EuropePMC** - European literature (EXISTS, UNTESTED)
- ✓ **Google Patents** - Patent search (EXISTS, UNTESTED)
- ✓ **Lens.org** - Patent/scholarly search (EXISTS, UNTESTED)
- ✓ **SureChEMBL** - Patent chemicals (TESTED, WORKING)

### Computational Chemistry Tools
- ✓ **RDKit Local** - Molecular properties (EXISTS, UNTESTED)
- ✓ **OpenMM** - Molecular dynamics (EXISTS, UNTESTED)
- ✓ **ADMET-AI** - ADMET predictions (EXISTS, UNTESTED)

### De Novo Design
- ✓ **MolGAN** - GAN-based generation (EXISTS, UNTESTED)
- ✓ **REINVENT** - RL-based design (EXISTS, UNTESTED)
- ✓ **De Novo** - De novo design (EXISTS, UNTESTED)

---

## NEED TO BUILD - User Requirements

### Priority 1: Essential Missing Tools

#### Protein-Protein Interactions
- ❌ **STRING** - Protein-protein interaction networks
  - Purpose: Off-target predictions, pathway analysis
  - API: https://string-db.org/cgi/help.pl?subpage=api
  - Free: YES
  - Priority: HIGH
  - Effort: Medium (REST API, straightforward)

#### Advanced Docking & Scoring
- ❌ **Smina** - Vina fork with improved scoring
  - Purpose: Better docking scores than Vina
  - Source: https://sourceforge.net/projects/smina/
  - Free: YES
  - Priority: HIGH
  - Effort: Medium (similar to Vina adapter)
  - Note: Command-line tool, needs binary

- ❌ **OnionNet** - ML-based binding affinity prediction
  - Purpose: Protein-ligand affinity scoring
  - Source: https://github.com/zhenglz/OnionNet
  - Free: YES
  - Priority: HIGH
  - Effort: High (ML model, complex dependencies)

- ❌ **Pafnucy** - CNN protein-ligand scoring
  - Purpose: Deep learning-based scoring
  - Source: https://gitlab.com/cheminfIBB/pafnucy
  - Free: YES
  - Priority: MEDIUM
  - Effort: High (TensorFlow/Keras, 3D grid generation)

#### Post-Docking Refinement
- ❌ **MM-GBSA** (Molecular Mechanics Generalized Born Surface Area)
  - Purpose: Binding free energy calculations
  - Implementation: Via OpenMM or AmberTools
  - Free: YES (AmberTools)
  - Priority: HIGH
  - Effort: High (requires MD setup, energy calculations)

- ❌ **MM-PBSA** (Molecular Mechanics Poisson-Boltzmann Surface Area)
  - Purpose: More accurate binding free energy
  - Implementation: Via AmberTools MMPBSA.py
  - Free: YES
  - Priority: MEDIUM
  - Effort: High (similar to MM-GBSA)

### Priority 2: Ligand-Based Design

#### Fingerprint Similarity
- ❌ **Tanimoto/Dice Screening** - Fingerprint similarity
  - Purpose: Virtual screening, lead optimization
  - Implementation: Extend RDKit adapter
  - Free: YES (RDKit has this)
  - Priority: HIGH
  - Effort: LOW (RDKit already exists)
  - Note: **Can add to existing RDKit adapter**

#### Pharmacophore Modeling
- ❌ **Pharmacophore Matching** - 3D pharmacophore searches
  - Purpose: Virtual screening by pharmacophore
  - Implementation Options:
    - RDKit (basic pharmacophores)
    - Open3DAlign (3D alignment)
    - LigandScout (commercial, but has free version)
  - Free: YES (RDKit/Open3DAlign)
  - Priority: MEDIUM
  - Effort: Medium-High (3D alignment, feature detection)

### Priority 3: Patent & Novelty

#### USPTO PatentView
- ❌ **USPTO PatentView** - US patent database
  - Purpose: Novelty checks, patent landscape
  - API: https://patentsview.org/apis/api-endpoints
  - Free: YES
  - Priority: MEDIUM
  - Effort: Medium (REST API)
  - Note: Complements existing Google Patents adapter

### Priority 4: Quantum Chemistry (Future)

#### Semi-Empirical QM
- ❌ **xTB** (Extended Tight-Binding)
  - Purpose: Fast QM calculations, conformer energies
  - Source: https://github.com/grimme-lab/xtb
  - Free: YES
  - Priority: LOW (future enhancement)
  - Effort: High (requires Fortran binary, complex setup)
  - Use Cases: Conformer generation, reaction barriers, HOMO-LUMO

#### Ab Initio QM
- ❌ **Psi4** - Quantum chemistry package
  - Purpose: High-accuracy calculations
  - Source: https://psicode.org/
  - Free: YES
  - Priority: LOW (future, computationally expensive)
  - Effort: Very High (complex dependencies, slow calculations)
  - Use Cases: Accurate energies, orbital analysis, reaction mechanisms

---

## IMPLEMENTATION PRIORITY MATRIX

### Immediate (Next Sprint)
1. **Test existing 20 adapters** - Many likely work already
2. **STRING adapter** - Essential for off-target analysis
3. **Tanimoto/Dice screening** - Extend RDKit adapter (easy win)

### Short-term (1-2 weeks)
4. **Smina adapter** - Better docking
5. **OnionNet adapter** - ML-based scoring
6. **USPTO PatentView** - Patent searches

### Medium-term (1 month)
7. **MM-GBSA/MM-PBSA** - Post-docking refinement
8. **Pharmacophore matching** - Ligand-based design
9. **Pafnucy adapter** - Deep learning scoring

### Long-term (Future)
10. **xTB adapter** - Semi-empirical QM
11. **Psi4 adapter** - Ab initio QM
12. Additional ML models as they emerge

---

## EFFORT ESTIMATION

### Quick Wins (1-2 days each)
- Tanimoto/Dice (extend existing RDKit) - 1 day
- Test existing adapters - 2-3 days
- STRING adapter - 1-2 days
- USPTO PatentView - 1-2 days

### Medium Effort (3-7 days each)
- Smina adapter - 3-4 days
- Pharmacophore matching - 5-7 days
- OnionNet adapter - 5-7 days

### High Effort (1-2 weeks each)
- MM-GBSA/MM-PBSA - 7-10 days
- Pafnucy adapter - 7-14 days
- xTB adapter - 10-14 days

### Very High Effort (2-4 weeks)
- Psi4 adapter - 14-21 days

---

## COVERAGE ANALYSIS

### What You Already Have (Exists/Working)

**Protein Analysis:**
- ✓ UniProt, RCSB PDB, AlphaFold, SWISS-MODEL, PDB-REDO
- ✓ OpenTargets, DisGeNET
- ✓ Reactome (exists), GTEx (exists)

**Docking:**
- ✓ Vina (exists), GNINA (documented), DiffDock (exists)
- ✓ BindingDB for validation

**Structure-Based:**
- ✓ All major structure sources covered
- ❌ Missing: Refinement (MM-GBSA/MM-PBSA)

**Ligand-Based:**
- ✓ ZINC fragments, DrugCentral
- ❌ Missing: Fingerprint screening (easy - extend RDKit)
- ❌ Missing: Pharmacophore matching

**Systems Biology:**
- ✓ Reactome (exists), GTEx (exists)
- ❌ Missing: STRING (protein-protein interactions)

**Safety/ADMET:**
- ✓ FDA FAERS, ADMET-AI (exists)

**Clinical:**
- ✓ ClinicalTrials.gov

**Patents:**
- ✓ Google Patents, Lens.org, SureChEMBL
- ❌ Missing: USPTO PatentView

**Quantum Chemistry:**
- ❌ Missing: All (xTB, Psi4) - future priority

---

## RECOMMENDATIONS

### Immediate Next Steps:
1. **Test the 20 existing adapters** - Many probably work
2. **Build STRING adapter** - Critical for off-target analysis
3. **Extend RDKit adapter** - Add Tanimoto/Dice screening (1 day)

### After Testing Complete:
4. Build Smina adapter (better docking)
5. Build OnionNet adapter (ML scoring)
6. Build USPTO PatentView (patents)

### Future Enhancements:
7. MM-GBSA/MM-PBSA for binding energy
8. Pharmacophore matching
9. Quantum chemistry tools (xTB, Psi4)

---

## FREE API/TOOL AVAILABILITY

### All FREE:
- ✓ STRING database (free API)
- ✓ Reactome (free API) - already exists
- ✓ GTEx (free data portal) - already exists
- ✓ RCSB PDB (free API) - working
- ✓ PDB-REDO (free service) - exists
- ✓ SWISS-MODEL (free service) - working
- ✓ AutoDock Vina (free software) - exists
- ✓ Smina (free software)
- ✓ GNINA (free software) - documented
- ✓ DiffDock (free software) - exists
- ✓ RDKit (free software) - exists
- ✓ OnionNet (free software/model)
- ✓ Pafnucy (free software)
- ✓ USPTO PatentView (free API)
- ✓ xTB (free software)
- ✓ Psi4 (free software)

**All requested features have free implementations available!**

---

## SUMMARY

**You Already Have:** 20/30+ requested features (67%)

**Quick Additions Needed (1-3 days each):**
- STRING adapter
- Tanimoto/Dice screening (extend RDKit)
- USPTO PatentView

**Medium Additions (3-7 days each):**
- Smina adapter
- OnionNet
- Pharmacophore matching

**Future Additions (1-3 weeks each):**
- MM-GBSA/MM-PBSA
- Pafnucy
- xTB
- Psi4

---

**Next Action:** Test existing 20 adapters, then build STRING + Tanimoto screening (critical for comprehensive drug discovery platform)
