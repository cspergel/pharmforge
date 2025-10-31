# PharmForge Adapter Expansion Plan
**Date:** 2025-10-30
**Goal:** Expand from 57 → 63+ adapters with free/easy-to-install tools from awesome-drug-discovery

---

## Current Status
- **Existing adapters:** 57
- **Source:** [awesome-drug-discovery](https://github.com/yboulaamane/awesome-drug-discovery)
- **Docker/Backend:** ✅ Ready (auto-installs from requirements.txt)

---

## Priority 1: Easy Pip-Installable (Target: +6 adapters = 63 total)

### Batch 1: Molecular Dynamics & Visualization (2 adapters)
1. **gmx_MMPBSA**
   - Category: MD Free Energy Analysis
   - Install: `pip install gmx_MMPBSA`
   - Use case: Calculate binding free energies from MD trajectories
   - Priority: HIGH - fills gap in MD analysis

2. **PyMOL (open-source version)**
   - Category: Molecular Visualization
   - Install: `pip install pymol-open-source`
   - Use case: Generate publication-quality molecular images
   - Priority: HIGH - essential visualization tool

### Batch 2: AutoML & Conformer Generation (2 adapters)
3. **Auto-sklearn**
   - Category: AutoML
   - Install: `pip install auto-sklearn`
   - Use case: Automated ML pipeline optimization (complements TPOT)
   - Priority: MEDIUM - provides alternative to TPOT

4. **MolScrub**
   - Category: Conformer Generation
   - Install: `pip install molscrub`
   - Use case: Generate clean, standardized conformers
   - Priority: MEDIUM - enhances docking prep

### Batch 3: Advanced ML & Web APIs (2 adapters)
5. **Oloren ChemEngine**
   - Category: Property Prediction
   - Install: `pip install olorenchemengine`
   - Use case: State-of-the-art molecular property prediction
   - Priority: HIGH - cutting-edge ML models

6. **SwissADME**
   - Category: ADMET Prediction (Web API)
   - Install: No install (web scraping or API)
   - Use case: Free web-based ADMET predictions
   - Priority: LOW - web scraping may be fragile

---

## Priority 2: Specialized Tools (Consider for +4 more = 67 total)

7. **RDKit Contrib** (additional features)
   - Already have rdkit_local, could expand with contrib modules

8. **Psi4** (Quantum Chemistry)
   - Install: `conda install psi4` (conda only)
   - Use case: QM calculations for small molecules
   - Priority: LOW - conda dependency

9. **KNIME Analytics Platform** (workflow tool)
   - Requires GUI/binary - SKIP

10. **ChemAxon** (proprietary)
    - Requires license - SKIP

---

## Priority 3: Database/API Adapters (Free APIs)

11. **DrugBank** (partial free API)
12. **BindingDB** (already have!)
13. **ZINC** (already have zinc_fragments!)
14. **ChemSpider** (API available)
15. **PDBe** (European PDB)

---

## Implementation Order

### Phase 1: Quick Wins (TODAY)
1. gmx_MMPBSA
2. PyMOL
3. MolScrub

**Estimated time:** 2-3 hours

### Phase 2: Advanced (NEXT)
4. Oloren ChemEngine
5. Auto-sklearn
6. SwissADME (if time permits)

**Estimated time:** 2-3 hours

---

## Exclusions (Complex Setup / Not Free)

❌ **Smina** - Requires compilation
❌ **AutoDock-GPU** - Requires CUDA compilation
❌ **LAMMPS** - Source compilation only
❌ **GROMACS** - Binary/source installation
❌ **Avogadro** - GUI binary download
❌ **VMD** - Proprietary binary
❌ **ChemAxon** - Proprietary license required

---

## Success Criteria

- [x] Requirements.txt updated with new dependencies
- [ ] All 6 adapters created with adapter.py, README.md, examples
- [ ] All adapters registered in adapter_registry.py
- [ ] Integration test passes for new adapters
- [ ] Count reaches 63 adapters

---

## Notes

- Docker will auto-install new dependencies on rebuild
- Each adapter should follow the protocol: `AdapterProtocol`
- Test with simple molecules (aspirin, caffeine) before complex workflows
- Document API keys/config requirements in each adapter's README


  📊 Current Status: 65 Adapters

  🎯 Major Gaps Found:

  Missing Categories:

  1. 🌿 Natural Products - No databases (COCONUT has 400k compounds!)
  2. 🧬 Metabolomics - No metabolite databases (HMDB essential)
  3. ☠️ Toxicity-Specific - Only general ADMET (Tox21, CompTox needed)
  4. 💉 Antibodies/Biologics - Growing field, zero coverage
  5. 🧬 RNA Databases - Only protein/DNA covered
  6. 🧪 Forward Synthesis - Only retrosynthesis covered
  7. 🔗 Peptides - No peptide-specific tools
  8. ⚛️ Quantum Chemistry - No QM calculations
  9. 💊💊 Drug-Drug Interactions - Critical for safety
  10. 📊 Mass Spec Prediction - Useful for metabolite ID

  🚀 Recommended Next Batch (Get to 71 adapters)

  Batch 5: Natural Products & Metabolomics (2 adapters)

  1. COCONUT - 400k natural products, CC0 license, free API
  2. HMDB - Human Metabolome Database, essential for metabolomics

  Batch 6: Toxicity & Safety (2 adapters)

  3. Tox21 - NIH toxicity dataset, public domain
  4. CompTox - EPA Chemistry Dashboard, comprehensive tox data

  Batch 7: Biologics (2 adapters)

  5. SAbDab - Antibody structures database, free API
  6. ImmuneBuilder - Antibody structure prediction, pip install

  All 6 are:
  - ✅ Fully free
  - ✅ Commercially friendly licenses
  - ✅ Easy to implement (API or pip install)
  - ✅ Fill critical gaps
  - ✅ High user value

  ---
  💡 Quick Wins Beyond 71:

  Batch 8 (Optional, → 75 adapters):
  7. RNAcentral - RNA database (RNA-targeting drugs)
  8. ORD - Open Reaction Database (forward synthesis)
  9. IntAct - Protein-protein interactions (complement STRING)
  10. xTB - Fast semiempirical QM (pip install)

  ---
  ❓ What would you like to do?

  Option A: Build Batch 5+6+7 now (6 adapters → 71 total) ⭐ RECOMMENDED