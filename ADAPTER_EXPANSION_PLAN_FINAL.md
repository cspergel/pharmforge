# PharmForge Adapter Expansion Plan - FINAL
**Date:** 2025-10-30
**Current:** 59 adapters
**Target:** 65-75 adapters (6-16 new)
**Focus:** Fully free, unique capabilities, commercial-friendly licenses

---

## Current Status
✅ **59 adapters** (57 original + gmx_MMPBSA + PyMOL)
🎯 **Target:** 65-75 adapters
📊 **Need:** 6-16 more adapters

---

## Priority Batches (All Fully Free & Unique)

### Batch 2: Conformer & Advanced ML (2 adapters) - NEXT
**Target: 61 adapters**

1. **MolScrub** ⭐⭐⭐
   - Category: Conformer Generation & Cleaning
   - Install: `pip install molscrub`
   - License: LGPL (commercial OK with dynamic linking)
   - Unique: Automated conformer cleaning for docking
   - Status: NOT in awesome-drug-discovery (found separately)

2. **Oloren ChemEngine** ⭐⭐⭐
   - Category: Advanced ML Property Prediction
   - Install: `pip install olorenchemengine`
   - License: Open source (check specific)
   - Unique: State-of-the-art molecular property models
   - Status: IN awesome-drug-discovery

---

### Batch 3: AutoML & European DB (2 adapters)
**Target: 63 adapters**

3. **Auto-sklearn** ⭐⭐
   - Category: AutoML
   - Install: `pip install auto-sklearn`
   - License: BSD (commercial OK)
   - Unique: Complements TPOT with different approach
   - Status: IN awesome-drug-discovery

4. **PDBe (European PDB)** ⭐⭐⭐
   - Category: Protein Structure Database (European)
   - Install: API only (requests)
   - License: Free API
   - Unique: Alternative to RCSB, European data
   - Status: NOT in awesome-drug-discovery

---

### Batch 4: Chemical Databases (2 adapters)
**Target: 65 adapters** ✅ **MINIMUM GOAL**

5. **ChemSpider** ⭐⭐⭐
   - Category: Chemical Database (Royal Society of Chemistry)
   - Install: API (requests)
   - License: Free API (RSC)
   - Unique: Aggregates 100+ databases
   - Status: NOT in awesome-drug-discovery

6. **BRENDA Enzymes** ⭐⭐
   - Category: Enzyme Database
   - Install: API/web scraping
   - License: Free academic/commercial
   - Unique: Comprehensive enzyme kinetics data
   - Status: NOT in awesome-drug-discovery

---

### Batch 5: Advanced Analysis (2 adapters) - OPTIONAL
**Target: 67 adapters**

7. **QM-Descriptors (RDKit extension)** ⭐⭐
   - Category: Quantum Mechanical Descriptors
   - Install: Uses RDKit + small extensions
   - License: BSD (RDKit)
   - Unique: Quantum chemistry features without full QM
   - Status: RDKit feature

8. **COCONUT (Natural Products)** ⭐⭐⭐
   - Category: Natural Product Database
   - Install: API (requests)
   - License: Free/Open
   - Unique: Largest natural products database
   - Status: NOT in awesome-drug-discovery

---

### Batch 6: Specialized Tools (2 adapters) - OPTIONAL
**Target: 69 adapters**

9. **ChemBERTA** ⭐⭐⭐
   - Category: Transformer-based molecular embeddings
   - Install: `pip install transformers`
   - License: Apache 2.0 (commercial OK)
   - Unique: Modern NLP-style embeddings for molecules
   - Status: IN awesome-drug-discovery (under transformers)

10. **ImmuneBuilder** ⭐⭐
    - Category: Antibody Structure Prediction
    - Install: `pip install immunebuilder`
    - License: Open source
    - Unique: Specialized for antibodies (complements AlphaFold)
    - Status: NOT in awesome-drug-discovery

---

### Batch 7: Pathway & Systems (2 adapters) - OPTIONAL
**Target: 71 adapters**

11. **WikiPathways** ⭐⭐⭐
    - Category: Biological Pathways Database
    - Install: API (requests)
    - License: Free/Open
    - Unique: Community-curated pathways
    - Status: NOT in awesome-drug-discovery

12. **STITCH (Protein-Chemical)** ⭐⭐
    - Category: Protein-Chemical Interaction Database
    - Install: API (requests)
    - License: Free
    - Unique: Links chemicals to proteins
    - Status: NOT in awesome-drug-discovery

---

### Batch 8: Visualization & Format (2 adapters) - OPTIONAL
**Target: 73 adapters**

13. **ChemDoodle Web** ⭐
    - Category: Web-based molecular drawing
    - Install: API/JS library
    - License: GPLv3 (commercial license available)
    - Unique: Interactive drawing in web
    - Status: NOT in awesome-drug-discovery
    - NOTE: May need commercial license check

14. **Open Babel Python** (enhanced) ⭐⭐
    - Category: Format conversion (enhanced beyond current)
    - Install: Already have, extend features
    - License: GPL (commercial OK)
    - Unique: Add more format conversions
    - Status: Already exists, can enhance

---

### Batch 9: Toxicity & Safety (2 adapters) - OPTIONAL
**Target: 75 adapters** ✅ **MAXIMUM GOAL**

15. **Tox21 Data** ⭐⭐⭐
    - Category: Toxicity Database
    - Install: API/dataset download
    - License: Public domain
    - Unique: NIH toxicity screening data
    - Status: IN awesome-drug-discovery (as dataset)

16. **CompTox Chemistry Dashboard** ⭐⭐⭐
    - Category: Toxicity & Chemistry Database (EPA)
    - Install: API (requests)
    - License: Free (US EPA)
    - Unique: Comprehensive toxicity data
    - Status: NOT in awesome-drug-discovery

---

## License Analysis (Commercial Friendliness)

### ✅ Definitely Commercial-Friendly
- **MIT/BSD/Apache 2.0**: Oloren, PyMOL, Auto-sklearn, ChemBERTA
- **Public APIs**: ChemSpider, PDBe, BRENDA, COCONUT, WikiPathways, Tox21, CompTox
- **LGPL (dynamic linking)**: MolScrub, OpenBabel

### ⚠️ Need Verification
- **GPL/GPLv3**: May require source disclosure or commercial license
- **Custom licenses**: Check each tool's specific terms

### ❌ Skip (Not Free/Commercial)
- SwissADME (web only, rate limited)
- Proprietary tools requiring licenses

---

## Implementation Priority

### TODAY (6-8 hours)
1. ✅ gmx_MMPBSA (DONE)
2. ✅ PyMOL (DONE)
3. **MolScrub** (1 hour)
4. **Oloren ChemEngine** (1.5 hours)
5. **Auto-sklearn** (1 hour)
6. **PDBe** (1 hour)

**Result: 63 adapters**

### TOMORROW (4-6 hours)
7. **ChemSpider** (1.5 hours)
8. **BRENDA** (1.5 hours)
9. **COCONUT** (1 hour)
10. **ChemBERTA** (1.5 hours)

**Result: 67 adapters**

### IF TIME PERMITS
11-16: WikiPathways, STITCH, ImmuneBuilder, Tox21, CompTox, QM-Descriptors

**Max: 73-75 adapters**

---

## Exclusions

### Too Complex / Not Free
❌ AutoDock-GPU (compilation)
❌ GROMACS (system install)
❌ Smina (compilation)
❌ VMD (proprietary binary)
❌ ChemAxon (commercial license)
❌ Schrodinger (commercial)

### Duplicates / Already Covered
❌ PubChem (already have)
❌ ChEMBL (already have)
❌ ZINC (have zinc_fragments)
❌ RDKit (have rdkit_local)

---

## Success Metrics

- [ ] 65+ adapters (minimum goal)
- [ ] All fully free with no API keys required (optional is OK)
- [ ] Unique capabilities (no duplicates)
- [ ] Commercial-friendly licenses verified
- [ ] All registered and tested
- [ ] Documentation complete

---

## Next Steps

**NOW**: Start Batch 2 (MolScrub + Oloren)
