# PharmForge Adapter Gap Analysis
**Date:** 2025-10-30
**Current Adapters:** 65
**Goal:** Identify missing categories and free/open-source opportunities

---

## Current Coverage (11 Categories)

### ✅ Well Covered
1. **Chemical Databases** (6): PubChem, ChEMBL, DrugCentral, ZINC, SureChEMBL, ChemSpider
2. **Protein Structures** (5): AlphaFold, RCSB PDB, PDBe, SWISS-MODEL, PDB-REDO
3. **Docking/MD** (6): Vina, GNINA, DiffDock, BindingDB, OpenMM, gmx_MMPBSA
4. **Literature** (4): PubMed, Europe PMC, Google Patents, Lens
5. **ML/Featurization** (18): Mordred, DeepChem, Chemprop, TorchDrug, etc.

### ⚠️ Partial Coverage
6. **ADMET** (6): ADMET-ai, pkCSM, SwissTarget, TargetNet, Oloren + RDKit
7. **Target/Disease** (5): OpenTargets, DisGeNET, UniProt, STRING, BioGRID
8. **Genomics** (4): GEO, GTEx, Reactome, KEGG

### 🔴 Minimal Coverage
9. **Clinical/Safety** (2): ClinicalTrials, FDA FAERS
10. **Retrosynthesis** (5): AiZynthFinder, LLM, REINVENT, MolGAN, DeNovo
11. **Enzymes** (1): BRENDA

---

## Missing Categories (High Value)

### 1. Natural Products 🌿
**Gap:** No natural product databases
**Free Options:**
- ✅ **COCONUT** (Natural products DB) - API/download
- ✅ **NP-MRD** (Marine natural products) - Web/API
- ✅ **StreptomeDB** (Streptomyces natural products) - Web
- ✅ **Super Natural II** - Download/API

**Recommendation:** Add COCONUT (largest, ~400k compounds)

---

### 2. Drug-Drug Interactions 💊💊
**Gap:** No DDI databases
**Free Options:**
- ✅ **DrugBank** (partial free) - Has DDI data
- ✅ **TWOSIDES** (side effects) - Open dataset
- ✅ **DrugCentral** (already have!) - Has DDI data, need to extend
- ⚠️ **Drugs.com API** - May have usage limits

**Recommendation:** Extend DrugCentral adapter with DDI queries

---

### 3. Metabolomics 🧬
**Gap:** No metabolite databases
**Free Options:**
- ✅ **HMDB** (Human Metabolome Database) - Free API
- ✅ **MetaboLights** (EMBL-EBI) - Free API
- ✅ **METLIN** (Metabolite database) - Free registration

**Recommendation:** Add HMDB (most comprehensive)

---

### 4. Toxicity Databases ☠️
**Gap:** Only general ADMET, no specific tox DBs
**Free Options:**
- ✅ **Tox21** (NIH) - Public dataset/API
- ✅ **CompTox** (EPA Chemistry Dashboard) - Free API
- ✅ **TOXNET/ChemIDplus** - Free
- ✅ **eTox** - Academic access

**Recommendation:** Add Tox21 + CompTox (2 adapters)

---

### 5. Pharmacophore & Fragment Libraries 🧩
**Gap:** No pharmacophore tools or fragment libraries
**Free Options:**
- ✅ **Pharmer** (pharmacophore search) - Open source
- ✅ **ZINC Fragments** (already have!)
- ✅ **ChEMBL Fragments** - Can query ChEMBL for fragments
- ✅ **PKIDB** (Kinase inhibitor fragments) - Free download

**Recommendation:** Add pharmacophore search to existing tools

---

### 6. RNA/Nucleotide Databases 🧬
**Gap:** Only protein databases, no RNA/DNA
**Free Options:**
- ✅ **RCSB PDB** (already have!) - Has RNA structures
- ✅ **RNAcentral** - Free API
- ✅ **miRBase** (microRNA) - Free
- ✅ **Rfam** (RNA families) - Free

**Recommendation:** Add RNAcentral for RNA-targeting drugs

---

### 7. Carbohydrate/Glycan Databases 🍬
**Gap:** No glycan/carbohydrate tools
**Free Options:**
- ✅ **GlyTouCan** (Glycan repository) - Free API
- ✅ **GlyConnect** - Free
- ✅ **CSDB** (Carbohydrate Structure DB) - Free

**Recommendation:** Add GlyTouCan if working on glycobiology

---

### 8. Mass Spectrometry Prediction 📊
**Gap:** No MS/NMR prediction tools
**Free Options:**
- ✅ **CFM-ID** (MS fragmentation) - Open source
- ✅ **MetFrag** - Open source
- ✅ **SIRIUS** - Free academic
- ✅ **nmrshiftdb2** (NMR prediction) - Free

**Recommendation:** Add CFM-ID for metabolite ID

---

### 9. Protein-Protein Interactions 🤝
**Gap:** Only STRING, could expand
**Free Options:**
- ✅ **STRING** (already have!)
- ✅ **BioGRID** (already have!)
- ✅ **IntAct** (EMBL-EBI) - Free API
- ✅ **MINT** - Free
- ✅ **DIP** - Free academic

**Recommendation:** Add IntAct (complementary to STRING)

---

### 10. Antibody/Biologics Tools 💉
**Gap:** No antibody-specific tools
**Free Options:**
- ✅ **SAbDab** (Antibody structures) - Free
- ✅ **IMGT** (Immunogenetics) - Free
- ✅ **ImmuneBuilder** (structure prediction) - Open source pip install
- ✅ **AbLang** (antibody language model) - Open source

**Recommendation:** Add SAbDab + ImmuneBuilder

---

### 11. Crystallization & Protein Design 💎
**Gap:** No crystallization or protein design tools
**Free Options:**
- ✅ **ProteinMPNN** (protein design) - Open source
- ✅ **ESMFold** (structure prediction) - Open source
- ✅ **OmegaFold** - Open source
- ⚠️ **RosettaFold** - Academic license

**Recommendation:** Add ProteinMPNN for protein design

---

### 12. Chemical Reactions & Synthesis 🧪
**Gap:** Only retrosynthesis, no forward synthesis
**Free Options:**
- ✅ **Reaxys** - ❌ Commercial
- ✅ **ORD** (Open Reaction Database) - Free
- ✅ **USPTO** (Patent reactions) - Free
- ✅ **RXN4Chemistry** (IBM) - Free API

**Recommendation:** Add ORD + RXN4Chemistry

---

### 13. Peptides & Cyclic Peptides 🔗
**Gap:** No peptide-specific tools
**Free Options:**
- ✅ **PepDB** (Peptide database) - Free
- ✅ **CyBase** (Cyclic peptides) - Free
- ✅ **BIOPEP** - Free
- ✅ **CAMP** (Antimicrobial peptides) - Free

**Recommendation:** Add CAMP for antimicrobial peptides

---

### 14. Quantum Chemistry (Lightweight) ⚛️
**Gap:** No QM calculations
**Free Options:**
- ✅ **Psi4** - Free (conda install)
- ✅ **PySCF** - pip install, BSD license
- ✅ **xTB** (semiempirical) - Free, very fast
- ⚠️ **ORCA** - Free for academic

**Recommendation:** Add xTB (fast semiempirical QM)

---

### 15. Lipophilicity & Membrane Transport 🧈
**Gap:** General ADMET but no transport-specific
**Free Options:**
- ✅ **MemSTAT** (membrane statistics) - In PDB
- ✅ **Transporter database** - Part of DrugBank
- ✅ **BBB prediction** - Part of pkCSM (already have!)

**Recommendation:** Extend existing ADMET adapters

---

## Priority Recommendations for 65 → 75 Adapters

### High Priority (Next 4 adapters → 69 total)
1. **COCONUT** - Natural products (unique niche)
2. **HMDB** - Human metabolome (metabolomics)
3. **Tox21** - Toxicity data (safety)
4. **CompTox** - EPA chemistry dashboard (tox/environmental)

### Medium Priority (Next 3 adapters → 72 total)
5. **ImmuneBuilder** - Antibody structure prediction
6. **SAbDab** - Antibody database
7. **RNAcentral** - RNA database

### Optional (Next 3 adapters → 75 total)
8. **ORD** - Open Reaction Database
9. **IntAct** - Protein interactions
10. **xTB** - Fast quantum chemistry

---

## Summary: What We're Missing

### By Impact:
1. 🌿 **Natural Products** - COCONUT
2. 🧬 **Metabolomics** - HMDB
3. ☠️ **Toxicity** - Tox21, CompTox
4. 💉 **Antibodies** - SAbDab, ImmuneBuilder
5. 🧬 **RNA** - RNAcentral
6. 🧪 **Reactions** - ORD, RXN4Chemistry
7. 🔗 **Peptides** - CAMP
8. ⚛️ **QM** - xTB

### By Ease of Implementation:
- **Easiest**: API-based (COCONUT, HMDB, Tox21, CompTox, SAbDab)
- **Medium**: pip install (ImmuneBuilder, xTB)
- **Harder**: Large datasets (ORD, RNAcentral)

---

## Recommended Next Batch (6 adapters to reach 71)

**Batch 5: Natural Products & Metabolomics**
1. COCONUT (natural products API)
2. HMDB (human metabolome API)

**Batch 6: Toxicity & Safety**
3. Tox21 (NIH toxicity dataset)
4. CompTox (EPA chemistry dashboard API)

**Batch 7: Biologics**
5. SAbDab (antibody structures API)
6. ImmuneBuilder (antibody prediction - pip install)

**Result:** 71 adapters with comprehensive coverage!

---

## License Check (All Recommended)

✅ **COCONUT** - CC0 (public domain)
✅ **HMDB** - Creative Commons (free for all)
✅ **Tox21** - Public domain (US govt)
✅ **CompTox** - Public domain (EPA)
✅ **SAbDab** - Academic free
✅ **ImmuneBuilder** - Open source (check license)
✅ **RNAcentral** - Open data
✅ **ORD** - CC-BY-SA
✅ **xTB** - LGPL (commercial OK)

All are commercially friendly! 🎉
