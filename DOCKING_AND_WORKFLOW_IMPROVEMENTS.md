# Docking Tools & Workflow Improvements Plan

## 🎯 Overall Goal

Transform PharmForge from a **compound screening tool** into a comprehensive **drug discovery platform** with full docking capabilities and improved workflow.

---

## 📦 Phase 1: Install Docking Tools ⏳ IN PROGRESS

### What's Being Installed:

**Molecular Dynamics:**
- ✅ OpenMM >= 8.1.0 - GPU-accelerated molecular dynamics
- ✅ PDBFixer >= 1.9 - Protein structure cleanup and preparation

**Docking:**
- ✅ Vina >= 1.2.5 - AutoDock Vina Python bindings (fast, reliable docking)
- ✅ PyTorch Geometric - For DiffDock (AI-powered docking)
  - torch-scatter
  - torch-sparse
  - torch-cluster

**Build Status:** Running in background (10-20 minutes)

### Files Modified:
1. `requirements.txt` - Added openmm, pdbfixer, vina
2. `Dockerfile.backend` - Added PyTorch Geometric installation

---

## 🎨 Phase 2: Add Protein Input to UI

### Current State:
```
SMILES Input Only
↓
Select Adapters
↓
Run
```

### New State:
```
SMILES Input
↓
Protein Input (conditional - shows for docking adapters)
  - Option A: Upload PDB file
  - Option B: Enter PDB ID (fetch from RCSB)
  - Option C: Enter UniProt ID (fetch + predict with AlphaFold)
↓
Select Adapters
↓
Run (passes both compound + protein to docking adapters)
```

### UI Changes Needed:

**1. Add Protein Input Section to Compound Testing:**
```python
# After SMILES input, before adapter selection
if any_docking_adapter_selected():
    st.divider()
    st.markdown("**Protein Target (Required for Docking)**")

    protein_input_mode = st.radio(
        "How would you like to provide the protein?",
        ["Upload PDB File", "PDB ID", "UniProt ID"],
        horizontal=True
    )

    if protein_input_mode == "Upload PDB File":
        protein_file = st.file_uploader("Upload PDB file", type=['pdb'])

    elif protein_input_mode == "PDB ID":
        pdb_id = st.text_input("PDB ID", placeholder="1HSG")
        if pdb_id:
            # Fetch from RCSB automatically
            protein_data = fetch_pdb(pdb_id)

    elif protein_input_mode == "UniProt ID":
        uniprot_id = st.text_input("UniProt ID", placeholder="P04637")
        if uniprot_id:
            # Fetch sequence + AlphaFold structure
            protein_data = fetch_alphafold_structure(uniprot_id)
```

**2. Conditional Display:**
- Show protein section ONLY when docking adapter is selected
- Clear, helpful messaging
- Auto-fetch from databases when possible

**3. Adapter Categorization:**
```python
DOCKING_ADAPTERS = {
    'vina', 'gnina', 'diffdock', 'smina'
}

MD_ADAPTERS = {
    'openmm'
}

STRUCTURE_ADAPTERS = {
    'alphafold', 'swissmodel', 'rcsb_pdb'
}
```

---

## 🔧 Phase 3: Connect Protein Data to Adapters

### Current Adapter Execution:
```python
response = api_client.test_adapter(adapter_name, smiles)
```

### New Adapter Execution:
```python
# For docking adapters
if adapter_name in DOCKING_ADAPTERS:
    response = api_client.test_adapter(
        adapter_name,
        smiles,
        protein_data=protein_data  # Pass protein
    )
else:
    # Non-docking adapters don't need protein
    response = api_client.test_adapter(adapter_name, smiles)
```

### Backend Changes:

**1. Update API Client (`frontend/components/api_client.py`):**
```python
def test_adapter(
    self,
    adapter_name: str,
    smiles: str,
    protein_data: Optional[Dict] = None
) -> AdapterResponse:
    data = {"input_data": {"smiles": smiles}}

    if protein_data:
        data["input_data"]["protein_data"] = protein_data

    return self._make_request(
        "POST",
        f"/api/adapters/{adapter_name}/execute",
        data=data
    )
```

**2. Update Adapter Implementations:**

Docking adapters will accept protein_data in their execute() method:
```python
async def execute(self, smiles: str, protein_data: Dict = None, **params):
    if not protein_data:
        return AdapterResult(
            success=False,
            error="Protein data required for docking"
        )

    # Run docking with protein
    results = self.dock(smiles, protein_data)
    ...
```

---

## 📊 Phase 4: Improve Visualizations & Results

### Docking Results Display:

**Current:**
```json
{
  "binding_affinity": -8.5,
  "pose_data": "..."
}
```

**Improved:**
```
Binding Affinity: -8.5 kcal/mol ⭐⭐⭐ (Excellent)
RMSD: 1.2 Å
Key Interactions:
  - TYR123 (H-bond)
  - ASP456 (Salt bridge)
  - LEU789 (Hydrophobic)

[3D Visualization of binding pose]
[Download PDB]
```

### Comparison View:

When testing multiple compounds against same protein:
```
Compound A: -8.5 kcal/mol ⭐⭐⭐
Compound B: -7.2 kcal/mol ⭐⭐
Compound C: -6.1 kcal/mol ⭐

[Scatter plot: ADMET Score vs Binding Affinity]
[Pareto frontier: binding + ADMET + synthesis]
```

---

## 🎯 Phase 5: Workflow Enhancements

### Current Workflow:
1. Enter compound
2. Select adapters
3. Run
4. View results

### Enhanced Workflows:

**Workflow A: Target-Based Discovery**
```
1. Select Target Protein (PDB/UniProt)
2. Upload/Generate Compound Library
3. Auto-select relevant adapters:
   - Docking (Vina, DiffDock)
   - ADMET (ADMET AI)
   - Novelty (ChEMBL)
4. Run parallel screening
5. Results table with multi-objective scoring
6. Export top candidates
```

**Workflow B: ADMET + Docking Pipeline**
```
1. SMILES input
2. Run ADMET AI (filter for druglike properties)
3. If passes → Run docking
4. If passes → Run synthesis check
5. Output: Compounds that are druglike, bind well, and synthesizable
```

**Workflow C: Comparative Analysis**
```
1. Upload 10-100 compounds (CSV)
2. Select target protein
3. Run comprehensive screening:
   - Parallel docking
   - ADMET profiling
   - Database searches
4. Rank by composite score
5. Generate comparison report
```

### UI for Workflows:

```
🏠 Home
🧪 Compound Testing ✅ (current - single compound)
🔬 Batch Screening 🆕 (upload CSV, screen library)
🎯 Workflow Builder 🆕 (drag-drop adapter chains)
📊 Results Dashboard 🆕 (visualize, compare, export)
```

---

## 🎨 Design Improvements

### Current Issues:
- ❌ All adapters shown at once (overwhelming)
- ❌ No categorization
- ❌ Unclear which adapters need additional inputs
- ❌ No saved workflows/templates

### Proposed Improvements:

**1. Adapter Categorization with Expandable Groups:**
```
▼ ADMET & Properties (6 adapters)
  ✓ ADMET AI - Comprehensive ADMET predictions 🔥
  □ RDKit Local - Basic molecular properties
  □ pkCSM - PK/safety predictions
  ...

▼ Docking & Binding (4 adapters)
  □ AutoDock Vina - Fast, reliable docking ⚡
  □ GNINA - CNN scoring docking
  □ DiffDock - AI-powered docking 🤖
  □ SMINA - Minimization docking
  [Requires protein target 🧬]

▼ Database Search (8 adapters)
  □ PubChem - Chemical database
  □ ChEMBL - Bioactivity database
  ...

▼ Literature & Patents (4 adapters)
  ...

▼ Protein Structure (3 adapters)
  ...

▼ Synthesis Planning (2 adapters)
  ...
```

**2. Smart Templates:**
```
Quick Start Templates:
📋 Basic Screening (ADMET + Databases)
🎯 Target-Based (Docking + ADMET + Synthesis)
🔬 Comprehensive (All relevant adapters)
💊 Drug Repurposing (Databases + Literature)

[Create Custom Template]
```

**3. Adapter Cards with Badges:**
```
┌────────────────────────────────┐
│ 🧪 ADMET AI               🔥   │
│ ML-based ADMET predictions     │
│ ⚡ Fast · 🎯 GPU · ⭐⭐⭐      │
│ 49 properties · 1-3s runtime   │
└────────────────────────────────┘

┌────────────────────────────────┐
│ 🎯 AutoDock Vina          ⚡   │
│ Protein-ligand docking         │
│ 🧬 Needs protein · ⭐⭐⭐      │
│ Gold standard · 10-30s         │
└────────────────────────────────┘
```

**4. Progress Visualization:**
```
Running 8 adapters in parallel...

✅ RDKit Local (0.5s)        █████████ DONE
✅ PubChem (2.1s)            █████████ DONE
⏳ ADMET AI (8/10s)          ████████░ 80%
⏳ Vina Docking (12/30s)     ████░░░░░ 40%
🔄 ChEMBL (querying...)      ░░░░░░░░░
⏸️  DiffDock (waiting...)    ░░░░░░░░░
```

---

## 📈 Success Metrics

**Before:**
- Single compound testing only
- No docking capabilities
- Sequential execution
- Limited to property prediction

**After:**
- ✅ Full docking workflow (compound + protein)
- ✅ 4 docking methods (Vina, GNINA, DiffDock, SMINA)
- ✅ Parallel execution (3-5x faster)
- ✅ Batch screening (CSV upload)
- ✅ Workflow templates
- ✅ Improved UX with categorization
- ✅ Multi-objective optimization
- ✅ Export & visualization

**Value Increase:**
- Current: 6/10
- With docking: 8/10
- With workflows: 9/10
- With batch + viz: 9.5/10

---

## 🚀 Implementation Timeline

**Phase 1** (IN PROGRESS): Install docking tools - 10-20 min
**Phase 2** (NEXT): Add protein input UI - 1-2 hours
**Phase 3**: Connect protein data to adapters - 1 hour
**Phase 4**: Improve visualizations - 2-3 hours
**Phase 5**: Workflow enhancements - 3-4 hours

**Total**: ~8-12 hours for complete transformation

---

## 📝 Next Actions

1. ⏳ Wait for build to complete (~10-20 min)
2. ✅ Verify docking tools installed
3. 🎨 Add protein input to compound testing UI
4. 🔧 Update adapter execution logic
5. 🧪 Test end-to-end docking workflow
6. 📊 Add result visualizations
7. 🎯 Build workflow templates

---

**Current Status:** Building docking tools in background...

Check build progress: `docker logs -f <container_id>`
