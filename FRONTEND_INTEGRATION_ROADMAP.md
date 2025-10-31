# PharmForge Frontend Integration Roadmap

**Date:** October 26, 2025
**Status:** Planning Phase
**Goal:** Integrate all 33 adapters into Streamlit/React frontend

---

## Current Status

### Backend Infrastructure ✅
- **23/33 adapters production-ready** (70%)
- **10/33 adapters need minor setup** (30%)
- OpenMM installed successfully
- RDKit, PyTorch, OpenBabel all working
- All ML adapters functional in fallback mode

### Enhanced Packages Status

| Package | Status | Notes |
|---------|--------|-------|
| **OpenMM** | ✅ Installed (8.3.1) | Molecular dynamics ready |
| **RDKit** | ✅ Installed (2023.09.6) | Core chemistry library |
| **PyTorch** | ✅ Installed (2.5.0) | ML framework |
| **OpenBabel** | ✅ Installed (3.1.1) | File conversion |
| **REINVENT** | ❌ Not on PyPI | Requires GitHub source install |

### REINVENT Installation Note
REINVENT is not available via pip. To install:
```bash
# Clone from GitHub
git clone https://github.com/MolecularAI/REINVENT4.git
cd REINVENT4
pip install -e .
```

However, **REINVENT adapter already works in RDKit fallback mode**, providing basic molecular generation functionality.

### Frontend Infrastructure

**Current:** Streamlit app exists at `frontend/streamlit_app.py`

**Design Spec:** Comprehensive React/Next.js design in `PharmForge_Frontend_Design_Spec.md`

**Gap:** Need to choose:
1. Enhance existing Streamlit app (faster MVP)
2. Build React/Next.js app per design spec (better UX, longer timeline)

---

## Adapter Integration Requirements

### Adapter Categories

#### 1. Core Drug Discovery (12 adapters)
**Ready:**
- UniProt, RCSB PDB, AlphaFold, SWISS-MODEL
- OpenTargets, BindingDB, DrugCentral, ZINC
- ChEMBL, PubChem, RDKit, ADMET-AI

**UI Components Needed:**
- Protein search widget
- Structure viewer (3D molecular visualization)
- Bioactivity data tables
- ADMET prediction cards

#### 2. Clinical & Safety (2 adapters)
**Ready:**
- ClinicalTrials.gov
- FDA FAERS

**UI Components Needed:**
- Clinical trial search
- Adverse event dashboard

#### 3. Literature & Patents (3 adapters)
**Ready:**
- PubMed, EuropePMC, SureChEMBL

**UI Components Needed:**
- Literature search interface
- Patent search widget
- Citation management

#### 4. ML/Computational (6 adapters)
**Ready:**
- REINVENT (fallback), MolGAN (fallback), De Novo (fragment)
- AiZynthFinder, Vina, OpenMM

**UI Components Needed:**
- Molecular generation interface
- Docking visualization
- MD simulation controls
- Retrosynthesis tree viewer

#### 5. Systems Biology (2 adapters)
**Ready:**
- Reactome, GTEx

**UI Components Needed:**
- Pathway visualization
- Expression heatmaps

#### 6. Needs Setup (8 adapters)
**Requires Configuration:**
- DiffDock (GPU + 5GB models)
- LLM Retrosynthesis (API keys)
- Lens.org, Google Patents (API keys)

---

## Integration Phases

### Phase 1: Core Infrastructure (Week 1)
**Goal:** Basic adapter connection to frontend

**Tasks:**
1. Create adapter registry endpoint (`/api/adapters`)
2. Create adapter status dashboard
3. Build unified search interface
4. Implement basic result display

**Deliverables:**
- List all 33 adapters in UI
- Show status (ready/needs-setup)
- Execute simple queries
- Display results in tables

**Estimated Time:** 3-5 days

---

### Phase 2: Essential Workflows (Week 2-3)
**Goal:** Enable core drug discovery workflows

**Tasks:**
1. **Protein Structure Workflow**
   - UniProt → AlphaFold → PDB
   - 3D structure viewer (Mol*, NGL Viewer)

2. **Compound Search Workflow**
   - PubChem/ChEMBL search
   - Property calculator (RDKit)
   - ADMET prediction (TDC)

3. **Docking Workflow**
   - Vina docking interface
   - Result visualization
   - Score interpretation

4. **Retrosynthesis Workflow**
   - AiZynthFinder interface
   - Route tree visualization
   - Stock availability checks

**Deliverables:**
- 4 complete workflows functional
- Progress tracking for long-running jobs
- Result export (CSV, JSON, SDF)

**Estimated Time:** 8-12 days

---

### Phase 3: Advanced Features (Week 4)
**Goal:** ML generation, batch processing

**Tasks:**
1. **ML Generation Interface**
   - REINVENT molecular generation
   - MolGAN/De Novo interfaces
   - Property-guided optimization

2. **Batch Processing**
   - CSV upload
   - Parallel execution
   - Progress dashboard

3. **Pipeline Builder**
   - Drag-and-drop interface
   - Save/load pipelines
   - Share pipelines

**Deliverables:**
- ML generation working
- Batch mode for 100+ compounds
- Custom pipeline builder

**Estimated Time:** 5-7 days

---

### Phase 4: Polish & Deploy (Week 5-6)
**Goal:** Production-ready UI

**Tasks:**
1. **Visualization**
   - 3D molecule viewer
   - Plotly charts
   - Heatmaps

2. **Documentation**
   - In-app tutorials
   - Example workflows
   - API documentation

3. **Performance**
   - Result caching
   - Lazy loading
   - WebSocket progress updates

**Deliverables:**
- Beautiful visualizations
- Complete documentation
- Performance optimized
- User testing complete

**Estimated Time:** 10-14 days

---

## Technology Stack Recommendation

### Option A: Streamlit MVP (Recommended for Speed)
**Pros:**
- Already started
- Python-native
- Rapid development
- Easy adapter integration

**Cons:**
- Limited customization
- Not as polished as React
- Performance limitations for complex UIs

**Timeline:** 3-4 weeks to production

### Option B: React/Next.js (Per Design Spec)
**Pros:**
- Full design spec ready
- Best user experience
- Highly customizable
- Better performance

**Cons:**
- Need to build frontend from scratch
- Longer development time
- Requires frontend expertise

**Timeline:** 6-8 weeks to production

### Recommendation
**Start with Streamlit, migrate to React later**

1. **Now (Weeks 1-4):** Build functional Streamlit app
2. **Later (Month 2):** Begin React migration
3. **Future:** Run both (Streamlit for power users, React for general users)

---

## Immediate Next Steps

### This Week (Days 1-3)

**Day 1: Adapter Registry API**
```python
# backend/api/adapters.py
@router.get("/adapters")
def list_adapters():
    return {
        "adapters": [
            {
                "name": "uniprot",
                "type": "database",
                "status": "ready",
                "description": "Protein sequence and annotation",
                "capabilities": ["protein_search", "sequence_fetch"]
            },
            # ... all 33 adapters
        ],
        "stats": {
            "total": 33,
            "ready": 23,
            "needs_setup": 10
        }
    }
```

**Day 2: Streamlit Dashboard**
```python
# frontend/streamlit_app.py
import streamlit as st
import requests

st.title("PharmForge - Drug Discovery Platform")

# Sidebar: Adapter Status
st.sidebar.header("Adapters")
response = requests.get("http://backend:8000/api/adapters")
adapters = response.json()

for adapter in adapters["adapters"]:
    status_icon = "✅" if adapter["status"] == "ready" else "⚠️"
    st.sidebar.write(f"{status_icon} {adapter['name']}")

# Main: Query Interface
query_type = st.selectbox("Select Workflow", [
    "Protein Structure Search",
    "Compound Property Prediction",
    "Molecular Docking",
    "Retrosynthesis Planning"
])

if query_type == "Protein Structure Search":
    protein_id = st.text_input("UniProt ID")
    if st.button("Search"):
        # Call UniProt → AlphaFold → PDB pipeline
        result = run_protein_pipeline(protein_id)
        display_protein_results(result)
```

**Day 3: Basic Workflows**
- Implement 2-3 simple workflows
- Test with real data
- Document usage

---

## API Endpoints Needed

### Adapter Management
```
GET  /api/adapters              # List all adapters
GET  /api/adapters/{name}       # Get adapter details
POST /api/adapters/{name}/test  # Test adapter connection
```

### Workflows
```
POST /api/workflows/protein     # Protein structure pipeline
POST /api/workflows/compound    # Compound search + properties
POST /api/workflows/docking     # Docking pipeline
POST /api/workflows/synthesis   # Retrosynthesis pipeline
```

### Jobs
```
POST /api/jobs                  # Submit job
GET  /api/jobs/{id}             # Get job status
GET  /api/jobs/{id}/results     # Get results
WS   /ws/jobs/{id}              # WebSocket progress updates
```

---

## Frontend Component Breakdown

### Core Components (Week 1)
1. **AdapterList** - Show all adapters with status
2. **SearchBar** - Unified search interface
3. **ResultsTable** - Display tabular results
4. **StatusBadge** - Show adapter/job status

### Workflow Components (Week 2-3)
5. **ProteinViewer** - 3D protein structures (Mol*)
6. **CompoundCard** - Molecule + properties display
7. **DockingResults** - Docking pose visualization
8. **SynthesisTree** - Retrosynthesis route tree
9. **ProgressTracker** - Job progress visualization

### Advanced Components (Week 4)
10. **MoleculeEditor** - SMILES/structure input
11. **PipelineBuilder** - Drag-and-drop workflow creator
12. **BatchUploader** - CSV/SDF file upload
13. **ExportPanel** - Result export options

---

## Data Flow Architecture

```
┌─────────────────────────────────────────────────┐
│         Streamlit Frontend (Port 8501)          │
│  • AdapterList                                  │
│  • Workflow Selector                            │
│  • Results Viewer                               │
└─────────────────┬───────────────────────────────┘
                  │ HTTP/WebSocket
┌─────────────────▼───────────────────────────────┐
│          FastAPI Backend (Port 8000)            │
│  • /api/adapters - Registry                     │
│  • /api/workflows - Pipelines                   │
│  • /api/jobs - Job Queue                        │
└─────────────────┬───────────────────────────────┘
                  │
┌─────────────────▼───────────────────────────────┐
│              Adapter Layer                      │
│  • 23 production-ready adapters                 │
│  • 10 setup-pending adapters                    │
└──────────────────────────────────────────────────┘
```

---

## Testing Strategy

### Week 1
- [ ] All 23 adapters accessible via API
- [ ] Adapter status dashboard functional
- [ ] Basic search works for 3 adapters

### Week 2-3
- [ ] Protein workflow: UniProt → AlphaFold
- [ ] Compound workflow: PubChem → RDKit → ADMET
- [ ] Docking workflow: Vina execution
- [ ] Retrosynthesis workflow: AiZynthFinder

### Week 4
- [ ] ML generation: 10 molecules generated
- [ ] Batch processing: 100 compounds
- [ ] Pipeline builder: Save/load works

### Week 5-6
- [ ] All visualizations working
- [ ] Performance: <3s response time
- [ ] Documentation complete
- [ ] 5 user tests passed

---

## Success Metrics

**Week 4 Checkpoint:**
- [ ] 4 workflows fully functional
- [ ] 23 adapters integrated
- [ ] 10+ example runs documented
- [ ] Ready for user testing

**Week 6 Checkpoint:**
- [ ] All 23 adapters integrated
- [ ] Beautiful visualizations
- [ ] Complete documentation
- [ ] 5 beta users onboarded

---

## Resources Needed

### Development Tools
- Streamlit (installed)
- Plotly (for charts)
- Mol*/NGL Viewer (3D molecules)
- React-Flow (pipeline builder, future)

### Documentation
- API endpoint documentation
- Workflow tutorials
- Video walkthroughs

### Testing
- Sample datasets
- Benchmark queries
- User feedback sessions

---

## Decision Points

### This Week
**Question:** Streamlit or React?
**Recommendation:** Start with Streamlit, plan React migration
**Reason:** Faster MVP, proven with adapters, easier Python integration

### Next Week
**Question:** Which workflows first?
**Recommendation:** Protein + Compound search
**Reason:** Highest user value, demonstrate core capabilities

### Week 4
**Question:** Continue Streamlit or start React?
**Recommendation:** Evaluate user feedback
**Reason:** Let usage data guide decision

---

## Conclusion

**Current Status:**
- ✅ Backend fully functional (23/33 adapters ready)
- ✅ OpenMM installed
- ✅ All core packages ready
- ⚠️ REINVENT works in fallback mode (no full install needed yet)
- 🔄 Frontend integration needed

**Recommended Path:**
1. **Week 1:** Basic adapter dashboard (Streamlit)
2. **Week 2-3:** Core workflows (protein, compound, docking, synthesis)
3. **Week 4:** ML generation + batch processing
4. **Week 5-6:** Polish, visualizations, documentation
5. **Month 2:** Consider React migration

**Timeline:** 4-6 weeks to fully functional Streamlit MVP

**Next Action:** Start with adapter registry API endpoint (Day 1 task above)

---

**Document Created:** October 26, 2025
**Next Review:** After Week 1 completion
**Owner:** Development Team

