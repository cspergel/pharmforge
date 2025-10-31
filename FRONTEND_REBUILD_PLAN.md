# Frontend Rebuild Plan - Functional UX
**Date:** October 26, 2025
**Status:** CRITICAL - Current UI is non-functional
**Priority:** HIGH

---

## Problem Statement

Current frontend issues identified:
1. ❌ NL Run fails - no proper backend integration
2. ❌ No adapter selection - can't choose which tools to use
3. ❌ No compound testing - can't just test a SMILES string
4. ❌ No design mode - only evaluation, no generation
5. ❌ No benchmark visualization - can't see adapter performance
6. ❌ No adapter comparison - can't compare similar tools
7. ❌ No LLM summarization - no AI explanations
8. ❌ Mock data everywhere - not using real adapters

---

## User Requirements (from feedback)

### 1. Simple Compound Testing Interface
**What user wants:**
- Put in a SMILES string directly
- Search for compounds (PubChem/ChEMBL lookup)
- Store favorite compounds
- Pick which adapters to run (manual selection)
- Run tests and get structured output
- Export results easily

**Current state:** Doesn't exist
**Priority:** CRITICAL

### 2. Adapter Management & Visualization
**What user wants:**
- See all 39 adapters available
- View benchmarks for each adapter
- Compare adapters that do the same thing
- Test individual adapters easily
- Customize adapter parameters
- Let orchestrator auto-pick OR manual selection

**Current state:** Only shows 4 mock adapters
**Priority:** CRITICAL

### 3. Design Mode (not just Evaluation)
**What user wants:**
- Generate new molecules (not just test existing)
- Optimize molecules iteratively
- Evolution/genetic algorithms
- Structure-based design
- Scaffold hopping

**Current state:** Doesn't exist (only evaluation)
**Priority:** HIGH

### 4. LLM Result Summarization
**What user wants:**
- AI summary of results in plain English
- Explain what the scores mean
- Suggest next steps
- Flag concerns (toxicity, synthesis difficulty)

**Current state:** Doesn't exist
**Priority:** MEDIUM

### 5. Orchestrator Integration
**What user wants:**
Two modes:
- **Auto mode:** LLM picks best adapters for task
- **Manual mode:** User selects exact adapters to run

**Current state:** Only has broken auto mode
**Priority:** CRITICAL

---

## Implementation Plan

### Phase 1: Core Functionality (This Week)

#### 1.1 Compound Testing Page ⚡ CRITICAL
**File:** `frontend/pages/compound_testing.py`

**Features:**
- Direct SMILES input field
- Compound search (PubChem/ChEMBL API)
- Adapter multi-select checkboxes (all 39)
- Parameter customization per adapter
- Run button → structured results
- Export to JSON/CSV/PDF

**UI Layout:**
```
┌────────────────────────────────────────┐
│  Compound Input                        │
│  [SMILES] [Search] [From Library]     │
│                                        │
│  Select Adapters to Run:               │
│  ☐ Docking (3 adapters)               │
│  ☐ ADMET (2 adapters)                 │
│  ☐ Retrosynthesis (2 adapters)        │
│  ... (expand categories)               │
│                                        │
│  [Run Selected Adapters]               │
│                                        │
│  Results:                              │
│  ├─ Vina Docking: -8.2 kcal/mol       │
│  ├─ TDC ADMET: 0.87/1.0               │
│  └─ AiZynthFinder: 3 steps            │
│                                        │
│  [Export] [Save to Library] [AI Summary]│
└────────────────────────────────────────┘
```

#### 1.2 Real Adapter Integration ⚡ CRITICAL
**File:** `backend/api/adapters.py`

**New endpoints:**
```python
GET  /api/adapters/list           # Get all 39 adapters with metadata
GET  /api/adapters/{name}/info    # Get adapter details
POST /api/adapters/{name}/test    # Test single adapter
POST /api/adapters/batch          # Run multiple adapters
GET  /api/adapters/{name}/benchmarks  # Get benchmark results
```

**Response format:**
```json
{
  "adapters": [
    {
      "name": "vina_docking",
      "display_name": "AutoDock Vina",
      "category": "Docking & Scoring",
      "description": "Molecular docking with AutoDock Vina",
      "version": "1.2.3",
      "status": "healthy",
      "required_inputs": ["smiles", "target_protein"],
      "optional_params": {
        "exhaustiveness": {"type": "int", "default": 8},
        "num_modes": {"type": "int", "default": 9}
      },
      "output_schema": {
        "binding_affinity": "float",
        "pose_rmsd": "float"
      },
      "benchmarks": {
        "dud_e_auc": 0.82,
        "avg_runtime": "45s"
      }
    },
    ...
  ]
}
```

#### 1.3 Adapter Browser Page ⚡ CRITICAL
**File:** `frontend/pages/adapter_browser.py`

**Features:**
- Grid/list view of all adapters
- Filter by category, status, benchmark score
- Click adapter → see details + benchmarks
- Test adapter with sample compound
- Compare 2-3 adapters side-by-side

**UI Layout:**
```
┌────────────────────────────────────────┐
│  Filters: [Category ▼] [Status ▼]     │
│  Search: [________________]            │
│                                        │
│  📊 Docking & Scoring (3 adapters)     │
│  ┌─────────┬─────────┬─────────┐      │
│  │ Vina    │ Smina   │ DiffDock│      │
│  │ ✅ 0.82 │ ✅ 0.79 │ ✅ 0.91 │      │
│  │ [Test]  │ [Test]  │ [Test]  │      │
│  └─────────┴─────────┴─────────┘      │
│                                        │
│  💊 ADMET & Toxicity (2 adapters)     │
│  ...                                   │
└────────────────────────────────────────┘
```

### Phase 2: Advanced Features (Next Week)

#### 2.1 Design Mode Page
**File:** `frontend/pages/design_mode.py`

**Features:**
- **Generate Tab:**
  - Starting scaffold input
  - Generative models (REINVENT, ChemGPT)
  - Constraints (MW, LogP, etc.)
  - Generate N candidates

- **Optimize Tab:**
  - Input molecule
  - Objectives (binding, ADMET, synthesis)
  - Genetic algorithm parameters
  - Run evolution

- **Scaffold Hop Tab:**
  - Input molecule
  - Similarity threshold
  - Find bioisosteres

#### 2.2 LLM Result Summarizer
**File:** `backend/api/llm_summary.py`

**Endpoint:**
```python
POST /api/results/{run_id}/summarize
```

**Integration:**
- Use GPT-4 or Claude API
- Template: "Summarize these drug discovery results..."
- Returns: Plain English summary + recommendations

**Example output:**
```
Summary:
Your pipeline identified 15 promising EGFR inhibitors. The top candidate
(SMILES: ...) shows excellent binding affinity (-11.2 kcal/mol) and good
ADMET properties (score: 0.89), but synthesis may be challenging (7 steps).

Recommendations:
1. Prioritize compounds ranked 2-4 for better synthesis feasibility
2. Consider structural modifications to reduce hERG liability
3. Test top 5 in vitro for validation

Concerns:
- 3 compounds flagged for potential hepatotoxicity
- BBB penetration lower than target for 8 compounds
```

#### 2.3 Orchestrator Mode Selector
**Enhancement to New Run page**

**Add mode selector:**
```
Run Mode:
○ Auto (AI picks adapters)
● Manual (I choose adapters)
○ Preset Pipeline
```

**Manual mode UI:**
```
Pipeline Builder:
┌─────────────────────────┐
│ 1. Input Source         │
│    ☑ ChEMBL Search     │
├─────────────────────────┤
│ 2. Property Prediction  │
│    ☑ RDKit Properties  │
│    ☑ TDC ADMET         │
├─────────────────────────┤
│ 3. Docking              │
│    ☑ Vina Docking      │
│    ☐ Smina             │
├─────────────────────────┤
│ 4. Retrosynthesis       │
│    ☑ AiZynthFinder     │
├─────────────────────────┤
│ 5. Ranking              │
│    Weights: [customize]│
└─────────────────────────┘
```

### Phase 3: Polish & UX (Week After)

#### 3.1 Benchmark Visualization Dashboard
**Features:**
- ROC curves for docking adapters
- MAE/RMSE comparison for ADMET
- Runtime vs accuracy tradeoffs
- Historical performance tracking

#### 3.2 Compound Library Management
**Features:**
- Save favorite compounds
- Tag/categorize molecules
- Import from files
- Export collections

#### 3.3 Result Comparison View
**Features:**
- Compare multiple runs side-by-side
- Diff adapter outputs
- Statistical significance tests

---

## Technical Architecture

### Backend API Endpoints (New)

```python
# Adapter endpoints
GET    /api/adapters/list
GET    /api/adapters/{name}/info
POST   /api/adapters/{name}/execute
POST   /api/adapters/batch
GET    /api/adapters/{name}/benchmarks

# Compound endpoints
POST   /api/compounds/search           # PubChem/ChEMBL lookup
POST   /api/compounds/validate         # Validate SMILES
GET    /api/compounds/library          # User's saved compounds
POST   /api/compounds/library/add

# Run endpoints (enhanced)
POST   /api/runs/create/manual        # Manual adapter selection
POST   /api/runs/create/auto          # Orchestrator mode
GET    /api/runs/{id}/results
POST   /api/runs/{id}/summarize       # LLM summary

# Design endpoints
POST   /api/design/generate           # Generate new molecules
POST   /api/design/optimize           # Optimize existing
POST   /api/design/scaffold-hop       # Find bioisosteres
```

### Frontend Pages Structure

```
frontend/
├── streamlit_app.py           # Main navigation
├── pages/
│   ├── compound_testing.py    # ⚡ NEW - Simple testing
│   ├── adapter_browser.py     # ⚡ NEW - All adapters
│   ├── design_mode.py         # ⚡ NEW - Generation
│   ├── pipeline_builder.py    # ⚡ NEW - Manual mode
│   ├── benchmarks.py          # ⚡ NEW - Comparisons
│   └── library.py             # ⚡ NEW - Saved compounds
├── components/
│   ├── adapter_card.py        # ⚡ NEW - Adapter display
│   ├── compound_input.py      # ⚡ NEW - SMILES input
│   ├── llm_summary.py         # ⚡ NEW - AI summary
│   └── pipeline_visualizer.py # ⚡ NEW - Flow diagram
```

---

## Implementation Order (Critical Path)

### Week 1: Core Functionality
**Days 1-2:**
1. ✅ Create `/api/adapters/list` endpoint (backend)
2. ✅ Create `/api/adapters/{name}/execute` endpoint (backend)
3. ✅ Build Adapter Browser page (frontend)
4. ✅ Build Compound Testing page (frontend)

**Days 3-4:**
5. ✅ Add adapter selection to pipeline
6. ✅ Fix NL Run → Orchestrator integration
7. ✅ Add manual mode to New Run page
8. ✅ Test end-to-end compound → adapters → results

**Day 5:**
9. ✅ Polish UX
10. ✅ Add error handling
11. ✅ Write user guide

### Week 2: Advanced Features
**Days 1-2:**
12. ⏳ LLM summarization endpoint + UI
13. ⏳ Benchmark visualization dashboard
14. ⏳ Adapter comparison view

**Days 3-4:**
15. ⏳ Design mode (generate molecules)
16. ⏳ Compound library management
17. ⏳ Pipeline builder visual editor

**Day 5:**
18. ⏳ Integration testing
19. ⏳ Documentation updates
20. ⏳ User acceptance testing

---

## Success Criteria

### Must Have (Week 1)
- ✅ User can input SMILES and test any adapter
- ✅ User can see all 39 real adapters (not mocks)
- ✅ User can manually select adapters for a run
- ✅ Results display correctly with structured output
- ✅ Export works (JSON, CSV)

### Should Have (Week 2)
- ⏳ LLM summarization working
- ⏳ Benchmark comparisons visible
- ⏳ Design mode (generate) functional
- ⏳ Compound library saves/loads

### Nice to Have (Week 3)
- ⏳ Pipeline visual editor
- ⏳ Advanced filtering
- ⏳ Collaboration features

---

## Next Steps (Immediate)

**RIGHT NOW:**
1. Create backend `/api/adapters/list` endpoint
2. Create backend `/api/adapters/{name}/execute` endpoint
3. Build new Compound Testing page
4. Build new Adapter Browser page
5. Test with real adapters

**User can then:**
- Put in any SMILES
- Pick any adapters
- Get real results
- See all 39 adapters
- Compare adapter performance

---

**Created:** October 26, 2025
**Owner:** PharmForge Development Team
**Status:** Ready to implement
**Timeline:** 2 weeks to full functionality
