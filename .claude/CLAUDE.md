# 🛠 Update Notes v2 — Runtime Fixes Baked In

**What changed (drop‑in safe):**
1) **Score normalization unified:** All objectives must be **0–1, higher=better** before ranking.
   - Docking (Vina): convert kcal/mol with `vina_affinity_to01()` (−12 → 1.0, −4 → 0.0).
   - Retrosynthesis steps: map fewer steps → higher score via `synthesis_steps_to01()`.
2) **Health endpoint status:** `/health` now returns **503** when DB/Redis fail; 200 when healthy.
3) **Background tasks:** Async pipelines wrapped so FastAPI `BackgroundTasks` actually execute them.
4) **Adapter registry hygiene:** Hard checks to ensure required adapters are registered by exact name.
5) **Frontend API URL:** UI reads `PHARMFORGE_API_URL` env, defaulting to `http://backend:8000`.
6) **Terraform completeness:** Notes added for ALB listener, target groups, security groups, and unique S3 names.
7) **Tests added:** Docking/retrosynthesis normalization tests, pipeline smoke, and health degradation path.

**Drop‑in helper snippets (reference):**

```python
# backend/core/scoring_utils.py
def to01(val: float, lo: float, hi: float) -> float:
    if lo == hi: return 0.0
    v = max(min(val, hi), lo)
    return (v - lo) / (hi - lo)

def vina_affinity_to01(kcal: float) -> float:   # more negative is better
    return to01(-kcal, 4.0, 12.0)               # maps [-12,-4] -> [1,0]

def synthesis_steps_to01(steps: int) -> float:  # fewer is better
    steps = max(1, min(int(steps), 10))
    return to01(11 - steps, 1.0, 10.0)
```

```python
# backend/app/health.py
from fastapi import APIRouter, Response
import json
router = APIRouter()
@router.get("/health")
def health(db_ok: bool = True, redis_ok: bool = True):
    checks = {"db": db_ok, "redis": redis_ok}
    healthy = all(checks.values())
    code = 200 if healthy else 503
    return Response(
        content=json.dumps({"status": "ok" if healthy else "degraded", **checks}),
        media_type="application/json",
        status_code=code,
    )
```

```python
# backend/app/run_tasks.py
import asyncio
def run_async_in_background(coro_fn, *args, **kwargs):
    def _runner():
        asyncio.run(coro_fn(*args, **kwargs))
    return _runner
# usage: background_tasks.add_task(run_async_in_background, pipeline.execute, smiles, run_id)
```

```python
# Registry name assertions
required = ["vina_docking", "tdc_admet", "aizynthfinder"]
for name in required:
    assert registry.get(name) is not None, f"Adapter not registered: {name}"
```

```python
# Frontend API URL
// Streamlit or React should read from env, e.g.:
API_URL = os.getenv("PHARMFORGE_API_URL", "http://backend:8000")
```


# PharmForge Build Guide - Executive Summary

**Target:** AI Coding Agent (Claude Code)  
**Timeline:** 12 weeks (90 days)  
**Team Size:** 2 founders working full-time  
**Goal:** Launch open-source MVP + cloud beta with first 10 paying customers

---

## 🎯 Project Overview

**PharmForge** is an open-source drug discovery workflow orchestrator that connects 10+ public APIs (PubChem, ChEMBL, OpenTargets) with local compute tools (docking, ADMET prediction, retrosynthesis) into complete *in silico* pipelines.

**Core Value Proposition:**
- Natural language → automated pipeline → ranked drug candidates
- Open-source core (MIT license)
- Web-first deployment (React + FastAPI)
- Three modes: Natural Language, Batch Processing, Evolution

---

## 🏗️ Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    User Interface Layer                      │
│  React/Next.js Web App + Streamlit (MVP) + CLI (optional)   │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│                  Control Plane (FastAPI)                     │
│  • Authentication (JWT)                                      │
│  • Arcana Orchestrator (NL → Pipeline)                      │
│  • Job Queue Management (Celery + Redis)                    │
│  • Results Storage (PostgreSQL)                              │
└─────────────────────┬───────────────────────────────────────┘
                      │
┌─────────────────────▼───────────────────────────────────────┐
│                    Adapter Layer                             │
│  • Public APIs: PubChem, ChEMBL, OpenTargets                │
│  • Local Compute: DiffDock, TDC ADMET, AiZynthFinder        │
│  • Caching: Redis (hot) + Disk (warm) + S3 (cold)          │
└──────────────────────────────────────────────────────────────┘
```

**Key Technical Decisions:**
- **Backend:** FastAPI (async, OpenAPI auto-gen, Pydantic validation)
- **Queue:** Celery + Redis (async task execution, progress tracking)
- **Database:** PostgreSQL (JSONB for flexible schema)
- **Caching:** Redis + disk for reproducibility
- **Container:** Docker + Docker Compose for development
- **Deployment:** AWS/GCP (cloud), self-hosted option

---

## 📁 Repository Structure

```
pharmforge/
├── backend/                 # FastAPI application
│   ├── api/                 # REST endpoints
│   ├── core/                # Business logic
│   │   ├── orchestrator.py  # Arcana NL planning
│   │   ├── pipeline.py      # Pipeline execution engine
│   │   └── adapters/        # Adapter implementations
│   ├── models/              # Pydantic schemas
│   ├── db/                  # Database models (SQLAlchemy)
│   └── tests/               # Backend tests
├── frontend/                # React/Next.js or Streamlit
│   ├── src/                 # React components (or streamlit_app.py)
│   ├── public/              # Static assets
│   └── tests/               # Frontend tests
├── adapters/                # Adapter plugin system
│   ├── pubchem/
│   ├── chembl/
│   ├── rdkit_local/
│   ├── tdc_admet/
│   ├── diffdock/
│   └── askcos/
├── config/                  # Configuration files
│   ├── presets.yaml         # Pipeline presets
│   ├── adapters.yaml        # Adapter registry
│   └── docker-compose.yml   # Local dev environment
├── docs/                    # Documentation
│   ├── tutorials/
│   ├── api/
│   └── architecture/
├── scripts/                 # Utility scripts
│   ├── setup.sh
│   └── benchmark.py
└── tests/                   # Integration tests
    ├── e2e/
    └── benchmarks/
```

---

## 🗓️ Phased Build Plan

### **Phase 1: Core Infrastructure (Weeks 1-4)**
**Focus:** Backend foundation, adapter system, first real adapters

**Deliverables:**
- ✅ Docker dev environment (`docker-compose up` works)
- ✅ FastAPI with health checks, auth middleware
- ✅ PostgreSQL schema + Redis caching
- ✅ Adapter protocol implementation
- ✅ Working adapters: PubChem, ChEMBL, RDKit, TDC ADMET, Vina
- ✅ Celery job queue with progress tracking
- ✅ Caching layer (deterministic keys)

**Success Criteria:**
- Query 100 compounds through full property + bioactivity + ADMET pipeline
- Results cached and reproducible
- Job queue handles concurrent requests

📄 **Build Document:** `phase1_weeks1-4.md`

---

### **Phase 2: Pipeline Completion (Weeks 5-8)**
**Focus:** Retrosynthesis, ranking, NL orchestration, UI

**Deliverables:**
- ✅ AiZynthFinder retrosynthesis adapter
- ✅ Multi-objective ranking (Pareto + weighted)
- ✅ Arcana orchestrator (GPT-4 integration)
- ✅ Streamlit UI with run wizard + results visualization
- ✅ Lockfile generation + export formats
- ✅ Documentation (README, tutorials, API docs)

**Success Criteria:**
- Natural language query → JSON pipeline plan → ranked results
- UI usable by non-technical user
- Lockfile enables bit-for-bit reproduction
- 80%+ test coverage on core modules

📄 **Build Document:** `phase2_weeks5-8.md`

---

### **Phase 3: Polish & Launch (Weeks 9-12)**
**Focus:** Preprint, cloud deployment, community launch

**Deliverables:**
- ✅ Validation benchmarks (DUD-E, TDC)
- ✅ Preprint submission (ChemRxiv)
- ✅ AWS cloud infrastructure
- ✅ Beta signup flow + Stripe billing
- ✅ Content marketing (3 blog posts)
- ✅ GitHub launch (public repo)
- ✅ First 10-20 beta users onboarded

**Success Criteria:**
- Break-even at ~35 customers ($8.7k MRR)
- 500+ GitHub stars
- 3 academic partnerships signed
- Published validation data

📄 **Build Document:** `phase3_weeks9-12.md`

---

## 🔧 Development Workflow

### **Daily Workflow**
```bash
# 1. Start services
docker-compose up -d

# 2. Run backend tests
cd backend && pytest

# 3. Start development server
uvicorn backend.main:app --reload

# 4. Start frontend
cd frontend && npm run dev  # or streamlit run app.py

# 5. Run integration test
python tests/e2e/test_full_pipeline.py
```

### **Git Workflow**
- Main branch: `main` (always deployable)
- Feature branches: `feature/adapter-pubchem`
- Weekly tags: `v0.1.0`, `v0.2.0`, etc.
- Daily commits even if incomplete

### **Testing Strategy**
- **Unit tests:** Each adapter, core logic module
- **Integration tests:** Full pipeline runs
- **End-to-end tests:** UI → backend → results
- **Benchmarks:** DUD-E, TDC (run weekly)

---

## 🎯 Success Metrics (Day 90)

| Metric | Target | Measurement |
|--------|--------|-------------|
| **Product** | MVP deployed (local + cloud) | `docker-compose up` works; cloud beta live |
| **GitHub** | 500 stars | Star count + 50 unique clones |
| **Users** | 50 active (10 paying) | MAU tracked in PostHog/Mixpanel |
| **MRR** | $2k-$5k | Stripe revenue |
| **Academic Partners** | 3 commitments | Signed DUAs |
| **Community Adapters** | 10 total | 5 team + 5 external |
| **Documentation** | Complete | README, API docs, 5 tutorials |
| **Validation** | Preprint submitted | ChemRxiv submission |

---

## 🚨 Critical Decision Points

### **Week 4 Checkpoint**
**Question:** Can we handle 100 compounds through full pipeline?  
**If No:** Simplify ADMET/docking; defer retrosynthesis  
**If Yes:** Proceed to Phase 2

### **Week 8 Checkpoint**
**Question:** Is Arcana NL planning accurate for common queries?  
**If No:** Ship with template-based planning; improve NL post-launch  
**If Yes:** Proceed to Phase 3

### **Week 12 Checkpoint**
**Question:** Do we have 3+ beta users providing feedback?  
**If No:** Delay public launch; focus on user recruitment  
**If Yes:** Execute GitHub launch

---

## 🔥 Failure Mode Prevention

### **Scope Creep Protection**
❌ **Don't Build:**
- Custom ML model training (use TDC, published models)
- Perfect UI (Streamlit MVP is enough)
- Advanced visualization (Plotly basics only)
- Social features (comments, sharing)
- Mobile app (web-responsive is enough)

✅ **Do Build:**
- Working end-to-end pipeline
- Reproducible results (lockfiles)
- Clear documentation
- Real user validation

### **Technical Debt Management**
- Reserve **20% of each week** for refactoring
- Write tests **before** complexity grows
- Document architectural decisions (ADRs)
- Review code quality weekly

### **Burnout Prevention**
- Take **one full day off per week**
- Cap work at **10 hours/day max**
- Switch tasks when stuck (don't grind)
- Celebrate small wins daily

---

## 📚 Key Resources

### **Documentation**
- [Master Documentation](PharmForge_Master_v4.6.2.md) - Full technical spec
- [90-Day Plan](PharmForge_90Day_Execution_Plan.md) - Week-by-week breakdown
- [Investor Q&A](PharmForge_Investor_QA_Supplement.md) - Business context
- [Deep Dive Q&A](PharmForge_Appendix_A11_DeepDive_QA.md) - Academic partnerships, unit economics

### **External APIs**
- PubChem: <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>
- ChEMBL: <https://chembl.gitbook.io/chembl-interface-documentation/web-services>
- OpenTargets: <https://platform-docs.opentargets.org/>
- TDC: <https://tdcommons.ai/> (Therapeutics Data Commons)

### **Tools & Libraries**
- FastAPI: <https://fastapi.tiangolo.com/>
- Celery: <https://docs.celeryproject.org/>
- RDKit: <https://www.rdkit.org/docs/>
- DiffDock: <https://github.com/gcorso/DiffDock>
- AiZynthFinder: <https://github.com/MolecularAI/aizynthfinder>

---

## 🤖 AI Agent Instructions

**When building PharmForge, you should:**

1. **Read phase documents sequentially** - Don't skip ahead
2. **Test after every major step** - Verify functionality before moving on
3. **Commit frequently** - Even incomplete work (mark with TODO)
4. **Ask for clarification** - If requirements are ambiguous
5. **Document decisions** - Create ADRs for non-obvious choices
6. **Flag blockers early** - Don't spend >2 hours stuck

**When encountering issues:**
- Check the troubleshooting section in each phase document
- Review existing adapters for patterns
- Search documentation for similar problems
- Ask human for guidance on business logic decisions

**Code quality standards:**
- Type hints for all Python functions
- Docstrings for public APIs
- Unit tests for core logic
- Integration tests for pipelines
- Error handling with specific exceptions
- Logging at appropriate levels

---

## 📞 Communication Protocol

**For Claude Code Agent:**
- Use `# TODO: [BLOCKED]` for blockers requiring human input
- Use `# DECISION NEEDED:` for business logic questions
- Use `# PERFORMANCE:` for optimization concerns
- Commit messages: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`

**Example blockers:**
```python
# TODO: [BLOCKED] Need ASKCOS API token
# DECISION NEEDED: Should we filter compounds with MW > 500 by default?
# PERFORMANCE: This loop takes 10min for 1000 compounds - optimize?
```

---

## 🎬 Getting Started

**For AI Agent (Claude Code):**

1. **Read this document first** - Understand project context
2. **Review Phase 1 document** - Week 1-4 detailed instructions
3. **Set up development environment** - Follow setup steps
4. **Start with Week 1, Day 1** - Project setup + tooling
5. **Execute step-by-step** - Don't skip validation steps

**Ready to begin?** → Open `phase1_weeks1-4.md`

---

**Version:** 1.0  
**Last Updated:** 2025-10-23  
**Next Review:** End of Phase 1 (Day 28)
