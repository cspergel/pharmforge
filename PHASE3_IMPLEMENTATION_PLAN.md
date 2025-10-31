# Phase 3 Implementation Plan

**Date:** October 26, 2025
**Working Directory:** `C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2`
**Duration:** Weeks 9-12 (Days 57-84)
**Status:** Starting Implementation

---

## Session Status

### ✅ Completed Today
- PDB-REDO adapter fixed (URL pattern updated)
- Adapter inventory complete (39 total adapters)
- 5 new FREE/OPEN adapters created (BioGRID, STRING-DB, GEO, pkCSM, KEGG)
- 3 adapters validated (OpenTargets, PDB-REDO, LLM Retrosynthesis)

### 🎯 Current Focus
- Phase 3: Polish & Launch
- Frontend integration per design spec
- Backend runtime fixes from phase3_weeks9-12_updated.md

---

## Phase 3 Objectives (from phase3_weeks9-12_updated.md)

By end of Week 12:
- ✅ Published validation benchmarks (DUD-E, TDC)
- ✅ Preprint submitted (ChemRxiv)
- ✅ AWS cloud infrastructure deployed
- ✅ Stripe billing integrated
- ✅ Beta signup flow live
- ✅ GitHub publicly launched (500+ stars)
- ✅ 3 academic partnerships signed
- ✅ First 10-20 paying customers

---

## Immediate Tasks (This Week)

### 1. Backend Runtime Fixes ⚡ HIGH PRIORITY

**From phase3_weeks9-12_updated.md:**

#### A. Score Normalization (backend/core/scoring_utils.py)
```python
def to01(val: float, lo: float, hi: float) -> float:
    if lo == hi: return 0.0
    v = max(min(val, hi), lo)
    return (v - lo) / (hi - lo)

def vina_affinity_to01(kcal: float) -> float:
    return to01(-kcal, 4.0, 12.0)  # maps [-12,-4] -> [1,0]

def synthesis_steps_to01(steps: int) -> float:
    steps = max(1, min(int(steps), 10))
    return to01(11 - steps, 1.0, 10.0)
```

#### B. Health Endpoint (backend/app/health.py)
```python
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

#### C. Async Background Tasks (backend/app/run_tasks.py)
```python
import asyncio
def run_async_in_background(coro_fn, *args, **kwargs):
    def _runner():
        asyncio.run(coro_fn(*args, **kwargs))
    return _runner
```

#### D. Registry Assertions
```python
required = ["vina_docking", "tdc_admet", "aizynthfinder"]
for name in required:
    assert registry.get(name) is not None, f"Adapter not registered: {name}"
```

### 2. Frontend Integration 🎨 HIGH PRIORITY

**Based on PharmForge_Frontend_Design_Spec.md:**

#### Phase 1: Infrastructure Setup
- [ ] Evaluate Streamlit vs React approach
- [ ] Set up frontend dev environment
- [ ] Configure API URL (`PHARMFORGE_API_URL`)
- [ ] Create design system (colors, typography)

#### Phase 2: Core Components
- [ ] Natural language input (ChatGPT-style)
- [ ] Command palette (⌘K functionality)
- [ ] Molecule viewer (3D visualization)
- [ ] Results dashboard
- [ ] Pipeline builder

#### Phase 3: Data Visualization
- [ ] Plotly/D3.js integration
- [ ] Interactive charts (docking scores, ADMET)
- [ ] Retrosynthesis tree viewer
- [ ] Progress indicators

### 3. Benchmark Suite 🧪 MEDIUM PRIORITY

**From phase3_weeks9-12_updated.md:**

#### A. DUD-E Benchmark
- [ ] Create `tests/benchmarks/dud_e_benchmark.py`
- [ ] Test docking enrichment
- [ ] Rank actives vs decoys
- [ ] Generate ROC curves

#### B. TDC Benchmark
- [ ] Test ADMET prediction accuracy
- [ ] Compare vs published baselines
- [ ] Generate performance metrics

---

## Week-by-Week Breakdown

### Week 9: Validation & Preprint (Days 57-63)

**Days 57-59: Benchmark Implementation**
- Implement DUD-E benchmark suite
- Implement TDC ADMET benchmark
- Generate validation data

**Days 60-61: Results Analysis**
- Analyze benchmark results
- Create figures and tables
- Document methodology

**Days 62-63: Preprint Writing**
- Write preprint draft
- Create supplementary materials
- Submit to ChemRxiv

### Week 10: Cloud Deployment (Days 64-70)

**Days 64-66: AWS Infrastructure**
- Set up Terraform configuration
- Deploy PostgreSQL RDS
- Configure S3 storage
- Set up ALB + target groups

**Days 67-68: Docker Deployment**
- Build production containers
- Configure secrets management
- Set up monitoring (CloudWatch)

**Days 69-70: Testing & Validation**
- Run smoke tests in cloud
- Load testing
- Security audit

### Week 11: Billing & Beta (Days 71-77)

**Days 71-73: Stripe Integration**
- Implement Stripe checkout
- Create subscription plans
- Add webhook handlers
- Test billing flow

**Days 74-76: Beta Signup**
- Create signup landing page
- Implement waitlist
- Set up email notifications
- Create onboarding flow

**Day 77: Launch Prep**
- Final QA testing
- Update documentation
- Prepare launch materials

### Week 12: Launch & Growth (Days 78-84)

**Days 78-79: GitHub Launch**
- Make repository public
- Submit to Show HN, r/bioinformatics
- Tweet announcement
- Email academic contacts

**Days 80-82: User Onboarding**
- Onboard first 10 beta users
- Collect feedback
- Fix critical bugs
- Create tutorials

**Days 83-84: Iteration**
- Implement user feedback
- Polish documentation
- Plan next sprint

---

## Technical Architecture Updates

### Adapter Count: 39 Total

| Category | Count | Status |
|----------|-------|--------|
| Molecular Databases | 5 | ✅ |
| Docking & Scoring | 3 | ✅ |
| Molecular Generation | 4 | ✅ |
| Retrosynthesis | 2 | ✅ |
| ADMET & Toxicity | 2 | ✅ |
| Target Prediction | 2 | ✅ |
| Protein Structure | 4 | ✅ |
| Molecular Dynamics | 1 | ✅ |
| Literature & Patents | 5 | ✅ |
| Clinical & Adverse Events | 2 | ✅ |
| Pathway & Systems Biology | 2 | ✅ |
| Gene Expression | 2 | ✅ |
| Protein Interactions | 2 | ✅ NEW |
| Target-Disease Associations | 1 | ✅ |

### Infrastructure

**Current:**
- Docker Compose dev environment
- FastAPI backend
- 39 adapters (38 production-ready)
- PostgreSQL database
- Redis cache
- Celery worker
- GPU support (RTX 5080)

**Needed:**
- Frontend (React/Streamlit)
- AWS cloud deployment
- Stripe billing
- Monitoring (CloudWatch)
- CI/CD pipeline

---

## Agent Task Assignments

### Agent 1: Backend Runtime Fixes
**Tasks:**
- Implement scoring_utils.py
- Update health endpoint
- Add async background task wrapper
- Add registry assertions
- Write tests for all fixes

### Agent 2: Frontend Infrastructure
**Tasks:**
- Evaluate Streamlit vs React
- Set up frontend dev environment
- Implement design system
- Create API client
- Build command palette (⌘K)

### Agent 3: Benchmark Suite
**Tasks:**
- Implement DUD-E benchmark
- Implement TDC benchmark
- Create benchmark runner
- Generate validation reports

### Agent 4: Phase 3 Documentation
**Tasks:**
- Update README with new adapters
- Create deployment guide
- Write API documentation
- Create user tutorials

---

## Success Metrics

### Technical Metrics
- [ ] All 39 adapters integrated into frontend
- [ ] Health endpoint returns proper status codes
- [ ] Score normalization working (0-1, higher=better)
- [ ] Benchmarks pass (DUD-E, TDC)
- [ ] Frontend loads in <2 seconds
- [ ] API response time <500ms (cached)

### Product Metrics
- [ ] 500+ GitHub stars
- [ ] 10-20 paying customers
- [ ] $2k-$5k MRR
- [ ] 3 academic partnerships
- [ ] Preprint submitted
- [ ] AWS deployment live

---

## Risk Mitigation

### Technical Risks
1. **Frontend Complexity** - Mitigate: Start with Streamlit MVP
2. **AWS Costs** - Mitigate: Use free tier, monitor costs
3. **Benchmark Accuracy** - Mitigate: Compare to published baselines

### Product Risks
1. **User Acquisition** - Mitigate: Academic partnerships first
2. **Billing Integration** - Mitigate: Use Stripe test mode extensively
3. **Scale Issues** - Mitigate: Load testing before launch

---

## Next Actions (Immediate)

1. ✅ **Create this plan** - Done
2. 🚀 **Launch 4 parallel agents:**
   - Backend runtime fixes
   - Frontend infrastructure
   - Benchmark suite
   - Documentation

3. 📊 **Track progress:**
   - Use TodoWrite tool
   - Update this plan daily
   - Review with user weekly

---

**Plan Created:** October 26, 2025
**Working Directory:** claude-code-agents-wizard-v2
**GPU:** ✅ Enabled (RTX 5080)
**Adapters:** 39 (38 production-ready)
**Ready to Execute:** YES
