# PharmForge - AI Orchestrated Build

This is the **PharmForge** project integrated with the **Claude Code Agent Orchestration System v2**.

## 🎯 What This Is

A 12-week (90-day) build plan for **PharmForge**, an open-source drug discovery workflow orchestrator, managed by specialized AI agents:

- **Claude (Orchestrator)** - You, managing the big picture with 200k context and TodoWrite
- **Coder Subagent** - Implements one todo at a time in clean context
- **Tester Subagent** - Verifies implementations using Playwright screenshots
- **Stuck Subagent** - Escalates to you when ANY problem occurs (no fallbacks!)

## 📁 Project Structure

```
claude-code-agents-wizard-v2/
├── .claude/
│   ├── CLAUDE.md              # PharmForge orchestrator instructions
│   └── agents/
│       ├── coder.md          # Implementation specialist
│       ├── tester.md         # Visual testing with Playwright
│       └── stuck.md          # Human escalation point
├── docs/phases/
│   ├── BUILD_GUIDE_INDEX.md         # Navigation guide
│   ├── phase1_weeks1-4_part1.md     # Week 1: Core infrastructure
│   ├── phase1_weeks1-4_part2.md     # Weeks 2-4: Adapters & queue
│   ├── phase2_weeks5-8_updated.md   # Weeks 5-8: Pipelines & UI
│   ├── phase3_weeks9-12_updated.md  # Weeks 9-12: Launch & validation
│   ├── PharmForge_Frontend_Design_Spec.md
│   ├── CommonMistakes.md
│   └── QuickStart.md
├── backend/                   # FastAPI application
├── frontend/                  # React/Streamlit UI
├── config/                    # Configuration files
├── scripts/                   # Utility scripts
├── tests/                     # Integration tests
├── .env.example              # Environment variables template
└── docker-compose.yml        # Development environment

```

## 🚀 Quick Start

### 1. Review the Build Plan
```bash
# Start here - understand the project
cat docs/phases/claude.md

# Navigation guide
cat docs/phases/BUILD_GUIDE_INDEX.md
```

### 2. Understanding the Agent System

**Orchestrator (Claude/You):**
- Reads phase documents
- Creates detailed todo lists with TodoWrite
- Delegates tasks to subagents using Task tool
- Maintains overall project context (200k window)
- Tracks progress across all phases

**Coder Subagent:**
- Gets ONE specific todo
- Implements it completely
- NEVER uses fallbacks
- Invokes stuck agent on ANY error
- Reports completion

**Tester Subagent:**
- Gets completed implementation
- Uses Playwright MCP for visual verification
- Takes screenshots as proof
- NEVER marks failing tests as passing
- Invokes stuck agent on ANY failure

**Stuck Subagent:**
- ONLY agent that can ask you questions
- Presents clear options
- Blocks progress until you respond
- Returns your decision to calling agent
- Ensures no blind fallbacks

## 🔄 The Workflow

```
YOU: "Start Phase 1, Week 1"
    ↓
CLAUDE (orchestrator): Creates detailed todos from phase1_weeks1-4_part1.md
    ↓
CLAUDE: Invokes coder subagent for todo #1
    ↓
CODER: Implements in clean context
    ↓
    ├─→ Problem? → Invokes STUCK → You decide → Continue
    ↓
CODER: Reports completion
    ↓
CLAUDE: Invokes tester subagent
    ↓
TESTER: Playwright screenshots & verification
    ↓
    ├─→ Test fails? → Invokes STUCK → You decide → Continue
    ↓
TESTER: Reports success
    ↓
CLAUDE: Marks todo complete, moves to next
    ↓
Repeat until phase complete ✅
```

## 📖 Phase Overview

### Phase 1: Core Infrastructure (Weeks 1-4)
**Location:** `docs/phases/phase1_weeks1-4_part1.md` & `phase1_weeks1-4_part2.md`

**Deliverables:**
- Docker dev environment
- FastAPI with PostgreSQL
- Adapter protocol
- Working adapters: PubChem, ChEMBL, RDKit, TDC ADMET, Vina
- Celery job queue
- Caching layer

**Success Criteria:** Process 100 compounds through full pipeline

### Phase 2: Pipeline Completion (Weeks 5-8)
**Location:** `docs/phases/phase2_weeks5-8_updated.md`

**Deliverables:**
- AiZynthFinder retrosynthesis
- Multi-objective ranking (Pareto + weighted)
- Arcana NL orchestrator (GPT-4)
- Streamlit UI
- Lockfile generation
- Documentation

**Success Criteria:** NL query → pipeline → ranked results

### Phase 3: Launch & Validation (Weeks 9-12)
**Location:** `docs/phases/phase3_weeks9-12_updated.md`

**Deliverables:**
- DUD-E & TDC validation benchmarks
- Preprint submission
- AWS cloud deployment
- Stripe billing
- GitHub public launch
- First 10 paying customers

**Success Criteria:** 500+ stars, $2-5k MRR, 3 academic partnerships

## 🎯 How to Begin

### Step 1: Tell Claude to Start
```
You: "Let's start building PharmForge. Begin with Phase 1, Week 1."
```

### Step 2: Claude Will Automatically:
1. Read `docs/phases/phase1_weeks1-4_part1.md`
2. Create detailed todos using TodoWrite
3. Delegate first todo to coder subagent
4. Verify with tester subagent
5. Ask you via stuck agent if problems occur
6. Move through todos systematically

### Step 3: You Monitor Progress
- Check todo list status anytime
- Respond to stuck agent questions
- Review screenshots from tester
- Validate at checkpoint milestones

## 🚨 The "No Fallbacks" Rule

**Key Differentiator:**
- Traditional AI: Error → tries workaround → might fail silently
- **This system**: Error → asks YOU → you decide → proceeds correctly

Every agent is hardwired to invoke stuck agent rather than using fallbacks. You stay in control.

## 💡 Pro Tips

1. **Trust the system** - Let Claude manage the todo list
2. **Review screenshots** - Tester provides visual proof
3. **Make decisions quickly** - Stuck agent presents clear options
4. **Don't interrupt flow** - Let subagents complete their work
5. **Check todos often** - Always know where you are

## 🔧 Environment Setup

```bash
# 1. Copy environment template
cp .env.example .env

# 2. Fill in your API keys
# Edit .env with your values

# 3. Start development environment
docker-compose up -d

# 4. Verify services are running
docker-compose ps
```

## 📚 Key Resources

- **Master Documentation**: `docs/phases/claude.md`
- **Quick Start**: `docs/phases/QuickStart.md`
- **Common Mistakes**: `docs/phases/CommonMistakes.md`
- **Frontend Spec**: `docs/phases/PharmForge_Frontend_Design_Spec.md`

## 🎓 Success Metrics (Day 90)

| Metric | Target |
|--------|--------|
| Product | MVP deployed (local + cloud) |
| GitHub Stars | 500+ |
| Active Users | 50 (10 paying) |
| MRR | $2k-$5k |
| Academic Partners | 3 commitments |
| Community Adapters | 10 total |
| Documentation | Complete |
| Validation | Preprint submitted |

## 🤝 Agent Communication

**Use these markers when needed:**
```python
# TODO: [BLOCKED] Need API key
# DECISION NEEDED: Should we use X or Y?
# PERFORMANCE: This is slow, optimize?
```

## ⚡ Ready to Build?

Simply say:
```
"Let's start building PharmForge. Begin with Phase 1, Week 1."
```

Claude will take it from there, creating todos and delegating to specialized agents!

---

**Built with:** Claude Code Agent Orchestration System v2
**Project:** PharmForge - Open-source drug discovery workflow orchestrator
**Timeline:** 12 weeks (90 days)
**Let's build the future of drug discovery! 🧬🚀**
