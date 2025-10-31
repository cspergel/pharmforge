# PharmForge - AI-Powered Drug Discovery Platform 🧬

An open-source workflow orchestrator for computational drug discovery, powered by Claude Code Agent Orchestration System v2.

## 🎯 What Is PharmForge?

**PharmForge** is a comprehensive platform that connects 39+ data sources and computational tools into complete *in silico* drug discovery pipelines. From target identification to lead optimization, PharmForge automates the entire computational workflow.

### Key Features

- **39 Integrated Adapters** - PubChem, ChEMBL, OpenTargets, UniProt, KEGG, STRING-DB, BioGRID, and 32 more
- **Natural Language Orchestration** - Describe your pipeline in plain English
- **Multi-Objective Optimization** - Pareto ranking for binding, ADMET, synthesis, and novelty
- **Reproducible Workflows** - Lockfiles with version control and DOIs
- **GPU-Accelerated** - AutoDock Vina, GNINA, DiffDock, OpenMM support
- **Production Ready** - 38/39 adapters tested and validated
- **Open Source** - MIT licensed, self-hostable

### AI Agent Architecture

PharmForge is built using the **Claude Code Agent Orchestration System v2**:

- **🧠 Claude (Orchestrator)** - Manages the 200k context window, creates todos, delegates tasks
- **✍️ Coder Subagent** - Implements features in clean context
- **👁️ Tester Subagent** - Validates implementations with Playwright
- **🆘 Stuck Subagent** - Escalates issues requiring human decision

## 🏗️ PharmForge Adapters (39 Total)

### Molecular Databases (5)
- **PubChem** - 110M+ compounds, properties, bioactivity
- **ChEMBL** - 2.3M+ compounds with bioactivity data
- **BindingDB** - 2.5M+ binding affinities
- **ZINC Fragments** - Drug-like fragment library
- **DrugCentral** - FDA-approved drugs and clinical data

### Docking & Scoring (3)
- **AutoDock Vina** - Fast molecular docking
- **GNINA** - CNN-based scoring
- **DiffDock** - ML-powered docking (GPU)

### Molecular Generation (4)
- **REINVENT** - RL-based generation
- **MolGAN** - GAN-based generation
- **De Novo** - Fragment-based design
- **RDKit Local** - Chemistry toolkit

### Retrosynthesis (2)
- **AiZynthFinder** - Retrosynthesis planning
- **LLM Retrosynthesis** - GPT-4-powered routes

### ADMET & Toxicity (2)
- **ADMET-AI** - MIT ADMET predictor
- **pkCSM** - 28 ADMET properties

### Target Prediction (2)
- **SwissTargetPrediction** - Target identification
- **TargetNet** - ML target prediction

### Protein Structure (4)
- **AlphaFold** - Structure prediction
- **RCSB PDB** - Experimental structures
- **PDB-REDO** - Re-refined structures
- **SWISS-MODEL** - Homology modeling

### Molecular Dynamics (1)
- **OpenMM** - MD simulations (GPU)

### Literature & Patents (5)
- **PubMed** - 35M+ biomedical articles
- **Europe PMC** - Life science literature
- **SureChEMBL** - 17M+ patent chemicals
- **Google Patents** - Patent search
- **Lens.org** - Patent analytics

### Clinical & Adverse Events (2)
- **ClinicalTrials.gov** - 450k+ clinical trials
- **FDA FAERS** - Adverse event reports

### Pathway & Systems Biology (2)
- **Reactome** - Biological pathways
- **KEGG** - Pathway database

### Gene Expression (2)
- **GTEx** - Tissue expression
- **GEO** - Gene expression datasets

### Protein Interactions (2)
- **BioGRID** - Protein-protein interactions
- **STRING-DB** - Interaction networks

### Target-Disease Associations (1)
- **OpenTargets** - Target validation

### Protein Information (1)
- **UniProt** - Protein sequences and functions

### Disease Information (1)
- **DisGeNET** - Gene-disease associations

## 🚀 Quick Start

### Prerequisites

- **Docker & Docker Compose** (for containers)
- **Python 3.9+** (for local development)
- **NVIDIA GPU** (optional, for GPU-accelerated adapters)
- **16GB+ RAM** recommended

### Installation

```bash
# Clone this repository
git clone https://github.com/your-org/pharmforge.git
cd pharmforge/claude-code-agents-wizard-v2

# Copy environment template
cp .env.example .env

# Edit .env with your API keys (optional - most adapters work without keys)
# Required: OPENAI_API_KEY (for LLM retrosynthesis)
# Optional: BIOGRID_ACCESS_KEY (free registration)

# Start all services
docker-compose up -d

# Check service health
curl http://localhost:8000/health

# Access the UI
open http://localhost:8501
```

### Your First Pipeline

```bash
# Example: Find drug candidates for EGFR
python -c "
from backend.core.pipeline import Pipeline

pipeline = Pipeline()
results = pipeline.execute(
    query='Find EGFR inhibitors with good ADMET',
    limit=10
)

print(f'Found {len(results)} candidates')
for r in results[:3]:
    print(f'  {r.smiles}: Score {r.score:.2f}')
"
```

### Running Tests

```bash
# Backend tests
docker-compose exec backend pytest

# Full integration test
python tests/e2e/test_full_pipeline.py
```

## 📊 Phase 3 Status (Current)

**Timeline:** Weeks 9-12 (Days 57-84)
**Focus:** Polish, Validation & Launch
**Status:** In Progress

### Completed ✅
- 39 adapters implemented (38 production-ready)
- 5 new FREE adapters added (BioGRID, STRING-DB, GEO, pkCSM, KEGG)
- All adapters validated and tested
- GPU support enabled (RTX 5080)
- Docker development environment ready
- Comprehensive documentation (8,000+ lines)

### In Progress 🔄
- Backend runtime fixes (score normalization, health endpoint)
- Frontend integration (Streamlit/React decision)
- Benchmark suite (DUD-E, TDC)
- AWS cloud deployment preparation
- Phase 3 documentation updates

### Planned 📋
- Validation benchmarks published
- Preprint submitted to ChemRxiv
- AWS infrastructure deployed
- Beta signup flow live
- GitHub public launch (target: 500+ stars)
- First 10-20 paying customers

### Metrics (Target by Day 84)
- **Adapters:** 39/39 ✅ (100% complete)
- **Production Ready:** 38/39 (97%)
- **Test Coverage:** 95%+
- **Documentation:** Comprehensive ✅
- **API Keys Required:** 2 (OpenAI, BioGRID - both free)
- **Monthly Cost:** $0-5 (OpenAI usage only)

## 📖 Documentation

### Quick Links
- **[Deployment Guide](DEPLOYMENT_GUIDE.md)** - Docker Compose and AWS setup
- **[User Guide](USER_GUIDE.md)** - Getting started and troubleshooting
- **[Adapter Inventory](FINAL_ADAPTER_INVENTORY.md)** - All 39 adapters documented
- **[Changelog](CHANGELOG.md)** - Version history and updates
- **[Phase 3 Plan](PHASE3_IMPLEMENTATION_PLAN.md)** - Current phase roadmap

### Architecture Docs
- **[Phase 1](docs/phases/phase1_weeks1-4_part1.md)** - Core infrastructure
- **[Phase 2](docs/phases/phase2_weeks5-8_updated.md)** - Pipeline completion
- **[Phase 3](docs/phases/phase3_weeks9-12_updated.md)** - Polish & launch
- **[Frontend Design](docs/phases/PharmForge_Frontend_Design_Spec.md)** - UI design system

## 📖 How to Use (Agent Workflow)

### Starting a Project

When you want to build something, just tell Claude your requirements:

```
You: "Build a todo app with React and TypeScript"
```

Claude will automatically:
1. Create a detailed todo list using TodoWrite
2. Delegate the first todo to the **coder** subagent
3. The coder implements in its own clean context window
4. Delegate verification to the **tester** subagent (Playwright screenshots)
5. If ANY problem occurs, the **stuck** subagent asks you what to do
6. Mark todo complete and move to the next one
7. Repeat until project complete

### The Workflow

```
USER: "Build X"
    ↓
CLAUDE: Creates detailed todos with TodoWrite
    ↓
CLAUDE: Invokes coder subagent for todo #1
    ↓
CODER (own context): Implements feature
    ↓
    ├─→ Problem? → Invokes STUCK → You decide → Continue
    ↓
CODER: Reports completion
    ↓
CLAUDE: Invokes tester subagent
    ↓
TESTER (own context): Playwright screenshots & verification
    ↓
    ├─→ Test fails? → Invokes STUCK → You decide → Continue
    ↓
TESTER: Reports success
    ↓
CLAUDE: Marks todo complete, moves to next
    ↓
Repeat until all todos done ✅
```

## 🛠️ How It Works

### Claude (The Orchestrator)
**Your 200k Context Window**

- Creates and maintains comprehensive todo lists
- Sees the complete project from A-Z
- Delegates individual todos to specialized subagents
- Tracks overall progress across all tasks
- Maintains project state and context

**How it works**: Claude IS the orchestrator - it uses its 200k context to manage everything

### Coder Subagent
**Fresh Context Per Task**

- Gets invoked with ONE specific todo item
- Works in its own clean context window
- Writes clean, functional code
- **Never uses fallbacks** - invokes stuck agent immediately
- Reports completion back to Claude

**When it's used**: Claude delegates each coding todo to this subagent

### Tester Subagent
**Fresh Context Per Verification**

- Gets invoked after each coder completion
- Works in its own clean context window
- Uses **Playwright MCP** to see rendered output
- Takes screenshots to verify layouts
- Tests interactions (clicks, forms, navigation)
- **Never marks failing tests as passing**
- Reports pass/fail back to Claude

**When it's used**: Claude delegates testing after every implementation

### Stuck Subagent
**Fresh Context Per Problem**

- Gets invoked when coder or tester hits a problem
- Works in its own clean context window
- **ONLY subagent** that can ask you questions
- Presents clear options for you to choose
- Blocks progress until you respond
- Returns your decision to the calling agent
- Ensures no blind fallbacks or workarounds

**When it's used**: Whenever ANY subagent encounters ANY problem

## 🚨 The "No Fallbacks" Rule

**This is the key differentiator:**

Traditional AI: Hits error → tries workaround → might fail silently
**This system**: Hits error → asks you → you decide → proceeds correctly

Every agent is **hardwired** to invoke the stuck agent rather than use fallbacks. You stay in control.

## 💡 Example Session

```
You: "Build a landing page with a contact form"

Claude creates todos:
  [ ] Set up HTML structure
  [ ] Create hero section
  [ ] Add contact form with validation
  [ ] Style with CSS
  [ ] Test form submission

Claude invokes coder(todo #1: "Set up HTML structure")

Coder (own context): Creates index.html
Coder: Reports completion to Claude

Claude invokes tester("Verify HTML structure loads")

Tester (own context): Uses Playwright to navigate
Tester: Takes screenshot
Tester: Verifies HTML structure visible
Tester: Reports success to Claude

Claude: Marks todo #1 complete ✓

Claude invokes coder(todo #2: "Create hero section")

Coder (own context): Implements hero section
Coder: ERROR - image file not found
Coder: Invokes stuck subagent

Stuck (own context): Asks YOU:
  "Hero image 'hero.jpg' not found. How to proceed?"
  Options:
  - Use placeholder image
  - Download from Unsplash
  - Skip image for now

You choose: "Download from Unsplash"

Stuck: Returns your decision to coder
Coder: Proceeds with Unsplash download
Coder: Reports completion to Claude

... and so on until all todos done
```

## 📁 Repository Structure

```
.
├── .claude/
│   ├── CLAUDE.md              # Orchestration instructions for main Claude
│   └── agents/
│       ├── coder.md          # Coder subagent definition
│       ├── tester.md         # Tester subagent definition
│       └── stuck.md          # Stuck subagent definition
├── .mcp.json                  # Playwright MCP configuration
├── .gitignore
└── README.md
```

## 🎓 Learn More

### Resources

- **[SEO Grove](https://seogrove.ai)** - AI-powered SEO automation platform
- **[ISS AI Automation School](https://www.skool.com/iss-ai-automation-school-6342/about)** - Join our community to learn AI automation
- **[Income Stream Surfers YouTube](https://www.youtube.com/incomestreamsurfers)** - Tutorials, breakdowns, and AI automation content

### Support

Have questions or want to share what you built?
- Join the [ISS AI Automation School community](https://www.skool.com/iss-ai-automation-school-6342/about)
- Subscribe to [Income Stream Surfers on YouTube](https://www.youtube.com/incomestreamsurfers)
- Check out [SEO Grove](https://seogrove.ai) for automated SEO solutions

## 🤝 Contributing

This is an open system! Feel free to:
- Add new specialized agents
- Improve existing agent prompts
- Share your agent configurations
- Submit PRs with enhancements

## 📝 How It Works Under the Hood

This system leverages Claude Code's [subagent system](https://docs.claude.com/en/docs/claude-code/sub-agents):

1. **CLAUDE.md** instructs main Claude to be the orchestrator
2. **Subagents** are defined in `.claude/agents/*.md` files
3. **Each subagent** gets its own fresh context window
4. **Main Claude** maintains the 200k context with todos and project state
5. **Playwright MCP** is configured in `.mcp.json` for visual testing

The magic happens because:
- **Claude (200k context)** = Maintains big picture, manages todos
- **Coder (fresh context)** = Implements one task at a time
- **Tester (fresh context)** = Verifies one implementation at a time
- **Stuck (fresh context)** = Handles one problem at a time with human input
- **Each subagent** has specific tools and hardwired escalation rules

## 🎯 Best Practices

1. **Trust Claude** - Let it create and manage the todo list
2. **Review screenshots** - The tester provides visual proof of every implementation
3. **Make decisions when asked** - The stuck agent needs your guidance
4. **Don't interrupt the flow** - Let subagents complete their work
5. **Check the todo list** - Always visible, tracks real progress

## 🔥 Pro Tips

- Use `/agents` command to see all available subagents
- Claude maintains the todo list in its 200k context - check anytime
- Screenshots from tester are saved and can be reviewed
- Each subagent has specific tools - check their `.md` files
- Subagents get fresh contexts - no context pollution!

## 📜 License

MIT - Use it, modify it, share it!

## 🙏 Credits

Built by [Income Stream Surfer](https://www.youtube.com/incomestreamsurfers)

Powered by Claude Code's agent system and Playwright MCP.

---

**Ready to build something amazing?** Just run `claude` in this directory and tell it what you want to create! 🚀
