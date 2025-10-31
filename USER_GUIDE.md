# PharmForge User Guide

**Version:** 0.3.0
**Last Updated:** October 26, 2025
**Target Audience:** Researchers, computational chemists, drug discovery scientists

---

## Table of Contents

1. [Getting Started](#getting-started)
2. [Running Your First Pipeline](#running-your-first-pipeline)
3. [Understanding Results](#understanding-results)
4. [Working with Adapters](#working-with-adapters)
5. [Pipeline Modes](#pipeline-modes)
6. [Advanced Features](#advanced-features)
7. [Troubleshooting](#troubleshooting)
8. [Best Practices](#best-practices)
9. [FAQ](#faq)

---

## Getting Started

### What is PharmForge?

PharmForge is a workflow orchestrator that connects **39 data sources and computational tools** into complete drug discovery pipelines. Instead of manually querying databases, running docking, and analyzing results, PharmForge automates the entire process.

### Key Concepts

- **Adapter** - A connector to a data source or tool (e.g., PubChem, AutoDock Vina)
- **Pipeline** - A sequence of adapters that work together
- **Run** - A single execution of a pipeline with specific inputs
- **Result** - The output from a pipeline (compounds, scores, rankings)
- **Lockfile** - A snapshot of a run that enables reproduction

### System Requirements

- **Browser:** Chrome, Firefox, Safari (for web UI)
- **Connection:** Internet (for API adapters)
- **Optional:** GPU for accelerated docking

### Accessing PharmForge

**Local Installation:**
```bash
# Start services
docker-compose up -d

# Open web UI
open http://localhost:8501
```

**Cloud Access:**
```
https://pharmforge.org  (coming soon)
```

---

## Running Your First Pipeline

### Example 1: Find Drug Candidates for EGFR

This example demonstrates the complete workflow from target to ranked candidates.

#### Step 1: Natural Language Query

Navigate to **New Run** → **Natural Language** mode.

```
Find EGFR inhibitors with:
- Good binding affinity (< -8 kcal/mol)
- Favorable ADMET properties
- Simple synthesis (< 5 steps)
- Limit to 20 compounds
```

#### Step 2: Review Pipeline Plan

PharmForge will generate a pipeline plan:

```yaml
Pipeline Plan:
  1. OpenTargets: Validate EGFR as target
  2. ChEMBL: Fetch known EGFR inhibitors (500 compounds)
  3. RDKit: Filter by Lipinski rules
  4. AutoDock Vina: Dock to EGFR structure (PDB: 1M17)
  5. ADMET-AI: Predict ADMET properties
  6. AiZynthFinder: Plan synthesis routes
  7. Ranking: Multi-objective Pareto optimization
```

Click **Approve & Run** to start execution.

#### Step 3: Monitor Progress

The UI will show real-time progress:

```
[=====>              ] 35% Complete
Current Step: Running ADMET predictions (175/500 compounds)
Estimated Time: 3 minutes remaining
```

#### Step 4: View Results

After completion (~6 minutes), results will display:

**Top 3 Candidates:**

| Rank | SMILES | Binding (kcal/mol) | ADMET Score | Synthesis Steps | Overall Score |
|------|--------|--------------------|-------------|-----------------|---------------|
| 1 | CC(C)Cc1ccc... | -9.2 | 0.88 | 3 | 0.92 |
| 2 | COc1ccccc1... | -8.7 | 0.91 | 4 | 0.89 |
| 3 | Cc1ccc(NC(... | -8.5 | 0.85 | 3 | 0.87 |

**Visual Results:**
- Pareto frontier plot (binding vs ADMET)
- 2D structure viewer
- ADMET property radar chart
- Synthesis route tree

#### Step 5: Export Results

```bash
# Download as CSV
curl http://localhost:8000/api/v1/runs/<run_id>/export?format=csv > results.csv

# Download lockfile (for reproduction)
curl http://localhost:8000/api/v1/runs/<run_id>/lockfile > lockfile.json

# Download as SDF (for docking)
curl http://localhost:8000/api/v1/runs/<run_id>/export?format=sdf > compounds.sdf
```

---

### Example 2: Batch Processing of Compounds

Upload a list of SMILES and analyze them.

#### Step 1: Prepare Input File

Create `compounds.csv`:

```csv
smiles,name
CC(C)Cc1ccc(cc1)C(C)C(=O)O,Ibuprofen
CC(=O)Oc1ccccc1C(=O)O,Aspirin
CN1C=NC2=C1C(=O)N(C(=O)N2C)C,Caffeine
```

#### Step 2: Upload to PharmForge

Navigate to **New Run** → **Batch Upload**.

- Upload `compounds.csv`
- Select pipeline preset: **ADMET Analysis**
- Click **Run**

#### Step 3: View Results

Results table:

| Name | MW | LogP | TPSA | BBB | hERG | AMES | Pass/Fail |
|------|-----|------|------|-----|------|------|-----------|
| Ibuprofen | 206 | 3.5 | 37 | High | Low | Negative | ✅ Pass |
| Aspirin | 180 | 1.2 | 63 | Low | Low | Negative | ✅ Pass |
| Caffeine | 194 | -0.1 | 58 | High | Low | Negative | ✅ Pass |

---

### Example 3: Target Validation Workflow

Validate a protein target before starting drug discovery.

#### Natural Language Query:

```
Validate TP53 as a cancer target:
- Check disease associations
- Find expression in cancer tissues
- Identify protein interaction partners
- Review literature evidence
```

#### PharmForge Pipeline:

1. **OpenTargets** - Disease associations (20+ cancer types)
2. **GTEx** - Expression levels in 53 tissues
3. **STRING-DB** - Protein interaction network (27 partners)
4. **PubMed** - Literature search (19,917 articles)

#### Results Summary:

```
Target Validation Score: 0.94/1.0 (Excellent)

Evidence:
✅ Strong disease associations (20 cancer types)
✅ Overexpressed in tumors (GTEx)
✅ Central hub in apoptosis pathway (STRING-DB)
✅ Extensive literature support (19,917 articles)
⚠️ Challenge: Undruggable protein (no binding pocket)

Recommendation: Focus on protein-protein interaction inhibitors
```

---

## Understanding Results

### Result Components

Every PharmForge run produces:

1. **Ranked Compounds** - Sorted by overall score
2. **Individual Scores** - Binding, ADMET, synthesis, novelty
3. **Visualizations** - Plots, charts, structures
4. **Metadata** - Runtime, versions, parameters
5. **Lockfile** - Complete reproduction info

### Score Interpretation

All scores are normalized to **0-1 scale** (higher = better):

| Score Range | Interpretation | Action |
|-------------|----------------|--------|
| 0.9 - 1.0 | Excellent | Prioritize for synthesis |
| 0.7 - 0.9 | Good | Consider for follow-up |
| 0.5 - 0.7 | Moderate | Review carefully |
| 0.0 - 0.5 | Poor | Likely not viable |

### Individual Objective Scores

#### Binding Affinity Score

- **Source:** AutoDock Vina docking
- **Range:** 0-1 (higher = better binding)
- **Conversion:** -12 kcal/mol → 1.0, -4 kcal/mol → 0.0
- **Target:** > 0.7 (corresponds to < -8 kcal/mol)

#### ADMET Score

- **Source:** ADMET-AI predictor
- **Components:** Absorption, Distribution, Metabolism, Excretion, Toxicity
- **Range:** 0-1 (higher = better drug-like properties)
- **Target:** > 0.6 for lead optimization

Key ADMET properties:
- **Absorption:** Caco-2 permeability, intestinal absorption
- **Distribution:** Blood-brain barrier, VDss
- **Metabolism:** CYP450 interactions
- **Excretion:** Renal clearance
- **Toxicity:** AMES mutagenicity, hERG inhibition

#### Synthesis Score

- **Source:** AiZynthFinder retrosynthesis
- **Range:** 0-1 (higher = simpler synthesis)
- **Conversion:** 1 step → 1.0, 10+ steps → 0.0
- **Target:** > 0.5 (< 6 steps)

#### Novelty Score

- **Source:** Tanimoto similarity to ChEMBL
- **Range:** 0-1 (higher = more novel)
- **Calculation:** 1 - max_similarity
- **Target:** > 0.3 for new chemical matter

### Overall Score Calculation

PharmForge uses **Pareto ranking** by default:

```
1. Find Pareto frontier (non-dominated solutions)
2. Assign front ranks (1 = best, 2 = second best, ...)
3. Within each front, sort by crowding distance
4. Normalize to 0-1 scale
```

**Alternative:** Weighted scoring

```python
overall_score = (
    w_binding * binding_score +
    w_admet * admet_score +
    w_synthesis * synthesis_score +
    w_novelty * novelty_score
)

# Default weights: [0.4, 0.3, 0.2, 0.1]
```

### Pareto Frontier Plot

The Pareto plot shows trade-offs between objectives:

```
     ADMET Score
          ↑
      1.0 │        ●  ← Pareto optimal
          │      ●   ●
      0.8 │    ●       ●
          │  ●     ○  ○  ← Dominated
      0.6 │○   ○  ○
          │  ○
      0.4 └──────────────────→
        -4  -6  -8  -10  -12
           Binding (kcal/mol)
```

**Interpretation:**
- **● (filled)** - Pareto optimal (no other compound is better in all objectives)
- **○ (hollow)** - Dominated (at least one other compound is strictly better)

---

## Working with Adapters

### Adapter Categories (39 Total)

PharmForge provides 39 adapters across 14 categories. See [README.md](README.md) for complete list.

### Checking Adapter Status

```bash
# List all adapters
curl http://localhost:8000/api/v1/adapters

# Check specific adapter
curl http://localhost:8000/api/v1/adapters/pubchem

# Test adapter
curl http://localhost:8000/api/v1/adapters/pubchem/test
```

Example response:

```json
{
  "name": "pubchem",
  "version": "1.0.0",
  "status": "active",
  "requires_auth": false,
  "rate_limit": "5 req/sec",
  "cache_ttl": 86400,
  "last_success": "2025-10-26T10:30:00Z"
}
```

### Using Individual Adapters

You can call adapters directly via API:

```bash
# Query PubChem for aspirin
curl -X POST http://localhost:8000/api/v1/adapters/pubchem/execute \
  -H "Content-Type: application/json" \
  -d '{"input": "CC(=O)Oc1ccccc1C(=O)O"}'

# Run docking
curl -X POST http://localhost:8000/api/v1/adapters/vina_docking/execute \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "target_pdb": "1M17",
    "box_center": [15.0, 20.0, 10.0],
    "box_size": [20.0, 20.0, 20.0]
  }'
```

### Adapter Dependencies

Some adapters require others:

| Adapter | Requires | Reason |
|---------|----------|--------|
| Vina Docking | RDKit | SMILES → 3D structure |
| REINVENT | RDKit | Molecule validation |
| AiZynthFinder | RDKit | SMILES canonicalization |
| DiffDock | AlphaFold/PDB | Protein structure |

PharmForge handles these dependencies automatically.

---

## Pipeline Modes

PharmForge supports three pipeline modes:

### 1. Natural Language Mode

**Best for:** Ad-hoc queries, exploratory research

```
User Query:
"Find kinase inhibitors with BBB penetration"

PharmForge interprets:
- Target family: kinases
- Constraint: BBB permeability > 0.5
- Output: Ranked compounds
```

**Advantages:**
- No coding required
- Flexible queries
- Automatic pipeline generation

**Limitations:**
- May misinterpret complex queries
- Limited fine-grained control

### 2. Batch Processing Mode

**Best for:** Analyzing existing compound libraries

**Input formats:**
- CSV (SMILES + metadata)
- SDF (structure files)
- SMILES list (one per line)

**Example CSV:**

```csv
smiles,name,source
CC(C)Cc1ccc(cc1)C(C)C(=O)O,Ibuprofen,ChEMBL
CC(=O)Oc1ccccc1C(=O)O,Aspirin,Literature
```

**Pipeline presets:**
- Property calculation (MW, LogP, TPSA)
- ADMET prediction
- Docking analysis
- Retrosynthesis planning
- Full analysis (all of the above)

### 3. Evolution Mode (Experimental)

**Best for:** Lead optimization, structure-activity relationship exploration

```python
# Start with seed compound
seed_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

# Define optimization objectives
objectives = {
    "binding": "maximize",
    "admet": "maximize",
    "similarity": "maintain > 0.7"  # Stay close to seed
}

# Run genetic algorithm
results = pharmforge.evolve(
    seed=seed_smiles,
    objectives=objectives,
    generations=10,
    population_size=100
)
```

**Evolution strategies:**
- Genetic algorithm (crossover + mutation)
- Reinforcement learning (REINVENT)
- Bayesian optimization

---

## Advanced Features

### Custom Pipelines

Define pipelines programmatically:

```python
from backend.core.pipeline import Pipeline, PipelineStep

pipeline = Pipeline(name="My Custom Pipeline")

# Add steps
pipeline.add_step(PipelineStep(
    adapter="chembl",
    params={"target": "EGFR", "limit": 1000}
))

pipeline.add_step(PipelineStep(
    adapter="rdkit_filter",
    params={"rules": "lipinski", "strict": True}
))

pipeline.add_step(PipelineStep(
    adapter="vina_docking",
    params={"target_pdb": "1M17", "exhaustiveness": 16}
))

# Execute
results = pipeline.execute()
```

### Lockfile Reproduction

Every run generates a lockfile that captures:
- Adapter versions
- Parameters
- Random seeds
- Timestamps
- Data source DOIs

```bash
# Export lockfile
curl http://localhost:8000/api/v1/runs/<run_id>/lockfile > lockfile.json

# Reproduce run
curl -X POST http://localhost:8000/api/v1/runs/reproduce \
  -H "Content-Type: application/json" \
  -d @lockfile.json
```

Lockfile structure:

```json
{
  "run_id": "abc123",
  "created_at": "2025-10-26T10:00:00Z",
  "pipeline": {
    "version": "0.3.0",
    "steps": [
      {
        "adapter": "chembl",
        "version": "1.2.0",
        "params": {"target": "EGFR"},
        "seed": 42
      }
    ]
  },
  "data_sources": [
    {
      "name": "ChEMBL",
      "version": "33",
      "doi": "10.6019/CHEMBL.database.33"
    }
  ]
}
```

### Caching Behavior

PharmForge uses 3-tier caching:

1. **Redis (hot)** - Fast, in-memory, 24h TTL
2. **Disk (warm)** - Persistent, 7-day TTL
3. **S3 (cold)** - Long-term, 30-day TTL (cloud only)

```bash
# Check cache statistics
curl http://localhost:8000/api/v1/stats/cache

# Clear cache (if needed)
curl -X DELETE http://localhost:8000/api/v1/cache
```

**Cache key generation:**
```python
cache_key = sha256(
    adapter_name +
    adapter_version +
    json_serialize(params) +
    input_data
)
```

### API Integration

PharmForge provides a REST API for programmatic access:

```python
import requests

# Start a new run
response = requests.post("http://localhost:8000/api/v1/runs", json={
    "query": "Find EGFR inhibitors",
    "mode": "natural_language"
})

run_id = response.json()["run_id"]

# Poll for completion
while True:
    status = requests.get(f"http://localhost:8000/api/v1/runs/{run_id}/status").json()
    if status["state"] == "completed":
        break
    time.sleep(5)

# Get results
results = requests.get(f"http://localhost:8000/api/v1/runs/{run_id}/results").json()
```

Full API documentation: http://localhost:8000/docs

---

## Troubleshooting

### Common Issues

#### Issue 1: "Adapter not available"

**Symptom:** Error message "Adapter 'xyz' not found in registry"

**Causes:**
- Adapter not installed
- Docker container not running
- Adapter initialization failed

**Solutions:**
```bash
# Check adapter list
curl http://localhost:8000/api/v1/adapters

# Restart backend
docker-compose restart backend

# Check logs
docker-compose logs backend | grep "adapter"
```

#### Issue 2: "API key required"

**Symptom:** `401 Unauthorized` or "Missing API key"

**Affected adapters:**
- LLM Retrosynthesis (OpenAI)
- BioGRID (free registration)
- Google Patents (optional)

**Solution:**
```bash
# Add keys to .env
echo "OPENAI_API_KEY=sk-your-key" >> .env
echo "BIOGRID_ACCESS_KEY=your-key" >> .env

# Restart services
docker-compose restart
```

#### Issue 3: "Pipeline timeout"

**Symptom:** Run stuck at "Running..." for > 10 minutes

**Causes:**
- Large compound set (> 1000 molecules)
- Expensive adapter (DiffDock, OpenMM)
- Network issues (API adapters)

**Solutions:**
```bash
# Check worker status
docker-compose logs worker

# Scale workers
docker-compose up -d --scale worker=4

# Reduce compound count
# Edit query to add "limit to 100 compounds"
```

#### Issue 4: "Out of memory"

**Symptom:** Docker container killed, "OOMKilled" status

**Solution:**
```bash
# Increase Docker memory limit
# Docker Desktop → Settings → Resources → Memory → 8GB

# Or edit docker-compose.yml
services:
  backend:
    deploy:
      resources:
        limits:
          memory: 4G
```

#### Issue 5: "Invalid SMILES"

**Symptom:** "Failed to parse SMILES string"

**Causes:**
- Typo in SMILES
- Non-standard representation
- Special characters

**Solution:**
```python
# Validate SMILES with RDKit
from rdkit import Chem

mol = Chem.MolFromSmiles("your_smiles_here")
if mol is None:
    print("Invalid SMILES")
else:
    canonical_smiles = Chem.MolToSmiles(mol)
    print(f"Valid SMILES: {canonical_smiles}")
```

#### Issue 6: "GPU not found"

**Symptom:** "CUDA not available" for GPU adapters

**Solution:**
```bash
# Check GPU availability
nvidia-smi

# Verify Docker GPU support
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi

# If fails, see DEPLOYMENT_GUIDE.md for GPU setup
```

---

## Best Practices

### 1. Start Small

- Begin with 10-100 compounds
- Test pipeline on small set first
- Scale up after validation

### 2. Use Filters Early

```python
# ❌ Bad: Process all 100,000 compounds
pipeline.add_step(chembl, limit=100000)
pipeline.add_step(vina_docking)  # Takes hours!

# ✅ Good: Filter first
pipeline.add_step(chembl, limit=100000)
pipeline.add_step(rdkit_filter, rules="lipinski")  # Reduce to ~10,000
pipeline.add_step(vina_docking)  # Manageable
```

### 3. Leverage Caching

```bash
# Run 1: Initial analysis (slow)
run_id_1 = pharmforge.execute(query="Find EGFR inhibitors")  # 10 minutes

# Run 2: Modify ranking only (fast)
run_id_2 = pharmforge.execute(
    query="Find EGFR inhibitors",
    ranking="weighted",  # Changed parameter
    weights=[0.5, 0.3, 0.1, 0.1]
)  # 30 seconds (cache hit on adapters)
```

### 4. Monitor Resource Usage

```bash
# Check resource consumption
docker stats

# Watch logs
docker-compose logs -f backend

# Check cache hit rate (target: > 50%)
curl http://localhost:8000/api/v1/stats/cache
```

### 5. Save Lockfiles

```bash
# After every successful run
curl http://localhost:8000/api/v1/runs/<run_id>/lockfile > \
  results/egfr_screen_$(date +%Y%m%d).json
```

### 6. Validate Critical Results

```bash
# Don't trust scores blindly
# Always review:
# - Top 10 structures (visual inspection)
# - Binding poses (if docking)
# - Literature precedent (PubMed search)
# - Experimental validation (if resources allow)
```

### 7. Understand Limitations

| Adapter | Accuracy | Limitations |
|---------|----------|-------------|
| Vina Docking | ~2 kcal/mol RMSD | No flexibility, protein rigid |
| ADMET-AI | ~80% accuracy | Training data bias |
| AiZynthFinder | Path exists ≠ feasible | No yield predictions |
| Novelty | Tanimoto distance | Bioisosteres not detected |

**Always validate in vitro/in vivo!**

---

## FAQ

### General Questions

**Q: Is PharmForge free?**

A: Yes! The core platform is MIT licensed. Most adapters (36/39) use free APIs. Optional:
- OpenAI API ($0.01-0.05 per query for LLM retrosynthesis)
- BioGRID (free registration required)

**Q: Can I run PharmForge offline?**

A: Partially. Local adapters work offline (RDKit, Vina, GNINA, OpenMM). API adapters (PubChem, ChEMBL, etc.) require internet.

**Q: How accurate are the predictions?**

A: Accuracy varies by adapter:
- Docking: ±2 kcal/mol (experimental validation recommended)
- ADMET: ~80% classification accuracy
- Retrosynthesis: Route exists ≠ practical synthesis

Always validate critical results experimentally.

**Q: Can I add my own adapters?**

A: Yes! See the [Adapter Development Guide](docs/adapters/README.md).

### Technical Questions

**Q: What's the maximum compound library size?**

A: Tested up to 10,000 compounds per run. For larger screens (100k+), use batch processing with filters.

**Q: How long does a typical run take?**

A: Depends on pipeline:
- Property calculation: 1-2 minutes (1000 compounds)
- ADMET prediction: 3-5 minutes
- Docking: 5-10 minutes (100 compounds)
- Full pipeline: 10-30 minutes

**Q: Can I run multiple pipelines in parallel?**

A: Yes! Each run is independent. Scale Celery workers for concurrency:
```bash
docker-compose up -d --scale worker=4
```

**Q: How much does AWS deployment cost?**

A: ~$94/month for basic setup (eligible for free tier first year). See [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) for details.

### Scientific Questions

**Q: Which docking program is most accurate?**

A: Depends on target:
- **Vina:** Fast, reliable for virtual screening
- **GNINA:** Better scoring with CNN, slower
- **DiffDock:** ML-based, best for complex targets (requires GPU)

**Q: How should I set compound filters?**

A: Common rules:
- **Lipinski (drug-like):** MW < 500, LogP < 5, HBD < 5, HBA < 10
- **Lead-like:** MW < 350, LogP < 3.5, HBD < 3, HBA < 6
- **Fragment-like:** MW < 250, LogP < 3, HBD < 3, HBA < 3

**Q: What's a good ADMET score threshold?**

A: Context-dependent:
- **CNS drugs:** BBB > 0.7, hERG < 0.3
- **Oral drugs:** Caco-2 > 0.6, intestinal absorption > 0.7
- **Safety:** AMES negative, hepatotoxicity < 0.3

---

## Support & Resources

### Documentation
- [Deployment Guide](DEPLOYMENT_GUIDE.md) - Installation and cloud setup
- [API Docs](http://localhost:8000/docs) - REST API reference
- [Adapter Inventory](FINAL_ADAPTER_INVENTORY.md) - All 39 adapters
- [Phase 3 Plan](PHASE3_IMPLEMENTATION_PLAN.md) - Roadmap

### Community
- **GitHub Issues:** https://github.com/your-org/pharmforge/issues
- **Email Support:** support@pharmforge.org
- **Twitter:** @PharmForge

### Academic Resources
- **Publications:** (preprint coming soon)
- **Benchmark Data:** [DUD-E results](benchmarks/dud_e/)
- **Example Workflows:** [examples/](examples/)

---

**Version:** 0.3.0
**Last Updated:** October 26, 2025
**Maintained by:** PharmForge Team

**Next Steps:**
1. Try the [Quick Start example](#running-your-first-pipeline)
2. Read the [Deployment Guide](DEPLOYMENT_GUIDE.md)
3. Explore [example workflows](examples/)
4. Join our community on GitHub!
