# PharmForge Marketplace - Community Hub Design
**Date:** October 26, 2025
**Status:** CONCEPT / POST-MVP
**Inspiration:** Hugging Face, Docker Hub, VS Code Marketplace

---

## Vision

**"The Hugging Face of Drug Discovery"**

A community-driven marketplace where researchers can:
- Share trained models (docking, ADMET, retrosynthesis)
- Publish datasets (binding affinity, toxicity, synthesis routes)
- Contribute custom adapters
- Benchmark and compare tools
- Collaborate on research

---

## Marketplace Categories

### 1. **Models** 🤖
Pre-trained ML models ready to use

**Examples:**
- Docking models (trained on specific protein families)
- ADMET predictors (organ-specific toxicity)
- Retrosynthesis models (reaction templates)
- Property predictors (LogP, solubility, BBB penetration)
- Generative models (molecule design, optimization)

**Metadata:**
- Model architecture (Transformer, GNN, CNN)
- Training dataset size
- Performance metrics (AUC, MAE, R²)
- License (MIT, Apache, Commercial)
- Hardware requirements (CPU, GPU, TPU)

---

### 2. **Datasets** 📊
Curated datasets for training and benchmarking

**Examples:**
- Binding affinity datasets (DUD-E, PDBbind, BindingDB)
- Toxicity datasets (Tox21, SIDER, hERG)
- Reaction datasets (USPTO, Reaxys)
- Protein structure datasets (AlphaFold DB, PDB subsets)
- Clinical trial outcomes

**Metadata:**
- Number of compounds/datapoints
- Data format (SMILES, SDF, PDB, CSV)
- Data split (train/val/test)
- Quality metrics (duplicate rate, coverage)
- License and usage restrictions

---

### 3. **Adapters** 🔌
Custom adapters for new data sources or tools

**Examples:**
- Private database connectors
- Custom ML model wrappers
- Internal LIMS integrations
- Commercial software adapters (Schrödinger, MOE)
- Experimental method adapters (NMR, MS, X-ray)

**Metadata:**
- Adapter protocol version
- Input/output schemas
- Dependencies
- Rate limits
- API credentials required?

---

### 4. **Pipelines** 🔄
Pre-configured workflows for common tasks

**Examples:**
- "Hit Identification" (screen → dock → ADMET → rank)
- "Lead Optimization" (generate → predict → synthesize → test)
- "Target Validation" (structure → binding site → design)
- "Toxicity Profiling" (multi-endpoint prediction)
- "Synthesis Planning" (retrosynthesis → route ranking)

**Metadata:**
- Adapters used
- Configuration parameters
- Expected runtime
- Success rate on benchmarks
- Use case description

---

### 5. **Benchmarks** 📈
Standardized tests for reproducible comparisons

**Examples:**
- DUD-E docking enrichment
- TDC ADMET prediction
- USPTO retrosynthesis success rate
- Molecule generation diversity metrics
- Virtual screening recall@1000

**Metadata:**
- Test set size
- Evaluation metrics
- Baseline scores
- Leaderboard rankings
- Submission guidelines

---

## Marketplace Features

### Core Features (Phase 1)

**1. Browse & Search**
```
┌─────────────────────────────────────────┐
│ 🔍 Search: "EGFR docking model"         │
│                                         │
│ Filters:                                │
│ ☑ Models  ☐ Datasets  ☐ Adapters       │
│ Category: [Docking ▼]                   │
│ License: [Any ▼]                        │
│ Sort by: [Most Downloads ▼]            │
│                                         │
│ Results: 47 models                      │
│                                         │
│ ┌───────────────────────────────────┐  │
│ │ 🏆 EGFR-Docking-v2 by @pharmalab  │  │
│ │ AutoDock Vina model trained on    │  │
│ │ 12,000 EGFR-ligand complexes      │  │
│ │                                   │  │
│ │ ⭐ 4.8/5 (234 ratings)            │  │
│ │ ⬇️ 12.3K downloads                │  │
│ │ 📊 AUC: 0.89 on DUD-E EGFR       │  │
│ │ 🔒 MIT License                    │  │
│ │                                   │  │
│ │ [View Details] [Download] [Fork]  │  │
│ └───────────────────────────────────┘  │
│                                         │
│ ┌───────────────────────────────────┐  │
│ │ EGFR-Focused-Library by @medchem  │  │
│ │ Dataset: 45K EGFR inhibitors with │  │
│ │ IC50 values from ChEMBL + patents │  │
│ │ ...                               │  │
│ └───────────────────────────────────┘  │
└─────────────────────────────────────────┘
```

**2. Item Detail Pages**
```
┌─────────────────────────────────────────┐
│ EGFR-Docking-v2                         │
│ by @pharmalab  •  Updated 2 weeks ago   │
│                                         │
│ [Download] [Fork] [Star ⭐ 234]        │
│                                         │
│ 📋 Description                          │
│ AutoDock Vina docking model specifically│
│ trained for EGFR kinase inhibitors...   │
│                                         │
│ 📊 Performance                          │
│ • DUD-E EGFR AUC: 0.89                 │
│ • Avg. docking time: 3.2s per ligand   │
│ • Tested on 5,000 compounds            │
│                                         │
│ 💾 Files                                │
│ • egfr_vina.pdbqt (2.3 MB)             │
│ • config.json (1.2 KB)                 │
│ • README.md                             │
│                                         │
│ 🔧 Usage                                │
│ ```python                               │
│ from pharmforge import download_model   │
│ model = download_model("pharmalab/egfr")│
│ results = model.dock(smiles, protein)   │
│ ```                                     │
│                                         │
│ 📈 Benchmarks                           │
│ [View full benchmark results →]        │
│                                         │
│ 💬 Community (47 comments)              │
│ @user123: "Great model, improved our..." │
│ @researcher: "How does this compare..." │
└─────────────────────────────────────────┘
```

**3. Upload / Contribute**
```
┌─────────────────────────────────────────┐
│ Upload New Model                        │
│                                         │
│ Model Name: [_________________]         │
│ Category: [Docking ▼]                   │
│ Description:                            │
│ [________________________________]      │
│ [________________________________]      │
│                                         │
│ Upload Files:                           │
│ [📁 Drag files here or click to browse]│
│                                         │
│ Metadata:                               │
│ Training Dataset Size: [____] compounds │
│ Performance Metrics:                    │
│   AUC: [____]  MAE: [____]             │
│                                         │
│ License: [MIT ▼]                        │
│ Tags: [+] docking [+] EGFR [+] kinase  │
│                                         │
│ [Cancel] [Submit for Review]            │
└─────────────────────────────────────────┘
```

---

### Advanced Features (Phase 2)

**4. Leaderboards**
Track top-performing models per task
```
┌─────────────────────────────────────────┐
│ 🏆 Docking Leaderboard - DUD-E EGFR     │
│                                         │
│ Rank  Model              AUC    Author  │
│ ─────────────────────────────────────── │
│  1.   DiffDock-EGFR    0.91   @deepdock│
│  2.   EGFR-Docking-v2  0.89   @pharmalab│
│  3.   VinaEGFR-Fine    0.87   @unichem │
│  4.   AutoDock-Base    0.82   @opensource│
│                                         │
│ [Submit Your Model]                     │
└─────────────────────────────────────────┘
```

**5. Spaces (Interactive Demos)**
Like Hugging Face Spaces - run models in browser
```
┌─────────────────────────────────────────┐
│ 🎮 Try EGFR-Docking-v2 Live             │
│                                         │
│ Input SMILES: [CCO________________]     │
│ Protein: [EGFR (PDB: 1M17) ▼]          │
│                                         │
│ [Run Docking]                           │
│                                         │
│ Results:                                │
│ • Binding Affinity: -8.2 kcal/mol      │
│ • Pose RMSD: 1.4 Å                     │
│ • [View 3D Structure]                   │
└─────────────────────────────────────────┘
```

**6. Collections / Lists**
Curate related items
```
┌─────────────────────────────────────────┐
│ 📚 Collection: "Complete EGFR Pipeline" │
│ by @pharmalab  •  ⭐ 89 stars           │
│                                         │
│ A full workflow for EGFR inhibitor      │
│ discovery, from screening to synthesis. │
│                                         │
│ Includes:                               │
│ 🤖 EGFR-Docking-v2 (model)             │
│ 📊 EGFR-Inhibitor-Dataset (dataset)    │
│ 🔄 EGFR-Opt-Pipeline (pipeline)        │
│ 📈 DUD-E-EGFR-Benchmark (benchmark)    │
│                                         │
│ [Use This Collection]                   │
└─────────────────────────────────────────┘
```

**7. Organizations**
Team accounts for labs/companies
```
┌─────────────────────────────────────────┐
│ 🏢 pharma-lab-stanford                  │
│                                         │
│ Stanford Pharma AI Lab                  │
│ 45 public models • 12 datasets          │
│ 234 followers                           │
│                                         │
│ [Follow] [Sponsor]                      │
│                                         │
│ Recent Uploads:                         │
│ • EGFR-Docking-v2 (2 weeks ago)        │
│ • Tox21-Predictor (1 month ago)        │
│ • USPTO-Templates (2 months ago)       │
└─────────────────────────────────────────┘
```

---

## Technical Architecture

### Database Schema

```sql
-- Models
CREATE TABLE marketplace_models (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(255) UNIQUE NOT NULL,  -- e.g., "pharmalab/egfr-docking-v2"
    author_id UUID REFERENCES users(id),
    category VARCHAR(50),  -- 'docking', 'admet', 'retrosynthesis', etc.
    description TEXT,
    model_type VARCHAR(50),  -- 'pytorch', 'tensorflow', 'onnx', 'vina_config', etc.
    license VARCHAR(50),
    downloads INT DEFAULT 0,
    stars INT DEFAULT 0,
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    metadata JSONB  -- performance metrics, hardware reqs, etc.
);

-- Datasets
CREATE TABLE marketplace_datasets (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(255) UNIQUE NOT NULL,
    author_id UUID REFERENCES users(id),
    category VARCHAR(50),
    description TEXT,
    num_datapoints INT,
    data_format VARCHAR(50),  -- 'smiles_csv', 'sdf', 'parquet', etc.
    license VARCHAR(50),
    downloads INT DEFAULT 0,
    stars INT DEFAULT 0,
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    metadata JSONB
);

-- Adapters
CREATE TABLE marketplace_adapters (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(255) UNIQUE NOT NULL,
    author_id UUID REFERENCES users(id),
    adapter_type VARCHAR(50),  -- 'api', 'local', 'ml', etc.
    description TEXT,
    protocol_version VARCHAR(20),
    input_schema JSONB,
    output_schema JSONB,
    license VARCHAR(50),
    downloads INT DEFAULT 0,
    stars INT DEFAULT 0,
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    metadata JSONB
);

-- Pipelines
CREATE TABLE marketplace_pipelines (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    slug VARCHAR(255) UNIQUE NOT NULL,
    author_id UUID REFERENCES users(id),
    description TEXT,
    adapter_ids UUID[],  -- Array of adapter IDs
    config JSONB,  -- Pipeline configuration
    license VARCHAR(50),
    downloads INT DEFAULT 0,
    stars INT DEFAULT 0,
    created_at TIMESTAMP,
    updated_at TIMESTAMP,
    metadata JSONB
);

-- Files (for models/datasets)
CREATE TABLE marketplace_files (
    id UUID PRIMARY KEY,
    item_id UUID,  -- References model, dataset, or adapter
    item_type VARCHAR(20),  -- 'model', 'dataset', 'adapter', 'pipeline'
    filename VARCHAR(255),
    size_bytes BIGINT,
    checksum VARCHAR(64),
    storage_url TEXT,  -- S3 URL, etc.
    created_at TIMESTAMP
);

-- Reviews/Ratings
CREATE TABLE marketplace_reviews (
    id UUID PRIMARY KEY,
    item_id UUID,
    item_type VARCHAR(20),
    user_id UUID REFERENCES users(id),
    rating INT CHECK (rating >= 1 AND rating <= 5),
    comment TEXT,
    created_at TIMESTAMP
);

-- Collections
CREATE TABLE marketplace_collections (
    id UUID PRIMARY KEY,
    name VARCHAR(255),
    description TEXT,
    author_id UUID REFERENCES users(id),
    items JSONB,  -- Array of {type, id} objects
    stars INT DEFAULT 0,
    created_at TIMESTAMP
);

-- Benchmarks
CREATE TABLE marketplace_benchmarks (
    id UUID PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    category VARCHAR(50),
    description TEXT,
    test_set_url TEXT,
    metrics JSONB,  -- Metric definitions
    created_at TIMESTAMP
);

-- Benchmark Results (leaderboard entries)
CREATE TABLE marketplace_benchmark_results (
    id UUID PRIMARY KEY,
    benchmark_id UUID REFERENCES marketplace_benchmarks(id),
    model_id UUID REFERENCES marketplace_models(id),
    scores JSONB,  -- {'auc': 0.89, 'mae': 0.12, ...}
    submitted_by UUID REFERENCES users(id),
    submitted_at TIMESTAMP,
    verified BOOLEAN DEFAULT FALSE
);
```

---

### API Endpoints

```python
# Browse & Search
GET  /api/marketplace/models
GET  /api/marketplace/datasets
GET  /api/marketplace/adapters
GET  /api/marketplace/pipelines
GET  /api/marketplace/search?q=egfr&type=model&category=docking

# Item Details
GET  /api/marketplace/models/{slug}
GET  /api/marketplace/datasets/{slug}
GET  /api/marketplace/adapters/{slug}

# Download
GET  /api/marketplace/models/{slug}/download
GET  /api/marketplace/datasets/{slug}/download

# Upload
POST /api/marketplace/models
POST /api/marketplace/datasets
POST /api/marketplace/adapters

# Interactions
POST /api/marketplace/models/{slug}/star
POST /api/marketplace/models/{slug}/review
POST /api/marketplace/models/{slug}/fork

# Benchmarks
GET  /api/marketplace/benchmarks
GET  /api/marketplace/benchmarks/{id}/leaderboard
POST /api/marketplace/benchmarks/{id}/submit

# Collections
GET  /api/marketplace/collections
POST /api/marketplace/collections
GET  /api/marketplace/collections/{id}

# Organizations
GET  /api/marketplace/orgs/{slug}
GET  /api/marketplace/orgs/{slug}/models
```

---

## Frontend Pages

### 1. Marketplace Home
```
frontend/pages/marketplace_home.py
```

**Layout:**
- Hero section: "Discover, Share, Collaborate"
- Trending models/datasets this week
- Featured collections
- Top contributors
- Recent uploads
- Search bar prominent

---

### 2. Browse Page
```
frontend/pages/marketplace_browse.py
```

**Features:**
- Filter sidebar (category, license, rating, downloads)
- Grid/list view toggle
- Sort options (trending, newest, most downloads, highest rated)
- Pagination
- Quick preview on hover

---

### 3. Item Detail Page
```
frontend/pages/marketplace_item.py
```

**Tabs:**
- Overview (description, metrics, usage)
- Files (download, preview)
- Benchmarks (performance comparison)
- Community (reviews, comments, forks)
- Version history

---

### 4. Upload Page
```
frontend/pages/marketplace_upload.py
```

**Form fields:**
- Name, description, category
- File upload (with progress bar)
- Metadata (metrics, dependencies, hardware)
- License selection
- Tags
- README editor

---

### 5. Leaderboard Page
```
frontend/pages/marketplace_leaderboards.py
```

**Features:**
- Select benchmark (DUD-E, TDC, USPTO, etc.)
- Select metric (AUC, MAE, success rate)
- Sortable table
- Filter by model type, date range
- Submit your model button

---

### 6. Spaces (Interactive Demos)
```
frontend/pages/marketplace_spaces.py
```

**Features:**
- Embedded model inference
- Input form (SMILES, protein, parameters)
- Run button
- Results visualization
- Share results link

---

## Integration with PharmForge

### 1. Adapter Registry Integration

**Auto-install marketplace adapters:**
```python
# User clicks "Install" on marketplace adapter
from backend.core.adapter_registry import registry

adapter_slug = "pharmalab/egfr-docking-v2"
marketplace_api.download_adapter(adapter_slug)
registry.register_from_marketplace(adapter_slug)

# Now available in PharmForge UI
```

**Adapter manifest:**
```json
{
  "name": "egfr_docking_v2",
  "version": "2.0.1",
  "author": "pharmalab",
  "marketplace_url": "https://pharmforge.ai/marketplace/pharmalab/egfr-docking-v2",
  "protocol_version": "1.0",
  "dependencies": ["numpy", "rdkit", "openbabel"],
  "input_schema": { "smiles": "string", "protein": "string" },
  "output_schema": { "affinity": "float", "pose": "string" }
}
```

---

### 2. Pipeline Builder Integration

**Import marketplace pipelines:**
```
┌─────────────────────────────────────────┐
│ New Run                                 │
│                                         │
│ Start from:                             │
│ ○ Blank pipeline                        │
│ ● Import from Marketplace               │
│                                         │
│ Search marketplace:                     │
│ [EGFR hit identification____________]   │
│                                         │
│ Results:                                │
│ 📦 EGFR-Complete-Pipeline by @pharmalab │
│    ⭐ 4.8 • ⬇️ 2.3K                     │
│    [Use This Pipeline]                  │
│                                         │
│ 📦 Kinase-Inhibitor-Screen by @medchem │
│    ⭐ 4.5 • ⬇️ 1.8K                     │
│    [Use This Pipeline]                  │
└─────────────────────────────────────────┘
```

---

### 3. Dataset Integration

**Use marketplace datasets in runs:**
```python
# Frontend: User selects dataset from marketplace
dataset = marketplace_api.get_dataset("chembl/egfr-inhibitors")

# Backend: Stream dataset into pipeline
for compound in dataset.iter_smiles():
    pipeline.run(compound)
```

---

### 4. Benchmark Integration

**Auto-submit results to leaderboards:**
```
┌─────────────────────────────────────────┐
│ Run Complete! ✅                        │
│                                         │
│ Your pipeline scored:                   │
│ • DUD-E EGFR AUC: 0.88                 │
│ • Runtime: 45 mins                      │
│                                         │
│ This qualifies for the leaderboard!     │
│                                         │
│ [Submit to DUD-E Leaderboard]           │
└─────────────────────────────────────────┘
```

---

## Monetization (Post-MVP)

### Free Tier
- Upload unlimited public items
- Download unlimited public items
- Star, review, fork items
- Basic analytics

### Pro Tier ($19/month)
- Private models/datasets
- Team collaboration (5 members)
- Advanced analytics (download stats, user demographics)
- Priority support

### Enterprise Tier ($199/month)
- Organization accounts
- SSO integration
- On-premise deployment
- Custom SLAs
- Unlimited team members

### Marketplace Revenue Share
- PharmForge takes 15% of paid model/dataset sales
- Authors keep 85%
- Free items remain free

---

## Community Governance

### Review Process
1. User uploads item
2. Automated checks (file size, malware scan, license validation)
3. Community review period (7 days)
4. Moderator approval
5. Published to marketplace

### Moderation
- Report abusive content
- Moderator team reviews reports
- Takedown policy for violations
- Appeal process

### Quality Standards
- README required
- License required
- Benchmarks encouraged
- Code of conduct

---

## Launch Strategy

### Phase 1: MVP (Week 10-11)
- Basic browse/search
- Upload models/datasets
- Download functionality
- Simple detail pages
- No authentication (view-only)

### Phase 2: Community (Week 12-13)
- User accounts
- Star/review system
- Comments
- Fork functionality
- Upload authentication

### Phase 3: Advanced (Week 14-16)
- Leaderboards
- Spaces (interactive demos)
- Collections
- Organizations
- Benchmark submission

### Phase 4: Monetization (Week 17+)
- Private items
- Paid models/datasets
- Pro/Enterprise tiers
- Analytics dashboard

---

## Success Metrics

### Engagement
- Monthly active contributors
- Upload rate (items/week)
- Download rate (downloads/week)
- Search queries

### Quality
- Average rating per category
- Review participation rate
- Benchmark submission rate
- Fork rate (indicates usefulness)

### Growth
- New users/week
- Items published/week
- Community size
- Social media mentions

---

## Competitive Analysis

### Hugging Face
**Strengths:** Massive model hub, Spaces, community
**Weaknesses:** Not drug discovery focused

### TorchDrug / TDC
**Strengths:** Drug-specific benchmarks
**Weaknesses:** No marketplace, limited models

### Docker Hub
**Strengths:** Simple UX, versioning
**Weaknesses:** Not ML-focused

### GitHub
**Strengths:** Version control, collaboration
**Weaknesses:** Not model-focused

**PharmForge Marketplace Advantage:**
- Drug discovery-specific
- Integrated with full platform (not just model hosting)
- Standardized adapter protocol
- Built-in benchmarking
- One-click deployment

---

## Technical Challenges

### 1. File Storage
**Challenge:** Large model files (GBs)
**Solution:** S3 + CloudFront CDN, chunked downloads

### 2. Version Control
**Challenge:** Tracking model versions
**Solution:** Git LFS-style versioning, semantic versioning

### 3. Reproducibility
**Challenge:** Ensuring models run correctly
**Solution:** Containerization (Docker), dependency pinning

### 4. Security
**Challenge:** Malicious models (data exfiltration)
**Solution:** Sandboxed execution, code review, malware scanning

### 5. Attribution
**Challenge:** Proper credit for community contributions
**Solution:** Citation.cff files, DOIs for datasets/models

---

## Next Steps (After MVP Launch)

1. **Build marketplace database schema**
2. **Create API endpoints**
3. **Design frontend pages**
4. **Implement file upload/download**
5. **Add authentication**
6. **Create initial seed content (official models/datasets)**
7. **Launch private beta with select users**
8. **Iterate based on feedback**
9. **Public launch**

---

**Timeline:** 4-6 weeks after MVP launch
**Team Size:** 2 backend, 1 frontend, 1 designer
**Infrastructure:** S3, CloudFront, PostgreSQL, Redis

---

**Key Insight:**
The marketplace becomes the moat. Once researchers publish models/datasets here, they'll keep coming back. Network effects = competitive advantage.

This is how PharmForge becomes the **default platform** for computational drug discovery.

---

**Created:** October 26, 2025
**Status:** Design phase, ready to implement post-MVP
**Owner:** PharmForge Core Team
