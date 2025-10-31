# Scientific Literature & Patent Database Adapters

## Overview

Four new adapters have been built for searching scientific literature and patent databases, enabling comprehensive pharmaceutical research and competitive intelligence within the PharmForge platform.

---

## Adapters Built

### 1. PubMed Adapter
**Location**: `claude-code-agents-wizard-v2/adapters/pubmed/`

**Purpose**: Search NCBI PubMed for peer-reviewed scientific literature

**Key Features**:
- ✅ E-utilities API integration
- ✅ Boolean search operators (AND, OR, NOT)
- ✅ MeSH term extraction and analysis
- ✅ Abstract retrieval
- ✅ Publication date filtering
- ✅ Citation information
- ✅ Customizable rate limiting (3 req/s, 10 with API key)
- ✅ Complex query support
- ✅ Author and journal metadata
- ✅ DOI cross-referencing

**API Documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25501/

**Example**:
```python
from adapters.pubmed import PubMedAdapter

adapter = PubMedAdapter(email="user@example.com", api_key="OPTIONAL_KEY")

result = await adapter(
    "aspirin AND cancer",
    max_results=50,
    sort_by="pub_date",
    filters={"date_from": "2020/01/01"}
)

# Access results
articles = result.data['articles']
mesh_terms = result.data['summary']['common_mesh_terms']
```

---

### 2. Google Patents Adapter
**Location**: `claude-code-agents-wizard-v2/adapters/google_patents/`

**Purpose**: Search and analyze patents for prior art and competitive intelligence

**Key Features**:
- ✅ Keyword-based patent search
- ✅ CPC (Cooperative Patent Classification) filtering
- ✅ Assignee/inventor search
- ✅ Patent metadata extraction
- ✅ Chemical entity identification from text
- ✅ Date range filtering
- ✅ Patent family tracking
- ✅ Automated chemical extraction using regex patterns
- ✅ Patent URL generation

**API Documentation**: Public Google Patents interface (for production: use BigQuery)

**Example**:
```python
from adapters.google_patents import GooglePatentsAdapter

adapter = GooglePatentsAdapter()

result = await adapter({
    "query": "kinase inhibitor",
    "cpc_codes": ["C07D", "A61P35"],  # Heterocyclics + Cancer
    "assignee": "Pfizer",
    "date_after": "20200101",
    "max_results": 30
})

# Access results
patents = result.data['patents']
chemicals = result.data['summary']['chemical_entities']
top_assignees = result.data['summary']['top_assignees']
```

**Common CPC Codes**:
- `A61K` - Preparations for medical purposes
- `A61K31` - Organic active ingredients
- `C07D` - Heterocyclic compounds
- `A61P35` - Cancer treatment

---

### 3. Lens.org Adapter
**Location**: `claude-code-agents-wizard-v2/adapters/lens/`

**Purpose**: Unified search across patents and scholarly works with citation analysis

**Key Features**:
- ✅ Combined patent + scholarly search
- ✅ Citation network analysis
- ✅ Patent family information
- ✅ Cross-referencing between patents and papers
- ✅ Advanced filtering (date, jurisdiction, open access)
- ✅ Chemical entity extraction (when available)
- ✅ Field of study classification
- ✅ DOI/PMID cross-references
- ✅ Citation counting and metrics

**API Documentation**: https://docs.api.lens.org/

**Authentication**: **Required** - Free academic API token available at https://www.lens.org/

**Example**:
```python
from adapters.lens import LensAdapter

adapter = LensAdapter(api_token="YOUR_TOKEN")

result = await adapter(
    {
        "query": "CRISPR gene editing",
        "search_type": "both",  # Search patents AND papers
        "filters": {
            "date_from": "2020-01-01",
            "jurisdiction": "US",
            "open_access": True
        }
    },
    build_citation_network=True
)

# Access results
patents = result.data['patents']
papers = result.data['scholarly']
network = result.data['citation_network']
summary = result.data['summary']
```

---

### 4. Europe PMC Adapter
**Location**: `claude-code-agents-wizard-v2/adapters/europepmc/`

**Purpose**: Full-text article search with advanced chemical entity extraction

**Key Features**:
- ✅ Full-text article search
- ✅ **Chemical entity extraction** (drugs, compounds, metabolites)
- ✅ Chemical annotations API integration
- ✅ MeSH term support
- ✅ Open access filtering
- ✅ Citation retrieval
- ✅ Advanced query syntax
- ✅ Full-text XML retrieval
- ✅ Database cross-references (PubMed, DOI, PMC)
- ✅ Chemical database IDs (ChEBI, ChEMBL, etc.)

**API Documentation**: https://europepmc.org/RestfulWebService

**Example**:
```python
from adapters.europepmc import EuropePMCAdapter

adapter = EuropePMCAdapter(email="user@example.com")

result = await adapter(
    {
        "query": "immunotherapy checkpoint inhibitor",
        "filters": {
            "has_chemicals": True,
            "has_fulltext": True,
            "open_access": True,
            "date_from": "2021-01-01"
        }
    },
    fetch_annotations=True  # Extract chemical entities
)

# Access results
articles = result.data['articles']
chemicals = result.data['summary']['common_chemicals']

# Chemical annotations per article
for article in articles:
    for chemical in article['chemicals']:
        print(f"{chemical['name']} - {chemical['type']}")
```

**Advanced Filters**:
- `HAS_CHEM:Y` - Articles with chemical annotations
- `HAS_FT:Y` - Articles with full text
- `OPEN_ACCESS:Y` - Open access only
- `SRC:PMC` - PubMed Central articles

---

## Example Scripts

### 1. PubMed Search Examples
**File**: `examples/literature_search/pubmed_search_example.py`

**Demonstrates**:
- Basic literature search for drug-disease associations
- Advanced search with date filters and boolean logic
- Drug-disease relationship analysis
- Gene-drug interaction searches
- Comparative drug study searches
- MeSH term analysis

**Run**: `python examples/literature_search/pubmed_search_example.py`

---

### 2. Patent Search Examples
**File**: `examples/literature_search/patent_search_example.py`

**Demonstrates**:
- Basic patent keyword search
- CPC classification filtering
- Assignee-specific searches (company portfolios)
- Chemistry-focused patent analysis
- Temporal trend analysis
- Chemical entity extraction from patent text

**Run**: `python examples/literature_search/patent_search_example.py`

---

### 3. Citation Network & Cross-Database Examples
**File**: `examples/literature_search/citation_network_example.py`

**Demonstrates**:
- Unified patent + paper search (Lens.org)
- Chemical entity extraction (Europe PMC)
- Full-text article search
- Cross-database drug comparison
- MeSH term analysis for topic clustering
- Citation-based influential paper discovery

**Run**: `python examples/literature_search/citation_network_example.py`

---

## Use Cases

### 1. Drug Discovery Literature Review

**Scenario**: Find all research on a specific drug target

```python
# Search PubMed for mechanism studies
pubmed = PubMedAdapter()
result = await pubmed("EGFR kinase inhibitor AND mechanism of action")

# Search patents for prior art
patents = GooglePatentsAdapter()
result = await patents({
    "query": "EGFR inhibitor",
    "cpc_codes": ["A61K31", "C07D"]
})

# Get chemical entities from full-text articles
europepmc = EuropePMCAdapter()
result = await europepmc(
    "EGFR inhibitor",
    filters={"has_chemicals": True},
    fetch_annotations=True
)
```

---

### 2. Competitive Intelligence

**Scenario**: Analyze competitor patent portfolio

```python
adapter = GooglePatentsAdapter()

competitors = ["Pfizer", "Merck", "Novartis"]
portfolio = {}

for company in competitors:
    result = await adapter({
        "query": "antibody drug conjugate",
        "assignee": company,
        "date_after": "20200101"
    })

    portfolio[company] = {
        "patent_count": len(result.data['patents']),
        "chemicals": result.data['summary']['chemical_entities'],
        "cpc_codes": result.data['summary']['top_cpc_codes']
    }
```

---

### 3. Systematic Literature Review with Chemical Extraction

**Scenario**: Extract all chemicals mentioned in cancer immunotherapy papers

```python
europepmc = EuropePMCAdapter()

result = await europepmc(
    {
        "query": "cancer immunotherapy",
        "filters": {
            "has_chemicals": True,
            "open_access": True,
            "date_from": "2020-01-01"
        }
    },
    fetch_annotations=True
)

# Extract unique chemicals
all_chemicals = {}
for article in result.data['articles']:
    for chem in article['chemicals']:
        name = chem['name']
        all_chemicals[name] = all_chemicals.get(name, 0) + 1

# Top 20 most mentioned chemicals
top_chemicals = sorted(all_chemicals.items(), key=lambda x: x[1], reverse=True)[:20]
```

---

### 4. Research Trend Analysis

**Scenario**: Track publication trends over time

```python
pubmed = PubMedAdapter()

trend_data = {}
for year in range(2018, 2025):
    result = await pubmed({
        "query": "mRNA vaccine",
        "filters": {
            "date_from": f"{year}/01/01",
            "date_to": f"{year}/12/31"
        }
    })
    trend_data[year] = len(result.data['pmids'])

# Visualize trend
import matplotlib.pyplot as plt
plt.plot(trend_data.keys(), trend_data.values())
plt.title("mRNA Vaccine Publications Over Time")
plt.savefig("trend_analysis.png")
```

---

### 5. Citation Network Analysis

**Scenario**: Build citation network for influential papers

```python
lens = LensAdapter(api_token="YOUR_TOKEN")

result = await lens(
    "CAR-T cell therapy",
    search_type="scholarly",
    build_citation_network=True
)

network = result.data['citation_network']

# Export to network analysis tools
import networkx as nx
G = nx.DiGraph()

for node in network['nodes']:
    G.add_node(node['id'], title=node['title'])

for edge in network['edges']:
    G.add_edge(edge['source'], edge['target'])

# Analyze network
centrality = nx.degree_centrality(G)
influential_papers = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:10]
```

---

## Integration with PharmForge

### Workflow Example: Literature-Guided Drug Design

```python
from adapters.pubmed import PubMedAdapter
from adapters.europepmc import EuropePMCAdapter
from adapters.rdkit_local import RDKitAdapter
from adapters.admet_ai import ADMETAdapter

# 1. Search literature for disease target
pubmed = PubMedAdapter()
lit_result = await pubmed("Alzheimer's disease AND drug target")

# 2. Extract chemical entities from full-text
europepmc = EuropePMCAdapter()
chem_result = await europepmc(
    "Alzheimer therapeutic compound",
    filters={"has_chemicals": True},
    fetch_annotations=True
)

# 3. Collect SMILES from literature (simplified)
compounds = []
for article in chem_result.data['articles']:
    for chem in article['chemicals']:
        # In practice, resolve chemical names to SMILES using PubChem
        compounds.append(chem['name'])

# 4. Analyze properties of literature compounds
rdkit = RDKitAdapter()
admet = ADMETAdapter()

for smiles in compounds:
    props = await rdkit(smiles)
    admet_pred = await admet(smiles)

    print(f"Compound: {smiles}")
    print(f"MW: {props.data['molecular_weight']}")
    print(f"LogP: {props.data['logp']}")
    print(f"Solubility: {admet_pred.data['solubility']}")
```

---

## Architecture

All adapters follow the **AdapterProtocol** pattern:

```python
class AdapterProtocol(ABC):
    async def execute(self, input_data: Any, **kwargs) -> AdapterResult
    def validate_input(self, input_data: Any) -> bool
    def generate_cache_key(self, input_data: Any, **kwargs) -> str
```

**Benefits**:
- ✅ Consistent interface across all adapters
- ✅ Automatic caching support
- ✅ Error handling with AdapterResult
- ✅ Metadata tracking
- ✅ Rate limiting built-in

---

## Rate Limits

| Adapter | Rate Limit | Authentication |
|---------|-----------|----------------|
| PubMed | 3 req/s (10 with key) | Optional (email + API key) |
| Google Patents | 1 req/s (conservative) | None |
| Lens.org | Varies by tier | **Required** (API token) |
| Europe PMC | 5 req/s | None (email recommended) |

---

## Caching

All adapters support automatic caching:

```python
# First call - fetches from API
result1 = await adapter(query)  # cache_hit = False

# Second call - returns cached result
result2 = await adapter(query)  # cache_hit = True

# Disable caching
result3 = await adapter(query, use_cache=False)
```

Cache keys are generated based on:
- Adapter name and version
- Input query/parameters
- Additional kwargs

---

## Error Handling

All adapters return `AdapterResult` with consistent error handling:

```python
result = await adapter(query)

if result.success:
    data = result.data
    print(f"Success! Found {len(data)} results")
else:
    print(f"Error: {result.error}")
    print(f"Metadata: {result.metadata}")
```

---

## Authentication Setup

### PubMed (Optional)
1. Email: Recommended for NCBI contact
2. API Key: Optional, increases rate limit
3. Get key: https://www.ncbi.nlm.nih.gov/account/

### Google Patents
- No authentication required
- For production: Use Google Patents Public Datasets on BigQuery

### Lens.org (Required)
1. Register: https://www.lens.org/
2. Free academic access available
3. Get token: https://www.lens.org/lens/user/subscriptions#scholar

### Europe PMC
- No authentication required
- Email recommended

---

## File Structure

```
claude-code-agents-wizard-v2/
├── adapters/
│   ├── pubmed/
│   │   ├── __init__.py
│   │   └── adapter.py          # PubMed adapter implementation
│   ├── google_patents/
│   │   ├── __init__.py
│   │   └── adapter.py          # Google Patents adapter
│   ├── lens/
│   │   ├── __init__.py
│   │   └── adapter.py          # Lens.org adapter
│   └── europepmc/
│       ├── __init__.py
│       └── adapter.py          # Europe PMC adapter
├── examples/
│   └── literature_search/
│       ├── README.md           # Comprehensive documentation
│       ├── pubmed_search_example.py
│       ├── patent_search_example.py
│       └── citation_network_example.py
└── LITERATURE_ADAPTERS_SUMMARY.md  # This file
```

---

## Testing

To test the adapters:

```bash
# Run PubMed examples
python examples/literature_search/pubmed_search_example.py

# Run patent search examples
python examples/literature_search/patent_search_example.py

# Run citation network examples (requires Lens.org token)
python examples/literature_search/citation_network_example.py
```

---

## Future Enhancements

### Potential Additions:
1. **Semantic Scholar Adapter** - AI-powered paper search
2. **Dimensions.ai Adapter** - Research metrics and citations
3. **USPTO PatentsView** - More robust patent API
4. **EPO OPS Adapter** - European Patent Office
5. **ChemDataExtractor Integration** - Better chemical entity extraction
6. **Network Visualization** - Interactive citation networks
7. **Auto-SMILES Resolution** - Convert chemical names to SMILES via PubChem
8. **PDF Parsing** - Extract data from full-text PDFs

---

## Resources

### API Documentation
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [Europe PMC REST API](https://europepmc.org/RestfulWebService)
- [Lens.org API](https://docs.api.lens.org/)
- [Google Patents](https://support.google.com/patents/)

### Alternative APIs
- [USPTO PatentsView](https://patentsview.org/apis/api-endpoints)
- [EPO Open Patent Services](https://www.epo.org/searching-for-patents/data/web-services/ops.html)
- [Semantic Scholar](https://www.semanticscholar.org/product/api)
- [Dimensions.ai](https://www.dimensions.ai/)

### Python Libraries
- **biopython** - NCBI integration
- **ChemDataExtractor** - Chemical NER from text
- **pypatent** - USPTO patent search
- **metapub** - PubMed wrapper

---

## License & Terms of Service

Ensure compliance with each API's terms of service:
- **PubMed**: Free for academic and commercial use
- **Google Patents**: Public data, respect robots.txt
- **Lens.org**: Free academic tier, commercial licenses available
- **Europe PMC**: Free for academic and commercial use

---

## Support

For issues or questions:
1. Check adapter-specific logs
2. Review API documentation
3. Test with simple queries first
4. Verify API tokens and authentication
5. Check rate limits and delays

---

## Summary

✅ **4 adapters built** for scientific literature and patent databases
✅ **AdapterProtocol compliance** - consistent interface
✅ **Comprehensive examples** - 3 example scripts with 15+ use cases
✅ **Chemical entity extraction** - automated from patents and papers
✅ **Citation analysis** - network building and metrics
✅ **Production-ready** - rate limiting, caching, error handling
✅ **Well-documented** - README with setup, usage, and troubleshooting

These adapters enable PharmForge users to:
- Search scientific literature for drug targets and mechanisms
- Analyze patent portfolios for competitive intelligence
- Extract chemical entities from full-text articles
- Build citation networks to identify influential research
- Track research trends over time
- Integrate literature findings with computational predictions
