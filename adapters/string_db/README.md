# STRING-DB Adapter

Adapter for the STRING Database (Search Tool for Retrieval of Interacting Genes/Proteins) - a comprehensive database of known and predicted protein-protein interactions.

## Overview

This adapter provides access to the FREE STRING-DB API to query:
- Protein-protein interaction networks
- Interaction confidence scores across multiple evidence channels
- Functional enrichment analysis
- Network statistics and PPI enrichment

## Features

- **Protein Interaction Networks**: Get networks of interacting proteins with confidence scores
- **Multiple Evidence Channels**: Scores from experimental data, databases, text mining, coexpression, etc.
- **Flexible Query Options**: Query by protein name, gene symbol, or STRING ID
- **Multi-organism Support**: Built-in support for common model organisms
- **Functional Enrichment**: GO terms, KEGG pathways, disease associations, tissue expression
- **Network Enrichment**: Statistical assessment of network connectivity (PPI enrichment p-values)
- **Interaction Partners**: Find proteins that interact with your query proteins

## API Endpoints Used

The adapter uses the following FREE STRING-DB API endpoints (JSON format):

1. **`/api/json/get_string_ids`** - Map gene/protein names to STRING identifiers
2. **`/api/json/network`** - Retrieve protein interaction networks
3. **`/api/json/interaction_partners`** - Get interaction partners for proteins
4. **`/api/json/enrichment`** - Functional enrichment analysis (GO, KEGG, etc.)
5. **`/api/json/ppi_enrichment`** - Network enrichment statistics

## Installation

No additional dependencies required beyond standard PharmForge requirements:
- `aiohttp` - for async HTTP requests
- `asyncio` - for async operations

## Usage

### Basic Usage

```python
from adapters.string_db import StringDBAdapter

# Initialize adapter
adapter = StringDBAdapter()

# Query single protein
result = await adapter.execute("TP53", species="human")

# Query multiple proteins
result = await adapter.execute(
    ["TP53", "MDM2", "BRCA1"],
    species="human",
    required_score=400  # Medium confidence (0-1000 scale)
)
```

### With Enrichment Analysis

```python
# Get network with functional enrichment
result = await adapter.execute(
    ["TP53", "BRCA1", "BRCA2", "ATM", "ATR"],
    species="human",
    required_score=400,
    include_enrichment=True,
    include_ppi_enrichment=True
)

# Access enrichment data
enrichment = result.data["enrichment"]
for category, terms in enrichment.items():
    print(f"{category}: {len(terms)} terms")
    for term in terms[:3]:
        print(f"  - {term['description']} (p={term['p_value']:.2e})")
```

### With Interaction Partners

```python
# Find proteins that interact with BRCA1
result = await adapter.execute(
    "BRCA1",
    species="human",
    required_score=600,  # High confidence
    include_partners=True,
    limit_partners=10
)

# Access partner data
partners = result.data["interaction_partners"]
for partner in partners:
    print(f"{partner['protein_b']}: {partner['combined_score']:.3f}")
```

### Dictionary Input Format

```python
# Using dict input with explicit parameters
result = await adapter.execute(
    {
        "identifiers": ["KRAS", "RAF1", "MAP2K1"],
        "species": "human"
    },
    required_score=700,
    network_type="physical"  # or "functional"
)
```

### Different Organisms

```python
# Common organism names are supported
result = await adapter.execute("Tp53", species="mouse")
result = await adapter.execute("Tp53", species="rat")
result = await adapter.execute("TP53", species="zebrafish")

# Or use NCBI taxonomy IDs directly
result = await adapter.execute("TP53", species=9606)  # Human
result = await adapter.execute("Tp53", species=10090)  # Mouse
```

## Parameters

### Required Parameters

- **`input_data`**: Protein identifier(s)
  - String: single protein name/ID
  - List: multiple protein names/IDs
  - Dict: `{"identifiers": [...], "species": "..."}`

### Optional Parameters

- **`species`** (str/int): Organism name or NCBI taxonomy ID [default: "human"/9606]
  - Supported names: "human", "mouse", "rat", "zebrafish", "yeast", "e_coli", etc.
  - Or use NCBI taxonomy ID: 9606 (human), 10090 (mouse), etc.

- **`required_score`** (int): Minimum confidence score, 0-1000 [default: 400]
  - 0-150: Low confidence
  - 150-400: Medium confidence
  - 400-700: High confidence
  - 700-1000: Very high confidence

- **`network_type`** (str): Type of interactions [default: "functional"]
  - "functional": All evidence types (recommended)
  - "physical": Only physical binding evidence

- **`add_nodes`** (int): Number of additional connected proteins to add [default: 0]

- **`include_partners`** (bool): Include interaction partners [default: False]

- **`include_enrichment`** (bool): Include functional enrichment [default: False]

- **`include_ppi_enrichment`** (bool): Include PPI enrichment p-value [default: False]

- **`limit_partners`** (int): Max partners per protein [default: 10]

- **`include_evidence_scores`** (bool): Include detailed evidence scores [default: True]

## Response Format

```python
{
    "interactions": [
        {
            "protein_a": "TP53",
            "protein_b": "MDM2",
            "combined_score": 0.998,  # 0-1 scale
            "string_id_a": "9606.ENSP00000269305",
            "string_id_b": "9606.ENSP00000258149",
            "evidence_scores": {
                "neighborhood": 0.0,
                "fusion": 0.0,
                "cooccurrence": 0.124,
                "coexpression": 0.092,
                "experimental": 0.993,
                "database": 0.900,
                "textmining": 0.993
            }
        }
    ],
    "num_interactions": 21,
    "protein_mapping": {
        "TP53": {
            "string_id": "9606.ENSP00000269305",
            "preferred_name": "TP53",
            "annotation": "Cellular tumor antigen p53..."
        }
    },
    "query_proteins": ["TP53", "MDM2"],
    "string_ids": ["9606.ENSP00000269305", "9606.ENSP00000258149"],
    "parameters": {
        "species": "9606",
        "required_score": 400,
        "network_type": "functional"
    },
    "enrichment": {  # if include_enrichment=True
        "Process": [
            {
                "term": "GO:0006974",
                "description": "Signal transduction in response to DNA damage",
                "number_of_genes": 7,
                "fdr": 1.36e-11,
                "p_value": 8.69e-16,
                "genes": ["TP53", "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2"]
            }
        ],
        "KEGG": [...],
        "DISEASES": [...],
        ...
    },
    "ppi_enrichment": {  # if include_ppi_enrichment=True
        "p_value": 2.33e-05,
        "number_of_nodes": 7,
        "number_of_edges": 21,
        "expected_edges": 7.0
    }
}
```

## Evidence Score Types

STRING combines evidence from multiple channels:

- **Neighborhood**: Gene co-occurrence across genomes
- **Fusion**: Gene fusion events
- **Cooccurrence**: Phylogenetic co-occurrence
- **Coexpression**: Co-expression patterns
- **Experimental**: Direct experimental evidence (e.g., yeast two-hybrid)
- **Database**: Curated interaction databases (e.g., KEGG, BioCyc)
- **Textmining**: Co-mentions in scientific literature

## Supported Organisms

Built-in support for common organisms (by name or taxonomy ID):

| Organism | Common Name | Taxonomy ID |
|----------|-------------|-------------|
| Homo sapiens | human | 9606 |
| Mus musculus | mouse | 10090 |
| Rattus norvegicus | rat | 10116 |
| Danio rerio | zebrafish | 7955 |
| Drosophila melanogaster | fruit_fly | 7227 |
| Caenorhabditis elegans | c_elegans | 6239 |
| Saccharomyces cerevisiae | yeast | 4932 |
| Escherichia coli | e_coli | 511145 |
| Arabidopsis thaliana | arabidopsis | 3702 |

Any NCBI taxonomy ID can be used directly.

## Rate Limiting

**IMPORTANT**: STRING-DB requires a **1-second delay between API requests**. This adapter automatically handles rate limiting with a default 1-second delay between calls.

## Example Use Cases

### 1. Cancer Pathway Analysis

```python
# Analyze p53 pathway proteins
cancer_genes = ["TP53", "MDM2", "MDM4", "ATM", "CHEK2", "CDKN1A"]
result = await adapter.execute(
    cancer_genes,
    species="human",
    required_score=500,
    include_enrichment=True
)
```

### 2. Finding Drug Targets

```python
# Find proteins interacting with a known drug target
result = await adapter.execute(
    "EGFR",
    species="human",
    required_score=700,
    include_partners=True,
    limit_partners=20
)
```

### 3. Cross-species Comparison

```python
# Compare networks across species
human_result = await adapter.execute("TP53", species="human")
mouse_result = await adapter.execute("Tp53", species="mouse")
```

### 4. Physical Binding Partners

```python
# Get only physical interactions
result = await adapter.execute(
    ["BRCA1", "BRCA2"],
    species="human",
    required_score=600,
    network_type="physical"
)
```

## Testing

Run the test suite:

```bash
# Quick test
python test_string_db.py --quick

# Run all tests
python test_string_db.py

# Run specific test
python test_string_db.py --test 1
```

Test coverage:
1. Single protein query
2. Multiple protein network
3. Interaction partners
4. Enrichment analysis
5. Mouse proteins
6. Physical network type
7. Dictionary input format
8. Comprehensive query (all features)
9. Error handling

## Limitations

1. **Rate Limiting**: 1 second between requests (enforced automatically)
2. **API Restrictions**: Free API has no authentication but fair use is expected
3. **Network Size**: Very large networks may be slow or timeout
4. **Enrichment**: Enrichment analysis requires at least 3-4 proteins for meaningful results
5. **Protein Names**: Ambiguous names may not map correctly; use UniProt IDs for precision

## References

- **STRING Database**: https://string-db.org/
- **API Documentation**: https://string-db.org/help/api/
- **Publication**: Szklarczyk et al. (2023) "The STRING database in 2023: protein-protein association networks and functional enrichment analyses for any sequenced genome of interest" Nucleic Acids Research

## Version

- **Version**: 1.0.0
- **Adapter Type**: API
- **Status**: Production Ready

## Contact

For issues or questions about this adapter, please refer to the PharmForge documentation or STRING-DB help resources.
