# BioGRID Adapter

Adapter for querying protein-protein interactions from the BioGRID (Biological General Repository for Interaction Datasets) REST API.

## Overview

BioGRID is a comprehensive database of protein, genetic, and chemical interactions. This adapter provides programmatic access to:
- Protein-protein physical interactions
- Genetic interactions
- Post-translational modifications
- Cross-referenced data from multiple organisms

**Database Size (v5.0.250):**
- 2,898,895+ interactions
- 87,320+ publications
- Multiple model organisms

## Installation Requirements

```bash
pip install aiohttp
```

## Authentication

BioGRID requires a free access key for API usage.

1. Visit: https://webservice.thebiogrid.org/
2. Fill out the registration form
3. Receive your access key immediately
4. Set the key when initializing the adapter

## Usage Examples

### Basic Usage

```python
import asyncio
from adapters.biogrid import BioGRIDAdapter

async def main():
    # Initialize with your access key
    adapter = BioGRIDAdapter(access_key="your_access_key_here")

    # Query interactions for a single gene
    result = await adapter.execute("TP53", organism="human")

    if result.success:
        data = result.data
        print(f"Found {data['num_interactions']} interactions")
        print(f"Unique genes: {data['summary']['num_unique_genes']}")

asyncio.run(main())
```

### Query Multiple Genes

```python
result = await adapter.execute(
    ["BRCA1", "BRCA2"],
    organism="human",
    max_results=100
)
```

### Filter by Evidence Type (Physical Interactions Only)

```python
result = await adapter.execute(
    "EGFR",
    organism="human",
    evidence_types=[
        "Affinity Capture-MS",
        "Affinity Capture-Western",
        "Co-crystal Structure"
    ]
)
```

### Network Expansion (Include Interactors)

```python
result = await adapter.execute(
    "TP53",
    organism="human",
    include_interactors=True,  # Get first-order interactors
    max_results=200
)
```

### Query Different Organisms

```python
# Common organisms supported
result_human = await adapter.execute("TP53", organism="human")
result_mouse = await adapter.execute("Trp53", organism="mouse")
result_yeast = await adapter.execute("CDC28", organism="yeast")

# Using taxonomy IDs
result = await adapter.execute("TP53", organism="9606")  # Human
```

### Get Available Organisms

```python
result = await adapter.execute("", query_type="organisms")
organisms = result.data["organisms"]
```

### Get Evidence Types

```python
result = await adapter.execute("", query_type="evidence")
evidence_types = result.data["evidence_types"]
```

## API Parameters

### `execute()` Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_data` | str or list | Required | Gene name(s) to query |
| `organism` | str | None | Organism name or taxonomy ID |
| `evidence_types` | list | None | Filter by experimental evidence types |
| `include_interactors` | bool | False | Include first-order interactors |
| `exclude_interspecies` | bool | True | Exclude cross-species interactions |
| `max_results` | int | 10000 | Maximum results (BioGRID limit: 10000) |
| `query_type` | str | "interactions" | Query type: "interactions", "organisms", "evidence" |

## Supported Organisms

Common organism shortcuts:
- `"human"` - Homo sapiens (9606)
- `"mouse"` - Mus musculus (10090)
- `"rat"` - Rattus norvegicus (10116)
- `"yeast"` - Saccharomyces cerevisiae (559292)
- `"fly"` - Drosophila melanogaster (7227)
- `"worm"` - Caenorhabditis elegans (6239)
- `"zebrafish"` - Danio rerio (7955)
- `"arabidopsis"` - Arabidopsis thaliana (3702)
- `"ecoli"` - Escherichia coli (83333)

You can also use any taxonomy ID directly.

## Response Format

```python
{
    "num_interactions": 150,
    "interactions": [
        {
            "interaction_id": "12345",
            "gene_a": "TP53",
            "gene_b": "MDM2",
            "entrez_gene_a": "7157",
            "entrez_gene_b": "4193",
            "organism_a": "9606",
            "organism_b": "9606",
            "experimental_system": "Affinity Capture-Western",
            "experimental_system_type": "physical",
            "pubmed_id": "12345678",
            "score": "",
            "uniprot_a": "P04637",
            "uniprot_b": "Q00987"
        },
        ...
    ],
    "summary": {
        "num_interactions": 150,
        "num_unique_genes": 85,
        "num_publications": 95,
        "experimental_systems": ["Affinity Capture-MS", "Co-IP", ...],
        "interaction_types": ["physical", "genetic"],
        "organisms": ["9606"],
        "unique_genes": ["TP53", "MDM2", "ATM", ...]
    }
}
```

## Common Evidence Types

**Physical Interactions:**
- Affinity Capture-MS
- Affinity Capture-Western
- Co-crystal Structure
- Co-fractionation
- Co-immunoprecipitation
- Far Western
- FRET
- Two-hybrid
- Protein-peptide

**Genetic Interactions:**
- Synthetic Growth Defect
- Synthetic Lethality
- Synthetic Rescue
- Dosage Growth Defect
- Dosage Rescue
- Phenotypic Enhancement
- Phenotypic Suppression

## Rate Limiting

The adapter implements rate limiting (5 requests/second) to be respectful of BioGRID resources.

## Caching

Results are automatically cached to reduce API calls:

```python
# First call hits the API
result1 = await adapter("TP53", organism="human")
print(result1.cache_hit)  # False

# Second identical call uses cache
result2 = await adapter("TP53", organism="human")
print(result2.cache_hit)  # True
```

## Error Handling

```python
result = await adapter.execute("INVALID_GENE")

if not result.success:
    print(f"Error: {result.error}")
    print(f"Metadata: {result.metadata}")
```

## Testing

Run the test suite:

```bash
# Set your access key
export BIOGRID_ACCESS_KEY="your_key_here"

# Run tests
python adapters/biogrid/test_adapter.py
```

## Integration with PharmForge

The adapter follows the PharmForge adapter protocol:

```python
from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

# Inherits from AdapterProtocol
# Returns AdapterResult objects
# Supports caching and metadata
```

## Limitations

- Maximum 10,000 results per query (BioGRID API limit)
- Requires internet connection
- Requires free access key
- Rate limited to 5 requests/second

## API Documentation

- BioGRID REST API: https://wiki.thebiogrid.org/doku.php/biogridrest
- BioGRID Homepage: https://thebiogrid.org/
- Access Key Registration: https://webservice.thebiogrid.org/

## References

- Oughtred R, et al. (2021) The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions. *Protein Science*, 30(1):187-200.
- Stark C, et al. (2011) The BioGRID Interaction Database: 2011 update. *Nucleic Acids Research*, 39:D698-D704.

## Example Workflows

### Build Protein Interaction Network

```python
async def build_network(gene, max_depth=2):
    adapter = BioGRIDAdapter(access_key=access_key)

    # Get direct interactions
    result = await adapter.execute(gene, organism="human")

    # Get interactions of interactors
    if result.success:
        partners = result.data['summary']['unique_genes'][:10]
        partner_result = await adapter.execute(
            partners,
            organism="human",
            max_results=100
        )
        return partner_result.data
```

### Find Common Interactors

```python
async def find_common_interactors(gene1, gene2):
    adapter = BioGRIDAdapter(access_key=access_key)

    result1 = await adapter.execute(gene1, organism="human")
    result2 = await adapter.execute(gene2, organism="human")

    if result1.success and result2.success:
        genes1 = set(result1.data['summary']['unique_genes'])
        genes2 = set(result2.data['summary']['unique_genes'])
        common = genes1 & genes2
        return list(common)
```

### Analyze Interaction Quality

```python
async def high_quality_interactions(gene):
    adapter = BioGRIDAdapter(access_key=access_key)

    # Query high-quality evidence types only
    result = await adapter.execute(
        gene,
        organism="human",
        evidence_types=[
            "Co-crystal Structure",
            "Affinity Capture-MS",
            "Affinity Capture-Western"
        ]
    )

    if result.success:
        return result.data['interactions']
```

## Version

1.0.0

## License

Follows PharmForge licensing. BioGRID data is freely available for research use.
