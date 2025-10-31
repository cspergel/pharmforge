# BioGRID Adapter - Quick Reference Card

## Setup (One Time)

1. Get free access key: https://webservice.thebiogrid.org/
2. Set environment variable:
   ```bash
   export BIOGRID_ACCESS_KEY="your_key_here"
   ```

## Import

```python
from adapters.biogrid import BioGRIDAdapter
import asyncio

adapter = BioGRIDAdapter(access_key="your_key")
```

## Basic Queries

### Single Gene
```python
result = await adapter.execute("TP53", organism="human")
```

### Multiple Genes
```python
result = await adapter.execute(["BRCA1", "BRCA2"], organism="human")
```

### With Evidence Filter
```python
result = await adapter.execute(
    "EGFR",
    organism="human",
    evidence_types=["Affinity Capture-MS", "Co-crystal Structure"]
)
```

### Network Expansion
```python
result = await adapter.execute(
    "TP53",
    organism="human",
    include_interactors=True
)
```

## Common Organisms

| Shortcut | Organism | Tax ID |
|----------|----------|--------|
| "human" | Homo sapiens | 9606 |
| "mouse" | Mus musculus | 10090 |
| "rat" | Rattus norvegicus | 10116 |
| "yeast" | S. cerevisiae | 559292 |
| "fly" | D. melanogaster | 7227 |

## Response Structure

```python
{
    "num_interactions": 150,
    "interactions": [...],  # List of interaction objects
    "summary": {
        "num_unique_genes": 85,
        "num_publications": 95,
        "experimental_systems": [...],
        "unique_genes": [...]
    }
}
```

## Common Evidence Types

**Physical:**
- Affinity Capture-MS
- Affinity Capture-Western
- Co-crystal Structure
- Two-hybrid
- Co-immunoprecipitation

**Genetic:**
- Synthetic Lethality
- Synthetic Growth Defect
- Phenotypic Enhancement

## Utility Queries

### List Organisms
```python
result = await adapter.execute("", query_type="organisms")
```

### List Evidence Types
```python
result = await adapter.execute("", query_type="evidence")
```

## Error Handling

```python
result = await adapter.execute("TP53", organism="human")

if result.success:
    data = result.data
    print(f"Found {data['num_interactions']} interactions")
else:
    print(f"Error: {result.error}")
```

## Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| organism | str | None | "human", "mouse", or tax ID |
| evidence_types | list | None | Filter by evidence |
| include_interactors | bool | False | Get network |
| max_results | int | 10000 | Result limit |
| exclude_interspecies | bool | True | Filter cross-species |

## Running Tests

```bash
export BIOGRID_ACCESS_KEY="your_key"
python adapters/biogrid/test_adapter.py
```

## Quick Examples

### Example 1: Cancer Gene Network
```python
async def main():
    adapter = BioGRIDAdapter(access_key="key")
    result = await adapter.execute(
        ["TP53", "BRCA1", "EGFR"],
        organism="human",
        max_results=200
    )
    print(f"Network size: {result.data['summary']['num_unique_genes']} genes")

asyncio.run(main())
```

### Example 2: Find Common Partners
```python
async def find_common(gene1, gene2):
    adapter = BioGRIDAdapter(access_key="key")

    r1 = await adapter.execute(gene1, organism="human")
    r2 = await adapter.execute(gene2, organism="human")

    genes1 = set(r1.data['summary']['unique_genes'])
    genes2 = set(r2.data['summary']['unique_genes'])

    return genes1 & genes2
```

### Example 3: High-Quality Only
```python
result = await adapter.execute(
    "EGFR",
    organism="human",
    evidence_types=[
        "Co-crystal Structure",
        "Affinity Capture-MS"
    ]
)
```

## Need Help?

- Full docs: `adapters/biogrid/README.md`
- Tests: `adapters/biogrid/test_adapter.py`
- Examples: `adapters/biogrid/example_usage.py`
- API docs: https://wiki.thebiogrid.org/doku.php/biogridrest
