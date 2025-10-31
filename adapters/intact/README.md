# IntAct Adapter

## Overview

The IntAct adapter provides access to the **IntAct Molecular Interaction Database**, a freely available resource for molecular interaction data derived from literature curation and direct user submissions.

**Database:** IntAct
**API Endpoint:** https://www.ebi.ac.uk/intact/ws/
**Type:** Public API (no authentication required)
**Rate Limit:** ~2 requests/second (be polite to EBI servers)

## Features

- **Curated Interactions**: Manually curated protein-protein interactions from published literature
- **Experimental Evidence**: Detailed detection methods and experimental conditions
- **Confidence Scores**: MI scores (molecular interaction scores) from 0-1
- **Cross-references**: Links to UniProt, PubMed, and other databases
- **Flexible Queries**: Search by UniProt ID, gene name, organism, or custom PSICQUIC queries
- **Batch Support**: Query multiple proteins simultaneously
- **Evidence Types**: Multiple detection methods and interaction types

## Installation

The adapter is included in PharmForge. No additional installation required.

```python
from adapters.intact import IntActAdapter

adapter = IntActAdapter()
```

## Basic Usage

### 1. Query Single Protein by UniProt ID

```python
import asyncio
from adapters.intact import IntActAdapter

adapter = IntActAdapter()

async def get_interactions():
    # Query TP53 (UniProt ID: P04637)
    result = await adapter.execute(
        "P04637",
        organism="human",
        min_mi_score=0.4
    )

    if result.success:
        data = result.data
        print(f"Found {data['num_interactions']} interactions")

        # Print top interactions
        for interaction in data['interactions'][:5]:
            print(f"{interaction['protein_a_name']} <-> {interaction['protein_b_name']}")
            print(f"  MI Score: {interaction['mi_score']:.3f}")
            print(f"  Method: {interaction['detection_method']}")
            print(f"  Type: {interaction['interaction_type']}")
            print()

asyncio.run(get_interactions())
```

### 2. Query Multiple Proteins

```python
async def get_multiple_interactions():
    # Query multiple proteins
    result = await adapter.execute(
        ["P04637", "Q00987", "P38398"],  # TP53, MDM2, BRCA1
        organism="human",
        min_mi_score=0.5
    )

    if result.success:
        data = result.data
        print(f"Proteins analyzed: {', '.join(data['proteins'])}")
        print(f"Total interactions: {data['num_interactions']}")
        print(f"Detection methods used: {', '.join(data['detection_methods'][:5])}")

asyncio.run(get_multiple_interactions())
```

### 3. Query by Gene Name

```python
async def query_by_gene():
    result = await adapter.execute(
        {"identifiers": ["TP53", "MDM2"], "organism": "human"},
        min_mi_score=0.6
    )

    if result.success:
        print(f"Found {result.data['num_interactions']} high-confidence interactions")

asyncio.run(query_by_gene())
```

### 4. Advanced Query with Filters

```python
async def advanced_query():
    result = await adapter.execute(
        "P04637",
        organism="human",
        min_mi_score=0.7,
        interaction_type="phosphorylation",
        max_results=100
    )

    if result.success:
        # Filter phosphorylation interactions
        for interaction in result.data['interactions']:
            if 'phosphorylation' in interaction['interaction_type'].lower():
                print(f"{interaction['protein_a_name']} phosphorylates {interaction['protein_b_name']}")
                print(f"  Evidence: {interaction['publication']}")

asyncio.run(advanced_query())
```

### 5. Custom PSICQUIC Query

```python
async def custom_query():
    # Raw PSICQUIC query syntax
    result = await adapter.execute(
        {"query": "id:P04637 AND taxid:9606 AND type:phosphorylation"},
        min_mi_score=0.5
    )

    if result.success:
        print(f"Custom query found {result.data['num_interactions']} interactions")

asyncio.run(custom_query())
```

## Input Formats

The adapter accepts multiple input formats:

### String (Single Identifier)
```python
result = await adapter.execute("P04637")  # UniProt ID
result = await adapter.execute("TP53")    # Gene name
```

### List (Multiple Identifiers)
```python
result = await adapter.execute(["P04637", "Q00987", "P38398"])
```

### Dictionary (Structured Query)
```python
result = await adapter.execute({
    "identifiers": ["TP53", "MDM2"],
    "organism": "human"
})
```

### PSICQUIC Query String
```python
result = await adapter.execute({
    "query": "id:P04637 AND taxid:9606"
})
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `organism` | str/int | None | Organism name or NCBI taxonomy ID |
| `min_mi_score` | float | 0.4 | Minimum MI score threshold (0-1) |
| `max_results` | int | 200 | Maximum number of results |
| `interaction_type` | str | None | Filter by interaction type |
| `detection_method` | str | None | Filter by detection method |
| `negative_interactions` | bool | False | Include negative results |
| `format` | str | "json" | Output format (json, tab25, tab27) |

## Output Format

```json
{
    "interactions": [
        {
            "protein_a": "P04637",
            "protein_b": "Q00987",
            "protein_a_name": "P53_HUMAN",
            "protein_b_name": "MDM2_HUMAN",
            "mi_score": 0.75,
            "detection_method": "anti tag coimmunoprecipitation",
            "interaction_type": "physical association",
            "publication": "pubmed:12345678",
            "interaction_ac": "EBI-123456",
            "source": "intact",
            "organism_a": "9606",
            "organism_b": "9606"
        }
    ],
    "num_interactions": 42,
    "num_proteins": 15,
    "proteins": ["P04637", "Q00987", ...],
    "detection_methods": ["anti tag coimmunoprecipitation", ...],
    "interaction_types": ["physical association", ...],
    "query": "id:P04637",
    "parameters": {
        "organism": "9606",
        "min_mi_score": 0.4,
        "max_results": 200
    }
}
```

## Supported Organisms

Common organisms with shortcuts:

| Organism | Name | NCBI Taxon ID |
|----------|------|---------------|
| Human | `human` or `homo_sapiens` | 9606 |
| Mouse | `mouse` or `mus_musculus` | 10090 |
| Rat | `rat` or `rattus_norvegicus` | 10116 |
| Zebrafish | `zebrafish` or `danio_rerio` | 7955 |
| Fruit Fly | `fruit_fly` or `drosophila_melanogaster` | 7227 |
| C. elegans | `c_elegans` or `caenorhabditis_elegans` | 6239 |
| Yeast | `yeast` or `saccharomyces_cerevisiae` | 4932 |
| E. coli | `e_coli` or `escherichia_coli` | 562 |

## Understanding MI Scores

**MI Score** (Molecular Interaction Score) ranges from 0 to 1:

- **0.0 - 0.4**: Low confidence (use with caution)
- **0.4 - 0.6**: Medium confidence (default threshold)
- **0.6 - 0.8**: High confidence (recommended for most analyses)
- **0.8 - 1.0**: Very high confidence (gold standard)

## Common Interaction Types

- `physical association` - General physical interaction
- `direct interaction` - Direct binding
- `colocalization` - Same cellular location
- `phosphorylation` - Post-translational modification
- `ubiquitination` - Protein degradation signal
- `methylation` - Epigenetic modification
- `enzymatic reaction` - Enzyme-substrate relationship

## Detection Methods

Common experimental methods:

- `anti tag coimmunoprecipitation` - Antibody-based pull-down
- `two hybrid` - Yeast two-hybrid screening
- `pull down` - Direct protein pull-down
- `fluorescence technology` - FRET, BRET, etc.
- `x-ray crystallography` - Structural evidence
- `mass spectrometry` - MS-based detection

## PSICQUIC Query Syntax

For advanced queries, use PSICQUIC syntax:

```python
# Query by ID
"id:P04637"

# Query by organism
"taxid:9606"

# Query by interaction type
'type:"phosphorylation"'

# Query by detection method
'detmethod:"two hybrid"'

# Combine with AND/OR
"id:P04637 AND taxid:9606 AND type:phosphorylation"
```

## Comparison with STRING-DB

| Feature | IntAct | STRING-DB |
|---------|--------|-----------|
| **Data Source** | Literature curation | Predicted + experimental |
| **Evidence Type** | Experimental only | Multi-evidence |
| **Curation** | Manual | Automated + manual |
| **Confidence** | MI scores | Combined scores |
| **Publications** | Direct links | Aggregated |
| **Best For** | Validated interactions | Network analysis |

**Recommendation**: Use IntAct for experimentally validated interactions and STRING-DB for comprehensive network analysis. Combine both for robust results.

## Integration with PharmForge

### Target Validation Pipeline

```python
async def validate_target(uniprot_id: str):
    # Get interactions from IntAct
    intact_result = await intact_adapter.execute(
        uniprot_id,
        organism="human",
        min_mi_score=0.6
    )

    # Get network from STRING-DB
    string_result = await string_adapter.execute(
        uniprot_id,
        species="human",
        required_score=600
    )

    # Compare results
    intact_partners = set(i['protein_b'] for i in intact_result.data['interactions'])
    string_partners = set(i['protein_b'] for i in string_result.data['interactions'])

    validated = intact_partners.intersection(string_partners)
    print(f"Validated interactions: {len(validated)}")

    return {
        "intact_only": intact_partners - string_partners,
        "string_only": string_partners - intact_partners,
        "validated": validated
    }
```

## Caching

Results are automatically cached using the PharmForge caching system. Cache keys include:

- Adapter name
- Query parameters
- Organism
- MI score threshold
- All filters

To bypass cache:
```python
result = await adapter(input_data, use_cache=False)
```

## Error Handling

```python
result = await adapter.execute("INVALID_ID")

if not result.success:
    print(f"Error: {result.error}")
    # Handle error gracefully
else:
    # Process results
    pass
```

## Performance Tips

1. **Use UniProt IDs** when possible (faster than gene name resolution)
2. **Set appropriate MI score threshold** (default 0.4 is balanced)
3. **Limit results** with `max_results` for faster queries
4. **Enable caching** for repeated queries (enabled by default)
5. **Batch queries** when possible (single request for multiple proteins)

## Troubleshooting

### No Results Found

- Check if protein exists in IntAct database
- Try lowering `min_mi_score` threshold
- Verify organism is correct
- Try alternative identifiers (gene name vs UniProt ID)

### Timeout Errors

- Reduce `max_results`
- Simplify query (remove filters)
- Check EBI server status: https://www.ebi.ac.uk/intact/

### Low MI Scores

- IntAct uses different scoring than STRING-DB
- Default threshold (0.4) is reasonable for most analyses
- Consider interaction type and detection method, not just score

## API Documentation

Full IntAct API documentation:
- REST API: https://www.ebi.ac.uk/intact/ws/
- PSICQUIC: https://psicquic.github.io/
- Database: https://www.ebi.ac.uk/intact/

## Citation

If you use IntAct data in publications:

> Orchard S, et al. (2014) The MIntAct project--IntAct as a common curation platform for 11 molecular interaction databases. Nucleic Acids Res. 42:D358-63.

## License

IntAct data is freely available under CC BY 4.0 license.

## Version History

- **1.0.0** (2025-10-31): Initial release
  - PSICQUIC query support
  - MI score filtering
  - Batch queries
  - Multiple organism support
  - Flexible input formats
