# IntAct Adapter - Quick Start Guide

## 5-Minute Setup

```python
from adapters.intact import IntActAdapter
import asyncio

adapter = IntActAdapter()
```

## Common Use Cases

### 1. Find Protein Interactions (UniProt ID)

```python
async def find_interactions():
    result = await adapter.execute("P04637", organism="human")
    print(f"Found {result.data['num_interactions']} interactions")

asyncio.run(find_interactions())
```

### 2. High-Confidence Interactions Only

```python
async def high_confidence():
    result = await adapter.execute(
        "P04637",
        organism="human",
        min_mi_score=0.7  # Only high-confidence (0.7-1.0)
    )

    for interaction in result.data['interactions'][:5]:
        print(f"{interaction['protein_a_name']} <-> {interaction['protein_b_name']}")
        print(f"  Score: {interaction['mi_score']:.3f}")

asyncio.run(high_confidence())
```

### 3. Multiple Proteins (Batch Query)

```python
async def batch_query():
    proteins = ["P04637", "Q00987", "P38398"]  # TP53, MDM2, BRCA1

    result = await adapter.execute(
        proteins,
        organism="human",
        min_mi_score=0.5
    )

    print(f"Network: {result.data['num_proteins']} proteins")
    print(f"Interactions: {result.data['num_interactions']}")

asyncio.run(batch_query())
```

### 4. Search by Gene Name

```python
async def gene_search():
    result = await adapter.execute(
        {"identifiers": ["TP53", "MDM2"], "organism": "human"},
        min_mi_score=0.6
    )

    print(f"Found {result.data['num_interactions']} interactions")

asyncio.run(gene_search())
```

### 5. Filter by Interaction Type

```python
async def phosphorylation():
    result = await adapter.execute(
        "P04637",
        organism="human",
        interaction_type="phosphorylation",
        min_mi_score=0.5
    )

    for interaction in result.data['interactions']:
        print(f"{interaction['protein_a_name']} -> {interaction['protein_b_name']}")
        print(f"  Method: {interaction['detection_method']}")

asyncio.run(phosphorylation())
```

## Quick Reference

### Input Formats

| Format | Example | Use Case |
|--------|---------|----------|
| String (UniProt) | `"P04637"` | Single protein |
| List | `["P04637", "Q00987"]` | Multiple proteins |
| Dict | `{"identifiers": ["TP53"], "organism": "human"}` | Structured query |
| PSICQUIC | `{"query": "id:P04637 AND taxid:9606"}` | Advanced queries |

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `organism` | None | `"human"`, `"mouse"`, or NCBI taxon ID |
| `min_mi_score` | 0.4 | Minimum confidence (0-1) |
| `max_results` | 200 | Max interactions to return |
| `interaction_type` | None | Filter by type (e.g., "phosphorylation") |

### MI Score Guide

| Score Range | Quality | Recommendation |
|-------------|---------|----------------|
| 0.0 - 0.4 | Low | Use with caution |
| 0.4 - 0.6 | Medium | Default threshold |
| 0.6 - 0.8 | High | Recommended |
| 0.8 - 1.0 | Very High | Gold standard |

### Common Organisms

| Name | Shortcut | NCBI ID |
|------|----------|---------|
| Human | `"human"` | 9606 |
| Mouse | `"mouse"` | 10090 |
| Rat | `"rat"` | 10116 |
| Yeast | `"yeast"` | 4932 |

## Output Structure

```python
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
            "publication": "pubmed:12345678"
        }
    ],
    "num_interactions": 42,
    "num_proteins": 15,
    "proteins": ["P04637", "Q00987", ...],
    "detection_methods": ["anti tag coimmunoprecipitation", ...],
    "interaction_types": ["physical association", ...]
}
```

## Common Patterns

### Get Top Interactors

```python
result = await adapter.execute("P04637", organism="human", min_mi_score=0.6)

top_5 = result.data['interactions'][:5]
for i in top_5:
    print(f"{i['protein_b_name']}: {i['mi_score']:.3f}")
```

### Filter Post-Translational Modifications

```python
ptm_types = ["phosphorylation", "ubiquitination", "methylation"]

for ptm in ptm_types:
    result = await adapter.execute(
        "P04637",
        organism="human",
        interaction_type=ptm,
        min_mi_score=0.5
    )
    print(f"{ptm}: {result.data['num_interactions']} events")
```

### Export for Cytoscape

```python
result = await adapter.execute(["P04637", "Q00987"], organism="human")

# Create edge list
edges = []
for i in result.data['interactions']:
    edges.append({
        "source": i['protein_a'],
        "target": i['protein_b'],
        "score": i['mi_score'],
        "type": i['interaction_type']
    })

import json
with open("network.json", "w") as f:
    json.dump({"edges": edges}, f, indent=2)
```

## Troubleshooting

### No Results?

1. Lower `min_mi_score` threshold (try 0.3)
2. Check protein exists in IntAct
3. Try alternative identifier (gene name vs UniProt)
4. Verify organism is correct

### Timeout?

1. Reduce `max_results`
2. Simplify query (remove filters)
3. Check EBI server status

### Different Results Than Expected?

- IntAct uses **experimental evidence only**
- Compare with STRING-DB for predicted interactions
- Check publication dates (new data added regularly)

## Next Steps

- 📖 Full documentation: [README.md](README.md)
- 💻 More examples: [example_usage.py](example_usage.py)
- 🌐 API docs: https://www.ebi.ac.uk/intact/ws/
- 🔬 Compare with: STRING-DB adapter

## One-Liners

```python
# Quick check for interactions
result = await IntActAdapter().execute("P04637", organism="human")

# High-confidence network
result = await IntActAdapter().execute(["P04637", "Q00987"], min_mi_score=0.7)

# Export to JSON
import json; json.dump(result.data, open("intact.json", "w"), indent=2)
```

## Integration Example

```python
# Combine with STRING-DB for comprehensive analysis
from adapters.intact import IntActAdapter
from adapters.string_db import StringDBAdapter

intact = IntActAdapter()
string = StringDBAdapter()

# Get experimental evidence from IntAct
intact_result = await intact.execute("P04637", organism="human", min_mi_score=0.6)

# Get predicted network from STRING
string_result = await string.execute("P04637", species="human", required_score=600)

# Find validated interactions (both databases)
intact_partners = set(i['protein_b'] for i in intact_result.data['interactions'])
string_partners = set(i['protein_b'] for i in string_result.data['interactions'])
validated = intact_partners & string_partners

print(f"IntAct only: {len(intact_partners - string_partners)}")
print(f"STRING only: {len(string_partners - intact_partners)}")
print(f"Validated (both): {len(validated)}")
```

---

**Ready to start?** Copy any example above and run it!

For comprehensive documentation, see [README.md](README.md)
