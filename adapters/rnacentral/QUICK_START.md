# RNAcentral Adapter - Quick Start Guide

Get started with the RNAcentral adapter in 5 minutes.

## What is RNAcentral?

RNAcentral is the world's largest RNA sequence database, integrating data from 43 expert databases. It's essential for:

- **RNA therapeutics** (antisense, siRNA, miRNA drugs)
- **Target identification** (disease-associated RNAs)
- **Biomarker discovery** (diagnostic RNAs)
- **RNA-protein interactions** (druggable targets)

## Installation

No special installation needed! The adapter uses standard libraries:

```bash
pip install aiohttp
```

## Quick Examples

### 1. Look Up an RNA by ID (10 seconds)

```python
from adapters.rnacentral.adapter import RNAcentralAdapter
import asyncio

async def quick_lookup():
    adapter = RNAcentralAdapter()
    result = await adapter.execute("URS0000527F89_9606")  # let-7a microRNA

    if result.success:
        print(f"Type: {result.data['rna_type']}")
        print(f"Length: {result.data['length']} nt")
        print(f"Sequence: {result.data['sequence']}")

asyncio.run(quick_lookup())
```

### 2. Search for Disease-Associated RNAs (30 seconds)

```python
async def find_targets():
    adapter = RNAcentralAdapter()

    # Find microRNAs linked to cancer
    result = await adapter.execute(
        "cancer",
        rna_type="miRNA",
        organism="Homo sapiens"
    )

    if result.success:
        print(f"Found {result.data['num_results']} microRNAs")
        for entry in result.data['entries'][:5]:
            print(f"- {entry['rnacentral_id']}: {entry['description'][:60]}")

asyncio.run(find_targets())
```

### 3. Get Database Cross-References (1 minute)

```python
async def get_xrefs():
    adapter = RNAcentralAdapter()

    result = await adapter.execute(
        "URS0000527F89_9606",
        get_xrefs=True
    )

    if result.success:
        print(f"Cross-references for {result.data['rnacentral_id']}:")
        for xref in result.data['cross_references'][:10]:
            print(f"  {xref['database']}: {xref['accession']}")

asyncio.run(get_xrefs())
```

## Common Use Cases

### Antisense Drug Design

```python
# Find target miRNA and get its sequence
result = await adapter.execute("URS0000527F89_9606")

if result.success:
    target_sequence = result.data['sequence']
    # Use sequence to design complementary antisense oligonucleotide
    antisense = reverse_complement(target_sequence)  # Your design function
    print(f"Target: {target_sequence}")
    print(f"ASO: {antisense}")
```

### RNA Biomarker Discovery

```python
# Survey all miRNAs in a disease context
result = await adapter.execute("Alzheimer", rna_type="miRNA")

if result.success:
    # Filter by length (typical mature miRNA: 18-25 nt)
    candidates = [
        e for e in result.data['entries']
        if 18 <= e['length'] <= 25
    ]
    print(f"Found {len(candidates)} candidate biomarkers")
```

### Target Validation Across Species

```python
# Check if RNA is conserved (important for target validation)
human_result = await adapter.execute("MALAT1", organism="Homo sapiens")
mouse_result = await adapter.execute("MALAT1", organism="Mus musculus")

if human_result.success and mouse_result.success:
    print(f"Human hits: {human_result.data['num_results']}")
    print(f"Mouse hits: {mouse_result.data['num_results']}")
    print("RNA is conserved!" if both found)
```

## Input Formats

| Input Type | Example | What It Does |
|------------|---------|--------------|
| RNA ID | `URS0000000001` | Direct lookup |
| RNA ID + Species | `URS0000000001_9606` | Species-specific |
| Keyword | `"cancer"` | Search all fields |
| Gene Name | `"HOTAIR"` | Find by gene name |

## Filter Options

```python
result = await adapter.execute(
    "query",
    rna_type="miRNA",      # Filter by RNA type
    organism="Homo sapiens", # Filter by organism
    get_xrefs=True         # Include cross-references
)
```

**Available RNA Types:**
- `miRNA` - MicroRNA (gene regulation)
- `lncRNA` - Long non-coding RNA (regulation)
- `rRNA` - Ribosomal RNA (translation)
- `tRNA` - Transfer RNA (translation)
- `snoRNA` - Small nucleolar RNA (modification)
- `piRNA` - PIWI-interacting RNA (transposon silencing)
- `ribozyme` - Catalytic RNA

## Output Structure

```python
{
    # For ID lookups:
    "rnacentral_id": "URS0000000001",
    "sequence": "ACGUACGU...",
    "length": 1234,
    "rna_type": "lncRNA",
    "description": "Long non-coding RNA...",
    "species": [...],
    "url": "https://rnacentral.org/...",

    # For searches:
    "num_results": 42,
    "rna_types": [{"type": "miRNA", "count": 25}],
    "entries": [...]
}
```

## Error Handling

```python
result = await adapter.execute("invalid_query")

if not result.success:
    print(f"Error: {result.error}")
    # Handle error appropriately
else:
    # Process data
    data = result.data
```

## Performance Tips

1. **Use RNA IDs when possible** - 10x faster than keyword search
2. **Enable caching** - Automatic, repeats are instant
3. **Add filters** - Reduces result size and processing time
4. **Batch queries** - Group similar requests together

## Real-World Workflow

```python
async def rna_drug_pipeline(disease_term):
    """Complete RNA target identification pipeline"""
    adapter = RNAcentralAdapter()

    # Step 1: Find disease-associated miRNAs
    print(f"Searching for {disease_term}-associated RNAs...")
    result = await adapter.execute(
        disease_term,
        rna_type="miRNA",
        organism="Homo sapiens"
    )

    if not result.success:
        return None

    # Step 2: Filter by ideal length for therapeutics
    candidates = [
        e for e in result.data['entries']
        if 18 <= e['length'] <= 25
    ]

    print(f"Found {len(candidates)} therapeutic candidates")

    # Step 3: Get detailed info for top candidates
    targets = []
    for candidate in candidates[:5]:
        detail = await adapter.execute(
            candidate['rnacentral_id'],
            get_xrefs=True
        )

        if detail.success:
            targets.append({
                'id': detail.data['rnacentral_id'],
                'sequence': detail.data['sequence'],
                'length': detail.data['length'],
                'description': detail.data['description'],
                'xrefs': len(detail.data.get('cross_references', []))
            })

    return targets

# Run the pipeline
targets = asyncio.run(rna_drug_pipeline("cancer"))

for t in targets:
    print(f"\nTarget: {t['id']}")
    print(f"  Length: {t['length']} nt")
    print(f"  Databases: {t['xrefs']} cross-refs")
```

## Integration with PharmForge

### In a Pipeline

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()

# Add RNAcentral step
pipeline.add_step("rnacentral", {
    "rna_type": "miRNA",
    "organism": "Homo sapiens"
})

# Execute
result = await pipeline.execute("Alzheimer biomarkers")
```

### With Ranking

```python
# Rank RNA targets by multiple criteria
results = await adapter.execute("cancer", rna_type="lncRNA")

if results.success:
    entries = results.data['entries']

    # Score by length and annotation quality
    scored = [
        {
            **e,
            'score': (
                (1.0 if e['length'] > 200 else 0.5) +  # lncRNAs are long
                (0.5 if 'cancer' in e['description'].lower() else 0)
            )
        }
        for e in entries
    ]

    # Sort by score
    ranked = sorted(scored, key=lambda x: x['score'], reverse=True)
```

## Next Steps

1. **Read the full README** - `adapters/rnacentral/README.md`
2. **Run examples** - `python adapters/rnacentral/example_usage.py`
3. **Explore RNAcentral** - https://rnacentral.org/
4. **Check API docs** - https://rnacentral.org/api

## Common Issues

**Q: No results found**
- Check RNA ID format (starts with "URS")
- Try broader search terms
- Remove filters

**Q: Timeout errors**
- Increase timeout: `adapter.config["timeout"] = 120`
- Try during off-peak hours
- Check network connection

**Q: Which RNA ID to use?**
- RNAcentral IDs: `URS0000000001` (database-wide)
- With species: `URS0000000001_9606` (human-specific)
- Search by gene name first if unsure

## Support Resources

- [RNAcentral Help](https://rnacentral.org/help)
- [RNA Drug Discovery Review](https://www.nature.com/articles/nrd.2016.117)
- [PharmForge Documentation](../../docs/)

---

**Ready to start?** Copy one of the examples above and run it!

**Need help?** Check the full README or file an issue on GitHub.
