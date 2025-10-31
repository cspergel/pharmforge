# COCONUT Adapter - Quick Start Guide

## 5-Minute Setup

### Import and Initialize
```python
import asyncio
from adapters.coconut.adapter import COCONUTAdapter

adapter = COCONUTAdapter()
```

### Simple SMILES Search
```python
async def search_caffeine():
    # Search for natural products similar to caffeine
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    result = await adapter(smiles)

    if result.success:
        print(f"Found {result.data['num_results']} natural products")
        for compound in result.data['results']:
            print(f"  {compound['name']} from {compound['natural_source']['organism']}")

asyncio.run(search_caffeine())
```

### Search by Organism
```python
async def find_streptomyces_compounds():
    query = {
        "query_type": "organism",
        "query": "Streptomyces",
        "max_results": 20
    }
    result = await adapter(query)

    if result.success:
        print(f"Found {result.data['num_results']} compounds from Streptomyces")

asyncio.run(find_streptomyces_compounds())
```

### Search with Property Filters
```python
async def find_drug_like_natural_products():
    query = {
        "query_type": "smiles",
        "query": "CC(C)=CCC1C(=O)CCC1C",  # Terpene scaffold
        "filters": {
            "molecular_weight": {"min": 200, "max": 500},
            "logp": {"min": 0, "max": 5}
        },
        "max_results": 50
    }
    result = await adapter(query)

    if result.success:
        # Filter for Lipinski's Rule of Five
        drug_like = [
            c for c in result.data['results']
            if c['molecular_weight'] <= 500
            and c['properties']['logp'] <= 5
            and c['properties']['hbd'] <= 5
            and c['properties']['hba'] <= 10
        ]
        print(f"Found {len(drug_like)} drug-like natural products")

asyncio.run(find_drug_like_natural_products())
```

## Common Query Types

### 1. SMILES Search
```python
{"query_type": "smiles", "query": "CC(C)O"}
```

### 2. Name Search
```python
{"query_type": "name", "query": "Taxol"}
```

### 3. Organism Search
```python
{"query_type": "organism", "query": "Streptomyces"}
```

### 4. ID Lookup
```python
{"query_type": "id", "query": "CNP0123456"}
```

## Property Filters

Available filters:
- `molecular_weight`: `{"min": 200, "max": 500}`
- `logp`: `{"min": 0, "max": 5}`
- `organism`: `"Streptomyces"` (substring match)

## Result Fields

Every result contains:
```python
{
    "coconut_id": "CNP0123456",
    "name": "Compound Name",
    "smiles": "CC(C)O",
    "molecular_weight": 250.0,
    "natural_source": {
        "organism": "Taxus brevifolia",
        "common_name": "Pacific Yew",
        "taxonomy": "Plantae;Pinophyta;..."
    },
    "properties": {
        "logp": 2.5,
        "hba": 3,
        "hbd": 1,
        "rotatable_bonds": 2
    },
    "biological_activities": ["Anticancer"],
    "references": ["PMID:12345678"]
}
```

## Integration with Pipeline

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()

# Find natural products
pipeline.add_step("coconut", {
    "query_type": "organism",
    "query": "Streptomyces",
    "filters": {"molecular_weight": {"min": 300, "max": 700}},
    "max_results": 100
})

# Predict ADMET
pipeline.add_step("tdc_admet", {
    "properties": ["Caco2_Wang", "Solubility_AqSolDB"]
})

# Dock to target
pipeline.add_step("vina_docking", {
    "receptor": "protein.pdb"
})

results = await pipeline.execute()
```

## Error Handling

```python
result = await adapter(query)

if result.success:
    # Process results
    for compound in result.data['results']:
        print(compound['name'])
else:
    # Handle error
    print(f"Error: {result.error}")

    # Check warnings
    if result.data and result.data.get('warnings'):
        for warning in result.data['warnings']:
            print(f"Warning: {warning}")
```

## Tips

1. **Start broad, filter narrow**: Search by organism first, then apply property filters
2. **Cache results**: Enable caching with `use_cache=True` (default)
3. **Batch processing**: Process multiple compounds with property filters
4. **Respect rate limits**: Default 0.5s delay between requests
5. **Validate structures**: Cross-check SMILES with RDKit before using

## Full Examples

See `example_usage.py` for 7 complete examples including:
- Natural product library screening
- Scaffold mining from nature
- Taxonomic analysis
- Traditional medicine validation

## Documentation

- **Full README**: `adapters/coconut/README.md`
- **Tests**: `backend/tests/test_coconut_adapter.py`
- **Protocol**: `backend/core/adapters/protocol.py`

## Quick Links

- COCONUT Database: https://coconut.naturalproducts.net/
- API Documentation: https://coconut.naturalproducts.net/api/
- PharmForge Docs: See main README

---

**Ready to use!** No API key required, no additional installation needed.
