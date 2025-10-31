# GEO Adapter Quick Start Guide

## Installation
No additional dependencies needed - uses standard PharmForge libraries.

## Basic Usage

```python
from adapters.geo.adapter import GEOAdapter

# Initialize
adapter = GEOAdapter(email="your_email@example.com")

# Run a query
result = await adapter.execute("TP53")

# Check results
if result.success:
    print(f"Found {result.data['total_count']} datasets")
```

## Common Query Patterns

### 1. Search by Gene
```python
result = await adapter.execute("BRCA1")
```

### 2. Search by Disease
```python
result = await adapter.execute({
    "query": "diabetes",
    "filters": {"organism": "Homo sapiens"}
})
```

### 3. Get Specific Dataset
```python
result = await adapter.execute("GSE1234")
```

### 4. Complex Query
```python
result = await adapter.execute({
    "query": "breast cancer AND TP53",
    "max_results": 20,
    "filters": {
        "organism": "Homo sapiens",
        "entry_type": "GSE",
        "date_from": "2020/01/01"
    }
})
```

## Working with Results

```python
result = await adapter.execute("TP53")

if result.success:
    # Get basic info
    total = result.data['total_count']
    datasets = result.data['datasets']

    # Iterate through datasets
    for ds in datasets:
        print(f"{ds['accession']}: {ds['title']}")
        print(f"Samples: {ds['n_samples']}")
        print(f"Organism: {ds['taxon']}")

        # Get download URLs
        print(f"Download: {ds['download_urls']['matrix']}")
```

## Quick Examples

### Example 1: Find Cancer Datasets
```python
adapter = GEOAdapter(email="test@example.com")
result = await adapter.execute({
    "query": "lung cancer",
    "max_results": 10
})
```

### Example 2: Get Gene Expression Profiles
```python
result = await adapter.execute({
    "query": "APOE",
    "database": "geoprofiles"
})
```

### Example 3: Drug Response Studies
```python
result = await adapter.execute({
    "query": "doxorubicin AND cancer",
    "filters": {"organism": "Homo sapiens"}
})
```

## Error Handling

```python
result = await adapter.execute("GSE999999999")

if not result.success:
    print(f"Error: {result.error}")
else:
    print(f"Success!")
```

## Rate Limiting

The adapter automatically handles rate limiting:
- 3 requests/second without API key
- 10 requests/second with API key

To use API key:
```python
adapter = GEOAdapter(
    email="your_email@example.com",
    api_key="your_ncbi_api_key"
)
```

Get API key: https://www.ncbi.nlm.nih.gov/account/

## Testing

Run tests:
```bash
python test_geo_adapter.py
```

Run examples:
```bash
python example_usage.py
```

## More Information

- Full documentation: `README.md`
- Implementation details: `IMPLEMENTATION_SUMMARY.md`
- Test cases: `test_geo_adapter.py`
- Usage examples: `example_usage.py`
