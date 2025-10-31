# GEO (Gene Expression Omnibus) Adapter

## Overview

The GEO adapter provides programmatic access to NCBI's Gene Expression Omnibus (GEO), the largest public repository of gene expression data. It uses the free NCBI E-utilities API to search and retrieve gene expression datasets without requiring authentication or API keys.

## Features

- **Multiple Query Types:**
  - Query by GEO accession (GSE, GDS, GPL, GSM)
  - Search by gene name/symbol
  - Search by disease/condition/keyword
  - Complex boolean queries (AND, OR, NOT)

- **Database Support:**
  - GEO DataSets (`gds`): Curated gene expression datasets
  - GEO Profiles (`geoprofiles`): Individual gene expression profiles

- **Filtering Capabilities:**
  - Filter by organism
  - Filter by entry type
  - Filter by platform
  - Filter by date range

- **Rich Metadata:**
  - Dataset titles and summaries
  - Sample counts
  - Platform information
  - PubMed links
  - Download URLs for raw data

## Installation

No additional dependencies required beyond the base PharmForge requirements.

## Usage

### Basic Gene Query

```python
from adapters.geo.adapter import GEOAdapter

# Initialize adapter
adapter = GEOAdapter(email="your_email@example.com")

# Search for datasets related to a gene
result = await adapter.execute("TP53")

if result.success:
    print(f"Found {result.data['total_count']} datasets")
    for dataset in result.data['datasets']:
        print(f"{dataset['accession']}: {dataset['title']}")
```

### Query by GEO Accession

```python
# Fetch specific dataset by accession
result = await adapter.execute("GSE1234")

if result.success:
    dataset = result.data['dataset']
    print(f"Title: {dataset['title']}")
    print(f"Samples: {dataset['n_samples']}")
    print(f"Download URL: {dataset['download_urls']['matrix']}")
```

### Disease/Condition Query with Filters

```python
# Search for breast cancer datasets in humans
result = await adapter.execute({
    "query": "breast cancer",
    "max_results": 20,
    "filters": {
        "organism": "Homo sapiens",
        "entry_type": "GSE"
    }
})

if result.success:
    summary = result.data['summary']
    print(f"Platform types: {summary['platform_types']}")
    print(f"Total datasets: {summary['total_datasets']}")
```

### Advanced Boolean Query

```python
# Complex query with boolean operators
result = await adapter.execute({
    "query": "diabetes AND (insulin OR glucose)",
    "max_results": 15,
    "filters": {
        "organism": "Homo sapiens",
        "date_from": "2020/01/01",
        "date_to": "2025/12/31"
    }
})
```

### Search GEO Profiles

```python
# Search gene expression profiles instead of datasets
result = await adapter.execute({
    "query": "BRCA1",
    "database": "geoprofiles",
    "max_results": 10
})

if result.success:
    for profile in result.data['datasets']:
        print(f"{profile['genesymbol']}: {profile['title']}")
```

## API Endpoints Used

The adapter uses the following NCBI E-utilities endpoints:

- **ESearch**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi`
  - Searches GEO databases and returns UIDs

- **ESummary**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi`
  - Retrieves dataset summaries and metadata

- **FTP Downloads**: `https://ftp.ncbi.nlm.nih.gov/geo/`
  - Direct download URLs for dataset files

## GEO Accession Types

- **GSE**: Series (a collection of related samples)
- **GDS**: Dataset (curated dataset with experimental design)
- **GPL**: Platform (array or sequencing platform)
- **GSM**: Sample (individual biological sample)

## Query Parameters

### Main Parameters

- `query` (str or dict): Search query or GEO accession
- `database` (str): Database to search (`gds` or `geoprofiles`)
- `max_results` (int): Maximum number of results to return (default: 100)
- `fetch_details` (bool): Whether to fetch full dataset details (default: True)

### Filter Options

- `organism` (str): Organism name (e.g., "Homo sapiens", "Mus musculus")
- `entry_type` (str): Entry type (e.g., "GSE", "GDS")
- `platform` (str): Platform identifier
- `date_from` (str): Start date in YYYY/MM/DD format
- `date_to` (str): End date in YYYY/MM/DD format

## Response Format

### Successful Response

```python
{
    "success": True,
    "data": {
        "query": "TP53",
        "database": "gds",
        "total_count": 1523,
        "returned_count": 100,
        "uids": ["200001234", "200001235", ...],
        "datasets": [
            {
                "uid": "200001234",
                "accession": "GDS1234",
                "title": "Gene expression in TP53 mutant cells",
                "summary": "This study examines...",
                "gpl": "GPL570",
                "gse": "GSE5678",
                "taxon": "Homo sapiens",
                "n_samples": 24,
                "pubmedids": [12345678],
                "download_urls": {
                    "soft": "https://ftp.ncbi.nlm.nih.gov/geo/...",
                    "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/...",
                    "supplementary": "https://ftp.ncbi.nlm.nih.gov/geo/..."
                }
            },
            ...
        ],
        "summary": {
            "total_datasets": 100,
            "datasets_with_pubmed": 85,
            "organisms": [
                {"organism": "Homo sapiens", "count": 75},
                {"organism": "Mus musculus", "count": 25}
            ],
            "platform_types": [
                {"type": "Expression profiling by array", "count": 60},
                {"type": "Expression profiling by high throughput sequencing", "count": 40}
            ],
            "entry_types": [
                {"type": "GSE", "count": 100}
            ]
        }
    },
    "metadata": {
        "source": "geo",
        "query": "TP53",
        "adapter_version": "1.0.0",
        "total_results": 1523,
        "datasets_fetched": 100
    }
}
```

## Rate Limits

- **Without API Key**: 3 requests per second
- **With API Key**: 10 requests per second

To use an API key:

```python
adapter = GEOAdapter(
    email="your_email@example.com",
    api_key="your_ncbi_api_key"
)
```

Get a free API key from: https://www.ncbi.nlm.nih.gov/account/

## Testing

Run the test suite:

```bash
cd adapters/geo
python test_geo_adapter.py
```

The test suite includes:
1. Query by GEO accession
2. Query by gene name
3. Query by disease/condition
4. Advanced query with filters
5. GEO Profiles search
6. Error handling

## Example Queries That Work

### 1. Cancer Research
```python
await adapter.execute({
    "query": "lung cancer AND smoking",
    "filters": {"organism": "Homo sapiens"}
})
```

### 2. Drug Response
```python
await adapter.execute({
    "query": "aspirin AND platelets",
    "max_results": 50
})
```

### 3. Gene Expression in Disease
```python
await adapter.execute({
    "query": "Alzheimer disease AND APOE",
    "filters": {
        "organism": "Homo sapiens",
        "date_from": "2020/01/01"
    }
})
```

### 4. Specific Dataset
```python
await adapter.execute("GSE123456")
```

### 5. RNA-seq Studies
```python
await adapter.execute({
    "query": "RNA-seq[Filter]",
    "max_results": 100
})
```

## Limitations

- **Rate Limits**: Respect NCBI's rate limits (3 req/s without key, 10 req/s with key)
- **Result Size**: Maximum 100,000 results per query (practical limit ~10,000)
- **Timeout**: Queries timeout after 30 seconds by default
- **Download**: The adapter provides URLs but does not download large data files
- **GPL/GSM Accessions**: Direct fetch by GPL or GSM accessions requires search-based retrieval

## Best Practices

1. **Always provide email**: NCBI requests an email for tracking purposes
2. **Use API key**: For production use, obtain and use an NCBI API key
3. **Cache results**: The adapter supports caching to reduce API calls
4. **Specific queries**: More specific queries return more relevant results
5. **Filter appropriately**: Use organism and date filters to narrow results
6. **Rate limiting**: Add delays between rapid successive queries

## Integration with PharmForge

The GEO adapter follows the PharmForge adapter protocol and can be used with:

- **Workflow Builder**: Chain with other adapters for multi-step analysis
- **Cache System**: Automatic caching of query results
- **Task Manager**: Integrate into automated pipelines
- **Data Export**: Results can be exported to various formats

## Common Use Cases

1. **Target Validation**: Find expression data for drug targets
2. **Disease Mechanism**: Explore gene expression in disease models
3. **Drug Response**: Identify datasets studying drug effects
4. **Biomarker Discovery**: Search for expression signatures
5. **Cross-Study Analysis**: Aggregate data across multiple studies

## Resources

- [GEO Homepage](https://www.ncbi.nlm.nih.gov/geo/)
- [E-utilities Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [GEO Programmatic Access](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html)
- [NCBI API Key](https://www.ncbi.nlm.nih.gov/account/)

## Support

For issues or questions:
- Check NCBI E-utilities documentation
- Review test_geo_adapter.py for examples
- Consult PharmForge adapter documentation
