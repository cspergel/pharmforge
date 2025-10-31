# DisGeNET Adapter

Adapter for the DisGeNET REST API, providing access to gene-disease associations and variant-disease associations based on curated databases and literature.

## Overview

DisGeNET is a comprehensive database of gene-disease and variant-disease associations from various expert-curated databases and scientific literature. This adapter queries the DisGeNET REST API to retrieve:

- Gene-disease associations (GDAs)
- Variant-disease associations (VDAs)
- Evidence types and sources
- Association scores and literature references

## API Documentation

- Base URL: `https://www.disgenet.org/api`
- API Docs: https://www.disgenet.org/api/
- Authentication: Optional API key for full access (public endpoints available)
- Rate Limit: ~1 request/second (self-imposed for public API)

## Installation

No additional dependencies required beyond the base adapter requirements:
- `aiohttp` for async HTTP requests
- `asyncio` for async/await support

## API Key (Optional)

For full access to all DisGeNET features, obtain a free API key:
1. Sign up at https://www.disgenet.org/signup
2. Verify your email
3. Get your API key from your profile
4. Pass it to the adapter: `DisGeNETAdapter(api_key="your_key_here")`

**Note:** Some endpoints may require authentication. The public API has limited functionality.

## Usage

### Basic Example - Gene-Disease Associations

```python
import asyncio
from adapters.disgenet import DisGeNETAdapter

async def main():
    adapter = DisGeNETAdapter()  # Or DisGeNETAdapter(api_key="your_key")

    # Query gene-disease associations
    result = await adapter("BRCA2", query_type="gene")

    if result.success:
        data = result.data
        print(f"Gene: {data['query']}")
        print(f"Associations found: {data['summary']['num_associations']}")
        print(f"Average score: {data['summary']['avg_score']:.3f}")

        print(f"\nTop 5 Disease Associations:")
        for assoc in data['associations'][:5]:
            print(f"  - {assoc['disease_name']}")
            print(f"    Score: {assoc['score']}")
            print(f"    Evidence: {assoc['evidence_type']}")
            print(f"    Source: {assoc['source']}")
    else:
        print(f"Error: {result.error}")

asyncio.run(main())
```

### Disease-Gene Associations

```python
# Find genes associated with a disease
result = await adapter("alzheimer disease", query_type="disease", max_results=50)

if result.success:
    print(f"Genes associated with {result.data['query']}:")
    for assoc in result.data['associations'][:10]:
        print(f"  {assoc['gene_symbol']}: {assoc['score']}")
```

### Variant-Disease Associations

```python
# Query variant-disease associations
result = await adapter("rs429358", query_type="variant")

if result.success:
    for assoc in result.data['associations']:
        print(f"Disease: {assoc['disease_name']}")
        print(f"Variant: {assoc['variant_id']}")
        print(f"Score: {assoc['score']}")
```

## Input Formats

The adapter accepts:
- **Gene symbols**: e.g., "BRCA2", "TP53", "APOE"
- **Disease names**: e.g., "breast cancer", "alzheimer disease"
- **Variant IDs**: e.g., "rs429358", "rs7412"

## Output Format

### AdapterResult Structure

```python
{
    "success": bool,
    "data": {
        "query": str,                    # Original query
        "query_type": str,               # "gene", "disease", or "variant"
        "summary": {
            "num_associations": int,
            "sources": List[str],        # Data sources
            "evidence_types": List[str], # Types of evidence
            "avg_score": float,          # Average association score
            "max_score": float,          # Highest score
            "min_score": float           # Lowest score
        },
        "associations": [
            {
                "gene_symbol": str,
                "disease_name": str,
                "disease_id": str,       # Disease identifier (UMLS, OMIM, etc.)
                "score": float,          # Association score (0-1)
                "evidence_type": str,    # Type of evidence
                "source": str,           # Data source
                "pmid": str,             # PubMed ID(s)
                # For variant queries only:
                "variant_id": str,
                "chromosome": str,
                "position": int
            }
        ]
    },
    "error": Optional[str],
    "cache_hit": bool,
    "metadata": {
        "source": "disgenet",
        "query": str,
        "query_type": str,
        "adapter_version": str,
        "authenticated": bool            # Whether API key was used
    }
}
```

## Configuration

Default configuration:
```python
{
    "rate_limit_delay": 1.0,    # 1 second between requests
    "timeout": 60,               # 60 second timeout
    "max_results": 100           # Maximum results per query
}
```

## Query Types

### 1. Gene-Disease Associations (query_type="gene")

Find diseases associated with a gene:
```python
result = await adapter("TP53", query_type="gene", max_results=50)
```

### 2. Disease-Gene Associations (query_type="disease")

Find genes associated with a disease:
```python
result = await adapter("parkinson disease", query_type="disease", max_results=50)
```

**Note:** May require API key for full results.

### 3. Variant-Disease Associations (query_type="variant")

Find diseases associated with a genetic variant:
```python
result = await adapter("rs429358", query_type="variant")
```

## Evidence Types

DisGeNET aggregates evidence from multiple sources:
- **Literature**: Text mining from PubMed
- **Animal Models**: Model organism databases
- **Genetic Variation**: GWAS, mutation databases
- **Expert Curated**: Manual curation

## Association Scores

- Score range: 0.0 to 1.0
- Higher scores indicate stronger evidence
- Scores combine multiple evidence sources
- Consider both score and evidence type

## Error Handling

The adapter handles:
- Invalid gene/disease/variant names (returns empty results or error)
- API timeouts (60 second timeout)
- Rate limiting (1 second delay between requests)
- Authentication errors (returns error if API key required)
- Network errors (returns error with details)

## Caching

Results are automatically cached using the AdapterProtocol caching system:
- Cache key: SHA256 hash of (adapter name, version, input, parameters)
- Cache duration: Configurable in backend cache settings
- Disable cache: `await adapter(query, use_cache=False)`

## Examples

### Finding Cancer-Related Genes

```python
# Find genes associated with lung cancer
result = await adapter("lung cancer", query_type="disease", max_results=100)

if result.success:
    # Filter by high-confidence associations
    high_conf = [
        assoc for assoc in result.data['associations']
        if assoc['score'] > 0.5
    ]

    print(f"High-confidence genes for lung cancer: {len(high_conf)}")
    for assoc in sorted(high_conf, key=lambda x: x['score'], reverse=True)[:20]:
        print(f"  {assoc['gene_symbol']}: {assoc['score']:.3f} ({assoc['evidence_type']})")
```

### Analyzing Gene Pleiotropy

```python
# Find all diseases associated with a gene
result = await adapter("APOE", query_type="gene", max_results=200)

if result.success:
    diseases = result.data['associations']
    print(f"APOE is associated with {len(diseases)} diseases")

    # Group by evidence type
    by_evidence = {}
    for assoc in diseases:
        ev_type = assoc['evidence_type']
        if ev_type not in by_evidence:
            by_evidence[ev_type] = []
        by_evidence[ev_type].append(assoc)

    print("\nBy evidence type:")
    for ev_type, assocs in by_evidence.items():
        print(f"  {ev_type}: {len(assocs)} associations")
```

### Variant Pathogenicity Analysis

```python
# Analyze a variant's disease associations
result = await adapter("rs429358", query_type="variant")

if result.success:
    print(f"Variant {result.data['query']} associations:")
    for assoc in result.data['associations']:
        print(f"\nDisease: {assoc['disease_name']}")
        print(f"Score: {assoc['score']:.3f}")
        print(f"Evidence: {assoc['evidence_type']}")
        if assoc['pmid']:
            print(f"References: PMID:{assoc['pmid']}")
```

## Limitations

- **Authentication**: Some endpoints require API key
- **Rate Limiting**: Public API limited to ~1 request/second
- **Data Sources**: Limited to curated databases and literature
- **Disease Names**: Must match DisGeNET nomenclature
- **Gene Symbols**: HGNC standard symbols preferred
- **Variant IDs**: dbSNP rs IDs supported

## Data Sources

DisGeNET integrates data from:
- CURATED: UniProt, CGI, ClinGen, Orphanet, PsyGeNET
- ANIMAL MODELS: MGD, RGD
- LITERATURE: BEFREE (text mining)
- GWAS: GWAS Catalog, GWAS db
- VARIANTS: ClinVar, UniProt

## Version History

- **1.0.0** (2025-10-25): Initial release
  - Gene-disease associations
  - Disease-gene associations
  - Variant-disease associations
  - Evidence type filtering
  - API key support
  - Caching support
  - Rate limiting
