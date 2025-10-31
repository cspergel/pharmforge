# OpenTargets Adapter

Adapter for the OpenTargets Platform GraphQL API, providing access to target-disease associations, drug mechanisms, and genetic evidence.

## Overview

The OpenTargets Platform integrates genetic, genomic, transcriptomic, and chemical data to identify and prioritize targets for drug discovery. This adapter queries the OpenTargets GraphQL API to retrieve:

- Target gene information and annotations
- Target-disease associations with evidence scores
- Known drugs and their mechanisms of action
- Genetic evidence linking genes to diseases

## API Documentation

- Base URL: `https://api.platform.opentargets.org/api/v4/graphql` (Updated as of v1.2.0)
- API Docs: https://platform-docs.opentargets.org/data-access/graphql-api
- GraphQL Playground: https://api.platform.opentargets.org/api/v4/graphql/browser
- Rate Limit: ~10 requests/second (self-imposed)

**Note:** The old endpoint `https://platform-api.opentargets.io` has been deprecated. Use the new endpoint above.

## Installation

No additional dependencies required beyond the base adapter requirements:
- `aiohttp` for async HTTP requests
- `asyncio` for async/await support

## Usage

### Basic Example

```python
import asyncio
from adapters.opentargets import OpenTargetsAdapter

async def main():
    adapter = OpenTargetsAdapter()

    # Query by gene symbol
    result = await adapter("BRCA2")

    if result.success:
        data = result.data
        print(f"Gene: {data['target_info']['symbol']}")
        print(f"Name: {data['target_info']['name']}")
        print(f"\nTop Disease Associations:")
        for assoc in data['disease_associations'][:5]:
            print(f"  - {assoc['disease_name']}: {assoc['overall_score']:.3f}")

        print(f"\nKnown Drugs:")
        for drug in data['drug_mechanisms'][:5]:
            print(f"  - {drug['drug_name']} (Phase {drug['max_clinical_phase']})")
    else:
        print(f"Error: {result.error}")

asyncio.run(main())
```

### Query by Ensembl Gene ID

```python
# Direct Ensembl ID query
result = await adapter("ENSG00000139618")  # BRCA2
```

### Specific Query Types

```python
# Get only target information
result = await adapter("BRCA2", query_type="target_info")

# Get only disease associations
result = await adapter("BRCA2", query_type="disease_associations", max_results=100)

# Get only drug mechanisms
result = await adapter("BRCA2", query_type="drug_mechanisms")
```

## Input Formats

The adapter accepts:
- **Gene symbols**: e.g., "BRCA2", "TP53", "EGFR"
- **Ensembl gene IDs**: e.g., "ENSG00000139618"

## Output Format

### AdapterResult Structure

```python
{
    "success": bool,
    "data": {
        "gene_id": str,              # Ensembl gene ID
        "query_type": str,           # Type of query executed
        "target_info": {             # Present if query_type includes "target_info"
            "symbol": str,
            "name": str,
            "biotype": str,
            "function": List[str],
            "go_terms": List[str]
        },
        "disease_associations": [    # Present if query_type includes "disease_associations"
            {
                "disease_id": str,
                "disease_name": str,
                "overall_score": float,
                "datatype_scores": {
                    "genetic_association": float,
                    "somatic_mutation": float,
                    "known_drug": float,
                    # ... other datatypes
                },
                "therapeutic_areas": List[str]
            }
        ],
        "drug_mechanisms": [         # Present if query_type includes "drug_mechanisms"
            {
                "drug_id": str,
                "drug_name": str,
                "drug_type": str,
                "mechanism": str,
                "phase": int,
                "max_clinical_phase": int,
                "status": str
            }
        ]
    },
    "error": Optional[str],
    "cache_hit": bool,
    "metadata": {
        "source": "opentargets",
        "gene_id": str,
        "original_query": str,
        "adapter_version": str
    }
}
```

## Configuration

Default configuration:
```python
{
    "rate_limit_delay": 0.1,   # 100ms between requests
    "timeout": 60,              # 60 second timeout
    "max_results": 50,          # Maximum results per query
    "max_retries": 3,           # Number of retry attempts (v1.2.0+)
    "retry_delay": 1.0,         # Initial retry delay in seconds (v1.2.0+)
    "base_url": None            # Custom API endpoint (v1.2.0+)
}
```

### Custom Configuration Example

```python
# Use custom configuration
adapter = OpenTargetsAdapter(config={
    "timeout": 120,
    "max_retries": 5,
    "retry_delay": 2.0,
    "max_results": 100
})
```

## Error Handling

The adapter handles:
- Invalid gene symbols/IDs (returns error)
- API timeouts (60 second timeout, configurable)
- Rate limiting (self-imposed 100ms delay)
- Network errors (automatic retry with exponential backoff - v1.2.0+)
- Connection errors (tries multiple endpoints - v1.2.0+)
- GraphQL errors (logged and returned as error)
- DNS resolution failures (tries fallback endpoints - v1.2.0+)

### Troubleshooting

If you experience connectivity issues:

1. **Check network connection:**
   ```bash
   curl https://api.platform.opentargets.org/api/v4/graphql
   ```

2. **Enable debug logging:**
   ```python
   import logging
   logging.basicConfig(level=logging.DEBUG)
   ```

3. **Increase timeout and retries:**
   ```python
   adapter = OpenTargetsAdapter(config={
       "timeout": 120,
       "max_retries": 5
   })
   ```

4. **Use HTTP fallback (if HTTPS blocked):**
   ```python
   adapter = OpenTargetsAdapter(config={
       "base_url": "http://api.platform.opentargets.org/api/v4/graphql"
   })
   ```

## Caching

Results are automatically cached using the AdapterProtocol caching system:
- Cache key: SHA256 hash of (adapter name, version, input, parameters)
- Cache duration: Configurable in backend cache settings
- Disable cache: `await adapter(gene_id, use_cache=False)`

## Examples

### Finding Drug Targets for a Disease

```python
# Find genes associated with breast cancer
result = await adapter("BRCA2", query_type="disease_associations")

if result.success:
    cancer_assocs = [
        assoc for assoc in result.data['disease_associations']
        if 'cancer' in assoc['disease_name'].lower()
    ]

    for assoc in sorted(cancer_assocs, key=lambda x: x['overall_score'], reverse=True)[:10]:
        print(f"{assoc['disease_name']}: {assoc['overall_score']:.3f}")
        print(f"  Evidence types: {list(assoc['datatype_scores'].keys())}")
```

### Finding Druggable Targets

```python
# Find targets with approved drugs
result = await adapter("EGFR", query_type="drug_mechanisms")

if result.success:
    approved_drugs = [
        drug for drug in result.data['drug_mechanisms']
        if drug['max_clinical_phase'] == 4  # Phase 4 = Approved
    ]

    print(f"Approved drugs targeting EGFR: {len(approved_drugs)}")
    for drug in approved_drugs:
        print(f"  - {drug['drug_name']}: {drug['mechanism']}")
```

## Limitations

- Requires internet connection
- Rate limited to ~10 requests/second (self-imposed)
- GraphQL queries may be complex and slow for large datasets
- Maximum 10,000 results per query (API limitation)
- Some data may require authentication (not supported in this version)

## Version History

- **1.2.0** (2025-10-25): Connectivity and reliability improvements
  - FIXED: Updated API endpoint from deprecated `platform-api.opentargets.io` to `api.platform.opentargets.org`
  - ADDED: Automatic retry logic with exponential backoff (3 retries by default)
  - ADDED: Fallback endpoint URLs for improved resilience
  - ADDED: Configurable endpoint URL, retry attempts, and delays
  - FIXED: GraphQL query schema issues (pagination, field names)
  - IMPROVED: Error messages with troubleshooting hints
  - IMPROVED: Connection handling with proper DNS and timeout settings
  - IMPROVED: Empty result handling (empty lists are now valid results)

- **1.0.0** (2025-10-25): Initial release
  - Target information queries
  - Disease association queries
  - Drug mechanism queries
  - Gene symbol search
  - Caching support
  - Rate limiting
