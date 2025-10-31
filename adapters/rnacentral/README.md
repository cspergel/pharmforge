# RNAcentral Adapter

Query RNA sequences and annotations from RNAcentral for RNA-targeting drug discovery.

## Overview

**RNAcentral** is a comprehensive database of non-coding RNA sequences from a broad range of organisms. This adapter provides access to:

- RNA sequence data (miRNA, lncRNA, rRNA, tRNA, etc.)
- Functional annotations
- Cross-references to other databases
- Species and taxonomic information
- RNA structural information

## Features

- **Search by RNA ID**: Direct lookup using RNAcentral URS identifiers
- **Keyword Search**: Find RNAs by gene name, function, or description
- **Filter by RNA Type**: Focus on specific RNA classes (miRNA, lncRNA, etc.)
- **Filter by Organism**: Target specific species
- **Cross-References**: Get links to related databases
- **Batch Support**: Process multiple queries with automatic caching

## Installation

The adapter uses only standard HTTP libraries and requires no API key:

```bash
pip install aiohttp
```

## Usage

### Basic RNA ID Lookup

```python
from adapters.rnacentral.adapter import RNAcentralAdapter

adapter = RNAcentralAdapter()

# Look up a specific RNA by ID
result = await adapter.execute("URS0000000001")

if result.success:
    rna_data = result.data
    print(f"RNA ID: {rna_data['rnacentral_id']}")
    print(f"Type: {rna_data['rna_type']}")
    print(f"Length: {rna_data['length']}")
    print(f"Sequence: {rna_data['sequence']}")
```

### Keyword Search

```python
# Search for microRNAs related to cancer
result = await adapter.execute(
    "cancer",
    rna_type="miRNA",
    organism="Homo sapiens"
)

if result.success:
    print(f"Found {result.data['num_results']} results")

    # Show RNA types distribution
    for rna_type in result.data['rna_types']:
        print(f"  {rna_type['type']}: {rna_type['count']}")

    # Show top results
    for entry in result.data['entries'][:5]:
        print(f"\n{entry['rnacentral_id']}")
        print(f"  {entry['description']}")
        print(f"  Type: {entry['rna_type']}, Length: {entry['length']}")
```

### Get Cross-References

```python
# Get database cross-references for an RNA
result = await adapter.execute(
    "URS0000000001",
    get_xrefs=True
)

if result.success and "cross_references" in result.data:
    print("Cross-references:")
    for xref in result.data['cross_references'][:10]:
        print(f"  {xref.get('database')}: {xref.get('accession')}")
```

### Filter by RNA Type

Supported RNA types include:

- `miRNA` - microRNA
- `lncRNA` - long non-coding RNA
- `rRNA` - ribosomal RNA
- `tRNA` - transfer RNA
- `snRNA` - small nuclear RNA
- `snoRNA` - small nucleolar RNA
- `piRNA` - PIWI-interacting RNA
- `ribozyme` - catalytic RNA

```python
# Search for long non-coding RNAs
result = await adapter.execute(
    "gene regulation",
    rna_type="lncRNA"
)
```

## Input Formats

### RNAcentral ID Format

- Format: `URS` + 10 digits (e.g., `URS0000000001`)
- Case-insensitive
- Automatically detected by adapter

### Search Queries

- Gene names (e.g., "HOTAIR", "MALAT1")
- Functions (e.g., "regulation", "cancer")
- Organisms (e.g., "human", "mouse")
- Any descriptive text

## Output Format

### For ID Lookups

```json
{
  "rnacentral_id": "URS0000000001",
  "sequence": "ACGUACGUACGU...",
  "length": 1234,
  "rna_type": "lncRNA",
  "description": "Long non-coding RNA involved in...",
  "species": [
    {
      "taxid": 9606,
      "name": "Homo sapiens",
      "common_name": "human"
    }
  ],
  "url": "https://rnacentral.org/rna/URS0000000001",
  "cross_references": [...],  // If get_xrefs=True
  "raw_data": {...}  // Full API response
}
```

### For Searches

```json
{
  "num_results": 42,
  "rna_types": [
    {"type": "miRNA", "count": 25},
    {"type": "lncRNA", "count": 17}
  ],
  "organisms": [
    {"name": "Homo sapiens", "count": 30}
  ],
  "entries": [
    {
      "rnacentral_id": "URS0000000001",
      "description": "...",
      "rna_type": "miRNA",
      "length": 22,
      "url": "https://rnacentral.org/rna/URS0000000001"
    }
  ],
  "raw_results": [...]  // Full search results
}
```

## Use Cases for Drug Discovery

### 1. Target Identification

Find RNAs involved in disease pathways:

```python
# Find RNAs associated with Alzheimer's disease
result = await adapter.execute("Alzheimer", organism="Homo sapiens")
```

### 2. Antisense Oligonucleotide Design

Get sequences for antisense drug design:

```python
# Get specific miRNA sequence for ASO design
result = await adapter.execute("URS00001C5F1A", get_xrefs=True)
sequence = result.data['sequence']
```

### 3. RNA-Targeting Small Molecules

Identify RNAs with structural features:

```python
# Search for riboswitches (RNA drug targets)
result = await adapter.execute("riboswitch", rna_type="ribozyme")
```

### 4. Cross-Species Validation

Find conserved RNAs across species:

```python
# Search for specific RNA across organisms
results_human = await adapter.execute("MALAT1", organism="Homo sapiens")
results_mouse = await adapter.execute("MALAT1", organism="Mus musculus")
```

## Configuration

The adapter supports these configuration options:

```python
adapter = RNAcentralAdapter()
adapter.config = {
    "rate_limit_delay": 0.3,  # Seconds between requests
    "timeout": 60,            # Request timeout in seconds
    "max_results": 100        # Maximum results per query
}
```

## Error Handling

```python
result = await adapter.execute("INVALID_ID")

if not result.success:
    print(f"Error: {result.error}")
    # Handle the error appropriately
```

Common errors:

- `Invalid input` - Input is not a string or is empty
- `Failed to retrieve RNA with ID: {id}` - RNA ID not found
- `Failed to search RNAcentral` - Search query failed

## API Rate Limits

RNAcentral has no published rate limits, but the adapter implements:

- Default 0.3 second delay between requests
- Automatic caching to minimize API calls
- Configurable timeout and retry handling

## Caching

The adapter automatically caches:

- RNA ID lookups
- Search queries (with all parameters)
- Cross-reference data

Cache keys include:

- Adapter name and version
- Input query
- All filter parameters (rna_type, organism, etc.)

## Performance Tips

1. **Use specific RNA IDs when possible** - Direct lookups are faster than searches
2. **Enable caching** - Repeated queries return instantly from cache
3. **Filter searches** - Use rna_type and organism filters to reduce result size
4. **Batch similar queries** - Group queries by organism or type

## Integration Examples

### Pipeline Integration

```python
from backend.core.pipeline import Pipeline

# Create a pipeline that finds target RNAs
pipeline = Pipeline()
pipeline.add_step("rnacentral", {
    "rna_type": "miRNA",
    "organism": "Homo sapiens"
})

result = await pipeline.execute("cancer biomarkers")
```

### Ranking RNA Targets

```python
# Rank RNAs by length and functional relevance
results = await adapter.execute("oncogene", rna_type="lncRNA")

if results.success:
    entries = results.data['entries']

    # Sort by length (longer lncRNAs often more functional)
    ranked = sorted(entries, key=lambda x: x['length'], reverse=True)

    for rna in ranked[:10]:
        print(f"{rna['rnacentral_id']}: {rna['length']} nt")
```

## Related Resources

- [RNAcentral Website](https://rnacentral.org/)
- [RNAcentral API Documentation](https://rnacentral.org/api)
- [RNAcentral Help](https://rnacentral.org/help)
- [RNA Biology Review](https://rnacentral.org/expert-databases)

## References

1. RNAcentral Consortium (2021). "RNAcentral 2021: secondary structure integration, improved sequence search and new member databases." *Nucleic Acids Research*.

2. For RNA therapeutics review: Matsui & Corey (2017). "Non-coding RNAs as drug targets." *Nature Reviews Drug Discovery*.

## Technical Details

- **Adapter Type**: `api`
- **Base URL**: `https://rnacentral.org/api/v1`
- **Authentication**: None required
- **Rate Limiting**: 0.3s default delay
- **Caching**: Automatic via AdapterProtocol
- **Version**: 1.0.0

## Troubleshooting

### No results found

- Check RNA ID format (should start with "URS")
- Try broader search terms
- Remove filters to see all results

### Timeout errors

- Increase timeout in config: `adapter.config["timeout"] = 120`
- Check network connection
- Try query during off-peak hours

### Invalid RNA ID

- Verify ID format: `URS` + 10 digits
- Check ID exists on [RNAcentral website](https://rnacentral.org/)
- Use search if unsure of exact ID

## Future Enhancements

Planned features:

- [ ] RNA structure prediction integration
- [ ] Batch ID lookup endpoint
- [ ] FASTA sequence export
- [ ] Advanced filtering (e.g., by length, GO terms)
- [ ] Integration with RNA folding tools

## Support

For issues or questions:

- Check [RNAcentral documentation](https://rnacentral.org/help)
- Review PharmForge adapter protocol
- File issues on PharmForge GitHub

---

**License**: MIT
**Maintainer**: PharmForge Team
**Last Updated**: 2025-10-31
