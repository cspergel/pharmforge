# UniProt Adapter

Adapter for the UniProt REST API, providing access to comprehensive protein sequence and annotation data from the Universal Protein Resource.

## Overview

UniProt (Universal Protein Resource) is the world's leading high-quality, comprehensive resource for protein sequence and annotation data. This adapter queries the UniProt REST API to retrieve:

- Protein sequences and properties
- Protein functions and annotations
- Post-translational modifications (PTMs)
- Structural information (domains, active sites, binding sites)
- Disease associations
- Cross-references to structure databases (PDB, AlphaFold)
- Subcellular localization

## API Documentation

- Base URL: `https://rest.uniprot.org`
- API Docs: https://www.uniprot.org/help/api
- Rate Limit: ~10 requests/second (self-imposed)
- No authentication required

## Installation

No additional dependencies required beyond the base adapter requirements:
- `aiohttp` for async HTTP requests
- `asyncio` for async/await support

## Usage

### Basic Example

```python
import asyncio
from adapters.uniprot import UniProtAdapter

async def main():
    adapter = UniProtAdapter()

    # Query by gene name
    result = await adapter("BRCA2")

    if result.success:
        data = result.data
        print(f"Protein: {data['protein_name']}")
        print(f"Gene: {', '.join(data['gene_names'])}")
        print(f"Organism: {data['organism']}")
        print(f"Length: {data['sequence_length']} amino acids")
        print(f"MW: {data['molecular_weight']/1000:.1f} kDa")

        print(f"\nFunction:")
        for func in data['function'][:2]:
            print(f"  {func[:200]}...")

        print(f"\nDomains: {len(data['domains'])}")
        for domain in data['domains'][:5]:
            print(f"  - {domain}")

        if data['pdb_structures']:
            print(f"\nPDB Structures: {', '.join(data['pdb_structures'][:5])}")

    else:
        print(f"Error: {result.error}")

asyncio.run(main())
```

### Query by UniProt Accession

```python
# Direct accession query
result = await adapter("P51587")  # BRCA2_HUMAN
```

### Query Without Full Sequence

```python
# Get protein info without sequence (faster, smaller response)
result = await adapter("TP53", include_sequence=False)
```

### Search for Multiple Proteins

```python
# Search returns multiple results
from adapters.uniprot.adapter import UniProtAdapter

adapter = UniProtAdapter()
search_results = await adapter._search_protein("kinase", max_results=10)
print(f"Found {len(search_results)} kinases")
```

## Input Formats

The adapter accepts:
- **Gene names**: e.g., "BRCA2", "TP53", "EGFR"
- **UniProt accessions**: e.g., "P51587", "P04637"
- **Protein names**: e.g., "insulin", "hemoglobin"

## Output Format

### AdapterResult Structure

```python
{
    "success": bool,
    "data": {
        "accession": str,                    # UniProt accession (e.g., "P51587")
        "entry_name": str,                   # Entry name (e.g., "BRCA2_HUMAN")
        "gene_names": List[str],             # Gene symbols
        "protein_name": str,                 # Recommended protein name
        "organism": str,                     # Scientific name
        "sequence": str,                     # Amino acid sequence (or None)
        "sequence_length": int,              # Sequence length
        "molecular_weight": float,           # Molecular weight in Da
        "function": List[str],               # Functional descriptions
        "subcellular_location": List[str],   # Where protein is located
        "disease_associations": List[str],   # Associated diseases
        "post_translational_modifications": List[str],  # PTMs
        "domains": List[str],                # Protein domains
        "active_sites": List[str],           # Catalytic residues
        "binding_sites": List[str],          # Ligand binding sites
        "modifications": List[str],          # Modified residues
        "pdb_structures": List[str],         # PDB IDs
        "alphafold_structures": List[str],   # AlphaFold IDs
        "keywords": List[str]                # Annotation keywords
    },
    "error": Optional[str],
    "cache_hit": bool,
    "metadata": {
        "source": "uniprot",
        "accession": str,
        "original_query": str,
        "adapter_version": str
    }
}
```

## Configuration

Default configuration:
```python
{
    "rate_limit_delay": 0.1,    # 100ms between requests
    "timeout": 60,               # 60 second timeout
    "format": "json"             # Response format
}
```

## Features

### Protein Sequence

```python
result = await adapter("INS")  # Insulin

if result.success:
    seq = result.data['sequence']
    print(f"Sequence length: {len(seq)} aa")
    print(f"Sequence: {seq}")
```

### Protein Function

```python
result = await adapter("EGFR")

if result.success:
    print("Functions:")
    for func in result.data['function']:
        print(f"  - {func}")
```

### Post-Translational Modifications

```python
result = await adapter("P04637")  # TP53

if result.success:
    print("PTMs:")
    for ptm in result.data['post_translational_modifications']:
        print(f"  - {ptm}")

    print("\nModifications:")
    for mod in result.data['modifications'][:10]:
        print(f"  - {mod}")
```

### Structural Information

```python
result = await adapter("BRCA1")

if result.success:
    print("Domains:")
    for domain in result.data['domains']:
        print(f"  - {domain}")

    print("\nPDB Structures:")
    for pdb in result.data['pdb_structures']:
        print(f"  - {pdb}")

    print("\nAlphaFold Models:")
    for af in result.data['alphafold_structures']:
        print(f"  - {af}")
```

### Disease Associations

```python
result = await adapter("APOE")

if result.success:
    print("Disease associations:")
    for disease in result.data['disease_associations']:
        print(f"  - {disease}")
```

## Error Handling

The adapter handles:
- Invalid protein IDs/names (returns error)
- API timeouts (60 second timeout)
- Rate limiting (100ms delay between requests)
- Network errors (returns error with details)
- Protein not found (returns error)

## Caching

Results are automatically cached using the AdapterProtocol caching system:
- Cache key: SHA256 hash of (adapter name, version, input, parameters)
- Cache duration: Configurable in backend cache settings
- Disable cache: `await adapter(protein_id, use_cache=False)`

## Examples

### Analyzing Enzyme Active Sites

```python
result = await adapter("P00698")  # Lysozyme

if result.success:
    print(f"Protein: {result.data['protein_name']}")
    print(f"\nActive sites:")
    for site in result.data['active_sites']:
        print(f"  - {site}")

    print(f"\nBinding sites:")
    for site in result.data['binding_sites']:
        print(f"  - {site}")
```

### Finding Drug Targets

```python
# Look for membrane proteins with known structures
result = await adapter("EGFR")

if result.success:
    data = result.data

    # Check if it's a good drug target
    is_membrane = any("membrane" in loc.lower() for loc in data['subcellular_location'])
    has_structure = len(data['pdb_structures']) > 0
    has_binding_sites = len(data['binding_sites']) > 0

    print(f"Membrane protein: {is_membrane}")
    print(f"Has 3D structure: {has_structure}")
    print(f"Has binding sites: {has_binding_sites}")

    if is_membrane and has_structure and has_binding_sites:
        print("\nGood drug target candidate!")
```

### Sequence Analysis

```python
result = await adapter("INS")  # Insulin

if result.success:
    seq = result.data['sequence']

    # Calculate basic properties
    cys_count = seq.count('C')
    aromatic_count = sum(seq.count(aa) for aa in 'FYW')
    charged_count = sum(seq.count(aa) for aa in 'DEKR')

    print(f"Length: {len(seq)} aa")
    print(f"Cysteine residues: {cys_count}")
    print(f"Aromatic residues: {aromatic_count}")
    print(f"Charged residues: {charged_count}")

    # Check for disulfide bonds
    disulfides = [mod for mod in result.data['modifications'] if 'Disulfide' in mod]
    print(f"Disulfide bonds: {len(disulfides)}")
```

### Comparing Protein Isoforms

```python
# Get different isoforms
result1 = await adapter("P51587")  # BRCA2 canonical
result2 = await adapter("P51587-2")  # BRCA2 isoform 2

if result1.success and result2.success:
    len1 = result1.data['sequence_length']
    len2 = result2.data['sequence_length']

    print(f"Canonical: {len1} aa")
    print(f"Isoform 2: {len2} aa")
    print(f"Difference: {abs(len1 - len2)} aa")
```

## Limitations

- No authentication required (public API)
- Rate limited to ~10 requests/second (self-imposed)
- Large proteins may have long response times
- Some annotations may be incomplete
- Historical data limited to current version

## Data Quality

UniProt entries are classified by review status:
- **Swiss-Prot**: Manually reviewed, high quality
- **TrEMBL**: Computationally analyzed, not reviewed

This adapter retrieves data from both sources but prioritizes Swiss-Prot entries.

## Cross-References

The adapter extracts cross-references to:
- **PDB**: Experimental 3D structures
- **AlphaFoldDB**: Predicted 3D structures
- Other databases available in raw data

## Version History

- **1.0.0** (2025-10-25): Initial release
  - Protein information retrieval
  - Sequence extraction
  - Function annotations
  - Post-translational modifications
  - Structural information
  - Disease associations
  - Caching support
  - Rate limiting
