# SAbDab Adapter

Structural Antibody Database (SAbDab) adapter for retrieving antibody structures, sequences, and annotations.

## Overview

SAbDab is a comprehensive database of all antibody structures available in the Protein Data Bank (PDB), curated and maintained by the Oxford Protein Informatics Group (OPIG). This adapter provides programmatic access to:

- **Antibody Structures**: All antibody-containing PDB entries
- **CDR Definitions**: Complementarity-Determining Region annotations
- **Antigen Information**: Antibody-antigen complex data
- **Sequence Data**: Heavy and light chain sequences
- **Quality Metrics**: Resolution, R-factor, experimental method

## Features

- ✅ Search antibodies by antigen name
- ✅ Retrieve antibody data by PDB ID
- ✅ Filter by antibody type (IgG, Fab, scFv, VHH, etc.)
- ✅ Filter by species (human, mouse, etc.)
- ✅ Quality filtering (resolution, R-factor)
- ✅ CDR sequence extraction
- ✅ Antibody classification information
- ✅ Chain assignment (heavy, light, antigen)

## Installation

No additional dependencies beyond PharmForge core requirements.

```bash
# Already included in PharmForge
pip install -r requirements.txt
```

## Usage

### Basic Usage - PDB ID Lookup

```python
from adapters.sabdab import SAbDabAdapter

# Initialize adapter
adapter = SAbDabAdapter()

# Get antibody information by PDB ID
result = await adapter.execute("7BWJ")

if result.success:
    data = result.data
    print(f"PDB ID: {data['pdb_id']}")
    print(f"Resolution: {data['antibody_data']['resolution']} Å")
    print(f"Antigen: {data['antibody_data']['antigen']}")
    print(f"Antibody Type: {data['antibody_data']['antibody_type']}")
    print(f"Species: {data['antibody_data']['species']}")
```

### Search by Antigen

```python
# Search for COVID-19 spike protein antibodies
search_params = {
    "antigen": "spike",
    "species": "human",
    "resolution_max": 3.0  # Only high-quality structures
}

result = await adapter.execute(search_params, max_results=50)

if result.success:
    antibodies = result.data['antibodies']
    print(f"Found {len(antibodies)} antibodies")

    for ab in antibodies[:5]:
        print(f"\nPDB: {ab['pdb_id']}")
        print(f"  Resolution: {ab['resolution']} Å")
        print(f"  Type: {ab['antibody_type']}")
        print(f"  Antigen: {ab['antigen']}")
```

### Filter by Antibody Type

```python
# Search for nanobodies (VHH) against a target
search_params = {
    "ab_type": "VHH",  # Nanobody/single-domain antibody
    "antigen": "SARS-CoV-2",
    "resolution_max": 2.5
}

result = await adapter.execute(search_params)

if result.success:
    nanobodies = result.data['antibodies']
    print(f"Found {len(nanobodies)} nanobodies")
```

### Get CDR Sequences

```python
# Retrieve antibody with CDR sequences
result = await adapter.execute("7BWJ", include_cdrs=True)

if result.success:
    cdrs = result.data['antibody_data'].get('cdrs', {})

    # Heavy chain CDRs
    heavy = cdrs.get('heavy_chain', {})
    print("Heavy Chain CDRs:")
    print(f"  H1: {heavy.get('h1')}")
    print(f"  H2: {heavy.get('h2')}")
    print(f"  H3: {heavy.get('h3')}")

    # Light chain CDRs
    light = cdrs.get('light_chain', {})
    print("Light Chain CDRs:")
    print(f"  L1: {light.get('l1')}")
    print(f"  L2: {light.get('l2')}")
    print(f"  L3: {light.get('l3')}")
```

### Quality Filtering

```python
# Find high-quality antibody structures
search_params = {
    "antigen": "Her2",
    "species": "human",
    "resolution_max": 2.0,  # High resolution
    "rfactor_max": 0.25     # Good R-factor
}

result = await adapter.execute(search_params)

if result.success:
    high_quality = result.data['antibodies']
    print(f"Found {len(high_quality)} high-quality structures")
```

## API Reference

### Input Types

The adapter accepts two types of input:

#### 1. PDB ID String
```python
pdb_id = "7BWJ"  # Direct PDB lookup
```

#### 2. Search Parameters Dictionary
```python
search_params = {
    "antigen": str,           # Antigen name (e.g., "spike", "Her2")
    "pdb_id": str,           # Specific PDB ID
    "species": str,          # Species (e.g., "human", "mouse")
    "ab_type": str,          # Antibody type (e.g., "IgG", "Fab", "VHH")
    "resolution_max": float, # Maximum resolution in Å
    "rfactor_max": float     # Maximum R-factor
}
```

### Output Format

#### PDB Lookup Result
```python
{
    "pdb_id": "7BWJ",
    "antibody_data": {
        "pdb_id": "7BWJ",
        "resolution": 2.45,
        "method": "X-ray diffraction",
        "rfactor": 0.23,
        "antibody_type": "Fab",
        "antigen": "SARS-CoV-2 spike",
        "heavy_chain": "H",
        "light_chain": "L",
        "antigen_chain": "A",
        "species": "Homo sapiens",
        "scfv": False,
        "engineered": True,
        "cdrs": {
            "heavy_chain": {"h1": ..., "h2": ..., "h3": ...},
            "light_chain": {"l1": ..., "l2": ..., "l3": ...}
        }
    },
    "source": "SAbDab",
    "reference": "http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/"
}
```

#### Search Result
```python
{
    "search_params": {...},
    "total_results": 42,
    "antibodies": [
        {
            "pdb_id": "7BWJ",
            "resolution": 2.45,
            "method": "X-ray diffraction",
            "antibody_type": "Fab",
            "antigen": "SARS-CoV-2 spike",
            "species": "Homo sapiens",
            "heavy_chain": "H",
            "light_chain": "L"
        },
        ...
    ],
    "source": "SAbDab"
}
```

## Antibody Types

Supported antibody types:

- **IgG**: Immunoglobulin G (full antibody)
- **IgA**: Immunoglobulin A
- **IgM**: Immunoglobulin M
- **IgE**: Immunoglobulin E
- **IgD**: Immunoglobulin D
- **Fab**: Fragment antigen-binding
- **scFv**: Single-chain variable fragment
- **VHH**: Single-domain antibody (nanobody)

## Quality Metrics

### Resolution
- **< 2.0 Å**: High quality, atomic details visible
- **2.0-3.0 Å**: Medium quality, good for most analyses
- **> 3.0 Å**: Lower quality, limited detail

### R-factor
- **< 0.25**: Good quality refinement
- **0.25-0.30**: Acceptable
- **> 0.30**: Poor quality

### Experimental Methods
- **X-ray diffraction**: Most common, typically high resolution
- **Cryo-EM**: Increasingly common, good for large complexes
- **NMR**: Solution structures, typically lower resolution

## Use Cases

### 1. Therapeutic Antibody Discovery
```python
# Find all human antibodies against a cancer target
search_params = {
    "antigen": "Her2",
    "species": "human",
    "ab_type": "IgG",
    "resolution_max": 3.0
}
result = await adapter.execute(search_params)
```

### 2. Nanobody Engineering
```python
# Find VHH nanobodies with high-quality structures
search_params = {
    "ab_type": "VHH",
    "resolution_max": 2.5,
    "rfactor_max": 0.25
}
result = await adapter.execute(search_params)
```

### 3. Structure-Based Design
```python
# Get specific antibody structure for modeling
result = await adapter.execute("7BWJ", include_cdrs=True)
# Use CDR sequences for design
```

### 4. Database Mining
```python
# Find all COVID-19 antibodies
search_params = {
    "antigen": "SARS-CoV-2"
}
result = await adapter.execute(search_params, max_results=1000)
```

## Configuration

```python
adapter = SAbDabAdapter(
    config={
        "cache_dir": "./cache/sabdab",  # Cache directory
        "timeout": 60,                   # Request timeout (seconds)
        "max_results": 100,              # Max search results
        "rate_limit_delay": 0.5          # Delay between requests
    }
)
```

## Caching

Results are automatically cached using PharmForge's caching system:

- **Cache key**: Generated from input parameters
- **TTL**: 24 hours (default)
- **Storage**: Redis (hot) + disk (warm)

## Error Handling

```python
result = await adapter.execute("INVALID")

if not result.success:
    print(f"Error: {result.error}")
    # Handle error appropriately
```

## Limitations

1. **API Availability**: Requires internet connection to SAbDab server
2. **Rate Limiting**: Respect 0.5s delay between requests
3. **Data Updates**: Database updated weekly
4. **CDR Parsing**: Full CDR sequence extraction requires structure file parsing

## References

- **Database**: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/
- **Paper**: Dunbar et al. (2014) "SAbDab: the structural antibody database" *Nucleic Acids Research*
- **DOI**: 10.1093/nar/gkt1043
- **OPIG Group**: http://opig.stats.ox.ac.uk/

## Support

For issues specific to this adapter:
- Check PharmForge documentation
- Review example usage above
- Check SAbDab website for database status

For SAbDab database questions:
- Visit: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/about/
- Contact: OPIG group at Oxford

## License

This adapter is part of PharmForge (MIT License).
SAbDab database is freely available for academic use.
