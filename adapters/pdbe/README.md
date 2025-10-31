# PDBe Adapter for PharmForge

European mirror of the Protein Data Bank with additional annotations and faster access for European users.

## Overview

The **Protein Data Bank in Europe (PDBe)** is the European node of the worldwide Protein Data Bank. It provides the same experimental structures as RCSB PDB but may have:
- Additional European annotations
- Faster access for European users
- FunPDBe functional annotations
- Alternative access for redundancy

## Installation

No additional dependencies required - uses standard Python libraries:

```bash
pip install aiohttp
```

## Features

- Download protein structures by PDB ID
- Search structures by keywords, organism, or compound
- Extract ligand and binding site information
- Get quality metrics (resolution, R-factors)
- Fetch functional annotations
- Support for multiple file formats (PDB, mmCIF)
- Automatic caching for performance
- Rate limiting to respect API limits

## Usage

### Basic Usage - Retrieve Structure by PDB ID

```python
import asyncio
from adapters.pdbe.adapter import PDBEAdapter

async def main():
    adapter = PDBEAdapter()

    # Fetch structure for COVID-19 main protease
    result = await adapter.execute("6LU7")

    if result.success:
        data = result.data
        print(f"PDB ID: {data['pdb_id']}")
        print(f"Title: {data['quality_metrics']['title']}")
        print(f"Method: {data['experimental_method']}")
        print(f"Resolution: {data['resolution']} Å")
        print(f"Organism: {data['organism']}")
        print(f"Number of ligands: {len(data['ligands'])}")
        print(f"Number of binding sites: {len(data['binding_sites'])}")

        # Access structure file
        if 'pdb' in data['structure_files']:
            pdb_content = data['structure_files']['pdb']
            print(f"PDB file size: {len(pdb_content)} bytes")

asyncio.run(main())
```

### Search for Structures

```python
async def search_example():
    adapter = PDBEAdapter()

    # Search for kinase structures
    query = {
        "query_type": "keyword",
        "query": "kinase",
        "max_results": 5
    }

    result = await adapter.execute(query)

    if result.success:
        for entry in result.data['entries']:
            print(f"{entry['pdb_id']}: {entry['title']}")
            print(f"  Method: {entry['method']}, Resolution: {entry['resolution']} Å")

asyncio.run(search_example())
```

### Get Ligand Information

```python
async def ligand_example():
    adapter = PDBEAdapter(config={
        "include_ligands": True,
        "include_binding_sites": True
    })

    result = await adapter.execute("1HSG")  # HIV-1 protease

    if result.success:
        data = result.data

        print("Ligands:")
        for ligand in data['ligands']:
            print(f"  {ligand['id']}: {ligand['name']}")
            print(f"    Formula: {ligand['formula']}")
            print(f"    Weight: {ligand['weight']}")

        print("\nBinding Sites:")
        for site in data['binding_sites']:
            print(f"  Site {site['site_id']}: {site['num_residues']} residues")
            print(f"    Ligand: {site['ligand_id']}")

asyncio.run(ligand_example())
```

### Download Structure Files

```python
async def download_example():
    adapter = PDBEAdapter(config={
        "download_pdb": True,
        "download_cif": True,
        "cache_dir": "./my_structures"
    })

    result = await adapter.execute("3CL0")

    if result.success:
        data = result.data

        # Files are automatically cached
        print(f"PDB file: {data['file_paths']['pdb']}")
        print(f"mmCIF file: {data['file_paths']['cif']}")

asyncio.run(download_example())
```

### Get Functional Annotations

```python
async def annotations_example():
    adapter = PDBEAdapter(config={
        "include_annotations": True
    })

    result = await adapter.execute("1A2B")

    if result.success:
        data = result.data

        if data['annotations']:
            print("Functional Annotations:")
            print(data['annotations'])

asyncio.run(annotations_example())
```

## Configuration Options

```python
config = {
    "cache_dir": "./cache/pdbe",           # Cache directory (default: ./cache/pdbe)
    "download_pdb": True,                  # Download PDB format (default: True)
    "download_cif": False,                 # Download mmCIF format (default: False)
    "include_ligands": True,               # Extract ligand info (default: True)
    "include_binding_sites": True,         # Extract binding sites (default: True)
    "include_annotations": False,          # Fetch functional annotations (default: False)
    "timeout": 30,                         # Request timeout in seconds (default: 30)
    "rate_limit_delay": 0.5                # Delay between requests (default: 0.5s)
}

adapter = PDBEAdapter(config=config)
```

## Input Formats

### PDB ID (String)

```python
result = await adapter.execute("1ABC")
```

### Query Dictionary

```python
query = {
    "query_type": "pdb_id",        # or "keyword", "organism", "compound"
    "query": "1A2B",               # PDB ID or search term
    "max_results": 10              # Maximum number of results (for searches)
}

result = await adapter.execute(query)
```

## Output Format

### For PDB ID Queries

```python
{
    "pdb_id": "6LU7",
    "summary": {...},                      # Full PDBe summary data
    "structure_files": {
        "pdb": "PDB file content...",
        "cif": "mmCIF file content..."
    },
    "file_paths": {
        "pdb": "/path/to/cached/pdb6lu7.ent",
        "cif": "/path/to/cached/6lu7.cif"
    },
    "quality_metrics": {
        "experimental_method": "X-RAY DIFFRACTION",
        "resolution": 2.16,
        "r_value": 0.195,
        "r_free": 0.235,
        "deposit_date": "2020-01-28",
        "release_date": "2020-02-05",
        "title": "Structure of COVID-19 main protease..."
    },
    "molecules": {...},                    # Molecule/chain information
    "ligands": [
        {
            "id": "PRD",
            "name": "N-[(5-METHYLISOXAZOL-3-YL)CARBONYL]ALANYL-L-VALYL-N~1~-((1R,2Z)-4-(BENZYLOXY)-4-OXO-1-{[(3R)-2-OXOPYRROLIDIN-3-YL]METHYL}BUT-2-ENYL)-L-LEUCINAMIDE",
            "formula": "C33 H45 N5 O7",
            "type": "BOUND",
            "weight": 623.747
        }
    ],
    "binding_sites": [
        {
            "site_id": "AC1",
            "ligand_id": "PRD",
            "residues": [...],
            "num_residues": 15
        }
    ],
    "annotations": {...},                  # Functional annotations (if requested)
    "experimental_method": "X-RAY DIFFRACTION",
    "resolution": 2.16,
    "organism": "Severe acute respiratory syndrome coronavirus 2",
    "num_chains": 2,
    "num_residues": 306,
    "reference": "Protein Data Bank in Europe (PDBe)"
}
```

### For Search Queries

```python
{
    "query_type": "keyword",
    "query": "kinase",
    "entries": [
        {
            "pdb_id": "1ATP",
            "title": "Structure of human protein kinase...",
            "organism": "Homo sapiens",
            "method": "X-RAY DIFFRACTION",
            "resolution": 2.0
        },
        ...
    ],
    "num_results": 5,
    "total_available": 1500
}
```

## API Endpoints

The adapter uses the following PDBe API endpoints:

- **Summary**: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}`
- **Molecules**: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}`
- **Ligands**: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id}`
- **Binding Sites**: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb_id}`
- **Annotations**: `https://www.ebi.ac.uk/pdbe/api/graph-api/pdb/funpdbe_annotation/{pdb_id}`
- **Search**: `https://www.ebi.ac.uk/pdbe/api/search/pdb/select`
- **Files**: `https://www.ebi.ac.uk/pdbe/entry-files/download/`

## Quality Metrics

### Resolution

- **High quality**: < 2.0 Å
- **Medium quality**: 2.0-3.0 Å
- **Low quality**: > 3.0 Å

### R-factors

- **R-value (good)**: < 0.25
- **R-free (good)**: < 0.28

Lower values indicate better agreement between the model and experimental data.

## Use Cases

### 1. Protein Structure Retrieval for Docking

```python
# Get protein structure for molecular docking
adapter = PDBEAdapter(config={"download_pdb": True})
result = await adapter.execute("6LU7")

if result.success:
    pdb_file = result.data['file_paths']['pdb']
    # Use pdb_file with AutoDock Vina adapter
```

### 2. Ligand-Protein Complex Analysis

```python
# Analyze existing ligand-protein complexes
adapter = PDBEAdapter(config={
    "include_ligands": True,
    "include_binding_sites": True
})

result = await adapter.execute("1HSG")  # HIV-1 protease with inhibitor

if result.success:
    for ligand in result.data['ligands']:
        print(f"Ligand: {ligand['name']}")
        print(f"Formula: {ligand['formula']}")
```

### 3. Alternative to RCSB for Redundancy

```python
# Use PDBe as backup if RCSB is unavailable
try:
    rcsb_result = await rcsb_adapter.execute("1ABC")
except Exception as e:
    # Fallback to PDBe
    pdbe_result = await pdbe_adapter.execute("1ABC")
```

### 4. European Data Center for Faster Access

```python
# Faster access for European users
adapter = PDBEAdapter()  # Uses European servers
result = await adapter.execute("3CL0")
```

### 5. Cross-Validation with RCSB Data

```python
# Compare data from both sources
rcsb_result = await rcsb_adapter.execute("6LU7")
pdbe_result = await pdbe_adapter.execute("6LU7")

# Verify consistency
assert rcsb_result.data['resolution'] == pdbe_result.data['resolution']
```

## Error Handling

```python
async def safe_retrieval(pdb_id: str):
    adapter = PDBEAdapter()
    result = await adapter.execute(pdb_id)

    if result.success:
        return result.data
    else:
        print(f"Error: {result.error}")
        return None
```

## Caching

The adapter automatically caches:
- Structure files (PDB, mmCIF)
- API responses (via PharmForge cache system)

Cached files are stored in the configured cache directory (default: `./cache/pdbe`).

## Rate Limiting

The adapter implements rate limiting to respect PDBe API limits:
- Default delay: 0.5 seconds between requests
- Configurable via `rate_limit_delay` parameter

## Comparison with RCSB PDB

| Feature | PDBe | RCSB PDB |
|---------|------|----------|
| Data Coverage | Same | Same |
| Location | Europe (UK) | USA |
| API Style | REST API | REST API |
| Annotations | FunPDBe included | Different annotations |
| Search | Solr-based | GraphQL available |
| Best For | European users, functional annotations | US users, advanced search |

## API Reference

### Main Methods

- `execute(input_data, **params)` - Main execution method
- `validate_input(input_data)` - Validate input format
- `generate_cache_key(input_data, **kwargs)` - Generate cache key
- `get_metadata()` - Get adapter metadata

### Internal Methods

- `_execute_pdb_id(pdb_id, **params)` - Execute PDB ID query
- `_execute_search(query_dict, **params)` - Execute search query
- `_get_summary(pdb_id)` - Get structure summary
- `_get_molecules(pdb_id)` - Get molecule information
- `_get_ligands(pdb_id)` - Get ligand information
- `_get_binding_sites(pdb_id)` - Get binding sites
- `_get_annotations(pdb_id)` - Get functional annotations
- `_download_structure(pdb_id, format)` - Download structure file

## Resources

- **Website**: https://www.ebi.ac.uk/pdbe/
- **API Documentation**: https://www.ebi.ac.uk/pdbe/api/doc/
- **Organization**: EMBL-EBI (European Bioinformatics Institute)
- **Citation**: Protein Data Bank in Europe (PDBe)

## License

Free API - No authentication required for most endpoints.

## Support

For issues specific to this adapter, please file an issue on the PharmForge GitHub repository.

For PDBe API questions, consult the official documentation: https://www.ebi.ac.uk/pdbe/api/doc/
