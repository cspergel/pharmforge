# Protein Structure Adapters for PharmForge

This directory contains three comprehensive adapters for accessing protein structure data from major databases. These adapters enable structure-based drug discovery workflows in PharmForge.

## Overview

### 1. AlphaFold DB Adapter (`adapters/alphafold/`)
Access AI-predicted protein structures from the AlphaFold Protein Structure Database.

**Key Features:**
- Download predicted structures by UniProt ID
- Access per-residue confidence scores (pLDDT)
- Predicted Aligned Error (PAE) data for domain analysis
- Structure quality metrics (pTM scores)
- Automatic file caching for offline access
- 200M+ predicted structures available

**Quality Metrics:**
- **pLDDT** (predicted LDDT): Per-residue confidence (0-100)
  - >90: Very high confidence
  - 70-90: Confident
  - <70: Low confidence
- **PAE**: Confidence in relative positions between residues
- **pTM**: Overall model confidence score

**Reference:** Jumper et al., Nature 2021, 596, 583-589

---

### 2. RCSB PDB Adapter (`adapters/rcsb_pdb/`)
Access experimental protein structures from the Research Collaboratory for Structural Bioinformatics Protein Data Bank.

**Key Features:**
- Download structures by PDB ID
- Full-text search by protein name or keywords
- Extract ligand information (including SMILES)
- Binding site identification
- Quality metrics (resolution, R-factors)
- Multiple file formats (PDB, mmCIF, PDBML)
- 200,000+ experimental structures

**Quality Metrics:**
- **Resolution**: Lower is better
  - <2.0 Å: High quality
  - 2.0-3.0 Å: Medium quality
  - >3.0 Å: Low quality
- **R-value/R-free**: Model quality (<0.25 is good)
- **Experimental methods**: X-ray, NMR, Cryo-EM, etc.

**Databases:** https://www.rcsb.org/

---

### 3. SWISS-MODEL Adapter (`adapters/swissmodel/`)
Access homology models from the SWISS-MODEL Repository.

**Key Features:**
- Homology model retrieval by UniProt ID
- Structure quality assessment (QMEAN, GMQE)
- Template structure information
- Sequence coverage analysis
- Quality filtering
- 800,000+ homology models

**Quality Metrics:**
- **QMEAN**: Quality Model Energy Analysis (-4 to 0, higher is better)
  - >-1.5: Very high quality
  - -2.5 to -1.5: High quality
  - <-2.5: Medium/Low quality
- **GMQE**: Global Model Quality Estimation (0-1, higher is better)
  - >0.7: High reliability
  - 0.4-0.7: Medium reliability
  - <0.4: Low reliability

**Reference:** Waterhouse et al., Nucleic Acids Res. 2018, 46, D296-D303

---

## Installation

All adapters are already included in PharmForge. Dependencies:

```bash
pip install aiohttp  # Async HTTP requests
```

No additional installation required!

---

## Quick Start

### AlphaFold: Get AI-predicted structure

```python
import asyncio
from adapters.alphafold import AlphaFoldAdapter

async def get_alphafold_structure():
    adapter = AlphaFoldAdapter(config={
        "download_pdb": True,
        "download_pae": True,
        "cache_dir": "./cache/alphafold"
    })

    # Get structure for P53 tumor suppressor
    result = await adapter.execute("P04637")

    if result.success:
        data = result.data
        print(f"Mean pLDDT: {data['plddt_scores']['mean']}")
        print(f"PDB file: {data['file_paths']['pdb']}")
        print(f"High confidence: {data['plddt_scores']['high_confidence_pct']}%")

    return result

asyncio.run(get_alphafold_structure())
```

### RCSB PDB: Get experimental structure

```python
import asyncio
from adapters.rcsb_pdb import RCSBPDBAdapter

async def get_experimental_structure():
    adapter = RCSBPDBAdapter(config={
        "download_pdb": True,
        "include_ligands": True,
        "include_binding_sites": True
    })

    # Get SARS-CoV-2 main protease structure
    result = await adapter.execute("6LU7")

    if result.success:
        data = result.data
        print(f"Resolution: {data['resolution']} Å")
        print(f"Ligands: {len(data['ligands'])}")
        print(f"Binding sites: {len(data['binding_sites'])}")

        # Get ligand SMILES
        for ligand in data['ligands']:
            if 'smiles' in ligand:
                print(f"Ligand: {ligand['name']}")
                print(f"SMILES: {ligand['smiles']}")

    return result

asyncio.run(get_experimental_structure())
```

### SWISS-MODEL: Get homology model

```python
import asyncio
from adapters.swissmodel import SwissModelAdapter

async def get_homology_model():
    adapter = SwissModelAdapter(config={
        "download_pdb": True,
        "min_qmean": -2.0,  # High quality only
        "min_gmqe": 0.5     # Reliable models only
    })

    # Get insulin homology model
    result = await adapter.execute("P01308")

    if result.success:
        data = result.data
        print(f"Models found: {data['model_count']}")
        print(f"Best QMEAN: {data['quality_metrics']['qmean']}")
        print(f"Best GMQE: {data['quality_metrics']['gmqe']}")
        print(f"Sequence identity: {data['quality_metrics']['sequence_identity']}%")

    return result

asyncio.run(get_homology_model())
```

---

## Advanced Usage

### Search for structures

```python
from adapters.rcsb_pdb import RCSBPDBAdapter

async def search_structures():
    adapter = RCSBPDBAdapter()

    # Search for COVID-19 related structures
    result = await adapter.search_structures(
        query="SARS-CoV-2 main protease",
        max_results=10
    )

    if result.success:
        pdb_ids = result.data['pdb_ids']
        print(f"Found {len(pdb_ids)} structures:")
        for pdb_id in pdb_ids:
            print(f"  - {pdb_id}")
```

### Compare predicted vs experimental structures

```python
async def compare_structures(uniprot_id, pdb_id):
    # Get AlphaFold prediction
    alphafold = AlphaFoldAdapter()
    af_result = await alphafold.execute(uniprot_id)

    # Get experimental structure
    pdb = RCSBPDBAdapter()
    pdb_result = await pdb.execute(pdb_id)

    if af_result.success and pdb_result.success:
        print(f"\nAlphaFold (predicted):")
        print(f"  Mean pLDDT: {af_result.data['plddt_scores']['mean']}")

        print(f"\nRCSB PDB (experimental):")
        print(f"  Resolution: {pdb_result.data['resolution']} Å")
        print(f"  R-value: {pdb_result.data['quality_metrics'].get('r_value', 'N/A')}")
```

### Batch processing

```python
async def batch_process_proteins(uniprot_ids):
    adapter = AlphaFoldAdapter()

    results = []
    for uniprot_id in uniprot_ids:
        result = await adapter.execute(uniprot_id)
        results.append((uniprot_id, result))

    # Filter by quality
    high_quality = [
        (uid, r) for uid, r in results
        if r.success and r.data['plddt_scores']['mean'] > 80
    ]

    print(f"High quality predictions: {len(high_quality)}/{len(uniprot_ids)}")
    return high_quality
```

---

## Choosing the Right Adapter

### Use AlphaFold when:
- You need structures for proteins without experimental data
- You want full-proteome coverage
- You need per-residue confidence scores
- You're working with recent proteins or isoforms
- Fast turnaround is important

### Use RCSB PDB when:
- You need experimentally validated structures
- You need ligand/binding site information
- You're working with drug targets
- You need the highest accuracy
- You need to cite experimental evidence

### Use SWISS-MODEL when:
- AlphaFold predictions aren't available
- You need template-based models
- You want to assess modeling feasibility
- You need intermediate quality between AlphaFold and experiments
- You need sequence similarity to templates

---

## Configuration Options

### AlphaFold Adapter

```python
config = {
    "cache_dir": "./cache/alphafold",     # Cache directory
    "download_pdb": True,                  # Download PDB format
    "download_cif": False,                 # Download mmCIF format
    "download_pae": False,                 # Download PAE data
    "timeout": 30                          # Request timeout (seconds)
}
```

### RCSB PDB Adapter

```python
config = {
    "cache_dir": "./cache/rcsb_pdb",       # Cache directory
    "download_pdb": True,                   # Download PDB format
    "download_cif": False,                  # Download mmCIF format
    "include_ligands": True,                # Extract ligand info
    "include_binding_sites": True,          # Extract binding sites
    "timeout": 30                           # Request timeout (seconds)
}
```

### SWISS-MODEL Adapter

```python
config = {
    "cache_dir": "./cache/swissmodel",     # Cache directory
    "download_pdb": True,                   # Download PDB models
    "include_templates": True,              # Include template info
    "min_qmean": -4.0,                     # Minimum QMEAN score
    "min_gmqe": 0.0,                       # Minimum GMQE score
    "timeout": 30                           # Request timeout (seconds)
}
```

---

## Testing

Each adapter includes comprehensive test suites:

```bash
# Test AlphaFold adapter
cd adapters/alphafold
python test_adapter.py

# Test RCSB PDB adapter
cd adapters/rcsb_pdb
python test_adapter.py

# Test SWISS-MODEL adapter
cd adapters/swissmodel
python test_adapter.py
```

---

## Caching

All adapters implement automatic file caching:

- **Downloaded structures** are cached to avoid re-downloading
- **API responses** are cached using PharmForge's global cache
- Cache keys are deterministic (same input = same cache key)
- Files are stored in configurable cache directories

**Benefits:**
- Faster repeated access
- Reduced API load
- Offline access to previously downloaded structures
- Significant speedup for batch processing

---

## Integration with PharmForge

These adapters integrate seamlessly with PharmForge workflows:

```python
from backend.core.adapters.protocol import registry

# Register adapters
from adapters.alphafold import AlphaFoldAdapter
from adapters.rcsb_pdb import RCSBPDBAdapter
from adapters.swissmodel import SwissModelAdapter

registry.register(AlphaFoldAdapter())
registry.register(RCSBPDBAdapter())
registry.register(SwissModelAdapter())

# Access from registry
alphafold = registry.get("alphafold")
result = await alphafold("P04637")
```

---

## Error Handling

All adapters follow PharmForge's `AdapterResult` pattern:

```python
result = await adapter.execute(input_data)

if result.success:
    # Process successful result
    data = result.data
else:
    # Handle error
    print(f"Error: {result.error}")

# Check cache hit
if result.cache_hit:
    print("Retrieved from cache")
```

---

## Performance Tips

1. **Enable caching**: Always configure cache directories
2. **Batch requests**: Group multiple requests to amortize overhead
3. **Filter early**: Use quality filters to reduce data transfer
4. **Choose format wisely**: PDB is smaller than mmCIF
5. **Download selectively**: Only download what you need (e.g., skip PAE if not needed)

---

## API Rate Limits

- **AlphaFold**: No documented rate limits, but be respectful
- **RCSB PDB**: No strict limits, but avoid excessive parallel requests
- **SWISS-MODEL**: Use responsibly; consider local caching

All adapters implement reasonable delays and caching to minimize API load.

---

## Examples

See test files for comprehensive examples:
- `adapters/alphafold/test_adapter.py`
- `adapters/rcsb_pdb/test_adapter.py`
- `adapters/swissmodel/test_adapter.py`

---

## Contributing

To add new protein structure sources:

1. Create new adapter directory
2. Implement `AdapterProtocol`
3. Add comprehensive tests
4. Update this README

---

## References

- **AlphaFold**: Jumper et al., Nature 2021, 596, 583-589
- **RCSB PDB**: Berman et al., Nucleic Acids Res. 2000, 28, 235-242
- **SWISS-MODEL**: Waterhouse et al., Nucleic Acids Res. 2018, 46, D296-D303

---

## License

These adapters are part of PharmForge and follow the same license.

---

## Support

For issues or questions:
1. Check test files for usage examples
2. Review adapter docstrings
3. Check PharmForge documentation
4. Open an issue on GitHub

---

## Future Enhancements

Planned improvements:
- [ ] Structural alignment between sources
- [ ] Quality comparison metrics
- [ ] Batch download optimization
- [ ] Additional file formats
- [ ] Integration with molecular dynamics
- [ ] Protein-ligand docking integration
