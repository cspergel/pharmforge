# Protein Structure Adapters - Quick Start Guide

## 🚀 Installation

```bash
# Already included in PharmForge!
# Only dependency: aiohttp (already installed)
```

---

## 📦 Three Adapters Available

| Adapter | Source | Coverage | Best For |
|---------|--------|----------|----------|
| **AlphaFold** | AI-predicted | 200M+ structures | Novel proteins, full proteome |
| **RCSB PDB** | Experimental | 200K+ structures | Drug targets, validated structures |
| **SWISS-MODEL** | Homology | 800K+ models | When AlphaFold unavailable |

---

## ⚡ Quick Examples

### 1. AlphaFold: Get AI prediction

```python
from adapters.alphafold import AlphaFoldAdapter
import asyncio

async def get_structure():
    adapter = AlphaFoldAdapter()
    result = await adapter.execute("P04637")  # TP53

    if result.success:
        print(f"Mean pLDDT: {result.data['plddt_scores']['mean']}")
        print(f"File: {result.data['file_paths']['pdb']}")

asyncio.run(get_structure())
```

### 2. RCSB PDB: Get experimental structure

```python
from adapters.rcsb_pdb import RCSBPDBAdapter
import asyncio

async def get_structure():
    adapter = RCSBPDBAdapter()
    result = await adapter.execute("6LU7")  # COVID-19 Mpro

    if result.success:
        print(f"Resolution: {result.data['resolution']} Å")
        print(f"Ligands: {len(result.data['ligands'])}")

asyncio.run(get_structure())
```

### 3. SWISS-MODEL: Get homology model

```python
from adapters.swissmodel import SwissModelAdapter
import asyncio

async def get_model():
    adapter = SwissModelAdapter()
    result = await adapter.execute("P01308")  # Insulin

    if result.success:
        print(f"QMEAN: {result.data['quality_metrics']['qmean']}")
        print(f"Models: {result.data['model_count']}")

asyncio.run(get_model())
```

---

## 🎯 Common Use Cases

### Search for structures

```python
adapter = RCSBPDBAdapter()
result = await adapter.search_structures("SARS-CoV-2 protease")

if result.success:
    print(f"Found: {result.data['pdb_ids']}")
```

### Download with quality filtering

```python
adapter = SwissModelAdapter(config={
    "min_qmean": -2.0,  # High quality only
    "min_gmqe": 0.5
})
result = await adapter.execute("P04637")
```

### Get multiple formats

```python
adapter = RCSBPDBAdapter(config={
    "download_pdb": True,
    "download_cif": True
})
result = await adapter.execute("1CRN")
```

### Extract ligand information

```python
adapter = RCSBPDBAdapter(config={
    "include_ligands": True,
    "include_binding_sites": True
})
result = await adapter.execute("6LU7")

if result.success:
    for ligand in result.data['ligands']:
        print(f"Ligand: {ligand['name']}")
        if 'smiles' in ligand:
            print(f"SMILES: {ligand['smiles']}")
```

---

## 📊 Quality Metrics

### AlphaFold pLDDT
- `>90` = Very high confidence ✅
- `70-90` = Confident ⚠️
- `<70` = Low confidence ❌

### RCSB Resolution
- `<2.0 Å` = High quality ✅
- `2.0-3.0 Å` = Medium quality ⚠️
- `>3.0 Å` = Low quality ❌

### SWISS-MODEL QMEAN
- `>-1.5` = Very high quality ✅
- `-2.5 to -1.5` = High quality ⚠️
- `<-2.5` = Medium/Low ❌

---

## ⚙️ Configuration

### AlphaFold
```python
config = {
    "download_pdb": True,      # Download PDB format
    "download_pae": True,      # Get PAE data
    "cache_dir": "./cache"     # Cache location
}
adapter = AlphaFoldAdapter(config=config)
```

### RCSB PDB
```python
config = {
    "download_pdb": True,          # PDB format
    "download_cif": False,         # mmCIF format
    "include_ligands": True,       # Extract ligands
    "include_binding_sites": True  # Extract binding sites
}
adapter = RCSBPDBAdapter(config=config)
```

### SWISS-MODEL
```python
config = {
    "download_pdb": True,    # Download models
    "min_qmean": -2.0,       # Quality filter
    "min_gmqe": 0.5          # Reliability filter
}
adapter = SwissModelAdapter(config=config)
```

---

## 🧪 Testing

```bash
# Test individual adapters
cd claude-code-agents-wizard-v2/adapters/alphafold
python test_adapter.py

cd ../rcsb_pdb
python test_adapter.py

cd ../swissmodel
python test_adapter.py

# Run integration examples
cd ..
python PROTEIN_STRUCTURE_INTEGRATION_EXAMPLE.py
```

---

## 📁 File Locations

```
claude-code-agents-wizard-v2/adapters/
├── alphafold/
│   ├── adapter.py           # AlphaFold adapter
│   ├── test_adapter.py      # Tests
│   └── __init__.py
├── rcsb_pdb/
│   ├── adapter.py           # RCSB PDB adapter
│   ├── test_adapter.py      # Tests
│   └── __init__.py
├── swissmodel/
│   ├── adapter.py           # SWISS-MODEL adapter
│   ├── test_adapter.py      # Tests
│   └── __init__.py
├── PROTEIN_STRUCTURE_ADAPTERS_README.md          # Full docs
├── PROTEIN_STRUCTURE_INTEGRATION_EXAMPLE.py      # Examples
└── PROTEIN_STRUCTURES_QUICK_START.md             # This file
```

---

## 🔍 Decision Tree

```
Do you know the structure ID?
│
├─ YES → Use RCSB PDB adapter
│         (experimental structures)
│
└─ NO → Do you have UniProt ID?
         │
         ├─ YES → Try AlphaFold first
         │         (AI-predicted, high coverage)
         │
         └─ NO → Search RCSB PDB by name
                  adapter.search_structures("protein name")
```

---

## 💡 Pro Tips

1. **Always check quality metrics** before using structures
2. **Use caching** for repeated access (significant speedup)
3. **Validate with multiple sources** when possible
4. **Experimental > AlphaFold > Homology** (generally)
5. **Check confidence regions** in AlphaFold predictions
6. **Look for ligands** in RCSB PDB for drug design

---

## 🆘 Common Issues

### "No structure found"
- Check UniProt/PDB ID format
- Try searching by protein name
- Check if structure exists in database

### "Low quality score"
- Consider using experimental structure
- Check alternative models
- Focus on high-confidence regions

### "Download failed"
- Check internet connection
- Verify API endpoints are accessible
- Check cache directory permissions

---

## 📚 More Information

- **Full Documentation**: `PROTEIN_STRUCTURE_ADAPTERS_README.md`
- **Examples**: `PROTEIN_STRUCTURE_INTEGRATION_EXAMPLE.py`
- **Summary**: `PROTEIN_STRUCTURE_ADAPTERS_SUMMARY.md` (in project root)

---

## 🎓 Example Workflow

```python
import asyncio
from adapters.alphafold import AlphaFoldAdapter
from adapters.rcsb_pdb import RCSBPDBAdapter

async def analyze_protein(uniprot_id, protein_name):
    # Step 1: Search for experimental structures
    pdb = RCSBPDBAdapter()
    search = await pdb.search_structures(protein_name, max_results=5)

    if search.success and search.data['pdb_ids']:
        # Use experimental structure
        pdb_id = search.data['pdb_ids'][0]
        result = await pdb.execute(pdb_id)
        print(f"Using experimental: {pdb_id}")
        print(f"Resolution: {result.data['resolution']} Å")
    else:
        # Fall back to AlphaFold
        alphafold = AlphaFoldAdapter()
        result = await alphafold.execute(uniprot_id)
        print(f"Using AlphaFold prediction")
        print(f"Mean pLDDT: {result.data['plddt_scores']['mean']}")

    return result

# Run analysis
asyncio.run(analyze_protein("P04637", "TP53 tumor suppressor"))
```

---

**Ready to use!** All adapters are production-ready and fully tested. 🚀
