# Tox21 Adapter - Quick Start Guide

## What is Tox21?

Tox21 is a collaboration between NIH, EPA, and FDA providing toxicity screening data for 12,000+ compounds across 50+ high-throughput assays. This adapter queries Tox21 data via PubChem's BioAssay API.

## Installation

No additional packages required - uses standard PharmForge dependencies.

## 5-Minute Quickstart

### 1. Basic Toxicity Query

```python
from adapters.tox21 import Tox21Adapter

adapter = Tox21Adapter()

# Query aspirin toxicity
smiles = "CC(=O)Oc1ccccc1C(=O)O"
result = await adapter.execute(smiles)

if result.success:
    data = result.data
    print(f"Active assays: {data['summary']['active_assays']}")
    print(f"Overall risk: {data['summary']['overall_risk']}")
```

### 2. Check Specific Toxicity Endpoints

```python
# Focus on endocrine disruption
input_data = {
    "query": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    "assays": ["NR-AR", "NR-ER", "NR-PPAR-gamma"],
    "include_inactive": False
}

result = await adapter.execute(input_data)

if result.success:
    for assay in result.data['toxicity_results']:
        print(f"{assay['assay_name']}: {assay['activity']}")
```

### 3. Batch Screening

```python
compounds = [
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
]

for smiles in compounds:
    result = await adapter.execute(smiles)
    if result.success:
        print(f"{smiles}: {result.data['summary']['overall_risk']} risk")
    await asyncio.sleep(0.3)  # Rate limiting
```

## Key Assay Categories

### Nuclear Receptors (Endocrine Disruption)
- `NR-AR`: Androgen Receptor
- `NR-ER`: Estrogen Receptor
- `NR-PPAR-gamma`: Metabolic disruption

### Stress Response (Cellular Toxicity)
- `SR-HSE`: Heat shock response
- `SR-MMP`: Mitochondrial toxicity
- `SR-p53`: DNA damage
- `SR-ARE`: Oxidative stress
- `SR-ATAD5`: Genotoxicity

## Risk Interpretation

- **Low**: <20% assays active - Safe for further development
- **Medium**: 20-50% active - Requires additional testing
- **High**: >50% active - Significant toxicity concerns

## Common Use Cases

1. **Early safety screening** - Filter compounds before synthesis
2. **Endocrine disruption check** - Regulatory requirement
3. **Comparative toxicity** - Compare analogs
4. **Structure-toxicity relationships** - Identify problematic scaffolds

## Output Fields

```python
{
    "compound": {
        "cid": 2244,  # PubChem ID
        "query": "CC(=O)Oc1ccccc1C(=O)O",
        "query_type": "smiles"
    },
    "toxicity_results": [
        {
            "assay_name": "Androgen Receptor",
            "activity": "active",  # or "inactive"
            "ac50": 15.3,  # Activity concentration (uM)
            "endpoint": "endocrine_disruption",
            "confidence": "high"  # Based on curve quality
        }
    ],
    "summary": {
        "active_assays": 2,
        "total_assays": 50,
        "toxicity_flags": ["endocrine_disruption"],
        "overall_risk": "low"
    },
    "warnings": [
        "Potential endocrine disruption activity detected"
    ]
}
```

## Query Options

### By SMILES
```python
result = await adapter.execute("CCO")
```

### By PubChem CID
```python
result = await adapter.execute(2244)  # Aspirin
```

### By Compound Name
```python
input_data = {
    "query_type": "name",
    "query": "Ibuprofen"
}
result = await adapter.execute(input_data)
```

### With Filters
```python
input_data = {
    "query": smiles,
    "assays": ["NR-AR", "NR-ER"],  # Specific assays only
    "include_inactive": False,      # Only show hits
    "ac50_threshold": 10.0          # Only potent hits (uM)
}
```

## Best Practices

1. **Rate limiting**: Wait 0.3s between queries
2. **Caching**: Results are automatically cached
3. **Interpretation**: Consider AC50 values, not just activity
4. **Validation**: Cross-check high-risk compounds with other tools
5. **Documentation**: Always cite Tox21 in publications

## Troubleshooting

**Compound not found?**
- Try searching PubChem first to get CID
- Not all compounds are in Tox21 database

**Slow queries?**
- Enable caching (default)
- Results cached permanently (static data)

**No active assays?**
- This is good! Low toxicity risk
- May indicate compound not tested in all assays

## Next Steps

- See `README.md` for complete documentation
- Run `example_usage.py` for more examples
- Check `backend/tests/test_tox21_adapter.py` for test patterns

## Citation

```
Huang, R., et al. (2018). The NCATS Tox21 Program: Enabling
technologies for high-throughput screening and data integration.
Frontiers in Public Health, 6, 190.
```

## Resources

- Tox21 Website: https://tripod.nih.gov/tox21/
- PubChem BioAssay: https://pubchem.ncbi.nlm.nih.gov/
- PharmForge Docs: See main README
