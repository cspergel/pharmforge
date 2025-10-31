# Tox21 Adapter

Provides access to NIH Toxicology in the 21st Century (Tox21) screening data for safety profiling in drug discovery.

## Overview

Tox21 is a collaboration between NIH, EPA, and FDA to develop predictive models for toxicity screening using high-throughput assays. This adapter queries Tox21 data via PubChem's BioAssay database.

**Key Features:**
- 12,000+ compounds tested
- 50+ toxicity assays
- Nuclear receptor binding assays
- Stress response pathway assays
- Cytotoxicity and genotoxicity endpoints
- AC50 values (activity concentrations)
- Public domain data

## Installation

No additional installation required. Uses standard Python libraries:
- `aiohttp` (already in PharmForge requirements)
- `asyncio` (standard library)

## Tox21 Assay Categories

### Nuclear Receptors (NR)

Test for endocrine disruption and metabolic effects:

- **NR-AR**: Androgen Receptor (agonist/antagonist)
- **NR-ER**: Estrogen Receptor Alpha (agonist/antagonist)
- **NR-ER-LBD**: Estrogen Receptor Ligand Binding Domain
- **NR-AR-LBD**: Androgen Receptor Ligand Binding Domain
- **NR-AhR**: Aryl Hydrocarbon Receptor (metabolic activation)
- **NR-Aromatase**: CYP19A1 inhibition
- **NR-PPAR-gamma**: Peroxisome Proliferator-Activated Receptor Gamma

### Stress Response (SR)

Test for cellular stress and damage:

- **SR-ARE**: Antioxidant Response Element (oxidative stress)
- **SR-ATAD5**: DNA damage response (genotoxicity)
- **SR-HSE**: Heat Shock Response Element
- **SR-MMP**: Mitochondrial Membrane Potential (mitochondrial toxicity)
- **SR-p53**: p53 tumor suppressor activation (DNA damage)

## Usage Examples

### Basic Query by SMILES

```python
from adapters.tox21 import Tox21Adapter

adapter = Tox21Adapter()

# Query toxicity data for aspirin
smiles = "CC(=O)Oc1ccccc1C(=O)O"
result = await adapter.execute(smiles)

if result.success:
    data = result.data
    print(f"Compound CID: {data['compound']['cid']}")
    print(f"Total assays: {data['summary']['total_assays']}")
    print(f"Active assays: {data['summary']['active_assays']}")
    print(f"Overall risk: {data['summary']['overall_risk']}")

    # Show toxicity flags
    if data['summary']['toxicity_flags']:
        print(f"Toxicity flags: {', '.join(data['summary']['toxicity_flags'])}")

    # Show active assays
    for assay in data['toxicity_results']:
        if assay['activity'] == 'active':
            print(f"  {assay['assay_name']}: {assay['endpoint']}")
            if assay['ac50']:
                print(f"    AC50: {assay['ac50']:.2f} uM")
```

### Query with Dictionary Input

```python
# Query by compound name
input_data = {
    "query_type": "name",
    "query": "Bisphenol A",
    "assays": ["NR-ER", "NR-AR", "NR-PPAR-gamma"],  # Specific assays
    "include_inactive": False,  # Only active results
    "ac50_threshold": 50.0  # Filter by potency (uM)
}

result = await adapter.execute(input_data)
```

### Query by PubChem CID

```python
# Direct query by CID
cid = 2244  # Aspirin
result = await adapter.execute(cid)
```

### Filter by Assay Category

```python
# Query only nuclear receptor assays
nuclear_receptor_assays = [
    "NR-AR", "NR-ER", "NR-ER-LBD",
    "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-PPAR-gamma"
]

input_data = {
    "query": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    "assays": nuclear_receptor_assays
}

result = await adapter.execute(input_data)
```

### Batch Screening

```python
import asyncio

async def screen_compounds(smiles_list):
    adapter = Tox21Adapter()
    results = []

    for smiles in smiles_list:
        result = await adapter.execute(smiles)
        if result.success:
            results.append({
                'smiles': smiles,
                'cid': result.data['compound']['cid'],
                'active_assays': result.data['summary']['active_assays'],
                'overall_risk': result.data['summary']['overall_risk'],
                'toxicity_flags': result.data['summary']['toxicity_flags']
            })
        # Rate limiting
        await asyncio.sleep(0.3)

    return results

# Screen multiple compounds
compounds = [
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
]

results = await screen_compounds(compounds)
```

## Output Format

### Complete Result Structure

```python
{
    "compound": {
        "cid": 2244,
        "query": "CC(=O)Oc1ccccc1C(=O)O",
        "query_type": "smiles"
    },
    "toxicity_results": [
        {
            "assay_id": "NR-AR",
            "assay_name": "Androgen Receptor",
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor",
            "pubchem_aid": 743040,
            "activity": "active",
            "ac50": 15.3,  # uM
            "efficacy": 78.5,  # % of positive control
            "curve_class": 1.1,  # Quality metric
            "confidence": "high"
        },
        {
            "assay_id": "NR-ER",
            "assay_name": "Estrogen Receptor Alpha",
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor",
            "pubchem_aid": 743077,
            "activity": "inactive",
            "ac50": null,
            "efficacy": null,
            "curve_class": null,
            "confidence": "unknown"
        }
    ],
    "summary": {
        "total_assays": 50,
        "active_assays": 2,
        "inactive_assays": 48,
        "toxicity_flags": ["endocrine_disruption"],
        "overall_risk": "low"  # low/medium/high/unknown
    },
    "warnings": [
        "Potential endocrine disruption activity detected"
    ]
}
```

### Curve Class Interpretation

Curve classes indicate the quality/reliability of the dose-response curve:

- **1.1, 1.2**: High quality, complete curves (high confidence)
- **2.1, 2.2**: Incomplete curves (high confidence)
- **3, 4**: Ambiguous curves (medium confidence)
- **Other**: Low confidence

## Toxicity Endpoints

### Endocrine Disruption
- Androgen/estrogen receptor activation or inhibition
- Can affect hormone signaling and development

### Metabolic Activation
- Aryl hydrocarbon receptor activation
- Can lead to metabolism of procarcinogens

### Metabolic Disruption
- PPAR-gamma activation
- Can affect glucose/lipid metabolism

### Oxidative Stress
- Antioxidant response element activation
- Indicates cellular stress response

### Genotoxicity
- DNA damage response pathway activation
- Potential for genetic damage

### Mitochondrial Toxicity
- Mitochondrial membrane potential disruption
- Can lead to cell death

### DNA Damage
- p53 tumor suppressor activation
- Indicates DNA damage response

## Use Cases

### 1. Early Safety Profiling

Screen drug candidates early to identify potential toxicity issues:

```python
async def safety_screen(smiles):
    adapter = Tox21Adapter()
    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        risk = data['summary']['overall_risk']

        if risk == "high":
            print("WARNING: High toxicity risk detected")
            print(f"Flags: {data['summary']['toxicity_flags']}")
            return False
        elif risk == "medium":
            print("CAUTION: Medium toxicity risk")
            return True
        else:
            print("PASS: Low toxicity risk")
            return True
```

### 2. Endocrine Disruption Screening

Focus on endocrine-related assays:

```python
endocrine_assays = ["NR-AR", "NR-ER", "NR-ER-LBD", "NR-AR-LBD", "NR-Aromatase"]

input_data = {
    "query": smiles,
    "assays": endocrine_assays,
    "include_inactive": False
}

result = await adapter.execute(input_data)
```

### 3. Comparative Toxicity Analysis

Compare toxicity profiles of similar compounds:

```python
async def compare_toxicity(compound_dict):
    adapter = Tox21Adapter()
    comparison = {}

    for name, smiles in compound_dict.items():
        result = await adapter.execute(smiles)
        if result.success:
            comparison[name] = {
                'active_assays': result.data['summary']['active_assays'],
                'risk': result.data['summary']['overall_risk'],
                'flags': result.data['summary']['toxicity_flags']
            }
        await asyncio.sleep(0.3)

    return comparison

compounds = {
    "Compound A": "CC(=O)Oc1ccccc1C(=O)O",
    "Compound B": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
}

comparison = await compare_toxicity(compounds)
```

## Assay Information

Get details about specific assays:

```python
adapter = Tox21Adapter()

# Get info about specific assay
info = adapter.get_assay_info("NR-AR")
print(f"Assay: {info['name']}")
print(f"Endpoint: {info['endpoint']}")
print(f"PubChem AIDs: {info['aids']}")

# List all assays
all_assays = adapter.list_assays()
print(f"Total assays: {len(all_assays)}")

# List by category
nr_assays = adapter.list_assays(category="nuclear_receptor")
sr_assays = adapter.list_assays(category="stress_response")
```

## Interpretation Guidelines

### Overall Risk Levels

- **Low**: 0-20% of assays active - Compound shows minimal toxicity signals
- **Medium**: 20-50% of assays active - Some toxicity concerns, requires further evaluation
- **High**: >50% of assays active - Significant toxicity signals, likely not suitable for drug development
- **Unknown**: No data available

### AC50 Values

AC50 (Activity Concentration 50%) indicates potency:
- **< 1 uM**: Very potent - high concern
- **1-10 uM**: Potent - moderate concern
- **10-100 uM**: Moderate potency - low concern
- **> 100 uM**: Weak activity - minimal concern

### Key Toxicity Flags

1. **Endocrine Disruption**: Can interfere with hormone systems
2. **Genotoxicity**: May cause DNA damage
3. **Mitochondrial Toxicity**: Can lead to cell death
4. **Oxidative Stress**: Cellular stress response

## Data Limitations

- Not all compounds are in Tox21 (12,000+ tested)
- In vitro assays may not reflect in vivo toxicity
- Some assays have high false positive rates
- Data is static (not updated in real-time)
- Results should be confirmed with additional testing

## Caching

Results are automatically cached by PharmForge's caching system:
- Cache key based on query and parameters
- Static data (Tox21 data doesn't change)
- Improves performance for repeated queries

## Rate Limiting

- Default: 0.3 seconds between requests (3 req/sec)
- PubChem allows ~5 req/sec without API key
- Adjust via config: `adapter.config['rate_limit_delay'] = 0.5`

## Error Handling

```python
result = await adapter.execute(smiles)

if not result.success:
    print(f"Error: {result.error}")
    # Common errors:
    # - "Could not find PubChem CID for query"
    # - "Invalid input data for Tox21 query"
    # - API timeout or connection issues
```

## References

- **Tox21 Homepage**: https://tripod.nih.gov/tox21/
- **PubChem BioAssay**: https://pubchem.ncbi.nlm.nih.gov/
- **Tox21 Data Paper**: Huang et al. (2018) Nucleic Acids Res. doi:10.1093/nar/gky1033
- **Tox21 Overview**: https://ntp.niehs.nih.gov/whatwestudy/tox21/index.html

## Citation

If using Tox21 data in publications:

```
Huang, R., et al. (2018). The NCATS Tox21 Program: Enabling technologies
for high-throughput screening and data integration. Frontiers in Public
Health, 6, 190. doi:10.3389/fpubh.2018.00190
```

## License

Tox21 data is public domain (US government data). No restrictions on use.

## Support

For issues or questions:
- Check PharmForge documentation
- Review PubChem API documentation
- Open an issue on GitHub
