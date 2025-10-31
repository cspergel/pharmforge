# CompTox Chemistry Dashboard Adapter

Access EPA chemistry and toxicity data through the CompTox Chemistry Dashboard API.

## Overview

The CompTox Chemistry Dashboard is a comprehensive database from the U.S. Environmental Protection Agency (EPA) containing over 900,000 chemicals with predicted toxicity endpoints, exposure pathways, physicochemical properties, bioactivity assays, and hazard classifications.

**Database:** https://comptox.epa.gov/dashboard/
**API Docs:** https://api-ccte.epa.gov/docs/
**License:** Public domain (U.S. EPA)
**No API Key Required:** Basic access is free

## Features

- **Chemical Search**: Search by DTXSID, CAS number, name, SMILES, or InChIKey
- **QSAR Predictions**: OPERA model toxicity predictions
- **Physicochemical Properties**: LogP, solubility, vapor pressure, melting/boiling points
- **Toxicity Endpoints**: Oral LD50, LC50, bioconcentration, biodegradation, mutagenicity
- **Bioactivity Data**: ToxCast/Tox21 assay results and gene targets
- **Exposure Analysis**: Consumer products, use categories, exposure pathways
- **Hazard Classifications**: GHS hazard codes and environmental fate
- **Environmental Fate**: Biodegradation, persistence, ecological hazard

## Installation

No additional dependencies required beyond the core PharmForge stack (aiohttp, asyncio).

```bash
pip install aiohttp
```

## Usage

### Basic Usage

```python
from adapters.comptox import CompToxAdapter

# Initialize adapter
adapter = CompToxAdapter()

# Search by chemical name
result = await adapter.execute("aspirin")

# Search by CAS number
result = await adapter.execute({
    "query": "50-78-2",
    "query_type": "cas"
})

# Search by DTXSID
result = await adapter.execute({
    "query": "DTXSID7020182",
    "query_type": "dtxsid"
})

# Search by SMILES
result = await adapter.execute({
    "query": "CC(=O)Oc1ccccc1C(=O)O",
    "query_type": "smiles"
})
```

### Advanced Usage

```python
# Selective data retrieval
result = await adapter.execute(
    "caffeine",
    include_toxicity=True,
    include_bioactivity=True,
    include_properties=True,
    include_exposure=True,
    include_hazards=True
)

# Access specific data
if result.success:
    data = result.data

    # Chemical identifiers
    chemical = data["chemical"]
    print(f"DTXSID: {chemical['dtxsid']}")
    print(f"CAS: {chemical['casrn']}")
    print(f"SMILES: {chemical['smiles']}")

    # Physicochemical properties
    props = data["properties"]
    print(f"LogP: {props['logp']}")
    print(f"Molecular Weight: {props['molecular_weight']}")

    # Toxicity predictions (OPERA models)
    toxicity = data["toxicity"]
    for endpoint in toxicity["predicted_endpoints"]:
        print(f"{endpoint['endpoint']}: {endpoint['value']} ({endpoint['source']})")

    # QSAR predictions
    qsar = toxicity["qsar_predictions"]
    print(f"Mutagenicity: {qsar.get('mutagenicity')}")
    print(f"Developmental Toxicity: {qsar.get('developmental_toxicity')}")

    # Bioactivity
    bioactivity = data["bioactivity"]
    print(f"Total Assays: {bioactivity['num_assays']}")
    print(f"Active Assays: {bioactivity['active_assays']}")
    print(f"Gene Targets: {', '.join(bioactivity['gene_targets'][:5])}")

    # Exposure pathways
    exposure = data["exposure"]
    print(f"Consumer Products: {', '.join(exposure['consumer_products'])}")
    print(f"Use Categories: {', '.join(exposure['use_categories'])}")

    # Hazard information
    hazard = data["hazard"]
    print(f"GHS Codes: {', '.join(hazard['ghs_classification'])}")
    print(f"Biodegradation: {hazard['environmental_fate'].get('biodegradation')}")
```

### With Caching

```python
# First call - fetches from API
result1 = await adapter(
    "benzene",
    use_cache=True,
    include_toxicity=True
)
print(f"Cache hit: {result1.cache_hit}")  # False

# Second call - retrieves from cache
result2 = await adapter(
    "benzene",
    use_cache=True,
    include_toxicity=True
)
print(f"Cache hit: {result2.cache_hit}")  # True
```

## Input Formats

### Simple Query (String)
```python
result = await adapter.execute("aspirin")
```
Defaults to name search.

### Structured Query (Dict)
```python
result = await adapter.execute({
    "query": "50-78-2",
    "query_type": "cas",  # or "dtxsid", "name", "smiles", "inchikey"
    "include_toxicity": True,
    "include_bioactivity": True,
    "include_exposure": True,
    "include_hazards": True
})
```

## Output Format

```python
{
    "chemical": {
        "dtxsid": "DTXSID7020182",
        "preferred_name": "Aspirin",
        "casrn": "50-78-2",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "molecular_formula": "C9H8O4",
        "molecular_weight": 180.157
    },
    "properties": {
        "logp": 1.19,
        "water_solubility": 4.6,  # g/L
        "vapor_pressure": 1.6e-7,  # mmHg
        "melting_point": 135,  # °C
        "boiling_point": 140  # °C
    },
    "toxicity": {
        "predicted_endpoints": [
            {
                "endpoint": "oral_rat_ld50",
                "value": 1750,
                "unit": "mg/kg",
                "source": "OPERA model"
            },
            {
                "endpoint": "bioconcentration_factor",
                "value": 3.16,
                "source": "OPERA model"
            }
        ],
        "qsar_predictions": {
            "mutagenicity": "negative",
            "carcinogenicity": "negative",
            "developmental_toxicity": "positive"
        }
    },
    "bioactivity": {
        "num_assays": 150,
        "active_assays": 45,
        "gene_targets": ["PTGS1", "PTGS2", "ALOX5"]
    },
    "exposure": {
        "consumer_products": ["pharmaceuticals", "cosmetics"],
        "exposure_pathways": ["oral", "dermal"],
        "use_categories": ["pharmaceutical", "food_additive"]
    },
    "hazard": {
        "ghs_classification": ["H302", "H319"],  # Harmful if swallowed, eye irritant
        "ecological_hazard": "low",
        "environmental_fate": {
            "biodegradation": "readily biodegradable",
            "persistence": "low"
        }
    },
    "warnings": []
}
```

## OPERA Models

CompTox provides QSAR predictions using OPERA (OPEn structure-activity Relationship App) models:

### Physicochemical Properties
- **LogP**: Octanol-water partition coefficient
- **MP**: Melting point (°C)
- **BP**: Boiling point (°C)
- **LogOH**: Atmospheric hydroxylation rate constant
- **KOA**: Octanol-air partition coefficient

### Environmental Fate
- **BCF**: Bioconcentration factor (log units)
- **BioDeg**: Ready biodegradability (probability)
- **HL**: Half-life in water/soil/sediment

### Toxicity Endpoints
- **Ames**: Mutagenicity (Ames test prediction)
- **RT**: Developmental toxicity
- **ER/AR**: Endocrine disruption (estrogen/androgen receptor)
- **Fish LC50**: Fish acute toxicity
- **Daphnia LC50**: Daphnia acute toxicity

## ToxCast/Tox21 Data

CompTox includes bioactivity data from EPA's high-throughput screening programs:

- **ToxCast**: ~1,000 assays across multiple biological targets
- **Tox21**: ~50 high-priority assays for toxicity pathways
- **Gene Targets**: Specific proteins, receptors, enzymes affected
- **Activity Flags**: Active/inactive/inconclusive

## GHS Hazard Classifications

Global Harmonized System (GHS) codes included:

- **H300-H3xx**: Health hazards (acute toxicity, corrosion, etc.)
- **H400-H4xx**: Environmental hazards (aquatic toxicity, etc.)
- **H200-H2xx**: Physical hazards (flammability, explosiveness)

## API Rate Limits

- **Rate Limit**: 10 requests/second (configured with 0.1s delay)
- **Timeout**: 45 seconds per request
- **Retries**: 3 attempts with exponential backoff
- **No Authentication**: Public API, no key required

## Use Cases

### 1. Comprehensive Toxicity Profiling
```python
result = await adapter.execute(
    "benzene",
    include_toxicity=True,
    include_hazards=True
)

toxicity = result.data["toxicity"]
print("Toxicity Endpoints:")
for endpoint in toxicity["predicted_endpoints"]:
    print(f"  {endpoint['endpoint']}: {endpoint['value']}")
```

### 2. Environmental Risk Assessment
```python
result = await adapter.execute("pesticide_name")

hazard = result.data["hazard"]
print(f"Biodegradation: {hazard['environmental_fate']['biodegradation']}")
print(f"Persistence: {hazard['environmental_fate']['persistence']}")
print(f"GHS Codes: {', '.join(hazard['ghs_classification'])}")
```

### 3. Green Chemistry Screening
```python
# Check if compound is environmentally friendly
result = await adapter.execute(compound_smiles)

bioactivity = result.data["bioactivity"]
hazard = result.data["hazard"]

is_green = (
    bioactivity["active_assays"] < 10 and
    hazard["environmental_fate"].get("biodegradation") == "readily biodegradable" and
    len(hazard["ghs_classification"]) == 0
)

print(f"Green chemistry candidate: {is_green}")
```

### 4. Consumer Product Safety
```python
result = await adapter.execute("ingredient_name")

exposure = result.data["exposure"]
print(f"Found in products: {', '.join(exposure['consumer_products'])}")
print(f"Exposure pathways: {', '.join(exposure['exposure_pathways'])}")

toxicity = result.data["toxicity"]
qsar = toxicity["qsar_predictions"]
print(f"Mutagenicity: {qsar.get('mutagenicity')}")
```

### 5. Drug Safety Screening
```python
# Screen drug candidate for toxicity flags
result = await adapter.execute(drug_smiles)

toxicity = result.data["toxicity"]
bioactivity = result.data["bioactivity"]

# Check for red flags
red_flags = []
if toxicity["qsar_predictions"].get("mutagenicity") == "positive":
    red_flags.append("Mutagenicity predicted")
if toxicity["qsar_predictions"].get("developmental_toxicity") == "positive":
    red_flags.append("Developmental toxicity predicted")
if bioactivity["active_assays"] > 100:
    red_flags.append("High promiscuity (many assay hits)")

print(f"Safety flags: {', '.join(red_flags) if red_flags else 'None'}")
```

## Error Handling

```python
result = await adapter.execute("nonexistent_chemical_xyz123")

if not result.success:
    print(f"Error: {result.error}")
    print(f"Metadata: {result.metadata}")
else:
    # Check warnings
    if result.data["warnings"]:
        print(f"Warnings: {', '.join(result.data['warnings'])}")
```

## Comparison with Other Toxicity Databases

| Feature | CompTox | TDC ADMET | pkCSM |
|---------|---------|-----------|-------|
| Chemical Coverage | 900,000+ | 100,000+ | Unlimited (predictions) |
| Toxicity Endpoints | 50+ (OPERA) | 22 | 30+ |
| Bioactivity Assays | Yes (ToxCast) | Limited | No |
| Environmental Fate | Yes | No | Limited |
| GHS Classifications | Yes | No | No |
| API Access | Free | Free | Free (rate limited) |
| QSAR Models | OPERA | Multiple | pkCSM |

## Integration Example

```python
from adapters.comptox import CompToxAdapter
from adapters.pubchem import PubChemAdapter

# Get compound from PubChem
pubchem = PubChemAdapter()
pc_result = await pubchem.execute("caffeine")

# Get toxicity from CompTox
comptox = CompToxAdapter()
ct_result = await comptox.execute({
    "query": pc_result.data["canonical_smiles"],
    "query_type": "smiles",
    "include_toxicity": True,
    "include_bioactivity": True
})

# Combine data
combined = {
    "smiles": pc_result.data["canonical_smiles"],
    "properties": pc_result.data,
    "toxicity": ct_result.data["toxicity"],
    "bioactivity": ct_result.data["bioactivity"]
}
```

## References

- **CompTox Dashboard**: https://comptox.epa.gov/dashboard/
- **API Documentation**: https://api-ccte.epa.gov/docs/
- **OPERA Models**: https://github.com/kmansouri/OPERA
- **ToxCast**: https://www.epa.gov/chemical-research/toxicity-forecasting
- **Tox21**: https://tripod.nih.gov/tox21/

## Support

For issues with this adapter, open an issue on the PharmForge GitHub repository.

For questions about the EPA CompTox Dashboard data or API, contact the EPA:
- Email: comptox@epa.gov
- Documentation: https://comptox.epa.gov/dashboard/help

## License

This adapter is part of PharmForge and is licensed under the MIT License.

The CompTox Chemistry Dashboard data is in the public domain (U.S. EPA).
