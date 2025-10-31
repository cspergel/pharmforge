# COCONUT Adapter

**COlleCtion of Open Natural prodUcTs** - Natural Products Database Adapter for PharmForge

## Overview

The COCONUT adapter provides access to the world's largest open database of natural products, containing over 400,000 natural product structures with biological source information, molecular properties, and bioactivity annotations.

### Key Features

- **400,000+ Natural Products**: Comprehensive coverage of natural product chemical space
- **Biological Source Information**: Organism, taxonomy, and biosource tracking
- **Multiple Search Modes**: By structure, name, properties, or organism
- **Molecular Descriptors**: Pre-calculated physicochemical properties
- **Open Access**: CC0 public domain license
- **No API Key Required**: Free public access

## Database Information

- **URL**: https://coconut.naturalproducts.net/
- **API**: https://coconut.naturalproducts.net/api/
- **License**: CC0 (public domain)
- **Citation**: See [COCONUT website](https://coconut.naturalproducts.net/)

## Installation

No additional dependencies required - uses `requests` and `aiohttp` which are already in PharmForge.

```bash
# No installation needed - adapter is ready to use
```

## Usage

### Basic Usage

```python
import asyncio
from adapters.coconut.adapter import COCONUTAdapter

async def search_natural_products():
    adapter = COCONUTAdapter()

    # Simple SMILES search
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    result = await adapter(smiles)

    if result.success:
        print(f"Found {result.data['num_results']} natural products")
        for compound in result.data['results'][:5]:
            print(f"  - {compound['name']}")
            print(f"    From: {compound['natural_source']['organism']}")
    else:
        print(f"Error: {result.error}")

asyncio.run(search_natural_products())
```

### Advanced Queries

#### Search with Property Filters

```python
query = {
    "query_type": "smiles",
    "query": "CC(C)=CCC1C(=O)CCC1C",  # Terpene scaffold
    "filters": {
        "molecular_weight": {"min": 200, "max": 500},
        "logp": {"min": 1, "max": 5},
        "organism": "Streptomyces"  # Optional organism filter
    },
    "include_taxonomy": True,
    "include_activities": True,
    "max_results": 50
}

result = await adapter(query)
```

#### Search by Organism

```python
query = {
    "query_type": "organism",
    "query": "Streptomyces",
    "max_results": 100
}

result = await adapter(query)
```

#### Search by Name

```python
query = {
    "query_type": "name",
    "query": "Taxol",
    "include_taxonomy": True,
    "include_activities": True
}

result = await adapter(query)
```

#### Get Compound by ID

```python
query = {
    "query_type": "id",
    "query": "CNP0123456"
}

result = await adapter(query)
```

## Input Format

### Simple String Input
```python
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin-like scaffold
result = await adapter(smiles)
```

### Structured Query Input
```python
{
    "query_type": "smiles",      # Required: "smiles", "name", "inchi", "inchikey", "properties", "organism", "id"
    "query": "your_query_here",  # Required: The search query
    "filters": {                 # Optional: Property filters
        "molecular_weight": {"min": 200, "max": 500},
        "logp": {"min": 0, "max": 5},
        "organism": "plant_name"  # Optional: filter by organism
    },
    "include_taxonomy": true,     # Optional: Include taxonomy info (default: true)
    "include_activities": true,   # Optional: Include biological activities (default: true)
    "max_results": 50            # Optional: Maximum results (default: 50)
}
```

## Output Format

```python
{
    "results": [
        {
            "coconut_id": "CNP0123456",
            "name": "Taxol",
            "smiles": "CC1=C2[C@@H](C(=O)...",
            "inchi": "InChI=1S/C47H51NO14/c1-25-31...",
            "inchikey": "RCINICONZNJXQF-UHFFFAOYSA-N",
            "molecular_formula": "C47H51NO14",
            "molecular_weight": 853.91,
            "natural_source": {
                "organism": "Taxus brevifolia",
                "common_name": "Pacific Yew",
                "taxonomy": "Plantae;Pinophyta;Pinopsida;Pinales;Taxaceae;Taxus"
            },
            "properties": {
                "logp": 4.2,
                "hba": 14,        # H-bond acceptors
                "hbd": 4,         # H-bond donors
                "rotatable_bonds": 12,
                "tpsa": 221.3     # Topological polar surface area
            },
            "biological_activities": [
                "Anticancer",
                "Microtubule stabilizer"
            ],
            "references": ["PMID:12345678"]
        }
    ],
    "num_results": 1,
    "query_type": "name",
    "query": "Taxol",
    "warnings": []
}
```

## Use Cases

### 1. Natural Product Library Screening

Build a virtual screening library from natural products:

```python
# Search for drug-like natural products
query = {
    "query_type": "properties",
    "filters": {
        "molecular_weight": {"min": 200, "max": 600},
        "logp": {"min": 0, "max": 5}
    },
    "max_results": 1000
}

result = await adapter(query)

# Filter for Lipinski's Rule of Five
drug_like = [
    c for c in result.data['results']
    if c['molecular_weight'] <= 500
    and c['properties']['logp'] <= 5
    and c['properties']['hbd'] <= 5
    and c['properties']['hba'] <= 10
]

print(f"Found {len(drug_like)} drug-like natural products")
```

### 2. Scaffold Mining from Nature

Find natural scaffolds similar to a known bioactive:

```python
# Find natural products similar to aspirin
aspirin = "CC(=O)Oc1ccccc1C(=O)O"
result = await adapter(aspirin)

# Analyze scaffold diversity
scaffolds = {}
for compound in result.data['results']:
    source = compound['natural_source']
    kingdom = source['taxonomy'].split(';')[0]
    scaffolds.setdefault(kingdom, []).append(compound)

for kingdom, compounds in scaffolds.items():
    print(f"{kingdom}: {len(compounds)} similar compounds")
```

### 3. Biosource Tracking

Find which organisms produce compounds with specific properties:

```python
query = {
    "query_type": "organism",
    "query": "Streptomyces",
    "filters": {
        "molecular_weight": {"min": 300, "max": 800}
    },
    "max_results": 200
}

result = await adapter(query)

# Analyze molecular weight distribution
weights = [c['molecular_weight'] for c in result.data['results']]
print(f"MW range: {min(weights):.1f} - {max(weights):.1f}")
print(f"Average MW: {sum(weights)/len(weights):.1f}")
```

### 4. Traditional Medicine Validation

Validate traditional medicine claims with natural product data:

```python
# Search for compounds from a medicinal plant
query = {
    "query_type": "organism",
    "query": "Artemisia annua",  # Sweet wormwood (antimalarial)
    "include_activities": True,
    "max_results": 50
}

result = await adapter(query)

# Check for antimalarial activities
antimalarials = [
    c for c in result.data['results']
    if 'biological_activities' in c
    and any('malaria' in act.lower() for act in c['biological_activities'])
]

print(f"Found {len(antimalarials)} compounds with antimalarial activity")
```

### 5. Natural Product-Inspired Drug Design

Generate analogs inspired by natural scaffolds:

```python
# Find natural products in a specific chemical space
result = await adapter("CC(C)=CCC1C(=O)CCC1C")  # Terpene scaffold

natural_scaffolds = []
for compound in result.data['results']:
    natural_scaffolds.append({
        'smiles': compound['smiles'],
        'source': compound['natural_source']['organism'],
        'activities': compound.get('biological_activities', [])
    })

# These scaffolds can now be used as starting points for:
# - Fragment-based drug design
# - De novo generation
# - Scaffold hopping
```

### 6. Taxonomic Analysis

Analyze which taxonomic groups produce specific compound classes:

```python
# Search for alkaloids (nitrogen-containing natural products)
# This would require filtering by SMILES patterns or molecular formula

results_by_kingdom = {}

for organism in ['Bacteria', 'Fungi', 'Plantae', 'Animalia']:
    query = {
        "query_type": "organism",
        "query": organism,
        "max_results": 100
    }
    result = await adapter(query)
    results_by_kingdom[organism] = result.data['num_results']

print("Natural products by kingdom:")
for kingdom, count in results_by_kingdom.items():
    print(f"  {kingdom}: {count}")
```

## Integration with PharmForge Pipeline

### Example Pipeline: Natural Product Screening

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()

# Step 1: Find natural products from a specific source
pipeline.add_step("coconut", {
    "query_type": "organism",
    "query": "Streptomyces",
    "filters": {
        "molecular_weight": {"min": 300, "max": 700}
    },
    "max_results": 100
})

# Step 2: Calculate ADMET properties
pipeline.add_step("tdc_admet", {
    "properties": ["Caco2_Wang", "Solubility_AqSolDB", "Lipophilicity_AstraZeneca"]
})

# Step 3: Predict targets
pipeline.add_step("swisstarget", {
    "organism": "Homo sapiens"
})

# Step 4: Dock to target
pipeline.add_step("vina_docking", {
    "receptor": "protein.pdb",
    "center": [10.0, 15.0, 20.0],
    "size": [20.0, 20.0, 20.0]
})

# Execute pipeline
results = await pipeline.execute()
```

## Configuration

The adapter can be configured through the `config` dictionary:

```python
adapter = COCONUTAdapter()
adapter.config.update({
    "rate_limit_delay": 1.0,  # Increase delay between requests
    "timeout": 120,            # Longer timeout for large searches
    "max_results": 100,        # Default max results
    "default_similarity_threshold": 0.85  # Stricter similarity
})
```

## API Endpoints

The adapter uses the following COCONUT API endpoints:

- `/search/structure` - Structure similarity search
- `/search/name` - Name/text search
- `/search/organism` - Organism/taxonomy search
- `/compound/{id}` - Get compound details by ID
- `/taxonomy` - Browse taxonomy tree

## Rate Limiting

The adapter implements respectful rate limiting:

- Default delay: 0.5 seconds between requests
- Configurable via `config["rate_limit_delay"]`
- Public API - no key required but be respectful

## Error Handling

The adapter handles common errors gracefully:

```python
result = await adapter(query)

if not result.success:
    if "timeout" in result.error.lower():
        print("Request timed out - try reducing max_results")
    elif "not found" in result.error.lower():
        print("No compounds matched your query")
    else:
        print(f"Error: {result.error}")
```

## Caching

Results are automatically cached using PharmForge's caching system:

```python
# First call - hits API
result1 = await adapter(smiles, use_cache=True)

# Second call - returns cached result
result2 = await adapter(smiles, use_cache=True)

print(f"Cache hit: {result2.cache_hit}")  # True
```

## Testing

Run the example usage script:

```bash
cd adapters/coconut
python example_usage.py
```

## Limitations

1. **API Availability**: Requires internet connection and COCONUT API availability
2. **Search Accuracy**: Structure similarity depends on COCONUT's fingerprinting method
3. **Rate Limits**: Public API may have undocumented rate limits
4. **Data Coverage**: Not all compounds have complete metadata
5. **Real-time Updates**: Database may not include the latest published natural products

## Complementary Adapters

The COCONUT adapter works well with:

- **ChEMBL**: Cross-reference bioactivity data
- **PubChem**: Validate structures and get additional properties
- **RDKit**: Calculate additional descriptors
- **TDC ADMET**: Predict ADMET properties
- **Vina/DiffDock**: Dock natural products to targets
- **AiZynthFinder**: Plan synthesis routes for natural products

## Best Practices

1. **Filter Early**: Use property filters to reduce result set size
2. **Cache Results**: Enable caching for repeated queries
3. **Batch Processing**: For large screens, process in batches
4. **Validate Structures**: Cross-check SMILES with RDKit
5. **Respect Rate Limits**: Don't hammer the public API

## Troubleshooting

### "No results found"
- Check SMILES validity with RDKit first
- Try relaxing property filters
- Try searching by name instead of structure

### "Timeout error"
- Reduce `max_results`
- Increase `timeout` in config
- Check internet connection

### "Invalid query type"
- Ensure query_type is one of: smiles, name, inchi, inchikey, properties, organism, id
- Check query dictionary structure

## Contributing

To improve this adapter:

1. Add support for additional search parameters
2. Implement batch search functionality
3. Add more sophisticated property filters
4. Improve error messages
5. Add unit tests

## References

- COCONUT Database: https://coconut.naturalproducts.net/
- Natural Products Atlas: https://www.npatlas.org/
- DNP (Dictionary of Natural Products)
- MarinLit (Marine Natural Products Database)

## License

This adapter is part of PharmForge and follows the PharmForge license.

The COCONUT database is licensed under CC0 (public domain).

## Version History

- **1.0.0** (2025-10-30): Initial implementation
  - Structure similarity search
  - Name search
  - Organism/taxonomy search
  - Property filtering
  - Biological activity annotations

## Contact

For issues or questions about this adapter, please open an issue on the PharmForge GitHub repository.

---

**Note**: This adapter provides access to natural product data for research purposes. Always validate computational predictions with experimental data before making any claims about biological activity or drug-likeness.
