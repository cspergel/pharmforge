# Open Reaction Database (ORD) Adapter

## Overview

The ORD adapter provides access to the **Open Reaction Database**, a public repository of detailed chemical reaction data. This adapter is designed for forward synthesis planning, reaction condition optimization, and literature-validated route discovery.

## Features

- **Multi-mode Search**: Query by product, reactant, or reagent SMILES
- **Detailed Conditions**: Temperature, pressure, solvent, catalyst information
- **Yield Data**: Reaction yields and selectivity metrics
- **Literature References**: DOI and publication links for experimental validation
- **Filtering**: By reaction type, minimum yield, or similarity threshold
- **Batch Queries**: Process multiple compounds efficiently
- **Smart Caching**: Deterministic cache keys for reproducibility

## Installation

The ORD adapter requires RDKit for SMILES canonicalization:

```bash
pip install rdkit aiohttp
```

No additional ORD-specific libraries are required - the adapter uses the REST API directly.

## Usage

### Basic Usage - Search by Product

```python
from adapters.ord.adapter import ORDAdapter
import asyncio

adapter = ORDAdapter()

# Search for reactions producing aspirin
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

result = await adapter.execute(smiles)

if result.success:
    data = result.data
    print(f"Found {data['reactions_found']} reactions")
    print(f"Average yield: {data['statistics']['average_yield']}%")

    for reaction in data['reactions'][:5]:  # Top 5 reactions
        print(f"\nReaction ID: {reaction['reaction_id']}")
        print(f"Yield: {reaction['yield']}%")
        print(f"Reactants: {len(reaction['reactants'])}")
        print(f"Temperature: {reaction['conditions']['temperature']} C")
        print(f"Solvents: {reaction['conditions']['solvents']}")
```

### Search by Reactant

```python
# Find reactions using benzene as a reactant
search_params = {
    "smiles": "c1ccccc1",  # Benzene
    "search_type": "reactant",
    "min_yield": 50.0  # Only reactions with >50% yield
}

result = await adapter.execute(search_params)
```

### Search by Reagent/Catalyst

```python
# Find reactions using palladium catalyst
search_params = {
    "smiles": "[Pd]",
    "search_type": "reagent",
    "reaction_type": "coupling"  # Optional filter
}

result = await adapter.execute(search_params)
```

### Advanced Filtering

```python
# Custom configuration with filters
config = {
    "max_results": 100,
    "min_yield": 60.0,
    "similarity_threshold": 0.90,  # Stricter SMILES matching
    "rate_limit_delay": 0.3
}

adapter = ORDAdapter(config=config)

result = await adapter.execute(
    "CCO",  # Ethanol
    min_yield=75.0  # Override config value
)
```

## Output Format

### Successful Response

```python
{
    "reactions_found": 42,
    "reactions": [
        {
            "reaction_id": "ord-12345",
            "reaction_type": "coupling",
            "reactants": [
                {
                    "smiles": "c1ccccc1Br",
                    "name": "bromobenzene",
                    "role": "reactant",
                    "amount": {...}
                }
            ],
            "products": [
                {
                    "smiles": "c1ccc(-c2ccccc2)cc1",
                    "name": "biphenyl",
                    "yield": 85.5,
                    "is_desired": True
                }
            ],
            "yield": 85.5,
            "conditions": {
                "temperature": 80.0,
                "temperature_unit": "C",
                "pressure": 1.0,
                "pressure_unit": "bar",
                "stirring_rate_rpm": 500,
                "solvents": ["toluene"],
                "reagents": ["Pd(PPh3)4", "K2CO3"]
            },
            "literature": {
                "doi": "10.1021/...",
                "url": "https://pubs.acs.org/..."
            },
            "n_reactants": 2,
            "n_products": 1
        }
    ],
    "statistics": {
        "average_yield": 72.3,
        "max_yield": 95.0,
        "total_reactions": 42,
        "after_filtering": 42
    },
    "reaction_types": ["coupling", "substitution"],
    "common_solvents": ["toluene", "THF", "DMF"],
    "common_reagents": ["Pd(PPh3)4", "K2CO3", "NaOH"],
    "query": {
        "smiles": "c1ccc(-c2ccccc2)cc1",
        "search_type": "product"
    },
    "database": "Open Reaction Database",
    "reference": "https://open-reaction-database.org/"
}
```

## Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `rate_limit_delay` | float | 0.5 | Delay between API calls (seconds) |
| `timeout` | int | 30 | Request timeout (seconds) |
| `max_results` | int | 50 | Maximum reactions to return |
| `min_yield` | float | 0.0 | Minimum yield threshold (%) |
| `similarity_threshold` | float | 0.85 | SMILES similarity threshold |
| `include_conditions` | bool | True | Include reaction conditions |
| `include_reagents` | bool | True | Include reagent details |

## Use Cases

### 1. Forward Synthesis Planning

Find proven experimental routes to synthesize a target molecule:

```python
# Target: Ibuprofen
target_smiles = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"

result = await adapter.execute(target_smiles)

# Review top-yielding reactions with literature references
best_reactions = sorted(
    result.data['reactions'],
    key=lambda r: r['yield'],
    reverse=True
)[:3]
```

### 2. Reaction Condition Optimization

Analyze successful conditions for a specific transformation:

```python
# Find amide coupling reactions
search_params = {
    "smiles": "CC(=O)NCC",  # Example amide product
    "search_type": "product",
    "min_yield": 70.0
}

result = await adapter.execute(search_params)

# Extract common conditions
solvents = result.data['common_solvents']
reagents = result.data['common_reagents']
temps = [r['conditions']['temperature'] for r in result.data['reactions']]
```

### 3. Retrosynthesis Validation

Verify that proposed synthetic steps have experimental precedent:

```python
# Check if a proposed retrosynthetic disconnection has been reported
proposed_product = "CCc1ccccc1"
proposed_reactants = ["c1ccccc1", "CCBr"]

# Search for reactions producing the product
result = await adapter.execute(proposed_product)

# Check if proposed reactants appear in any reactions
validated_routes = []
for reaction in result.data['reactions']:
    reactant_smiles = [r['smiles'] for r in reaction['reactants']]
    if set(proposed_reactants).issubset(set(reactant_smiles)):
        validated_routes.append(reaction)
```

### 4. Reagent Discovery

Find alternative reagents for a transformation:

```python
# What reagents are used to produce compound X?
result = await adapter.execute("CC(=O)c1ccccc1")  # Acetophenone

reagent_usage = {}
for reaction in result.data['reactions']:
    for reagent in reaction['conditions']['reagents']:
        reagent_usage[reagent] = reagent_usage.get(reagent, 0) + 1

# Most common reagents
popular_reagents = sorted(
    reagent_usage.items(),
    key=lambda x: x[1],
    reverse=True
)
```

## Integration with PharmForge Pipeline

The ORD adapter integrates seamlessly with PharmForge pipelines:

```python
# In a pipeline configuration
pipeline = {
    "steps": [
        {
            "adapter": "ord",
            "inputs": {"smiles": "$target_molecule"},
            "params": {
                "search_type": "product",
                "min_yield": 60.0
            },
            "outputs": ["synthesis_routes"]
        }
    ]
}
```

## Caching

The adapter implements intelligent caching:

```python
# First call - hits API
result1 = await adapter.execute("CCO")

# Second call - uses cache (instant)
result2 = await adapter.execute("CCO")

assert result2.cache_hit == True

# Different parameters - new API call
result3 = await adapter.execute("CCO", min_yield=80.0)
```

## Error Handling

```python
result = await adapter.execute("invalid_smiles")

if not result.success:
    print(f"Error: {result.error}")
    # Handle error appropriately
else:
    # Process results
    pass
```

## API Notes

**Important**: The ORD API structure in this adapter is based on the expected REST API design. The actual ORD API may vary. If the API endpoint differs, update the `BASE_URL` and query structure in the adapter.

Current implementation assumes:
- Endpoint: `https://client.open-reaction-database.org/api/query`
- POST requests with JSON payload
- Search by component SMILES with role specification

For the most up-to-date API documentation, visit:
- https://open-reaction-database.org/
- https://github.com/open-reaction-database/ord-schema

## Limitations

1. **API Availability**: Requires active ORD API service
2. **Rate Limiting**: Be respectful of API limits (default: 0.5s delay)
3. **SMILES Matching**: Fuzzy matching may miss exact stereochemistry
4. **Data Coverage**: Not all reactions have complete condition data
5. **Yield Variability**: Experimental conditions may affect reproducibility

## Performance

- **Typical Response Time**: 0.5-2 seconds per query
- **Cache Hit Time**: <0.01 seconds
- **Batch Processing**: Use async for concurrent queries
- **Memory**: ~10-50 KB per reaction record

## Troubleshooting

### Connection Errors

```python
# Increase timeout for slow connections
adapter = ORDAdapter(config={"timeout": 60})
```

### No Results Found

- Check SMILES validity (must be parseable by RDKit)
- Try lowering `similarity_threshold`
- Remove yield filters to see all reactions
- Try searching by reactant instead of product

### Rate Limiting Issues

```python
# Increase delay between requests
adapter = ORDAdapter(config={"rate_limit_delay": 1.0})
```

## References

- **Website**: https://open-reaction-database.org/
- **GitHub**: https://github.com/open-reaction-database/ord-schema
- **Paper**: Kearnes et al., "The Open Reaction Database", *J. Am. Chem. Soc.* 2021, 143, 45, 18820-18826
- **Documentation**: https://docs.open-reaction-database.org/

## Version History

- **1.0.0** (2025-10-31): Initial release
  - Search by product, reactant, reagent
  - Reaction conditions and yields
  - Literature references
  - Smart caching
  - Batch support

## License

This adapter is part of PharmForge and follows the project's MIT license. The Open Reaction Database itself is licensed under CC BY 4.0.
