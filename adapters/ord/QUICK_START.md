# ORD Adapter Quick Start Guide

## Installation

```bash
pip install rdkit aiohttp
```

## 30-Second Usage

```python
from adapters.ord.adapter import ORDAdapter
import asyncio

adapter = ORDAdapter()

# Search for reactions producing aspirin
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

print(f"Found {result.data['reactions_found']} reactions")
print(f"Average yield: {result.data['statistics']['average_yield']}%")
```

## Common Patterns

### 1. Product Search (Default)

```python
# Search: "What reactions make this molecule?"
result = await adapter.execute("CCO")  # Ethanol
```

### 2. Reactant Search

```python
# Search: "What can I make with benzene?"
result = await adapter.execute({
    "smiles": "c1ccccc1",
    "search_type": "reactant"
})
```

### 3. High-Yield Reactions Only

```python
# Only show reactions with >70% yield
result = await adapter.execute("CCO", min_yield=70.0)
```

### 4. Get Reaction Conditions

```python
result = await adapter.execute("CC(=O)O")

for reaction in result.data['reactions'][:5]:
    conditions = reaction['conditions']
    print(f"Temp: {conditions['temperature']} °C")
    print(f"Solvent: {conditions['solvents']}")
    print(f"Reagents: {conditions['reagents']}")
```

### 5. Find Literature References

```python
result = await adapter.execute("c1ccccc1")

for reaction in result.data['reactions']:
    if reaction['literature']['doi']:
        print(f"DOI: {reaction['literature']['doi']}")
```

### 6. Batch Processing

```python
compounds = ["CCO", "CC(=O)C", "c1ccccc1"]

results = await asyncio.gather(*[
    adapter.execute(smiles) for smiles in compounds
])
```

## Configuration

```python
# Custom settings
adapter = ORDAdapter(config={
    "max_results": 100,      # Return more reactions
    "min_yield": 50.0,       # Filter low-yield reactions
    "timeout": 60,           # Longer timeout
    "rate_limit_delay": 0.3  # Faster queries
})
```

## Output Structure

```python
result.data = {
    "reactions_found": 42,
    "reactions": [
        {
            "reaction_id": "ord-12345",
            "yield": 85.5,
            "reactants": [...],
            "products": [...],
            "conditions": {
                "temperature": 80.0,
                "solvents": ["toluene"],
                "reagents": ["Pd(PPh3)4"]
            },
            "literature": {"doi": "..."}
        }
    ],
    "statistics": {
        "average_yield": 72.3,
        "max_yield": 95.0
    },
    "common_solvents": ["toluene", "THF"],
    "common_reagents": ["Pd(PPh3)4", "K2CO3"]
}
```

## Error Handling

```python
result = await adapter.execute("invalid_smiles")

if result.success:
    # Process results
    data = result.data
else:
    # Handle error
    print(f"Error: {result.error}")
```

## Caching

```python
# Automatic caching (default)
result1 = await adapter(smiles)  # API call
result2 = await adapter(smiles)  # Cache hit (instant)

# Disable cache
result = await adapter(smiles, use_cache=False)
```

## Integration with PharmForge

```python
# In pipeline configuration
{
    "adapter": "ord",
    "inputs": {"smiles": "$compound"},
    "params": {
        "search_type": "product",
        "min_yield": 60.0
    }
}
```

## Tips

1. **Canonicalize SMILES**: The adapter auto-canonicalizes for consistent caching
2. **Filter by yield**: Use `min_yield` to focus on practical reactions
3. **Batch processing**: Use `asyncio.gather()` for multiple compounds
4. **Check literature**: Most reactions have DOI references for validation
5. **Analyze conditions**: Extract common solvents/reagents for optimization

## Troubleshooting

| Problem | Solution |
|---------|----------|
| No results | Try lower `similarity_threshold` or `min_yield` |
| Timeout | Increase `config["timeout"]` |
| Rate limiting | Increase `rate_limit_delay` |
| Invalid SMILES | Check with RDKit: `Chem.MolFromSmiles(smiles)` |

## Examples

See `example_usage.py` for 7 detailed examples including:
- Product/reactant/reagent search
- Condition optimization
- Catalyst discovery
- Retrosynthesis validation
- Batch processing

## Resources

- **Documentation**: [README.md](README.md)
- **API**: https://open-reaction-database.org/
- **Paper**: Kearnes et al., JACS 2021
- **Support**: File issues on PharmForge GitHub

---

**Need help?** Check the full [README.md](README.md) or run `python example_usage.py`
