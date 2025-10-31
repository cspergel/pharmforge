# BRENDA Adapter - Quick Reference Card

## Import
```python
from adapters.brenda import BRENDAAdapter
import asyncio
```

## Basic Usage
```python
adapter = BRENDAAdapter()
result = await adapter.execute("1.1.1.1")  # EC number

if result.success:
    enzyme = result.data['enzyme']
    kinetics = result.data['kinetics']
    inhibitors = result.data['inhibitors']
```

## Query Formats

### Simple (String)
```python
result = await adapter.execute("1.1.1.1")
```

### Advanced (Dict)
```python
query = {
    "query_type": "ec_number",       # Required
    "query": "1.1.1.1",              # Required
    "parameters": ["km", "kcat"],    # Optional (default: ["km", "kcat", "ki"])
    "organisms": ["Homo sapiens"],   # Optional (default: all)
    "max_results": 50                # Optional (default: 50)
}
result = await adapter.execute(query)
```

## Output Structure
```python
{
    "enzyme": {
        "ec_number": "1.1.1.1",
        "systematic_name": "Alcohol:NAD+ oxidoreductase",
        "common_name": "Alcohol dehydrogenase",
        "reaction": "Alcohol + NAD+ = Aldehyde + NADH",
        "substrates": ["Ethanol", "Methanol"],
        "products": ["Acetaldehyde"],
        "cofactors": ["NAD+", "Zn2+"]
    },
    "kinetics": [
        {
            "parameter": "km",
            "value": 0.5,
            "unit": "mM",
            "substrate": "Ethanol",
            "organism": "Homo sapiens",
            "conditions": {"pH": 7.4, "temperature": 37},
            "reference": "PMID:12345678"
        }
    ],
    "inhibitors": [
        {
            "compound": "Disulfiram",
            "ki_value": 10.0,
            "unit": "nM",
            "inhibition_type": "competitive"
        }
    ]
}
```

## Common Patterns

### Filter by Organism
```python
result = await adapter.execute({
    "query": "1.1.1.1",
    "organisms": ["Homo sapiens"]
})
```

### Get Specific Parameters
```python
result = await adapter.execute({
    "query": "1.1.1.1",
    "parameters": ["km", "ki"]  # Only Km and Ki
})
```

### Find Potent Inhibitors
```python
result = await adapter.execute("1.1.1.1")
inhibitors = result.data['inhibitors']
potent = [i for i in inhibitors if i['ki_value'] < 100]  # < 100 nM
```

### Calculate Catalytic Efficiency
```python
kinetics = result.data['kinetics']
km_values = [k for k in kinetics if k['parameter'] == 'km']
kcat_values = [k for k in kinetics if k['parameter'] == 'kcat']

for km in km_values:
    kcat = next((k for k in kcat_values if k['substrate'] == km['substrate']), None)
    if kcat:
        efficiency = kcat['value'] / km['value']
        print(f"{km['substrate']}: {efficiency:.2f} mM⁻¹s⁻¹")
```

## Error Handling
```python
result = await adapter.execute("1.1.1.1")

if not result.success:
    print(f"Error: {result.error}")
else:
    # Process data
    pass
```

## Configuration
```python
adapter = BRENDAAdapter()

# Change rate limit (seconds between requests)
adapter.config['rate_limit_delay'] = 0.5

# Change timeout (seconds)
adapter.config['timeout'] = 120

# Change max results
adapter.config['max_results'] = 100
```

## Caching
```python
# With cache (default)
result = await adapter("1.1.1.1", use_cache=True)

# Without cache
result = await adapter("1.1.1.1", use_cache=False)

# Check if cached
if result.cache_hit:
    print("Data from cache")
```

## Registry Access
```python
from backend.core.adapters.protocol import registry

brenda = registry.get('brenda')
result = await brenda.execute("1.1.1.1")
```

## Common EC Numbers
```python
"1.1.1.1"      # Alcohol dehydrogenase
"3.4.21.5"     # Thrombin
"1.14.13.39"   # CYP2D6
"2.7.1.1"      # Hexokinase
"5.3.1.9"      # Glucose-6-phosphate isomerase
```

## Parameters
- `km` - Michaelis constant (substrate affinity)
- `kcat` - Turnover number (catalytic rate)
- `ki` - Inhibition constant (inhibitor potency)
- `kd` - Dissociation constant (binding affinity)

## Units
- Km: mM (millimolar)
- Kcat: s⁻¹ (per second)
- Ki: nM (nanomolar) or µM (micromolar)

## Tips
1. Use organism filtering for species-specific data
2. Cache is enabled by default - use it!
3. Lower Ki = more potent inhibitor
4. Competitive inhibitors are best for drug design
5. Kcat/Km = catalytic efficiency

## Files
- **Adapter:** `adapters/brenda/adapter.py`
- **Tests:** `backend/tests/test_brenda_adapter.py`
- **Docs:** `adapters/brenda/README.md`
- **Examples:** `adapters/brenda/example_usage.py`

## Support
- GitHub Issues: [pharmforge/issues]
- BRENDA Support: brenda@tu-braunschweig.de
- API Docs: https://www.brenda-enzymes.org/soap.php
