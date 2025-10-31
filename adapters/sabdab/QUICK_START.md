# SAbDab Adapter - Quick Start

Fast reference for using the SAbDab (Structural Antibody Database) adapter.

## Installation

```bash
# No additional dependencies needed
# SAbDab uses standard HTTP requests
```

## 5-Minute Quick Start

### 1. Basic Import

```python
from adapters.sabdab import SAbDabAdapter

adapter = SAbDabAdapter()
```

### 2. Look Up by PDB ID

```python
# Get antibody structure information
result = await adapter("7BWJ")

if result.success:
    ab = result.data['antibody_data']
    print(f"Antigen: {ab['antigen']}")
    print(f"Resolution: {ab['resolution']} Å")
```

### 3. Search by Antigen

```python
# Find antibodies against specific antigen
result = await adapter({
    "antigen": "spike",
    "species": "human"
})

antibodies = result.data['antibodies']
```

### 4. Filter by Quality

```python
# High-quality structures only
result = await adapter({
    "antigen": "Her2",
    "resolution_max": 2.5,
    "rfactor_max": 0.25
})
```

## Common Use Cases

### COVID-19 Antibodies

```python
result = await adapter({
    "antigen": "SARS-CoV-2",
    "species": "human",
    "resolution_max": 3.0
}, max_results=50)
```

### Nanobodies (VHH)

```python
result = await adapter({
    "ab_type": "VHH",
    "resolution_max": 2.5
})
```

### Therapeutic Antibodies (IgG)

```python
result = await adapter({
    "ab_type": "IgG",
    "antigen": "Her2",
    "species": "human"
})
```

### Fragment Structures (Fab)

```python
result = await adapter({
    "ab_type": "Fab",
    "resolution_max": 2.0
})
```

## Input Formats

### PDB ID (String)
```python
pdb_id = "7BWJ"
result = await adapter(pdb_id)
```

### Search Parameters (Dictionary)
```python
params = {
    "antigen": "spike",           # Antigen name
    "species": "human",           # Species filter
    "ab_type": "IgG",            # Antibody type
    "resolution_max": 3.0,        # Max resolution (Å)
    "rfactor_max": 0.30          # Max R-factor
}
result = await adapter(params)
```

## Output Format

### PDB Lookup
```python
{
    "pdb_id": "7BWJ",
    "antibody_data": {
        "resolution": 2.45,
        "antibody_type": "Fab",
        "antigen": "SARS-CoV-2 spike",
        "species": "Homo sapiens",
        "heavy_chain": "H",
        "light_chain": "L",
        "cdrs": {...}
    }
}
```

### Search Results
```python
{
    "total_results": 42,
    "antibodies": [
        {
            "pdb_id": "7BWJ",
            "resolution": 2.45,
            "antibody_type": "Fab",
            "antigen": "SARS-CoV-2 spike"
        },
        ...
    ]
}
```

## Antibody Types

| Type | Description |
|------|-------------|
| `IgG` | Full immunoglobulin G |
| `Fab` | Fragment antigen-binding |
| `scFv` | Single-chain variable fragment |
| `VHH` | Nanobody/single-domain |

## Quality Guidelines

### Resolution (Å)
- ⭐⭐⭐ `< 2.0` - Excellent
- ⭐⭐ `2.0-3.0` - Good
- ⭐ `> 3.0` - Acceptable

### R-factor
- ⭐⭐⭐ `< 0.25` - Excellent
- ⭐⭐ `0.25-0.30` - Acceptable
- ⭐ `> 0.30` - Poor

## Advanced Options

### Include CDR Sequences
```python
result = await adapter("7BWJ", include_cdrs=True)
cdrs = result.data['antibody_data']['cdrs']
```

### Limit Results
```python
result = await adapter(params, max_results=20)
```

### Disable Caching
```python
result = await adapter(pdb_id, use_cache=False)
```

## Error Handling

```python
result = await adapter("INVALID_ID")

if result.success:
    # Process data
    data = result.data
else:
    # Handle error
    print(f"Error: {result.error}")
```

## Complete Example

```python
import asyncio
from adapters.sabdab import SAbDabAdapter

async def main():
    adapter = SAbDabAdapter()

    # Search for COVID antibodies
    result = await adapter({
        "antigen": "spike",
        "species": "human",
        "resolution_max": 3.0
    }, max_results=10)

    if result.success:
        print(f"Found {result.data['total_results']} antibodies")

        for ab in result.data['antibodies'][:5]:
            print(f"{ab['pdb_id']}: {ab['antigen']} "
                  f"({ab['resolution']} Å)")

asyncio.run(main())
```

## Metadata

```python
# Get adapter information
metadata = adapter.get_metadata()

print(f"Database: {metadata['database']['name']}")
print(f"URL: {metadata['database']['url']}")
```

## Tips

1. **Start broad**: Begin with loose filters, then narrow down
2. **Quality first**: Use resolution < 3.0 Å for reliable structures
3. **Cache results**: Default caching speeds up repeated queries
4. **Batch lookups**: Process multiple PDB IDs efficiently
5. **Species matters**: Filter by species for therapeutic applications

## Common Antigens

- `spike` - SARS-CoV-2 spike protein
- `Her2` - Human epidermal growth factor receptor 2
- `EGFR` - Epidermal growth factor receptor
- `PD-1` - Programmed cell death protein 1
- `CD20` - B-lymphocyte antigen

## Resources

- **Database**: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/
- **Paper**: Dunbar et al. (2014) Nucleic Acids Research
- **Full README**: See `README.md` in this directory
- **Examples**: See `example_usage.py` for detailed examples

## Troubleshooting

### Connection Timeout
```python
adapter = SAbDabAdapter(config={"timeout": 120})
```

### Too Many Results
```python
result = await adapter(params, max_results=50)
```

### No Results Found
- Check spelling of antigen name
- Broaden quality filters
- Try different antibody types
- Check if data exists in SAbDab

## Next Steps

- 📖 Read full documentation in `README.md`
- 💻 Run `example_usage.py` for detailed examples
- 🔬 Explore SAbDab website for available data
- 🧪 Integrate with your PharmForge pipeline

---

**Quick Reference Card**

```python
# Import
from adapters.sabdab import SAbDabAdapter
adapter = SAbDabAdapter()

# PDB lookup
await adapter("7BWJ")

# Antigen search
await adapter({"antigen": "spike"})

# Quality filter
await adapter({"antigen": "Her2", "resolution_max": 2.5})

# Nanobodies
await adapter({"ab_type": "VHH"})
```
