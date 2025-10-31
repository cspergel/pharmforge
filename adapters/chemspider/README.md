# ChemSpider Adapter

**Royal Society of Chemistry Chemical Database**

ChemSpider is a free chemical structure database aggregating data from over 100 sources including PubChem, ChEMBL, DrugBank, and more. It provides comprehensive chemical information, identifiers, properties, and safety data.

## Features

- **Multi-format search**: SMILES, InChI, InChIKey, name, or molecular formula
- **Aggregated data**: Pulls information from 100+ chemical databases
- **Cross-database validation**: Verify compounds across multiple sources
- **Comprehensive identifiers**: Get all synonyms and trade names
- **Property lookup**: Molecular weight, LogP, and calculated properties
- **Safety information**: Access regulatory and safety data (when available)

## Installation

No additional packages required - uses standard `aiohttp` for API calls.

```bash
pip install aiohttp
```

## API Key Setup

ChemSpider requires a free API key for full functionality:

1. **Register for free**: https://developer.rsc.org/
2. **Get your API key**: Navigate to your profile → API Keys
3. **Set environment variable**:

```bash
# Linux/Mac
export CHEMSPIDER_API_KEY="your-api-key-here"

# Windows
set CHEMSPIDER_API_KEY=your-api-key-here

# Or in .env file
CHEMSPIDER_API_KEY=your-api-key-here
```

4. **Verify setup**:

```python
import os
print(os.getenv("CHEMSPIDER_API_KEY"))  # Should print your key
```

## Usage

### Basic SMILES Search

```python
from adapters.chemspider import ChemSpiderAdapter
import asyncio

async def search_aspirin():
    adapter = ChemSpiderAdapter()

    # Search by SMILES
    result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

    if result.success:
        print(f"Found {result.data['num_results']} compounds")
        for compound in result.data['results']:
            print(f"ChemSpider ID: {compound['chemspider_id']}")
            print(f"Name: {compound['common_name']}")
            print(f"Formula: {compound['formula']}")
            print(f"Molecular Weight: {compound['molecular_weight']}")
            print(f"Synonyms: {', '.join(compound['synonyms'][:5])}")
            print(f"Data Sources: {', '.join(compound['data_sources'][:5])}")
    else:
        print(f"Error: {result.error}")

asyncio.run(search_aspirin())
```

### Advanced Search Options

```python
async def advanced_search():
    adapter = ChemSpiderAdapter()

    # Search by different identifier types
    search_params = {
        "query": "aspirin",
        "query_type": "name",  # name, smiles, inchi, inchikey, formula
        "include_properties": True,
        "include_synonyms": True,
        "max_results": 5
    }

    result = await adapter.execute(search_params)

    if result.success:
        for compound in result.data['results']:
            print(f"\nCompound: {compound['common_name']}")
            print(f"SMILES: {compound['smiles']}")
            print(f"InChIKey: {compound['inchikey']}")

            # Properties
            if compound['properties']:
                props = compound['properties']
                print(f"ALogP: {props.get('alogp')}")
                print(f"XLogP: {props.get('xlogp')}")

asyncio.run(advanced_search())
```

### Search by Molecular Formula

```python
async def search_by_formula():
    adapter = ChemSpiderAdapter()

    result = await adapter.execute({
        "query": "C9H8O4",  # Aspirin formula
        "query_type": "formula"
    })

    if result.success:
        print(f"Found {result.data['num_results']} compounds with formula C9H8O4")
        for compound in result.data['results']:
            print(f"- {compound['common_name']} (ID: {compound['chemspider_id']})")

asyncio.run(search_by_formula())
```

### Search by InChIKey

```python
async def search_by_inchikey():
    adapter = ChemSpiderAdapter()

    # Aspirin InChIKey
    result = await adapter.execute({
        "query": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "query_type": "inchikey"
    })

    if result.success and result.data['results']:
        compound = result.data['results'][0]
        print(f"Compound: {compound['common_name']}")
        print(f"Aggregated from {len(compound['data_sources'])} databases")
        print(f"Databases: {compound['data_sources']}")

asyncio.run(search_by_inchikey())
```

### Chemical Validation

```python
async def validate_structure():
    """
    Use ChemSpider to validate a chemical structure across multiple databases
    """
    adapter = ChemSpiderAdapter()

    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    result = await adapter.execute(smiles)

    if result.success and result.data['results']:
        compound = result.data['results'][0]
        num_sources = len(compound['data_sources'])

        print(f"Structure validated across {num_sources} databases")

        if num_sources >= 5:
            print("✓ High confidence - found in multiple authoritative sources")
        elif num_sources >= 2:
            print("○ Medium confidence - found in some sources")
        else:
            print("! Low confidence - limited sources")

        # Check for name consistency
        synonyms = compound['synonyms']
        print(f"\nKnown as: {', '.join(synonyms[:10])}")

asyncio.run(validate_structure())
```

### Batch Processing

```python
async def batch_search():
    adapter = ChemSpiderAdapter()

    compounds = [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    ]

    results = []
    for smiles in compounds:
        result = await adapter.execute(smiles)
        if result.success and result.data['results']:
            compound = result.data['results'][0]
            results.append({
                "smiles": smiles,
                "name": compound['common_name'],
                "chemspider_id": compound['chemspider_id'],
                "sources": len(compound['data_sources'])
            })
        await asyncio.sleep(1)  # Respect rate limits

    for r in results:
        print(f"{r['name']}: {r['sources']} sources (ID: {r['chemspider_id']})")

asyncio.run(batch_search())
```

## Output Format

### Successful Search Result

```python
{
    "results": [
        {
            "chemspider_id": "2157",
            "common_name": "Aspirin",
            "formula": "C9H8O4",
            "molecular_weight": 180.157,
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "synonyms": [
                "Aspirin",
                "Acetylsalicylic acid",
                "2-Acetoxybenzoic acid",
                "ASA",
                "Ecotrin"
            ],
            "properties": {
                "alogp": 1.23,
                "xlogp": 1.19,
                "molecular_formula": "C9H8O4",
                "nominal_mass": 180,
                "average_mass": 180.157,
                "monoisotopic_mass": 180.042259
            },
            "data_sources": [
                "PubChem",
                "ChEMBL",
                "DrugBank",
                "KEGG",
                "HMDB"
            ]
        }
    ],
    "num_results": 1,
    "total_found": 1,
    "warnings": []
}
```

### No Results Found (Not an Error)

```python
{
    "results": [],
    "num_results": 0,
    "warnings": ["No ChemSpider entries found for smiles: INVALID_SMILES"]
}
```

## Configuration

Default configuration:

```python
{
    "rate_limit_delay": 0.5,  # 2 requests/second
    "timeout": 60,            # 60 second timeout
    "max_results": 10         # Return up to 10 results
}
```

Override in initialization:

```python
adapter = ChemSpiderAdapter()
adapter.config["max_results"] = 20
adapter.config["rate_limit_delay"] = 1.0  # More conservative
```

## Rate Limiting

ChemSpider has rate limits:
- **Free tier**: ~2 requests/second recommended
- **Built-in protection**: Adapter includes 0.5s delay between requests
- **Batch processing**: Add additional sleep between compounds

```python
# Good practice for batch processing
for smiles in compound_list:
    result = await adapter.execute(smiles)
    # Process result...
    await asyncio.sleep(1)  # Extra safety margin
```

## Error Handling

```python
async def robust_search():
    adapter = ChemSpiderAdapter()

    try:
        result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

        if not result.success:
            print(f"Search failed: {result.error}")
            return None

        if result.data['num_results'] == 0:
            print("No results found")
            return None

        if result.data['warnings']:
            print(f"Warnings: {result.data['warnings']}")

        return result.data['results'][0]

    except Exception as e:
        print(f"Unexpected error: {e}")
        return None
```

## Common Use Cases

### 1. Chemical Validation
Verify a structure exists in authoritative databases

```python
result = await adapter.execute("YOUR_SMILES_HERE")
if result.success and len(result.data['results'][0]['data_sources']) >= 3:
    print("Structure validated!")
```

### 2. Name to Structure Conversion

```python
result = await adapter.execute({
    "query": "ibuprofen",
    "query_type": "name"
})
smiles = result.data['results'][0]['smiles']
```

### 3. Cross-Database Lookup
Find all identifiers and synonyms for a compound

```python
result = await adapter.execute(smiles)
compound = result.data['results'][0]
print(f"Found in: {compound['data_sources']}")
print(f"Also known as: {compound['synonyms']}")
```

### 4. Property Aggregation
Get properties from multiple sources

```python
result = await adapter.execute(smiles)
props = result.data['results'][0]['properties']
print(f"LogP: {props['alogp']} (ChemSpider aggregate)")
```

## Caching

Results are automatically cached by the PharmForge adapter framework:

```python
# First call - hits API
result1 = await adapter(smiles, use_cache=True)
print(result1.cache_hit)  # False

# Second call - from cache
result2 = await adapter(smiles, use_cache=True)
print(result2.cache_hit)  # True
```

## Integration with PharmForge Pipeline

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()
pipeline.add_step("chemspider_lookup", {
    "adapter": "chemspider",
    "params": {
        "query_type": "smiles",
        "include_properties": True
    }
})

result = await pipeline.execute("CC(=O)Oc1ccccc1C(=O)O")
```

## Advantages Over Other Databases

1. **Aggregation**: Pulls data from 100+ sources (PubChem only has PubChem data)
2. **Validation**: Cross-reference across multiple authoritative databases
3. **Completeness**: More synonyms and identifiers than single-source databases
4. **Regulatory data**: Includes safety and regulatory information
5. **Free access**: No cost for API access (just requires registration)

## Limitations

1. **Rate limits**: Free tier has modest rate limits
2. **API complexity**: Async polling pattern for search results
3. **Latency**: Multiple database queries can be slower than single-source
4. **API key required**: Need to register (though it's free)

## Troubleshooting

### Authentication Errors

```
Error: ChemSpider: Authentication failed - check API key
```

**Solution**: Verify your API key is set correctly:

```python
import os
print(os.getenv("CHEMSPIDER_API_KEY"))
```

### No Results Found

If searches return 0 results:
- Try different query types (name vs SMILES vs InChI)
- Check for typos in chemical names
- Verify SMILES string is valid

### Timeout Errors

If requests timeout:
- Increase timeout in config
- Check your internet connection
- ChemSpider may be experiencing issues

```python
adapter.config["timeout"] = 120  # Increase to 2 minutes
```

## API Reference

### ChemSpider API Documentation
- **Developer Portal**: https://developer.rsc.org/
- **API Docs**: https://developer.rsc.org/compounds-v1/apis
- **Terms of Use**: https://www.rsc.org/terms-conditions/

### Adapter Methods

**`execute(input_data, **kwargs)`**
- Primary method for searching ChemSpider
- Input: SMILES string or dict with query parameters
- Returns: AdapterResult with compound data

**`validate_input(input_data)`**
- Validates input format
- Returns: Boolean

## Contributing

To improve this adapter:

1. Add support for more search types (molecular weight range, etc.)
2. Implement safety/regulatory data extraction
3. Add structure similarity search
4. Improve error messages and retry logic

## License

This adapter is part of PharmForge (MIT License).
ChemSpider data is provided by the Royal Society of Chemistry.
See: https://www.rsc.org/terms-conditions/

## Support

- **ChemSpider Issues**: https://www.chemspider.com/Contact.aspx
- **PharmForge Issues**: [GitHub Issues](https://github.com/yourusername/pharmforge/issues)
- **API Status**: https://status.rsc.org/

## Version History

- **v1.0.0** (2025-10-30): Initial implementation
  - SMILES, InChI, InChIKey, name, and formula search
  - Property and synonym lookup
  - Multi-source data aggregation
  - Automatic caching support
