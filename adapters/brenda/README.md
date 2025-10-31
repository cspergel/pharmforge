# BRENDA Adapter

**Comprehensive Enzyme Information Database**

The BRENDA adapter provides access to the BRaunschweig ENzyme DAtabase, the world's most comprehensive enzyme information system with over 7 million kinetic parameters.

## Overview

BRENDA contains:
- Enzyme kinetic parameters (Km, Kcat, Ki, Kd)
- Substrate specificity data
- Inhibitor information
- Organism-specific enzyme data
- Reaction mechanisms
- Cofactor requirements
- Post-translational modifications
- pH and temperature optima

## Installation

No additional installation required. The adapter uses standard Python libraries (`aiohttp`, `xml`).

```bash
# All dependencies are in requirements.txt
pip install aiohttp
```

## API Access

BRENDA offers several access methods:

### 1. Web Interface (Free)
- URL: https://www.brenda-enzymes.org/
- Access: Unrestricted browsing
- Best for: Manual queries and exploration

### 2. REST API (Limited)
- URL: https://www.brenda-enzymes.org/enzyme.php
- Access: Free, rate-limited
- Best for: Simple queries

### 3. SOAP API (Full Access)
- URL: https://www.brenda-enzymes.org/soap.php
- Access: Requires free registration
- Best for: Programmatic batch queries
- Registration: https://www.brenda-enzymes.org/soap.php

### Registration (Recommended)

For full API access:

1. Visit https://www.brenda-enzymes.org/soap.php
2. Register for a free account (academic or commercial)
3. Receive API credentials (email + password)
4. Set environment variables:

```bash
export BRENDA_EMAIL="your.email@university.edu"
export BRENDA_PASSWORD="your_password"
```

## Usage

### Basic Query (EC Number)

```python
from adapters.brenda import BRENDAAdapter
import asyncio

async def query_enzyme():
    adapter = BRENDAAdapter()

    # Query by EC number
    result = await adapter.execute("1.1.1.1")

    if result.success:
        data = result.data
        print(f"Enzyme: {data['enzyme']['common_name']}")
        print(f"EC: {data['enzyme']['ec_number']}")
        print(f"Reaction: {data['enzyme']['reaction']}")
        print(f"\nKinetic parameters: {len(data['kinetics'])}")
        print(f"Inhibitors: {len(data['inhibitors'])}")

        # Print Km values
        for kinetic in data['kinetics']:
            if kinetic['parameter'] == 'km':
                print(f"  Km = {kinetic['value']} {kinetic['unit']} ({kinetic['substrate']})")

asyncio.run(query_enzyme())
```

### Advanced Query (Full Parameters)

```python
async def advanced_query():
    adapter = BRENDAAdapter()

    query = {
        "query_type": "ec_number",
        "query": "1.1.1.1",
        "parameters": ["km", "kcat", "ki", "kd"],
        "organisms": ["Homo sapiens", "Saccharomyces cerevisiae"],
        "include_references": True,
        "max_results": 100
    }

    result = await adapter.execute(query)

    if result.success:
        enzyme = result.data['enzyme']
        kinetics = result.data['kinetics']
        inhibitors = result.data['inhibitors']

        # Filter human-specific data
        human_kinetics = [k for k in kinetics if k['organism'] == 'Homo sapiens']

        print(f"Human {enzyme['common_name']} kinetics:")
        for k in human_kinetics:
            print(f"  {k['parameter'].upper()}: {k['value']} {k['unit']}")
            print(f"    Substrate: {k['substrate']}")
            print(f"    Conditions: pH {k['conditions']['pH']}, {k['conditions']['temperature']}°C")
            print(f"    Reference: {k['reference']}\n")

asyncio.run(advanced_query())
```

### Query by Enzyme Name

```python
async def search_by_name():
    adapter = BRENDAAdapter()

    query = {
        "query_type": "enzyme_name",
        "query": "alcohol dehydrogenase",
        "parameters": ["km", "kcat"],
        "max_results": 50
    }

    result = await adapter.execute(query)

    if result.success:
        print(f"Found: {result.data['enzyme']['ec_number']}")
```

## Output Format

### Complete Response Structure

```python
{
    "enzyme": {
        "ec_number": "1.1.1.1",
        "systematic_name": "Alcohol:NAD+ oxidoreductase",
        "common_name": "Alcohol dehydrogenase",
        "reaction": "Alcohol + NAD+ = Aldehyde + NADH + H+",
        "substrates": ["Ethanol", "Methanol", "Propanol"],
        "products": ["Acetaldehyde", "Formaldehyde"],
        "cofactors": ["NAD+", "Zn2+"]
    },
    "kinetics": [
        {
            "parameter": "km",
            "value": 0.5,
            "unit": "mM",
            "substrate": "Ethanol",
            "organism": "Homo sapiens",
            "conditions": {
                "pH": 7.4,
                "temperature": 37
            },
            "reference": "PMID:12345678"
        }
    ],
    "inhibitors": [
        {
            "compound": "Disulfiram",
            "ki_value": 10.0,
            "unit": "nM",
            "inhibition_type": "competitive",
            "organism": "Homo sapiens",
            "reference": "PMID:34567890"
        }
    ],
    "num_kinetic_records": 25,
    "num_inhibitor_records": 5,
    "warnings": []
}
```

## Use Cases

### 1. Drug Target Validation

```python
async def validate_enzyme_target():
    """Check if enzyme is druggable based on known inhibitors"""
    adapter = BRENDAAdapter()

    result = await adapter.execute({
        "query": "3.4.23.15",  # Renin
        "parameters": ["ki"],
        "organisms": ["Homo sapiens"]
    })

    if result.success:
        inhibitors = result.data['inhibitors']
        potent_inhibitors = [i for i in inhibitors if i['ki_value'] < 100]  # < 100 nM

        print(f"Found {len(potent_inhibitors)} potent inhibitors")
        if potent_inhibitors:
            print("Target is druggable!")
```

### 2. Enzyme Kinetics for Drug Design

```python
async def analyze_enzyme_kinetics():
    """Get kinetic parameters for inhibitor design"""
    adapter = BRENDAAdapter()

    result = await adapter.execute({
        "query": "1.14.13.39",  # CYP2D6
        "parameters": ["km", "kcat"],
        "organisms": ["Homo sapiens"]
    })

    if result.success:
        kinetics = result.data['kinetics']

        # Calculate catalytic efficiency (Kcat/Km)
        for k in kinetics:
            if k['parameter'] == 'km':
                km = k['value']
                # Find matching Kcat
                kcat_record = next(
                    (x for x in kinetics
                     if x['parameter'] == 'kcat'
                     and x['substrate'] == k['substrate']),
                    None
                )

                if kcat_record:
                    kcat = kcat_record['value']
                    efficiency = kcat / km
                    print(f"{k['substrate']}: Efficiency = {efficiency:.2f} mM⁻¹s⁻¹")
```

### 3. Species-Specific Enzyme Differences

```python
async def compare_species():
    """Compare enzyme kinetics across species"""
    adapter = BRENDAAdapter()

    organisms = ["Homo sapiens", "Mus musculus", "Rattus norvegicus"]

    results = {}
    for org in organisms:
        result = await adapter.execute({
            "query": "1.1.1.1",
            "parameters": ["km"],
            "organisms": [org]
        })

        if result.success:
            km_values = [k['value'] for k in result.data['kinetics'] if k['parameter'] == 'km']
            if km_values:
                results[org] = sum(km_values) / len(km_values)

    print("Average Km by species:")
    for org, km in results.items():
        print(f"  {org}: {km:.2f} mM")
```

### 4. Inhibitor Screening

```python
async def screen_inhibitors():
    """Find competitive inhibitors for drug repurposing"""
    adapter = BRENDAAdapter()

    result = await adapter.execute({
        "query": "3.4.21.5",  # Thrombin
        "parameters": ["ki"],
        "organisms": ["Homo sapiens"]
    })

    if result.success:
        inhibitors = result.data['inhibitors']

        # Filter competitive inhibitors
        competitive = [i for i in inhibitors if i['inhibition_type'] == 'competitive']

        # Sort by potency
        competitive.sort(key=lambda x: x['ki_value'])

        print("Top 5 competitive inhibitors:")
        for i, inh in enumerate(competitive[:5], 1):
            print(f"{i}. {inh['compound']}: Ki = {inh['ki_value']} {inh['unit']}")
```

### 5. Metabolic Pathway Analysis

```python
async def analyze_pathway():
    """Analyze enzymes in a metabolic pathway"""
    adapter = BRENDAAdapter()

    # Glycolysis pathway enzymes
    glycolysis_enzymes = [
        "2.7.1.1",   # Hexokinase
        "5.3.1.9",   # Glucose-6-phosphate isomerase
        "2.7.1.11",  # Phosphofructokinase
        "4.1.2.13",  # Fructose-bisphosphate aldolase
    ]

    for ec in glycolysis_enzymes:
        result = await adapter.execute({
            "query": ec,
            "parameters": ["km", "kcat"],
            "organisms": ["Homo sapiens"]
        })

        if result.success:
            enzyme = result.data['enzyme']
            print(f"\n{enzyme['common_name']} (EC {ec})")
            print(f"Reaction: {enzyme['reaction']}")

            # Show rate-limiting information
            kinetics = result.data['kinetics']
            if kinetics:
                avg_km = sum(k['value'] for k in kinetics if k['parameter'] == 'km') / len([k for k in kinetics if k['parameter'] == 'km'])
                print(f"Average Km: {avg_km:.2f} mM")
```

## API Reference

### Input Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query_type` | str | No | Query type: "ec_number", "enzyme_name", "substrate", "organism" (default: "ec_number") |
| `query` | str | Yes | EC number or search term |
| `parameters` | list | No | Kinetic parameters to fetch: ["km", "kcat", "ki", "kd"] (default: ["km", "kcat", "ki"]) |
| `organisms` | list | No | Filter by organisms (e.g., ["Homo sapiens"]) |
| `include_references` | bool | No | Include PubMed references (default: true) |
| `max_results` | int | No | Maximum results to return (default: 50) |

### Output Fields

#### Enzyme Object
- `ec_number`: EC classification number
- `systematic_name`: IUPAC systematic name
- `common_name`: Recommended common name
- `reaction`: Reaction equation
- `substrates`: List of known substrates
- `products`: List of products
- `cofactors`: Required cofactors

#### Kinetics Array
- `parameter`: Parameter type (km, kcat, ki, kd)
- `value`: Numeric value
- `unit`: Unit of measurement
- `substrate`: Substrate name
- `organism`: Organism source
- `conditions`: pH, temperature, etc.
- `reference`: PubMed ID

#### Inhibitors Array
- `compound`: Inhibitor name/structure
- `ki_value`: Inhibition constant
- `unit`: Unit (nM, µM, mM)
- `inhibition_type`: competitive, non-competitive, uncompetitive
- `organism`: Organism tested
- `reference`: PubMed ID

## Configuration

### Rate Limiting

BRENDA requests are rate-limited to be respectful:

```python
adapter = BRENDAAdapter()
adapter.config['rate_limit_delay'] = 1.0  # 1 second between requests
```

### Timeout

Adjust timeout for slow connections:

```python
adapter.config['timeout'] = 120  # 2 minutes
```

### Result Limits

Control maximum results:

```python
adapter.config['max_results'] = 100
```

## Error Handling

```python
result = await adapter.execute("invalid_ec")

if not result.success:
    print(f"Error: {result.error}")
    print(f"Metadata: {result.metadata}")
```

Common errors:
- `"Enzyme not found for EC number"` - Invalid or non-existent EC number
- `"Query type not yet implemented"` - Feature under development
- `"Timeout"` - Network/server timeout
- `"Invalid input"` - Malformed query parameters

## Performance Tips

1. **Use Caching**: Results are automatically cached
2. **Batch Queries**: Query multiple enzymes sequentially with rate limiting
3. **Filter by Organism**: Reduces data volume
4. **Limit Parameters**: Only request needed kinetic parameters
5. **Set max_results**: Prevent large response sizes

## Integration Examples

### With PharmForge Pipeline

```python
from backend.core.pipeline import Pipeline
from adapters.brenda import BRENDAAdapter

async def enzyme_pipeline():
    pipeline = Pipeline()

    # Add BRENDA adapter to pipeline
    brenda = BRENDAAdapter()
    pipeline.add_step("enzyme_kinetics", brenda)

    # Run pipeline
    result = await pipeline.execute(
        input_data="1.1.1.1",
        steps=["enzyme_kinetics"]
    )

    return result
```

### With Target Validation

```python
async def validate_target_with_brenda():
    from adapters.uniprot import UniProtAdapter
    from adapters.brenda import BRENDAAdapter

    # Get protein info
    uniprot = UniProtAdapter()
    protein_result = await uniprot.execute("P00326")  # ADH1A_HUMAN

    # Get enzyme kinetics
    brenda = BRENDAAdapter()
    enzyme_result = await brenda.execute({
        "query": "1.1.1.1",
        "parameters": ["ki"],
        "organisms": ["Homo sapiens"]
    })

    # Combine data
    if protein_result.success and enzyme_result.success:
        print(f"Protein: {protein_result.data['protein_name']}")
        print(f"Enzyme: {enzyme_result.data['enzyme']['common_name']}")
        print(f"Known inhibitors: {len(enzyme_result.data['inhibitors'])}")
```

## Limitations

### Current Implementation

This version provides:
- ✅ EC number queries
- ✅ Basic enzyme information
- ✅ Kinetic parameters (Km, Kcat, Ki)
- ✅ Inhibitor data
- ✅ Organism filtering

Not yet implemented:
- ⏳ Enzyme name search
- ⏳ Substrate-based search
- ⏳ Advanced SOAP API integration
- ⏳ Full reference extraction
- ⏳ pH/temperature optimization curves

### Data Coverage

- Example data is used for demonstration
- Production deployment should integrate with BRENDA SOAP API
- Some kinetic parameters may not be available for all enzymes

## Resources

- **BRENDA Homepage**: https://www.brenda-enzymes.org/
- **SOAP API Documentation**: https://www.brenda-enzymes.org/soap.php
- **EC Number Database**: https://enzyme.expasy.org/
- **Kinetics Tutorial**: https://www.brenda-enzymes.org/kinetics.php

## Citation

If you use BRENDA data in publications, please cite:

> Chang A, Jeske L, Ulbrich S, Hofmann J, Koblitz J, Schomburg I, Neumann-Schaal M, Jahn D, Schomburg D (2021)
> BRENDA, the ELIXIR core data resource in 2021: new developments and updates.
> Nucleic Acids Res 49:D498-D508

## License

BRENDA data is free for academic use. Commercial use requires a license.
Contact: brenda@tu-braunschweig.de

## Support

For issues or questions:
- PharmForge Issues: [GitHub Issues](https://github.com/your-repo/pharmforge/issues)
- BRENDA Support: brenda@tu-braunschweig.de
