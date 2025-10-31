# Database Adapters - Quick Reference Guide

## Quick Start

```python
import asyncio
from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter
from adapters.surechembl.adapter import SureChEMBLAdapter
from adapters.pubchem.adapter_enhanced import PubChemEnhancedAdapter
from adapters.chembl.adapter_enhanced import ChEMBLEnhancedAdapter

# Initialize
zinc = ZINCFragmentsAdapter()
surechembl = SureChEMBLAdapter()
pubchem = PubChemEnhancedAdapter()
chembl = ChEMBLEnhancedAdapter()

# Use
result = await zinc(smiles, search_type="similarity")
```

## Common Tasks

### 1. Find Purchasable Building Blocks

```python
# ZINC search
result = await zinc(
    smiles,
    search_type="similarity",
    similarity_threshold=0.8,
    get_purchasability=True,
    mw_max=250
)

# Check vendors
if result.success:
    for frag in result.data['fragments']:
        purch = frag.get('purchasability', {})
        if purch.get('purchasable'):
            print(f"{frag['zinc_id']}: {purch['vendor_count']} vendors")
```

### 2. Check Patent Landscape

```python
# SureChEMBL search
result = await surechembl(
    smiles,
    mode="structure_search",
    search_type="similarity",
    similarity_threshold=0.85
)

if result.success:
    num_patents = result.data['num_patents']
    if num_patents > 10:
        print("⚠️  High patent risk - consider alternatives")
```

### 3. Get Compound Properties

```python
# PubChem lookup
result = await pubchem(smiles, mode="properties")

if result.success:
    data = result.data
    print(f"MW: {data['molecular_weight']}")
    print(f"LogP: {data['logp']}")
    print(f"TPSA: {data['tpsa']}")
```

### 4. Find Bioactive Compounds

```python
# ChEMBL search
result = await chembl(
    smiles,
    mode="bioactivity",
    activity_type="IC50",
    max_value=100.0  # IC50 ≤ 100 nM
)

if result.success:
    print(f"Best IC50: {result.data['best_ic50_nm']} nM")
    print(f"Targets: {result.data['num_targets']}")
```

### 5. Similarity Search

```python
# PubChem similarity
result = await pubchem(
    smiles,
    mode="similarity",
    similarity_threshold=85,
    max_results=50
)

# ChEMBL similarity (bioactive space)
result = await chembl(
    smiles,
    mode="similarity",
    similarity=80
)
```

### 6. Check Drug Status

```python
# ChEMBL drug lookup
result = await chembl(smiles, mode="drug")

if result.success:
    if result.data['is_drug']:
        drug_info = result.data['drug_info']
        print(f"Max Phase: {drug_info['max_phase']}")
        print(f"First Approval: {drug_info['first_approval']}")
```

## Search Type Reference

### ZINC Fragments

| Search Type | Description | Parameters |
|------------|-------------|------------|
| `similarity` | Find similar fragments | `similarity_threshold` (0-1) |
| `substructure` | Contains substructure | - |
| `exact` | Exact match | - |

### SureChEMBL

| Mode | Description | Parameters |
|------|-------------|------------|
| `structure_search` | Search patents by structure | `search_type`, `similarity_threshold` |
| `patent_lookup` | Get patent details | `patent_id` |
| `extract_compounds` | Extract compounds from patent | `patent_id`, `max_results` |

### PubChem Enhanced

| Mode | Description | Returns |
|------|-------------|---------|
| `properties` | Molecular properties | MW, LogP, TPSA, HBD, HBA, etc. |
| `similarity` | Similar compounds | List of CIDs |
| `substructure` | Substructure search | List of CIDs |
| `bioassays` | Bioassay data | Assay records |
| `vendors` | Commercial suppliers | Vendor list |
| `patents` | Patent mentions | Patent IDs |

### ChEMBL Enhanced

| Mode | Description | Returns |
|------|-------------|---------|
| `bioactivity` | Activity data | IC50, Ki, Kd, targets |
| `target` | Target info | Type, organism, components |
| `drug` | Drug status | Max phase, approval date |
| `mechanism` | Mechanism of action | MOA, action type |
| `similarity` | Similar bioactive compounds | Molecules |
| `substructure` | Substructure in bioactive space | Molecules |

## Property Filters

### ZINC Fragments

```python
result = await zinc(
    smiles,
    mw_max=300,      # Maximum molecular weight
    logp_max=3.0,    # Maximum LogP
    hbd_max=3,       # Maximum H-bond donors
    hba_max=3        # Maximum H-bond acceptors
)
```

### ChEMBL Activities

```python
result = await chembl(
    smiles,
    mode="bioactivity",
    activity_type="IC50",  # or "Ki", "Kd", "EC50"
    max_value=100.0,       # Maximum value (nM)
    min_value=0.1,         # Minimum value (nM)
    target_type="PROTEIN"  # Target type filter
)
```

## Retrosynthesis Workflows

### Building Block Validation

```python
async def validate_building_block(smiles):
    # Check ZINC
    zinc_result = await zinc(smiles, search_type="exact", get_purchasability=True)

    # Check PubChem
    pubchem_result = await pubchem(smiles, mode="vendors")

    # Aggregate
    zinc_avail = zinc_result.success and zinc_result.data['fragments']
    pubchem_avail = pubchem_result.success and pubchem_result.data['num_vendors'] > 0

    return zinc_avail or pubchem_avail
```

### IP Risk Assessment

```python
async def check_ip_risk(smiles):
    result = await surechembl(
        smiles,
        mode="structure_search",
        search_type="exact"
    )

    if result.success:
        num_patents = result.data['num_patents']

        if num_patents == 0:
            return "LOW"
        elif num_patents < 5:
            return "MEDIUM"
        else:
            return "HIGH"

    return "UNKNOWN"
```

### Bioactivity Check

```python
async def get_bioactivity_score(smiles):
    result = await chembl(smiles, mode="bioactivity")

    if result.success:
        return {
            'num_activities': result.data['num_activities'],
            'num_targets': result.data['num_targets'],
            'best_ic50': result.data['best_ic50_nm']
        }

    return None
```

### Comprehensive Evaluation

```python
async def evaluate_route(intermediates):
    scores = []

    for smiles in intermediates:
        # Parallel queries
        tasks = [
            validate_building_block(smiles),
            check_ip_risk(smiles),
            get_bioactivity_score(smiles)
        ]

        available, ip_risk, bioactivity = await asyncio.gather(*tasks)

        score = {
            'smiles': smiles,
            'available': available,
            'ip_risk': ip_risk,
            'bioactivity': bioactivity
        }

        scores.append(score)

    return scores
```

## Batch Operations

### Parallel Queries

```python
# Query multiple compounds at once
smiles_list = ["smiles1", "smiles2", "smiles3"]

tasks = [zinc(s, search_type="exact") for s in smiles_list]
results = await asyncio.gather(*tasks)

for result in results:
    if result.success:
        print(f"Found: {len(result.data['fragments'])} fragments")
```

### Multi-Database Profiling

```python
# Get comprehensive data from all databases
async def profile_compound(smiles):
    tasks = [
        zinc(smiles, search_type="similarity", get_purchasability=True),
        pubchem(smiles, mode="properties"),
        pubchem(smiles, mode="vendors"),
        chembl(smiles, mode="bioactivity"),
        surechembl(smiles, mode="structure_search")
    ]

    results = await asyncio.gather(*tasks)

    return {
        'zinc': results[0].data if results[0].success else None,
        'properties': results[1].data if results[1].success else None,
        'vendors': results[2].data if results[2].success else None,
        'bioactivity': results[3].data if results[3].success else None,
        'patents': results[4].data if results[4].success else None
    }
```

## Configuration Quick Reference

### Rate Limiting

```python
adapter.config['rate_limit_delay'] = 1.0  # seconds
```

### Timeouts

```python
adapter.config['timeout'] = 90  # seconds
```

### Result Limits

```python
adapter.config['max_results'] = 100
adapter.config['max_activities'] = 200  # ChEMBL only
```

### Similarity Thresholds

```python
# ZINC (Tanimoto coefficient)
similarity_threshold=0.7  # 70% similarity (0.0-1.0)

# PubChem
similarity_threshold=85   # 85% similarity (0-100)

# ChEMBL
similarity=80             # 80% similarity (0-100)

# SureChEMBL
similarity_threshold=0.8  # 80% similarity (0.0-1.0)
```

## Error Handling

```python
result = await adapter(smiles)

if result.success:
    # Use result.data
    data = result.data
else:
    # Handle error
    print(f"Error: {result.error}")

    # Check metadata
    if result.metadata:
        print(f"Metadata: {result.metadata}")
```

## Caching

```python
# Automatic caching (default)
result1 = await adapter(smiles)  # Hits API
result2 = await adapter(smiles)  # From cache

# Disable caching
result3 = await adapter(smiles, use_cache=False)

# Check if cached
if result.cache_hit:
    print("Result from cache")
```

## Common Patterns

### Fragment-Based Design

```python
# 1. Find fragments similar to target
zinc_result = await zinc(
    target_smiles,
    search_type="similarity",
    similarity_threshold=0.6,
    mw_max=200
)

# 2. Filter for purchasable
purchasable = [
    f for f in zinc_result.data['fragments']
    if f.get('purchasability', {}).get('purchasable')
]

# 3. Check bioactivity (privileged scaffolds)
for frag in purchasable:
    chembl_result = await chembl(frag['smiles'], mode="bioactivity")
    if chembl_result.data['num_activities'] > 10:
        print(f"Privileged scaffold: {frag['zinc_id']}")
```

### Lead Optimization

```python
# 1. Find analogs in PubChem
pubchem_result = await pubchem(
    lead_smiles,
    mode="similarity",
    similarity_threshold=85
)

# 2. Get bioactivity for analogs from ChEMBL
for cid in pubchem_result.data['cids'][:20]:
    # Convert CID to SMILES first, then:
    chembl_result = await chembl(
        smiles,
        mode="bioactivity",
        activity_type="IC50",
        max_value=1000.0
    )

    if chembl_result.success:
        print(f"CID {cid}: {chembl_result.data['best_ic50_nm']} nM")
```

### IP Landscape Analysis

```python
# 1. Search patent space
patent_result = await surechembl(
    compound_smiles,
    mode="structure_search",
    search_type="similarity",
    similarity_threshold=0.9
)

# 2. Analyze applicants
if patent_result.success:
    applicants = patent_result.data['top_applicants']

    # Check concentration
    total = patent_result.data['num_patents']
    top3 = sum(a['count'] for a in applicants[:3])
    concentration = (top3 / total * 100) if total > 0 else 0

    if concentration > 70:
        print("⚠️  High concentration - crowded IP space")
```

## Performance Tips

1. **Use caching** for repeated queries
2. **Batch operations** with `asyncio.gather()`
3. **Apply filters** early to reduce result sizes
4. **Respect rate limits** to avoid throttling
5. **Parallel databases** for comprehensive profiling
6. **Strategic queries** - start broad, narrow down

## Troubleshooting

### Timeout Errors

```python
adapter.config['timeout'] = 120  # Increase timeout
```

### Rate Limit Errors

```python
adapter.config['rate_limit_delay'] = 2.0  # Increase delay
```

### No Results

```python
# Lower similarity threshold
result = await adapter(smiles, similarity_threshold=0.6)

# Broaden search
result = await adapter(smiles, search_type="similarity", max_results=200)
```

### API Errors

```python
# Check result status
if not result.success:
    print(f"Error: {result.error}")

    # Retry with different parameters
    result = await adapter(smiles, search_type="substructure")
```

## Example Files

- `zinc_fragments/example_usage.py` - Fragment search examples
- `surechembl/example_usage.py` - Patent search examples
- `pubchem/example_enhanced.py` - PubChem features
- `chembl/example_enhanced.py` - Bioactivity queries
- `INTEGRATION_EXAMPLES.py` - Complete workflows

Run:
```bash
python adapters/zinc_fragments/example_usage.py
python adapters/INTEGRATION_EXAMPLES.py
```
