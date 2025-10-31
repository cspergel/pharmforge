# Chemical Database and Compound Library Adapters

Comprehensive adapters for searching chemical databases, finding purchasable compounds, analyzing bioactivity, and checking patent landscapes.

## Overview

These adapters integrate major chemical databases into PharmForge workflows, enabling:
- **Building block discovery** - Find commercially available starting materials
- **Bioactivity analysis** - Assess drug-likeness and target interactions
- **Patent landscape** - Check freedom-to-operate for synthesis routes
- **Vendor sourcing** - Identify suppliers and pricing
- **Fragment-based design** - Discover privileged scaffolds

## Adapters

### 1. ZINC Fragments Adapter (`zinc_fragments/`)

**Database**: ZINC15/ZINC20 fragment library
**API**: http://zinc15.docking.org/

#### Capabilities
- Fragment search by SMILES
- Similarity search (Tanimoto threshold)
- Substructure search
- Purchasability and vendor information
- Property filtering (MW, LogP, HBD, HBA)

#### Example Usage

```python
from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter

adapter = ZINCFragmentsAdapter()

# Similarity search for fragments
result = await adapter(
    "c1ccccc1",  # Query SMILES
    search_type="similarity",
    similarity_threshold=0.8,
    max_results=20,
    mw_max=250,  # Fragment size limit
    get_purchasability=True
)

if result.success:
    fragments = result.data['fragments']
    print(f"Found {len(fragments)} purchasable fragments")

    for frag in fragments:
        purch = frag.get('purchasability', {})
        print(f"  {frag['zinc_id']}: {purch['vendor_count']} vendors")
```

#### Search Types
- **`similarity`** - Find similar fragments (Tanimoto threshold)
- **`substructure`** - Find fragments containing substructure
- **`exact`** - Exact match search

#### Property Filters
- `mw_max` - Maximum molecular weight (Da)
- `logp_max` - Maximum LogP
- `hbd_max` - Maximum H-bond donors
- `hba_max` - Maximum H-bond acceptors

#### Purchasability Data
When `get_purchasability=True`:
```python
{
    'zinc_id': 'ZINC000001234567',
    'purchasable': True,
    'vendor_count': 5,
    'suppliers': ['Sigma-Aldrich', 'TCI', 'Acros'],
    'price_range': '$50-$200',
    'delivery_days': '2-3'
}
```

### 2. SureChEMBL Adapter (`surechembl/`)

**Database**: SureChEMBL patent chemistry database
**API**: https://www.surechembl.org/api/

#### Capabilities
- Patent compound search by structure
- Similarity search in patent literature
- Patent ID lookup and details
- Extract compounds from patents
- Publication dates and patent families

#### Example Usage

```python
from adapters.surechembl.adapter import SureChEMBLAdapter

adapter = SureChEMBLAdapter()

# Search patents for similar structures
result = await adapter(
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    mode="structure_search",
    search_type="similarity",
    similarity_threshold=0.85,
    max_results=20
)

if result.success:
    print(f"Found in {result.data['num_patents']} patents")
    print(f"Top applicants:")
    for app in result.data['top_applicants'][:5]:
        print(f"  {app['name']}: {app['count']} patents")
```

#### Modes

**`structure_search`** - Search by chemical structure
```python
result = await adapter(
    smiles,
    mode="structure_search",
    search_type="similarity",  # or "substructure", "exact"
    similarity_threshold=0.8
)
```

**`patent_lookup`** - Get patent details
```python
result = await adapter(
    "US20190123456",
    mode="patent_lookup"
)
# Returns: title, abstract, applicants, inventors, compounds
```

**`extract_compounds`** - Extract all compounds from patent
```python
result = await adapter(
    "US20190123456",
    mode="extract_compounds",
    max_results=50
)
# Returns: list of compound structures
```

#### IP Landscape Analysis
```python
{
    'num_patents': 45,
    'num_patent_families': 12,
    'date_range': {
        'earliest': '2010-03-15',
        'latest': '2023-11-20'
    },
    'top_applicants': [
        {'name': 'Pfizer', 'count': 15},
        {'name': 'Merck', 'count': 8}
    ]
}
```

### 3. Enhanced PubChem Adapter (`pubchem/`)

**Database**: PubChem
**API**: https://pubchem.ncbi.nlm.nih.gov/rest/pug/

#### Capabilities
- Compound property lookup (enhanced)
- Similarity search (2D Tanimoto)
- Substructure search
- Bioassay data retrieval
- Vendor and purchasability information
- Patent mentions
- Synonyms and chemical names

#### Example Usage

```python
from adapters.pubchem.adapter_enhanced import PubChemEnhancedAdapter

adapter = PubChemEnhancedAdapter()

# Enhanced property lookup
result = await adapter(smiles, mode="properties")
# Returns: MW, LogP, TPSA, HBD, HBA, complexity, IUPAC name, InChI

# Similarity search
result = await adapter(
    smiles,
    mode="similarity",
    similarity_threshold=90,  # 90% similarity
    max_results=100
)
# Returns: list of similar compound CIDs

# Bioassay data
result = await adapter(smiles, mode="bioassays")
# Returns: bioassay records with activity outcomes

# Vendor lookup
result = await adapter(smiles, mode="vendors")
# Returns: list of commercial suppliers

# Patent mentions
result = await adapter(smiles, mode="patents")
# Returns: patent IDs mentioning compound
```

#### Modes Summary
- **`properties`** - Molecular properties and identifiers
- **`similarity`** - Find similar compounds (CIDs)
- **`substructure`** - Find compounds with substructure (CIDs)
- **`bioassays`** - Retrieve bioassay data
- **`vendors`** - Get commercial suppliers
- **`patents`** - Find patent mentions

#### Comprehensive Profiling
```python
# Get all data types at once
tasks = [
    adapter(smiles, mode="properties"),
    adapter(smiles, mode="bioassays"),
    adapter(smiles, mode="vendors"),
    adapter(smiles, mode="patents")
]
results = await asyncio.gather(*tasks)
```

### 4. Enhanced ChEMBL Adapter (`chembl/`)

**Database**: ChEMBL bioactivity database
**API**: https://www.ebi.ac.uk/chembl/api/

#### Capabilities
- Bioactivity data search with filtering
- Target information and queries
- Drug/clinical candidate lookup
- Mechanism of action data
- Similarity search in bioactive space
- Substructure search
- Advanced activity filtering (IC50, Ki, Kd)

#### Example Usage

```python
from adapters.chembl.adapter_enhanced import ChEMBLEnhancedAdapter

adapter = ChEMBLEnhancedAdapter()

# Bioactivity search with filtering
result = await adapter(
    smiles,
    mode="bioactivity",
    activity_type="IC50",
    max_value=100.0,  # Only IC50 ≤ 100 nM
    target_type="PROTEIN"
)

if result.success:
    print(f"Activities: {result.data['num_activities']}")
    print(f"Targets: {result.data['num_targets']}")
    print(f"Best IC50: {result.data['best_ic50_nm']} nM")
```

#### Modes

**`bioactivity`** - Search bioactivity data
```python
result = await adapter(
    smiles,
    mode="bioactivity",
    activity_type="IC50",  # or "Ki", "Kd", "EC50"
    target_type="PROTEIN",
    max_value=100.0,  # Maximum activity value (nM)
    min_value=0.1     # Minimum activity value (nM)
)
```

**`target`** - Get target information
```python
result = await adapter(
    "TARGET:CHEMBL203",  # EGFR
    mode="target"
)
# Returns: target type, organism, name, components
```

**`drug`** - Check drug status
```python
result = await adapter(
    chembl_id,
    mode="drug"
)
# Returns: max_phase, first_approval, indication_class
```

**`mechanism`** - Get mechanism of action
```python
result = await adapter(
    chembl_id,
    mode="mechanism"
)
# Returns: mechanism_of_action, action_type, target
```

**`similarity`** - Find similar bioactive compounds
```python
result = await adapter(
    smiles,
    mode="similarity",
    similarity=80  # 80% similarity
)
```

**`substructure`** - Substructure search
```python
result = await adapter(
    substructure_smiles,
    mode="substructure"
)
```

#### Advanced Filtering
```python
# Find highly potent compounds
result = await adapter(
    smiles,
    mode="bioactivity",
    activity_type="IC50",
    max_value=10.0,  # Very potent (IC50 < 10 nM)
    relation="<"     # Relation type
)
```

## Retrosynthesis Integration

### Workflow 1: Building Block Validation

Validate commercial availability before synthesis:

```python
from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter
from adapters.pubchem.adapter_enhanced import PubChemEnhancedAdapter

zinc = ZINCFragmentsAdapter()
pubchem = PubChemEnhancedAdapter()

building_blocks = ["c1ccccc1Br", "CC(=O)Cl"]

for smiles in building_blocks:
    # Check ZINC
    zinc_result = await zinc(smiles, search_type="exact", get_purchasability=True)

    # Check PubChem
    pubchem_result = await pubchem(smiles, mode="vendors")

    # Decision
    if zinc_result.success or pubchem_result.data['num_vendors'] > 0:
        print(f"✓ {smiles} - COMMERCIALLY AVAILABLE")
    else:
        print(f"✗ {smiles} - NEEDS SYNTHESIS")
```

### Workflow 2: Fragment-Based Retrosynthesis

Use ZINC fragments to guide disconnections:

```python
from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter
from adapters.chembl.adapter_enhanced import ChEMBLEnhancedAdapter

zinc = ZINCFragmentsAdapter()
chembl = ChEMBLEnhancedAdapter()

# Find fragments similar to target
zinc_result = await zinc(
    target_smiles,
    search_type="similarity",
    similarity_threshold=0.6,
    mw_max=200  # Smaller than target
)

# Check which fragments are bioactive (privileged scaffolds)
for frag in zinc_result.data['fragments']:
    chembl_result = await chembl(frag['smiles'], mode="bioactivity")

    if chembl_result.data['num_activities'] > 10:
        print(f"Privileged scaffold: {frag['zinc_id']}")
        # Use as retrosynthesis starting point
```

### Workflow 3: IP-Aware Retrosynthesis

Check patent landscape for intermediates:

```python
from adapters.surechembl.adapter import SureChEMBLAdapter

surechembl = SureChEMBLAdapter()

intermediates = ["c1ccc(cc1)C(=O)Cl", "CC(C)Cc1ccccc1"]

for smiles in intermediates:
    result = await surechembl(
        smiles,
        mode="structure_search",
        search_type="exact"
    )

    num_patents = result.data['num_patents']

    if num_patents == 0:
        print(f"✓ Clear for use")
    elif num_patents < 5:
        print(f"⚠️  {num_patents} patents - Review recommended")
    else:
        print(f"✗ {num_patents} patents - HIGH RISK")
        # Consider alternative route
```

### Workflow 4: Bioactivity-Guided Synthesis

Prioritize routes with drug-like intermediates:

```python
from adapters.chembl.adapter_enhanced import ChEMBLEnhancedAdapter

chembl = ChEMBLEnhancedAdapter()

routes = {
    "Route A": ["intermediate1_smiles", "intermediate2_smiles"],
    "Route B": ["intermediate3_smiles", "intermediate4_smiles"]
}

for route_name, intermediates in routes.items():
    total_activities = 0

    for smiles in intermediates:
        result = await chembl(smiles, mode="bioactivity")
        total_activities += result.data['num_activities']

    print(f"{route_name}: {total_activities} activities")

# Choose route with more bioactive intermediates
```

### Workflow 5: Comprehensive Route Evaluation

Score routes by multiple criteria:

```python
# For each intermediate, check:
# 1. Purchasability (ZINC + PubChem)
# 2. Patent risk (SureChEMBL)
# 3. Bioactivity (ChEMBL)
# 4. Structural complexity

score = 0

# Purchasability
zinc_result = await zinc(smiles, get_purchasability=True)
if zinc_result.data['fragments'][0]['purchasability']['purchasable']:
    score += 1

# Patent risk
patent_result = await surechembl(smiles, mode="structure_search")
if patent_result.data['num_patents'] < 5:
    score += 1

# Low bioactivity (safer intermediate)
chembl_result = await chembl(smiles, mode="bioactivity")
if chembl_result.data['num_activities'] < 10:
    score += 1

# Simple structure
if len(smiles) < 20:
    score += 1

# Score: 4/4 = Excellent route
```

## Configuration

### Rate Limiting

All adapters respect API rate limits:

```python
adapter = ZINCFragmentsAdapter()
adapter.config['rate_limit_delay'] = 1.0  # 1 second between requests
```

### Timeouts

```python
adapter.config['timeout'] = 90  # 90 seconds
```

### Result Limits

```python
# ZINC
result = await adapter(smiles, max_results=100)

# PubChem
result = await adapter(smiles, mode="similarity", max_results=200)

# ChEMBL
adapter.config['max_activities'] = 500
```

## Caching

All adapters support automatic caching via `AdapterProtocol`:

```python
# First call - hits API
result1 = await adapter(smiles)  # cache_hit=False

# Second call - from cache
result2 = await adapter(smiles)  # cache_hit=True

# Disable caching
result3 = await adapter(smiles, use_cache=False)
```

## Error Handling

```python
result = await adapter(smiles, mode="bioactivity")

if result.success:
    data = result.data
    print(f"Success: {data}")
else:
    print(f"Error: {result.error}")
    # Handle gracefully
```

## Example Scripts

Each adapter includes comprehensive examples:

- `adapters/zinc_fragments/example_usage.py` - Fragment search examples
- `adapters/surechembl/example_usage.py` - Patent search examples
- `adapters/pubchem/example_enhanced.py` - PubChem features
- `adapters/chembl/example_enhanced.py` - Bioactivity queries
- `adapters/INTEGRATION_EXAMPLES.py` - Retrosynthesis workflows

Run examples:
```bash
python adapters/zinc_fragments/example_usage.py
python adapters/INTEGRATION_EXAMPLES.py
```

## Use Cases

### Drug Discovery
- **Lead identification** - Find bioactive compounds in ChEMBL
- **Analog searching** - Similarity search in PubChem/ChEMBL
- **Target validation** - ChEMBL target queries
- **Drug repurposing** - Mechanism of action analysis

### Synthesis Planning
- **Building block sourcing** - ZINC/PubChem vendor lookup
- **Route optimization** - Purchasability-guided retrosynthesis
- **IP analysis** - SureChEMBL freedom-to-operate
- **Fragment-based design** - ZINC fragment library

### Competitive Intelligence
- **Patent landscape** - SureChEMBL competitor analysis
- **Clinical pipeline** - ChEMBL drug status
- **Market analysis** - Patent applicant trends
- **Compound extraction** - Mine competitor patents

## API Endpoints

### ZINC15
- Base: `http://zinc15.docking.org`
- Search: `/substances/search/`
- Details: `/substances/{zinc_id}.json`

### SureChEMBL
- Base: `https://www.surechembl.org/api`
- Structure search: `/search/similarity`
- Patent lookup: `/patent/{patent_id}`

### PubChem
- Base: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`
- Properties: `/compound/smiles/{smiles}/property/{props}/JSON`
- Similarity: `/compound/fastsimilarity_2d/smiles/{smiles}/cids/JSON`
- Vendors: `/compound/cid/{cid}/synonyms/JSON`

### ChEMBL
- Base: `https://www.ebi.ac.uk/chembl/api/data`
- Molecule: `/molecule.json` or `/molecule/{chembl_id}.json`
- Activity: `/activity.json`
- Target: `/target/{target_id}.json`
- Drug: `/drug.json`
- Mechanism: `/mechanism.json`

## Best Practices

1. **Rate limiting** - Respect API limits, use delays
2. **Caching** - Enable caching for repeated queries
3. **Error handling** - Always check `result.success`
4. **Batch operations** - Use `asyncio.gather()` for parallel queries
5. **Filtering** - Apply filters to reduce result sizes
6. **Validation** - Validate building blocks before synthesis
7. **IP checks** - Check patents early in planning

## Dependencies

```python
# Required packages
aiohttp>=3.8.0
asyncio
logging
```

## License

These adapters are part of PharmForge and follow the same license.

## Support

For issues or questions:
- Check example scripts for usage patterns
- Review API documentation for each database
- Test with small queries before scaling up

## Future Enhancements

- [ ] eMolecules integration for additional vendors
- [ ] ChemSpider for additional compound data
- [ ] Reaxys integration for reaction data
- [ ] SciFinder integration for literature
- [ ] Automated route scoring and ranking
- [ ] ML-based purchasability prediction
