# KEGG Adapter for PharmForge

## Overview

The KEGG (Kyoto Encyclopedia of Genes and Genomes) adapter provides access to comprehensive pathway, gene, compound, disease, and drug information using the FREE KEGG REST API.

## Features

- **Pathway Queries**: Search and retrieve metabolic and signaling pathway information
- **Gene Queries**: Access gene information across organisms with pathway associations
- **Compound Queries**: Query chemical compounds by ID, name, or formula
- **Disease Queries**: Get disease information with associated genes and pathways
- **Drug Queries**: Access pharmaceutical compound data with target information
- **Auto-Detection**: Automatically detects query type based on input format
- **Cross-References**: Links between databases (pathways, genes, compounds, diseases)

## Installation

The adapter is already integrated into PharmForge. No additional dependencies required.

```python
from adapters.kegg.adapter import KEGGAdapter

adapter = KEGGAdapter()
```

## API Endpoints Used

The adapter uses the FREE KEGG REST API (https://rest.kegg.jp):

- `/info/<database>` - Database information
- `/list/<database>` - List entries in database
- `/find/<database>/<query>` - Search entries
- `/get/<entry_ids>` - Get full entry details
- `/link/<target_db>/<source>` - Cross-reference between databases
- `/conv/<target_db>/<source_db>` - Convert between ID systems

## Usage Examples

### 1. Query Pathways

#### Get Specific Pathway
```python
result = await adapter.execute({
    "operation": "query_pathway",
    "pathway_id": "hsa00010"  # Glycolysis
})

if result.success:
    data = result.data
    print(f"Pathway: {data['data']['NAME']}")
    print(f"Description: {data['data']['DESCRIPTION']}")
```

#### Search Pathways by Keyword
```python
result = await adapter.execute({
    "operation": "query_pathway",
    "search_term": "apoptosis"
})

for pathway in result.data['results']:
    print(f"{pathway['id']}: {pathway['description']}")
```

#### List Organism-Specific Pathways
```python
result = await adapter.execute({
    "operation": "query_pathway",
    "organism": "hsa"  # Human
})

print(f"Found {result.data['count']} human pathways")
```

### 2. Query Genes

#### Get Specific Gene
```python
result = await adapter.execute({
    "operation": "query_gene",
    "gene_id": "hsa:7157"  # TP53
})

if result.success:
    gene = result.data
    print(f"Gene: {gene['data']['NAME']}")
    print(f"Pathways: {len(gene['pathways'])}")
```

#### Search Genes by Keyword
```python
result = await adapter.execute({
    "operation": "query_gene",
    "search_term": "insulin receptor"
})

for gene in result.data['results'][:5]:
    print(f"{gene['id']}: {gene['description']}")
```

#### Auto-Detection (Gene ID)
```python
result = await adapter.execute("hsa:7157")  # Automatically detects as gene query
```

### 3. Query Compounds

#### Get Specific Compound
```python
result = await adapter.execute({
    "operation": "query_compound",
    "compound_id": "C00031"  # Glucose
})

if result.success:
    compound = result.data
    print(f"Compound: {compound['data']['NAME']}")
    print(f"Formula: {compound['data']['FORMULA']}")
    print(f"Pathways: {len(compound['pathways'])}")
```

#### Search Compounds by Name
```python
result = await adapter.execute({
    "operation": "query_compound",
    "search_term": "aspirin"
})
```

#### Search by Chemical Formula
```python
result = await adapter.execute({
    "operation": "query_compound",
    "formula": "C6H12O6"
})

for compound in result.data['results']:
    print(f"{compound['id']}: {compound['description']}")
```

#### Auto-Detection (Compound ID)
```python
result = await adapter.execute("C00031")  # Automatically detects as compound query
```

### 4. Query Diseases

#### Get Specific Disease
```python
result = await adapter.execute({
    "operation": "query_disease",
    "disease_id": "H00056"  # Alzheimer's
})

if result.success:
    disease = result.data
    print(f"Disease: {disease['data']['NAME']}")
    print(f"Associated genes: {len(disease['genes'])}")
    print(f"Associated pathways: {len(disease['pathways'])}")
```

#### Search Diseases by Keyword
```python
result = await adapter.execute({
    "operation": "query_disease",
    "search_term": "cancer"
})

for disease in result.data['results'][:5]:
    print(f"{disease['id']}: {disease['description']}")
```

#### Auto-Detection (Disease ID)
```python
result = await adapter.execute("H00056")  # Automatically detects as disease query
```

### 5. Query Drugs

#### Get Specific Drug
```python
result = await adapter.execute({
    "operation": "query_drug",
    "drug_id": "D00109"  # Aspirin
})

if result.success:
    drug = result.data
    print(f"Drug: {drug['data']['NAME']}")
    print(f"Targets: {len(drug['targets'])}")
```

#### Search Drugs by Keyword
```python
result = await adapter.execute({
    "operation": "query_drug",
    "search_term": "metformin"
})
```

#### Auto-Detection (Drug ID)
```python
result = await adapter.execute("D00109")  # Automatically detects as drug query
```

## ID Format Reference

The adapter auto-detects query types based on KEGG ID formats:

- **Gene IDs**: `<org>:<id>` (e.g., `hsa:7157`)
- **Compound IDs**: `C#####` (e.g., `C00031`)
- **Drug IDs**: `D#####` (e.g., `D00109`)
- **Disease IDs**: `H#####` (e.g., `H00056`)
- **Pathway IDs**: `<org>##### or map#####` (e.g., `hsa00010`, `map00010`)

## Organism Codes

Common organism codes:
- `hsa` - Homo sapiens (human)
- `mmu` - Mus musculus (mouse)
- `rno` - Rattus norvegicus (rat)
- `dme` - Drosophila melanogaster (fruit fly)
- `cel` - Caenorhabditis elegans (worm)
- `sce` - Saccharomyces cerevisiae (yeast)
- `eco` - Escherichia coli K-12 MG1655

Full list: https://www.genome.jp/kegg/catalog/org_list.html

## Advanced Integration Examples

### Trace Disease → Genes → Pathways
```python
# Search for disease
result = await adapter.execute({
    "operation": "query_disease",
    "search_term": "diabetes"
})

disease_id = result.data['results'][0]['id']

# Get disease details with genes and pathways
result = await adapter.execute({
    "operation": "query_disease",
    "disease_id": disease_id
})

print(f"Disease: {result.data['data']['NAME']}")
print(f"Associated genes: {len(result.data['genes'])}")
print(f"Associated pathways: {len(result.data['pathways'])}")

for pathway in result.data['pathways'][:5]:
    print(f"  - {pathway['target']}")
```

### Compound Pathway Analysis
```python
# Query compound
result = await adapter.execute("C00031")  # Glucose

if result.success:
    print(f"Compound: {result.data['data']['NAME']}")
    print(f"Found in {len(result.data['pathways'])} pathways")

    # Get pathway details for top pathways
    for pathway_link in result.data['pathways'][:3]:
        pathway_id = pathway_link['target'].replace('path:', '')

        pathway_result = await adapter.execute({
            "operation": "query_pathway",
            "pathway_id": pathway_id
        })

        if pathway_result.success:
            print(f"\nPathway: {pathway_result.data['data']['NAME']}")
```

### Multi-Gene Pathway Enrichment
```python
# Search multiple genes
genes = ["TP53", "BRCA1", "PTEN"]
gene_pathways = {}

for gene in genes:
    result = await adapter.execute({
        "operation": "query_gene",
        "search_term": gene
    })

    if result.success and result.data['results']:
        gene_id = result.data['results'][0]['id']

        # Get gene details with pathways
        gene_result = await adapter.execute({
            "operation": "query_gene",
            "gene_id": gene_id
        })

        if gene_result.success:
            gene_pathways[gene] = [p['target'] for p in gene_result.data['pathways']]

# Find common pathways
common_pathways = set(gene_pathways[genes[0]])
for gene in genes[1:]:
    common_pathways = common_pathways.intersection(gene_pathways[gene])

print(f"Common pathways for {', '.join(genes)}: {len(common_pathways)}")
```

## Response Format

All queries return an `AdapterResult` object:

```python
@dataclass
class AdapterResult:
    success: bool              # Query success status
    data: Any                  # Query results
    error: Optional[str]       # Error message if failed
    cache_hit: bool           # Whether result was cached
    metadata: Dict[str, Any]  # Additional metadata
```

### Successful Response Example
```python
result = await adapter.execute("C00031")

# result.success = True
# result.data = {
#     "compound_id": "C00031",
#     "data": {
#         "NAME": "D-Glucose; Grape sugar; Dextrose...",
#         "FORMULA": "C6H12O6",
#         "MOL_WEIGHT": "180.0634",
#         ...
#     },
#     "pathways": [
#         {"source": "cpd:C00031", "target": "path:map00010"},
#         ...
#     ],
#     "raw": "..."  # Raw KEGG flat file
# }
```

## Error Handling

```python
result = await adapter.execute({"operation": "query_compound", "compound_id": "INVALID"})

if not result.success:
    print(f"Error: {result.error}")
    # Output: "Error: Compound INVALID not found"
```

## Rate Limiting

The adapter implements automatic rate limiting (0.3 seconds between requests) to be respectful to the free KEGG API. This is configurable:

```python
adapter = KEGGAdapter()
adapter.config['rate_limit_delay'] = 0.5  # Increase delay to 0.5 seconds
```

## Limitations

1. **Entry Limits**: Some operations are limited to 10 entries per request (KEGG API limitation)
2. **No Authentication**: Uses FREE public API (no API key required)
3. **Rate Limiting**: Recommended delay between requests to avoid overwhelming the service
4. **Data Format**: Returns KEGG flat file format which is parsed into dictionaries
5. **Link Operations**: Some cross-database links may not be available for all entity types

## API Limitations and Notes

- Maximum 10 entries per GET request
- Some /link operations return 400 for incompatible database pairs (expected behavior)
- Search results are automatically limited (e.g., 50-100 entries) to prevent overwhelming responses
- Images and KGML files can be retrieved by using `option` parameter in `_get_entry` method

## Troubleshooting

### Common Issues

1. **"Resource not found" error**
   - Verify the ID format is correct (e.g., `C00031` not `c00031`)
   - Check that the ID exists in KEGG database

2. **Timeout errors**
   - Increase timeout in config: `adapter.config['timeout'] = 60`
   - Check network connectivity

3. **Empty results**
   - Try broadening search terms
   - Verify organism code is correct for gene/pathway queries

4. **400 errors for link operations**
   - Not all entity types can be linked (e.g., drugs → genes may not work)
   - This is expected KEGG API behavior

## Testing

Run the comprehensive test suite:

```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"
python adapters/kegg/test_adapter.py
```

Tests cover:
- Pathway queries (by ID, search, organism)
- Gene queries (by ID, search, organism)
- Compound queries (by ID, name, formula)
- Disease queries (by ID, search)
- Drug queries (by ID, search)
- Auto-detection of ID formats
- Integration scenarios

## References

- **KEGG REST API**: https://www.kegg.jp/kegg/rest/keggapi.html
- **KEGG Database**: https://www.kegg.jp/kegg/
- **Organism Codes**: https://www.genome.jp/kegg/catalog/org_list.html
- **KEGG Identifiers**: https://www.kegg.jp/kegg/kegg1.html

## Version

Current version: 1.0.0

## License

Part of PharmForge project. Uses the FREE KEGG REST API.
