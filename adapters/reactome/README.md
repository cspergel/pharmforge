# Reactome Adapter

Comprehensive adapter for pathway analysis and biological process enrichment using the Reactome Pathway Database.

## Features

- **Pathway Search**: Search pathways by name, keyword, or biological process
- **Gene/Protein Analysis**: Map genes/proteins to biological pathways
- **Enrichment Analysis**: Identify overrepresented pathways in gene sets
- **Pathway Details**: Get comprehensive pathway information including:
  - Hierarchical structure
  - Participating molecules
  - Literature references
  - Disease associations
- **Multi-Species Support**: Human, mouse, rat, and other model organisms

## API Documentation

- Base URL: https://reactome.org/ContentService
- Documentation: https://reactome.org/dev/content-service
- Database: Curated by EMBL-EBI and OICR

## Installation

No additional dependencies beyond the base requirements.

## Usage

### Search Pathways by Name

```python
from adapters.reactome.adapter import ReactomeAdapter

adapter = ReactomeAdapter()

# Search for pathways
result = await adapter.execute("apoptosis")

if result.success:
    pathways = result.data['pathways']
    print(f"Found {len(pathways)} pathways")

    for pathway in pathways[:5]:
        print(f"{pathway['pathway_id']}: {pathway['name']}")
```

### Gene Set Enrichment Analysis

```python
# Analyze a list of genes
genes = ["TP53", "BRCA1", "PTEN", "RB1", "MYC"]

result = await adapter.execute(genes, enrichment=True)

if result.success:
    enriched = result.data['enriched_pathways']
    print(f"Found {len(enriched)} enriched pathways")

    for pathway in enriched[:10]:
        print(f"\n{pathway['name']}")
        print(f"  ID: {pathway['pathway_id']}")
        print(f"  Sub-pathways: {pathway.get('sub_pathways', 0)}")
```

### Get Detailed Pathway Information

```python
# Get details for a specific pathway
result = await adapter.execute(
    {"pathway_id": "R-HSA-69278"},  # Cell Cycle
    include_hierarchy=True,
    include_participants=True
)

if result.success:
    pathway = result.data['pathway']

    print(f"Name: {pathway['name']}")
    print(f"Description: {pathway.get('description', 'N/A')}")

    if pathway.get('references'):
        print("\nReferences:")
        for ref in pathway['references'][:3]:
            print(f"  PMID {ref['pubmed_id']}: {ref['title']}")

    if pathway.get('participant_types'):
        print("\nParticipants:")
        for ptype, count in pathway['participant_types'].items():
            print(f"  {ptype}: {count}")
```

## Input Parameters

### Simple String Input
- **query** (str): Search term for pathways or single gene/protein

### List Input
- **genes** (list): List of gene symbols or protein identifiers

### Dictionary Input
- **pathway_id** (str): Reactome pathway stable identifier (R-HSA-xxxxx)
- **genes** (list): List of gene symbols
- **proteins** (list): List of protein identifiers
- **species** (str): Species name (default: "Homo sapiens")

### Keyword Arguments
- **species** (str): Organism (default: "Homo sapiens")
- **enrichment** (bool): Perform enrichment analysis for gene lists (default: True)
- **include_hierarchy** (bool): Include pathway hierarchy (default: True)
- **include_participants** (bool): Include participating molecules (default: True)

## Output Format

### Search Results
```python
{
    "pathways": [
        {
            "pathway_id": "R-HSA-109581",
            "name": "Apoptosis",
            "species": "Homo sapiens",
            "type": "Pathway",
            "matched_gene": "TP53"  # If searching by gene
        }
    ],
    "total_found": 15
}
```

### Enrichment Results
```python
{
    "enriched_pathways": [
        {
            "pathway_id": "R-HSA-69278",
            "name": "Cell Cycle",
            "species": "Homo sapiens",
            "type": "Pathway",
            "sub_pathways": 15,
            "sub_pathway_list": [
                {
                    "id": "R-HSA-69242",
                    "name": "S Phase",
                    "type": "Pathway"
                }
            ],
            "num_participants": 500,
            "participant_types": {
                "Protein": 300,
                "Complex": 150,
                "Small Molecule": 50
            }
        }
    ],
    "input_genes": ["TP53", "BRCA1", "PTEN"],
    "species": "Homo sapiens",
    "total_pathways": 25
}
```

### Detailed Pathway
```python
{
    "pathway": {
        "pathway_id": "R-HSA-69278",
        "name": "Cell Cycle",
        "species": "Homo sapiens",
        "type": "Pathway",
        "description": "The replication of the genome and the subsequent segregation...",
        "disease_related": false,
        "inferred": false,
        "references": [
            {
                "pubmed_id": "12345678",
                "title": "Cell cycle regulation..."
            }
        ],
        "sub_pathways": 15,
        "sub_pathway_list": [...],
        "num_participants": 500,
        "participant_types": {...}
    }
}
```

## Reactome Pathway Hierarchy

Reactome organizes pathways hierarchically:

```
Top-level pathway
├── Sub-pathway 1
│   ├── Sub-sub-pathway 1.1
│   └── Sub-sub-pathway 1.2
└── Sub-pathway 2
```

### Example Hierarchy

```python
result = await adapter.execute(
    {"pathway_id": "R-HSA-1640170"},  # Cell Cycle
    include_hierarchy=True
)

if result.success:
    pathway = result.data['pathway']

    print(f"{pathway['name']} contains {pathway['sub_pathways']} sub-pathways:")
    for sub in pathway['sub_pathway_list']:
        print(f"  - {sub['name']} ({sub['id']})")
```

## Example Workflows

### 1. Cancer Pathway Analysis

```python
adapter = ReactomeAdapter()

# Oncogenes and tumor suppressors
cancer_genes = ["TP53", "BRCA1", "KRAS", "MYC", "PTEN", "RB1"]

result = await adapter.execute(cancer_genes, enrichment=False)

if result.success:
    pathways = result.data['pathways']

    # Find cancer-related pathways
    cancer_pathways = [
        p for p in pathways
        if any(term in p['name'].lower() for term in ['cancer', 'apoptosis', 'cell cycle', 'dna repair'])
    ]

    print(f"Found {len(cancer_pathways)} cancer-related pathways:")
    for pathway in cancer_pathways[:10]:
        print(f"  - {pathway['name']}")
```

### 2. Drug Target Pathway Mapping

```python
# Map drug targets to pathways
drug_targets = ["EGFR", "ERBB2", "VEGFA", "KDR"]

result = await adapter.execute(drug_targets, enrichment=False)

if result.success:
    pathways = result.data['pathways']

    # Group by target
    target_pathways = {}
    for pathway in pathways:
        gene = pathway.get('matched_gene')
        if gene:
            if gene not in target_pathways:
                target_pathways[gene] = []
            target_pathways[gene].append(pathway['name'])

    for target, pathway_list in target_pathways.items():
        print(f"\n{target} participates in {len(pathway_list)} pathways:")
        for p in pathway_list[:5]:
            print(f"  - {p}")
```

### 3. Pathway Comparison Between Gene Sets

```python
# Compare two gene sets
geneset1 = ["TP53", "BRCA1", "PTEN"]
geneset2 = ["KRAS", "MYC", "EGFR"]

result1 = await adapter.execute(geneset1, enrichment=False)
result2 = await adapter.execute(geneset2, enrichment=False)

if result1.success and result2.success:
    pathways1 = set(p['pathway_id'] for p in result1.data['pathways'])
    pathways2 = set(p['pathway_id'] for p in result2.data['pathways'])

    common = pathways1 & pathways2
    unique1 = pathways1 - pathways2
    unique2 = pathways2 - pathways1

    print(f"Common pathways: {len(common)}")
    print(f"Unique to set 1: {len(unique1)}")
    print(f"Unique to set 2: {len(unique2)}")
```

### 4. Immune System Pathway Analysis

```python
# Analyze immune system pathways
result = await adapter.execute("immune system")

if result.success:
    pathways = result.data['pathways']

    # Categorize immune pathways
    categories = {
        'Innate': [],
        'Adaptive': [],
        'Cytokine': [],
        'Other': []
    }

    for pathway in pathways:
        name = pathway['name'].lower()
        if 'innate' in name:
            categories['Innate'].append(pathway['name'])
        elif 'adaptive' in name or 'b cell' in name or 't cell' in name:
            categories['Adaptive'].append(pathway['name'])
        elif 'cytokine' in name or 'interleukin' in name:
            categories['Cytokine'].append(pathway['name'])
        else:
            categories['Other'].append(pathway['name'])

    for category, pathway_list in categories.items():
        if pathway_list:
            print(f"\n{category} Immunity: {len(pathway_list)} pathways")
            for p in pathway_list[:3]:
                print(f"  - {p}")
```

### 5. Disease Pathway Discovery

```python
# Search for disease-associated pathways
diseases = ["diabetes", "alzheimer", "parkinson", "cancer"]

for disease in diseases:
    result = await adapter.execute(disease)

    if result.success:
        pathways = result.data['pathways']
        print(f"\n{disease.title()}: {len(pathways)} pathways")

        for pathway in pathways[:3]:
            print(f"  - {pathway['name']}")
```

## Multi-Species Support

Reactome supports multiple organisms:

```python
# Analyze mouse genes
mouse_genes = ["Tp53", "Brca1", "Pten"]

result = await adapter.execute(
    mouse_genes,
    species="Mus musculus",
    enrichment=False
)

# Supported species
species_list = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Caenorhabditis elegans",
    "Drosophila melanogaster",
    "Saccharomyces cerevisiae",
    "Danio rerio",
    # ... and many more
]
```

## Understanding Reactome Data

### Pathway Types
- **Pathway**: Biological process with multiple steps
- **TopLevelPathway**: High-level category
- **Reaction**: Single molecular event

### Participant Types
- **Protein**: Individual protein
- **Complex**: Protein complex
- **EntitySet**: Group of similar entities
- **SimpleEntity**: Small molecule
- **Polymer**: DNA, RNA sequences

### Evidence Codes
- **IDA**: Inferred from Direct Assay
- **TAS**: Traceable Author Statement
- **IEA**: Inferred from Electronic Annotation

## Integration with Other Adapters

### With DrugCentral
```python
# Get drug targets, then analyze pathways
drug_result = await drugcentral.execute("imatinib", include_targets=True)
targets = [t['gene'] for t in drug_result.data['drugs'][0]['targets']]

pathway_result = await reactome.execute(targets, enrichment=False)
```

### With Clinical Trials
```python
# Find pathways, then search for trials targeting those pathways
pathway_result = await reactome.execute("cancer")
# Use pathway genes to find relevant trials
```

## Visualization

Reactome provides pathway diagrams:

```python
# Get diagram URL
pathway_id = "R-HSA-69278"
diagram_url = f"https://reactome.org/PathwayBrowser/#{pathway_id}"
print(f"View diagram: {diagram_url}")

# Export diagram
export_url = f"https://reactome.org/ContentService/exporter/diagram/{pathway_id}.png"
```

## Error Handling

```python
result = await adapter.execute("invalid_input")

if not result.success:
    print(f"Error: {result.error}")
elif result.data.get('total_found', 0) == 0:
    print("No pathways found")
else:
    # Process data
    pass
```

## Rate Limiting

```python
adapter = ReactomeAdapter()
adapter.config['rate_limit_delay'] = 0.5  # 500ms between requests
```

## Best Practices

1. **Use Enrichment**: For gene lists >5, use enrichment analysis
2. **Species Specific**: Always specify correct species
3. **Hierarchy**: Explore pathway hierarchy for context
4. **References**: Check literature references for validation
5. **Multiple Queries**: Compare results from different queries

## Notes

- Reactome is manually curated by experts
- Updated quarterly with latest research
- Focus on human pathways (most comprehensive)
- Cross-references to other databases (UniProt, ChEBI, PubMed)

## Limitations

- Not all genes have pathway annotations
- Some pathways are incomplete
- Inference across species may have gaps
- API rate limits apply

## Future Enhancements

- [ ] Statistical enrichment with p-values
- [ ] Pathway diagram integration
- [ ] Cross-species pathway comparison
- [ ] Drug target pathway prioritization
- [ ] Disease pathway networks
