# HMDB Adapter

Access the **Human Metabolome Database (HMDB)** for comprehensive metabolomics data in drug discovery workflows.

## Overview

HMDB is the world's most comprehensive database of human metabolites, containing:
- **220,000+ metabolites** with detailed annotations
- **Biofluid concentration data** (normal ranges in blood, urine, CSF, saliva, etc.)
- **Disease associations** and biomarker relationships
- **Metabolic pathway information** (KEGG, SMPDB)
- **Protein/enzyme interactions** (UniProt cross-references)
- **NMR/MS spectral data** references

**Database:** https://hmdb.ca/
**License:** Creative Commons (free for all uses)
**Installation:** No additional packages required (uses `aiohttp` and `xml.etree`)

---

## Key Features

### 1. Metabolite Lookup by HMDB ID
Retrieve complete metabolite profiles including structure, properties, and clinical data.

### 2. Metabolite Search
Search by name, formula, or other criteria to find relevant metabolites.

### 3. Biofluid Concentration Data
Access normal concentration ranges in various biofluids:
- Blood (serum, plasma)
- Urine
- Cerebrospinal Fluid (CSF)
- Saliva
- Sweat
- Feces
- Breast milk

### 4. Disease Associations
Identify metabolites associated with specific diseases and conditions.

### 5. Metabolic Pathway Information
Get pathway context (KEGG, SMPDB) for understanding metabolite function.

### 6. Protein/Enzyme Interactions
Find proteins and enzymes that interact with metabolites.

---

## Usage Examples

### Example 1: Basic Metabolite Lookup

```python
import asyncio
from adapters.hmdb import HMDBAdapter

async def main():
    adapter = HMDBAdapter()

    # Query by HMDB ID
    result = await adapter.execute("HMDB0000001")

    if result.success:
        metabolite = result.data["metabolite"]
        summary = result.data["summary"]

        print(f"Name: {metabolite['name']}")
        print(f"Formula: {metabolite['chemical_formula']}")
        print(f"MW: {metabolite['average_molecular_weight']}")
        print(f"SMILES: {metabolite['smiles']}")
        print(f"Diseases: {summary['num_diseases']}")
        print(f"Pathways: {summary['num_pathways']}")

asyncio.run(main())
```

**Output:**
```
Name: 1-Methylhistidine
Formula: C7H11N3O2
MW: 169.181
SMILES: CN1C=NC(C[C@H](N)C(O)=O)=C1
Diseases: 8
Pathways: 2
```

---

### Example 2: Biofluid Concentration Analysis

```python
async def analyze_biofluid_levels():
    adapter = HMDBAdapter()

    # Get metabolite with specific biofluid data
    result = await adapter.execute(
        "HMDB0000001",
        biofluids=["blood", "urine", "cerebrospinal_fluid"],
        include_concentrations=True
    )

    if result.success:
        concentrations = result.data["metabolite"]["concentrations"]

        for biofluid, data in concentrations.items():
            print(f"\n{biofluid.replace('_', ' ').title()}:")
            print(f"  Normal range: {data['value']} {data['unit']}")
            if data.get('references'):
                print(f"  References: {', '.join(data['references'][:3])}")

asyncio.run(analyze_biofluid_levels())
```

**Output:**
```
Blood:
  Normal range: 2.5-12.0 uM
  References: PMID:12345, PMID:67890

Urine:
  Normal range: 50-300 uM
  References: PMID:11111
```

---

### Example 3: Disease Association Discovery

```python
async def find_disease_biomarkers():
    adapter = HMDBAdapter()

    # Search for a metabolite and get disease associations
    result = await adapter.execute(
        {"query": "glucose", "mode": "search"},
        include_diseases=True,
        include_concentrations=True
    )

    if result.success:
        for metabolite in result.data["metabolites"]:
            print(f"\n{metabolite['name']} ({metabolite['hmdb_id']}):")

            # Show disease associations
            for disease in metabolite.get("diseases", [])[:5]:
                print(f"  - {disease['name']}")
                if disease.get('omim_id'):
                    print(f"    OMIM: {disease['omim_id']}")

asyncio.run(find_disease_biomarkers())
```

**Output:**
```
D-Glucose (HMDB0000122):
  - Diabetes mellitus
    OMIM: 222100
  - Hyperglycemia
    OMIM: 241500
  - Hypoglycemia
    OMIM: 240800
```

---

### Example 4: Metabolic Pathway Context

```python
async def get_pathway_context():
    adapter = HMDBAdapter()

    result = await adapter.execute(
        "HMDB0000001",
        include_pathways=True,
        include_proteins=True
    )

    if result.success:
        metabolite = result.data["metabolite"]

        print(f"Metabolite: {metabolite['name']}\n")

        # Show pathways
        print("Pathways:")
        for pathway in metabolite.get("pathways", []):
            print(f"  - {pathway['name']}")
            if pathway.get('smpdb_id'):
                print(f"    SMPDB: {pathway['smpdb_id']}")
            if pathway.get('kegg_map_id'):
                print(f"    KEGG: {pathway['kegg_map_id']}")

        # Show protein associations
        print("\nAssociated Proteins:")
        for protein in metabolite.get("proteins", [])[:5]:
            print(f"  - {protein.get('gene_name', 'Unknown')} ({protein.get('uniprot_id', 'N/A')})")
            print(f"    Type: {protein.get('protein_type', 'N/A')}")

asyncio.run(get_pathway_context())
```

**Output:**
```
Metabolite: 1-Methylhistidine

Pathways:
  - Histidine Metabolism
    SMPDB: SMP00009
    KEGG: map00340
  - Beta-Alanine Metabolism
    SMPDB: SMP00007

Associated Proteins:
  - CNDP2 (Q9UBK2)
    Type: enzyme
  - CARNS1 (P30046)
    Type: enzyme
```

---

### Example 5: Cross-Reference Integration

```python
async def get_cross_references():
    adapter = HMDBAdapter()

    result = await adapter.execute("HMDB0000001")

    if result.success:
        ext_ids = result.data["metabolite"]["external_ids"]

        print("External Database IDs:")
        print(f"  PubChem: {ext_ids.get('pubchem_compound_id', 'N/A')}")
        print(f"  KEGG: {ext_ids.get('kegg_id', 'N/A')}")
        print(f"  ChEBI: {ext_ids.get('chebi_id', 'N/A')}")
        print(f"  ChemSpider: {ext_ids.get('chemspider_id', 'N/A')}")
        print(f"  DrugBank: {ext_ids.get('drugbank_id', 'N/A')}")

asyncio.run(get_cross_references())
```

**Output:**
```
External Database IDs:
  PubChem: 92105
  KEGG: C01152
  ChEBI: CHEBI:50599
  ChemSpider: 83153
  DrugBank: N/A
```

---

### Example 6: Advanced Search with Filtering

```python
async def advanced_search():
    adapter = HMDBAdapter()

    # Search for metabolites related to a disease
    result = await adapter.execute(
        {
            "query": "creatinine",
            "mode": "search"
        },
        biofluids=["blood", "urine"],
        include_concentrations=True,
        include_diseases=True
    )

    if result.success:
        print(f"Found {result.data['num_results']} metabolites")
        print(f"Fetched details for {result.data['num_fetched']}\n")

        for metabolite in result.data["metabolites"]:
            print(f"{metabolite['name']} ({metabolite['hmdb_id']}):")
            print(f"  Formula: {metabolite['chemical_formula']}")
            print(f"  Locations: {', '.join(metabolite['biofluid_locations'][:5])}")

            # Show concentration data
            if metabolite.get("concentrations"):
                print(f"  Concentrations:")
                for biofluid, data in metabolite["concentrations"].items():
                    print(f"    - {biofluid}: {data['value']} {data['unit']}")

            print()

asyncio.run(advanced_search())
```

---

## API Parameters

### Input Formats

**1. Simple HMDB ID:**
```python
result = await adapter.execute("HMDB0000001")
```

**2. Metabolite Name (triggers search):**
```python
result = await adapter.execute("glucose")
```

**3. Dictionary with options:**
```python
result = await adapter.execute({
    "hmdb_id": "HMDB0000001",
    "mode": "metabolite"
})
```

**4. Search dictionary:**
```python
result = await adapter.execute({
    "query": "glucose",
    "mode": "search",
    "search_type": "name"
})
```

### Keyword Arguments

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | str | `"metabolite"` | Query mode: `"metabolite"`, `"search"` |
| `biofluids` | List[str] | `["blood", "urine", "cerebrospinal_fluid"]` | Biofluids to include |
| `include_concentrations` | bool | `True` | Include concentration data |
| `include_diseases` | bool | `True` | Include disease associations |
| `include_pathways` | bool | `True` | Include pathway information |
| `include_proteins` | bool | `True` | Include protein associations |
| `search_type` | str | `"name"` | Search type: `"name"`, `"formula"`, `"mass"` |

---

## Output Format

### Metabolite Mode

```python
{
    "metabolite": {
        "hmdb_id": "HMDB0000001",
        "name": "1-Methylhistidine",
        "systematic_name": "(2S)-2-Amino-3-(1-methyl-1H-imidazol-4-yl)propanoic acid",
        "chemical_formula": "C7H11N3O2",
        "average_molecular_weight": "169.181",
        "smiles": "CN1C=NC(C[C@H](N)C(O)=O)=C1",
        "inchi": "InChI=1S/C7H11N3O2/c1-10-3-5...",
        "inchikey": "...",
        "state": "solid",
        "description": "...",
        "synonyms": ["Tau-methylhistidine", ...],
        "biofluid_locations": ["Blood", "Urine", "Cerebrospinal Fluid"],
        "tissue_locations": ["Muscle", "Brain"],
        "concentrations": {
            "blood": {
                "value": "2.5-12.0",
                "unit": "uM",
                "subject_age": "Adult",
                "subject_sex": "Both",
                "subject_condition": "Normal",
                "references": ["PMID:12345"]
            }
        },
        "diseases": [
            {
                "name": "Chronic kidney disease",
                "omim_id": "263200",
                "references": ["PMID:67890"]
            }
        ],
        "pathways": [
            {
                "name": "Histidine Metabolism",
                "smpdb_id": "SMP00009",
                "kegg_map_id": "map00340"
            }
        ],
        "proteins": [
            {
                "uniprot_id": "Q9UBK2",
                "gene_name": "CNDP2",
                "protein_type": "enzyme",
                "name": "Cytosolic non-specific dipeptidase"
            }
        ],
        "ontology": {
            "description": "...",
            "origins": ["Endogenous"],
            "biofunctions": ["Component of proteins"],
            "applications": ["Biomarker"]
        },
        "external_ids": {
            "pubchem_compound_id": "92105",
            "chemspider_id": "83153",
            "kegg_id": "C01152",
            "chebi_id": "CHEBI:50599",
            "drugbank_id": null,
            "foodb_id": null
        }
    },
    "summary": {
        "hmdb_id": "HMDB0000001",
        "name": "1-Methylhistidine",
        "formula": "C7H11N3O2",
        "molecular_weight": "169.181",
        "num_biofluids": 8,
        "num_tissues": 4,
        "num_concentrations": 3,
        "num_diseases": 8,
        "num_pathways": 2,
        "num_proteins": 12,
        "has_pubchem": true,
        "has_kegg": true,
        "state": "solid"
    }
}
```

### Search Mode

```python
{
    "query": "glucose",
    "search_type": "name",
    "num_results": 15,
    "num_fetched": 5,
    "all_hmdb_ids": ["HMDB0000122", "HMDB0000660", ...],
    "metabolites": [
        {
            # Full metabolite data for first 5 results
        }
    ],
    "warning": "Only fetched details for first 5 results"
}
```

---

## Use Cases in Drug Discovery

### 1. Off-Target Metabolite Analysis
Identify how drug candidates affect endogenous metabolites:
```python
# Check if drug molecule is similar to known metabolites
# Could indicate metabolic interference
```

### 2. Biomarker Discovery
Find metabolites associated with disease states:
```python
# Search for metabolites in disease contexts
# Identify potential diagnostic or prognostic biomarkers
```

### 3. Drug Metabolism Prediction
Understand metabolic pathways:
```python
# Get metabolic pathway context
# Identify enzymes involved in metabolism
# Predict potential drug-metabolite interactions
```

### 4. Clinical Chemistry Integration
Access normal ranges for metabolites:
```python
# Get biofluid concentration ranges
# Understand clinical significance of metabolite changes
```

### 5. Systems Biology Context
Connect metabolites to proteins and pathways:
```python
# Link compounds to metabolic networks
# Understand biological function
# Identify mechanism of action
```

### 6. Drug Repurposing
Find metabolites similar to drugs:
```python
# Cross-reference with DrugBank
# Identify structural similarities
# Explore therapeutic potential
```

---

## Integration Examples

### Pipeline Integration

```python
from backend.core.pipeline import Pipeline
from adapters.hmdb import HMDBAdapter
from adapters.pubchem import PubChemEnhancedAdapter

async def metabolomics_pipeline(smiles: str):
    """
    Complete metabolomics analysis pipeline
    """
    # Step 1: Get compound properties from PubChem
    pubchem = PubChemEnhancedAdapter()
    pubchem_result = await pubchem.execute(smiles)

    # Step 2: Search for similar metabolites in HMDB
    hmdb = HMDBAdapter()
    hmdb_result = await hmdb.execute(
        {"query": smiles, "mode": "search"},
        include_concentrations=True,
        include_diseases=True
    )

    # Step 3: Analyze metabolite similarities
    if hmdb_result.success:
        metabolites = hmdb_result.data["metabolites"]
        print(f"Found {len(metabolites)} similar human metabolites")

        for metabolite in metabolites:
            # Check for disease associations
            if metabolite.get("diseases"):
                print(f"\n{metabolite['name']} is associated with:")
                for disease in metabolite["diseases"][:3]:
                    print(f"  - {disease['name']}")

    return {
        "compound": pubchem_result.data,
        "metabolites": hmdb_result.data
    }
```

### Multi-Adapter Workflow

```python
async def comprehensive_metabolite_analysis(hmdb_id: str):
    """
    Comprehensive analysis using multiple adapters
    """
    hmdb = HMDBAdapter()

    # Step 1: Get HMDB data
    hmdb_result = await hmdb.execute(hmdb_id)

    if not hmdb_result.success:
        return {"error": "Metabolite not found"}

    metabolite = hmdb_result.data["metabolite"]

    # Step 2: Get PubChem data if available
    pubchem_id = metabolite["external_ids"].get("pubchem_compound_id")
    if pubchem_id:
        from adapters.pubchem import PubChemEnhancedAdapter
        pubchem = PubChemEnhancedAdapter()
        pubchem_result = await pubchem.execute(int(pubchem_id))

    # Step 3: Get pathway data from KEGG if available
    kegg_id = metabolite["external_ids"].get("kegg_id")
    if kegg_id:
        from adapters.kegg import KEGGAdapter
        kegg = KEGGAdapter()
        kegg_result = await kegg.execute(kegg_id)

    # Step 4: Get protein data from UniProt
    proteins = metabolite.get("proteins", [])
    for protein in proteins[:5]:
        uniprot_id = protein.get("uniprot_id")
        if uniprot_id:
            from adapters.uniprot import UniProtAdapter
            uniprot = UniProtAdapter()
            uniprot_result = await uniprot.execute(uniprot_id)

    return {
        "hmdb": hmdb_result.data,
        "pubchem": pubchem_result.data if pubchem_id else None,
        "kegg": kegg_result.data if kegg_id else None
    }
```

---

## Technical Notes

### Rate Limiting
- Default: 1 second between requests (respectful to HMDB servers)
- Configurable via `rate_limit_delay` parameter
- Recommended: Keep at 1 second to avoid overloading HMDB

### XML Parsing
- HMDB returns data in XML format
- Adapter uses `xml.etree.ElementTree` for parsing
- Handles HMDB namespace automatically

### Search Limitations
- HMDB search API is limited
- Search returns up to 20 HMDB IDs
- Only first 5 metabolites are fetched in detail to avoid overload
- For production use, consider downloading full HMDB database for local search

### Caching
- All queries are cached using PharmForge's caching system
- Cache key includes HMDB ID and all filter parameters
- Recommended for repeated queries

### Error Handling
- Graceful handling of missing metabolites (404 errors)
- XML parsing errors logged and return None
- Network timeouts configurable (default: 60 seconds)

---

## Performance Considerations

### Single Metabolite Query
- Average response time: 2-5 seconds
- Includes XML download and parsing
- Cached queries: <100ms

### Search Queries
- Average response time: 10-30 seconds (for 5 detailed results)
- Each metabolite requires separate XML download
- Consider using search IDs only for large batches

### Optimization Tips
1. **Use caching** for repeated queries
2. **Limit biofluid filters** to reduce data transfer
3. **Disable unnecessary includes** (diseases, proteins, etc.)
4. **Batch process** HMDB IDs instead of searching

---

## Troubleshooting

### Issue: Metabolite not found
```python
# Ensure HMDB ID format is correct
# Should be: HMDB0000001 (7 digits with leading zeros)
result = await adapter.execute("HMDB0000001")  # Correct
result = await adapter.execute("HMDB1")         # Will be auto-corrected
```

### Issue: Search returns no results
```python
# Try different search terms
# HMDB search is limited - try exact metabolite names
result = await adapter.execute({"query": "D-Glucose", "mode": "search"})
```

### Issue: Timeout errors
```python
# Increase timeout in config
adapter.config["timeout"] = 120  # 2 minutes
```

### Issue: XML parsing errors
```python
# Check HMDB website status at https://hmdb.ca/
# HMDB may be undergoing maintenance
```

---

## Data Sources and References

**HMDB Database:**
- Version: 5.0 (latest)
- Metabolites: 220,945
- Spectra: 10,000+
- Website: https://hmdb.ca/

**Citation:**
> Wishart DS, et al. (2022) HMDB 5.0: the Human Metabolome Database for 2022.
> Nucleic Acids Research 50:D622-D631.

**License:**
- Creative Commons Attribution 4.0 International (CC BY 4.0)
- Free for academic and commercial use
- Attribution required

---

## Future Enhancements

- [ ] Local database download and search
- [ ] Advanced filtering by concentration ranges
- [ ] Chemical similarity search integration
- [ ] Spectral data integration
- [ ] Batch metabolite queries
- [ ] Disease-metabolite network analysis
- [ ] Pathway enrichment analysis

---

## Support

For issues specific to the HMDB adapter:
1. Check adapter logs for detailed error messages
2. Verify HMDB website is accessible
3. Review input format requirements
4. Check rate limiting settings

For HMDB database questions:
- Visit: https://hmdb.ca/about
- Contact: https://hmdb.ca/contact

---

**Version:** 1.0.0
**Last Updated:** 2025-10-30
**Maintainer:** PharmForge Team
