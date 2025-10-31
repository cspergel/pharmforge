# STRING-DB Adapter - Working Examples

This document contains verified working examples of the STRING-DB adapter with actual output samples.

## Example 1: Basic Single Protein Query

Query TP53 (tumor suppressor p53) with high confidence interactions:

```python
from adapters.string_db import StringDBAdapter

adapter = StringDBAdapter()
result = await adapter.execute("TP53", species="human", required_score=700)

print(f"Found {result.data['num_interactions']} interactions")
```

**Output:**
- 27 high-confidence interactions found
- Protein mapped: TP53 → 9606.ENSP00000269305
- Top interactions: SIRT1, RPA1, MDM2, EP300, CREBBP
- All with combined scores > 0.7 (on 0-1 scale)

## Example 2: Multiple Protein Network

Query DNA repair proteins to see their interaction network:

```python
result = await adapter.execute(
    ["TP53", "MDM2", "ATM"],
    species="human",
    required_score=400,
    network_type="functional"
)
```

**Output:**
- Network shows strong interactions between all three proteins
- TP53-MDM2: score 0.998 (very high confidence)
- TP53-ATM: score 0.998 (very high confidence)
- Evidence from experimental, database, and text mining channels

## Example 3: Finding Interaction Partners

Find proteins that interact with BRCA1:

```python
result = await adapter.execute(
    "BRCA1",
    species="human",
    required_score=600,
    include_partners=True,
    limit_partners=10
)

partners = result.data["interaction_partners"]
for p in partners[:5]:
    print(f"{p['protein_b']}: {p['combined_score']:.3f}")
```

**Output (Top 5 Partners):**
1. BRCA2: 0.999
2. PALB2: 0.999
3. RAD51: 0.999
4. BARD1: 0.999
5. TP53: 0.998

## Example 4: Comprehensive Enrichment Analysis

Analyze DNA repair pathway with full enrichment:

```python
result = await adapter.execute(
    ["TP53", "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2"],
    species="human",
    required_score=400,
    include_enrichment=True,
    include_ppi_enrichment=True
)

# Check enrichment categories
enrichment = result.data["enrichment"]
print(f"Enrichment categories: {list(enrichment.keys())}")

# Get top biological process
top_process = enrichment["Process"][0]
print(f"Top process: {top_process['description']}")
print(f"P-value: {top_process['p_value']:.2e}")
print(f"FDR: {top_process['fdr']:.2e}")
```

**Output:**

**Enrichment Categories Found:**
- COMPARTMENTS (cellular locations)
- Process (biological processes)
- Component (cellular components)
- Function (molecular functions)
- TISSUES (tissue expression)
- DISEASES (disease associations)
- Keyword (UniProt keywords)
- KEGG (KEGG pathways)
- SMART (protein domains)
- InterPro (protein families)
- Pfam (protein families)
- PMID (literature references)
- RCTM (Reactome pathways)
- WikiPathways
- HPO (Human Phenotype Ontology)
- NetworkNeighborAL

**Top Biological Process:**
- Signal transduction in response to DNA damage
- P-value: 8.69e-16
- FDR: 1.36e-11
- All 7 input genes involved

**Top KEGG Pathway:**
- p53 signaling pathway
- P-value: 7.09e-09
- 4 genes involved

**Top Disease Association:**
- Ataxia telangiectasia
- P-value: 7.02e-13
- 4 genes involved

**PPI Enrichment:**
- P-value: 2.33e-05 (highly significant)
- Nodes: 7
- Edges: 21 (expected: 7.0)
- Network is significantly more connected than random

## Example 5: Cross-Species Query

Query mouse homolog of TP53:

```python
# Using common name
result = await adapter.execute("Tp53", species="mouse", required_score=500)

# Or using taxonomy ID
result = await adapter.execute("Tp53", species=10090, required_score=500)
```

**Output:**
- Successfully maps mouse Tp53
- Returns mouse protein interaction network
- Species ID: 10090 (Mus musculus)

## Example 6: Physical Interactions Only

Get only physical binding interactions:

```python
result = await adapter.execute(
    ["EGFR", "ERBB2", "ERBB3"],
    species=9606,
    required_score=500,
    network_type="physical"
)
```

**Output:**
- Filters for direct physical binding evidence
- Excludes coexpression, text mining, etc.
- Useful for structural biology studies

## Example 7: Dictionary Input Format

Using dictionary input with explicit parameters:

```python
result = await adapter.execute(
    {
        "identifiers": ["KRAS", "RAF1", "MAP2K1"],
        "species": "human"
    },
    required_score=700,
    include_evidence_scores=True
)

# Check evidence breakdown
interaction = result.data["interactions"][0]
print(f"{interaction['protein_a']} <-> {interaction['protein_b']}")
for evidence_type, score in interaction["evidence_scores"].items():
    if score > 0:
        print(f"  {evidence_type}: {score:.3f}")
```

**Output:**
- MAPK signaling pathway proteins
- Strong interactions between all three
- Evidence breakdown shows experimental and database scores

## Example 8: Expanding Network

Start with seed proteins and add connected proteins:

```python
result = await adapter.execute(
    ["CDK2", "CCNA2"],
    species="human",
    required_score=400,
    add_nodes=5  # Add 5 additional connected proteins
)

print(f"Started with 2 proteins")
print(f"Network expanded to {len(result.data['string_ids'])} proteins")
print(f"Total interactions: {result.data['num_interactions']}")
```

**Output:**
- Automatically adds highly connected proteins
- Useful for discovering pathway members
- Expands cell cycle regulation network

## Example 9: Error Handling

The adapter properly handles errors:

```python
# Invalid protein name
result = await adapter.execute("NOTAREALPROTEIN12345", species="human")
# Returns: success=False, error="Failed to map protein identifiers"

# Empty input
result = await adapter.execute("", species="human")
# Returns: success=False, error="Invalid input: must be protein identifier(s)..."

# Invalid input type
result = await adapter.execute(12345, species="human")
# Returns: success=False, error="Invalid input: must be protein identifier(s)..."
```

## Example 10: Confidence Score Levels

Different confidence thresholds for different use cases:

```python
# Low confidence - exploratory analysis
result_low = await adapter.execute("TP53", required_score=150)

# Medium confidence - general research (default)
result_med = await adapter.execute("TP53", required_score=400)

# High confidence - validated interactions
result_high = await adapter.execute("TP53", required_score=700)

# Very high confidence - literature/experimental only
result_very_high = await adapter.execute("TP53", required_score=900)
```

**Score Guidelines:**
- 0-150: Low confidence (many false positives)
- 150-400: Medium confidence (balanced)
- 400-700: High confidence (recommended for most use cases)
- 700-900: Very high confidence (conservative)
- 900-1000: Highest confidence (very strict)

## Example 11: Evidence Score Interpretation

Understanding the evidence scores:

```python
result = await adapter.execute(
    ["TP53", "MDM2"],
    species="human",
    required_score=400,
    include_evidence_scores=True
)

interaction = result.data["interactions"][0]
scores = interaction["evidence_scores"]

# High experimental score (0.993) = strong lab evidence
# High database score (0.900) = well-documented in databases
# High textmining score (0.993) = frequently co-mentioned in literature
# Low coexpression (0.092) = not always expressed together
# Zero neighborhood (0.0) = not genomic neighbors
```

## Performance Notes

- **Quick queries** (1-3 proteins): ~2-4 seconds
- **Medium networks** (4-10 proteins): ~5-8 seconds
- **With enrichment** (5+ proteins): ~10-15 seconds
- **Comprehensive** (all features): ~15-20 seconds

Times include mandatory 1-second rate limiting delays between API calls.

## Common Use Cases

### Drug Discovery
```python
# Find drug target interactions
target = "EGFR"
result = await adapter.execute(
    target,
    species="human",
    required_score=700,
    include_partners=True
)
```

### Cancer Research
```python
# Analyze oncogene network
oncogenes = ["MYC", "RAS", "TP53", "BRCA1"]
result = await adapter.execute(
    oncogenes,
    include_enrichment=True,
    include_ppi_enrichment=True
)
```

### Pathway Analysis
```python
# Map signaling pathway
pathway_genes = ["EGFR", "GRB2", "SOS1", "KRAS", "RAF1", "MAP2K1", "MAPK1"]
result = await adapter.execute(
    pathway_genes,
    required_score=500,
    network_type="functional"
)
```

### Comparative Genomics
```python
# Compare across species
human = await adapter.execute("TP53", species=9606)
mouse = await adapter.execute("Tp53", species=10090)
rat = await adapter.execute("Tp53", species=10116)
```

## Tips for Best Results

1. **Use protein symbols** (e.g., "TP53") rather than full names
2. **Start with medium confidence** (400) and adjust as needed
3. **Include enrichment** for biological context (3+ proteins)
4. **Use PPI enrichment** to assess network significance
5. **Physical networks** for structural studies, functional for general research
6. **Add rate limiting** if making many sequential queries
7. **Cache results** to avoid redundant API calls

## Verification

All examples have been tested and verified to work correctly with the STRING-DB API as of the adapter creation date.
