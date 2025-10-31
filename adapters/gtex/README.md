# GTEx Adapter

Adapter for the GTEx (Genotype-Tissue Expression) Portal API, providing access to gene expression data across human tissues and eQTL (expression quantitative trait loci) information.

## Overview

The GTEx (Genotype-Tissue Expression) project is a comprehensive resource for studying tissue-specific gene expression and regulation. This adapter queries the GTEx Portal API v2 to retrieve:

- Median gene expression levels across 54 human tissues
- Tissue-specific expression patterns
- Expression quantitative trait loci (eQTLs)
- Gene annotations and metadata

## API Documentation

- Base URL: `https://gtexportal.org/api/v2`
- API Docs: https://gtexportal.org/api/v2/redoc
- Dataset: GTEx v8 (default, can be configured)
- Rate Limit: ~5 requests/second (self-imposed)
- No authentication required

## Installation

No additional dependencies required beyond the base adapter requirements:
- `aiohttp` for async HTTP requests
- `asyncio` for async/await support

## Usage

### Basic Example - Gene Expression Across Tissues

```python
import asyncio
from adapters.gtex import GTExAdapter

async def main():
    adapter = GTExAdapter()

    # Query gene expression across all tissues
    result = await adapter("BRCA2", query_type="expression")

    if result.success:
        data = result.data
        print(f"Gene: {data['gene_symbol']}")
        print(f"Tissues with expression: {data['num_tissues_expressed']}")
        print(f"Highest expression: {data['max_expression_tissue']}")
        print(f"  TPM: {data['max_expression_value']:.2f}")

        print(f"\nTop 10 Expressing Tissues:")
        sorted_tissues = sorted(
            data['median_expression'].items(),
            key=lambda x: x[1],
            reverse=True
        )
        for tissue, tpm in sorted_tissues[:10]:
            print(f"  {tissue}: {tpm:.2f} TPM")

    else:
        print(f"Error: {result.error}")

asyncio.run(main())
```

### Tissue-Specific Expression

```python
# Get expression in a specific tissue
result = await adapter(
    "TP53",
    query_type="tissue",
    tissue="Brain_Cortex"
)

if result.success:
    print(f"Expression in {result.data['tissue']}:")
    print(result.data['expression_data'])
```

### eQTL Analysis

```python
# Get eQTL data for a gene
result = await adapter("APOE", query_type="eqtl")

if result.success:
    print(f"eQTLs found: {result.data['num_eqtls']}")
    for eqtl in result.data['eqtls'][:5]:
        print(f"  Variant: {eqtl.get('variantId')}")
        print(f"  P-value: {eqtl.get('pValue')}")
        print(f"  Tissue: {eqtl.get('tissueSiteDetailId')}")
```

## Input Formats

The adapter accepts:
- **Gene symbols**: e.g., "BRCA2", "TP53", "APOE"
- **Ensembl gene IDs**: e.g., "ENSG00000139618"

## Output Format

### Expression Query (query_type="expression")

```python
{
    "success": bool,
    "data": {
        "gene_symbol": str,              # Gene symbol/ID
        "query_type": "expression",
        "dataset": str,                  # GTEx dataset version
        "tissues": List[str],            # List of tissue names
        "median_expression": {           # TPM values by tissue
            "Brain_Cortex": float,
            "Heart_Left_Ventricle": float,
            # ... more tissues
        },
        "max_expression_tissue": str,    # Tissue with highest expression
        "max_expression_value": float,   # Maximum TPM value
        "num_tissues_expressed": int,    # Count where TPM > 1
        "gene_info": {                   # Optional gene metadata
            "geneSymbol": str,
            "geneId": str,
            "chromosome": str,
            "start": int,
            "end": int
        }
    },
    "error": Optional[str],
    "cache_hit": bool,
    "metadata": {
        "source": "gtex",
        "gene": str,
        "query_type": str,
        "dataset": str,
        "adapter_version": str
    }
}
```

### Tissue Query (query_type="tissue")

```python
{
    "data": {
        "gene_symbol": str,
        "query_type": "tissue",
        "tissue": str,
        "expression_data": {
            # Tissue-specific expression details
        }
    }
}
```

### eQTL Query (query_type="eqtl")

```python
{
    "data": {
        "gene_symbol": str,
        "query_type": "eqtl",
        "num_eqtls": int,
        "eqtls": [
            {
                "variantId": str,
                "pValue": float,
                "tissueSiteDetailId": str,
                "slope": float,
                "geneSymbol": str
            }
        ]
    }
}
```

## Configuration

Default configuration:
```python
{
    "rate_limit_delay": 0.2,    # 200ms between requests
    "timeout": 60,               # 60 second timeout
    "max_tissues": 100           # Maximum tissues to retrieve
}
```

## GTEx Tissues

GTEx v8 includes 54 tissues:
- **Brain**: Amygdala, Anterior cingulate cortex, Caudate, Cerebellar Hemisphere, Cerebellum, Cortex, Frontal Cortex, Hippocampus, Hypothalamus, Nucleus accumbens, Putamen, Spinal cord, Substantia nigra
- **Heart**: Atrial Appendage, Left Ventricle
- **Muscle**: Skeletal
- **Blood**: Cells, EBV-transformed lymphocytes, Whole Blood
- **Adipose**: Subcutaneous, Visceral
- **Artery**: Aorta, Coronary, Tibial
- **Lung**
- **Skin**: Sun Exposed, Not Sun Exposed
- **Digestive**: Colon, Esophagus (Gastroesophageal Junction, Mucosa, Muscularis), Liver, Minor Salivary Gland, Pancreas, Small Intestine, Stomach
- **Reproductive**: Breast, Cervix, Fallopian Tube, Ovary, Prostate, Testis, Uterus, Vagina
- **Endocrine**: Adrenal Gland, Pituitary, Thyroid
- **Kidney**: Cortex, Medulla
- **Nerve**: Tibial
- **Spleen**

## Query Types

### 1. Expression Query (default)

Get median expression across all tissues:
```python
result = await adapter("BRCA2", query_type="expression")
```

### 2. Tissue-Specific Query

Get expression in a specific tissue:
```python
result = await adapter(
    "TP53",
    query_type="tissue",
    tissue="Liver"
)
```

### 3. eQTL Query

Get expression quantitative trait loci:
```python
result = await adapter("APOE", query_type="eqtl")

# eQTL for specific tissue
result = await adapter(
    "APOE",
    query_type="eqtl",
    tissue="Brain_Cortex"
)
```

## Understanding Expression Values

- **TPM (Transcripts Per Million)**: Normalized gene expression metric
  - TPM < 1: Low/no expression
  - TPM 1-10: Low expression
  - TPM 10-100: Moderate expression
  - TPM > 100: High expression

## Error Handling

The adapter handles:
- Invalid gene symbols/IDs (returns error)
- API timeouts (60 second timeout)
- Rate limiting (200ms delay between requests)
- Network errors (returns error with details)
- Missing data (returns appropriate error message)

## Caching

Results are automatically cached using the AdapterProtocol caching system:
- Cache key: SHA256 hash of (adapter name, version, input, parameters)
- Cache duration: Configurable in backend cache settings
- Disable cache: `await adapter(gene_id, use_cache=False)`

## Examples

### Finding Tissue-Specific Genes

```python
# Find genes highly expressed in brain
result = await adapter("SYN1", query_type="expression")  # Synapsin 1

if result.success:
    brain_tissues = {
        tissue: tpm
        for tissue, tpm in result.data['median_expression'].items()
        if 'Brain' in tissue
    }

    print("Brain expression:")
    for tissue, tpm in sorted(brain_tissues.items(), key=lambda x: x[1], reverse=True):
        print(f"  {tissue}: {tpm:.2f} TPM")
```

### Comparing Expression Patterns

```python
# Compare two genes
gene1_result = await adapter("BRCA1", query_type="expression")
gene2_result = await adapter("BRCA2", query_type="expression")

if gene1_result.success and gene2_result.success:
    # Find tissues where both are highly expressed
    common_high = []
    for tissue in gene1_result.data['tissues']:
        tpm1 = gene1_result.data['median_expression'].get(tissue, 0)
        tpm2 = gene2_result.data['median_expression'].get(tissue, 0)
        if tpm1 > 10 and tpm2 > 10:
            common_high.append((tissue, tpm1, tpm2))

    print("Tissues with high expression of both:")
    for tissue, tpm1, tpm2 in common_high:
        print(f"  {tissue}: BRCA1={tpm1:.1f}, BRCA2={tpm2:.1f}")
```

### Finding eQTLs

```python
# Find regulatory variants affecting gene expression
result = await adapter("APOE", query_type="eqtl")

if result.success:
    # Filter significant eQTLs (p < 1e-5)
    significant = [
        eqtl for eqtl in result.data['eqtls']
        if eqtl.get('pValue', 1) < 1e-5
    ]

    print(f"Significant eQTLs: {len(significant)}")
    for eqtl in sorted(significant, key=lambda x: x.get('pValue', 1))[:10]:
        print(f"\nVariant: {eqtl.get('variantId')}")
        print(f"P-value: {eqtl.get('pValue'):.2e}")
        print(f"Tissue: {eqtl.get('tissueSiteDetailId')}")
        print(f"Effect size: {eqtl.get('slope', 'N/A')}")
```

### Identifying Druggable Targets

```python
# Find genes highly expressed in disease-relevant tissues
result = await adapter("IL6", query_type="expression")

if result.success:
    # Focus on immune-related tissues
    immune_tissues = [
        "Whole_Blood",
        "Spleen",
        "Cells_EBV-transformed_lymphocytes"
    ]

    print("Expression in immune tissues:")
    for tissue in immune_tissues:
        tpm = result.data['median_expression'].get(tissue, 0)
        print(f"  {tissue}: {tpm:.2f} TPM")

    # High expression in immune tissues suggests it's a good target
    avg_immune = sum(
        result.data['median_expression'].get(t, 0)
        for t in immune_tissues
    ) / len(immune_tissues)

    if avg_immune > 10:
        print("\nGood immune-related drug target!")
```

### Tissue Specificity Score

```python
# Calculate tissue specificity (tau)
result = await adapter("INS", query_type="expression")  # Insulin

if result.success:
    expressions = list(result.data['median_expression'].values())
    max_exp = max(expressions)

    if max_exp > 0:
        tau = sum(1 - (exp / max_exp) for exp in expressions) / (len(expressions) - 1)
        print(f"Tissue specificity (tau): {tau:.3f}")
        # tau close to 1 = highly specific
        # tau close to 0 = ubiquitously expressed

        if tau > 0.8:
            print("Highly tissue-specific gene")
        elif tau < 0.3:
            print("Ubiquitously expressed gene")
```

## Limitations

- **GTEx Version**: Currently uses v8 (can be configured)
- **Data Updates**: GTEx is periodically updated; API may change
- **Rate Limiting**: Self-imposed ~5 requests/second
- **Sample Size**: Varies by tissue (30-800+ donors)
- **Age/Sex**: Data is aggregated; individual-level data not available via API
- **Disease States**: GTEx uses healthy donor tissues only

## Dataset Information

- **GTEx v8**: Current default
  - 17,382 samples
  - 54 tissues
  - 948 donors
  - Released: 2019

## Version History

- **1.0.0** (2025-10-25): Initial release
  - Gene expression queries
  - Tissue-specific expression
  - eQTL data retrieval
  - Multiple tissues support
  - Caching support
  - Rate limiting
