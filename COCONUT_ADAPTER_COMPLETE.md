# COCONUT Adapter - Implementation Complete

**Date:** 2025-10-30
**Adapter:** COCONUT (COlleCtion of Open Natural prodUcTs)
**Status:** ✅ COMPLETE AND TESTED

## Summary

Successfully built and integrated a comprehensive COCONUT adapter for PharmForge that provides access to the world's largest open natural products database (400,000+ compounds) with biological source information, taxonomy, and bioactivity annotations.

## Deliverables

### 1. Main Adapter Implementation
**File:** `adapters/coconut/adapter.py`

**Features:**
- ✅ Full AdapterProtocol compliance
- ✅ Multiple search modes (SMILES, name, organism, ID, InChI/InChIKey)
- ✅ Property-based filtering (MW, LogP, organism)
- ✅ Taxonomy and biological source tracking
- ✅ Biological activity annotations
- ✅ Rate limiting (0.5s default delay)
- ✅ Comprehensive error handling
- ✅ Async/await implementation
- ✅ Automatic caching support

**Search Capabilities:**
1. **Structure Search** - Find natural products similar to a given SMILES
2. **Name Search** - Search by compound name (common or systematic)
3. **Organism Search** - Find all compounds from specific organisms
4. **ID Search** - Retrieve detailed info by COCONUT ID
5. **Property Filters** - Filter by MW, LogP, organism simultaneously

### 2. Module Initialization
**File:** `adapters/coconut/__init__.py`

- Package exports for clean imports
- Version tracking (1.0.0)

### 3. Example Usage
**File:** `adapters/coconut/example_usage.py`

**7 Comprehensive Examples:**
1. Simple SMILES search
2. Advanced query with property filters
3. Search by organism
4. Search by compound name
5. Get compound by COCONUT ID
6. Build natural product library (Lipinski Ro5 analysis)
7. Scaffold mining from natural sources

### 4. Documentation
**File:** `adapters/coconut/README.md`

**Contents:**
- Overview and key features
- Installation (no additional deps needed)
- Complete API documentation
- Input/output format specifications
- 6 detailed use cases with code examples
- Integration with PharmForge pipelines
- Configuration options
- Error handling guide
- Best practices
- Troubleshooting guide

### 5. Comprehensive Tests
**File:** `backend/tests/test_coconut_adapter.py`

**Test Coverage:**
- ✅ Adapter initialization
- ✅ Input validation (SMILES and query dict)
- ✅ Input parsing
- ✅ Property filtering logic
- ✅ Result formatting
- ✅ Cache key generation
- ✅ Metadata extraction
- ✅ Invalid input handling
- ✅ Simple SMILES execution
- ✅ Filtered queries
- ✅ Caching functionality
- ✅ Integration tests (name search, organism search, rate limiting)

**Test Results:** 16/16 PASSED ✅

### 6. Registry Integration
**File:** `backend/core/adapter_registry.py`

- ✅ Imported COCONUTAdapter
- ✅ Added to adapter_classes list
- ✅ Updated adapter count (63 → 64)
- ✅ Successfully registered and accessible via `registry.get('coconut')`

## API Endpoints Used

The adapter integrates with COCONUT's REST API:

```
Base URL: https://coconut.naturalproducts.net/api

Endpoints:
- /search/structure  - Structure similarity search
- /search/name       - Name/text search
- /search/organism   - Organism/taxonomy search
- /compound/{id}     - Get compound details by ID
- /taxonomy          - Browse taxonomy tree
```

## Input Format

### Simple SMILES
```python
result = await adapter("CC(=O)Oc1ccccc1C(=O)O")
```

### Structured Query
```python
query = {
    "query_type": "organism",
    "query": "Streptomyces",
    "filters": {
        "molecular_weight": {"min": 300, "max": 600},
        "logp": {"min": 1, "max": 5}
    },
    "include_taxonomy": True,
    "include_activities": True,
    "max_results": 50
}
result = await adapter(query)
```

## Output Format

```python
{
    "results": [
        {
            "coconut_id": "CNP0123456",
            "name": "Taxol",
            "smiles": "CC1=C2[C@@H]...",
            "inchi": "InChI=1S/C47H51NO14/...",
            "molecular_formula": "C47H51NO14",
            "molecular_weight": 853.91,
            "natural_source": {
                "organism": "Taxus brevifolia",
                "common_name": "Pacific Yew",
                "taxonomy": "Plantae;Pinophyta;..."
            },
            "properties": {
                "logp": 4.2,
                "hba": 14,
                "hbd": 4,
                "rotatable_bonds": 12,
                "tpsa": 221.3
            },
            "biological_activities": ["Anticancer", "Microtubule stabilizer"],
            "references": ["PMID:12345678"]
        }
    ],
    "num_results": 1,
    "query_type": "name",
    "query": "Taxol",
    "warnings": []
}
```

## Use Cases

### 1. Natural Product Library Screening
Build drug-like natural product libraries for virtual screening.

### 2. Scaffold Mining from Nature
Find natural scaffolds similar to known bioactives for scaffold hopping.

### 3. Biosource Tracking
Track which organisms produce specific compound classes.

### 4. Taxonomic Analysis
Analyze natural product distribution across taxonomic groups.

### 5. Traditional Medicine Validation
Validate traditional medicine claims with natural product data.

### 6. Natural Product-Inspired Drug Design
Generate analogs inspired by natural scaffolds.

## Integration with PharmForge Pipeline

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()

# Step 1: Find natural products from Streptomyces
pipeline.add_step("coconut", {
    "query_type": "organism",
    "query": "Streptomyces",
    "filters": {"molecular_weight": {"min": 300, "max": 700}},
    "max_results": 100
})

# Step 2: Calculate ADMET properties
pipeline.add_step("tdc_admet", {
    "properties": ["Caco2_Wang", "Solubility_AqSolDB"]
})

# Step 3: Predict targets
pipeline.add_step("swisstarget", {"organism": "Homo sapiens"})

# Step 4: Dock to target
pipeline.add_step("vina_docking", {
    "receptor": "protein.pdb",
    "center": [10.0, 15.0, 20.0],
    "size": [20.0, 20.0, 20.0]
})

results = await pipeline.execute()
```

## Configuration

```python
adapter = COCONUTAdapter()
adapter.config.update({
    "rate_limit_delay": 1.0,               # Delay between requests
    "timeout": 120,                         # Request timeout
    "max_results": 100,                     # Default max results
    "default_similarity_threshold": 0.85    # Similarity cutoff
})
```

## Complementary Adapters

Works well with:
- **ChEMBL** - Cross-reference bioactivity data
- **PubChem** - Validate structures and get additional properties
- **RDKit** - Calculate additional descriptors
- **TDC ADMET** - Predict ADMET properties
- **Vina/DiffDock** - Dock natural products to targets
- **AiZynthFinder** - Plan synthesis routes

## Key Technical Features

### Rate Limiting
- Default: 0.5 seconds between requests
- Configurable via `config["rate_limit_delay"]`
- Respectful to public API

### Caching
- Automatic caching via PharmForge's cache system
- Deterministic cache keys (SHA256)
- Cache hit tracking

### Error Handling
- Graceful timeout handling
- API error detection
- Input validation with clear error messages
- Warning system for partial results

### Async Performance
- Full async/await implementation
- Non-blocking I/O operations
- Efficient concurrent execution

## Testing Results

```
======================== test session starts =========================
platform win32 -- Python 3.12.7, pytest-8.4.2, pluggy-1.6.0
collected 16 items

backend/tests/test_coconut_adapter.py::TestCOCONUTAdapter::
  test_adapter_initialization                   PASSED [  6%]
  test_validate_input_simple_smiles            PASSED [ 12%]
  test_validate_input_query_dict               PASSED [ 18%]
  test_parse_input_simple_smiles               PASSED [ 25%]
  test_parse_input_query_dict                  PASSED [ 31%]
  test_apply_property_filters                  PASSED [ 37%]
  test_format_compound_result                  PASSED [ 43%]
  test_cache_key_generation                    PASSED [ 50%]
  test_metadata                                PASSED [ 56%]
  test_execute_invalid_input                   PASSED [ 62%]
  test_execute_simple_smiles                   PASSED [ 68%]
  test_execute_with_filters                    PASSED [ 75%]
  test_caching                                 PASSED [ 81%]

backend/tests/test_coconut_adapter.py::TestCOCONUTIntegration::
  test_name_search_integration                 PASSED [ 87%]
  test_organism_search_integration             PASSED [ 93%]
  test_rate_limiting                           PASSED [100%]

=================== 16 passed, 1 warning in 10.80s ==================
```

## Files Created

```
adapters/coconut/
├── __init__.py              # Module initialization (483 bytes)
├── adapter.py               # Main adapter (17,467 bytes)
├── example_usage.py         # Usage examples (9,461 bytes)
└── README.md                # Documentation (13,714 bytes)

backend/tests/
└── test_coconut_adapter.py  # Test suite (10,887 bytes)

Total: 5 files, ~52 KB of production code
```

## Database Information

- **Database:** COCONUT (COlleCtion of Open Natural prodUcTs)
- **Size:** 400,000+ natural product structures
- **License:** CC0 (public domain)
- **URL:** https://coconut.naturalproducts.net/
- **API:** Public REST API (no key required)
- **Coverage:** Natural products from bacteria, fungi, plants, animals

## Quality Metrics

- **AdapterProtocol Compliance:** ✅ 100%
- **Test Coverage:** ✅ 16/16 tests passing
- **Documentation:** ✅ Complete (README + examples)
- **Error Handling:** ✅ Comprehensive
- **Code Quality:** ✅ Type hints, docstrings, logging
- **Performance:** ✅ Async, cached, rate-limited

## Success Criteria - ALL MET ✅

- ✅ Adapter follows AdapterProtocol exactly
- ✅ All methods implemented with proper signatures
- ✅ Tests pass with real data
- ✅ Error handling covers common cases
- ✅ Adapter registered and accessible via API
- ✅ Documentation added to adapter file
- ✅ Example usage provided
- ✅ Integration with PharmForge pipeline demonstrated

## Next Steps (Optional Enhancements)

While the adapter is complete and production-ready, future enhancements could include:

1. **Batch Search** - Support searching multiple compounds simultaneously
2. **Advanced Filters** - Add more property filters (HBA, HBD, TPSA, etc.)
3. **Substructure Search** - Add SMARTS-based substructure matching
4. **Export Formats** - Support SDF, MOL file downloads
5. **Taxonomy Browser** - Browse entire taxonomy tree
6. **Activity Enrichment** - Statistical analysis of activity distributions
7. **Source Map** - Geographic mapping of natural sources
8. **Fragment Analysis** - Natural product fragment library generation

## References

- **COCONUT Database:** https://coconut.naturalproducts.net/
- **COCONUT Paper:** Sorokina et al., J. Cheminformatics (2021)
- **Natural Products Atlas:** https://www.npatlas.org/
- **PharmForge AdapterProtocol:** `backend/core/adapters/protocol.py`

## Conclusion

The COCONUT adapter is a fully functional, well-tested, and documented integration that brings the power of the world's largest open natural products database to PharmForge. It enables researchers to:

- Screen natural product libraries
- Mine scaffolds from nature
- Track biosources and taxonomy
- Design natural product-inspired drugs
- Validate traditional medicine

The adapter is production-ready and can be used immediately in PharmForge pipelines for natural product discovery and drug design workflows.

---

**Adapter Name:** `coconut`
**Version:** 1.0.0
**Type:** API
**Status:** Production Ready ✅
**License:** PharmForge license (adapter), CC0 (database)

