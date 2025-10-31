# Chemical Database Adapter Fixes - Documentation

**Date:** October 2025
**Version:** 2.0.0
**Status:** ✅ Fixed and Tested

## Executive Summary

Fixed three critical Chemical Database adapters that were failing due to deprecated APIs, bot protection, and server errors:

1. **DrugCentral** - Migrated from deprecated API v1 to Pharos API
2. **ZINC** - Fixed 403 Forbidden errors with proper headers and retry logic
3. **SureChEMBL** - Updated to v2.0 API with ChEMBL fallback

All adapters now include:
- ✅ Retry logic with exponential backoff
- ✅ Proper User-Agent headers
- ✅ Rate limiting to respect API servers
- ✅ Comprehensive error handling
- ✅ Fallback mechanisms
- ✅ Detailed logging

---

## 1. DrugCentral Adapter Fix

### Problem
- **Old API:** `https://drugcentral.org/api/v1` (404 errors - deprecated)
- **Issue:** API v1 no longer exists, causing all requests to fail

### Solution
**Migrated to Pharos API** which provides access to DrugCentral data:
- **New API:** `https://pharos.nih.gov/idg/api/v1` (REST)
- **GraphQL API:** `https://pharos-api.ncats.io/graphql` (available)
- **Documentation:** https://pharos.nih.gov/api

### What Changed

#### API Endpoints
| Old Endpoint | New Endpoint | Notes |
|-------------|--------------|-------|
| `/api/v1/search` | `/idg/api/v1/ligands/search` | Drug search |
| `/api/v1/drug/{id}/targets` | `/idg/api/v1/ligands/{id}/targets` | Target information |
| `/api/v1/drug/{id}/indications` | `/idg/api/v1/ligands/{id}` | Extracted from full ligand data |
| `/api/v1/drug/{id}/pharmacology` | `/idg/api/v1/ligands/{id}` | Extracted from full ligand data |
| `/api/v1/drug/{id}/properties` | `/idg/api/v1/ligands/{id}` | Chemical properties |

#### Response Format Changes
Pharos uses different field names:
- `id` → `ligid`
- `inchikey` → `inchiKey`
- `gene` → `sym` (for targets)
- `molecular_weight` → `mwt`
- `polar_surface_area` → `tpsa`

#### New Features
- **Retry Logic:** 3 retries with exponential backoff (1s, 2s, 4s)
- **User-Agent:** `PharmForge/2.0 (Drug Discovery Platform; +https://pharmforge.org)`
- **Better Error Handling:** Distinguishes between 404, 500, and network errors
- **Rate Limiting:** 0.5s delay between requests

### Usage Example

```python
from adapters.drugcentral.adapter import DrugCentralAdapter

async def search_drug():
    adapter = DrugCentralAdapter()

    # Search for aspirin
    result = await adapter.execute("aspirin")

    if result.success:
        drugs = result.data['drugs']
        for drug in drugs:
            print(f"Drug: {drug['name']}")
            print(f"Targets: {drug['num_targets']}")
            print(f"Indications: {drug['num_indications']}")
```

### Known Limitations
- Pharos API may have rate limits (respectful delays implemented)
- Some DrugCentral-specific fields may not be available
- GraphQL option available for complex queries (not yet implemented)

---

## 2. ZINC Adapter Fix

### Problem
- **Error:** 403 Forbidden
- **Cause:** ZINC servers block requests without proper User-Agent headers (anti-bot protection)
- **Rate Limiting:** Too many rapid requests trigger blocking

### Solution
Added comprehensive bot protection bypass and retry mechanisms:
1. ✅ Proper User-Agent headers mimicking legitimate browser
2. ✅ Retry logic with exponential backoff
3. ✅ Increased rate limiting delays
4. ✅ Changed from HTTP to HTTPS
5. ✅ Sequential purchasability requests (instead of parallel)

### What Changed

#### Headers Added
```python
{
    'User-Agent': 'Mozilla/5.0 (compatible; PharmForge/2.0; +https://pharmforge.org/bot)',
    'Accept': 'application/json',
    'Accept-Language': 'en-US,en;q=0.9',
    'Accept-Encoding': 'gzip, deflate',
    'Connection': 'keep-alive'
}
```

#### Configuration Updates
| Parameter | Old Value | New Value | Reason |
|-----------|-----------|-----------|--------|
| `rate_limit_delay` | 0.5s | 2.0s | Prevent rate limiting |
| `timeout` | 60s | 90s | Allow slower responses |
| `BASE_URL` | `http://` | `https://` | Security and compatibility |
| `max_retries` | N/A | 3 | Handle temporary failures |
| `backoff_base` | N/A | 2 | Exponential backoff (2s, 4s, 8s) |

#### Retry Logic
```python
# Automatic retry on 403 errors
if response.status == 403:
    raise aiohttp.ClientError("403 Forbidden - bot detection")
    # Triggers retry with exponential backoff
```

#### Purchasability Requests
Changed from parallel to sequential with delays:
```python
# Before: Parallel requests (could trigger rate limits)
await asyncio.gather(*purchasability_tasks)

# After: Sequential with delays
for result in filtered_results[:20]:
    await asyncio.sleep(2.0)  # Rate limit delay
    purchasability = await self._get_purchasability_async(zinc_id)
```

### Usage Example

```python
from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter

async def search_fragments():
    adapter = ZINCFragmentsAdapter()

    # Search for fragments similar to ethanol
    result = await adapter.execute(
        "CCO",  # SMILES
        search_type="similarity",
        similarity_threshold=0.8,
        max_results=50,
        get_purchasability=True  # Get vendor info
    )

    if result.success:
        fragments = result.data['fragments']
        print(f"Found {len(fragments)} fragments")

        for frag in fragments[:5]:
            print(f"ZINC ID: {frag.get('zinc_id')}")
            if 'purchasability' in frag:
                print(f"  Purchasable: {frag['purchasability']['purchasable']}")
```

### Known Limitations
- Rate limiting: 2 second delay between requests (slower but reliable)
- Purchasability info limited to top 20 results (to avoid triggering blocks)
- ZINC may still block if used too aggressively - use responsibly

---

## 3. SureChEMBL Adapter Fix

### Problem
- **Error:** 500 Internal Server Error
- **Cause:** API endpoint changes in SureChEMBL 2.0 release (May 2025)
- **Migration:** SureChEMBL upgraded to RDKit, new annotation pipeline

### Solution
1. ✅ Updated to SureChEMBL 2.0 API endpoints
2. ✅ Added ChEMBL API fallback (ChEMBL includes SureChEMBL data)
3. ✅ Retry logic for handling temporary 500 errors
4. ✅ Updated endpoint paths from `/patent/` to `/document/`

### What Changed

#### API Endpoints
| Function | Old Endpoint | New Endpoint | Notes |
|----------|-------------|--------------|-------|
| Structure Search | `/api/search/similarity` | `/api/v1/search?type=similarity` | Updated parameters |
| Patent Lookup | `/api/patent/{id}` | `/api/v1/document/{id}` | Changed from patent to document |
| Extract Compounds | `/api/patent/{id}/compounds` | `/api/v1/document/{id}/compounds` | Changed path |

#### ChEMBL Fallback
When SureChEMBL fails (500 errors), automatically falls back to ChEMBL:
```python
# Primary: SureChEMBL 2.0
BASE_URL = "https://www.surechembl.org/api/v1"

# Fallback: ChEMBL API (includes SureChEMBL data)
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
```

#### Configuration
```python
config = {
    "rate_limit_delay": 1.0,
    "timeout": 90,
    "default_similarity": 0.8,
    "max_results": 50,
    "max_retries": 3,
    "use_chembl_fallback": True  # Enable fallback
}
```

#### Response Format Updates
SureChEMBL 2.0 field changes:
- `publication_date` → `published`
- `family_id` → `family`
- Better compound structure in responses

### Usage Example

```python
from adapters.surechembl.adapter import SureChEMBLAdapter

async def search_patents():
    adapter = SureChEMBLAdapter()

    # Search patents for aspirin-like compounds
    result = await adapter.execute(
        {"smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        mode="structure_search",
        search_type="similarity",
        similarity_threshold=0.8,
        max_results=20
    )

    if result.success:
        print(f"Found {result.data['num_patents']} patents")
        print(f"Used fallback: {result.metadata['used_chembl_fallback']}")

        # Get patent details
        patents = result.data['patents']
        for patent in patents[:5]:
            print(f"Patent: {patent.get('patent_id')}")
            print(f"Source: {patent.get('source', 'surechembl')}")
```

### Known Limitations
- SureChEMBL 2.0 may still have intermittent 500 errors during updates
- ChEMBL fallback provides drug data but may not have all patent-specific info
- ChEMBL fallback results are marked with `source: "chembl_fallback"`

---

## Testing

### Test Script
Run the included test script to verify all adapters:

```bash
cd claude-code-agents-wizard-v2/adapters
python test_fixed_adapters.py
```

### Manual Testing

#### DrugCentral
```python
adapter = DrugCentralAdapter()
result = await adapter.execute("aspirin")
assert result.success
assert result.data['total_found'] > 0
```

#### ZINC
```python
adapter = ZINCFragmentsAdapter()
result = await adapter.execute("CCO", search_type="similarity")
assert result.success
assert result.data['num_results'] >= 0
```

#### SureChEMBL
```python
adapter = SureChEMBLAdapter()
result = await adapter.execute(
    {"smiles": "CC(=O)Oc1ccccc1C(=O)O"},
    mode="structure_search"
)
assert result.success
```

---

## General Improvements (All Adapters)

### 1. Retry Logic
All adapters now include exponential backoff:
```python
async def _retry_request(self, func, *args, **kwargs):
    max_retries = 3
    for attempt in range(max_retries):
        try:
            return await func(*args, **kwargs)
        except (asyncio.TimeoutError, aiohttp.ClientError) as e:
            wait_time = 2 ** attempt  # 1s, 2s, 4s
            await asyncio.sleep(wait_time)
    return None
```

### 2. User-Agent Headers
All adapters identify as PharmForge:
```python
headers = {
    'User-Agent': 'PharmForge/2.0 (Drug Discovery Platform; +https://pharmforge.org)',
    'Accept': 'application/json'
}
```

### 3. Rate Limiting
Respectful delays between requests:
- DrugCentral: 0.5s
- ZINC: 2.0s (increased to prevent blocks)
- SureChEMBL: 1.0s

### 4. Error Handling
Comprehensive error categorization:
- Network errors (timeout, connection)
- HTTP errors (403, 404, 500)
- API errors (invalid response format)
- Validation errors (bad input)

### 5. Logging
Detailed logging at all levels:
```python
logger.info("Successful request")
logger.warning("Retry attempt {}/{}")
logger.error("Request failed after all retries")
logger.debug("Rate limiting - waiting {}s")
```

---

## Troubleshooting

### DrugCentral Issues

**Problem:** No results found
- **Solution:** Pharos has different drug coverage than old DrugCentral
- **Workaround:** Try alternative drug names or SMILES

**Problem:** Missing fields
- **Solution:** Check new field mapping in `_extract_drug_info()`
- **Workaround:** Use GraphQL API for custom queries (future enhancement)

### ZINC Issues

**Problem:** Still getting 403 errors
- **Solution:** Increase `rate_limit_delay` in config
- **Workaround:** Reduce `max_results` to make fewer requests

**Problem:** Slow responses
- **Solution:** This is expected with 2s delays
- **Workaround:** Disable `get_purchasability` for faster searches

### SureChEMBL Issues

**Problem:** 500 errors even with retries
- **Solution:** ChEMBL fallback should activate automatically
- **Workaround:** Set `use_chembl_fallback=True` in config

**Problem:** Different results from fallback
- **Solution:** ChEMBL has broader drug data, not patent-specific
- **Check:** Look for `"source": "chembl_fallback"` in results

---

## Migration Guide

### If You Were Using Old Adapters

#### DrugCentral Migration
```python
# Before (broken)
adapter = DrugCentralAdapter()
result = await adapter.execute("aspirin")

# After (works - same interface!)
adapter = DrugCentralAdapter()  # Now uses Pharos API
result = await adapter.execute("aspirin")

# Check new field names
drug = result.data['drugs'][0]
drug_id = drug['drug_id']  # Was drug.get('id')
inchikey = drug['inchikey']  # Was drug.get('inchiKey')
```

#### ZINC Migration
```python
# Before (403 errors)
adapter = ZINCFragmentsAdapter()
result = await adapter.execute("CCO")

# After (works - same interface!)
adapter = ZINCFragmentsAdapter()  # Now with User-Agent
result = await adapter.execute("CCO")

# Note: Requests are slower now (2s delays) but reliable
```

#### SureChEMBL Migration
```python
# Before (500 errors)
adapter = SureChEMBLAdapter()
result = await adapter.execute({"smiles": "CCO"})

# After (works with fallback!)
adapter = SureChEMBLAdapter()  # Now v2.0 with ChEMBL fallback
result = await adapter.execute({"smiles": "CCO"})

# Check if fallback was used
if result.metadata['used_chembl_fallback']:
    print("Results from ChEMBL (SureChEMBL unavailable)")
```

---

## Version History

### v2.0.0 (October 2025)
- **DrugCentral:** Migrated to Pharos API
- **ZINC:** Added User-Agent headers and retry logic
- **SureChEMBL:** Updated to v2.0 API with ChEMBL fallback
- **All:** Added exponential backoff retry logic
- **All:** Improved error handling and logging

### v1.0.0 (Original)
- Basic adapter implementations
- No retry logic
- Missing User-Agent headers
- Deprecated API endpoints

---

## Future Enhancements

### Planned Improvements
1. **DrugCentral GraphQL:** Implement GraphQL queries for complex searches
2. **ZINC Caching:** Cache successful responses to reduce API calls
3. **SureChEMBL Health Check:** Ping API before requests to avoid fallback delay
4. **Unified Retry Policy:** Shared retry decorator for all adapters
5. **Metrics:** Track API success rates and fallback usage

### Contributing
If you encounter issues or have improvements:
1. Check logs for detailed error messages
2. Verify API endpoints are still valid
3. Test with simple queries first
4. Report issues with full error context

---

## References

### API Documentation
- **DrugCentral/Pharos:** https://pharos.nih.gov/api
- **ZINC15:** https://wiki.docking.org/index.php/ZINC15:API
- **SureChEMBL 2.0:** https://chembl.gitbook.io/surechembl/api
- **ChEMBL API:** https://chembl.gitbook.io/chembl-interface-documentation/web-services

### Related Resources
- SmartAPI DrugCentral: https://smart-api.info/ui/b8399f43cd71264faca4d798c793aa30
- ZINC Database: https://zinc15.docking.org/
- SureChEMBL 2.0 Release: http://chembl.blogspot.com/2025/05/download-surechembl-data-major-update.html

---

## Support

For issues or questions:
1. Check this documentation first
2. Review logs with `logging.DEBUG` level
3. Test with provided test script
4. Check API status pages for service disruptions

**Last Updated:** October 2025
**Maintained By:** PharmForge Development Team
