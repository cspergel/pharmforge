# PDB-REDO Adapter: Required Code Fixes

## Quick Reference for Fixing the Adapter

This document provides the exact code changes needed to fix the PDB-REDO adapter.

---

## Fix 1: Update `_check_structure_exists` method (Lines 234-255)

**Location:** Line 244-246

**Current Code:**
```python
# PDB-REDO organizes files by middle two characters
mid = pdb_id[1:3]
url = f"{self.BASE_URL}/{mid}/{pdb_id}/{pdb_id}_final.pdb"
```

**Fixed Code:**
```python
# PDB-REDO uses /db/ prefix for all entries
url = f"{self.BASE_URL}/db/{pdb_id}/{pdb_id}_final.pdb"
```

---

## Fix 2: Update `_download_structure` method (Lines 257-318)

**Location:** Lines 272-285

**Current Code:**
```python
# PDB-REDO organizes files by middle two characters
mid = pdb_id[1:3]

# Construct filename and URL
if format == "pdb":
    filename = f"{pdb_id}_final.pdb"
elif format == "cif":
    filename = f"{pdb_id}_final.cif"
elif format == "mtz":
    filename = f"{pdb_id}_final.mtz"
else:
    return None, None

url = f"{self.BASE_URL}/{mid}/{pdb_id}/{filename}"
```

**Fixed Code:**
```python
# Construct filename and URL
if format == "pdb":
    filename = f"{pdb_id}_final.pdb"
elif format == "cif":
    filename = f"{pdb_id}_final.cif"
elif format == "mtz":
    filename = f"{pdb_id}_final.mtz"
else:
    return None, None

# PDB-REDO uses /db/ prefix for all entries
url = f"{self.BASE_URL}/db/{pdb_id}/{filename}"
```

---

## Fix 3: Update `_get_quality_metrics` method (Lines 320-345)

**Location:** Lines 330-331

**Current Code:**
```python
mid = pdb_id[1:3]
url = f"{self.BASE_URL}/{mid}/{pdb_id}/{pdb_id}_final.json"
```

**Fixed Code:**
```python
url = f"{self.BASE_URL}/db/{pdb_id}/{pdb_id}_final.json"
```

**⚠️ WARNING:** The JSON format returned by this endpoint may not match the expected format in `_parse_quality_metrics`. Additional investigation needed.

---

## Fix 4: Update `_get_validation_report` method (Lines 384-408)

**Location:** Lines 394-395

**Current Code:**
```python
mid = pdb_id[1:3]
url = f"{self.BASE_URL}/{mid}/{pdb_id}/{pdb_id}_validation.json"
```

**Fixed Code:**
```python
url = f"{self.BASE_URL}/db/{pdb_id}/{pdb_id}_validation.json"
```

---

## Fix 5: Improve Input Validation (Lines 86-112)

**Location:** Lines 96-111

**Current Code:**
```python
if not pdb_id or not isinstance(pdb_id, str):
    return False

# PDB IDs are exactly 4 characters: 1 digit + 3 alphanumeric
pdb_id = pdb_id.strip().lower()  # PDB-REDO uses lowercase
if len(pdb_id) != 4:
    return False
```

**Issue:** The strip() happens after the isinstance check but the length check happens after strip, which is correct. However, the method parameter should be stripped and reassigned.

**Better Approach:**
```python
if not pdb_id or not isinstance(pdb_id, str):
    return False

# PDB IDs are exactly 4 characters: 1 digit + 3 alphanumeric
# Strip whitespace and convert to lowercase (PDB-REDO uses lowercase)
pdb_id = pdb_id.strip()
if len(pdb_id) != 4:
    return False

pdb_id = pdb_id.lower()
```

**Note:** This is a minor issue since the execute() method also strips (line 143), but validation should be consistent.

---

## Summary of Changes

### Files to Modify
- `adapters/pdb_redo/adapter.py`

### Lines to Change
1. **Line 245-246** (remove `mid` calculation)
2. **Line 246** (update URL)
3. **Line 273** (remove `mid` calculation)
4. **Line 285** (update URL)
5. **Line 331** (remove `mid` and update URL)
6. **Line 395** (remove `mid` and update URL)

### Pattern to Find/Replace

**Search for:**
```python
mid = pdb_id[1:3]
```

**Action:** Delete these lines (appears 4 times)

**Search for:**
```python
f"{self.BASE_URL}/{mid}/{pdb_id}/
```

**Replace with:**
```python
f"{self.BASE_URL}/db/{pdb_id}/
```

---

## Testing After Fixes

Run the test script to verify:
```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"
python test_pdb_redo_adapter.py
```

**Expected Results:**
- Test 1 (Initialization): ✓ PASS
- Test 2 (Validation): ✓ PASS
- Test 3 (Existence): ✓ PASS (with correct True/False for structures)
- Test 4 (Full Retrieval): ✓ PASS (with actual file content for 6lu7)
- Test 5 (Multiple Structures): ✓ PASS
- Test 6 (Cache): ✓ PASS (files cached to disk)

---

## Alternative: Use Known Working Structure in Tests

If 1abc doesn't exist in PDB-REDO, update the test script to use 6lu7 instead:

**In `test_pdb_redo_adapter.py`:**

Change line 237:
```python
# OLD
results["full_retrieval_1abc"] = await test_full_retrieval(adapter, "1abc")

# NEW
results["full_retrieval_6lu7"] = await test_full_retrieval(adapter, "6lu7")
```

---

## Verification URLs

Test these URLs after fixes to verify they work:

```bash
# Should return 200 OK
curl -I https://pdb-redo.eu/db/6lu7/6lu7_final.pdb
curl -I https://pdb-redo.eu/db/6lu7/6lu7_final.cif
curl -I https://pdb-redo.eu/db/6lu7/6lu7_final.mtz

# Should return 404 (structure not in PDB-REDO)
curl -I https://pdb-redo.eu/db/1abc/1abc_final.pdb
curl -I https://pdb-redo.eu/db/1crn/1crn_final.pdb
```

---

**Last Updated:** 2025-10-26
**Created By:** Claude Code Agent
