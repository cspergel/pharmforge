# Download fpscores.pkl.gz for SAScore

**Required for:** `synthesis/sascore.py` (Synthetic Accessibility Score)

---

## Option 1: Direct Download (Recommended)

```bash
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SA_Score/fpscores.pkl.gz \
  -O adapters/custom/resources/fpscores.pkl.gz
```

Or with curl:
```bash
curl -L https://github.com/rdkit/rdkit/raw/master/Contrib/SA_Score/fpscores.pkl.gz \
  -o adapters/custom/resources/fpscores.pkl.gz
```

---

## Option 2: Alternative Source

```bash
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SAscore/data/fpscores.pkl.gz \
  -O adapters/custom/resources/fpscores.pkl.gz
```

---

## Verify Download

```bash
ls -lh adapters/custom/resources/fpscores.pkl.gz
```

**Expected output:**
```
-rw-r--r-- 1 user user 1.0M fpscores.pkl.gz
```

**File size should be ~1 MB**

---

## Test SAScore

```bash
python -c "
from adapters.custom.synthesis import SynthesisAccessibilityEvaluator

evaluator = SynthesisAccessibilityEvaluator()
result = evaluator.evaluate('CC(=O)Oc1ccccc1C(=O)O')  # Aspirin
print(f'SAScore: {result[\"score\"]} (1=easy, 10=hard)')
print(f'Difficulty: {result[\"summary\"][\"synthesis_difficulty\"]}')
"
```

**Expected output:**
```
[SAScore] Loaded 85479 fragment scores.
SAScore: 2.3 (1=easy, 10=hard)
Difficulty: low
```

---

## Troubleshooting

### Download fails
- Try the alternative source (Option 2)
- Manually download from browser: https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score
- Right-click fpscores.pkl.gz → Save As → place in `adapters/custom/resources/`

### File size wrong
- Should be ~1 MB compressed
- If much smaller, re-download (may have downloaded HTML error page)

### Permission issues
```bash
chmod 644 adapters/custom/resources/fpscores.pkl.gz
```

---

## What's in this file?

Fragment frequency data from millions of known compounds:
- 85,479 molecular fragments (Morgan fingerprints)
- Each fragment has a score based on occurrence frequency
- Used to estimate how "common" or "easy to make" a molecule is

**Reference:** Ertl & Schuffenhauer, J. Cheminform. 2009

---

## Quick Start Script

Save this as `download_fpscores.sh`:

```bash
#!/bin/bash
cd "$(dirname "$0")"
echo "Downloading fpscores.pkl.gz..."
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SA_Score/fpscores.pkl.gz || \
wget https://github.com/rdkit/rdkit/raw/master/Contrib/SAscore/data/fpscores.pkl.gz
echo "Done! File saved to: adapters/custom/resources/fpscores.pkl.gz"
ls -lh fpscores.pkl.gz
```

Run:
```bash
chmod +x adapters/custom/resources/download_fpscores.sh
./adapters/custom/resources/download_fpscores.sh
```
