# Compound Testing - Quick Start Guide

## What is Compound Testing?

Compound Testing is PharmForge's fastest way to analyze a single molecule. Enter a SMILES string, select adapters, and get results in seconds.

## 5-Minute Tutorial

### Step 1: Access the Page

**Option A:** Sidebar Navigation
```
Click "🧪 Compound Testing" in the left sidebar
```

**Option B:** Home Page
```
Click "Test Compound" button on home page
```

### Step 2: Enter a Compound

**Try the Example:**
```
1. Click "📋 Aspirin Example" button
2. Click "✓ Validate SMILES" button
3. See 2D structure appear below
```

**Or Enter Your Own:**
```
1. Paste SMILES in the text box
2. Click "✓ Validate SMILES"
3. If valid, molecule appears
```

### Step 3: Select Adapters

**Quick Selection:**
```
1. Expand a category (e.g., "Docking")
2. Click "Select All" to choose all adapters in category
3. Or check individual adapter boxes
```

**Recommended First Test:**
```
✓ TDC ADMET (fast, gives drug properties)
✓ ChEMBL Similarity (checks novelty)
✓ ZINC Database (checks availability)
```

### Step 4: Run Analysis

```
1. Click "🚀 Run Selected Adapters" button
2. Watch progress bar
3. Wait for completion (usually 5-30 seconds)
```

### Step 5: View Results

**Summary Tab:**
```
- See success/failure counts
- Check total execution time
- View quick results table
```

**Detailed Results:**
```
- Expand any adapter card
- See full JSON output
- Check execution time
```

**Export Data:**
```
- Click "📥 Download CSV" for spreadsheet
- Or "📥 Download JSON" for raw data
```

## Example SMILES Strings

Copy and paste these for quick testing:

| Drug | SMILES | Use Case |
|------|--------|----------|
| **Aspirin** | `CC(=O)Oc1ccccc1C(=O)O` | Basic testing |
| **Caffeine** | `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` | Small molecule |
| **Ibuprofen** | `CC(C)Cc1ccc(cc1)C(C)C(=O)O` | NSAID drug |
| **Metformin** | `CN(C)C(=N)NC(=N)N` | Diabetes drug |
| **Sildenafil** | `CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12` | Complex structure |

## Common Workflows

### 1. Basic Validation
**Goal:** Check if molecule is valid
```
Input → Validate → View structure → Done
Time: ~10 seconds
```

### 2. Drug-likeness Check
**Goal:** Is this compound drug-like?
```
Input → Validate → Select "ADMET" adapters → Run → Check scores
Time: ~30 seconds
```

### 3. Novelty Assessment
**Goal:** Is this compound novel?
```
Input → Validate → Select "ChEMBL Similarity" → Run → Check similarity %
Time: ~15 seconds
```

### 4. Full Analysis
**Goal:** Complete compound profile
```
Input → Validate → Select All adapters → Run → Export results
Time: ~60 seconds
```

## Tips & Tricks

### Speed Up Testing
- Start with fast adapters (ADMET, Similarity)
- Save slow adapters (Docking, Retrosynthesis) for later
- Use "Select All" by category

### Get Better Results
- Validate SMILES before running
- Check adapter status (✅ = reliable)
- Read adapter descriptions
- Export results for record-keeping

### Troubleshooting
- **Invalid SMILES**: Check for typos, use examples first
- **No adapters**: Backend may be offline, mock data will load
- **Slow execution**: Some adapters take 10-30s, be patient
- **Error results**: Check adapter status, try different adapters

## Advanced Features

### Custom Parameters
```
1. Select adapters
2. Expand "⚙️ Advanced: Customize Parameters"
3. Adjust values (e.g., exhaustiveness, thresholds)
4. Run analysis
```

### Batch Testing (Coming Soon)
```
Will support:
- CSV upload with multiple SMILES
- Parallel execution
- Comparison tables
```

## Result Interpretation

### ADMET Adapters
```json
{
  "caco2_permeability": 0.85,  // Higher = better absorption
  "herg_inhibition": 0.15,     // Lower = safer (cardiac)
  "cyp3a4_inhibition": 0.35,   // Lower = fewer interactions
  "lipinski_violations": 0      // 0 = drug-like
}
```

### Docking Adapters
```json
{
  "binding_affinity": -8.5,  // More negative = stronger binding
  "rmsd": 1.2,               // Lower = better prediction
  "binding_mode": "Active site"
}
```

### Similarity Adapters
```json
{
  "max_similarity": 0.78,    // 0-1 scale, higher = less novel
  "novelty_score": 0.22,     // 1 - similarity
  "similar_compounds": 23    // Count in database
}
```

## Next Steps

After testing a compound:

1. **Looks promising?** → Create a full pipeline run with "🚀 New Run"
2. **Want similar compounds?** → Use similarity search in batch mode
3. **Need more analysis?** → Select additional adapters
4. **Ready to synthesize?** → Check retrosynthesis results

## Support

- **Documentation**: See `COMPOUND_TESTING_README.md`
- **Examples**: Use built-in example buttons
- **API Docs**: Check `/api/docs` endpoint
- **Issues**: Report bugs via GitHub

## Performance Benchmarks

| Adapter Type | Typical Time | CPU/Memory |
|-------------|--------------|------------|
| ADMET | 2-5 seconds | Low |
| Similarity | 3-8 seconds | Low |
| Database | 5-10 seconds | Low |
| Docking | 15-30 seconds | High |
| Retrosynthesis | 20-60 seconds | High |

**Total for all adapters:** ~2-3 minutes

## Keyboard Shortcuts

- `Tab`: Navigate between fields
- `Enter`: Submit current field
- `Ctrl+Z`: Undo text input
- `Ctrl+V`: Paste SMILES

## Mobile Support

The page is responsive and works on mobile devices:
- Portrait mode: Stacked layout
- Landscape mode: Side-by-side layout
- Touch-friendly buttons
- Swipe to expand/collapse sections

---

**Ready to test your first compound?** Click "🧪 Compound Testing" now!
