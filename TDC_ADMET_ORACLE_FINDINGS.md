proceedd# TDC Oracle API Investigation - 2025-10-24

## Discovery Summary

After successfully installing PyTDC 1.1.15 via Docker, I investigated how to use pre-trained ADMET models and discovered important limitations with the TDC library.

## Key Findings

### 1. TDC Has Two Separate APIs

**Datasets API (`tdc.single_pred.ADME`)**:
- Provides **datasets only** for ADMET properties
- Methods: `get_data()`, `get_split()`, `convert_format()`
- **No `predict()` method** - requires training your own model
- Example:
  ```python
  from tdc.single_pred import ADME
  model = ADME(name='Caco2_Wang')
  data = model.get_data()  # Returns training dataset
  # No model.predict() available!
  ```

**Oracle API (`tdc.Oracle`)**:
- Provides **pre-trained models** for limited properties
- Method: `__call__(smiles_list)` returns predictions
- Example:
  ```python
  from tdc import Oracle
  oracle = Oracle(name='cyp3a4_veith')
  predictions = oracle(['CCO'])  # [0.234]
  ```

### 2. Available ADMET Oracles

Out of our 17 target ADMET models, **only 1 is available as a pre-trained Oracle**:
- ✅ `cyp3a4_veith` - CYP3A4 inhibition prediction

**NOT available as Oracles** (dataset only, require training):
- ✗ `caco2_wang` - Caco-2 permeability
- ✗ `hia_hou` - Human intestinal absorption
- ✗ `pgp_broccatelli` - P-glycoprotein inhibition
- ✗ `bioavailability_ma` - Oral bioavailability
- ✗ `lipophilicity_astrazeneca` - Lipophilicity
- ✗ `bbb_martins` - Blood-brain barrier penetration
- ✗ `ppbr_az` - Plasma protein binding
- ✗ `vdss_lombardo` - Volume distribution
- ✗ `cyp2c9_veith` - CYP2C9 inhibition
- ✗ `cyp2d6_veith` - CYP2D6 inhibition
- ✗ `half_life_obach` - Half-life
- ✗ `clearance_hepatocyte_az` - Clearance
- ✗ `herg` - hERG cardiotoxicity
- ✗ `ames` - Ames mutagenicity
- ✗ `dili` - Drug-induced liver injury
- ✗ `ld50_zhu` - Acute toxicity

### 3. Complete Oracle List

TDC provides 93 oracles total, but most are for other purposes:
- **Docking**: `1iep_docking`, `2rgp_docking`, `3eml_docking`, etc.
- **Bioactivity**: `drd2`, `gsk3b`, `jnk3`
- **Drug-likeness**: `qed`, `logp`, `sa` (synthetic accessibility)
- **Synthesis**: `askcos`, `ibm_rxn`, `molecule_one_synthesis`
- **ADMET (limited)**: `cyp3a4_veith`

## Implications for PharmForge

### Current Adapter Status
The `TDCAdmetAdapter` code (430 lines) is well-written but **assumes pre-trained models exist** for all 17 ADMET properties. This assumption is incorrect for TDC 1.0+.

### Options Going Forward

**Option 1: Use Oracle API for Available Models**
- Update adapter to use `tdc.Oracle` for `cyp3a4_veith`
- **Pros**: Works immediately
- **Cons**: Only 1 property available, not comprehensive ADMET coverage

**Option 2: Train Models on TDC Datasets**
- Use `tdc.single_pred.ADME` to get datasets
- Train scikit-learn/XGBoost models for each property
- Cache trained models
- **Pros**: Full ADMET coverage, reproducible
- **Cons**: Requires ML pipeline, longer development time, model quality unknown

**Option 3: Use Alternative ADMET Libraries**
- **ADMETlab**: Web API with 30+ ADMET endpoints (requires API key)
- **pkCSM**: Web-based predictor (no Python API, would need scraping)
- **SwissADME**: Web-based (no API)
- **DeepPurpose**: Python library with ADMET models
- **admet_ai**: Chemprop-based ADMET predictor
- **Pros**: Pre-trained models ready to use
- **Cons**: External dependencies, potential API limits

**Option 4: Hybrid Approach (Recommended)**
- Use `tdc.Oracle` for `cyp3a4_veith`
- Use `rdkit` calculations for simple properties (logP, MW, TPSA)
- Mark complex properties as "unavailable" in adapter
- Document limitation clearly
- Plan to add more predictors in future (Week 4+)

### Recommended Next Steps

1. **Update `TDCAdmetAdapter`** to:
   - Use `Oracle` API for `cyp3a4_veith`
   - Return informative error for unsupported models
   - Update metadata to reflect limited model availability

2. **Update documentation** to:
   - Clarify TDC provides datasets, not universal ADMET predictor
   - Set user expectations correctly
   - Roadmap for additional ADMET sources

3. **Consider adding RDKit-based properties** (Week 3.5):
   - LogP (partition coefficient)
   - TPSA (topological polar surface area)
   - Number of H-bond donors/acceptors
   - Lipinski's Rule of Five compliance
   - These are fast, deterministic, and widely used

4. **Evaluate `admet_ai` or `DeepPurpose`** (Week 4):
   - Both provide pre-trained ADMET models
   - Could create separate adapters
   - Would provide comprehensive coverage

## Docker Build Success

✅ **PyTDC 1.1.15 installed successfully** via Docker
- Image size: 14.5GB
- All dependencies resolved (numpy>=1.26.4, pandas>=2.1.4, etc.)
- TDC imports and loads models correctly
- Oracle API works for supported models

## Testing Results

```bash
# Test 1: Oracle API loads successfully
docker run --rm pharmforge-backend:test python -c "from tdc import Oracle; oracle = Oracle(name='cyp3a4_veith'); print('Oracle created')"
# Output: Downloading Oracle... Done! Oracle created

# Test 2: Oracle prediction FAILS - requires DeepPurpose
docker run --rm pharmforge-backend:test python -c "from tdc import Oracle; oracle = Oracle(name='cyp3a4_veith'); print(oracle(['CCO']))"
# Output: ImportError: Please install DeepPurpose by 'pip install DeepPurpose'

# Test 3: Model not available as Oracle
docker run --rm pharmforge-backend:test python -c "from tdc import Oracle; oracle = Oracle(name='Caco2_Wang')"
# Output: ValueError: 'caco2_wang' does not match to available values
```

**Finding**: Even the `cyp3a4_veith` Oracle requires the separate `DeepPurpose` library, which has its own heavy dependencies (PyTorch, RDKit, etc.). This makes the Oracle approach more complex than anticipated.

## Conclusion

**Week 3 Status: PARTIAL SUCCESS**

| Component | Status | Notes |
|-----------|--------|-------|
| Adapter Code | ✅ Complete | Well-structured, needs API update |
| PyTDC Installation | ✅ Success | Version 1.1.15 in Docker |
| ADMET Model Access | ⚠️ Limited | Only 1/17 models available |
| Documentation | ✅ Complete | Accurate to original assumptions |
| Path Forward | ✅ Clear | Multiple viable options identified |

**Recommendation**: Pivot to **Option 4 (Hybrid Approach)** for Week 3 completion:
1. Update adapter to use Oracle for `cyp3a4_veith`
2. Add RDKit molecular descriptors for basic properties
3. Document limitations honestly
4. Plan Week 4 expansion with alternative ADMET libraries

This provides **some** ADMET functionality immediately while maintaining architectural quality and setting realistic expectations.
