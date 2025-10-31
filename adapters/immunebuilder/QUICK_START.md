# ImmuneBuilder Quick Start Guide

Fast antibody structure prediction in seconds.

## Installation

```bash
pip install immunebuilder
```

## Basic Usage

```python
from adapters.immunebuilder import ImmuneBuilderAdapter

# Initialize
adapter = ImmuneBuilderAdapter()

# Predict Fab structure
result = await adapter.execute({
    "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
    "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
    "name": "my_antibody"
})

# Access results
pdb_file = result.data["file_path"]
confidence = result.data["confidence_scores"]["mean"]
```

## Common Patterns

### VHH Nanobody

```python
result = await adapter.execute({
    "heavy_chain": "EVQLVESGGGLVQAGGSLRLSCAASGRTFSTYAMGWFRQAPGKEREFVARITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPFDYWGQGTQVTVSS",
    "antibody_type": "vhh",
    "name": "nanobody_1"
})
```

### Batch Processing

```python
antibodies = [
    {"heavy_chain": "...", "light_chain": "...", "name": "ab1"},
    {"heavy_chain": "...", "light_chain": "...", "name": "ab2"}
]

results = await adapter.batch_predict(antibodies)
```

### Custom Output Directory

```python
adapter = ImmuneBuilderAdapter(config={
    "output_dir": "./my_structures",
    "save_structures": True
})
```

## Input Requirements

**Required:**
- `heavy_chain`: Heavy chain sequence (VH domain minimum)

**Optional:**
- `light_chain`: Light chain sequence (VL domain)
- `antibody_type`: "fab", "scfv", "vhh", or "auto" (default)
- `name`: Identifier for output files

**Sequence Format:**
- Single-letter amino acid code (ACDEFGHIKLMNPQRSTVWY)
- No gaps or special characters
- Typical VH length: 110-130 residues
- Typical VL length: 105-115 residues

## Output Structure

```python
{
    "name": "my_antibody",
    "antibody_type": "fab",
    "structure_pdb": "<PDB string>",
    "file_path": "/path/to/my_antibody_fab.pdb",
    "sequences": {
        "heavy_chain": "...",
        "light_chain": "...",
        "heavy_chain_length": 128,
        "light_chain_length": 112
    },
    "confidence_scores": {
        "mean": 87.3,
        "median": 89.1,
        "min": 65.2,
        "max": 95.8
    },
    "structure_stats": {
        "num_residues": 240,
        "num_chains": 2,
        "chain_ids": ["H", "L"]
    },
    "prediction_time_seconds": 2.3
}
```

## Key Features

| Feature | Description |
|---------|-------------|
| **Speed** | 1-5 seconds per structure |
| **Accuracy** | Comparable to AlphaFold for VH/VL |
| **GPU** | Optional (5-10x speedup) |
| **Formats** | Fab, scFv, VHH nanobodies |
| **Batch** | Multiple predictions supported |

## Confidence Interpretation

- **>90**: Very high confidence (well-structured regions)
- **70-90**: High confidence (most of structure)
- **50-70**: Medium confidence (some flexibility)
- **<50**: Low confidence (flexible loops, termini)

## Comparison: ImmuneBuilder vs AlphaFold

| Aspect | ImmuneBuilder | AlphaFold |
|--------|--------------|-----------|
| Speed | ~2 seconds | ~10 minutes |
| Training | SAbDab antibodies | All PDB proteins |
| Best for | High-throughput | Highest accuracy |
| GPU needed | Optional | Recommended |

## Troubleshooting

### Import Error
```bash
pip install immunebuilder
```

### Invalid Sequence
- Check for only valid amino acids (ACDEFGHIKLMNPQRSTVWY)
- Remove gaps/spaces
- Ensure minimum length (VH: ~110 aa, VL: ~105 aa)

### Low Confidence
- Verify full variable domains included
- Check CDR-H3 length (>20 residues may be challenging)
- Consider AlphaFold for difficult cases

## Integration Examples

### With SAbDab

```python
# Get sequence from SAbDab, predict structure
from adapters.sabdab import SAbDabAdapter
from adapters.immunebuilder import ImmuneBuilderAdapter

sabdab = SAbDabAdapter()
immunebuilder = ImmuneBuilderAdapter()

# Query SAbDab
sabdab_result = await sabdab.execute("7BWJ")
heavy = sabdab_result.data["antibody_data"]["heavy_sequence"]
light = sabdab_result.data["antibody_data"]["light_sequence"]

# Predict structure
pred = await immunebuilder.execute({
    "heavy_chain": heavy,
    "light_chain": light,
    "name": "7BWJ_predicted"
})
```

### In Pipeline

```python
from backend.core.pipeline import Pipeline

pipeline = Pipeline()
pipeline.add_step("predict", {
    "adapter": "immunebuilder",
    "input": {"heavy_chain": "...", "light_chain": "..."}
})

results = await pipeline.execute()
```

## Performance Tips

1. **Use GPU**: Enable for 5-10x speedup
2. **Batch Processing**: Use `batch_predict()` for multiple structures
3. **Cache Results**: Leverage PharmForge caching
4. **Pre-validate**: Check sequences before prediction

## Use Cases

- Therapeutic antibody design
- High-throughput screening
- Nanobody engineering
- CDR grafting/humanization
- Structure-based optimization
- Virtual screening

## Resources

- **Paper**: Abanades et al., Communications Biology (2023)
- **GitHub**: https://github.com/oxpig/ImmuneBuilder
- **OPIG**: http://opig.stats.ox.ac.uk/

## Next Steps

1. Try the examples: `python example_usage.py`
2. Read full documentation: `README.md`
3. Integrate with your pipeline
4. Explore advanced features

---

**Need Help?**
- Check README.md for detailed documentation
- See example_usage.py for code examples
- Visit GitHub for issues/questions
