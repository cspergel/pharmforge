# ImmuneBuilder Adapter for PharmForge

Fast antibody structure prediction from sequence using ImmuneBuilder, a deep learning method from the Oxford Protein Informatics Group.

## Overview

ImmuneBuilder provides rapid antibody structure prediction (~1-5 seconds per structure) using deep learning trained on the SAbDab database. It's significantly faster than AlphaFold while maintaining comparable accuracy for antibody variable domains.

## Features

- **Fast Prediction**: Structure prediction in 1-5 seconds
- **Multiple Formats**: Support for Fab, scFv, VHH nanobodies
- **High Quality**: Accuracy comparable to AlphaFold for antibody domains
- **Confidence Scores**: Per-residue confidence estimates
- **Batch Processing**: Predict multiple antibodies efficiently
- **PDB Output**: Standard PDB format structures

## Installation

```bash
pip install immunebuilder
```

### Dependencies

- Python 3.7+
- PyTorch
- Biopython
- NumPy

## Usage

### Basic Structure Prediction

```python
from adapters.immunebuilder import ImmuneBuilderAdapter

# Initialize adapter
adapter = ImmuneBuilderAdapter()

# Predict Fab structure
result = await adapter.execute({
    "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
    "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
    "name": "my_antibody"
})

# Access structure
pdb_structure = result.data["structure_pdb"]
file_path = result.data["file_path"]
confidence = result.data["confidence_scores"]
```

### VHH Nanobody Prediction

```python
# Predict VHH nanobody structure (heavy chain only)
result = await adapter.execute({
    "heavy_chain": "EVQLVESGGGLVQAGGSLRLSCAASGRTFSTYAMGWFRQAPGKEREFVARITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPFDYWGQGTQVTVSS",
    "antibody_type": "vhh",
    "name": "nanobody_1"
})
```

### Batch Prediction

```python
antibodies = [
    {
        "heavy_chain": "EVQLVES...",
        "light_chain": "DIQMTQS...",
        "name": "ab_1"
    },
    {
        "heavy_chain": "QVQLQES...",
        "light_chain": "DIVMTQS...",
        "name": "ab_2"
    }
]

results = await adapter.batch_predict(antibodies)

# Process results
for i, result in enumerate(results):
    if result.success:
        print(f"Antibody {i+1}: Success")
        print(f"  File: {result.data['file_path']}")
        print(f"  Confidence: {result.data['confidence_scores']['mean']}")
    else:
        print(f"Antibody {i+1}: Failed - {result.error}")
```

### Custom Configuration

```python
adapter = ImmuneBuilderAdapter(config={
    "output_dir": "./my_antibodies",
    "num_models": 1,
    "save_structures": True,
    "use_gpu": True
})
```

## Input Format

### Required Fields

- **heavy_chain**: Heavy chain amino acid sequence (single-letter code)

### Optional Fields

- **light_chain**: Light chain amino acid sequence (required for Fab/scFv)
- **antibody_type**: Specify type ("fab", "scfv", "vhh", or "auto")
- **name**: Identifier for the antibody (default: "antibody")

### Sequence Requirements

- Use single-letter amino acid code (ACDEFGHIKLMNPQRSTVWY)
- Remove any gaps or special characters
- Sequences should include at least the variable domain (VH/VL)
- Typical lengths:
  - VH: ~110-130 residues
  - VL: ~105-115 residues
  - VHH: ~110-130 residues

## Output Format

```python
{
    "name": "my_antibody",
    "antibody_type": "fab",
    "structure_pdb": "<PDB format structure string>",
    "file_path": "/path/to/output/my_antibody_fab.pdb",
    "sequences": {
        "heavy_chain": "EVQLVES...",
        "light_chain": "DIQMTQS...",
        "heavy_chain_length": 128,
        "light_chain_length": 112
    },
    "confidence_scores": {
        "mean": 87.3,
        "median": 89.1,
        "min": 65.2,
        "max": 95.8,
        "stdev": 8.4,
        "num_residues": 240
    },
    "structure_stats": {
        "num_atoms": 1850,
        "num_residues": 240,
        "num_chains": 2,
        "chain_ids": ["H", "L"]
    },
    "prediction_time_seconds": 2.3,
    "method": "ImmuneBuilder (ABodyBuilder2)",
    "reference": "Abanades et al., Communications Biology (2023)"
}
```

## Antibody Types

### Fab (Fragment Antigen-Binding)

- Heavy chain variable domain (VH)
- Light chain variable domain (VL)
- Most common format for therapeutic antibodies

### scFv (Single-Chain Variable Fragment)

- VH and VL connected by a linker peptide
- Typically provided as single sequence or separate chains

### VHH Nanobody

- Heavy chain variable domain only
- Found in camelids (camels, llamas, alpacas)
- Smaller size, higher stability

### Full Antibody

- If constant domains are included in sequences
- VH, VL, CH1, CL domains

## Confidence Scores

Confidence scores are stored in the B-factor column of the PDB file:

- **>90**: Very high confidence
- **70-90**: High confidence
- **50-70**: Medium confidence
- **<50**: Low confidence

Low confidence regions typically correspond to:
- CDR loops (especially CDR-H3)
- Flexible linker regions
- Terminal residues

## Integration with PharmForge

### With SAbDab Adapter

```python
from adapters.sabdab import SAbDabAdapter
from adapters.immunebuilder import ImmuneBuilderAdapter

# Get antibody sequences from SAbDab
sabdab = SAbDabAdapter()
sabdab_result = await sabdab.execute("7BWJ")

# Extract CDR sequences and predict structure
heavy_seq = sabdab_result.data["antibody_data"]["heavy_sequence"]
light_seq = sabdab_result.data["antibody_data"]["light_sequence"]

# Predict structure
immunebuilder = ImmuneBuilderAdapter()
structure_result = await immunebuilder.execute({
    "heavy_chain": heavy_seq,
    "light_chain": light_seq,
    "name": "7BWJ_predicted"
})
```

### In Pipeline Workflows

```python
# Example: Predict structures for antibody library
from backend.core.pipeline import Pipeline

pipeline = Pipeline()
pipeline.add_step("predict_structures", {
    "adapter": "immunebuilder",
    "input": antibody_sequences,
    "params": {
        "save_structures": True,
        "num_models": 1
    }
})

results = await pipeline.execute()
```

## Comparison with AlphaFold

| Feature | ImmuneBuilder | AlphaFold |
|---------|--------------|-----------|
| **Speed** | ~1-5 seconds | ~5-30 minutes |
| **Accuracy (VH/VL)** | High (comparable) | Very High |
| **Accuracy (CDR-H3)** | Good | Excellent |
| **Training Data** | SAbDab (antibodies) | PDB (all proteins) |
| **GPU Required** | Optional | Recommended |
| **Use Case** | High-throughput screening | High-accuracy modeling |

**When to use ImmuneBuilder:**
- High-throughput antibody screening
- Initial structure generation
- CDR grafting/humanization
- Speed is critical

**When to use AlphaFold:**
- Highest accuracy needed
- Complex antibody-antigen modeling
- Full-length antibody structures
- Novel CDR-H3 conformations

## Validation

ImmuneBuilder has been validated against:
- SAbDab test set
- Recent antibody crystal structures
- AlphaFold predictions

Typical performance:
- **Backbone RMSD**: 1-2 Å for VH/VL domains
- **CDR-H3 RMSD**: 2-4 Å (most challenging region)
- **Success rate**: >95% for standard antibodies

## Use Cases

1. **Therapeutic Antibody Design**
   - Rapid structure generation for drug candidates
   - Humanization and optimization
   - Epitope mapping

2. **Antibody Engineering**
   - CDR grafting
   - Affinity maturation
   - Stability optimization

3. **Nanobody Development**
   - VHH structure prediction
   - Library screening
   - Camelid antibody engineering

4. **High-Throughput Screening**
   - Process thousands of sequences
   - Structure-based filtering
   - Virtual screening workflows

5. **Structure-Based Design**
   - Docking studies
   - Interface analysis
   - Epitope prediction

## Troubleshooting

### Import Error

```python
ImportError: No module named 'ImmuneBuilder'
```

**Solution**: Install immunebuilder
```bash
pip install immunebuilder
```

### Invalid Sequence Error

```python
error: "Heavy chain contains invalid amino acid characters"
```

**Solution**:
- Use only standard amino acids (ACDEFGHIKLMNPQRSTVWY)
- Remove gaps, spaces, or special characters
- Check for lowercase letters (should be uppercase)

### Low Confidence Scores

**Common reasons**:
- Very long CDR-H3 loops (>20 residues)
- Unusual sequence features
- Missing/incomplete variable domains

**Solutions**:
- Verify sequence quality
- Check that full VH/VL domains are included
- Consider using AlphaFold for difficult cases

### Slow Predictions

**Solutions**:
- Ensure GPU is available and enabled
- Check PyTorch CUDA installation
- Use batch_predict() for multiple structures
- Consider CPU-only for small-scale work

## Performance Tips

1. **Use Batch Prediction**: Process multiple antibodies together
2. **Enable GPU**: 5-10x speedup with GPU acceleration
3. **Cache Results**: Use PharmForge's built-in caching
4. **Pre-validate Sequences**: Check sequences before prediction
5. **Parallel Processing**: Run multiple adapters concurrently

## API Reference

### ImmuneBuilderAdapter

#### Methods

- **execute(input_data, **params)**: Predict single antibody structure
- **batch_predict(antibodies, **params)**: Predict multiple structures
- **validate_input(input_data)**: Validate sequence input
- **generate_cache_key(input_data, **kwargs)**: Generate cache key
- **get_metadata()**: Get adapter information

#### Configuration

- **output_dir**: Output directory for PDB files
- **num_models**: Number of models to generate (default: 1)
- **save_structures**: Save PDB files to disk (default: True)
- **use_gpu**: Use GPU acceleration (default: True)

## References

1. **ImmuneBuilder Paper**:
   - Abanades et al., "ImmuneBuilder: Deep-Learning models for predicting the structures of immune proteins"
   - Communications Biology (2023)
   - DOI: 10.1038/s42003-023-04927-7

2. **SAbDab Database**:
   - Dunbar et al., "SAbDab: the structural antibody database"
   - Nucleic Acids Research (2014)

3. **GitHub Repository**:
   - https://github.com/oxpig/ImmuneBuilder

## License

ImmuneBuilder is developed by the Oxford Protein Informatics Group (OPIG) and is available under the BSD 3-Clause License.

## Support

- **PharmForge Issues**: Submit issues via GitHub
- **ImmuneBuilder Issues**: https://github.com/oxpig/ImmuneBuilder/issues
- **OPIG Website**: http://opig.stats.ox.ac.uk/
