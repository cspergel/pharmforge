# OpenMM Adapter for PharmForge

A high-performance molecular dynamics adapter for PharmForge that provides energy minimization, molecular dynamics simulation, and property calculation for small molecules using the OpenMM toolkit.

## Overview

OpenMM is a powerful molecular simulation toolkit that can leverage CPUs and GPUs (CUDA/OpenCL) for high-performance simulations. This adapter integrates OpenMM into PharmForge's workflow for:

- **Molecular stability assessment** - Validate synthesized molecules before lab work
- **Conformer generation** - Find low-energy conformations
- **Property prediction** - Calculate RMSD, radius of gyration, energies
- **Drug candidate screening** - Compare stability of multiple candidates
- **Post-synthesis validation** - Verify molecular properties after retrosynthesis planning

## What is OpenMM?

OpenMM is a toolkit for molecular simulation with the following features:

- **High Performance**: Optimized for modern CPUs and GPUs
- **GPU Acceleration**: 5-50x faster with CUDA/OpenCL support
- **Force Fields**: Support for AMBER, GAFF, CHARMM, and custom force fields
- **Python API**: Easy integration with scientific Python ecosystem

**Reference**: Eastman et al., *PLoS Comput. Biol.* 2017, 13(7), e1005659
**Website**: http://openmm.org/

## When to Use This Adapter

### Ideal Use Cases

1. **Pre-Synthesis Validation**
   - Validate molecule stability before expensive synthesis
   - Predict which synthesis targets are most feasible
   - Identify potential structural issues

2. **Drug Candidate Screening**
   - Compare stability of multiple drug candidates
   - Rank molecules by predicted stability
   - Filter unstable or infeasible molecules

3. **Conformer Analysis**
   - Generate low-energy conformations
   - Identify most stable molecular geometry
   - Compare conformational flexibility

4. **Integration with Retrosynthesis**
   - Use after AiZynthFinder or LLM retrosynthesis
   - Validate that synthesis targets are stable
   - Prioritize synthesis based on stability scores

### When NOT to Use

- **Large proteins** - OpenMM excels at small molecules; for proteins, consider specialized tools
- **Ultra-fast screening** - MD simulations take seconds to minutes; use faster descriptors for initial filtering
- **Batch processing without GPU** - Can be slow on CPU; consider GPU acceleration or caching

## Installation

### Recommended: Conda Installation (includes GPU support)

```bash
conda install -c conda-forge openmm
conda install -c conda-forge pdbfixer
```

### Alternative: Pip Installation (CPU only)

```bash
pip install openmm>=8.0.0
pip install pdbfixer>=1.9
```

### GPU Acceleration

For **CUDA** (NVIDIA GPUs):
- Install [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
- OpenMM will automatically detect and use CUDA
- Performance: 5-50x faster than CPU

For **OpenCL** (AMD/Intel GPUs):
- Install OpenCL drivers for your GPU
- OpenMM supports OpenCL for non-NVIDIA GPUs

### Verify Installation

```python
import openmm
print(f"OpenMM version: {openmm.version.version}")

# Check available platforms
from openmm import Platform
for i in range(Platform.getNumPlatforms()):
    platform = Platform.getPlatform(i)
    print(f"Platform {i}: {platform.getName()}")
```

## Configuration Options

The adapter supports extensive configuration:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `force_field` | str | `"gaff"` | Force field to use (gaff, amber14, amber99) |
| `minimize_steps` | int | `1000` | Energy minimization steps |
| `minimize_tolerance` | float | `10.0` | Convergence tolerance (kJ/mol) |
| `run_md` | bool | `False` | Whether to run MD simulation |
| `md_steps` | int | `10000` | MD simulation steps |
| `md_temperature` | float | `300.0` | Temperature in Kelvin |
| `md_timestep` | float | `2.0` | Timestep in femtoseconds |
| `save_trajectory` | bool | `False` | Save trajectory coordinates |
| `platform` | str | `None` | OpenMM platform (CPU, CUDA, OpenCL, Reference) |
| `friction_coefficient` | float | `1.0` | Langevin friction (1/ps) |
| `report_interval` | int | `100` | Report every N steps |

### Example Configuration

```python
from adapters.openmm.adapter import OpenMMAdapter

# Quick stability check (default - fast)
adapter_fast = OpenMMAdapter(config={
    "minimize_steps": 1000,
    "run_md": False
})

# Full MD simulation (comprehensive - slow)
adapter_full = OpenMMAdapter(config={
    "force_field": "gaff",
    "minimize_steps": 2000,
    "run_md": True,
    "md_steps": 50000,  # 100 ps at 2 fs timestep
    "md_temperature": 300.0,
    "platform": "CUDA"  # Use GPU
})

# Custom parameters
adapter_custom = OpenMMAdapter(config={
    "minimize_steps": 1500,
    "run_md": True,
    "md_steps": 20000,
    "md_temperature": 310.0,  # Body temperature
    "save_trajectory": True
})
```

## Usage Examples

### Basic Usage

```python
from adapters.openmm.adapter import OpenMMAdapter
import asyncio

async def main():
    # Initialize adapter
    adapter = OpenMMAdapter()

    # Run energy minimization (fast validation)
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        print(f"Energy: {data['minimization']['final_energy']:.2f} kJ/mol")
        print(f"Stability: {data['stability_score']:.3f}")
        print(f"Feasibility: {data['feasibility']}")

asyncio.run(main())
```

### With Molecular Dynamics

```python
async def analyze_with_md():
    adapter = OpenMMAdapter(config={
        "run_md": True,
        "md_steps": 20000
    })

    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        md = data['molecular_dynamics']

        print(f"MD Runtime: {md['time']:.2f} ps")
        print(f"RMSD: {md['rmsd']:.3f} Å")
        print(f"Radius of gyration: {md['radius_of_gyration']:.3f} Å")
        print(f"Average energy: {md['average_energy']:.2f} kJ/mol")

asyncio.run(analyze_with_md())
```

### Integration with Retrosynthesis

```python
from adapters.aizynthfinder.adapter import AiZynthFinderAdapter
from adapters.openmm.adapter import OpenMMAdapter

async def validate_synthesis_target():
    # Plan synthesis
    retro_adapter = AiZynthFinderAdapter()
    retro_result = await retro_adapter.execute("target_smiles")

    if retro_result.success and retro_result.data['routes_found'] > 0:
        # Validate stability before synthesis
        md_adapter = OpenMMAdapter()
        stability_result = await md_adapter.execute("target_smiles")

        if stability_result.success:
            synth_score = retro_result.data['synthesis_score']
            stability_score = stability_result.data['stability_score']

            # Combined feasibility
            combined_score = 0.6 * synth_score + 0.4 * stability_score

            print(f"Synthesis feasibility: {synth_score:.3f}")
            print(f"Stability score: {stability_score:.3f}")
            print(f"Combined score: {combined_score:.3f}")

            if combined_score > 0.7:
                print("✓ Recommended for synthesis")
            else:
                print("⚠ Consider alternative targets")

asyncio.run(validate_synthesis_target())
```

### Comparing Multiple Candidates

```python
async def compare_candidates():
    adapter = OpenMMAdapter()

    candidates = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O"
    }

    results = []
    for name, smiles in candidates.items():
        result = await adapter.execute(smiles)
        if result.success:
            results.append({
                "name": name,
                "stability": result.data['stability_score'],
                "energy": result.data['minimization']['final_energy']
            })

    # Sort by stability
    results.sort(key=lambda x: x['stability'], reverse=True)

    print("Stability Ranking:")
    for i, r in enumerate(results, 1):
        print(f"{i}. {r['name']}: {r['stability']:.3f} ({r['energy']:.2f} kJ/mol)")

asyncio.run(compare_candidates())
```

## Output Data Structure

### Energy Minimization Only

```python
{
    "smiles": str,  # Input SMILES
    "minimization": {
        "initial_energy": float,  # kJ/mol
        "final_energy": float,    # kJ/mol
        "energy_change": float,   # kJ/mol
        "converged": bool,
        "steps": int
    },
    "structure": {
        "num_atoms": int,
        "molecular_weight": float,  # g/mol
        "pdb_string": str  # 3D coordinates in PDB format
    },
    "stability_score": float,  # 0-1 (higher is more stable)
    "feasibility": str,  # "high", "medium", "low"
    "model": str,  # "OpenMM X.X.X"
    "force_field": str,
    "reference": str
}
```

### With Molecular Dynamics

Additional field added:

```python
{
    ...,  # All fields from minimization
    "molecular_dynamics": {
        "temperature": float,  # Kelvin
        "steps": int,
        "time": float,  # picoseconds
        "timestep": float,  # femtoseconds
        "final_energy": float,  # kJ/mol
        "average_energy": float,  # kJ/mol
        "energy_std": float,  # kJ/mol (standard deviation)
        "rmsd": float,  # Ångströms (vs. initial structure)
        "radius_of_gyration": float,  # Ångströms
        "energy_trajectory": List[float],  # If save_trajectory=True
        "time_trajectory": List[float]     # If save_trajectory=True
    }
}
```

## Interpreting Results

### Stability Score (0-1)

- **0.7-1.0** (High): Molecule is stable and feasible
  - Low energy configuration
  - Low RMSD during MD
  - Recommended for synthesis

- **0.4-0.7** (Medium): Moderate stability
  - Acceptable energy
  - Some structural flexibility
  - Consider for synthesis with caution

- **0.0-0.4** (Low): Poor stability
  - High energy or large conformational changes
  - May be unstable or difficult to synthesize
  - Consider alternative designs

### Energy Values

Energy ranges for small drug-like molecules (approximate):

- **< -100 kJ/mol**: Very stable
- **-100 to +100 kJ/mol**: Moderate stability
- **> +100 kJ/mol**: Potentially unstable

**Note**: Absolute energy values depend on force field, molecular size, and charge state. Always compare relative energies within the same system.

### RMSD (Root Mean Square Deviation)

Measures structural change during MD simulation:

- **< 1 Å**: Very stable structure
- **1-3 Å**: Moderate flexibility (typical for small molecules)
- **> 3 Å**: High flexibility or instability

### Radius of Gyration (Rg)

Indicates molecular compactness:

- Smaller Rg → More compact
- Larger Rg → More extended
- Changes over time indicate conformational transitions

## Computational Costs

### Performance Benchmarks

Approximate times for a typical drug-like molecule (~30 atoms):

| Operation | CPU Time | GPU Time (CUDA) | Speedup |
|-----------|----------|-----------------|---------|
| Energy minimization (1000 steps) | 2-5 sec | 0.2-0.5 sec | 10x |
| MD simulation (10 ps, 5000 steps) | 10-30 sec | 1-3 sec | 10-30x |
| MD simulation (100 ps, 50000 steps) | 100-300 sec | 10-30 sec | 10-30x |

**Factors affecting performance:**
- Molecular size (more atoms = slower)
- Force field complexity
- Hardware (CPU model, GPU model)
- Platform (CPU < OpenCL < CUDA)

### Recommendations

1. **For quick screening**: Use minimization only (`run_md=False`)
2. **For batch processing**: Use GPU acceleration if available
3. **For detailed analysis**: Run MD on filtered candidates only
4. **For production**: Cache results (adapter handles this automatically)

## System Requirements

### Minimum Requirements

- **CPU**: Any modern multi-core processor
- **RAM**: 2 GB (for small molecules)
- **Python**: 3.8 or higher
- **OS**: Linux, macOS, or Windows

### Recommended for Production

- **GPU**: NVIDIA GPU with CUDA support (GTX 1060 or better)
- **RAM**: 8 GB
- **Storage**: SSD for faster I/O
- **OS**: Linux (best performance and compatibility)

### GPU Support

| Platform | Hardware | Performance |
|----------|----------|-------------|
| **CUDA** | NVIDIA GPUs | Best (10-50x faster) |
| **OpenCL** | AMD, Intel GPUs | Good (5-20x faster) |
| **CPU** | All CPUs | Baseline |

## Troubleshooting

### OpenMM Not Found

```python
ImportError: No module named 'openmm'
```

**Solution**: Install OpenMM
```bash
conda install -c conda-forge openmm
# or
pip install openmm
```

### Force Field Error

```
Exception: Could not load force field 'gaff'
```

**Solution**: Use simple force field or install openmmforcefields
```bash
pip install openmmforcefields
```

Or use config:
```python
adapter = OpenMMAdapter(config={"force_field": "amber14"})
```

### GPU Not Detected

```
Available platforms: CPU, Reference (no CUDA/OpenCL)
```

**Solution**: Check GPU drivers
- CUDA: Install [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit)
- OpenCL: Install GPU vendor's OpenCL runtime

### Simulation Crashes

```
Exception during minimization/MD
```

**Possible causes**:
1. Invalid molecule structure → Check SMILES validity
2. Extreme initial energy → Increase minimize_steps
3. Incompatible force field → Try different force field
4. Memory limit → Reduce md_steps or use smaller molecules

**Debug mode**:
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## Advanced Usage

### Custom Platform Selection

```python
# Force CPU (even if GPU available)
adapter = OpenMMAdapter(config={"platform": "CPU"})

# Force CUDA
adapter = OpenMMAdapter(config={"platform": "CUDA"})

# Auto-select fastest
adapter = OpenMMAdapter(config={"platform": None})
```

### Temperature-Dependent Simulations

```python
# Room temperature
adapter_room = OpenMMAdapter(config={
    "run_md": True,
    "md_temperature": 298.0
})

# Body temperature
adapter_body = OpenMMAdapter(config={
    "run_md": True,
    "md_temperature": 310.0
})

# High temperature (conformational sampling)
adapter_hot = OpenMMAdapter(config={
    "run_md": True,
    "md_temperature": 400.0
})
```

### Trajectory Analysis

```python
adapter = OpenMMAdapter(config={
    "run_md": True,
    "md_steps": 50000,
    "save_trajectory": True
})

result = await adapter.execute(smiles)

if result.success:
    energies = result.data['molecular_dynamics']['energy_trajectory']
    times = result.data['molecular_dynamics']['time_trajectory']

    # Plot energy vs. time
    import matplotlib.pyplot as plt
    plt.plot(times, energies)
    plt.xlabel('Time (ps)')
    plt.ylabel('Energy (kJ/mol)')
    plt.show()
```

## Integration with PharmForge Workflow

### Recommended Workflow

1. **Target Selection** → Choose drug candidates
2. **Retrosynthesis** → Plan synthesis routes (AiZynthFinder/LLM)
3. **Stability Validation** → Run OpenMM minimization (this adapter)
4. **Ranking** → Combine synthesis + stability scores
5. **Detailed Analysis** → Run MD on top candidates
6. **Synthesis** → Proceed with validated, stable molecules

### Example Pipeline

```python
async def complete_drug_design_pipeline(target_smiles: str):
    # Step 1: Retrosynthesis
    retro = AiZynthFinderAdapter()
    synth_result = await retro.execute(target_smiles)

    if not synth_result.success or synth_result.data['routes_found'] == 0:
        print("No synthesis routes found")
        return

    # Step 2: Quick stability check
    md = OpenMMAdapter(config={"minimize_steps": 1000})
    stability_result = await md.execute(target_smiles)

    if not stability_result.success:
        print("Stability validation failed")
        return

    # Step 3: Combined scoring
    synth_score = synth_result.data['synthesis_score']
    stability_score = stability_result.data['stability_score']
    combined = 0.6 * synth_score + 0.4 * stability_score

    print(f"Synthesis: {synth_score:.3f}")
    print(f"Stability: {stability_score:.3f}")
    print(f"Combined: {combined:.3f}")

    # Step 4: Detailed MD if promising
    if combined > 0.7:
        md_detailed = OpenMMAdapter(config={
            "run_md": True,
            "md_steps": 50000
        })
        md_result = await md_detailed.execute(target_smiles)

        if md_result.success:
            md_data = md_result.data['molecular_dynamics']
            print(f"MD RMSD: {md_data['rmsd']:.3f} Å")
            print(f"MD Energy: {md_data['average_energy']:.2f} kJ/mol")
            print("✓ Recommended for synthesis")
    else:
        print("⚠ Score too low, consider alternatives")
```

## Known Limitations

1. **Small Molecules Only**: Optimized for drug-like molecules (~10-100 atoms)
   - For proteins, use specialized tools
   - For polymers, consider alternative methods

2. **Force Field Accuracy**: GAFF/AMBER are general force fields
   - May not be highly accurate for all molecule types
   - Consider QM calculations for critical applications

3. **Timescale**: MD simulations limited to ns-μs timescales
   - Cannot observe slow processes (ms-sec)
   - Folding/binding events may not be captured

4. **Solvent Effects**: Current implementation uses vacuum or implicit solvent
   - Explicit solvent increases computational cost significantly
   - May affect accuracy for charged molecules

5. **Computational Cost**: MD simulations are expensive
   - Use GPU acceleration when possible
   - Cache results to avoid re-computation
   - Consider pre-filtering with faster methods

## References

1. **OpenMM**: Eastman et al., "OpenMM 7: Rapid development of high performance algorithms for molecular dynamics", *PLoS Comput. Biol.* 2017, 13(7), e1005659
   - Website: http://openmm.org/
   - Documentation: http://docs.openmm.org/

2. **Force Fields**:
   - GAFF: Wang et al., *J. Comput. Chem.* 2004, 25, 1157-1174
   - AMBER: https://ambermd.org/

3. **RDKit** (3D structure generation): https://www.rdkit.org/

## Contributing

To improve this adapter:

1. Add support for more force fields
2. Implement explicit solvent simulations
3. Add enhanced sampling methods (replica exchange, etc.)
4. Optimize performance for batch processing
5. Add more property calculations

## License

Part of PharmForge - see main project license.

## Support

For issues or questions:
- Check troubleshooting section above
- Review examples in `examples/molecular_dynamics_workflow.py`
- Test with `test_openmm_adapter.py`
- Consult OpenMM documentation: http://docs.openmm.org/

---

**Version**: 1.0.0
**Last Updated**: 2024-10-25
**Maintainer**: PharmForge Team
