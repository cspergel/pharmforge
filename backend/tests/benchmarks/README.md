# PharmForge Validation Benchmark Suite

Comprehensive benchmarking system for validating PharmForge adapter accuracy and performance against established datasets.

## Overview

This benchmark suite provides standardized tests for:

1. **DUD-E (Database of Useful Decoys: Enhanced)** - Tests molecular docking enrichment
2. **TDC ADMET** - Tests ADMET prediction accuracy against Therapeutics Data Commons

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Benchmarks](#benchmarks)
  - [DUD-E Benchmark](#dud-e-benchmark)
  - [TDC ADMET Benchmark](#tdc-admet-benchmark)
- [Running Benchmarks](#running-benchmarks)
- [Interpreting Results](#interpreting-results)
- [Advanced Usage](#advanced-usage)
- [Adding Custom Benchmarks](#adding-custom-benchmarks)

## Installation

### Requirements

```bash
# Core dependencies
pip install numpy scikit-learn

# For TDC benchmarks
pip install PyTDC

# For real adapter testing
pip install rdkit admet-ai

# Optional: For visualization
pip install matplotlib seaborn
```

### Setup

```bash
cd backend/tests/benchmarks
```

## Quick Start

### Run All Benchmarks (Quick Mode)

```bash
python run_benchmarks.py --quick
```

This runs a quick validation with reduced sample sizes (~1 minute).

### Run Specific Benchmark

```bash
# DUD-E benchmark only
python run_benchmarks.py --benchmark dude

# TDC ADMET benchmark only
python run_benchmarks.py --benchmark tdc
```

### Full Benchmark Suite

```bash
# Run with full datasets (may take hours)
python run_benchmarks.py --benchmark all
```

## Benchmarks

### DUD-E Benchmark

Tests the ability of docking adapters to distinguish known active compounds from property-matched decoys.

**Purpose**: Validate virtual screening enrichment

**Dataset**:
- 102 protein targets
- ~22,000 actives
- ~50 decoys per active (matched by physicochemical properties)

**Metrics**:
- **AUC-ROC**: Area under ROC curve (0.5 = random, 1.0 = perfect)
- **AUC-PR**: Area under precision-recall curve
- **Enrichment Factors (EF)**: Fold-improvement over random at 1%, 5%, 10% of dataset
- **Separation Index**: Mean difference in binding affinity (in standard deviations)

**Usage**:

```python
from backend.tests.benchmarks.dud_e_benchmark import DUDEBenchmark
from adapters.vina.adapter import VinaAdapter

# Configure adapter
adapter = VinaAdapter(config={
    'receptor_path': '/path/to/receptor.pdbqt',
    'center_x': 10.0,
    'center_y': 20.0,
    'center_z': 15.0,
})

# Load DUD-E data
actives = [...]  # Load from DUD-E dataset
decoys = [...]

# Run benchmark
benchmark = DUDEBenchmark(
    adapter=adapter,
    target_name="EGFR",
    actives=actives,
    decoys=decoys
)

metrics = await benchmark.run()
print(f"AUC-ROC: {metrics.auc_roc:.3f}")
print(f"EF 1%: {metrics.enrichment_factor_1:.2f}x")
```

**Expected Results** (Literature Baselines):

| Method | AUC-ROC | EF 1% | EF 5% |
|--------|---------|-------|-------|
| Random | 0.50 | 1.0x | 1.0x |
| AutoDock Vina | 0.70-0.80 | 5-15x | 3-8x |
| GNINA (CNN) | 0.75-0.85 | 10-25x | 5-12x |
| DiffDock | 0.80-0.90 | 15-35x | 8-18x |

### TDC ADMET Benchmark

Tests ADMET property prediction accuracy against curated TDC datasets with published baselines.

**Purpose**: Validate ADMET prediction models

**Dataset**:
- 49 ADMET properties
- Absorption, Distribution, Metabolism, Excretion, Toxicity endpoints
- Regression and classification tasks

**Properties Tested**:

| Category | Properties |
|----------|-----------|
| **Absorption** | Caco2, HIA, Bioavailability, Lipophilicity |
| **Distribution** | BBB, PPBR, VDss |
| **Metabolism** | CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4, Clearance |
| **Excretion** | Half-life, Clearance |
| **Toxicity** | AMES, hERG, DILI, Carcinogenicity, Skin Sensitization |

**Metrics**:

*Regression Tasks*:
- **R²** (Coefficient of Determination): Variance explained (1.0 = perfect)
- **MAE** (Mean Absolute Error): Average prediction error
- **RMSE** (Root Mean Squared Error): RMS prediction error

*Classification Tasks*:
- **ROC-AUC**: Area under ROC curve
- **Accuracy**: Overall classification accuracy
- **F1 Score**: Harmonic mean of precision and recall

**Usage**:

```python
from backend.tests.benchmarks.tdc_admet_benchmark import TDCADMETBenchmark
from adapters.admet_ai.adapter import ADMETaiAdapter

# Configure adapter
adapter = ADMETaiAdapter()

# Run benchmark
benchmark = TDCADMETBenchmark(
    adapter=adapter,
    properties=['AMES', 'hERG', 'Caco2_Wang', 'BBB_Martins']
)

summary = await benchmark.run()
print(f"Mean Classification ROC-AUC: {summary.mean_classification_roc_auc:.3f}")
print(f"Mean Regression R²: {summary.mean_regression_r2:.3f}")
```

**Expected Results** (TDC Baselines):

| Property | Type | Baseline | Good | Excellent |
|----------|------|----------|------|-----------|
| AMES | Classification | 0.85 ROC-AUC | >0.88 | >0.92 |
| hERG | Classification | 0.82 ROC-AUC | >0.85 | >0.90 |
| Caco2 | Regression | 0.45 R² | >0.55 | >0.65 |
| BBB | Classification | 0.91 ROC-AUC | >0.93 | >0.95 |
| CYP3A4 | Classification | 0.85 ROC-AUC | >0.87 | >0.92 |

## Running Benchmarks

### Command-Line Interface

```bash
# Basic usage
python run_benchmarks.py

# Options
python run_benchmarks.py \
    --benchmark {dude|tdc|all} \
    --quick \
    --output benchmark_results \
    --format {json|csv|both} \
    --verbose
```

**Arguments**:
- `--benchmark`: Which benchmark to run (default: all)
- `--quick`: Quick mode with reduced samples (default: False)
- `--output`: Output directory (default: benchmark_results)
- `--format`: Output format (default: json)
- `--verbose`: Enable debug logging (default: False)

### Programmatic Usage

```python
from backend.tests.benchmarks.run_benchmarks import BenchmarkRunner

runner = BenchmarkRunner(
    output_dir="benchmark_results",
    quick_mode=False
)

# Run specific benchmarks
await runner.run_dude_benchmark(docking_adapter)
await runner.run_tdc_admet_benchmark(admet_adapter)

# Or run all
await runner.run_all_benchmarks(
    docking_adapters=[vina_adapter, gnina_adapter],
    admet_adapters=[admet_ai_adapter]
)
```

## Interpreting Results

### Output Files

After running benchmarks, you'll find:

```
benchmark_results/
├── benchmark_summary.json           # Overall summary
├── benchmark_summary.csv            # CSV format
├── dud_e_EGFR_vina_docking.json    # Detailed DUD-E results
├── dud_e_EGFR_vina_docking_summary.txt
├── tdc_admet_admet_ai.json          # Detailed TDC results
└── tdc_admet_admet_ai_summary.txt
```

### Understanding Metrics

#### DUD-E Metrics

**AUC-ROC (Area Under ROC Curve)**:
- Random: 0.5
- Acceptable: >0.7
- Good: >0.8
- Excellent: >0.9

**Enrichment Factor (EF)**:
- Measures fold-improvement over random selection
- EF 1% = How many more actives found in top 1% vs random
- Good: >5x
- Excellent: >10x

**Interpretation**:
```
AUC-ROC = 0.82, EF1% = 12.5x
→ "Good discrimination; top 1% contains 12.5x more actives than random"
```

#### TDC ADMET Metrics

**R² (Regression)**:
- 0.0 = No predictive power
- 0.5 = Moderate
- 0.8 = Good
- >0.9 = Excellent

**ROC-AUC (Classification)**:
- 0.5 = Random
- 0.7 = Acceptable
- 0.8 = Good
- >0.9 = Excellent

**Baseline Comparison**:
```
R² = 0.65, Baseline = 0.45
Improvement = +44% over baseline
→ "Significantly better than published methods"
```

### Performance Tiers

| Tier | DUD-E AUC-ROC | TDC Mean ROC-AUC | TDC Mean R² |
|------|---------------|------------------|-------------|
| Poor | <0.6 | <0.7 | <0.3 |
| Acceptable | 0.6-0.7 | 0.7-0.8 | 0.3-0.5 |
| Good | 0.7-0.8 | 0.8-0.85 | 0.5-0.7 |
| Excellent | 0.8-0.9 | 0.85-0.9 | 0.7-0.8 |
| Outstanding | >0.9 | >0.9 | >0.8 |

## Advanced Usage

### Custom DUD-E Data

```python
# Load your own actives/decoys
import pandas as pd

actives_df = pd.read_csv('my_actives.csv')
decoys_df = pd.read_csv('my_decoys.csv')

actives = actives_df['SMILES'].tolist()
decoys = decoys_df['SMILES'].tolist()

benchmark = DUDEBenchmark(
    adapter=adapter,
    target_name="MyTarget",
    actives=actives,
    decoys=decoys,
    max_actives=100,  # Limit for testing
    max_decoys=500
)
```

### Custom TDC Properties

```python
# Test specific ADMET properties
properties = [
    'AMES',           # Mutagenicity
    'hERG',           # Cardiotoxicity
    'CYP3A4_Veith',   # Metabolism
    'BBB_Martins',    # Blood-brain barrier
    'Caco2_Wang',     # Intestinal absorption
]

benchmark = TDCADMETBenchmark(
    adapter=adapter,
    properties=properties,
    max_samples_per_property=200
)
```

### Loading Real TDC Data

To use real TDC data instead of mock data:

```python
# Install TDC
# pip install PyTDC

# Modify tdc_admet_benchmark.py to use real data loader
# Uncomment the _load_tdc_property() method implementation

from tdc.single_pred import ADME, Tox

# Example: Load Caco2 data
data = ADME(name='Caco2_Wang')
df = data.get_data()

# Use in benchmark
samples = [
    {'smiles': row['Drug'], 'y': row['Y']}
    for _, row in df.iterrows()
]
```

### Batch Testing Multiple Adapters

```python
# Compare multiple adapters
adapters = [
    VinaAdapter(config={...}),
    GNINAAdapter(config={...}),
    DiffDockAdapter(config={...})
]

results = {}
for adapter in adapters:
    benchmark = DUDEBenchmark(adapter=adapter, ...)
    metrics = await benchmark.run()
    results[adapter.name] = metrics

# Compare
for name, metrics in results.items():
    print(f"{name}: AUC-ROC={metrics.auc_roc:.3f}")
```

## Adding Custom Benchmarks

To add a new benchmark:

1. Create new file: `backend/tests/benchmarks/my_benchmark.py`

2. Implement benchmark class:

```python
from typing import List
from dataclasses import dataclass

@dataclass
class MyMetrics:
    metric1: float
    metric2: float

class MyBenchmark:
    def __init__(self, adapter, **kwargs):
        self.adapter = adapter

    async def run(self) -> MyMetrics:
        # Implement benchmark logic
        results = []

        for sample in self.samples:
            result = await self.adapter.execute(sample)
            results.append(result)

        metrics = self._calculate_metrics(results)
        self._save_results(metrics)

        return metrics
```

3. Add to runner:

```python
# In run_benchmarks.py
async def run_my_benchmark(self, adapter):
    benchmark = MyBenchmark(adapter=adapter)
    return await benchmark.run()
```

## Troubleshooting

### Common Issues

**Issue**: `ImportError: No module named 'tdc'`
```bash
Solution: pip install PyTDC
```

**Issue**: `RDKit not available`
```bash
Solution: pip install rdkit
```

**Issue**: Benchmark runs but all predictions fail
```bash
Check:
1. Adapter configuration (receptor file, model path, etc.)
2. Input SMILES validity
3. Adapter logs for specific errors
```

**Issue**: Low performance metrics
```bash
Possible causes:
1. Adapter not properly configured
2. Training data mismatch
3. Model not loaded correctly
4. Using default/random predictions

Check adapter.execute() returns valid predictions
```

### Performance Optimization

For large benchmarks:

```python
# Use quick mode for testing
runner = BenchmarkRunner(quick_mode=True)

# Limit samples
benchmark = TDCADMETBenchmark(
    adapter=adapter,
    max_samples_per_property=100  # Instead of full dataset
)

# Test subset of properties
properties = ['AMES', 'hERG']  # Instead of all 49
```

## References

### DUD-E
- Mysinger et al. (2012) "Directory of Useful Decoys, Enhanced (DUD-E): Better Ligands and Decoys for Better Benchmarking" *J. Med. Chem.* 55(14), 6582-6594
- Website: http://dude.docking.org/

### TDC
- Huang et al. (2021) "Therapeutics Data Commons: Machine Learning Datasets and Tasks for Drug Discovery and Development" *NeurIPS Datasets and Benchmarks*
- Website: https://tdcommons.ai/
- GitHub: https://github.com/mims-harvard/TDC

### ADMET-AI
- Swanson et al. (2023) "ADMET-AI: A machine learning ADMET platform for evaluation of large-scale chemical libraries"
- GitHub: https://github.com/swansonk14/admet_ai

## License

This benchmark suite is part of PharmForge and follows the same license.

## Contributing

To contribute new benchmarks:

1. Follow the existing benchmark structure
2. Include proper documentation
3. Add example usage
4. Include baseline metrics from literature
5. Add tests

## Support

For issues or questions:
- Check the troubleshooting section
- Review adapter documentation
- Open an issue on GitHub
