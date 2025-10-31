"""
Example Usage of PharmForge Benchmark Suite

This script demonstrates how to use the benchmark suite with real PharmForge adapters.

Examples:
1. DUD-E benchmark with Vina adapter
2. TDC ADMET benchmark with ADMET-AI adapter
3. Comparing multiple adapters
4. Custom benchmark configurations
"""

import asyncio
import logging
import sys
from pathlib import Path

# Add parent directories to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from backend.tests.benchmarks.dud_e_benchmark import DUDEBenchmark
from backend.tests.benchmarks.tdc_admet_benchmark import TDCADMETBenchmark
from backend.tests.benchmarks.run_benchmarks import BenchmarkRunner

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Example 1: DUD-E Benchmark with Vina Adapter
# ============================================================================

async def example_1_dude_vina():
    """
    Example: Run DUD-E benchmark with AutoDock Vina adapter.

    This tests virtual screening enrichment for a specific protein target.
    """
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE 1: DUD-E Benchmark with Vina")
    logger.info("="*70 + "\n")

    # Import Vina adapter
    try:
        from adapters.vina.adapter import VinaAdapter
    except ImportError:
        logger.warning("Vina adapter not available, using mock adapter")
        from backend.tests.benchmarks.test_benchmarks import MockDockingAdapter
        adapter = MockDockingAdapter(name="vina_docking")
    else:
        # Configure Vina adapter
        adapter = VinaAdapter(config={
            'receptor_path': '/path/to/receptor.pdbqt',  # UPDATE THIS
            'center_x': 10.0,  # UPDATE THESE
            'center_y': 20.0,
            'center_z': 15.0,
            'size_x': 25,
            'size_y': 25,
            'size_z': 25,
        })

    # Sample actives (known binders) - replace with real DUD-E data
    actives = [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
    ]

    # Sample decoys (property-matched non-binders) - replace with real DUD-E data
    decoys = [
        "CCO",  # Ethanol
        "c1ccccc1",  # Benzene
        "CC(C)O",  # Isopropanol
        "CCCCCC",  # Hexane
        "CCC(C)O",  # 2-Butanol
    ]

    # Create and run benchmark
    benchmark = DUDEBenchmark(
        adapter=adapter,
        target_name="EGFR",  # Example target
        actives=actives,
        decoys=decoys,
        output_dir="benchmark_results"
    )

    metrics = await benchmark.run()

    # Print results
    logger.info("\nResults:")
    logger.info(f"  AUC-ROC: {metrics.auc_roc:.3f}")
    logger.info(f"  AUC-PR:  {metrics.auc_pr:.3f}")
    logger.info(f"  EF 1%:   {metrics.enrichment_factor_1:.2f}x")
    logger.info(f"  EF 5%:   {metrics.enrichment_factor_5:.2f}x")
    logger.info(f"  EF 10%:  {metrics.enrichment_factor_10:.2f}x")

    if metrics.mean_active_affinity:
        logger.info(f"\n  Mean Active Affinity: {metrics.mean_active_affinity:.2f} kcal/mol")
        logger.info(f"  Mean Decoy Affinity:  {metrics.mean_decoy_affinity:.2f} kcal/mol")
        logger.info(f"  Separation:           {metrics.separation_index:.2f} σ")


# ============================================================================
# Example 2: TDC ADMET Benchmark with ADMET-AI Adapter
# ============================================================================

async def example_2_tdc_admet():
    """
    Example: Run TDC ADMET benchmark with ADMET-AI adapter.

    This tests ADMET property prediction accuracy against published benchmarks.
    """
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE 2: TDC ADMET Benchmark with ADMET-AI")
    logger.info("="*70 + "\n")

    # Import ADMET-AI adapter
    try:
        from adapters.admet_ai.adapter import ADMETaiAdapter
    except ImportError:
        logger.warning("ADMET-AI adapter not available, using mock adapter")
        from backend.tests.benchmarks.test_benchmarks import MockADMETAdapter
        adapter = MockADMETAdapter()
    else:
        # Configure ADMET-AI adapter
        adapter = ADMETaiAdapter()

    # Select properties to test
    properties = [
        # Toxicity
        'AMES',              # Mutagenicity
        'hERG',              # Cardiotoxicity

        # Absorption
        'Caco2_Wang',        # Intestinal permeability

        # Distribution
        'BBB_Martins',       # Blood-brain barrier

        # Metabolism
        'CYP3A4_Veith',      # CYP3A4 inhibition
    ]

    # Create and run benchmark
    benchmark = TDCADMETBenchmark(
        adapter=adapter,
        properties=properties,
        max_samples_per_property=50,  # Limit for quick testing
        output_dir="benchmark_results"
    )

    summary = await benchmark.run()

    # Print results
    logger.info("\nResults:")
    logger.info(f"  Properties tested: {summary.num_properties}")
    logger.info(f"  Total samples:     {summary.total_samples}")
    logger.info(f"  Success rate:      {100*summary.total_successful/summary.total_samples:.1f}%")

    if summary.mean_classification_roc_auc:
        logger.info(f"\n  Mean Classification ROC-AUC: {summary.mean_classification_roc_auc:.3f}")

    if summary.mean_regression_r2:
        logger.info(f"  Mean Regression R²:          {summary.mean_regression_r2:.3f}")

    # Per-property results
    logger.info("\n  Per-Property Results:")
    for metrics in summary.property_metrics:
        logger.info(f"\n    {metrics.property_name}:")
        logger.info(f"      Type:    {metrics.task_type}")

        if metrics.task_type == 'classification' and metrics.roc_auc:
            logger.info(f"      ROC-AUC: {metrics.roc_auc:.3f}")
            if metrics.baseline_metric:
                logger.info(f"      Baseline: {metrics.baseline_metric:.3f}")
                logger.info(f"      Improvement: {metrics.improvement_over_baseline:+.1f}%")

        elif metrics.task_type == 'regression' and metrics.r2:
            logger.info(f"      R²:      {metrics.r2:.3f}")
            logger.info(f"      MAE:     {metrics.mae:.3f}")
            if metrics.baseline_metric:
                logger.info(f"      Baseline: {metrics.baseline_metric:.3f}")
                logger.info(f"      Improvement: {metrics.improvement_over_baseline:+.1f}%")


# ============================================================================
# Example 3: Comparing Multiple Adapters
# ============================================================================

async def example_3_compare_adapters():
    """
    Example: Compare multiple docking adapters.

    This demonstrates running the same benchmark on different adapters
    to compare their performance.
    """
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE 3: Comparing Multiple Adapters")
    logger.info("="*70 + "\n")

    # Import adapters (using mocks for demonstration)
    from backend.tests.benchmarks.test_benchmarks import MockDockingAdapter

    adapters = [
        MockDockingAdapter(name="vina", bias_actives=True),
        MockDockingAdapter(name="gnina", bias_actives=True),
        MockDockingAdapter(name="random", bias_actives=False),
    ]

    # Sample data
    actives = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    ]

    decoys = [
        "CCO", "c1ccccc1", "CC(C)O", "CCCCCC", "CCC(C)O"
    ]

    # Run benchmark for each adapter
    results = {}

    for adapter in adapters:
        logger.info(f"\nTesting {adapter.name}...")

        benchmark = DUDEBenchmark(
            adapter=adapter,
            target_name="TEST",
            actives=actives,
            decoys=decoys,
            output_dir="benchmark_results"
        )

        metrics = await benchmark.run()
        results[adapter.name] = metrics

    # Compare results
    logger.info("\n" + "="*70)
    logger.info("COMPARISON RESULTS")
    logger.info("="*70 + "\n")

    logger.info(f"{'Adapter':<15} {'AUC-ROC':<10} {'EF 1%':<10} {'EF 5%':<10}")
    logger.info("-" * 50)

    for name, metrics in results.items():
        logger.info(
            f"{name:<15} "
            f"{metrics.auc_roc:<10.3f} "
            f"{metrics.enrichment_factor_1:<10.2f} "
            f"{metrics.enrichment_factor_5:<10.2f}"
        )

    # Find best adapter
    best_adapter = max(results.items(), key=lambda x: x[1].auc_roc)
    logger.info(f"\nBest Adapter: {best_adapter[0]} (AUC-ROC: {best_adapter[1].auc_roc:.3f})")


# ============================================================================
# Example 4: Using the Benchmark Runner
# ============================================================================

async def example_4_benchmark_runner():
    """
    Example: Use the benchmark runner to run multiple benchmarks.

    This demonstrates the high-level API for running comprehensive benchmarks.
    """
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE 4: Using Benchmark Runner")
    logger.info("="*70 + "\n")

    # Create mock adapters
    from backend.tests.benchmarks.test_benchmarks import (
        MockDockingAdapter,
        MockADMETAdapter
    )

    docking_adapter = MockDockingAdapter()
    admet_adapter = MockADMETAdapter()

    # Create runner
    runner = BenchmarkRunner(
        output_dir="benchmark_results",
        quick_mode=True  # Use quick mode for demonstration
    )

    # Run DUD-E benchmark
    await runner.run_dude_benchmark(
        adapter=docking_adapter,
        target_name="EGFR"
    )

    # Run TDC ADMET benchmark
    await runner.run_tdc_admet_benchmark(
        adapter=admet_adapter,
        properties=['AMES', 'hERG']
    )

    # Save results in CSV format
    runner.save_results_csv()

    logger.info("\nBenchmarks complete! Check benchmark_results/ for outputs.")


# ============================================================================
# Example 5: Loading Real DUD-E Data
# ============================================================================

def example_5_load_dude_data():
    """
    Example: How to load real DUD-E data.

    This shows how to download and parse actual DUD-E datasets.
    """
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE 5: Loading Real DUD-E Data")
    logger.info("="*70 + "\n")

    logger.info("""
To use real DUD-E data:

1. Download DUD-E dataset for your target:
   Visit: http://dude.docking.org/targets

2. Extract actives and decoys:
   - actives_final.ism (active compounds)
   - decoys_final.ism (decoy compounds)

3. Load in Python:

   import pandas as pd

   # Load actives
   actives_df = pd.read_csv('actives_final.ism', sep=' ', header=None, names=['SMILES', 'ID'])
   actives = actives_df['SMILES'].tolist()

   # Load decoys
   decoys_df = pd.read_csv('decoys_final.ism', sep=' ', header=None, names=['SMILES', 'ID'])
   decoys = decoys_df['SMILES'].tolist()

   # Run benchmark
   benchmark = DUDEBenchmark(
       adapter=adapter,
       target_name="EGFR",
       actives=actives,
       decoys=decoys
   )

   metrics = await benchmark.run()

Available targets:
- EGFR (Epidermal Growth Factor Receptor)
- ESR1 (Estrogen Receptor)
- PPARG (Peroxisome Proliferator-Activated Receptor Gamma)
- ... and 99 more targets

For a complete list, visit: http://dude.docking.org/targets
    """)


# ============================================================================
# Example 6: Loading Real TDC Data
# ============================================================================

def example_6_load_tdc_data():
    """
    Example: How to load real TDC data.

    This shows how to use the TDC library to get benchmark datasets.
    """
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE 6: Loading Real TDC Data")
    logger.info("="*70 + "\n")

    logger.info("""
To use real TDC data:

1. Install TDC library:
   pip install PyTDC

2. Load benchmark data:

   from tdc.single_pred import ADME, Tox

   # Example: Load AMES toxicity data
   data = Tox(name='AMES')
   df = data.get_data()

   # Get train/test split
   split = data.get_split()
   train_df = split['train']
   test_df = split['test']

   # Use test set for benchmarking
   samples = [
       {'smiles': row['Drug'], 'y': row['Y']}
       for _, row in test_df.iterrows()
   ]

3. Available datasets:

   Absorption:
   - Caco2_Wang (permeability)
   - HIA_Hou (human intestinal absorption)
   - Bioavailability_Ma

   Distribution:
   - BBB_Martins (blood-brain barrier)
   - PPBR_AZ (plasma protein binding)
   - VDss_Lombardo (volume of distribution)

   Metabolism:
   - CYP1A2_Veith, CYP2C9_Veith, CYP2C19_Veith
   - CYP2D6_Veith, CYP3A4_Veith
   - Clearance_Hepatocyte_AZ

   Excretion:
   - Half_Life_Obach
   - Clearance_Microsome_AZ

   Toxicity:
   - AMES (mutagenicity)
   - hERG (cardiotoxicity)
   - DILI (drug-induced liver injury)

For full documentation: https://tdcommons.ai/
    """)


# ============================================================================
# Main
# ============================================================================

async def main():
    """Run all examples."""
    import argparse

    parser = argparse.ArgumentParser(description="PharmForge Benchmark Examples")
    parser.add_argument(
        '--example',
        type=int,
        choices=[1, 2, 3, 4, 5, 6],
        help="Run specific example (1-6)"
    )

    args = parser.parse_args()

    if args.example == 1 or args.example is None:
        await example_1_dude_vina()

    if args.example == 2 or args.example is None:
        await example_2_tdc_admet()

    if args.example == 3 or args.example is None:
        await example_3_compare_adapters()

    if args.example == 4 or args.example is None:
        await example_4_benchmark_runner()

    if args.example == 5:
        example_5_load_dude_data()

    if args.example == 6:
        example_6_load_tdc_data()


if __name__ == "__main__":
    asyncio.run(main())
