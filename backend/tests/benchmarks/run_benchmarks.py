"""
Benchmark Runner for PharmForge Adapters

Runs all validation benchmarks and generates comprehensive reports.

Usage:
    # Run all benchmarks
    python run_benchmarks.py

    # Run specific benchmark
    python run_benchmarks.py --benchmark dude

    # Run with specific adapter
    python run_benchmarks.py --adapter vina_docking

    # Quick test mode (fewer samples)
    python run_benchmarks.py --quick

    # Save results to CSV
    python run_benchmarks.py --format csv
"""

import asyncio
import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import csv

# Add parent directories to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from backend.tests.benchmarks.dud_e_benchmark import DUDEBenchmark, EnrichmentMetrics
from backend.tests.benchmarks.tdc_admet_benchmark import TDCADMETBenchmark, BenchmarkSummary

logger = logging.getLogger(__name__)


class BenchmarkRunner:
    """
    Main benchmark runner for PharmForge adapters.

    Orchestrates running multiple benchmarks and generating summary reports.
    """

    def __init__(
        self,
        output_dir: str = "benchmark_results",
        quick_mode: bool = False
    ):
        """
        Initialize benchmark runner.

        Args:
            output_dir: Directory to save results
            quick_mode: If True, run with fewer samples for quick testing
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.quick_mode = quick_mode

        self.results: Dict[str, Any] = {
            'timestamp': datetime.now().isoformat(),
            'quick_mode': quick_mode,
            'benchmarks': {}
        }

        logger.info(f"Benchmark Runner initialized (output: {output_dir})")

    async def run_dude_benchmark(
        self,
        adapter,
        target_name: str = "TEST",
        actives: Optional[List[str]] = None,
        decoys: Optional[List[str]] = None
    ) -> EnrichmentMetrics:
        """
        Run DUD-E benchmark for a docking adapter.

        Args:
            adapter: Docking adapter instance
            target_name: Protein target name
            actives: List of active compound SMILES
            decoys: List of decoy compound SMILES

        Returns:
            EnrichmentMetrics
        """
        logger.info(f"\n{'='*70}")
        logger.info(f"Running DUD-E Benchmark: {target_name}")
        logger.info(f"{'='*70}\n")

        # Use default test data if not provided
        if actives is None or decoys is None:
            actives, decoys = self._get_default_dude_data()

        # Limit samples in quick mode
        if self.quick_mode:
            actives = actives[:5]
            decoys = decoys[:10]

        benchmark = DUDEBenchmark(
            adapter=adapter,
            target_name=target_name,
            actives=actives,
            decoys=decoys,
            output_dir=str(self.output_dir)
        )

        metrics = await benchmark.run()

        # Store results
        self.results['benchmarks'][f'dude_{target_name}_{adapter.name}'] = {
            'benchmark': 'DUD-E',
            'target': target_name,
            'adapter': adapter.name,
            'metrics': {
                'auc_roc': metrics.auc_roc,
                'auc_pr': metrics.auc_pr,
                'ef_1': metrics.enrichment_factor_1,
                'ef_5': metrics.enrichment_factor_5,
                'ef_10': metrics.enrichment_factor_10,
                'num_actives': metrics.num_actives,
                'num_decoys': metrics.num_decoys,
                'num_successful': metrics.num_successful
            }
        }

        return metrics

    async def run_tdc_admet_benchmark(
        self,
        adapter,
        properties: Optional[List[str]] = None
    ) -> BenchmarkSummary:
        """
        Run TDC ADMET benchmark for an ADMET adapter.

        Args:
            adapter: ADMET adapter instance
            properties: List of properties to test (None = all)

        Returns:
            BenchmarkSummary
        """
        logger.info(f"\n{'='*70}")
        logger.info(f"Running TDC ADMET Benchmark")
        logger.info(f"{'='*70}\n")

        # Use default properties if not provided
        if properties is None:
            properties = ['AMES', 'hERG', 'Caco2_Wang', 'BBB_Martins']

        # Limit properties in quick mode
        if self.quick_mode:
            properties = properties[:2]

        max_samples = 10 if self.quick_mode else None

        benchmark = TDCADMETBenchmark(
            adapter=adapter,
            properties=properties,
            max_samples_per_property=max_samples,
            output_dir=str(self.output_dir)
        )

        summary = await benchmark.run()

        # Store results
        self.results['benchmarks'][f'tdc_admet_{adapter.name}'] = {
            'benchmark': 'TDC ADMET',
            'adapter': adapter.name,
            'metrics': {
                'num_properties': summary.num_properties,
                'total_samples': summary.total_samples,
                'total_successful': summary.total_successful,
                'mean_regression_r2': summary.mean_regression_r2,
                'mean_classification_roc_auc': summary.mean_classification_roc_auc
            }
        }

        return summary

    async def run_all_benchmarks(
        self,
        docking_adapters: Optional[List[Any]] = None,
        admet_adapters: Optional[List[Any]] = None
    ):
        """
        Run all available benchmarks.

        Args:
            docking_adapters: List of docking adapters to test
            admet_adapters: List of ADMET adapters to test
        """
        logger.info(f"\n{'#'*70}")
        logger.info("RUNNING ALL BENCHMARKS")
        logger.info(f"{'#'*70}\n")

        start_time = time.time()

        # Run DUD-E benchmarks for docking adapters
        if docking_adapters:
            for adapter in docking_adapters:
                try:
                    await self.run_dude_benchmark(adapter)
                except Exception as e:
                    logger.error(f"DUD-E benchmark failed for {adapter.name}: {e}")

        # Run TDC ADMET benchmarks for ADMET adapters
        if admet_adapters:
            for adapter in admet_adapters:
                try:
                    await self.run_tdc_admet_benchmark(adapter)
                except Exception as e:
                    logger.error(f"TDC ADMET benchmark failed for {adapter.name}: {e}")

        elapsed_time = time.time() - start_time

        # Generate summary report
        self._generate_summary_report(elapsed_time)

    def _generate_summary_report(self, elapsed_time: float):
        """Generate summary report of all benchmarks."""
        logger.info(f"\n{'='*70}")
        logger.info("BENCHMARK SUMMARY")
        logger.info(f"{'='*70}\n")

        logger.info(f"Total benchmarks run: {len(self.results['benchmarks'])}")
        logger.info(f"Total time: {elapsed_time:.1f}s")
        logger.info("")

        # Print summary for each benchmark
        for name, result in self.results['benchmarks'].items():
            logger.info(f"{result['benchmark']} - {result['adapter']}:")

            metrics = result['metrics']
            if result['benchmark'] == 'DUD-E':
                logger.info(f"  AUC-ROC: {metrics['auc_roc']:.3f}")
                logger.info(f"  EF 1%:   {metrics['ef_1']:.2f}x")
            elif result['benchmark'] == 'TDC ADMET':
                logger.info(f"  Properties: {metrics['num_properties']}")
                if metrics['mean_classification_roc_auc']:
                    logger.info(f"  Mean ROC-AUC: {metrics['mean_classification_roc_auc']:.3f}")
                if metrics['mean_regression_r2']:
                    logger.info(f"  Mean R²:      {metrics['mean_regression_r2']:.3f}")

            logger.info("")

        # Save summary to JSON
        self._save_summary_json()

    def _save_summary_json(self):
        """Save summary results to JSON file."""
        output_file = self.output_dir / "benchmark_summary.json"

        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        logger.info(f"Summary saved to: {output_file}")

    def save_results_csv(self):
        """Save results to CSV format."""
        output_file = self.output_dir / "benchmark_summary.csv"

        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)

            # Write header
            writer.writerow([
                'Benchmark',
                'Adapter',
                'Metric',
                'Value'
            ])

            # Write data
            for name, result in self.results['benchmarks'].items():
                benchmark = result['benchmark']
                adapter = result['adapter']

                for metric_name, value in result['metrics'].items():
                    if value is not None:
                        writer.writerow([benchmark, adapter, metric_name, value])

        logger.info(f"CSV saved to: {output_file}")

    @staticmethod
    def _get_default_dude_data() -> tuple:
        """Get default DUD-E test data."""
        # Sample actives (known drugs)
        actives = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "CN1CCC[C@H]1c2cccnc2",  # Nicotine
            "CC(C)NCC(O)c1ccc(O)c(CO)c1",  # Salbutamol
        ]

        # Sample decoys
        decoys = [
            "CCO",  # Ethanol
            "CCCCCC",  # Hexane
            "c1ccccc1",  # Benzene
            "CC(C)O",  # Isopropanol
            "CCC(C)O",  # 2-Butanol
            "CCCCCCCC",  # Octane
            "CC(C)C",  # Isobutane
            "CCCCC",  # Pentane
            "CCOC",  # Dimethyl ether
            "CC(C)CC(C)C",  # 2,4-Dimethylpentane
        ]

        return actives, decoys


async def main():
    """Main entry point for benchmark runner."""
    parser = argparse.ArgumentParser(
        description="Run PharmForge adapter benchmarks"
    )
    parser.add_argument(
        '--benchmark',
        choices=['dude', 'tdc', 'all'],
        default='all',
        help="Which benchmark to run (default: all)"
    )
    parser.add_argument(
        '--adapter',
        help="Specific adapter name to test"
    )
    parser.add_argument(
        '--quick',
        action='store_true',
        help="Run in quick mode (fewer samples)"
    )
    parser.add_argument(
        '--output',
        default='benchmark_results',
        help="Output directory for results"
    )
    parser.add_argument(
        '--format',
        choices=['json', 'csv', 'both'],
        default='json',
        help="Output format (default: json)"
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Initialize runner
    runner = BenchmarkRunner(
        output_dir=args.output,
        quick_mode=args.quick
    )

    # Create mock adapters for demonstration
    # In production, you would import and configure real adapters
    logger.info("NOTE: Using mock adapters for demonstration")
    logger.info("In production, configure real adapters from PharmForge\n")

    # Mock docking adapter
    class MockDockingAdapter:
        def __init__(self):
            self.name = "mock_docking"
            self.version = "1.0.0"

        async def execute(self, smiles):
            import random
            from backend.core.adapters.protocol import AdapterResult

            affinity = random.uniform(-10, -5)
            return AdapterResult(
                success=True,
                data={
                    'binding_affinity': affinity,
                    'binding_score': (affinity + 10) / 5
                },
                error=None,
                metadata={'adapter_name': self.name}
            )

    # Mock ADMET adapter
    class MockADMETAdapter:
        def __init__(self):
            self.name = "mock_admet"
            self.version = "1.0.0"

        async def execute(self, smiles):
            import random
            from backend.core.adapters.protocol import AdapterResult

            properties = {
                'AMES': random.choice([0, 1]),
                'hERG': random.choice([0, 1]),
                'Caco2_Wang': random.uniform(-7, -4),
                'BBB_Martins': random.choice([0, 1]),
            }

            return AdapterResult(
                success=True,
                data={'properties': properties},
                error=None,
                metadata={'adapter_name': self.name}
            )

    # Run benchmarks
    try:
        if args.benchmark == 'dude' or args.benchmark == 'all':
            docking_adapter = MockDockingAdapter()
            await runner.run_dude_benchmark(docking_adapter)

        if args.benchmark == 'tdc' or args.benchmark == 'all':
            admet_adapter = MockADMETAdapter()
            await runner.run_tdc_admet_benchmark(admet_adapter)

        # Save results
        if args.format in ['csv', 'both']:
            runner.save_results_csv()

        logger.info("\n" + "="*70)
        logger.info("BENCHMARKS COMPLETE")
        logger.info("="*70)
        logger.info(f"\nResults saved to: {args.output}/")

    except KeyboardInterrupt:
        logger.warning("\nBenchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\nBenchmark failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())
