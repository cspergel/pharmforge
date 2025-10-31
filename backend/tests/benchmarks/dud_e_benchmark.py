"""
DUD-E (Database of Useful Decoys: Enhanced) Benchmark

Tests molecular docking enrichment by ranking known actives versus decoys.
Evaluates the ability of docking adapters to distinguish true binders from non-binders.

DUD-E is a widely-used benchmarking dataset for virtual screening, containing
102 protein targets with ~22,000 actives and ~50 decoys per active (matched by
physical properties but topologically dissimilar).

Reference:
    Mysinger et al. (2012) J. Chem. Inf. Model. 52, 3103-3125
    http://dude.docking.org/
"""

import asyncio
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve, auc

logger = logging.getLogger(__name__)


@dataclass
class DockingResult:
    """Result of a single docking run"""
    smiles: str
    is_active: bool
    binding_affinity: Optional[float] = None
    binding_score: Optional[float] = None
    success: bool = False
    error: Optional[str] = None
    runtime_ms: Optional[float] = None


@dataclass
class EnrichmentMetrics:
    """Metrics for evaluating docking enrichment"""
    auc_roc: float
    auc_pr: float  # Area under precision-recall curve
    enrichment_factor_1: float  # EF at 1% of dataset
    enrichment_factor_5: float  # EF at 5% of dataset
    enrichment_factor_10: float  # EF at 10% of dataset
    num_actives: int
    num_decoys: int
    num_successful: int
    mean_active_affinity: Optional[float] = None
    mean_decoy_affinity: Optional[float] = None
    separation_index: Optional[float] = None  # (mean_active - mean_decoy) / std


class DUDEBenchmark:
    """
    DUD-E benchmark for docking adapter validation.

    Tests the ability to enrich known actives over decoys using molecular docking.
    """

    def __init__(
        self,
        adapter,
        target_name: str,
        actives: List[str],
        decoys: List[str],
        max_actives: Optional[int] = None,
        max_decoys: Optional[int] = None,
        output_dir: str = "benchmark_results",
    ):
        """
        Initialize DUD-E benchmark.

        Args:
            adapter: Docking adapter instance (e.g., VinaAdapter, GNINAAdapter)
            target_name: Name of protein target (e.g., "EGFR", "ESR1")
            actives: List of SMILES strings for active compounds
            decoys: List of SMILES strings for decoy compounds
            max_actives: Maximum number of actives to test (None = all)
            max_decoys: Maximum number of decoys to test (None = all)
            output_dir: Directory to save benchmark results
        """
        self.adapter = adapter
        self.target_name = target_name
        self.actives = actives[:max_actives] if max_actives else actives
        self.decoys = decoys[:max_decoys] if max_decoys else decoys
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.results: List[DockingResult] = []

        logger.info(
            f"DUD-E Benchmark initialized: {target_name} "
            f"({len(self.actives)} actives, {len(self.decoys)} decoys)"
        )

    async def run(self) -> EnrichmentMetrics:
        """
        Run the benchmark: dock all actives and decoys, calculate metrics.

        Returns:
            EnrichmentMetrics with AUC-ROC, enrichment factors, etc.
        """
        logger.info(f"Starting DUD-E benchmark for {self.target_name}...")
        start_time = time.time()

        # Dock actives
        logger.info(f"Docking {len(self.actives)} active compounds...")
        active_results = await self._dock_molecules(self.actives, is_active=True)
        self.results.extend(active_results)

        # Dock decoys
        logger.info(f"Docking {len(self.decoys)} decoy compounds...")
        decoy_results = await self._dock_molecules(self.decoys, is_active=False)
        self.results.extend(decoy_results)

        # Calculate metrics
        metrics = self._calculate_metrics()

        # Save results
        self._save_results(metrics)

        elapsed_time = time.time() - start_time
        logger.info(
            f"Benchmark complete in {elapsed_time:.1f}s: "
            f"AUC-ROC={metrics.auc_roc:.3f}, EF1%={metrics.enrichment_factor_1:.2f}"
        )

        return metrics

    async def _dock_molecules(
        self,
        smiles_list: List[str],
        is_active: bool
    ) -> List[DockingResult]:
        """
        Dock a list of molecules.

        Args:
            smiles_list: List of SMILES strings
            is_active: Whether these are actives (True) or decoys (False)

        Returns:
            List of DockingResult objects
        """
        results = []

        for i, smiles in enumerate(smiles_list):
            if (i + 1) % 10 == 0:
                logger.info(
                    f"  Progress: {i + 1}/{len(smiles_list)} "
                    f"{'actives' if is_active else 'decoys'}"
                )

            start_time = time.time()

            try:
                # Execute docking
                adapter_result = await self.adapter.execute(smiles)

                runtime_ms = (time.time() - start_time) * 1000

                if adapter_result.success:
                    result = DockingResult(
                        smiles=smiles,
                        is_active=is_active,
                        binding_affinity=adapter_result.data.get('binding_affinity'),
                        binding_score=adapter_result.data.get('binding_score'),
                        success=True,
                        runtime_ms=runtime_ms
                    )
                else:
                    result = DockingResult(
                        smiles=smiles,
                        is_active=is_active,
                        success=False,
                        error=adapter_result.error,
                        runtime_ms=runtime_ms
                    )

                results.append(result)

            except Exception as e:
                logger.error(f"Error docking {smiles}: {e}")
                results.append(
                    DockingResult(
                        smiles=smiles,
                        is_active=is_active,
                        success=False,
                        error=str(e)
                    )
                )

        return results

    def _calculate_metrics(self) -> EnrichmentMetrics:
        """
        Calculate enrichment metrics from docking results.

        Returns:
            EnrichmentMetrics object
        """
        # Filter successful dockings
        successful = [r for r in self.results if r.success and r.binding_affinity is not None]

        if len(successful) == 0:
            logger.error("No successful docking results!")
            return EnrichmentMetrics(
                auc_roc=0.0,
                auc_pr=0.0,
                enrichment_factor_1=0.0,
                enrichment_factor_5=0.0,
                enrichment_factor_10=0.0,
                num_actives=len([r for r in self.results if r.is_active]),
                num_decoys=len([r for r in self.results if not r.is_active]),
                num_successful=0
            )

        # Extract labels and scores
        y_true = np.array([1 if r.is_active else 0 for r in successful])
        # Use negative affinity as score (more negative = better)
        y_score = np.array([-r.binding_affinity for r in successful])

        # Calculate AUC-ROC
        auc_roc = roc_auc_score(y_true, y_score)

        # Calculate AUC-PR (Precision-Recall)
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        auc_pr = auc(recall, precision)

        # Calculate enrichment factors
        ef_1 = self._calculate_enrichment_factor(y_true, y_score, 0.01)
        ef_5 = self._calculate_enrichment_factor(y_true, y_score, 0.05)
        ef_10 = self._calculate_enrichment_factor(y_true, y_score, 0.10)

        # Calculate mean affinities
        active_affinities = [r.binding_affinity for r in successful if r.is_active]
        decoy_affinities = [r.binding_affinity for r in successful if not r.is_active]

        mean_active_affinity = np.mean(active_affinities) if active_affinities else None
        mean_decoy_affinity = np.mean(decoy_affinities) if decoy_affinities else None

        # Calculate separation index
        separation_index = None
        if active_affinities and decoy_affinities:
            all_affinities = active_affinities + decoy_affinities
            std = np.std(all_affinities)
            if std > 0:
                separation_index = (mean_active_affinity - mean_decoy_affinity) / std

        return EnrichmentMetrics(
            auc_roc=auc_roc,
            auc_pr=auc_pr,
            enrichment_factor_1=ef_1,
            enrichment_factor_5=ef_5,
            enrichment_factor_10=ef_10,
            num_actives=len([r for r in successful if r.is_active]),
            num_decoys=len([r for r in successful if not r.is_active]),
            num_successful=len(successful),
            mean_active_affinity=mean_active_affinity,
            mean_decoy_affinity=mean_decoy_affinity,
            separation_index=separation_index
        )

    @staticmethod
    def _calculate_enrichment_factor(
        y_true: np.ndarray,
        y_score: np.ndarray,
        fraction: float
    ) -> float:
        """
        Calculate enrichment factor at a given fraction of the dataset.

        EF = (Actives_found / Actives_total) / (N_screened / N_total)

        Args:
            y_true: True labels (1=active, 0=decoy)
            y_score: Predicted scores (higher=better)
            fraction: Fraction of dataset to consider (e.g., 0.01 for top 1%)

        Returns:
            Enrichment factor
        """
        n_total = len(y_true)
        n_screened = max(1, int(n_total * fraction))
        n_actives_total = np.sum(y_true)

        if n_actives_total == 0:
            return 0.0

        # Get top-scoring compounds
        top_indices = np.argsort(y_score)[::-1][:n_screened]
        n_actives_found = np.sum(y_true[top_indices])

        # EF = (actives in top X%) / (expected actives by random)
        ef = (n_actives_found / n_actives_total) / (n_screened / n_total)

        return ef

    def _save_results(self, metrics: EnrichmentMetrics):
        """Save benchmark results to JSON file."""
        output_file = self.output_dir / f"dud_e_{self.target_name}_{self.adapter.name}.json"

        # Prepare result data
        result_data = {
            "benchmark": "DUD-E",
            "target": self.target_name,
            "adapter": self.adapter.name,
            "adapter_version": getattr(self.adapter, 'version', 'unknown'),
            "metrics": asdict(metrics),
            "results": [
                {
                    "smiles": r.smiles,
                    "is_active": r.is_active,
                    "binding_affinity": r.binding_affinity,
                    "binding_score": r.binding_score,
                    "success": r.success,
                    "error": r.error,
                    "runtime_ms": r.runtime_ms
                }
                for r in self.results
            ]
        }

        # Save to file
        with open(output_file, 'w') as f:
            json.dump(result_data, f, indent=2)

        logger.info(f"Results saved to: {output_file}")

        # Also save summary
        self._save_summary(metrics)

    def _save_summary(self, metrics: EnrichmentMetrics):
        """Save human-readable summary."""
        summary_file = self.output_dir / f"dud_e_{self.target_name}_{self.adapter.name}_summary.txt"

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 70 + "\n")
            f.write("DUD-E BENCHMARK RESULTS\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Target:          {self.target_name}\n")
            f.write(f"Adapter:         {self.adapter.name} v{getattr(self.adapter, 'version', 'unknown')}\n")
            f.write(f"Actives:         {metrics.num_actives}\n")
            f.write(f"Decoys:          {metrics.num_decoys}\n")
            f.write(f"Successful:      {metrics.num_successful} / {len(self.results)}\n")
            f.write("\n")
            f.write("-" * 70 + "\n")
            f.write("ENRICHMENT METRICS\n")
            f.write("-" * 70 + "\n\n")
            f.write(f"AUC-ROC:         {metrics.auc_roc:.4f}\n")
            f.write(f"AUC-PR:          {metrics.auc_pr:.4f}\n")
            f.write(f"EF 1%:           {metrics.enrichment_factor_1:.2f}x\n")
            f.write(f"EF 5%:           {metrics.enrichment_factor_5:.2f}x\n")
            f.write(f"EF 10%:          {metrics.enrichment_factor_10:.2f}x\n")
            f.write("\n")

            if metrics.mean_active_affinity is not None:
                f.write("-" * 70 + "\n")
                f.write("BINDING AFFINITY ANALYSIS\n")
                f.write("-" * 70 + "\n\n")
                f.write(f"Mean Active:     {metrics.mean_active_affinity:.2f} kcal/mol\n")
                f.write(f"Mean Decoy:      {metrics.mean_decoy_affinity:.2f} kcal/mol\n")
                f.write(f"Separation:      {metrics.separation_index:.2f} σ\n")

            f.write("\n")
            f.write("=" * 70 + "\n")
            f.write("INTERPRETATION\n")
            f.write("=" * 70 + "\n\n")
            f.write("AUC-ROC:  Random = 0.5, Perfect = 1.0\n")
            f.write("          >0.7 = Good, >0.8 = Excellent, >0.9 = Outstanding\n\n")
            f.write("EF:       Enrichment factor compared to random selection\n")
            f.write("          >1.0 = Better than random\n")
            f.write("          >5.0 = Good enrichment\n")
            f.write("          >10.0 = Excellent enrichment\n\n")

        logger.info(f"Summary saved to: {summary_file}")


# Example usage / testing code
async def example_dud_e_benchmark():
    """
    Example of running DUD-E benchmark.

    This is a mock example - in practice, you would:
    1. Load real actives/decoys from DUD-E dataset
    2. Configure a real docking adapter with target protein
    3. Run the benchmark
    """
    # Mock actives and decoys (in practice, load from DUD-E)
    actives = [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
    ]

    decoys = [
        "CCO",  # Ethanol
        "CCCCCC",  # Hexane
        "c1ccccc1",  # Benzene
        "CC(C)O",  # Isopropanol
        "CCC(C)O",  # 2-Butanol
    ]

    # Mock adapter (you would use real VinaAdapter, GNINAAdapter, etc.)
    class MockDockingAdapter:
        def __init__(self):
            self.name = "mock_docking"
            self.version = "1.0.0"

        async def execute(self, smiles):
            # Mock docking result
            import random
            affinity = random.uniform(-10, -5)

            from backend.core.adapters.protocol import AdapterResult
            return AdapterResult(
                success=True,
                data={
                    'binding_affinity': affinity,
                    'binding_score': (affinity + 10) / 5
                },
                error=None,
                metadata={'adapter_name': self.name}
            )

    # Run benchmark
    adapter = MockDockingAdapter()
    benchmark = DUDEBenchmark(
        adapter=adapter,
        target_name="EGFR",
        actives=actives,
        decoys=decoys,
        output_dir="benchmark_results"
    )

    metrics = await benchmark.run()

    print(f"\nBenchmark Results:")
    print(f"AUC-ROC: {metrics.auc_roc:.3f}")
    print(f"EF 1%: {metrics.enrichment_factor_1:.2f}x")
    print(f"EF 5%: {metrics.enrichment_factor_5:.2f}x")


if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Run example
    asyncio.run(example_dud_e_benchmark())
