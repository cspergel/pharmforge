"""
TDC ADMET Benchmark

Tests ADMET prediction accuracy against Therapeutics Data Commons (TDC) benchmarks.
Evaluates ML-based ADMET adapters on standard datasets with published baselines.

TDC provides curated ADMET datasets across multiple categories:
- Absorption: Caco2, HIA, Bioavailability, etc.
- Distribution: BBB, PPBR, VDss, etc.
- Metabolism: CYP inhibition, clearance, etc.
- Excretion: Half-life, clearance, etc.
- Toxicity: hERG, AMES, DILI, etc.

Reference:
    Huang et al. (2021) Nature Chemical Biology
    https://tdcommons.ai/
"""

import asyncio
import json
import logging
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, asdict
import numpy as np
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
    accuracy_score,
    roc_auc_score,
    precision_score,
    recall_score,
    f1_score
)

logger = logging.getLogger(__name__)


@dataclass
class PredictionResult:
    """Result of a single ADMET prediction"""
    smiles: str
    property_name: str
    true_value: Union[float, int]
    predicted_value: Optional[Union[float, int]] = None
    success: bool = False
    error: Optional[str] = None
    runtime_ms: Optional[float] = None


@dataclass
class PropertyMetrics:
    """Metrics for a single ADMET property"""
    property_name: str
    task_type: str  # "regression" or "classification"
    num_samples: int
    num_successful: int

    # Regression metrics (if applicable)
    mae: Optional[float] = None
    rmse: Optional[float] = None
    r2: Optional[float] = None

    # Classification metrics (if applicable)
    accuracy: Optional[float] = None
    roc_auc: Optional[float] = None
    precision: Optional[float] = None
    recall: Optional[float] = None
    f1: Optional[float] = None

    # Baseline comparison
    baseline_metric: Optional[float] = None
    improvement_over_baseline: Optional[float] = None


@dataclass
class BenchmarkSummary:
    """Overall benchmark summary"""
    adapter_name: str
    adapter_version: str
    num_properties: int
    total_samples: int
    total_successful: int
    mean_regression_r2: Optional[float] = None
    mean_classification_roc_auc: Optional[float] = None
    property_metrics: List[PropertyMetrics] = None


class TDCADMETBenchmark:
    """
    TDC ADMET benchmark for ADMET adapter validation.

    Tests prediction accuracy on standard TDC datasets with published baselines.
    """

    # TDC ADMET benchmarks with task types
    TDC_BENCHMARKS = {
        # Absorption (regression)
        'Caco2_Wang': {'type': 'regression', 'unit': 'log Papp (cm/s)', 'baseline_r2': 0.45},
        'HIA_Hou': {'type': 'classification', 'baseline_roc_auc': 0.98},
        'Bioavailability_Ma': {'type': 'classification', 'baseline_roc_auc': 0.70},

        # Distribution (regression)
        'BBB_Martins': {'type': 'classification', 'baseline_roc_auc': 0.91},
        'PPBR_AZ': {'type': 'regression', 'unit': '% bound', 'baseline_r2': 0.59},
        'VDss_Lombardo': {'type': 'regression', 'unit': 'log L/kg', 'baseline_r2': 0.54},

        # Metabolism (classification)
        'CYP2C9_Veith': {'type': 'classification', 'baseline_roc_auc': 0.79},
        'CYP2D6_Veith': {'type': 'classification', 'baseline_roc_auc': 0.75},
        'CYP3A4_Veith': {'type': 'classification', 'baseline_roc_auc': 0.85},
        'CYP2C19_Veith': {'type': 'classification', 'baseline_roc_auc': 0.78},
        'CYP1A2_Veith': {'type': 'classification', 'baseline_roc_auc': 0.83},

        # Excretion (regression)
        'Clearance_Hepatocyte_AZ': {'type': 'regression', 'unit': 'log mL/min/kg', 'baseline_r2': 0.43},
        'Half_Life_Obach': {'type': 'regression', 'unit': 'log h', 'baseline_r2': 0.35},

        # Toxicity (classification)
        'AMES': {'type': 'classification', 'baseline_roc_auc': 0.85},
        'hERG': {'type': 'classification', 'baseline_roc_auc': 0.82},
        'DILI': {'type': 'classification', 'baseline_roc_auc': 0.90},
        'Carcinogens_Laetz': {'type': 'classification', 'baseline_roc_auc': 0.71},
    }

    def __init__(
        self,
        adapter,
        properties: Optional[List[str]] = None,
        max_samples_per_property: Optional[int] = None,
        output_dir: str = "benchmark_results",
    ):
        """
        Initialize TDC ADMET benchmark.

        Args:
            adapter: ADMET adapter instance (e.g., ADMETaiAdapter)
            properties: List of properties to test (None = all available)
            max_samples_per_property: Max samples per property (None = all)
            output_dir: Directory to save benchmark results
        """
        self.adapter = adapter
        self.properties = properties or list(self.TDC_BENCHMARKS.keys())
        self.max_samples = max_samples_per_property
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.results: List[PredictionResult] = []
        self.tdc_data: Dict[str, Any] = {}

        logger.info(
            f"TDC ADMET Benchmark initialized: {len(self.properties)} properties"
        )

    async def run(self) -> BenchmarkSummary:
        """
        Run the benchmark: test all properties, calculate metrics.

        Returns:
            BenchmarkSummary with overall and per-property metrics
        """
        logger.info("Starting TDC ADMET benchmark...")
        start_time = time.time()

        # Load TDC data
        self._load_tdc_data()

        # Test each property
        property_metrics = []
        for prop in self.properties:
            if prop in self.tdc_data:
                logger.info(f"\nTesting property: {prop}")
                metrics = await self._test_property(prop)
                property_metrics.append(metrics)
            else:
                logger.warning(f"Property {prop} not available in TDC data")

        # Calculate summary metrics
        summary = self._calculate_summary(property_metrics)

        # Save results
        self._save_results(summary)

        elapsed_time = time.time() - start_time
        logger.info(f"\nBenchmark complete in {elapsed_time:.1f}s")

        return summary

    def _load_tdc_data(self):
        """
        Load TDC benchmark datasets.

        In practice, this would load real TDC data using the TDC library:
            from tdc.benchmark_group import admet_group
            group = admet_group(path='data/')
            benchmark = group.get('Caco2_Wang')
            train, test = benchmark['train'], benchmark['test']

        For this implementation, we use mock data as a template.
        """
        logger.info("Loading TDC benchmark data...")

        # Mock implementation - replace with real TDC data loading
        for prop in self.properties:
            if prop in self.TDC_BENCHMARKS:
                # In production, load real data:
                # self.tdc_data[prop] = self._load_tdc_property(prop)

                # Mock data for demonstration
                self.tdc_data[prop] = self._create_mock_data(prop)

        logger.info(f"Loaded {len(self.tdc_data)} TDC datasets")

    def _create_mock_data(self, property_name: str) -> Dict[str, Any]:
        """Create mock data for demonstration (replace with real TDC loading)."""
        import random
        from rdkit import Chem

        # Sample SMILES
        sample_smiles = [
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "CCO",  # Ethanol
            "c1ccccc1",  # Benzene
        ]

        task_type = self.TDC_BENCHMARKS[property_name]['type']

        samples = []
        for smiles in sample_smiles:
            if task_type == 'classification':
                value = random.choice([0, 1])
            else:  # regression
                value = random.uniform(-2, 2)

            samples.append({
                'smiles': smiles,
                'y': value
            })

        return {
            'property': property_name,
            'type': task_type,
            'samples': samples
        }

    def _load_tdc_property(self, property_name: str) -> Dict[str, Any]:
        """
        Load real TDC data for a property.

        This method should be implemented to use the real TDC library.
        """
        try:
            from tdc.single_pred import ADME, Tox

            # Map property names to TDC classes
            if property_name in ['AMES', 'hERG', 'DILI', 'Carcinogens_Laetz']:
                data = Tox(name=property_name)
            else:
                data = ADME(name=property_name)

            df = data.get_data()

            # Convert to our format
            samples = [
                {'smiles': row['Drug'], 'y': row['Y']}
                for _, row in df.iterrows()
            ]

            if self.max_samples:
                samples = samples[:self.max_samples]

            task_type = self.TDC_BENCHMARKS[property_name]['type']

            return {
                'property': property_name,
                'type': task_type,
                'samples': samples
            }

        except Exception as e:
            logger.error(f"Failed to load TDC data for {property_name}: {e}")
            return None

    async def _test_property(self, property_name: str) -> PropertyMetrics:
        """
        Test adapter on a single ADMET property.

        Args:
            property_name: Name of the property to test

        Returns:
            PropertyMetrics with evaluation results
        """
        prop_data = self.tdc_data[property_name]
        samples = prop_data['samples']
        task_type = prop_data['type']

        logger.info(f"  Testing {len(samples)} samples...")

        # Run predictions
        y_true = []
        y_pred = []
        successful = 0

        for i, sample in enumerate(samples):
            if (i + 1) % 10 == 0:
                logger.info(f"    Progress: {i + 1}/{len(samples)}")

            smiles = sample['smiles']
            true_value = sample['y']

            start_time = time.time()

            try:
                # Execute adapter
                result = await self.adapter.execute(smiles)

                runtime_ms = (time.time() - start_time) * 1000

                if result.success:
                    # Extract predicted value for this property
                    pred_value = self._extract_prediction(result, property_name)

                    if pred_value is not None:
                        y_true.append(true_value)
                        y_pred.append(pred_value)
                        successful += 1

                        self.results.append(PredictionResult(
                            smiles=smiles,
                            property_name=property_name,
                            true_value=true_value,
                            predicted_value=pred_value,
                            success=True,
                            runtime_ms=runtime_ms
                        ))
                    else:
                        self.results.append(PredictionResult(
                            smiles=smiles,
                            property_name=property_name,
                            true_value=true_value,
                            success=False,
                            error="Property not found in result",
                            runtime_ms=runtime_ms
                        ))
                else:
                    self.results.append(PredictionResult(
                        smiles=smiles,
                        property_name=property_name,
                        true_value=true_value,
                        success=False,
                        error=result.error,
                        runtime_ms=runtime_ms
                    ))

            except Exception as e:
                logger.error(f"Error predicting {smiles}: {e}")
                self.results.append(PredictionResult(
                    smiles=smiles,
                    property_name=property_name,
                    true_value=true_value,
                    success=False,
                    error=str(e)
                ))

        # Calculate metrics
        metrics = self._calculate_property_metrics(
            property_name,
            task_type,
            y_true,
            y_pred,
            len(samples),
            successful
        )

        logger.info(f"  Results: {successful}/{len(samples)} successful")
        if task_type == 'regression' and metrics.r2 is not None:
            logger.info(f"  R² = {metrics.r2:.3f}, MAE = {metrics.mae:.3f}")
        elif task_type == 'classification' and metrics.roc_auc is not None:
            logger.info(f"  ROC-AUC = {metrics.roc_auc:.3f}, Accuracy = {metrics.accuracy:.3f}")

        return metrics

    def _extract_prediction(
        self,
        result: Any,
        property_name: str
    ) -> Optional[Union[float, int]]:
        """
        Extract prediction for a specific property from adapter result.

        Args:
            result: AdapterResult from adapter.execute()
            property_name: Name of the property

        Returns:
            Predicted value or None if not found
        """
        if not result.data or 'properties' not in result.data:
            return None

        properties = result.data['properties']

        # Try exact match
        if property_name in properties:
            return properties[property_name]

        # Try case-insensitive match
        for key, value in properties.items():
            if key.lower() == property_name.lower():
                return value

        return None

    def _calculate_property_metrics(
        self,
        property_name: str,
        task_type: str,
        y_true: List[Union[float, int]],
        y_pred: List[Union[float, int]],
        num_samples: int,
        num_successful: int
    ) -> PropertyMetrics:
        """Calculate metrics for a single property."""
        if len(y_true) == 0 or len(y_pred) == 0:
            return PropertyMetrics(
                property_name=property_name,
                task_type=task_type,
                num_samples=num_samples,
                num_successful=0
            )

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        baseline_info = self.TDC_BENCHMARKS.get(property_name, {})

        if task_type == 'regression':
            mae = mean_absolute_error(y_true, y_pred)
            rmse = np.sqrt(mean_squared_error(y_true, y_pred))
            r2 = r2_score(y_true, y_pred)

            baseline_r2 = baseline_info.get('baseline_r2')
            improvement = None
            if baseline_r2 is not None:
                improvement = ((r2 - baseline_r2) / baseline_r2) * 100

            return PropertyMetrics(
                property_name=property_name,
                task_type=task_type,
                num_samples=num_samples,
                num_successful=num_successful,
                mae=mae,
                rmse=rmse,
                r2=r2,
                baseline_metric=baseline_r2,
                improvement_over_baseline=improvement
            )

        else:  # classification
            # Convert to binary if needed
            y_true_binary = (y_true > 0.5).astype(int)
            y_pred_binary = (y_pred > 0.5).astype(int)

            accuracy = accuracy_score(y_true_binary, y_pred_binary)

            # ROC-AUC (if we have probabilities)
            roc_auc = None
            try:
                roc_auc = roc_auc_score(y_true_binary, y_pred)
            except:
                pass

            precision = precision_score(y_true_binary, y_pred_binary, zero_division=0)
            recall = recall_score(y_true_binary, y_pred_binary, zero_division=0)
            f1 = f1_score(y_true_binary, y_pred_binary, zero_division=0)

            baseline_auc = baseline_info.get('baseline_roc_auc')
            improvement = None
            if baseline_auc is not None and roc_auc is not None:
                improvement = ((roc_auc - baseline_auc) / baseline_auc) * 100

            return PropertyMetrics(
                property_name=property_name,
                task_type=task_type,
                num_samples=num_samples,
                num_successful=num_successful,
                accuracy=accuracy,
                roc_auc=roc_auc,
                precision=precision,
                recall=recall,
                f1=f1,
                baseline_metric=baseline_auc,
                improvement_over_baseline=improvement
            )

    def _calculate_summary(
        self,
        property_metrics: List[PropertyMetrics]
    ) -> BenchmarkSummary:
        """Calculate overall benchmark summary."""
        total_samples = sum(m.num_samples for m in property_metrics)
        total_successful = sum(m.num_successful for m in property_metrics)

        # Calculate mean R² for regression tasks
        regression_r2s = [m.r2 for m in property_metrics if m.task_type == 'regression' and m.r2 is not None]
        mean_regression_r2 = np.mean(regression_r2s) if regression_r2s else None

        # Calculate mean ROC-AUC for classification tasks
        classification_aucs = [m.roc_auc for m in property_metrics if m.task_type == 'classification' and m.roc_auc is not None]
        mean_classification_roc_auc = np.mean(classification_aucs) if classification_aucs else None

        return BenchmarkSummary(
            adapter_name=self.adapter.name,
            adapter_version=getattr(self.adapter, 'version', 'unknown'),
            num_properties=len(property_metrics),
            total_samples=total_samples,
            total_successful=total_successful,
            mean_regression_r2=mean_regression_r2,
            mean_classification_roc_auc=mean_classification_roc_auc,
            property_metrics=property_metrics
        )

    def _save_results(self, summary: BenchmarkSummary):
        """Save benchmark results to JSON file."""
        output_file = self.output_dir / f"tdc_admet_{self.adapter.name}.json"

        # Prepare result data
        result_data = {
            "benchmark": "TDC ADMET",
            "adapter": self.adapter.name,
            "adapter_version": summary.adapter_version,
            "summary": {
                "num_properties": summary.num_properties,
                "total_samples": summary.total_samples,
                "total_successful": summary.total_successful,
                "mean_regression_r2": summary.mean_regression_r2,
                "mean_classification_roc_auc": summary.mean_classification_roc_auc
            },
            "property_metrics": [asdict(m) for m in summary.property_metrics],
            "results": [
                {
                    "smiles": r.smiles,
                    "property_name": r.property_name,
                    "true_value": r.true_value,
                    "predicted_value": r.predicted_value,
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

        logger.info(f"\nResults saved to: {output_file}")

        # Also save summary
        self._save_summary(summary)

    def _save_summary(self, summary: BenchmarkSummary):
        """Save human-readable summary."""
        summary_file = self.output_dir / f"tdc_admet_{self.adapter.name}_summary.txt"

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("TDC ADMET BENCHMARK RESULTS\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Adapter:         {summary.adapter_name} v{summary.adapter_version}\n")
            f.write(f"Properties:      {summary.num_properties}\n")
            f.write(f"Total Samples:   {summary.total_samples}\n")
            f.write(f"Successful:      {summary.total_successful} / {summary.total_samples} ")
            f.write(f"({100*summary.total_successful/summary.total_samples:.1f}%)\n")
            f.write("\n")

            if summary.mean_regression_r2 is not None:
                f.write(f"Mean R² (regression):          {summary.mean_regression_r2:.3f}\n")
            if summary.mean_classification_roc_auc is not None:
                f.write(f"Mean ROC-AUC (classification): {summary.mean_classification_roc_auc:.3f}\n")

            f.write("\n")
            f.write("-" * 80 + "\n")
            f.write("PER-PROPERTY RESULTS\n")
            f.write("-" * 80 + "\n\n")

            for m in summary.property_metrics:
                f.write(f"{m.property_name}:\n")
                f.write(f"  Type:       {m.task_type}\n")
                f.write(f"  Samples:    {m.num_successful} / {m.num_samples}\n")

                if m.task_type == 'regression':
                    if m.r2 is not None:
                        f.write(f"  R²:         {m.r2:.3f}")
                        if m.baseline_metric:
                            f.write(f" (baseline: {m.baseline_metric:.3f})")
                        f.write("\n")
                    if m.mae is not None:
                        f.write(f"  MAE:        {m.mae:.3f}\n")
                    if m.rmse is not None:
                        f.write(f"  RMSE:       {m.rmse:.3f}\n")
                else:  # classification
                    if m.roc_auc is not None:
                        f.write(f"  ROC-AUC:    {m.roc_auc:.3f}")
                        if m.baseline_metric:
                            f.write(f" (baseline: {m.baseline_metric:.3f})")
                        f.write("\n")
                    if m.accuracy is not None:
                        f.write(f"  Accuracy:   {m.accuracy:.3f}\n")
                    if m.f1 is not None:
                        f.write(f"  F1:         {m.f1:.3f}\n")

                if m.improvement_over_baseline is not None:
                    f.write(f"  Improvement: {m.improvement_over_baseline:+.1f}%\n")

                f.write("\n")

        logger.info(f"Summary saved to: {summary_file}")


# Example usage
async def example_tdc_benchmark():
    """Example of running TDC ADMET benchmark."""
    # Mock adapter
    class MockADMETAdapter:
        def __init__(self):
            self.name = "mock_admet"
            self.version = "1.0.0"

        async def execute(self, smiles):
            import random
            from backend.core.adapters.protocol import AdapterResult

            # Mock predictions
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

    # Run benchmark
    adapter = MockADMETAdapter()
    benchmark = TDCADMETBenchmark(
        adapter=adapter,
        properties=['AMES', 'hERG', 'Caco2_Wang'],
        output_dir="benchmark_results"
    )

    summary = await benchmark.run()

    print(f"\nBenchmark Results:")
    print(f"Properties tested: {summary.num_properties}")
    print(f"Mean Classification ROC-AUC: {summary.mean_classification_roc_auc:.3f}")
    print(f"Mean Regression R²: {summary.mean_regression_r2:.3f}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    asyncio.run(example_tdc_benchmark())
