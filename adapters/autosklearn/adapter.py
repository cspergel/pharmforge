"""
Auto-sklearn Adapter - Automated Machine Learning with Bayesian Optimization
Uses Bayesian optimization and meta-learning for automated ML pipeline optimization
"""
from typing import Any, Dict, Optional, List
import logging
import json
import pickle
import tempfile
from pathlib import Path

try:
    import numpy as np
    import autosklearn.classification
    import autosklearn.regression
    import autosklearn.metrics
    from sklearn.model_selection import cross_val_score
    AUTOSKLEARN_AVAILABLE = True
except ImportError:
    AUTOSKLEARN_AVAILABLE = False
    logging.warning("Auto-sklearn not available - install with: pip install auto-sklearn")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class AutoSklearnAdapter(AdapterProtocol):
    """
    Adapter for automated ML pipeline optimization using Auto-sklearn

    Uses Bayesian optimization (SMAC3) and meta-learning to automatically:
    - Select ML algorithms
    - Configure hyperparameters
    - Perform feature preprocessing
    - Build ensembles of best models

    Auto-sklearn differs from TPOT by:
    - Using Bayesian optimization instead of genetic programming
    - Leveraging meta-learning from past runs
    - Building weighted ensembles instead of single pipelines
    - Supporting warm-starting from previous runs

    Requirements:
        - Auto-sklearn: pip install auto-sklearn
        - Works on Linux/Mac (Windows support limited)
    """

    def __init__(self):
        super().__init__(
            name="autosklearn",
            adapter_type="local",
            config={
                "timeout": 7200,  # 2 hour default timeout
                "default_time_limit": 300,  # 5 minutes default
                "default_per_run_time_limit": 30,  # 30 seconds per model
                "default_ensemble_size": 50,
                "default_cv": 5,
                "memory_limit": 3072,  # MB (3GB default)
            }
        )
        self.version = "1.0.0"

        if not AUTOSKLEARN_AVAILABLE:
            logger.error("Auto-sklearn is not installed! Install with: pip install auto-sklearn")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input contains required data for Auto-sklearn optimization

        Args:
            input_data: Dictionary with training data and configuration

        Returns:
            True if valid, False otherwise
        """
        if not AUTOSKLEARN_AVAILABLE:
            return False

        if not isinstance(input_data, dict):
            logger.error("Input must be a dictionary")
            return False

        # Check required fields
        required_fields = ["X_train", "y_train"]
        for field in required_fields:
            if field not in input_data:
                logger.error(f"Missing required field: {field}")
                return False

        # Validate X_train
        X_train = input_data["X_train"]
        if not isinstance(X_train, (list, tuple, np.ndarray)):
            logger.error("X_train must be a list, tuple, or numpy array")
            return False

        if len(X_train) == 0:
            logger.error("X_train cannot be empty")
            return False

        # Validate y_train
        y_train = input_data["y_train"]
        if not isinstance(y_train, (list, tuple, np.ndarray)):
            logger.error("y_train must be a list, tuple, or numpy array")
            return False

        if len(y_train) == 0:
            logger.error("y_train cannot be empty")
            return False

        # Check that X_train and y_train have same length
        if len(X_train) != len(y_train):
            logger.error(f"X_train and y_train must have same length: {len(X_train)} vs {len(y_train)}")
            return False

        # Validate task type
        task = input_data.get("task", "classification")
        if task not in ["classification", "regression"]:
            logger.error(f"Invalid task type: {task}. Must be 'classification' or 'regression'")
            return False

        return True

    def _convert_to_numpy(self, data: Any) -> 'np.ndarray':
        """
        Convert input data to numpy array

        Args:
            data: Input data (list, tuple, or numpy array)

        Returns:
            Numpy array
        """
        if isinstance(data, np.ndarray):
            return data
        return np.array(data)

    def _get_metric(self, task: str, metric_name: Optional[str] = None):
        """
        Get Auto-sklearn metric object for the task

        Args:
            task: Task type (classification or regression)
            metric_name: User-specified metric name (optional)

        Returns:
            Auto-sklearn metric object
        """
        if metric_name is not None:
            # Map common metric names to Auto-sklearn metrics
            metric_map = {
                # Classification metrics
                "accuracy": autosklearn.metrics.accuracy,
                "roc_auc": autosklearn.metrics.roc_auc,
                "f1": autosklearn.metrics.f1,
                "precision": autosklearn.metrics.precision,
                "recall": autosklearn.metrics.recall,
                "log_loss": autosklearn.metrics.log_loss,
                # Regression metrics
                "r2": autosklearn.metrics.r2,
                "mse": autosklearn.metrics.mean_squared_error,
                "mae": autosklearn.metrics.mean_absolute_error,
                "rmse": autosklearn.metrics.root_mean_squared_error,
            }

            if metric_name in metric_map:
                return metric_map[metric_name]
            else:
                logger.warning(f"Unknown metric '{metric_name}', using default")

        # Default metrics
        if task == "classification":
            return autosklearn.metrics.accuracy
        else:  # regression
            return autosklearn.metrics.r2

    def _parse_include_preprocessors(self, include_list: Optional[List[str]]) -> Optional[Dict]:
        """
        Parse include_preprocessors list into Auto-sklearn format

        Args:
            include_list: List of preprocessor names to include

        Returns:
            Dictionary for Auto-sklearn include parameter or None
        """
        if not include_list:
            return None

        # Map common names to Auto-sklearn preprocessor names
        preprocessor_map = {
            "PCA": "feature_agglomeration",
            "StandardScaler": "rescaling:standardize",
            "MinMaxScaler": "rescaling:minmax",
            "Normalizer": "rescaling:normalize",
            "RobustScaler": "rescaling:robust_scaler",
            "QuantileTransformer": "rescaling:quantile_transformer",
            "PowerTransformer": "rescaling:power_transformer",
            "SelectKBest": "feature_preprocessor:select_percentile",
            "VarianceThreshold": "feature_preprocessor:select_rates",
        }

        # Create include dict
        include_dict = {}
        for item in include_list:
            if item in preprocessor_map:
                autosklearn_name = preprocessor_map[item]
                if ":" in autosklearn_name:
                    component_type, component_name = autosklearn_name.split(":")
                    if component_type not in include_dict:
                        include_dict[component_type] = []
                    include_dict[component_type].append(component_name)

        return include_dict if include_dict else None

    def _get_model_statistics(self, automl) -> Dict[str, Any]:
        """
        Extract statistics from fitted Auto-sklearn model

        Args:
            automl: Fitted Auto-sklearn instance

        Returns:
            Dictionary with model statistics
        """
        stats = {
            "total_models_evaluated": 0,
            "ensemble_size": 0,
            "ensemble_models": [],
            "best_model": None,
            "meta_features_used": False,
        }

        try:
            # Get ensemble information
            if hasattr(automl, "get_models_with_weights"):
                ensemble = automl.get_models_with_weights()
                stats["ensemble_size"] = len(ensemble)

                ensemble_models = []
                for weight, model in ensemble:
                    model_info = {
                        "weight": float(weight),
                        "model_id": str(model),
                        "algorithm": type(model).__name__,
                    }
                    ensemble_models.append(model_info)

                stats["ensemble_models"] = ensemble_models

                # Best model is the one with highest weight
                if ensemble_models:
                    stats["best_model"] = max(ensemble_models, key=lambda x: x["weight"])

            # Get run history
            if hasattr(automl, "cv_results_"):
                stats["total_models_evaluated"] = len(automl.cv_results_["status"])

            # Check if meta-learning was used
            if hasattr(automl, "automl_") and hasattr(automl.automl_, "initial_configurations"):
                stats["meta_features_used"] = True

        except Exception as e:
            logger.warning(f"Could not extract full model statistics: {e}")

        return stats

    def _get_leaderboard(self, automl) -> List[Dict[str, Any]]:
        """
        Extract leaderboard (top models) from Auto-sklearn run

        Args:
            automl: Fitted Auto-sklearn instance

        Returns:
            List of top models with scores
        """
        leaderboard = []

        try:
            if hasattr(automl, "cv_results_"):
                cv_results = automl.cv_results_

                # Get model scores and status
                scores = cv_results.get("mean_test_score", [])
                statuses = cv_results.get("status", [])

                # Only include successful runs
                valid_results = []
                for i, (score, status) in enumerate(zip(scores, statuses)):
                    if status == "Success":
                        valid_results.append({
                            "model_id": i,
                            "score": float(score),
                            "rank": 0,  # Will be set after sorting
                        })

                # Sort by score (descending)
                valid_results.sort(key=lambda x: x["score"], reverse=True)

                # Add ranks
                for rank, model in enumerate(valid_results[:20], 1):  # Top 20
                    model["rank"] = rank
                    leaderboard.append(model)

        except Exception as e:
            logger.warning(f"Could not extract leaderboard: {e}")

        return leaderboard

    def _export_model(self, automl, model_path: str) -> bool:
        """
        Export fitted model to disk using pickle

        Args:
            automl: Fitted Auto-sklearn instance
            model_path: Path to save model

        Returns:
            True if successful, False otherwise
        """
        try:
            with open(model_path, 'wb') as f:
                pickle.dump(automl, f)
            logger.info(f"Model exported to: {model_path}")
            return True
        except Exception as e:
            logger.error(f"Failed to export model: {e}")
            return False

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Auto-sklearn automated ML pipeline optimization

        Args:
            input_data: Dictionary containing:
                - X_train: Feature matrix (required)
                - y_train: Target values (required)
                - task: "classification" or "regression" (default: "classification")
                - time_limit: Total time budget in seconds (default: 300)
                - per_run_time_limit: Time limit per model in seconds (default: 30)
                - ensemble_size: Number of models in ensemble (default: 50)
                - metric: Metric to optimize (optional, auto-selected based on task)
                - include_preprocessors: List of preprocessor names to include (optional)
                - memory_limit: Memory limit in MB (default: 3072)
                - n_jobs: Number of parallel jobs (default: 1, -1 for all cores)
                - seed: Random seed for reproducibility (optional)
                - resampling_strategy: CV strategy, e.g., "cv" or "holdout" (default: "holdout-iterative-fit")
                - resampling_strategy_arguments: Dict with CV args (optional)
                - X_test: Test features for evaluation (optional)
                - y_test: Test targets for evaluation (optional)
                - export_model_path: Path to export fitted model (optional)
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing optimized model and results
        """
        # Check if Auto-sklearn is available
        if not AUTOSKLEARN_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="Auto-sklearn is not installed. Install with: pip install auto-sklearn"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input data for Auto-sklearn optimization"
            )

        try:
            # Extract configuration
            X_train = self._convert_to_numpy(input_data["X_train"])
            y_train = self._convert_to_numpy(input_data["y_train"])
            task = input_data.get("task", "classification")
            time_limit = input_data.get("time_limit", self.config["default_time_limit"])
            per_run_time_limit = input_data.get("per_run_time_limit", self.config["default_per_run_time_limit"])
            ensemble_size = input_data.get("ensemble_size", self.config["default_ensemble_size"])
            memory_limit = input_data.get("memory_limit", self.config["memory_limit"])
            n_jobs = input_data.get("n_jobs", 1)
            seed = input_data.get("seed", None)
            resampling_strategy = input_data.get("resampling_strategy", "holdout-iterative-fit")
            resampling_strategy_arguments = input_data.get("resampling_strategy_arguments", None)

            # Get metric
            metric = self._get_metric(task, input_data.get("metric", None))

            # Parse include preprocessors
            include_preprocessors = self._parse_include_preprocessors(
                input_data.get("include_preprocessors", None)
            )

            logger.info(f"Starting Auto-sklearn {task} optimization:")
            logger.info(f"  - Time limit: {time_limit}s")
            logger.info(f"  - Per-run time limit: {per_run_time_limit}s")
            logger.info(f"  - Ensemble size: {ensemble_size}")
            logger.info(f"  - Metric: {metric.name}")
            logger.info(f"  - Training samples: {len(X_train)}")
            logger.info(f"  - Features: {X_train.shape[1] if len(X_train.shape) > 1 else 1}")
            logger.info(f"  - Resampling strategy: {resampling_strategy}")

            # Create Auto-sklearn instance
            if task == "classification":
                automl = autosklearn.classification.AutoSklearnClassifier(
                    time_left_for_this_task=time_limit,
                    per_run_time_limit=per_run_time_limit,
                    ensemble_size=ensemble_size,
                    memory_limit=memory_limit,
                    n_jobs=n_jobs,
                    seed=seed,
                    resampling_strategy=resampling_strategy,
                    resampling_strategy_arguments=resampling_strategy_arguments,
                    metric=metric,
                    include=include_preprocessors,
                )
            else:  # regression
                automl = autosklearn.regression.AutoSklearnRegressor(
                    time_left_for_this_task=time_limit,
                    per_run_time_limit=per_run_time_limit,
                    ensemble_size=ensemble_size,
                    memory_limit=memory_limit,
                    n_jobs=n_jobs,
                    seed=seed,
                    resampling_strategy=resampling_strategy,
                    resampling_strategy_arguments=resampling_strategy_arguments,
                    metric=metric,
                    include=include_preprocessors,
                )

            # Fit Auto-sklearn
            logger.info("Fitting Auto-sklearn model (this may take a while)...")
            automl.fit(X_train, y_train)

            # Get training score
            train_score = automl.score(X_train, y_train)
            logger.info(f"Auto-sklearn optimization completed! Training score: {train_score:.6f}")

            # Get model statistics
            model_stats = self._get_model_statistics(automl)

            # Get leaderboard
            leaderboard = self._get_leaderboard(automl)

            # Get predictions on training set
            train_predictions = automl.predict(X_train).tolist()

            # Prepare result data
            result_data = {
                "model_id": f"autosklearn_{hash(str(automl))}",
                "task": task,
                "best_model": model_stats.get("best_model"),
                "ensemble": {
                    "size": model_stats.get("ensemble_size", 0),
                    "models": model_stats.get("ensemble_models", []),
                },
                "predictions": {
                    "train": train_predictions,
                },
                "performance": {
                    "train_score": float(train_score),
                    "metric": metric.name,
                },
                "leaderboard": leaderboard,
                "optimization_history": {
                    "total_models_evaluated": model_stats.get("total_models_evaluated", 0),
                    "meta_learning_used": model_stats.get("meta_features_used", False),
                },
                "configuration": {
                    "time_limit": time_limit,
                    "per_run_time_limit": per_run_time_limit,
                    "ensemble_size": ensemble_size,
                    "resampling_strategy": resampling_strategy,
                },
                "warnings": [],
            }

            # Evaluate on test set if provided
            if "X_test" in input_data and "y_test" in input_data:
                X_test = self._convert_to_numpy(input_data["X_test"])
                y_test = self._convert_to_numpy(input_data["y_test"])
                test_score = automl.score(X_test, y_test)
                test_predictions = automl.predict(X_test).tolist()

                result_data["predictions"]["test"] = test_predictions
                result_data["performance"]["test_score"] = float(test_score)

                logger.info(f"Test set score: {test_score:.6f}")

            # Get feature importance if available
            try:
                if hasattr(automl, "feature_importances_"):
                    result_data["feature_importance"] = automl.feature_importances_.tolist()
            except Exception as e:
                logger.debug(f"Could not extract feature importance: {e}")

            # Export model if requested
            export_path = input_data.get("export_model_path")
            if export_path:
                if self._export_model(automl, export_path):
                    result_data["model_path"] = export_path

            # Add warnings for common issues
            if model_stats.get("ensemble_size", 0) == 0:
                result_data["warnings"].append("No ensemble was built - single best model used")

            if model_stats.get("total_models_evaluated", 0) < 10:
                result_data["warnings"].append(
                    "Very few models evaluated - consider increasing time_limit"
                )

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "autosklearn",
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "task": task,
                    "time_limit": time_limit,
                    "ensemble_size": ensemble_size,
                    "metric": metric.name,
                }
            )

        except Exception as e:
            logger.error(f"Auto-sklearn optimization failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"Auto-sklearn optimization failed: {str(e)}"
            )
