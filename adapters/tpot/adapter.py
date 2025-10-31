"""
TPOT Adapter - Tree-based Pipeline Optimization Tool
Automated ML pipeline optimization using genetic programming
"""
from typing import Any, Dict, Optional, List
import logging
import io
import json

try:
    import numpy as np
    from tpot import TPOTClassifier, TPOTRegressor
    from sklearn.model_selection import cross_val_score
    TPOT_AVAILABLE = True
except ImportError:
    TPOT_AVAILABLE = False
    logging.warning("TPOT not available - install with: pip install tpot")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class TPOTAdapter(AdapterProtocol):
    """
    Adapter for automated ML pipeline optimization using TPOT
    Uses genetic programming to evolve optimal ML pipelines
    """

    def __init__(self):
        super().__init__(
            name="tpot",
            adapter_type="local",
            config={
                "timeout": 7200,  # 2 hour default timeout for optimization
                "default_generations": 5,
                "default_population_size": 20,
                "default_cv": 5,
                "max_time_mins": 60,
                "max_eval_time_mins": 5
            }
        )
        self.version = "1.0.0"

        if not TPOT_AVAILABLE:
            logger.error("TPOT is not installed! Install with: pip install tpot")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input contains required data for TPOT optimization

        Args:
            input_data: Dictionary with training data and configuration

        Returns:
            True if valid, False otherwise
        """
        if not TPOT_AVAILABLE:
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
        if not isinstance(X_train, (list, tuple)):
            logger.error("X_train must be a list or tuple")
            return False

        if len(X_train) == 0:
            logger.error("X_train cannot be empty")
            return False

        # Validate y_train
        y_train = input_data["y_train"]
        if not isinstance(y_train, (list, tuple)):
            logger.error("y_train must be a list or tuple")
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

    def _get_scoring_metric(self, task: str, scoring: Optional[str] = None) -> str:
        """
        Get appropriate scoring metric for the task

        Args:
            task: Task type (classification or regression)
            scoring: User-specified scoring metric (optional)

        Returns:
            Scoring metric string
        """
        if scoring is not None:
            return scoring

        # Default scoring metrics
        if task == "classification":
            return "roc_auc"
        else:  # regression
            return "neg_mean_squared_error"

    def _create_tpot_instance(
        self,
        task: str,
        generations: int,
        population_size: int,
        cv: int,
        scoring: str,
        random_state: Optional[int] = None,
        verbosity: int = 2,
        max_time_mins: Optional[int] = None,
        max_eval_time_mins: Optional[int] = None,
        n_jobs: int = -1,
        config_dict: Optional[str] = None
    ):
        """
        Create TPOT instance based on task type

        Args:
            task: Task type (classification or regression)
            generations: Number of generations to run
            population_size: Population size for genetic algorithm
            cv: Number of cross-validation folds
            scoring: Scoring metric
            random_state: Random seed for reproducibility
            verbosity: Verbosity level (0-3)
            max_time_mins: Maximum time in minutes
            max_eval_time_mins: Maximum evaluation time per pipeline
            n_jobs: Number of parallel jobs (-1 for all cores)
            config_dict: Configuration dictionary name (e.g., 'TPOT light', 'TPOT MDR', 'TPOT sparse')

        Returns:
            TPOT instance (TPOTClassifier or TPOTRegressor)
        """
        common_params = {
            "generations": generations,
            "population_size": population_size,
            "cv": cv,
            "scoring": scoring,
            "random_state": random_state,
            "verbosity": verbosity,
            "n_jobs": n_jobs,
        }

        # Add time limits if specified
        if max_time_mins is not None:
            common_params["max_time_mins"] = max_time_mins
        if max_eval_time_mins is not None:
            common_params["max_eval_time_mins"] = max_eval_time_mins

        # Add config dict if specified
        if config_dict is not None:
            common_params["config_dict"] = config_dict

        if task == "classification":
            return TPOTClassifier(**common_params)
        else:  # regression
            return TPOTRegressor(**common_params)

    def _extract_pipeline_info(self, tpot_instance) -> Dict[str, Any]:
        """
        Extract information about the fitted pipeline

        Args:
            tpot_instance: Fitted TPOT instance

        Returns:
            Dictionary with pipeline information
        """
        pipeline_info = {
            "fitted_pipeline": tpot_instance.fitted_pipeline_,
            "pipeline_steps": []
        }

        # Extract pipeline steps
        if hasattr(tpot_instance.fitted_pipeline_, "steps"):
            for step_name, step_obj in tpot_instance.fitted_pipeline_.steps:
                step_info = {
                    "name": step_name,
                    "type": type(step_obj).__name__,
                    "parameters": step_obj.get_params() if hasattr(step_obj, "get_params") else {}
                }
                pipeline_info["pipeline_steps"].append(step_info)

        return pipeline_info

    def _export_pipeline_code(self, tpot_instance) -> str:
        """
        Export the optimized pipeline as Python code

        Args:
            tpot_instance: Fitted TPOT instance

        Returns:
            Python code as string
        """
        try:
            # Create a string buffer to capture the exported code
            code_buffer = io.StringIO()
            tpot_instance.export(code_buffer)
            code = code_buffer.getvalue()
            code_buffer.close()
            return code
        except Exception as e:
            logger.error(f"Failed to export pipeline code: {e}")
            return f"# Error exporting code: {str(e)}"

    def _get_optimization_history(self, tpot_instance) -> List[Dict[str, Any]]:
        """
        Extract optimization history from TPOT

        Args:
            tpot_instance: Fitted TPOT instance

        Returns:
            List of generation statistics
        """
        history = []

        try:
            # TPOT stores evaluation history
            if hasattr(tpot_instance, "evaluated_individuals_"):
                eval_individuals = tpot_instance.evaluated_individuals_

                # Group by generation
                generations_data = {}
                for pipeline, scores in eval_individuals.items():
                    if isinstance(scores, dict) and "generation" in scores:
                        gen = scores["generation"]
                        if gen not in generations_data:
                            generations_data[gen] = []
                        generations_data[gen].append({
                            "pipeline": pipeline,
                            "score": scores.get("internal_cv_score", None),
                            "operator_count": scores.get("operator_count", None),
                        })

                # Format as list
                for gen, pipelines in sorted(generations_data.items()):
                    scores = [p["score"] for p in pipelines if p["score"] is not None]
                    if scores:
                        history.append({
                            "generation": gen,
                            "pipelines_evaluated": len(pipelines),
                            "best_score": max(scores),
                            "mean_score": np.mean(scores),
                            "std_score": np.std(scores),
                        })

        except Exception as e:
            logger.warning(f"Could not extract optimization history: {e}")

        return history

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute TPOT automated ML pipeline optimization

        Args:
            input_data: Dictionary containing:
                - X_train: Feature matrix (required)
                - y_train: Target values (required)
                - task: "classification" or "regression" (default: "classification")
                - generations: Number of generations (default: 5)
                - population_size: Population size (default: 20)
                - scoring: Scoring metric (optional, auto-selected based on task)
                - cv: Cross-validation folds (default: 5)
                - random_state: Random seed (optional)
                - verbosity: Verbosity level 0-3 (default: 2)
                - max_time_mins: Maximum time in minutes (optional)
                - max_eval_time_mins: Max evaluation time per pipeline (optional)
                - n_jobs: Number of parallel jobs (default: -1)
                - config_dict: Configuration preset (e.g., 'TPOT light')
                - X_test: Test features for evaluation (optional)
                - y_test: Test targets for evaluation (optional)
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing optimized pipeline and results
        """
        # Check if TPOT is available
        if not TPOT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="TPOT is not installed. Install with: pip install tpot"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input data for TPOT optimization"
            )

        try:
            # Extract configuration
            X_train = self._convert_to_numpy(input_data["X_train"])
            y_train = self._convert_to_numpy(input_data["y_train"])
            task = input_data.get("task", "classification")
            generations = input_data.get("generations", self.config["default_generations"])
            population_size = input_data.get("population_size", self.config["default_population_size"])
            cv = input_data.get("cv", self.config["default_cv"])
            random_state = input_data.get("random_state", None)
            verbosity = input_data.get("verbosity", 2)
            max_time_mins = input_data.get("max_time_mins", None)
            max_eval_time_mins = input_data.get("max_eval_time_mins", None)
            n_jobs = input_data.get("n_jobs", -1)
            config_dict = input_data.get("config_dict", None)

            # Get scoring metric
            scoring = self._get_scoring_metric(task, input_data.get("scoring", None))

            logger.info(f"Starting TPOT {task} optimization:")
            logger.info(f"  - Generations: {generations}")
            logger.info(f"  - Population size: {population_size}")
            logger.info(f"  - CV folds: {cv}")
            logger.info(f"  - Scoring: {scoring}")
            logger.info(f"  - Training samples: {len(X_train)}")
            logger.info(f"  - Features: {X_train.shape[1] if len(X_train.shape) > 1 else 1}")

            # Create TPOT instance
            tpot = self._create_tpot_instance(
                task=task,
                generations=generations,
                population_size=population_size,
                cv=cv,
                scoring=scoring,
                random_state=random_state,
                verbosity=verbosity,
                max_time_mins=max_time_mins,
                max_eval_time_mins=max_eval_time_mins,
                n_jobs=n_jobs,
                config_dict=config_dict
            )

            # Fit TPOT
            logger.info("Fitting TPOT model (this may take a while)...")
            tpot.fit(X_train, y_train)

            # Get best pipeline
            best_pipeline = tpot.fitted_pipeline_
            best_score = tpot.score(X_train, y_train)

            logger.info(f"TPOT optimization completed! Best CV score: {best_score:.6f}")

            # Extract pipeline information
            pipeline_info = self._extract_pipeline_info(tpot)

            # Export pipeline code
            pipeline_code = self._export_pipeline_code(tpot)

            # Get optimization history
            optimization_history = self._get_optimization_history(tpot)

            # Prepare result data
            result_data = {
                "best_pipeline": str(best_pipeline),
                "best_score": float(best_score),
                "pipeline_code": pipeline_code,
                "fitted_pipeline": best_pipeline,
                "pipeline_steps": pipeline_info["pipeline_steps"],
                "optimization_history": optimization_history,
                "metadata": {
                    "task": task,
                    "generations_evaluated": generations,
                    "population_size": population_size,
                    "scoring": scoring,
                    "cv_folds": cv,
                    "n_samples": len(X_train),
                    "n_features": X_train.shape[1] if len(X_train.shape) > 1 else 1,
                }
            }

            # Evaluate on test set if provided
            if "X_test" in input_data and "y_test" in input_data:
                X_test = self._convert_to_numpy(input_data["X_test"])
                y_test = self._convert_to_numpy(input_data["y_test"])
                test_score = tpot.score(X_test, y_test)
                result_data["test_score"] = float(test_score)
                result_data["metadata"]["n_test_samples"] = len(X_test)
                logger.info(f"Test set score: {test_score:.6f}")

            # Get best generation info
            if optimization_history:
                best_gen = max(optimization_history, key=lambda x: x["best_score"])
                result_data["metadata"]["best_generation"] = best_gen["generation"]
                result_data["metadata"]["best_generation_score"] = best_gen["best_score"]
                result_data["metadata"]["pipelines_tested"] = sum(
                    h["pipelines_evaluated"] for h in optimization_history
                )

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "tpot",
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "task": task,
                    "generations": generations,
                    "population_size": population_size,
                }
            )

        except Exception as e:
            logger.error(f"TPOT optimization failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"TPOT optimization failed: {str(e)}"
            )
