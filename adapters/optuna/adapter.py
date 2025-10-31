"""
Optuna Adapter - Hyperparameter optimization framework
Provides automated hyperparameter tuning for other adapters and models
"""
from typing import Any, Dict, Optional, List, Callable
import logging
import json

try:
    import optuna
    from optuna.samplers import TPESampler, RandomSampler, GridSampler
    from optuna.pruners import MedianPruner
    OPTUNA_AVAILABLE = True
except ImportError:
    OPTUNA_AVAILABLE = False
    logging.warning("Optuna not available - install with: pip install optuna")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class OptunaAdapter(AdapterProtocol):
    """
    Adapter for hyperparameter optimization using Optuna
    Supports various sampling strategies and optimization algorithms
    """

    def __init__(self):
        super().__init__(
            name="optuna",
            adapter_type="local",
            config={
                "timeout": 3600,  # 1 hour default timeout for optimization
                "default_n_trials": 100,
                "default_sampler": "TPE"
            }
        )
        self.version = "1.0.0"

        if not OPTUNA_AVAILABLE:
            logger.error("Optuna is not installed! Install with: pip install optuna")

        # Store available samplers
        self.available_samplers = {
            "TPE": TPESampler if OPTUNA_AVAILABLE else None,
            "Random": RandomSampler if OPTUNA_AVAILABLE else None,
            "Grid": GridSampler if OPTUNA_AVAILABLE else None,
        }

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input contains required optimization configuration

        Args:
            input_data: Dictionary with optimization configuration

        Returns:
            True if valid, False otherwise
        """
        if not OPTUNA_AVAILABLE:
            return False

        if not isinstance(input_data, dict):
            logger.error("Input must be a dictionary")
            return False

        # Check required fields
        required_fields = ["search_space", "objective_function"]
        for field in required_fields:
            if field not in input_data:
                logger.error(f"Missing required field: {field}")
                return False

        # Validate search_space
        if not isinstance(input_data["search_space"], dict):
            logger.error("search_space must be a dictionary")
            return False

        if len(input_data["search_space"]) == 0:
            logger.error("search_space cannot be empty")
            return False

        # Validate search space parameters
        for param_name, param_config in input_data["search_space"].items():
            if not isinstance(param_config, dict):
                logger.error(f"Parameter config for {param_name} must be a dictionary")
                return False

            if "type" not in param_config:
                logger.error(f"Parameter {param_name} missing 'type' field")
                return False

            param_type = param_config["type"]
            if param_type not in ["float", "int", "categorical"]:
                logger.error(f"Invalid parameter type: {param_type}")
                return False

            # Validate based on type
            if param_type in ["float", "int"]:
                if "low" not in param_config or "high" not in param_config:
                    logger.error(f"Parameter {param_name} missing 'low' or 'high'")
                    return False
            elif param_type == "categorical":
                if "choices" not in param_config:
                    logger.error(f"Parameter {param_name} missing 'choices'")
                    return False
                if not isinstance(param_config["choices"], list):
                    logger.error(f"Parameter {param_name} choices must be a list")
                    return False

        return True

    def _create_sampler(self, sampler_name: str, sampler_kwargs: Optional[Dict] = None) -> Any:
        """
        Create an Optuna sampler based on name

        Args:
            sampler_name: Name of the sampler (TPE, Random, Grid)
            sampler_kwargs: Additional kwargs for the sampler

        Returns:
            Optuna sampler instance
        """
        sampler_kwargs = sampler_kwargs or {}

        if sampler_name not in self.available_samplers:
            logger.warning(f"Unknown sampler {sampler_name}, using TPE")
            sampler_name = "TPE"

        sampler_class = self.available_samplers[sampler_name]

        try:
            if sampler_name == "TPE":
                # TPE sampler with sensible defaults
                return sampler_class(
                    n_startup_trials=sampler_kwargs.get("n_startup_trials", 10),
                    n_ei_candidates=sampler_kwargs.get("n_ei_candidates", 24),
                    seed=sampler_kwargs.get("seed", None)
                )
            elif sampler_name == "Random":
                return sampler_class(seed=sampler_kwargs.get("seed", None))
            elif sampler_name == "Grid":
                # Grid sampler requires search_space
                if "search_space" not in sampler_kwargs:
                    logger.error("Grid sampler requires search_space in sampler_kwargs")
                    return TPESampler()  # Fallback to TPE
                return sampler_class(sampler_kwargs["search_space"])
            else:
                return TPESampler()  # Default fallback
        except Exception as e:
            logger.error(f"Error creating sampler {sampler_name}: {e}")
            return TPESampler()  # Fallback

    def _create_objective_function(
        self,
        objective_config: Any,
        search_space: Dict[str, Dict],
        custom_objective: Optional[Callable] = None
    ) -> Callable:
        """
        Create an objective function for Optuna

        Args:
            objective_config: Configuration for the objective
            search_space: Parameter search space configuration
            custom_objective: Optional custom objective function

        Returns:
            Objective function compatible with Optuna
        """
        def objective(trial: 'optuna.Trial') -> float:
            """
            Objective function for Optuna optimization

            Args:
                trial: Optuna trial object

            Returns:
                Objective value to minimize/maximize
            """
            # Suggest parameters based on search space
            params = {}
            for param_name, param_config in search_space.items():
                param_type = param_config["type"]

                if param_type == "float":
                    params[param_name] = trial.suggest_float(
                        param_name,
                        param_config["low"],
                        param_config["high"],
                        log=param_config.get("log", False)
                    )
                elif param_type == "int":
                    params[param_name] = trial.suggest_int(
                        param_name,
                        param_config["low"],
                        param_config["high"],
                        log=param_config.get("log", False)
                    )
                elif param_type == "categorical":
                    params[param_name] = trial.suggest_categorical(
                        param_name,
                        param_config["choices"]
                    )

            # Call custom objective if provided
            if custom_objective is not None:
                return custom_objective(params, trial)

            # Otherwise, use the objective_config
            if isinstance(objective_config, dict) and "function" in objective_config:
                # This is a placeholder - in real use, you would call the actual function
                # For now, we return a dummy value
                logger.warning("No custom objective provided, using dummy objective")
                # Simple dummy objective: sum of squared parameters
                return sum(v**2 if isinstance(v, (int, float)) else 0 for v in params.values())

            # Default fallback
            logger.error("Invalid objective configuration")
            return float('inf')

        return objective

    def _extract_trial_history(self, study: 'optuna.Study') -> List[Dict[str, Any]]:
        """
        Extract trial history from Optuna study

        Args:
            study: Completed Optuna study

        Returns:
            List of trial information dictionaries
        """
        trials = []
        for trial in study.trials:
            trial_info = {
                "number": trial.number,
                "value": trial.value,
                "params": trial.params,
                "state": trial.state.name,
                "datetime_start": trial.datetime_start.isoformat() if trial.datetime_start else None,
                "datetime_complete": trial.datetime_complete.isoformat() if trial.datetime_complete else None,
                "duration": trial.duration.total_seconds() if trial.duration else None,
            }
            trials.append(trial_info)

        return trials

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Optuna hyperparameter optimization

        Args:
            input_data: Dictionary containing:
                - search_space: Dict defining parameter search space
                - objective_function: Function to optimize or config dict
                - n_trials: Number of optimization trials (optional)
                - sampler: Sampler type (TPE, Random, Grid) (optional)
                - direction: "minimize" or "maximize" (optional, default: "minimize")
                - study_name: Name for the study (optional)
                - sampler_kwargs: Additional sampler arguments (optional)
            **kwargs: Additional parameters including:
                - custom_objective: Custom objective function (optional)
                - pruner: Pruner configuration (optional)

        Returns:
            AdapterResult containing best parameters and optimization history
        """
        # Check if Optuna is available
        if not OPTUNA_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="Optuna is not installed. Install with: pip install optuna"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid optimization configuration"
            )

        try:
            # Extract configuration
            search_space = input_data["search_space"]
            objective_config = input_data["objective_function"]
            n_trials = input_data.get("n_trials", self.config["default_n_trials"])
            sampler_name = input_data.get("sampler", self.config["default_sampler"])
            direction = input_data.get("direction", "minimize")
            study_name = input_data.get("study_name", "optuna_optimization")
            sampler_kwargs = input_data.get("sampler_kwargs", {})

            # Get custom objective from kwargs if provided
            custom_objective = kwargs.get("custom_objective", None)

            # Create sampler
            sampler = self._create_sampler(sampler_name, sampler_kwargs)

            # Create pruner if specified
            pruner = None
            if "pruner" in kwargs and kwargs["pruner"]:
                pruner = MedianPruner()

            # Create study
            study = optuna.create_study(
                direction=direction,
                sampler=sampler,
                pruner=pruner,
                study_name=study_name
            )

            # Create objective function
            objective_fn = self._create_objective_function(
                objective_config,
                search_space,
                custom_objective
            )

            # Run optimization
            logger.info(f"Starting Optuna optimization with {n_trials} trials using {sampler_name} sampler")
            study.optimize(objective_fn, n_trials=n_trials)

            # Extract results
            best_params = study.best_params
            best_value = study.best_value
            best_trial = study.best_trial.number

            # Get trial history
            trial_history = self._extract_trial_history(study)

            # Prepare result data
            result_data = {
                "best_params": best_params,
                "best_value": best_value,
                "best_trial": best_trial,
                "n_trials": len(study.trials),
                "trial_history": trial_history,
                "sampler": sampler_name,
                "direction": direction,
            }

            logger.info(f"Optimization completed: best_value={best_value:.6f} at trial {best_trial}")

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "optuna",
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "study_name": study_name,
                    "n_trials": n_trials,
                    "sampler": sampler_name,
                }
            )

        except Exception as e:
            logger.error(f"Optuna optimization failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"Optimization failed: {str(e)}"
            )
