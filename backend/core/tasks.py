"""
Celery Tasks for Pipeline Execution
Handles async execution of drug discovery pipelines
"""
from celery import Task
from sqlalchemy.orm import Session
from typing import Dict, Any, List
import logging
from datetime import datetime

from backend.celery_app import celery_app
from backend.db.database import SessionLocal
from backend.db.models import PipelineRun, CompoundResult
from backend.core.adapter_registry import register_all_adapters
from backend.core.adapters.protocol import registry

logger = logging.getLogger(__name__)


class DatabaseTask(Task):
    """
    Base task class that provides database session management
    All database tasks should inherit from this class
    """
    _db: Session = None

    def after_return(self, *args, **kwargs):
        """Close database session after task completes"""
        if self._db is not None:
            try:
                self._db.close()
            except Exception as e:
                logger.error(f"Error closing database session: {e}")
            finally:
                self._db = None

    @property
    def db(self) -> Session:
        """Get or create database session"""
        if self._db is None:
            self._db = SessionLocal()
        return self._db


@celery_app.task(bind=True, base=DatabaseTask, name="backend.core.tasks.execute_pipeline")
def execute_pipeline(self, run_id: str) -> Dict[str, Any]:
    """
    Execute a drug discovery pipeline asynchronously using the Pipeline orchestrator

    Args:
        run_id: Unique identifier for the pipeline run

    Returns:
        Dictionary containing execution results and status
    """
    logger.info(f"Starting pipeline execution for run_id: {run_id}")

    # Get database session from task
    db: Session = self.db

    try:
        # Ensure adapters are registered
        register_all_adapters()

        # Retrieve pipeline run from database
        run = db.query(PipelineRun).filter(PipelineRun.run_id == run_id).first()

        if not run:
            error_msg = f"Pipeline run not found: {run_id}"
            logger.error(error_msg)
            return {"success": False, "error": error_msg}

        # Get input SMILES
        input_smiles: List[str] = run.input_smiles
        if not input_smiles:
            error_msg = "No input SMILES provided"
            logger.error(error_msg)
            run.status = "failed"
            run.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}

        logger.info(f"Processing {len(input_smiles)} SMILES strings")

        # Import Pipeline orchestrator
        from backend.core.pipeline import create_pipeline, PipelineConfig

        # Create pipeline configuration from run config
        config = PipelineConfig()
        if run.pipeline_config:
            # Override defaults with user config
            if "adapter_sequence" in run.pipeline_config:
                config.adapter_sequence = run.pipeline_config["adapter_sequence"]
            if "ranking_method" in run.pipeline_config:
                config.ranking_method = run.pipeline_config["ranking_method"]
            if "ranking_weights" in run.pipeline_config:
                config.ranking_weights = run.pipeline_config["ranking_weights"]
            if "n_top_candidates" in run.pipeline_config:
                config.n_top_candidates = run.pipeline_config["n_top_candidates"]

        # Create pipeline
        pipeline = create_pipeline(db=db, config=config)

        # Execute pipeline (runs async)
        import asyncio
        result = asyncio.run(pipeline.execute(input_smiles, run_id))

        logger.info(f"Pipeline execution completed: {result}")

        return result

    except Exception as e:
        error_msg = f"Pipeline execution failed: {str(e)}"
        logger.error(error_msg, exc_info=True)

        # Update run status to FAILED
        try:
            run = db.query(PipelineRun).filter(PipelineRun.run_id == run_id).first()
            if run:
                run.status = "failed"
                run.error_message = error_msg
                db.commit()
        except Exception as db_error:
            logger.error(f"Failed to update run status: {db_error}")

        return {
            "success": False,
            "run_id": run_id,
            "error": error_msg
        }
