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
    Execute a drug discovery pipeline asynchronously

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

        # Update status to RUNNING
        run.status = "running"
        db.commit()
        logger.info(f"Pipeline run {run_id} status updated to RUNNING")

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

        # Define adapter execution order
        adapter_sequence = [
            "rdkit_local",    # Local molecular properties
            "admet_ai",       # ADMET predictions
            "pubchem",        # PubChem data enrichment
            "chembl"          # Bioactivity data
        ]

        # Verify all required adapters are registered
        missing_adapters = []
        for adapter_name in adapter_sequence:
            if registry.get(adapter_name) is None:
                missing_adapters.append(adapter_name)

        if missing_adapters:
            error_msg = f"Required adapters not registered: {', '.join(missing_adapters)}"
            logger.error(error_msg)
            run.status = "failed"
            run.error_message = error_msg
            db.commit()
            return {"success": False, "error": error_msg}

        # Process each compound
        all_results = []
        total_compounds = len(input_smiles)

        for idx, smiles in enumerate(input_smiles):
            logger.info(f"Processing compound {idx + 1}/{total_compounds}: {smiles}")

            compound_data = {
                "smiles": smiles,
                "properties": None,
                "admet": None,
                "pubchem": None,
                "bioactivity": None
            }

            # Execute adapters sequentially
            for adapter_idx, adapter_name in enumerate(adapter_sequence):
                adapter = registry.get(adapter_name)

                try:
                    logger.info(f"Running adapter: {adapter_name} for {smiles}")

                    # Execute adapter (adapters are async but Celery workers handle this)
                    import asyncio
                    result = asyncio.run(adapter.execute(smiles))

                    if result.success:
                        # Store result in appropriate field
                        if adapter_name == "rdkit_local":
                            compound_data["properties"] = result.data
                        elif adapter_name == "admet_ai":
                            compound_data["admet"] = result.data
                        elif adapter_name == "pubchem":
                            compound_data["pubchem"] = result.data
                        elif adapter_name == "chembl":
                            compound_data["bioactivity"] = result.data

                        logger.info(f"✓ Adapter {adapter_name} completed successfully")
                    else:
                        logger.warning(f"✗ Adapter {adapter_name} failed: {result.error}")
                        compound_data[f"{adapter_name}_error"] = result.error

                except Exception as e:
                    error_msg = f"Exception in adapter {adapter_name}: {str(e)}"
                    logger.error(error_msg, exc_info=True)
                    compound_data[f"{adapter_name}_error"] = error_msg

                # Update progress after each adapter
                # Progress = (compound_index * num_adapters + adapter_index + 1) / (total_compounds * num_adapters)
                progress = ((idx * len(adapter_sequence)) + adapter_idx + 1) / (total_compounds * len(adapter_sequence))
                self.update_state(
                    state="PROGRESS",
                    meta={
                        "current": idx + 1,
                        "total": total_compounds,
                        "status": f"Processing {smiles} with {adapter_name}",
                        "progress": progress
                    }
                )
                logger.info(f"Progress: {progress:.1%}")

            # Create CompoundResult entry
            compound_result = CompoundResult(
                compound_id=f"{run_id}_{idx}",
                run_id=run.id,
                smiles=smiles,
                properties=compound_data.get("properties"),
                admet=compound_data.get("admet"),
                bioactivity=compound_data.get("bioactivity")
            )
            db.add(compound_result)

            all_results.append(compound_data)
            logger.info(f"Compound {idx + 1}/{total_compounds} completed")

        # Store aggregated results
        run.results = {
            "compounds": all_results,
            "total_processed": len(all_results),
            "adapter_sequence": adapter_sequence,
            "completed_at": datetime.utcnow().isoformat()
        }

        # Update status to COMPLETED
        run.status = "completed"
        run.completed_at = datetime.utcnow()
        db.commit()

        logger.info(f"Pipeline run {run_id} COMPLETED successfully")
        logger.info(f"Processed {len(all_results)} compounds")

        return {
            "success": True,
            "run_id": run_id,
            "total_compounds": len(all_results),
            "status": "completed"
        }

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
