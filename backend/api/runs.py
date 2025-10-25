"""
Pipeline Runs API Endpoints
Handles creation, retrieval, and management of pipeline runs
"""
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from pydantic import BaseModel, Field
import logging
import uuid
from datetime import datetime

from backend.db.database import get_db
from backend.db.models import PipelineRun, CompoundResult
from backend.core.tasks import execute_pipeline

logger = logging.getLogger(__name__)

router = APIRouter()


# Pydantic models for request/response
class CreateRunRequest(BaseModel):
    """Request model for creating a new pipeline run"""
    input_smiles: List[str] = Field(..., description="List of SMILES strings to process", min_items=1)
    pipeline_config: Optional[dict] = Field(None, description="Optional pipeline configuration")
    user_id: Optional[str] = Field(None, description="Optional user identifier")

    class Config:
        json_schema_extra = {
            "example": {
                "input_smiles": ["CCO", "CC(=O)O", "c1ccccc1"],
                "pipeline_config": {"adapters": ["rdkit_local", "admet_ai"]},
                "user_id": "user123"
            }
        }


class RunResponse(BaseModel):
    """Response model for pipeline run information"""
    run_id: str
    status: str
    created_at: datetime
    updated_at: Optional[datetime]
    completed_at: Optional[datetime]
    input_smiles: List[str]
    pipeline_config: Optional[dict]
    results: Optional[dict]
    error_message: Optional[str]
    user_id: Optional[str]

    class Config:
        from_attributes = True


class RunListResponse(BaseModel):
    """Response model for list of runs"""
    runs: List[RunResponse]
    total: int
    page: int
    page_size: int


@router.post("/runs", response_model=RunResponse, status_code=201)
def create_run(
    request: CreateRunRequest,
    db: Session = Depends(get_db)
):
    """
    Create a new pipeline run and trigger async execution

    Args:
        request: CreateRunRequest containing input SMILES and config
        db: Database session

    Returns:
        RunResponse with run_id and initial status
    """
    logger.info(f"Creating new pipeline run with {len(request.input_smiles)} SMILES")

    # Generate unique run_id
    run_id = f"run_{uuid.uuid4().hex[:16]}"

    # Create pipeline run entry
    pipeline_run = PipelineRun(
        run_id=run_id,
        status="pending",
        input_smiles=request.input_smiles,
        pipeline_config=request.pipeline_config,
        user_id=request.user_id
    )

    # Save to database
    db.add(pipeline_run)
    db.commit()
    db.refresh(pipeline_run)

    logger.info(f"Created pipeline run: {run_id}")

    # Trigger async execution via Celery
    try:
        task = execute_pipeline.delay(run_id)
        logger.info(f"Triggered Celery task {task.id} for run {run_id}")
    except Exception as e:
        logger.error(f"Failed to trigger Celery task: {e}")
        pipeline_run.status = "failed"
        pipeline_run.error_message = f"Failed to start pipeline execution: {str(e)}"
        db.commit()
        db.refresh(pipeline_run)

    # Return response
    return RunResponse(
        run_id=pipeline_run.run_id,
        status=pipeline_run.status,
        created_at=pipeline_run.created_at,
        updated_at=pipeline_run.updated_at,
        completed_at=pipeline_run.completed_at,
        input_smiles=pipeline_run.input_smiles,
        pipeline_config=pipeline_run.pipeline_config,
        results=pipeline_run.results,
        error_message=pipeline_run.error_message,
        user_id=pipeline_run.user_id
    )


@router.get("/runs/{run_id}", response_model=RunResponse)
def get_run(
    run_id: str,
    db: Session = Depends(get_db)
):
    """
    Get status and results for a specific pipeline run

    Args:
        run_id: Unique run identifier
        db: Database session

    Returns:
        RunResponse with current status and results

    Raises:
        HTTPException: 404 if run not found
    """
    logger.info(f"Retrieving run: {run_id}")

    # Query database
    pipeline_run = db.query(PipelineRun).filter(PipelineRun.run_id == run_id).first()

    if not pipeline_run:
        logger.warning(f"Run not found: {run_id}")
        raise HTTPException(status_code=404, detail=f"Run not found: {run_id}")

    # Return response
    return RunResponse(
        run_id=pipeline_run.run_id,
        status=pipeline_run.status,
        created_at=pipeline_run.created_at,
        updated_at=pipeline_run.updated_at,
        completed_at=pipeline_run.completed_at,
        input_smiles=pipeline_run.input_smiles,
        pipeline_config=pipeline_run.pipeline_config,
        results=pipeline_run.results,
        error_message=pipeline_run.error_message,
        user_id=pipeline_run.user_id
    )


@router.get("/runs", response_model=RunListResponse)
def list_runs(
    page: int = Query(1, ge=1, description="Page number (1-indexed)"),
    page_size: int = Query(20, ge=1, le=100, description="Items per page"),
    status: Optional[str] = Query(None, description="Filter by status"),
    user_id: Optional[str] = Query(None, description="Filter by user_id"),
    db: Session = Depends(get_db)
):
    """
    List all pipeline runs with pagination

    Args:
        page: Page number (1-indexed)
        page_size: Number of items per page (max 100)
        status: Optional status filter (pending, running, completed, failed)
        user_id: Optional user_id filter
        db: Database session

    Returns:
        RunListResponse with paginated runs
    """
    logger.info(f"Listing runs: page={page}, page_size={page_size}, status={status}, user_id={user_id}")

    # Build query
    query = db.query(PipelineRun)

    # Apply filters
    if status:
        query = query.filter(PipelineRun.status == status)
    if user_id:
        query = query.filter(PipelineRun.user_id == user_id)

    # Get total count
    total = query.count()

    # Apply pagination
    offset = (page - 1) * page_size
    runs = query.order_by(PipelineRun.created_at.desc()).offset(offset).limit(page_size).all()

    logger.info(f"Found {total} runs, returning {len(runs)} for page {page}")

    # Convert to response models
    run_responses = [
        RunResponse(
            run_id=run.run_id,
            status=run.status,
            created_at=run.created_at,
            updated_at=run.updated_at,
            completed_at=run.completed_at,
            input_smiles=run.input_smiles,
            pipeline_config=run.pipeline_config,
            results=run.results,
            error_message=run.error_message,
            user_id=run.user_id
        )
        for run in runs
    ]

    return RunListResponse(
        runs=run_responses,
        total=total,
        page=page,
        page_size=page_size
    )


@router.delete("/runs/{run_id}", status_code=204)
def delete_run(
    run_id: str,
    db: Session = Depends(get_db)
):
    """
    Delete a pipeline run and all associated results

    Args:
        run_id: Unique run identifier
        db: Database session

    Raises:
        HTTPException: 404 if run not found
    """
    logger.info(f"Deleting run: {run_id}")

    # Query database
    pipeline_run = db.query(PipelineRun).filter(PipelineRun.run_id == run_id).first()

    if not pipeline_run:
        logger.warning(f"Run not found: {run_id}")
        raise HTTPException(status_code=404, detail=f"Run not found: {run_id}")

    # Delete (cascade will handle compound_results)
    db.delete(pipeline_run)
    db.commit()

    logger.info(f"Deleted run: {run_id}")
    return None
