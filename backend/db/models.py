"""
SQLAlchemy database models for PharmForge
"""
from sqlalchemy import Column, Integer, String, Float, DateTime, JSON, ForeignKey, Text, Boolean
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func
from .database import Base

class PipelineRun(Base):
    """
    Stores information about pipeline execution runs
    """
    __tablename__ = "pipeline_runs"

    id = Column(Integer, primary_key=True, index=True)
    run_id = Column(String(64), unique=True, index=True, nullable=False)
    status = Column(String(32), default="pending")  # pending, running, completed, failed
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
    completed_at = Column(DateTime(timezone=True), nullable=True)

    # Input data
    input_smiles = Column(JSON, nullable=False)  # List of SMILES strings
    pipeline_config = Column(JSON, nullable=True)  # Pipeline configuration

    # Results
    results = Column(JSON, nullable=True)  # Final ranked results
    lockfile = Column(JSON, nullable=True)  # Reproducibility lockfile

    # Metadata
    user_id = Column(String(128), nullable=True)
    error_message = Column(Text, nullable=True)

    # Relationships
    compound_results = relationship("CompoundResult", back_populates="pipeline_run", cascade="all, delete-orphan")


class CompoundResult(Base):
    """
    Stores individual compound results within a pipeline run
    """
    __tablename__ = "compound_results"

    id = Column(Integer, primary_key=True, index=True)
    compound_id = Column(String(128), index=True, nullable=False)
    run_id = Column(Integer, ForeignKey("pipeline_runs.id"), nullable=False)

    # Compound data
    smiles = Column(String(512), nullable=False)
    name = Column(String(256), nullable=True)

    # Adapter results
    properties = Column(JSON, nullable=True)  # Molecular properties (from PubChem/RDKit)
    bioactivity = Column(JSON, nullable=True)  # Bioactivity data (from ChEMBL)
    admet = Column(JSON, nullable=True)  # ADMET predictions (from TDC)
    docking = Column(JSON, nullable=True)  # Docking scores (from Vina/DiffDock)
    synthesis = Column(JSON, nullable=True)  # Retrosynthesis routes (from AiZynthFinder)

    # Scoring
    scores = Column(JSON, nullable=True)  # Normalized scores (0-1, higher=better)
    final_score = Column(Float, nullable=True)  # Weighted or Pareto rank
    rank = Column(Integer, nullable=True)  # Final ranking position

    # Timestamps
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())

    # Relationships
    pipeline_run = relationship("PipelineRun", back_populates="compound_results")


class CacheEntry(Base):
    """
    Stores cached adapter results for reproducibility
    """
    __tablename__ = "cache_entries"

    id = Column(Integer, primary_key=True, index=True)
    cache_key = Column(String(256), unique=True, index=True, nullable=False)
    adapter_name = Column(String(64), index=True, nullable=False)

    # Input and output
    input_data = Column(JSON, nullable=False)
    output_data = Column(JSON, nullable=False)

    # Metadata
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    accessed_at = Column(DateTime(timezone=True), server_default=func.now(), onupdate=func.now())
    access_count = Column(Integer, default=1)

    # Flags
    is_valid = Column(Boolean, default=True)
    version = Column(String(32), nullable=True)  # Adapter version


class Adapter(Base):
    """
    Registry of available adapters
    """
    __tablename__ = "adapters"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String(64), unique=True, index=True, nullable=False)
    adapter_type = Column(String(32), nullable=False)  # api, local, ml

    # Configuration
    config = Column(JSON, nullable=True)
    enabled = Column(Boolean, default=True)

    # Metadata
    description = Column(Text, nullable=True)
    version = Column(String(32), nullable=True)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
    updated_at = Column(DateTime(timezone=True), onupdate=func.now())
