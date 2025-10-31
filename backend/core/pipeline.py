"""
Pipeline Orchestration System

Chains adapters together to execute complete drug discovery workflows.
Handles progress tracking, error handling, and result aggregation.
"""

import logging
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, field
import asyncio
from datetime import datetime
from sqlalchemy.orm import Session

from backend.core.adapters.protocol import registry, AdapterResult
from backend.core.ranking import MultiObjectiveRanker, CompoundScore
from backend.db.models import PipelineRun, CompoundResult
from backend.core.scoring_utils import vina_affinity_to01, synthesis_steps_to01

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Pipeline execution configuration."""
    adapter_sequence: List[str] = field(default_factory=lambda: [
        "rdkit_local",
        "admet_ai",
        "pubchem",
        "chembl"
    ])
    use_caching: bool = True
    timeout_minutes: int = 60
    ranking_method: str = "pareto"  # "pareto", "weighted", or "hybrid"
    ranking_weights: Optional[Dict[str, float]] = None
    n_top_candidates: int = 10


class Pipeline:
    """
    Drug discovery pipeline orchestrator.

    Executes adapters in sequence:
    1. RDKit local properties
    2. ADMET predictions
    3. PubChem data enrichment
    4. ChEMBL bioactivity data
    5. Multi-objective ranking

    Tracks progress and handles errors gracefully.
    """

    def __init__(
        self,
        config: PipelineConfig,
        db: Session
    ):
        """
        Initialize pipeline.

        Args:
            config: Pipeline configuration
            db: Database session for progress tracking
        """
        self.config = config
        self.db = db
        logger.info(f"Initialized pipeline with adapters: {', '.join(config.adapter_sequence)}")

    async def execute(
        self,
        input_smiles: List[str],
        run_id: str
    ) -> Dict[str, Any]:
        """
        Execute complete pipeline.

        Args:
            input_smiles: List of input SMILES strings
            run_id: Database run ID for tracking

        Returns:
            Dictionary containing ranked results and metadata
        """
        logger.info(f"Starting pipeline execution for run {run_id}")
        logger.info(f"Processing {len(input_smiles)} compounds")

        # Retrieve pipeline run from database
        run = self.db.query(PipelineRun).filter(PipelineRun.run_id == run_id).first()

        if not run:
            error_msg = f"Pipeline run not found: {run_id}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        # Update status to RUNNING
        run.status = "running"
        run.updated_at = datetime.utcnow()
        self.db.commit()
        logger.info(f"Pipeline run {run_id} status updated to RUNNING")

        try:
            # Verify all required adapters are registered
            self._verify_adapters()

            # Process all compounds through adapter sequence
            all_results = await self._process_compounds(input_smiles, run)

            # Rank compounds using multi-objective ranking
            ranked_results = self._rank_compounds(all_results)

            # Store results in database
            self._store_results(run, all_results, ranked_results)

            # Update status to COMPLETED
            run.status = "completed"
            run.completed_at = datetime.utcnow()
            run.updated_at = datetime.utcnow()
            self.db.commit()

            logger.info(f"Pipeline run {run_id} COMPLETED successfully")
            logger.info(f"Processed {len(all_results)} compounds")

            return {
                "success": True,
                "run_id": run_id,
                "n_input": len(input_smiles),
                "n_processed": len(all_results),
                "n_top_candidates": min(self.config.n_top_candidates, len(ranked_results)),
                "top_compounds": [c.to_dict() for c in ranked_results[:self.config.n_top_candidates]],
                "status": "completed",
                "completed_at": run.completed_at.isoformat() if run.completed_at else None
            }

        except Exception as e:
            error_msg = f"Pipeline execution failed: {str(e)}"
            logger.error(error_msg, exc_info=True)

            # Update run status to FAILED
            run.status = "failed"
            run.error_message = error_msg
            run.updated_at = datetime.utcnow()
            self.db.commit()

            raise

    def _verify_adapters(self):
        """
        Verify all required adapters are registered.

        Raises:
            ValueError: If any required adapter is not registered
        """
        missing_adapters = []

        for adapter_name in self.config.adapter_sequence:
            if registry.get(adapter_name) is None:
                missing_adapters.append(adapter_name)

        if missing_adapters:
            error_msg = f"Required adapters not registered: {', '.join(missing_adapters)}"
            logger.error(error_msg)
            raise ValueError(error_msg)

        logger.info(f"✓ All {len(self.config.adapter_sequence)} required adapters are registered")

    async def _process_compounds(
        self,
        input_smiles: List[str],
        run: PipelineRun
    ) -> List[Dict[str, Any]]:
        """
        Process all compounds through adapter sequence.

        Args:
            input_smiles: List of SMILES strings
            run: Pipeline run database object

        Returns:
            List of compound data dictionaries
        """
        all_results = []
        total_compounds = len(input_smiles)
        total_steps = total_compounds * len(self.config.adapter_sequence)
        current_step = 0

        for idx, smiles in enumerate(input_smiles):
            logger.info(f"Processing compound {idx + 1}/{total_compounds}: {smiles}")

            compound_data = {
                "smiles": smiles,
                "properties": None,
                "admet": None,
                "pubchem": None,
                "bioactivity": None,
                "errors": {}
            }

            # Execute adapters sequentially
            for adapter_name in self.config.adapter_sequence:
                adapter = registry.get(adapter_name)

                try:
                    logger.debug(f"Running adapter: {adapter_name} for {smiles}")

                    # Execute adapter (call directly, caching is handled in protocol)
                    result: AdapterResult = await adapter(
                        smiles,
                        use_cache=self.config.use_caching
                    )

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

                        logger.debug(f"✓ Adapter {adapter_name} completed successfully")
                    else:
                        logger.warning(f"✗ Adapter {adapter_name} failed: {result.error}")
                        compound_data["errors"][adapter_name] = result.error

                except Exception as e:
                    error_msg = f"Exception in adapter {adapter_name}: {str(e)}"
                    logger.error(error_msg, exc_info=True)
                    compound_data["errors"][adapter_name] = error_msg

                # Update progress
                current_step += 1
                progress = current_step / total_steps
                run.updated_at = datetime.utcnow()
                self.db.commit()

                logger.debug(f"Progress: {progress:.1%} ({current_step}/{total_steps})")

            all_results.append(compound_data)
            logger.info(f"Compound {idx + 1}/{total_compounds} completed")

        return all_results

    def _rank_compounds(
        self,
        compound_data_list: List[Dict[str, Any]]
    ) -> List[CompoundScore]:
        """
        Rank compounds using multi-objective ranking.

        Args:
            compound_data_list: List of compound data dictionaries

        Returns:
            List of ranked CompoundScore objects
        """
        logger.info("Starting multi-objective ranking")

        # Prepare data for ranking
        ranking_data = []

        for compound_data in compound_data_list:
            # Extract and normalize scores
            smiles = compound_data["smiles"]

            # Binding score (from docking if available, else 0)
            binding_score = 0.0
            if compound_data.get("bioactivity"):
                # TODO: Extract actual binding affinity when docking adapter is available
                # For now, use placeholder
                binding_score = 0.5

            # ADMET score (composite from ADMET predictions)
            admet_score = 0.0
            if compound_data.get("admet"):
                admet_data = compound_data["admet"]
                # Average ADMET predictions if available
                if isinstance(admet_data, dict):
                    scores = [v for v in admet_data.values() if isinstance(v, (int, float))]
                    if scores:
                        admet_score = sum(scores) / len(scores)

            # Synthesis score (placeholder for now)
            synthesis_score = 0.5

            # Novelty score (placeholder for now)
            novelty_score = 0.5

            ranking_data.append({
                "smiles": smiles,
                "binding_score": binding_score,
                "admet_score": admet_score,
                "synthesis_score": synthesis_score,
                "novelty_score": novelty_score
            })

        # Create ranker
        ranker = MultiObjectiveRanker(weights=self.config.ranking_weights)

        # Rank compounds
        ranked_compounds = ranker.rank(
            ranking_data,
            method=self.config.ranking_method
        )

        logger.info(f"Ranked {len(ranked_compounds)} compounds")
        logger.info(f"Top compound: {ranked_compounds[0].smiles} (score: {ranked_compounds[0].composite_score:.3f})")

        return ranked_compounds

    def _store_results(
        self,
        run: PipelineRun,
        all_results: List[Dict[str, Any]],
        ranked_compounds: List[CompoundScore]
    ):
        """
        Store results in database.

        Args:
            run: Pipeline run database object
            all_results: List of compound data dictionaries
            ranked_compounds: List of ranked CompoundScore objects
        """
        logger.info("Storing results in database")

        # Create a mapping from SMILES to rank and scores
        rank_map = {
            compound.smiles: {
                "rank": idx + 1,
                "composite_score": compound.composite_score,
                "pareto_rank": compound.pareto_rank,
                "binding_score": compound.binding_score,
                "admet_score": compound.admet_score,
                "synthesis_score": compound.synthesis_score,
                "novelty_score": compound.novelty_score
            }
            for idx, compound in enumerate(ranked_compounds)
        }

        # Store individual compound results
        for idx, compound_data in enumerate(all_results):
            smiles = compound_data["smiles"]
            rank_info = rank_map.get(smiles, {})

            compound_result = CompoundResult(
                compound_id=f"{run.run_id}_{idx}",
                run_id=run.id,
                smiles=smiles,
                properties=compound_data.get("properties"),
                admet=compound_data.get("admet"),
                bioactivity=compound_data.get("bioactivity"),
                scores=rank_info,
                final_score=rank_info.get("composite_score"),
                rank=rank_info.get("rank")
            )

            self.db.add(compound_result)

        # Store aggregated results in run
        run.results = {
            "compounds": [c.to_dict() for c in ranked_compounds],
            "total_processed": len(all_results),
            "adapter_sequence": self.config.adapter_sequence,
            "ranking_method": self.config.ranking_method,
            "completed_at": datetime.utcnow().isoformat()
        }

        self.db.commit()
        logger.info(f"Stored {len(all_results)} compound results in database")


def create_pipeline(
    db: Session,
    config: Optional[PipelineConfig] = None
) -> Pipeline:
    """
    Factory function to create a pipeline instance.

    Args:
        db: Database session
        config: Optional pipeline configuration (uses defaults if not provided)

    Returns:
        Configured Pipeline instance
    """
    if config is None:
        config = PipelineConfig()

    return Pipeline(config=config, db=db)
