"""
Multi-objective ranking system.

Combines multiple scores (binding, ADMET, synthesis, novelty)
into a final ranking using Pareto optimization.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

@dataclass
class CompoundScore:
    """Compound with multi-objective scores."""
    smiles: str
    binding_score: float
    admet_score: float
    synthesis_score: float
    novelty_score: float
    composite_score: float = 0.0
    pareto_rank: int = 0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        return {
            "smiles": self.smiles,
            "binding_score": self.binding_score,
            "admet_score": self.admet_score,
            "synthesis_score": self.synthesis_score,
            "novelty_score": self.novelty_score,
            "composite_score": self.composite_score,
            "pareto_rank": self.pareto_rank
        }

class ParetoRanker:
    """
    Pareto frontier ranking.

    Ranks compounds by non-dominated sorting:
    - Rank 1 = Pareto frontier (not dominated by any other)
    - Rank 2 = Dominated only by Rank 1
    - etc.
    """

    @staticmethod
    def dominates(a: CompoundScore, b: CompoundScore) -> bool:
        """
        Check if compound a dominates compound b.

        a dominates b if:
        - a is >= b on all objectives
        - a is > b on at least one objective

        Args:
            a: First compound
            b: Second compound

        Returns:
            True if a dominates b, False otherwise
        """
        better_or_equal = (
            a.binding_score >= b.binding_score and
            a.admet_score >= b.admet_score and
            a.synthesis_score >= b.synthesis_score and
            a.novelty_score >= b.novelty_score
        )

        strictly_better = (
            a.binding_score > b.binding_score or
            a.admet_score > b.admet_score or
            a.synthesis_score > b.synthesis_score or
            a.novelty_score > b.novelty_score
        )

        return better_or_equal and strictly_better

    @staticmethod
    def compute_pareto_ranks(compounds: List[CompoundScore]) -> List[CompoundScore]:
        """
        Compute Pareto ranks for all compounds.

        Uses non-dominated sorting to assign ranks. Compounds on the
        Pareto frontier (not dominated by any other) get rank 1.
        After removing rank 1, the next frontier gets rank 2, etc.

        Args:
            compounds: List of compounds to rank

        Returns:
            Same list of compounds with pareto_rank assigned
        """
        remaining = compounds.copy()
        current_rank = 1

        while remaining:
            # Find Pareto frontier in remaining compounds
            frontier = []

            for comp in remaining:
                # Check if dominated by any other remaining compound
                is_dominated = any(
                    ParetoRanker.dominates(other, comp)
                    for other in remaining
                    if other != comp
                )

                if not is_dominated:
                    frontier.append(comp)

            # Assign rank to frontier
            for comp in frontier:
                comp.pareto_rank = current_rank

            # Remove frontier from remaining
            remaining = [c for c in remaining if c not in frontier]
            current_rank += 1

        return compounds

class MultiObjectiveRanker:
    """
    Multi-objective ranking system.

    Combines multiple adapter outputs into final ranking using
    different strategies (Pareto, weighted, or hybrid).
    """

    def __init__(
        self,
        weights: Optional[Dict[str, float]] = None
    ):
        """
        Initialize ranker.

        Args:
            weights: Objective weights for composite score
                     (binding, admet, synthesis, novelty)
                     Default: equal weights (0.25 each)
        """
        self.weights = weights or {
            "binding": 0.25,
            "admet": 0.25,
            "synthesis": 0.25,
            "novelty": 0.25
        }

        # Validate weights sum to 1.0
        weight_sum = sum(self.weights.values())
        if not (0.99 <= weight_sum <= 1.01):  # Allow small floating point error
            logger.warning(f"Weights sum to {weight_sum}, normalizing to 1.0")
            total = sum(self.weights.values())
            self.weights = {k: v/total for k, v in self.weights.items()}

    def rank(
        self,
        compounds_data: List[Dict[str, Any]],
        method: str = "pareto"
    ) -> List[CompoundScore]:
        """
        Rank compounds.

        Args:
            compounds_data: List of dicts with scores for each compound.
                           Each dict should have keys: smiles, binding_score,
                           admet_score, synthesis_score, novelty_score
            method: Ranking method ("pareto", "weighted", or "hybrid")
                   - "pareto": Pure Pareto ranking, ties broken by composite score
                   - "weighted": Sort by weighted composite score only
                   - "hybrid": Same as "pareto" (Pareto + composite tie-breaking)

        Returns:
            Sorted list of CompoundScore objects
        """
        logger.info(f"Ranking {len(compounds_data)} compounds using {method}")

        # Convert to CompoundScore objects
        compounds = []
        for data in compounds_data:
            # Normalize scores to 0-1 (defensive in case already normalized)
            compound = CompoundScore(
                smiles=data["smiles"],
                binding_score=self._normalize_score(data.get("binding_score", 0)),
                admet_score=self._normalize_score(data.get("admet_score", 0)),
                synthesis_score=self._normalize_score(data.get("synthesis_score", 0)),
                novelty_score=self._normalize_score(data.get("novelty_score", 0))
            )

            # Compute composite score (weighted average)
            compound.composite_score = (
                self.weights["binding"] * compound.binding_score +
                self.weights["admet"] * compound.admet_score +
                self.weights["synthesis"] * compound.synthesis_score +
                self.weights["novelty"] * compound.novelty_score
            )

            compounds.append(compound)

        # Rank based on method
        if method == "pareto":
            compounds = ParetoRanker.compute_pareto_ranks(compounds)
            # Sort by Pareto rank, then composite score
            compounds.sort(key=lambda c: (c.pareto_rank, -c.composite_score))

        elif method == "weighted":
            # Sort by composite score only
            compounds.sort(key=lambda c: -c.composite_score)

        elif method == "hybrid":
            # Pareto ranking, but break ties with composite score
            compounds = ParetoRanker.compute_pareto_ranks(compounds)
            compounds.sort(key=lambda c: (c.pareto_rank, -c.composite_score))

        else:
            raise ValueError(f"Unknown ranking method: {method}. Use 'pareto', 'weighted', or 'hybrid'")

        return compounds

    @staticmethod
    def _normalize_score(score: float) -> float:
        """
        Normalize score to 0-1 range.

        Args:
            score: Raw score value

        Returns:
            Score clamped to [0, 1]
        """
        return max(0.0, min(1.0, score))

def rank_pipeline_results(
    results: Dict[str, Any],
    method: str = "pareto",
    weights: Optional[Dict[str, float]] = None
) -> pd.DataFrame:
    """
    Convenience function to rank pipeline results.

    Args:
        results: Dict with adapter results for each compound.
                Keys are SMILES strings, values are dicts with adapter outputs.
        method: Ranking method ("pareto", "weighted", or "hybrid")
        weights: Optional custom weights for composite scoring

    Returns:
        DataFrame with ranked compounds, 1-indexed ranking
    """
    ranker = MultiObjectiveRanker(weights=weights)

    # Extract compound data
    compounds_data = []
    for smiles, data in results.items():
        compounds_data.append({
            "smiles": smiles,
            "binding_score": data.get("docking", {}).get("binding_affinity", 0),
            "admet_score": data.get("admet", {}).get("composite_score", 0),
            "synthesis_score": data.get("retrosynthesis", {}).get("synthesis_score", 0),
            "novelty_score": data.get("novelty", {}).get("novelty_score", 0)
        })

    # Rank
    ranked_compounds = ranker.rank(compounds_data, method=method)

    # Convert to DataFrame
    df = pd.DataFrame([c.to_dict() for c in ranked_compounds])
    df.index = range(1, len(df) + 1)  # 1-indexed ranking

    return df
