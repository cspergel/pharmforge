"""
Demo script for the multi-objective ranking system.

Shows how to use the ranking system with example compound data.
"""

from backend.core.ranking import (
    MultiObjectiveRanker,
    rank_pipeline_results
)

def demo_weighted_ranking():
    """Demonstrate weighted ranking."""
    print("=" * 80)
    print("WEIGHTED RANKING DEMO")
    print("=" * 80)

    # Create ranker with custom weights (prioritize binding and ADMET)
    ranker = MultiObjectiveRanker(weights={
        "binding": 0.4,
        "admet": 0.4,
        "synthesis": 0.1,
        "novelty": 0.1
    })

    # Example compounds
    compounds = [
        {
            "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen-like
            "binding_score": 0.85,
            "admet_score": 0.75,
            "synthesis_score": 0.90,
            "novelty_score": 0.20
        },
        {
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin-like
            "binding_score": 0.70,
            "admet_score": 0.80,
            "synthesis_score": 0.95,
            "novelty_score": 0.15
        },
        {
            "smiles": "Cc1ccc(cc1)S(=O)(=O)N",  # Novel candidate
            "binding_score": 0.90,
            "admet_score": 0.85,
            "synthesis_score": 0.60,
            "novelty_score": 0.95
        },
        {
            "smiles": "c1ccccc1",  # Benzene (poor candidate)
            "binding_score": 0.30,
            "admet_score": 0.40,
            "synthesis_score": 1.00,
            "novelty_score": 0.05
        }
    ]

    ranked = ranker.rank(compounds, method="weighted")

    print("\nRanked compounds (by weighted composite score):")
    print("-" * 80)
    for i, comp in enumerate(ranked, 1):
        print(f"\nRank {i}: {comp.smiles}")
        print(f"  Composite Score: {comp.composite_score:.3f}")
        print(f"  - Binding:   {comp.binding_score:.2f}")
        print(f"  - ADMET:     {comp.admet_score:.2f}")
        print(f"  - Synthesis: {comp.synthesis_score:.2f}")
        print(f"  - Novelty:   {comp.novelty_score:.2f}")

def demo_pareto_ranking():
    """Demonstrate Pareto ranking."""
    print("\n\n" + "=" * 80)
    print("PARETO RANKING DEMO")
    print("=" * 80)

    ranker = MultiObjectiveRanker()  # Equal weights

    # Example compounds with trade-offs
    compounds = [
        {
            "smiles": "Compound_A",
            "binding_score": 0.95,  # Excellent binding
            "admet_score": 0.60,    # Poor ADMET
            "synthesis_score": 0.50,
            "novelty_score": 0.40
        },
        {
            "smiles": "Compound_B",
            "binding_score": 0.70,
            "admet_score": 0.90,    # Excellent ADMET
            "synthesis_score": 0.70,
            "novelty_score": 0.60
        },
        {
            "smiles": "Compound_C",
            "binding_score": 0.85,  # Balanced
            "admet_score": 0.80,
            "synthesis_score": 0.75,
            "novelty_score": 0.70
        },
        {
            "smiles": "Compound_D",
            "binding_score": 0.50,  # Poor all around
            "admet_score": 0.55,
            "synthesis_score": 0.60,
            "novelty_score": 0.45
        }
    ]

    ranked = ranker.rank(compounds, method="pareto")

    print("\nRanked compounds (by Pareto frontier):")
    print("-" * 80)

    current_rank = None
    for i, comp in enumerate(ranked, 1):
        if comp.pareto_rank != current_rank:
            current_rank = comp.pareto_rank
            print(f"\n--- Pareto Rank {current_rank} (Frontier {current_rank}) ---")

        print(f"\n  {comp.smiles}")
        print(f"  Composite: {comp.composite_score:.3f}")
        print(f"  Scores: B={comp.binding_score:.2f}, A={comp.admet_score:.2f}, "
              f"S={comp.synthesis_score:.2f}, N={comp.novelty_score:.2f}")

def demo_pipeline_results():
    """Demonstrate ranking with pipeline results format."""
    print("\n\n" + "=" * 80)
    print("PIPELINE RESULTS RANKING DEMO")
    print("=" * 80)

    # Simulated pipeline results
    results = {
        "CCO": {
            "docking": {"binding_affinity": 0.75},
            "admet": {"composite_score": 0.80},
            "retrosynthesis": {"synthesis_score": 0.90},
            "novelty": {"novelty_score": 0.30}
        },
        "CC(C)O": {
            "docking": {"binding_affinity": 0.85},
            "admet": {"composite_score": 0.70},
            "retrosynthesis": {"synthesis_score": 0.85},
            "novelty": {"novelty_score": 0.50}
        },
        "c1ccccc1O": {
            "docking": {"binding_affinity": 0.90},
            "admet": {"composite_score": 0.75},
            "retrosynthesis": {"synthesis_score": 0.70},
            "novelty": {"novelty_score": 0.85}
        }
    }

    # Rank with custom weights
    df = rank_pipeline_results(
        results,
        method="pareto",
        weights={
            "binding": 0.35,
            "admet": 0.35,
            "synthesis": 0.15,
            "novelty": 0.15
        }
    )

    print("\nRanking DataFrame:")
    print("-" * 80)
    print(df.to_string())

    print("\n\nTop compound:")
    print("-" * 80)
    top = df.iloc[0]
    print(f"SMILES: {top['smiles']}")
    print(f"Pareto Rank: {top['pareto_rank']}")
    print(f"Composite Score: {top['composite_score']:.3f}")
    print(f"Individual Scores:")
    print(f"  - Binding:   {top['binding_score']:.2f}")
    print(f"  - ADMET:     {top['admet_score']:.2f}")
    print(f"  - Synthesis: {top['synthesis_score']:.2f}")
    print(f"  - Novelty:   {top['novelty_score']:.2f}")

if __name__ == "__main__":
    demo_weighted_ranking()
    demo_pareto_ranking()
    demo_pipeline_results()

    print("\n" + "=" * 80)
    print("DEMO COMPLETE")
    print("=" * 80)
