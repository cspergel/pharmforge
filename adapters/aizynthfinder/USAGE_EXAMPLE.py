"""
AiZynthFinder Adapter Usage Examples

Complete examples showing how to use the AiZynthFinder adapter
in various PharmForge scenarios.
"""

import asyncio
import sys
import os

# Add PharmForge to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from backend.core.adapter_registry import registry, register_all_adapters
from adapters.aizynthfinder.adapter import AiZynthFinderAdapter


async def example_1_basic_usage():
    """
    Example 1: Basic retrosynthesis analysis
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 1: Basic Retrosynthesis Analysis")
    print("=" * 70)

    # Create adapter
    adapter = AiZynthFinderAdapter()

    # Aspirin synthesis
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"

    print(f"\nAnalyzing: Aspirin")
    print(f"SMILES: {aspirin_smiles}")

    # Note: This will only work if AiZynthFinder is installed
    # Otherwise, it will return an error with installation instructions
    result = await adapter.execute(aspirin_smiles, max_routes=3, expansion_time=30)

    if result.success:
        if result.data["routes_found"] > 0:
            print(f"\n[OK] Found {result.data['routes_found']} synthesis routes")

            best = result.data["best_route"]
            print(f"\nBest Route:")
            print(f"  Steps: {best['n_steps']}")
            print(f"  Synthesis Score: {best['synthesis_score']:.3f} (0-1, higher=better)")
            print(f"  Feasibility: {best['feasibility']}")
            print(f"  Starting Materials: {len(best['starting_materials'])}")
        else:
            print("\n[WARN] No synthesis routes found")
            print("This may indicate the molecule is too complex or disconnected")
    else:
        print(f"\n[ERROR] {result.error}")
        if "not installed" in result.error:
            print("\nTo use this adapter, install AiZynthFinder:")
            print("  pip install aizynthfinder")


async def example_2_using_registry():
    """
    Example 2: Using adapter through registry (recommended)
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Using Adapter Through Registry")
    print("=" * 70)

    # Register all adapters
    register_all_adapters()

    # Get adapter from registry
    adapter = registry.get("aizynthfinder")

    if adapter is None:
        print("[ERROR] AiZynthFinder adapter not registered")
        return

    print(f"[OK] Adapter loaded: {adapter.name}")
    print(f"  Version: {adapter.version}")
    print(f"  Type: {adapter.adapter_type}")

    # Use it
    smiles = "c1ccccc1"  # Benzene
    print(f"\nAnalyzing: Benzene (should be simple/commercial)")
    print(f"SMILES: {smiles}")

    result = await adapter.execute(smiles, max_routes=1, expansion_time=20)

    if result.success:
        print(f"\n[OK] Analysis complete")
        print(f"  Routes found: {result.data['routes_found']}")
        if result.data['routes_found'] > 0:
            print(f"  Synthesis score: {result.data['synthesis_score']:.3f}")
    else:
        print(f"\n[WARN] {result.error}")


async def example_3_multiple_molecules():
    """
    Example 3: Analyzing multiple molecules
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Multi-Molecule Analysis")
    print("=" * 70)

    adapter = AiZynthFinderAdapter(config={
        "max_routes": 3,
        "expansion_time": 30
    })

    # Test molecules
    molecules = {
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    }

    results = {}

    print("\nAnalyzing molecules...")
    for name, smiles in molecules.items():
        print(f"\n  {name}...")
        result = await adapter.execute(smiles)

        if result.success:
            results[name] = result.data
            print(f"    [OK] {result.data['routes_found']} routes found")
        else:
            print(f"    [ERROR] {result.error}")

    # Rank by synthesis ease
    if results:
        print("\n" + "-" * 70)
        print("Ranking by Synthesis Ease:")
        print("-" * 70)

        ranked = sorted(
            results.items(),
            key=lambda x: x[1].get("synthesis_score", 0.0),
            reverse=True
        )

        for rank, (name, data) in enumerate(ranked, 1):
            if data['routes_found'] > 0:
                print(f"{rank}. {name}")
                print(f"   Score: {data['synthesis_score']:.3f}")
                print(f"   Steps: {data['n_steps']}")
                print(f"   Feasibility: {data['feasibility']}")
            else:
                print(f"{rank}. {name} - No routes found")


async def example_4_pipeline_integration():
    """
    Example 4: Integration with PharmForge pipeline
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 4: Pipeline Integration")
    print("=" * 70)

    # Simulating a pipeline that ranks drug candidates
    print("\nSimulating drug candidate ranking pipeline...")

    candidates = [
        ("Candidate A", "CN1CCN(CC1)C(=O)C2=C(C=CC(=C2)OC)OC"),
        ("Candidate B", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Candidate C", "c1ccccc1"),
    ]

    adapter = AiZynthFinderAdapter(config={"max_routes": 3, "expansion_time": 30})

    ranked_candidates = []

    for name, smiles in candidates:
        print(f"\n  Evaluating {name}...")

        # Run retrosynthesis
        retro_result = await adapter.execute(smiles)

        # Simulated other scores (in reality, these come from other adapters)
        simulated_docking_score = 0.8  # From Vina adapter
        simulated_admet_score = 0.7    # From ADMET-AI adapter

        if retro_result.success and retro_result.data['routes_found'] > 0:
            synthesis_score = retro_result.data['synthesis_score']
        else:
            synthesis_score = 0.0

        # Multi-objective score (equal weighting)
        combined_score = (
            simulated_docking_score * 0.4 +
            simulated_admet_score * 0.3 +
            synthesis_score * 0.3
        )

        ranked_candidates.append({
            "name": name,
            "smiles": smiles,
            "docking": simulated_docking_score,
            "admet": simulated_admet_score,
            "synthesis": synthesis_score,
            "combined": combined_score,
            "n_steps": retro_result.data.get('n_steps', None) if retro_result.success else None
        })

        print(f"    Synthesis: {synthesis_score:.3f}")
        print(f"    Combined:  {combined_score:.3f}")

    # Sort by combined score
    ranked_candidates.sort(key=lambda x: x['combined'], reverse=True)

    print("\n" + "-" * 70)
    print("Final Ranking:")
    print("-" * 70)
    print(f"{'Rank':<6}{'Name':<15}{'Docking':<10}{'ADMET':<10}{'Synthesis':<12}{'Combined':<10}{'Steps'}")
    print("-" * 70)

    for rank, cand in enumerate(ranked_candidates, 1):
        steps = cand['n_steps'] if cand['n_steps'] is not None else "N/A"
        print(
            f"{rank:<6}"
            f"{cand['name']:<15}"
            f"{cand['docking']:<10.3f}"
            f"{cand['admet']:<10.3f}"
            f"{cand['synthesis']:<12.3f}"
            f"{cand['combined']:<10.3f}"
            f"{steps}"
        )


async def example_5_custom_configuration():
    """
    Example 5: Using custom configuration
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Custom Configuration")
    print("=" * 70)

    # Create adapter with custom settings
    custom_adapter = AiZynthFinderAdapter(config={
        "max_routes": 10,          # Find more routes
        "expansion_time": 120,     # Longer search time
        "timeout": 180,            # Higher timeout
        "stock": "emolecules",     # Different stock database
    })

    print("\nCustom Configuration:")
    print(f"  Max routes: {custom_adapter.config['max_routes']}")
    print(f"  Expansion time: {custom_adapter.config['expansion_time']}s")
    print(f"  Stock: {custom_adapter.config['stock']}")

    # Get metadata
    metadata = custom_adapter.get_metadata()
    print(f"\nAdapter Metadata:")
    print(f"  Name: {metadata['name']}")
    print(f"  Description: {metadata['description']}")
    print(f"  Capabilities: {list(metadata['capabilities'].keys())}")


async def main():
    """Run all examples"""
    print("\n" + "=" * 70)
    print("AiZynthFinder Adapter - Usage Examples")
    print("=" * 70)

    print("\nNote: These examples demonstrate the adapter interface.")
    print("Full execution requires AiZynthFinder installation and configuration.")
    print("See adapters/aizynthfinder/README.md for setup instructions.")

    # Run examples
    await example_1_basic_usage()
    await example_2_using_registry()
    await example_3_multiple_molecules()
    await example_4_pipeline_integration()
    await example_5_custom_configuration()

    print("\n" + "=" * 70)
    print("Examples Complete!")
    print("=" * 70)
    print("\nNext Steps:")
    print("1. Install AiZynthFinder: pip install aizynthfinder")
    print("2. Configure models (see README.md)")
    print("3. Run your own retrosynthesis analyses!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
