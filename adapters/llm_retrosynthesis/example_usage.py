"""
Example Usage of LLM Retrosynthesis Adapter

This script demonstrates various use cases of the LLM retrosynthesis adapter.
"""

import asyncio
import os
from pathlib import Path
import sys

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from adapters.llm_retrosynthesis import LLMRetrosynthesisAdapter


async def example_basic_usage():
    """Example 1: Basic usage with default settings."""
    print("\n" + "="*60)
    print("EXAMPLE 1: Basic Usage")
    print("="*60)

    # Initialize with defaults (Claude)
    adapter = LLMRetrosynthesisAdapter()

    # Run retrosynthesis for aspirin
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    result = await adapter.execute(aspirin)

    if result.success:
        print(f"\n✓ Found {result.data['routes_found']} routes")
        print(f"Best route: {result.data['n_steps']} steps")
        print(f"Score: {result.data['synthesis_score']:.3f}")
        print(f"Feasibility: {result.data['feasibility']}")

        # Show starting materials
        print("\nStarting materials needed:")
        for material in result.data['requirements']['starting_materials']:
            print(f"  - {material}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_custom_config():
    """Example 2: Custom configuration."""
    print("\n" + "="*60)
    print("EXAMPLE 2: Custom Configuration")
    print("="*60)

    # Configure adapter with specific settings
    adapter = LLMRetrosynthesisAdapter(config={
        "provider": "claude",
        "model": "claude-3-5-sonnet-20241022",
        "num_routes": 5,        # Request 5 routes
        "temperature": 0.5,     # Lower temperature for more conservative routes
        "max_tokens": 4000
    })

    # Test with ibuprofen
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    print(f"\nTarget: Ibuprofen")
    print(f"SMILES: {ibuprofen}")

    result = await adapter.execute(ibuprofen)

    if result.success:
        print(f"\n✓ Generated {result.data['routes_found']} routes")

        # Show all routes
        for route in result.data['routes']:
            print(f"\nRoute {route['route_id']}:")
            print(f"  Steps: {route['n_steps']}")
            print(f"  Score: {route['synthesis_score']:.3f} ({route['feasibility']})")
            print(f"  Strategy: {route.get('strategy', 'N/A')[:60]}...")


async def example_openai_provider():
    """Example 3: Using OpenAI instead of Claude."""
    print("\n" + "="*60)
    print("EXAMPLE 3: OpenAI Provider")
    print("="*60)

    # Initialize with OpenAI
    adapter = LLMRetrosynthesisAdapter(config={
        "provider": "openai",
        "model": "gpt-4o",
        "num_routes": 3
    })

    # Check if API key is available
    if not adapter.api_key:
        print("\n⚠ OpenAI API key not found")
        print("Set OPENAI_API_KEY environment variable to run this example")
        return

    caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    print(f"\nTarget: Caffeine")
    print(f"SMILES: {caffeine}")

    result = await adapter.execute(caffeine)

    if result.success:
        print(f"\n✓ Found {result.data['routes_found']} routes using GPT-4")
        best = result.data['best_route']
        print(f"\nBest route details:")
        print(f"  Steps: {best['n_steps']}")
        print(f"  Score: {best['synthesis_score']:.3f}")
        print(f"  Strategy: {best.get('strategy', 'N/A')}")


async def example_error_handling():
    """Example 4: Error handling."""
    print("\n" + "="*60)
    print("EXAMPLE 4: Error Handling")
    print("="*60)

    adapter = LLMRetrosynthesisAdapter()

    # Test with invalid SMILES
    print("\nTest 1: Invalid SMILES")
    result = await adapter.execute("not_a_valid_smiles")
    print(f"  Success: {result.success}")
    print(f"  Error: {result.error}")

    # Test with empty SMILES
    print("\nTest 2: Empty SMILES")
    result = await adapter.execute("")
    print(f"  Success: {result.success}")
    print(f"  Error: {result.error}")

    # Test without API key
    print("\nTest 3: No API key")
    adapter_no_key = LLMRetrosynthesisAdapter(config={
        "provider": "claude",
        "api_key": None
    })
    # Clear environment variable temporarily
    old_key = os.environ.get("ANTHROPIC_API_KEY")
    if old_key:
        del os.environ["ANTHROPIC_API_KEY"]

    result = await adapter_no_key.execute("CC(=O)O")
    print(f"  Success: {result.success}")
    print(f"  Error: {result.error}")

    # Restore key
    if old_key:
        os.environ["ANTHROPIC_API_KEY"] = old_key


async def example_detailed_route_info():
    """Example 5: Accessing detailed route information."""
    print("\n" + "="*60)
    print("EXAMPLE 5: Detailed Route Information")
    print("="*60)

    adapter = LLMRetrosynthesisAdapter(config={"num_routes": 2})

    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    result = await adapter.execute(aspirin)

    if result.success and result.data['routes_found'] > 0:
        best_route = result.data['best_route']

        print(f"\n=== BEST ROUTE FOR ASPIRIN ===")
        print(f"\nOverall Strategy:")
        print(f"  {best_route.get('strategy', 'N/A')}")

        print(f"\nSynthesis Steps:")
        for step in best_route.get('steps', []):
            print(f"\n  Step {step['step_number']}: {step['description']}")
            print(f"    Reagents: {', '.join(step.get('reagents', []))}")
            print(f"    Conditions: {step.get('conditions', 'N/A')}")
            if 'product_smiles' in step:
                print(f"    Product: {step['product_smiles']}")

        print(f"\nStarting Materials:")
        for material in best_route.get('starting_materials', []):
            print(f"  - {material}")

        print(f"\nChallenges:")
        print(f"  {best_route.get('challenges', 'None specified')}")

        print(f"\nAdvantages:")
        print(f"  {best_route.get('advantages', 'None specified')}")

        print(f"\nMetrics:")
        print(f"  Number of steps: {best_route['n_steps']}")
        print(f"  Synthesis score: {best_route['synthesis_score']:.3f}")
        print(f"  Feasibility: {best_route['feasibility']}")


async def example_comparison_multiple_molecules():
    """Example 6: Compare routes for multiple molecules."""
    print("\n" + "="*60)
    print("EXAMPLE 6: Multiple Molecules Comparison")
    print("="*60)

    adapter = LLMRetrosynthesisAdapter(config={"num_routes": 2})

    molecules = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    }

    results = []

    for name, smiles in molecules.items():
        print(f"\nAnalyzing {name}...")
        result = await adapter.execute(smiles)

        if result.success:
            results.append({
                "name": name,
                "smiles": smiles,
                "routes": result.data['routes_found'],
                "steps": result.data['n_steps'],
                "score": result.data['synthesis_score'],
                "feasibility": result.data['feasibility']
            })

    # Display comparison table
    print("\n" + "="*60)
    print("COMPARISON TABLE")
    print("="*60)
    print(f"{'Molecule':<15} {'Routes':<8} {'Steps':<8} {'Score':<8} {'Feasibility'}")
    print("-" * 60)

    for r in results:
        print(f"{r['name']:<15} {r['routes']:<8} {r['steps']:<8} "
              f"{r['score']:<8.3f} {r['feasibility']}")


async def example_with_caching():
    """Example 7: Demonstrate caching behavior."""
    print("\n" + "="*60)
    print("EXAMPLE 7: Caching Demonstration")
    print("="*60)

    adapter = LLMRetrosynthesisAdapter()

    aspirin = "CC(=O)Oc1ccccc1C(=O)O"

    # First call - will hit API
    print("\nFirst call (should call API)...")
    import time
    start = time.time()
    result1 = await adapter(aspirin, use_cache=True)
    elapsed1 = time.time() - start

    print(f"  Cache hit: {result1.cache_hit}")
    print(f"  Time: {elapsed1:.2f}s")

    # Second call - should use cache
    print("\nSecond call (should use cache)...")
    start = time.time()
    result2 = await adapter(aspirin, use_cache=True)
    elapsed2 = time.time() - start

    print(f"  Cache hit: {result2.cache_hit}")
    print(f"  Time: {elapsed2:.2f}s")
    print(f"  Speedup: {elapsed1/elapsed2:.1f}x faster")

    # Call without cache
    print("\nThird call (cache disabled)...")
    start = time.time()
    result3 = await adapter(aspirin, use_cache=False)
    elapsed3 = time.time() - start

    print(f"  Cache hit: {result3.cache_hit}")
    print(f"  Time: {elapsed3:.2f}s")


async def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("LLM RETROSYNTHESIS ADAPTER - USAGE EXAMPLES")
    print("="*70)

    # Check if API key is available
    adapter = LLMRetrosynthesisAdapter()
    if not adapter.api_key:
        print("\n⚠ WARNING: No API key found!")
        print("\nTo run these examples, set an API key:")
        print("  export ANTHROPIC_API_KEY=your-key  # for Claude")
        print("  export OPENAI_API_KEY=your-key     # for OpenAI")
        print("\nRunning error handling example only...\n")

        await example_error_handling()
        return

    print("\nAPI key found - running all examples...")

    # Run examples
    try:
        await example_basic_usage()
        await example_custom_config()
        await example_openai_provider()
        await example_error_handling()
        await example_detailed_route_info()
        await example_comparison_multiple_molecules()
        await example_with_caching()

        print("\n" + "="*70)
        print("ALL EXAMPLES COMPLETED")
        print("="*70)

    except Exception as e:
        print(f"\n❌ Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(main())
