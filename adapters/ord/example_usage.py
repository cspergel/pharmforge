"""
Example Usage of the Open Reaction Database (ORD) Adapter

This script demonstrates various use cases for the ORD adapter in PharmForge.
"""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.ord.adapter import ORDAdapter


async def example_1_basic_product_search():
    """
    Example 1: Basic search for reactions producing a target molecule
    """
    print("\n" + "="*70)
    print("Example 1: Search for reactions producing aspirin")
    print("="*70)

    adapter = ORDAdapter()

    # Aspirin SMILES
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"

    result = await adapter.execute(aspirin)

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['reactions_found']} reactions")
        print(f"  Average yield: {data['statistics']['average_yield']:.1f}%")
        print(f"  Max yield: {data['statistics']['max_yield']:.1f}%")
        print(f"  Reaction types: {', '.join(data['reaction_types'][:3])}")

        # Show top 3 reactions
        print("\nTop 3 Reactions by Yield:")
        for i, reaction in enumerate(data['reactions'][:3], 1):
            print(f"\n  {i}. Reaction ID: {reaction['reaction_id']}")
            print(f"     Yield: {reaction['yield']}%")
            print(f"     Reactants: {reaction['n_reactants']}")
            print(f"     Temperature: {reaction['conditions']['temperature']} °C")
            print(f"     Solvents: {', '.join(reaction['conditions']['solvents'][:2])}")
            if reaction['literature']['doi']:
                print(f"     DOI: {reaction['literature']['doi']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_2_search_by_reactant():
    """
    Example 2: Find reactions using benzene as a reactant
    """
    print("\n" + "="*70)
    print("Example 2: Find reactions using benzene as reactant")
    print("="*70)

    adapter = ORDAdapter()

    search_params = {
        "smiles": "c1ccccc1",  # Benzene
        "search_type": "reactant",
        "min_yield": 50.0  # Only high-yield reactions
    }

    result = await adapter.execute(search_params)

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['reactions_found']} reactions using benzene")
        print(f"  Average yield: {data['statistics']['average_yield']:.1f}%")

        # Analyze products
        products = set()
        for reaction in data['reactions']:
            for product in reaction['products']:
                if product['is_desired']:
                    products.add(product['name'] or product['smiles'][:30])

        print(f"\n  Unique products formed: {len(products)}")
        print(f"  Example products: {list(products)[:5]}")
    else:
        print(f"✗ Error: {result.error}")


async def example_3_catalyst_search():
    """
    Example 3: Find reactions using palladium catalysts
    """
    print("\n" + "="*70)
    print("Example 3: Find coupling reactions with Pd catalyst")
    print("="*70)

    adapter = ORDAdapter(config={
        "max_results": 100,
        "min_yield": 60.0
    })

    # Note: This is a simplified example
    # In practice, you'd search for products and filter by catalyst
    search_params = {
        "smiles": "c1ccc(-c2ccccc2)cc1",  # Biphenyl (common coupling product)
        "search_type": "product"
    }

    result = await adapter.execute(search_params)

    if result.success:
        data = result.data

        # Filter reactions using Pd catalyst
        pd_reactions = []
        for reaction in data['reactions']:
            reagents = reaction['conditions']['reagents']
            if any('Pd' in r for r in reagents):
                pd_reactions.append(reaction)

        print(f"\n✓ Found {len(pd_reactions)} Pd-catalyzed coupling reactions")

        # Analyze common Pd catalysts
        pd_catalysts = {}
        for reaction in pd_reactions:
            for reagent in reaction['conditions']['reagents']:
                if 'Pd' in reagent:
                    pd_catalysts[reagent] = pd_catalysts.get(reagent, 0) + 1

        print("\n  Most common Pd catalysts:")
        for catalyst, count in sorted(pd_catalysts.items(), key=lambda x: x[1], reverse=True)[:5]:
            print(f"    - {catalyst}: {count} reactions")
    else:
        print(f"✗ Error: {result.error}")


async def example_4_condition_optimization():
    """
    Example 4: Analyze optimal conditions for amide bond formation
    """
    print("\n" + "="*70)
    print("Example 4: Analyze conditions for amide coupling")
    print("="*70)

    adapter = ORDAdapter()

    # Simple amide product
    amide_smiles = "CC(=O)NCC"

    result = await adapter.execute(amide_smiles, min_yield=70.0)

    if result.success:
        data = result.data

        if data['reactions_found'] > 0:
            print(f"\n✓ Analyzed {data['reactions_found']} high-yield amide formations")

            # Analyze temperatures
            temps = [
                r['conditions']['temperature']
                for r in data['reactions']
                if r['conditions']['temperature'] is not None
            ]

            if temps:
                avg_temp = sum(temps) / len(temps)
                print(f"\n  Average temperature: {avg_temp:.1f} °C")
                print(f"  Temperature range: {min(temps):.1f} - {max(temps):.1f} °C")

            # Common solvents
            print(f"\n  Most common solvents:")
            for solvent in data['common_solvents'][:5]:
                print(f"    - {solvent}")

            # Common coupling reagents
            print(f"\n  Common coupling reagents:")
            for reagent in data['common_reagents'][:5]:
                print(f"    - {reagent}")
        else:
            print(f"\n  No high-yield reactions found. Try lowering yield threshold.")
    else:
        print(f"✗ Error: {result.error}")


async def example_5_batch_processing():
    """
    Example 5: Batch process multiple compounds
    """
    print("\n" + "="*70)
    print("Example 5: Batch processing of multiple targets")
    print("="*70)

    adapter = ORDAdapter()

    # List of target molecules
    targets = [
        ("Ethanol", "CCO"),
        ("Acetone", "CC(=O)C"),
        ("Acetic Acid", "CC(=O)O"),
        ("Benzene", "c1ccccc1")
    ]

    print(f"\nProcessing {len(targets)} compounds...\n")

    # Process concurrently
    tasks = [
        adapter(smiles, use_cache=True)  # Use __call__ with caching
        for name, smiles in targets
    ]

    results = await asyncio.gather(*tasks)

    # Summary
    summary = []
    for (name, smiles), result in zip(targets, results):
        if result.success:
            n_reactions = result.data['reactions_found']
            avg_yield = result.data['statistics']['average_yield']
            summary.append({
                'name': name,
                'reactions': n_reactions,
                'avg_yield': avg_yield
            })

    print("Results Summary:")
    print(f"{'Compound':<15} {'Reactions':<12} {'Avg Yield':<12}")
    print("-" * 40)
    for item in summary:
        print(f"{item['name']:<15} {item['reactions']:<12} {item['avg_yield']:<12.1f}%")


async def example_6_retrosynthesis_validation():
    """
    Example 6: Validate a proposed retrosynthetic disconnection
    """
    print("\n" + "="*70)
    print("Example 6: Validate retrosynthetic proposal")
    print("="*70)

    adapter = ORDAdapter()

    # Proposed synthesis
    product = "c1ccc(-c2ccccc2)cc1"  # Biphenyl
    proposed_reactants = ["c1ccccc1Br", "c1ccccc1B(O)O"]  # Bromobenzene + phenylboronic acid

    print(f"\nProposed: Bromobenzene + Phenylboronic acid → Biphenyl")

    result = await adapter.execute(product)

    if result.success:
        data = result.data

        # Check if proposed reactants appear in any reactions
        validated_routes = []
        for reaction in data['reactions']:
            reactant_smiles = [r['smiles'] for r in reaction['reactants']]

            # Check for similarity (exact match may not work due to canonicalization)
            if any('Br' in s for s in reactant_smiles) and any('B' in s for s in reactant_smiles):
                validated_routes.append(reaction)

        if validated_routes:
            print(f"\n✓ Found {len(validated_routes)} literature precedents!")

            best_route = max(validated_routes, key=lambda r: r['yield'])
            print(f"\n  Best literature example:")
            print(f"    Yield: {best_route['yield']}%")
            print(f"    Temperature: {best_route['conditions']['temperature']} °C")
            print(f"    Catalyst: {', '.join(best_route['conditions']['reagents'][:2])}")
            if best_route['literature']['doi']:
                print(f"    Reference: {best_route['literature']['doi']}")
        else:
            print(f"\n  No exact match found, but {data['reactions_found']} routes exist")
            print(f"  Consider alternative disconnections")
    else:
        print(f"✗ Error: {result.error}")


async def example_7_caching_demonstration():
    """
    Example 7: Demonstrate caching functionality
    """
    print("\n" + "="*70)
    print("Example 7: Caching demonstration")
    print("="*70)

    adapter = ORDAdapter()

    smiles = "CCO"  # Ethanol

    # First call - API hit
    print("\nFirst call (API hit)...")
    import time
    start = time.time()
    result1 = await adapter(smiles, use_cache=True)
    time1 = time.time() - start

    print(f"  Time: {time1:.2f}s")
    print(f"  Cache hit: {result1.cache_hit}")
    print(f"  Reactions found: {result1.data['reactions_found']}")

    # Second call - cache hit
    print("\nSecond call (cache hit)...")
    start = time.time()
    result2 = await adapter(smiles, use_cache=True)
    time2 = time.time() - start

    print(f"  Time: {time2:.2f}s")
    print(f"  Cache hit: {result2.cache_hit}")
    print(f"  Speedup: {time1/time2:.1f}x faster")

    # Different parameters - new API call
    print("\nThird call with different parameters (API hit)...")
    start = time.time()
    result3 = await adapter(smiles, use_cache=True, min_yield=80.0)
    time3 = time.time() - start

    print(f"  Time: {time3:.2f}s")
    print(f"  Cache hit: {result3.cache_hit}")
    print(f"  Reactions found: {result3.data['reactions_found']}")


async def main():
    """
    Run all examples
    """
    print("\n" + "="*70)
    print("ORD ADAPTER EXAMPLES FOR PHARMFORGE")
    print("="*70)

    examples = [
        ("Basic Product Search", example_1_basic_product_search),
        ("Search by Reactant", example_2_search_by_reactant),
        ("Catalyst Search", example_3_catalyst_search),
        ("Condition Optimization", example_4_condition_optimization),
        ("Batch Processing", example_5_batch_processing),
        ("Retrosynthesis Validation", example_6_retrosynthesis_validation),
        ("Caching Demo", example_7_caching_demonstration),
    ]

    # Run all examples
    for name, example_func in examples:
        try:
            await example_func()
        except Exception as e:
            print(f"\n✗ Example '{name}' failed: {e}")

    print("\n" + "="*70)
    print("ALL EXAMPLES COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
