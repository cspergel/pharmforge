"""
BRENDA Adapter - Example Usage
Demonstrates various ways to query enzyme kinetics and inhibitor data
"""
import asyncio
from adapters.brenda import BRENDAAdapter


async def example_basic_query():
    """Example 1: Basic EC number query"""
    print("=" * 60)
    print("Example 1: Basic EC Number Query")
    print("=" * 60)

    adapter = BRENDAAdapter()

    # Simple query by EC number
    result = await adapter.execute("1.1.1.1")

    if result.success:
        data = result.data

        print(f"\nEnzyme: {data['enzyme']['common_name']}")
        print(f"EC Number: {data['enzyme']['ec_number']}")
        print(f"Systematic Name: {data['enzyme']['systematic_name']}")
        print(f"Reaction: {data['enzyme']['reaction']}")
        print(f"\nSubstrates: {', '.join(data['enzyme']['substrates'])}")
        print(f"Products: {', '.join(data['enzyme']['products'])}")
        print(f"Cofactors: {', '.join(data['enzyme']['cofactors'])}")

        print(f"\nKinetic Parameters Found: {data['num_kinetic_records']}")
        print(f"Inhibitors Found: {data['num_inhibitor_records']}")

        if data['warnings']:
            print(f"\nWarnings: {', '.join(data['warnings'])}")
    else:
        print(f"Error: {result.error}")


async def example_kinetic_parameters():
    """Example 2: Query specific kinetic parameters"""
    print("\n" + "=" * 60)
    print("Example 2: Kinetic Parameters Query")
    print("=" * 60)

    adapter = BRENDAAdapter()

    # Query with specific parameters
    query = {
        "query_type": "ec_number",
        "query": "1.1.1.1",
        "parameters": ["km", "kcat"],
        "organisms": ["Homo sapiens"],
        "max_results": 50
    }

    result = await adapter.execute(query)

    if result.success:
        kinetics = result.data['kinetics']

        print(f"\nFound {len(kinetics)} kinetic parameters for human enzyme\n")

        # Group by parameter type
        km_values = [k for k in kinetics if k['parameter'] == 'km']
        kcat_values = [k for k in kinetics if k['parameter'] == 'kcat']

        print("Km Values (Michaelis constant):")
        for km in km_values:
            print(f"  - {km['substrate']}: {km['value']} {km['unit']}")
            print(f"    pH {km['conditions']['pH']}, {km['conditions']['temperature']}°C")
            print(f"    Reference: {km['reference']}\n")

        print("Kcat Values (Turnover number):")
        for kcat in kcat_values:
            print(f"  - {kcat['substrate']}: {kcat['value']} {kcat['unit']}")
            print(f"    Reference: {kcat['reference']}\n")

        # Calculate catalytic efficiency
        if km_values and kcat_values:
            print("Catalytic Efficiency (Kcat/Km):")
            for km in km_values:
                kcat = next(
                    (k for k in kcat_values if k['substrate'] == km['substrate']),
                    None
                )
                if kcat:
                    efficiency = kcat['value'] / km['value']
                    print(f"  - {km['substrate']}: {efficiency:.2f} mM⁻¹s⁻¹")


async def example_inhibitor_data():
    """Example 3: Query inhibitor information"""
    print("\n" + "=" * 60)
    print("Example 3: Inhibitor Data Query")
    print("=" * 60)

    adapter = BRENDAAdapter()

    query = {
        "query": "1.1.1.1",
        "parameters": ["ki"],
        "organisms": ["Homo sapiens"]
    }

    result = await adapter.execute(query)

    if result.success:
        inhibitors = result.data['inhibitors']

        print(f"\nFound {len(inhibitors)} inhibitors\n")

        # Sort by potency (lower Ki = more potent)
        inhibitors_sorted = sorted(inhibitors, key=lambda x: x['ki_value'])

        print("Inhibitors (sorted by potency):")
        for i, inh in enumerate(inhibitors_sorted, 1):
            print(f"{i}. {inh['compound']}")
            print(f"   Ki: {inh['ki_value']} {inh['unit']}")
            print(f"   Type: {inh['inhibition_type']}")
            print(f"   Organism: {inh['organism']}")
            print(f"   Reference: {inh['reference']}\n")

        # Find competitive inhibitors with Ki < 100 nM
        potent_competitive = [
            i for i in inhibitors
            if i['inhibition_type'] == 'competitive' and i['ki_value'] < 100
        ]

        if potent_competitive:
            print(f"Found {len(potent_competitive)} potent competitive inhibitors (Ki < 100 nM)")
            print("These are excellent lead compounds for drug design!")


async def example_multi_organism():
    """Example 4: Compare enzyme kinetics across organisms"""
    print("\n" + "=" * 60)
    print("Example 4: Multi-Organism Comparison")
    print("=" * 60)

    adapter = BRENDAAdapter()

    organisms = [
        "Homo sapiens",
        "Mus musculus",
        "Saccharomyces cerevisiae"
    ]

    print(f"\nComparing enzyme 1.1.1.1 (Alcohol dehydrogenase) across species:\n")

    for organism in organisms:
        query = {
            "query": "1.1.1.1",
            "parameters": ["km"],
            "organisms": [organism],
            "max_results": 10
        }

        result = await adapter.execute(query)

        if result.success:
            kinetics = result.data['kinetics']

            if kinetics:
                # Average Km for this organism
                km_values = [k['value'] for k in kinetics if k['parameter'] == 'km']
                avg_km = sum(km_values) / len(km_values)

                print(f"{organism}:")
                print(f"  Data points: {len(km_values)}")
                print(f"  Average Km: {avg_km:.2f} mM")
                print(f"  Range: {min(km_values):.2f} - {max(km_values):.2f} mM\n")
            else:
                print(f"{organism}: No data available\n")


async def example_drug_target_validation():
    """Example 5: Drug target validation workflow"""
    print("\n" + "=" * 60)
    print("Example 5: Drug Target Validation")
    print("=" * 60)

    adapter = BRENDAAdapter()

    # Target enzyme: Renin (3.4.23.15)
    ec_number = "1.1.1.1"  # Using ADH as example
    print(f"\nValidating target: EC {ec_number}\n")

    query = {
        "query": ec_number,
        "parameters": ["km", "ki"],
        "organisms": ["Homo sapiens"]
    }

    result = await adapter.execute(query)

    if result.success:
        enzyme = result.data['enzyme']
        kinetics = result.data['kinetics']
        inhibitors = result.data['inhibitors']

        print(f"Target: {enzyme['common_name']}")
        print(f"Reaction: {enzyme['reaction']}")

        # Validation criteria
        print("\n--- Druggability Assessment ---")

        # 1. Known inhibitors
        print(f"\n1. Known Inhibitors: {len(inhibitors)}")
        if len(inhibitors) >= 5:
            print("   ✓ Good: Multiple known inhibitors exist")
        elif len(inhibitors) >= 1:
            print("   ~ Moderate: Some inhibitors known")
        else:
            print("   ✗ Poor: No known inhibitors")

        # 2. Potent inhibitors
        potent = [i for i in inhibitors if i['ki_value'] < 100]  # < 100 nM
        print(f"\n2. Potent Inhibitors (Ki < 100 nM): {len(potent)}")
        if len(potent) >= 3:
            print("   ✓ Excellent: Multiple potent inhibitors")
        elif len(potent) >= 1:
            print("   ✓ Good: Potent inhibitors exist")
        else:
            print("   ~ Moderate: No highly potent inhibitors yet")

        # 3. Competitive inhibitors (best for drug design)
        competitive = [i for i in inhibitors if i['inhibition_type'] == 'competitive']
        print(f"\n3. Competitive Inhibitors: {len(competitive)}")
        if len(competitive) >= 3:
            print("   ✓ Good: Multiple competitive inhibitors")
        elif len(competitive) >= 1:
            print("   ~ Moderate: Some competitive inhibitors")
        else:
            print("   Note: Non-competitive inhibitors may still be valuable")

        # 4. Kinetic data availability
        print(f"\n4. Kinetic Data: {len(kinetics)} parameters")
        if len(kinetics) >= 5:
            print("   ✓ Excellent: Rich kinetic data available")
        elif len(kinetics) >= 2:
            print("   ✓ Good: Sufficient kinetic data")
        else:
            print("   ~ Limited kinetic data available")

        # Overall assessment
        print("\n--- Overall Assessment ---")
        score = 0
        if len(inhibitors) >= 5:
            score += 1
        if len(potent) >= 1:
            score += 1
        if len(competitive) >= 1:
            score += 1
        if len(kinetics) >= 5:
            score += 1

        if score >= 3:
            print("✓ STRONG TARGET: Excellent druggability indicators")
            print("  Recommendation: Proceed with drug design")
        elif score >= 2:
            print("~ MODERATE TARGET: Some druggability evidence")
            print("  Recommendation: Further validation recommended")
        else:
            print("? WEAK TARGET: Limited druggability evidence")
            print("  Recommendation: Consider alternative targets")

        # Show top inhibitors
        if potent:
            print("\nTop Lead Compounds:")
            for i, inh in enumerate(sorted(potent, key=lambda x: x['ki_value'])[:3], 1):
                print(f"  {i}. {inh['compound']} (Ki: {inh['ki_value']} {inh['unit']})")


async def example_caching():
    """Example 6: Demonstrate caching behavior"""
    print("\n" + "=" * 60)
    print("Example 6: Caching Demonstration")
    print("=" * 60)

    adapter = BRENDAAdapter()

    print("\nFirst query (will fetch from BRENDA)...")
    import time
    start = time.time()
    result1 = await adapter("1.1.1.1")
    time1 = time.time() - start
    print(f"Time: {time1:.2f}s")
    print(f"Cache hit: {result1.cache_hit}")

    print("\nSecond query (should use cache)...")
    start = time.time()
    result2 = await adapter("1.1.1.1")
    time2 = time.time() - start
    print(f"Time: {time2:.2f}s")
    print(f"Cache hit: {result2.cache_hit}")

    if result2.cache_hit:
        print(f"\n✓ Cache working! Speedup: {time1/time2:.1f}x faster")
    else:
        print("\nNote: Cache may not be enabled in development mode")


async def example_error_handling():
    """Example 7: Error handling"""
    print("\n" + "=" * 60)
    print("Example 7: Error Handling")
    print("=" * 60)

    adapter = BRENDAAdapter()

    # Invalid EC number
    print("\n1. Invalid EC number:")
    result = await adapter.execute("999.999.999.999")
    if not result.success:
        print(f"   Error: {result.error}")

    # Invalid query type
    print("\n2. Invalid query format:")
    result = await adapter.execute({"invalid": "data"})
    if not result.success:
        print(f"   Error: {result.error}")

    # Unsupported query type
    print("\n3. Unsupported query type:")
    result = await adapter.execute({
        "query_type": "substrate",
        "query": "ethanol"
    })
    if not result.success:
        print(f"   Error: {result.error}")

    print("\n✓ Error handling working correctly")


async def main():
    """Run all examples"""
    print("\n" + "=" * 60)
    print("BRENDA Adapter - Comprehensive Examples")
    print("=" * 60)

    examples = [
        example_basic_query,
        example_kinetic_parameters,
        example_inhibitor_data,
        example_multi_organism,
        example_drug_target_validation,
        example_caching,
        example_error_handling
    ]

    for example in examples:
        try:
            await example()
        except Exception as e:
            print(f"\nError in {example.__name__}: {e}")

    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
