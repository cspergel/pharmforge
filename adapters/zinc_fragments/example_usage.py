"""
ZINC Fragments Adapter - Example Usage
Demonstrates fragment search, similarity search, substructure search, and purchasability lookups
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter


async def example_similarity_search():
    """
    Example 1: Find similar fragments to a query compound
    Use case: Fragment-based drug design
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Similarity Search for Fragments")
    print("="*80)

    adapter = ZINCFragmentsAdapter()

    # Benzene ring - a common fragment
    query_smiles = "c1ccccc1"

    print(f"\nSearching for fragments similar to: {query_smiles}")
    print("Use case: Finding purchasable fragments for fragment-based design")

    result = await adapter(
        query_smiles,
        search_type="similarity",
        similarity_threshold=0.8,
        max_results=20,
        mw_max=250  # Typical fragment size
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} fragments (filtered from {data['num_total']} total)")
        print(f"  Filters: MW ≤ {data['filters_applied']['mw_max']} Da")

        if data['fragments']:
            print("\nTop 5 Fragments:")
            for i, frag in enumerate(data['fragments'][:5], 1):
                print(f"\n  {i}. ZINC ID: {frag.get('zinc_id', 'N/A')}")
                print(f"     MW: {frag.get('mw', 'N/A')} Da")
                print(f"     LogP: {frag.get('logp', 'N/A')}")
                print(f"     SMILES: {frag.get('smiles', 'N/A')[:50]}...")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_substructure_search():
    """
    Example 2: Find fragments containing a specific substructure
    Use case: Scaffold hopping and bioisostere replacement
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Substructure Search")
    print("="*80)

    adapter = ZINCFragmentsAdapter()

    # Pyridine substructure - common in drugs
    substructure_smiles = "c1ccncc1"

    print(f"\nSearching for fragments containing: {substructure_smiles} (pyridine)")
    print("Use case: Finding pyridine-containing fragments for bioisostere replacement")

    result = await adapter(
        substructure_smiles,
        search_type="substructure",
        max_results=15,
        mw_max=300,
        logp_max=3.0
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} fragments")
        print(f"  Filters: MW ≤ 300 Da, LogP ≤ 3.0")

        if data['fragments']:
            print("\nTop 5 Results:")
            for i, frag in enumerate(data['fragments'][:5], 1):
                print(f"\n  {i}. {frag.get('zinc_id', 'N/A')}")
                print(f"     Properties: MW={frag.get('mw', 'N/A')}, LogP={frag.get('logp', 'N/A')}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_purchasability_search():
    """
    Example 3: Find purchasable fragments with vendor information
    Use case: Building block selection for synthesis planning
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Purchasability and Vendor Lookup")
    print("="*80)

    adapter = ZINCFragmentsAdapter()

    # Naphthalene - useful scaffold
    query_smiles = "c1ccc2ccccc2c1"

    print(f"\nSearching for purchasable fragments similar to: {query_smiles}")
    print("Use case: Identifying commercially available building blocks")

    result = await adapter(
        query_smiles,
        search_type="similarity",
        similarity_threshold=0.75,
        max_results=10,
        get_purchasability=True,  # Get vendor info
        mw_max=250
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} fragments")

        # Show purchasability info
        purchasable = [f for f in data['fragments'] if f.get('purchasability', {}).get('purchasable')]
        print(f"  Purchasable: {len(purchasable)}")

        if purchasable:
            print("\nTop Purchasable Fragments:")
            for i, frag in enumerate(purchasable[:5], 1):
                purch = frag.get('purchasability', {})
                print(f"\n  {i}. {frag.get('zinc_id', 'N/A')}")
                print(f"     Vendors: {purch.get('vendor_count', 0)}")
                if purch.get('suppliers'):
                    print(f"     Available from: {', '.join(purch['suppliers'][:3])}")
                if purch.get('price_range'):
                    print(f"     Price range: {purch['price_range']}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_property_filtering():
    """
    Example 4: Advanced property filtering for fragment selection
    Use case: Rule-of-Three compliance for fragment libraries
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Property-Based Fragment Filtering (Rule of Three)")
    print("="*80)

    adapter = ZINCFragmentsAdapter()

    query_smiles = "CCc1ccccc1"  # Ethylbenzene

    print(f"\nSearching for Rule-of-Three compliant fragments")
    print("Rule of Three: MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3")

    result = await adapter(
        query_smiles,
        search_type="similarity",
        similarity_threshold=0.7,
        max_results=50,
        mw_max=300,
        logp_max=3.0,
        hbd_max=3,
        hba_max=3
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} Rule-of-Three compliant fragments")
        print(f"  (from {data['num_total']} total hits)")

        print("\nApplied Filters:")
        for key, value in data['filters_applied'].items():
            if value is not None:
                print(f"  - {key}: ≤ {value}")

        if data['fragments']:
            print("\nSample Fragments:")
            for i, frag in enumerate(data['fragments'][:3], 1):
                print(f"\n  {i}. {frag.get('zinc_id', 'N/A')}")
                print(f"     MW: {frag.get('mw', 'N/A')} Da")
                print(f"     LogP: {frag.get('logp', 'N/A')}")
                print(f"     HBD/HBA: {frag.get('hbd', 'N/A')}/{frag.get('hba', 'N/A')}")
    else:
        print(f"\n✗ Error: {result.error}")


async def main():
    """Run all examples"""
    print("\n" + "#"*80)
    print("# ZINC Fragments Adapter - Example Usage")
    print("# Fragment-based drug discovery with purchasability information")
    print("#"*80)

    await example_similarity_search()
    await example_substructure_search()
    await example_purchasability_search()
    await example_property_filtering()

    print("\n" + "#"*80)
    print("# Examples Complete!")
    print("#"*80)
    print("\nKey Features Demonstrated:")
    print("  ✓ Similarity search for fragment discovery")
    print("  ✓ Substructure search for scaffold hopping")
    print("  ✓ Purchasability and vendor lookups")
    print("  ✓ Property-based filtering (Rule of Three)")
    print("\nUse Cases:")
    print("  • Fragment-based drug design")
    print("  • Building block selection for synthesis")
    print("  • Scaffold hopping and bioisostere replacement")
    print("  • Fragment library curation")


if __name__ == "__main__":
    asyncio.run(main())
