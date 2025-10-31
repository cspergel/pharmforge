"""
COCONUT Adapter - Example Usage

Demonstrates various ways to search and retrieve natural product data
from the COCONUT database.
"""
import asyncio
import json
from adapters.coconut.adapter import COCONUTAdapter


async def example_simple_smiles_search():
    """
    Example 1: Simple SMILES search
    Search for natural products similar to a given structure
    """
    print("\n" + "="*70)
    print("Example 1: Simple SMILES Search")
    print("="*70)

    adapter = COCONUTAdapter()

    # Caffeine SMILES
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    result = await adapter(smiles)

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} natural products similar to caffeine")

        if data['results']:
            first = data['results'][0]
            print(f"\nFirst result:")
            print(f"  ID: {first.get('coconut_id')}")
            print(f"  Name: {first.get('name')}")
            print(f"  Formula: {first.get('molecular_formula')}")
            print(f"  MW: {first.get('molecular_weight')}")

            if 'natural_source' in first:
                source = first['natural_source']
                print(f"  Source: {source.get('organism')} ({source.get('common_name')})")
    else:
        print(f"Error: {result.error}")


async def example_advanced_query():
    """
    Example 2: Advanced query with filters
    Search with molecular property constraints
    """
    print("\n" + "="*70)
    print("Example 2: Advanced Query with Property Filters")
    print("="*70)

    adapter = COCONUTAdapter()

    query = {
        "query_type": "smiles",
        "query": "CC(C)=CCC1C(=O)CCC1C",  # Terpene-like structure
        "filters": {
            "molecular_weight": {"min": 200, "max": 500},
            "logp": {"min": 1, "max": 5}
        },
        "include_taxonomy": True,
        "include_activities": True,
        "max_results": 10
    }

    result = await adapter(query)

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} terpenes matching criteria")

        for i, compound in enumerate(data['results'][:3], 1):
            print(f"\n{i}. {compound.get('name', 'Unknown')}")
            print(f"   MW: {compound.get('molecular_weight'):.2f}")
            props = compound.get('properties', {})
            print(f"   LogP: {props.get('logp')}")

            if 'natural_source' in compound:
                source = compound['natural_source']
                print(f"   From: {source.get('organism')}")
    else:
        print(f"Error: {result.error}")


async def example_organism_search():
    """
    Example 3: Search by organism
    Find all compounds from a specific source organism
    """
    print("\n" + "="*70)
    print("Example 3: Search by Organism")
    print("="*70)

    adapter = COCONUTAdapter()

    query = {
        "query_type": "organism",
        "query": "Streptomyces",
        "max_results": 20
    }

    result = await adapter(query)

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} compounds from Streptomyces")

        # Summarize by molecular weight
        weights = [c.get('molecular_weight', 0) for c in data['results'] if c.get('molecular_weight')]
        if weights:
            print(f"  MW range: {min(weights):.1f} - {max(weights):.1f}")
            print(f"  Average MW: {sum(weights)/len(weights):.1f}")

        # Show first few
        for compound in data['results'][:3]:
            print(f"\n  - {compound.get('name')}")
            print(f"    SMILES: {compound.get('smiles')}")
    else:
        print(f"Error: {result.error}")


async def example_name_search():
    """
    Example 4: Search by compound name
    Find natural products by common or systematic name
    """
    print("\n" + "="*70)
    print("Example 4: Search by Name")
    print("="*70)

    adapter = COCONUTAdapter()

    query = {
        "query_type": "name",
        "query": "Taxol",
        "include_taxonomy": True,
        "include_activities": True
    }

    result = await adapter(query)

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} compounds matching 'Taxol'")

        for compound in data['results']:
            print(f"\n  Name: {compound.get('name')}")
            print(f"  COCONUT ID: {compound.get('coconut_id')}")
            print(f"  Formula: {compound.get('molecular_formula')}")

            if 'natural_source' in compound:
                source = compound['natural_source']
                print(f"  Source: {source.get('organism')}")
                print(f"  Taxonomy: {source.get('taxonomy')}")

            if 'biological_activities' in compound:
                activities = compound['biological_activities']
                print(f"  Activities: {', '.join(activities[:3])}")
    else:
        print(f"Error: {result.error}")


async def example_get_by_id():
    """
    Example 5: Get compound by COCONUT ID
    Retrieve detailed information for a specific compound
    """
    print("\n" + "="*70)
    print("Example 5: Get Compound by ID")
    print("="*70)

    adapter = COCONUTAdapter()

    query = {
        "query_type": "id",
        "query": "CNP0123456"  # Example COCONUT ID
    }

    result = await adapter(query)

    if result.success:
        data = result.data
        if data['num_results'] > 0:
            compound = data['results'][0]
            print(f"\nCompound Details:")
            print(json.dumps(compound, indent=2))
        else:
            print("Compound not found")
    else:
        print(f"Error: {result.error}")


async def example_natural_product_library():
    """
    Example 6: Build a natural product library
    Extract compounds with specific characteristics
    """
    print("\n" + "="*70)
    print("Example 6: Build Natural Product Library")
    print("="*70)

    adapter = COCONUTAdapter()

    # Search for drug-like natural products from plants
    query = {
        "query_type": "organism",
        "query": "Plantae",
        "filters": {
            "molecular_weight": {"min": 200, "max": 600},
            "logp": {"min": 0, "max": 5}
        },
        "max_results": 50
    }

    result = await adapter(query)

    if result.success:
        data = result.data
        print(f"\nBuilt library of {data['num_results']} drug-like plant compounds")

        # Analyze Lipinski's Rule of Five compliance
        ro5_compliant = []
        for compound in data['results']:
            props = compound.get('properties', {})
            mw = compound.get('molecular_weight', 0)
            logp = props.get('logp', 0)
            hbd = props.get('hbd', 0)
            hba = props.get('hba', 0)

            if mw and logp is not None and hbd is not None and hba is not None:
                if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                    ro5_compliant.append(compound)

        print(f"  Lipinski Ro5 compliant: {len(ro5_compliant)}/{data['num_results']}")
        print(f"  Compliance rate: {len(ro5_compliant)/max(data['num_results'],1)*100:.1f}%")

        # Show top candidates
        print("\n  Top 3 drug-like natural products:")
        for i, compound in enumerate(ro5_compliant[:3], 1):
            print(f"\n  {i}. {compound.get('name')}")
            print(f"     MW: {compound.get('molecular_weight'):.1f}")
            props = compound.get('properties', {})
            print(f"     LogP: {props.get('logp'):.2f}")
            print(f"     HBD: {props.get('hbd')}, HBA: {props.get('hba')}")
    else:
        print(f"Error: {result.error}")


async def example_scaffold_mining():
    """
    Example 7: Scaffold mining from nature
    Find natural products with similar scaffolds
    """
    print("\n" + "="*70)
    print("Example 7: Scaffold Mining from Natural Products")
    print("="*70)

    adapter = COCONUTAdapter()

    # Search for compounds similar to aspirin (acetylsalicylic acid)
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"

    result = await adapter(aspirin_smiles)

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} natural products with aspirin-like scaffolds")

        # Group by source kingdom
        kingdoms = {}
        for compound in data['results']:
            source = compound.get('natural_source', {})
            taxonomy = source.get('taxonomy', 'Unknown')
            kingdom = taxonomy.split(';')[0] if ';' in taxonomy else taxonomy

            if kingdom not in kingdoms:
                kingdoms[kingdom] = []
            kingdoms[kingdom].append(compound)

        print("\n  Distribution by kingdom:")
        for kingdom, compounds in kingdoms.items():
            print(f"    {kingdom}: {len(compounds)} compounds")
    else:
        print(f"Error: {result.error}")


async def main():
    """Run all examples"""
    print("\n" + "="*70)
    print("COCONUT Adapter - Usage Examples")
    print("Natural Products Database Integration")
    print("="*70)

    # Run each example
    await example_simple_smiles_search()
    await example_advanced_query()
    await example_organism_search()
    await example_name_search()
    await example_get_by_id()
    await example_natural_product_library()
    await example_scaffold_mining()

    print("\n" + "="*70)
    print("All examples completed!")
    print("="*70 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
