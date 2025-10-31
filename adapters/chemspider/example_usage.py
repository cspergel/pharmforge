"""
ChemSpider Adapter - Example Usage
Demonstrates various ways to use the ChemSpider adapter
"""
import asyncio
import json
from adapters.chemspider import ChemSpiderAdapter


async def example_1_basic_smiles_search():
    """
    Example 1: Basic SMILES search
    Search for a compound by SMILES string
    """
    print("\n" + "="*60)
    print("Example 1: Basic SMILES Search (Aspirin)")
    print("="*60)

    adapter = ChemSpiderAdapter()

    # Aspirin SMILES
    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} compound(s)")

        if data['results']:
            compound = data['results'][0]
            print(f"\nChemSpider ID: {compound['chemspider_id']}")
            print(f"Common Name: {compound['common_name']}")
            print(f"Formula: {compound['formula']}")
            print(f"Molecular Weight: {compound['molecular_weight']}")
            print(f"\nSMILES: {compound['smiles']}")
            print(f"InChI: {compound['inchi'][:60]}...")
            print(f"InChIKey: {compound['inchikey']}")

            print(f"\nSynonyms ({len(compound['synonyms'])} total):")
            for syn in compound['synonyms'][:5]:
                print(f"  - {syn}")

            print(f"\nData Sources ({len(compound['data_sources'])} total):")
            for source in compound['data_sources'][:5]:
                print(f"  - {source}")

            if compound['properties']:
                print(f"\nProperties:")
                for key, value in compound['properties'].items():
                    if value is not None:
                        print(f"  {key}: {value}")
    else:
        print(f"Error: {result.error}")


async def example_2_name_search():
    """
    Example 2: Search by compound name
    """
    print("\n" + "="*60)
    print("Example 2: Search by Name (Ibuprofen)")
    print("="*60)

    adapter = ChemSpiderAdapter()

    result = await adapter.execute({
        "query": "ibuprofen",
        "query_type": "name",
        "max_results": 3
    })

    if result.success:
        print(f"\nFound {result.data['num_results']} compound(s)")

        for i, compound in enumerate(result.data['results'], 1):
            print(f"\nResult {i}:")
            print(f"  ChemSpider ID: {compound['chemspider_id']}")
            print(f"  Name: {compound['common_name']}")
            print(f"  Formula: {compound['formula']}")
            print(f"  SMILES: {compound['smiles']}")
            print(f"  Data sources: {len(compound['data_sources'])}")
    else:
        print(f"Error: {result.error}")


async def example_3_formula_search():
    """
    Example 3: Search by molecular formula
    """
    print("\n" + "="*60)
    print("Example 3: Search by Molecular Formula (C9H8O4)")
    print("="*60)

    adapter = ChemSpiderAdapter()

    result = await adapter.execute({
        "query": "C9H8O4",
        "query_type": "formula",
        "max_results": 5
    })

    if result.success:
        print(f"\nFound {result.data['num_results']} compounds with formula C9H8O4")
        print(f"Total in database: {result.data['total_found']}")

        for compound in result.data['results']:
            print(f"\n  - {compound['common_name']} (ID: {compound['chemspider_id']})")
            print(f"    SMILES: {compound['smiles']}")
            print(f"    MW: {compound['molecular_weight']}")
    else:
        print(f"Error: {result.error}")


async def example_4_inchikey_search():
    """
    Example 4: Search by InChIKey
    """
    print("\n" + "="*60)
    print("Example 4: Search by InChIKey (Caffeine)")
    print("="*60)

    adapter = ChemSpiderAdapter()

    # Caffeine InChIKey
    inchikey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"

    result = await adapter.execute({
        "query": inchikey,
        "query_type": "inchikey"
    })

    if result.success and result.data['results']:
        compound = result.data['results'][0]
        print(f"\nCompound: {compound['common_name']}")
        print(f"ChemSpider ID: {compound['chemspider_id']}")
        print(f"Formula: {compound['formula']}")
        print(f"SMILES: {compound['smiles']}")
        print(f"\nAggregated from {len(compound['data_sources'])} databases:")
        for source in compound['data_sources']:
            print(f"  - {source}")
    else:
        print(f"Error or no results: {result.error if not result.success else 'No results'}")


async def example_5_chemical_validation():
    """
    Example 5: Validate a chemical structure across multiple databases
    """
    print("\n" + "="*60)
    print("Example 5: Chemical Validation")
    print("="*60)

    adapter = ChemSpiderAdapter()

    # Test compound - Acetaminophen (Paracetamol)
    smiles = "CC(=O)Nc1ccc(O)cc1"

    result = await adapter.execute(smiles)

    if result.success and result.data['results']:
        compound = result.data['results'][0]
        num_sources = len(compound['data_sources'])

        print(f"\nValidating: {smiles}")
        print(f"Identified as: {compound['common_name']}")
        print(f"Found in {num_sources} database sources")

        # Validation confidence
        if num_sources >= 5:
            confidence = "HIGH"
            symbol = "✓"
        elif num_sources >= 3:
            confidence = "MEDIUM"
            symbol = "○"
        else:
            confidence = "LOW"
            symbol = "!"

        print(f"\nValidation Confidence: {symbol} {confidence}")
        print(f"Reasoning: Structure found in {num_sources} authoritative sources")

        # Check name consistency
        print(f"\nKnown synonyms ({len(compound['synonyms'])} total):")
        for syn in compound['synonyms'][:10]:
            print(f"  - {syn}")

        # Identifier verification
        print(f"\nStandardized Identifiers:")
        print(f"  SMILES: {compound['smiles']}")
        print(f"  InChIKey: {compound['inchikey']}")
        print(f"  Formula: {compound['formula']}")
    else:
        print(f"Validation failed: {result.error if not result.success else 'No results'}")


async def example_6_batch_processing():
    """
    Example 6: Batch processing multiple compounds
    """
    print("\n" + "="*60)
    print("Example 6: Batch Processing")
    print("="*60)

    adapter = ChemSpiderAdapter()

    compounds = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Acetaminophen", "CC(=O)Nc1ccc(O)cc1")
    ]

    results = []

    print(f"\nProcessing {len(compounds)} compounds...")

    for name, smiles in compounds:
        print(f"\n  Processing: {name}...")

        result = await adapter.execute(smiles)

        if result.success and result.data['results']:
            compound = result.data['results'][0]
            results.append({
                "input_name": name,
                "smiles": smiles,
                "chemspider_name": compound['common_name'],
                "chemspider_id": compound['chemspider_id'],
                "formula": compound['formula'],
                "molecular_weight": compound['molecular_weight'],
                "num_sources": len(compound['data_sources']),
                "num_synonyms": len(compound['synonyms'])
            })
            print(f"    ✓ Found: {compound['common_name']} (ID: {compound['chemspider_id']})")
        else:
            print(f"    ✗ Not found or error")

        # Respect rate limits
        await asyncio.sleep(1)

    # Summary table
    print("\n" + "="*60)
    print("BATCH PROCESSING SUMMARY")
    print("="*60)
    print(f"\n{'Name':<15} {'Formula':<10} {'MW':<8} {'Sources':<8} {'Synonyms':<10}")
    print("-" * 60)

    for r in results:
        print(f"{r['input_name']:<15} {r['formula']:<10} {r['molecular_weight']:<8.2f} "
              f"{r['num_sources']:<8} {r['num_synonyms']:<10}")

    print(f"\nTotal processed: {len(results)}/{len(compounds)}")


async def example_7_property_comparison():
    """
    Example 7: Compare properties from ChemSpider vs other sources
    """
    print("\n" + "="*60)
    print("Example 7: Property Comparison")
    print("="*60)

    adapter = ChemSpiderAdapter()

    # Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    result = await adapter.execute(smiles)

    if result.success and result.data['results']:
        compound = result.data['results'][0]

        print(f"\nCompound: {compound['common_name']}")
        print(f"ChemSpider aggregated properties:")

        props = compound['properties']
        print(f"\n  Molecular Weight: {compound['molecular_weight']}")
        print(f"  Formula: {compound['formula']}")

        if props.get('alogp') is not None:
            print(f"  ALogP: {props['alogp']:.2f}")
        if props.get('xlogp') is not None:
            print(f"  XLogP: {props['xlogp']:.2f}")
        if props.get('monoisotopic_mass') is not None:
            print(f"  Monoisotopic Mass: {props['monoisotopic_mass']:.4f}")

        print(f"\n  Aggregated from {len(compound['data_sources'])} sources:")
        for source in compound['data_sources'][:8]:
            print(f"    - {source}")


async def example_8_error_handling():
    """
    Example 8: Proper error handling
    """
    print("\n" + "="*60)
    print("Example 8: Error Handling")
    print("="*60)

    adapter = ChemSpiderAdapter()

    # Test with invalid SMILES
    invalid_smiles = "INVALID_SMILES_STRING_12345"

    print(f"\nTesting with invalid SMILES: {invalid_smiles}")

    try:
        result = await adapter.execute(invalid_smiles)

        if not result.success:
            print(f"\nSearch failed (expected):")
            print(f"  Error: {result.error}")
            print(f"  Metadata: {result.metadata}")
        elif result.data['num_results'] == 0:
            print(f"\nNo results found (expected):")
            print(f"  Warnings: {result.data['warnings']}")
        else:
            print(f"\nUnexpectedly found results")

    except Exception as e:
        print(f"\nException caught: {type(e).__name__}: {e}")

    # Test with valid but obscure compound
    print(f"\n\nTesting with valid but obscure SMILES:")
    obscure_smiles = "C1=CC=C2C(=C1)C=CC=N2"  # Quinoline

    result = await adapter.execute(obscure_smiles)

    if result.success:
        if result.data['num_results'] > 0:
            print(f"  Found {result.data['num_results']} result(s)")
            print(f"  Name: {result.data['results'][0]['common_name']}")
        else:
            print(f"  No results found")
            if result.data['warnings']:
                print(f"  Warnings: {result.data['warnings']}")


async def example_9_caching_demo():
    """
    Example 9: Demonstrate caching behavior
    """
    print("\n" + "="*60)
    print("Example 9: Caching Demonstration")
    print("="*60)

    adapter = ChemSpiderAdapter()

    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

    print(f"\nFirst call (should hit API):")
    result1 = await adapter(smiles, use_cache=True)
    print(f"  Success: {result1.success}")
    print(f"  Cache hit: {result1.cache_hit}")
    print(f"  Cache key: {result1.metadata.get('cache_key', 'N/A')[:32]}...")

    print(f"\nSecond call (should use cache):")
    result2 = await adapter(smiles, use_cache=True)
    print(f"  Success: {result2.success}")
    print(f"  Cache hit: {result2.cache_hit}")
    print(f"  Same data: {result1.data == result2.data}")

    print(f"\nThird call (bypass cache):")
    result3 = await adapter(smiles, use_cache=False)
    print(f"  Success: {result3.success}")
    print(f"  Cache hit: {result3.cache_hit}")


async def main():
    """
    Run all examples
    """
    print("\n" + "="*60)
    print("ChemSpider Adapter - Usage Examples")
    print("="*60)
    print("\nNOTE: These examples require a valid CHEMSPIDER_API_KEY")
    print("Set it via: export CHEMSPIDER_API_KEY='your-key-here'")
    print("Get your free key at: https://developer.rsc.org/")

    import os
    if not os.getenv("CHEMSPIDER_API_KEY"):
        print("\n⚠ WARNING: CHEMSPIDER_API_KEY not found in environment")
        print("Examples may fail or have limited functionality\n")
        response = input("Continue anyway? (y/n): ")
        if response.lower() != 'y':
            return

    # Run examples
    try:
        await example_1_basic_smiles_search()
        await asyncio.sleep(1)

        await example_2_name_search()
        await asyncio.sleep(1)

        await example_3_formula_search()
        await asyncio.sleep(1)

        await example_4_inchikey_search()
        await asyncio.sleep(1)

        await example_5_chemical_validation()
        await asyncio.sleep(1)

        await example_6_batch_processing()
        await asyncio.sleep(1)

        await example_7_property_comparison()
        await asyncio.sleep(1)

        await example_8_error_handling()
        await asyncio.sleep(1)

        await example_9_caching_demo()

        print("\n" + "="*60)
        print("All examples completed!")
        print("="*60)

    except KeyboardInterrupt:
        print("\n\nExamples interrupted by user")
    except Exception as e:
        print(f"\n\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(main())
