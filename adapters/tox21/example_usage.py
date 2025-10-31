"""
Example Usage of Tox21 Adapter

Demonstrates various use cases for querying Tox21 toxicity screening data.
"""

import asyncio
import json
from adapters.tox21 import Tox21Adapter


async def example_1_basic_query():
    """Example 1: Basic toxicity query by SMILES"""
    print("\n" + "="*60)
    print("Example 1: Basic Toxicity Query")
    print("="*60)

    adapter = Tox21Adapter()

    # Query aspirin toxicity
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"\nQuerying Tox21 data for Aspirin: {smiles}")

    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        print(f"\n✓ Success!")
        print(f"  Compound CID: {data['compound']['cid']}")
        print(f"  Total assays tested: {data['summary']['total_assays']}")
        print(f"  Active assays: {data['summary']['active_assays']}")
        print(f"  Inactive assays: {data['summary']['inactive_assays']}")
        print(f"  Overall risk: {data['summary']['overall_risk']}")

        if data['summary']['toxicity_flags']:
            print(f"  Toxicity flags: {', '.join(data['summary']['toxicity_flags'])}")

        if data['warnings']:
            print("\n  Warnings:")
            for warning in data['warnings']:
                print(f"    - {warning}")

        # Show first 3 active assays
        active_results = [r for r in data['toxicity_results'] if r['activity'] == 'active']
        if active_results:
            print("\n  Active Assays (first 3):")
            for assay in active_results[:3]:
                print(f"    - {assay['assay_name']} ({assay['endpoint']})")
                if assay['ac50']:
                    print(f"      AC50: {assay['ac50']:.2f} uM, Confidence: {assay['confidence']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_2_query_by_name():
    """Example 2: Query by compound name"""
    print("\n" + "="*60)
    print("Example 2: Query by Compound Name")
    print("="*60)

    adapter = Tox21Adapter()

    input_data = {
        "query_type": "name",
        "query": "Bisphenol A",
        "include_inactive": False  # Only show active results
    }

    print(f"\nQuerying Tox21 data for: {input_data['query']}")
    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\n✓ Found compound (CID: {data['compound']['cid']})")
        print(f"  Active assays: {data['summary']['active_assays']}")
        print(f"  Overall risk: {data['summary']['overall_risk']}")

        if data['toxicity_results']:
            print("\n  Active Toxicity Results:")
            for assay in data['toxicity_results']:
                print(f"    - {assay['assay_name']}")
                print(f"      Category: {assay['category']}")
                print(f"      Endpoint: {assay['endpoint']}")
                if assay['ac50']:
                    print(f"      AC50: {assay['ac50']:.2f} uM")
    else:
        print(f"✗ Error: {result.error}")


async def example_3_endocrine_screening():
    """Example 3: Focused endocrine disruption screening"""
    print("\n" + "="*60)
    print("Example 3: Endocrine Disruption Screening")
    print("="*60)

    adapter = Tox21Adapter()

    # Focus on nuclear receptor assays
    endocrine_assays = [
        "NR-AR", "NR-ER", "NR-ER-LBD",
        "NR-AR-LBD", "NR-Aromatase"
    ]

    input_data = {
        "query": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "assays": endocrine_assays,
        "include_inactive": False
    }

    print(f"\nScreening for endocrine disruption (Ibuprofen)")
    print(f"Testing {len(endocrine_assays)} nuclear receptor assays")

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        active = data['summary']['active_assays']

        if active == 0:
            print("\n✓ PASS: No endocrine disruption detected")
        else:
            print(f"\n⚠ WARNING: {active} endocrine assay(s) active")
            for assay in data['toxicity_results']:
                print(f"  - {assay['assay_name']}")
                if assay['ac50']:
                    print(f"    AC50: {assay['ac50']:.2f} uM")
    else:
        print(f"✗ Error: {result.error}")


async def example_4_batch_screening():
    """Example 4: Batch screening multiple compounds"""
    print("\n" + "="*60)
    print("Example 4: Batch Toxicity Screening")
    print("="*60)

    adapter = Tox21Adapter()

    compounds = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    print(f"\nScreening {len(compounds)} compounds...\n")

    results = []
    for name, smiles in compounds.items():
        print(f"  Querying {name}...", end=" ")
        result = await adapter.execute(smiles)

        if result.success:
            data = result.data
            results.append({
                'name': name,
                'smiles': smiles,
                'cid': data['compound']['cid'],
                'active_assays': data['summary']['active_assays'],
                'total_assays': data['summary']['total_assays'],
                'risk': data['summary']['overall_risk'],
                'flags': data['summary']['toxicity_flags']
            })
            print("✓")
        else:
            print(f"✗ ({result.error})")

        # Rate limiting
        await asyncio.sleep(0.3)

    # Display results table
    print("\n" + "-"*80)
    print(f"{'Compound':<15} {'CID':<10} {'Active/Total':<15} {'Risk':<10} {'Flags'}")
    print("-"*80)
    for r in results:
        flags = ', '.join(r['flags']) if r['flags'] else 'None'
        print(f"{r['name']:<15} {r['cid']:<10} {r['active_assays']}/{r['total_assays']:<13} {r['risk']:<10} {flags}")
    print("-"*80)


async def example_5_potency_filtering():
    """Example 5: Filter by potency (AC50 threshold)"""
    print("\n" + "="*60)
    print("Example 5: Potency-Based Filtering")
    print("="*60)

    adapter = Tox21Adapter()

    input_data = {
        "query": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "include_inactive": False,
        "ac50_threshold": 10.0  # Only show hits with AC50 < 10 uM
    }

    print(f"\nQuerying with AC50 threshold: {input_data['ac50_threshold']} uM")
    print("(Only showing potent hits)")

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\n✓ Found {len(data['toxicity_results'])} potent hits")

        if data['toxicity_results']:
            print("\n  Potent Toxicity Hits (AC50 < 10 uM):")
            for assay in data['toxicity_results']:
                print(f"    - {assay['assay_name']}")
                print(f"      AC50: {assay['ac50']:.2f} uM")
                print(f"      Endpoint: {assay['endpoint']}")
                print(f"      Confidence: {assay['confidence']}")
        else:
            print("\n  No potent hits found (all AC50 > 10 uM)")
    else:
        print(f"✗ Error: {result.error}")


async def example_6_assay_metadata():
    """Example 6: Explore assay information"""
    print("\n" + "="*60)
    print("Example 6: Assay Information and Metadata")
    print("="*60)

    adapter = Tox21Adapter()

    # Get adapter metadata
    metadata = adapter.get_metadata()
    print(f"\nAdapter Version: {metadata['version']}")
    print(f"Total Assays Available: {metadata['assays']['total']}")
    print(f"  - Nuclear Receptor: {metadata['assays']['nuclear_receptor']}")
    print(f"  - Stress Response: {metadata['assays']['stress_response']}")

    # List all assays by category
    print("\n" + "-"*60)
    print("Nuclear Receptor Assays:")
    print("-"*60)
    nr_assays = adapter.list_assays(category="nuclear_receptor")
    for assay_id in nr_assays:
        info = adapter.get_assay_info(assay_id)
        print(f"  {assay_id:<20} {info['name']}")
        print(f"  {'':20} Endpoint: {info['endpoint']}")
        print(f"  {'':20} PubChem AIDs: {info['aids']}")

    print("\n" + "-"*60)
    print("Stress Response Assays:")
    print("-"*60)
    sr_assays = adapter.list_assays(category="stress_response")
    for assay_id in sr_assays:
        info = adapter.get_assay_info(assay_id)
        print(f"  {assay_id:<20} {info['name']}")
        print(f"  {'':20} Endpoint: {info['endpoint']}")


async def example_7_comparative_analysis():
    """Example 7: Compare toxicity profiles"""
    print("\n" + "="*60)
    print("Example 7: Comparative Toxicity Analysis")
    print("="*60)

    adapter = Tox21Adapter()

    # Compare two similar compounds
    compounds = {
        "Compound A": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "Compound B": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
    }

    print("\nComparing toxicity profiles of two NSAIDs:\n")

    comparison = {}
    for name, smiles in compounds.items():
        print(f"  Analyzing {name}...", end=" ")
        result = await adapter.execute(smiles)

        if result.success:
            data = result.data
            comparison[name] = {
                'cid': data['compound']['cid'],
                'active': data['summary']['active_assays'],
                'total': data['summary']['total_assays'],
                'risk': data['summary']['overall_risk'],
                'flags': data['summary']['toxicity_flags'],
                'active_endpoints': set()
            }

            # Collect active endpoints
            for assay in data['toxicity_results']:
                if assay['activity'] == 'active':
                    comparison[name]['active_endpoints'].add(assay['endpoint'])

            print("✓")
        else:
            print(f"✗")

        await asyncio.sleep(0.3)

    # Display comparison
    print("\n" + "-"*80)
    print("Toxicity Profile Comparison:")
    print("-"*80)

    for name, data in comparison.items():
        print(f"\n{name} (CID: {data['cid']}):")
        print(f"  Active assays: {data['active']}/{data['total']}")
        print(f"  Overall risk: {data['risk']}")
        print(f"  Toxicity flags: {', '.join(data['flags']) if data['flags'] else 'None'}")

        if data['active_endpoints']:
            print(f"  Active endpoints: {', '.join(data['active_endpoints'])}")

    # Find common toxicity signals
    if len(comparison) == 2:
        names = list(comparison.keys())
        common = comparison[names[0]]['active_endpoints'] & comparison[names[1]]['active_endpoints']
        unique_a = comparison[names[0]]['active_endpoints'] - comparison[names[1]]['active_endpoints']
        unique_b = comparison[names[1]]['active_endpoints'] - comparison[names[0]]['active_endpoints']

        print("\n" + "-"*80)
        print("Comparison Summary:")
        print("-"*80)
        print(f"  Common toxicity signals: {', '.join(common) if common else 'None'}")
        print(f"  Unique to {names[0]}: {', '.join(unique_a) if unique_a else 'None'}")
        print(f"  Unique to {names[1]}: {', '.join(unique_b) if unique_b else 'None'}")


async def example_8_export_results():
    """Example 8: Export results to JSON"""
    print("\n" + "="*60)
    print("Example 8: Export Results to JSON")
    print("="*60)

    adapter = Tox21Adapter()

    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"\nQuerying and exporting results for: {smiles}")

    result = await adapter.execute(smiles)

    if result.success:
        # Export to JSON file
        output_file = "tox21_results.json"
        with open(output_file, 'w') as f:
            json.dump(result.data, f, indent=2)

        print(f"\n✓ Results exported to: {output_file}")
        print(f"  File size: {len(json.dumps(result.data))} bytes")
        print(f"  Contains {len(result.data['toxicity_results'])} assay results")

        # Also show metadata
        print(f"\n  Metadata:")
        for key, value in result.metadata.items():
            print(f"    {key}: {value}")
    else:
        print(f"✗ Error: {result.error}")


async def main():
    """Run all examples"""
    print("\n" + "="*80)
    print(" "*20 + "TOX21 ADAPTER USAGE EXAMPLES")
    print("="*80)

    examples = [
        example_1_basic_query,
        example_2_query_by_name,
        example_3_endocrine_screening,
        example_4_batch_screening,
        example_5_potency_filtering,
        example_6_assay_metadata,
        example_7_comparative_analysis,
        example_8_export_results
    ]

    for i, example in enumerate(examples, 1):
        try:
            await example()
        except Exception as e:
            print(f"\n✗ Error in example {i}: {e}")

        if i < len(examples):
            print("\n" + "."*80)
            await asyncio.sleep(1)  # Rate limiting between examples

    print("\n" + "="*80)
    print("All examples completed!")
    print("="*80 + "\n")


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
