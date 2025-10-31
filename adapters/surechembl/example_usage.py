"""
SureChEMBL Adapter - Example Usage
Demonstrates patent chemistry search, compound extraction, and IP landscape analysis
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from adapters.surechembl.adapter import SureChEMBLAdapter


async def example_structure_search():
    """
    Example 1: Search patents for compounds similar to a query structure
    Use case: IP landscape analysis and freedom-to-operate
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Patent Structure Search")
    print("="*80)

    adapter = SureChEMBLAdapter()

    # Aspirin-like structure
    query_smiles = "CC(=O)Oc1ccccc1C(=O)O"

    print(f"\nSearching patents for structures similar to: {query_smiles}")
    print("Use case: Freedom-to-operate analysis for aspirin analogs")

    result = await adapter(
        query_smiles,
        mode="structure_search",
        search_type="similarity",
        similarity_threshold=0.85,
        max_results=20
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_compounds']} compound mentions in {data['num_patents']} patents")
        print(f"  Patent families: {data['num_patent_families']}")

        if data.get('date_range'):
            print(f"  Date range: {data['date_range']['earliest']} to {data['date_range']['latest']}")

        if data.get('top_applicants'):
            print("\nTop Patent Applicants:")
            for i, app in enumerate(data['top_applicants'][:5], 1):
                print(f"  {i}. {app['name']}: {app['count']} patents")

        if data['patents']:
            print("\nSample Patents:")
            for i, patent in enumerate(data['patents'][:3], 1):
                print(f"\n  {i}. Patent: {patent.get('patent_id', 'N/A')}")
                print(f"     Publication: {patent.get('publication_date', 'N/A')}")
                print(f"     Applicant: {patent.get('applicant', 'N/A')}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_substructure_search():
    """
    Example 2: Find all patents containing a specific substructure
    Use case: Prior art search and novelty assessment
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Patent Substructure Search")
    print("="*80)

    adapter = SureChEMBLAdapter()

    # Thiazole ring - common in antibiotics
    substructure = "c1scnc1"

    print(f"\nSearching for thiazole-containing compounds in patents")
    print("Use case: Prior art search for thiazole-based antibiotics")

    result = await adapter(
        substructure,
        mode="structure_search",
        search_type="substructure",
        max_results=30
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_compounds']} thiazole-containing compounds")
        print(f"  Across {data['num_patents']} patents")

        if data.get('top_applicants'):
            print("\nKey Players in Thiazole Chemistry:")
            for app in data['top_applicants'][:5]:
                print(f"  • {app['name']}: {app['count']} patents")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_patent_lookup():
    """
    Example 3: Get detailed information about a specific patent
    Use case: Patent analysis and compound extraction
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Patent Details Lookup")
    print("="*80)

    adapter = SureChEMBLAdapter()

    # Example patent ID (format varies by jurisdiction)
    patent_id = "US20190123456"  # Example format

    print(f"\nRetrieving details for patent: {patent_id}")
    print("Use case: Detailed patent analysis and competitive intelligence")

    result = await adapter(
        patent_id,
        mode="patent_lookup"
    )

    if result.success:
        data = result.data
        print(f"\n✓ Patent Information:")
        print(f"  Title: {data.get('title', 'N/A')}")
        print(f"  Publication Date: {data.get('publication_date', 'N/A')}")
        print(f"  Number of Compounds: {data.get('num_compounds', 0)}")

        if data.get('applicants'):
            print(f"\n  Applicants:")
            for applicant in data['applicants'][:3]:
                print(f"    • {applicant}")

        if data.get('inventors'):
            print(f"\n  Inventors:")
            for inventor in data['inventors'][:3]:
                print(f"    • {inventor}")

        print(f"\n  Patent Family: {data.get('patent_family', 'N/A')}")

        if data.get('abstract'):
            print(f"\n  Abstract: {data['abstract'][:200]}...")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_extract_compounds():
    """
    Example 4: Extract all chemical structures from a patent
    Use case: Competitive intelligence and compound library building
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Extract Compounds from Patent")
    print("="*80)

    adapter = SureChEMBLAdapter()

    patent_id = "US20190123456"  # Example

    print(f"\nExtracting all compounds from patent: {patent_id}")
    print("Use case: Building a competitive intelligence compound library")

    result = await adapter(
        patent_id,
        mode="extract_compounds",
        max_results=50
    )

    if result.success:
        data = result.data
        print(f"\n✓ Extracted {data['num_compounds']} compounds from patent")

        if data['compounds']:
            print("\nSample Extracted Compounds:")
            for i, compound in enumerate(data['compounds'][:5], 1):
                print(f"\n  {i}. {compound.get('surechembl_id', 'N/A')}")
                print(f"     SMILES: {compound.get('smiles', 'N/A')[:60]}...")
                if compound.get('mol_weight'):
                    print(f"     MW: {compound['mol_weight']} Da")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_ip_landscape():
    """
    Example 5: Comprehensive IP landscape analysis
    Use case: Strategic patent analysis for drug discovery program
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: IP Landscape Analysis")
    print("="*80)

    adapter = SureChEMBLAdapter()

    # Novel kinase inhibitor scaffold
    scaffold_smiles = "c1ccc2c(c1)ncnc2N"

    print(f"\nAnalyzing IP landscape for kinase inhibitor scaffold")
    print("Use case: Freedom-to-operate assessment for new drug candidate")

    result = await adapter(
        scaffold_smiles,
        mode="structure_search",
        search_type="similarity",
        similarity_threshold=0.9,
        max_results=50
    )

    if result.success:
        data = result.data
        print(f"\n✓ IP Landscape Summary:")
        print(f"  Total compounds found: {data['num_compounds']}")
        print(f"  Unique patents: {data['num_patents']}")
        print(f"  Patent families: {data['num_patent_families']}")

        if data.get('date_range'):
            dr = data['date_range']
            print(f"\n  Patent Activity Timeline:")
            print(f"    First patent: {dr['earliest']}")
            print(f"    Most recent: {dr['latest']}")

        if data.get('top_applicants'):
            print(f"\n  Competitive Landscape (Top Applicants):")
            for i, app in enumerate(data['top_applicants'][:10], 1):
                print(f"    {i:2d}. {app['name']}: {app['count']} patents")

        # Analyze patent concentration
        if data.get('top_applicants'):
            total_patents = data['num_patents']
            top_3_count = sum(app['count'] for app in data['top_applicants'][:3])
            concentration = (top_3_count / total_patents * 100) if total_patents > 0 else 0

            print(f"\n  Market Concentration:")
            print(f"    Top 3 applicants hold {concentration:.1f}% of patents")

            if concentration > 70:
                print(f"    ⚠️  High concentration - crowded IP space")
            elif concentration < 30:
                print(f"    ✓ Low concentration - more freedom to operate")

    else:
        print(f"\n✗ Error: {result.error}")


async def main():
    """Run all examples"""
    print("\n" + "#"*80)
    print("# SureChEMBL Adapter - Example Usage")
    print("# Patent chemistry search and IP landscape analysis")
    print("#"*80)

    await example_structure_search()
    await example_substructure_search()
    await example_patent_lookup()
    await example_extract_compounds()
    await example_ip_landscape()

    print("\n" + "#"*80)
    print("# Examples Complete!")
    print("#"*80)
    print("\nKey Features Demonstrated:")
    print("  ✓ Patent structure search (similarity & substructure)")
    print("  ✓ Patent details and metadata retrieval")
    print("  ✓ Compound extraction from patents")
    print("  ✓ IP landscape analysis")
    print("\nUse Cases:")
    print("  • Freedom-to-operate analysis")
    print("  • Prior art searches")
    print("  • Competitive intelligence")
    print("  • Patent portfolio analysis")
    print("  • Strategic IP planning")


if __name__ == "__main__":
    asyncio.run(main())
