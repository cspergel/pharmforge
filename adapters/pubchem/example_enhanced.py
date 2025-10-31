"""
Enhanced PubChem Adapter - Example Usage
Demonstrates advanced features: similarity search, bioassays, vendors, and patents
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from adapters.pubchem.adapter_enhanced import PubChemEnhancedAdapter


async def example_basic_properties():
    """
    Example 1: Basic property lookup (enhanced with more properties)
    Use case: Compound characterization
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Enhanced Property Lookup")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # Caffeine
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    print(f"\nLookup properties for: {smiles} (Caffeine)")

    result = await adapter(smiles, mode="properties")

    if result.success:
        data = result.data
        print(f"\n✓ Compound Properties:")
        print(f"  PubChem CID: {data.get('cid')}")
        print(f"  IUPAC Name: {data.get('iupac_name', 'N/A')[:60]}...")
        print(f"  Molecular Weight: {data.get('molecular_weight')} Da")
        print(f"  LogP: {data.get('logp')}")
        print(f"  TPSA: {data.get('tpsa')} Ų")
        print(f"  H-Bond Donors: {data.get('h_bond_donors')}")
        print(f"  H-Bond Acceptors: {data.get('h_bond_acceptors')}")
        print(f"  Rotatable Bonds: {data.get('rotatable_bonds')}")
        print(f"  Complexity: {data.get('complexity')}")
        print(f"\n  InChIKey: {data.get('inchikey')}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_similarity_search():
    """
    Example 2: Find similar compounds in PubChem
    Use case: Analog searching and SAR exploration
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Similarity Search")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # Ibuprofen
    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

    print(f"\nSearching for compounds similar to: {smiles} (Ibuprofen)")
    print("Use case: Finding structural analogs for SAR analysis")

    result = await adapter(
        smiles,
        mode="similarity",
        similarity_threshold=85,  # 85% similarity
        max_results=50
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} similar compounds")
        print(f"  Similarity threshold: {data['similarity_threshold']}%")
        print(f"\n  Top 10 PubChem CIDs:")
        for i, cid in enumerate(data['cids'][:10], 1):
            print(f"    {i:2d}. CID {cid}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_substructure_search():
    """
    Example 3: Substructure search
    Use case: Scaffold-based searching
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Substructure Search")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # Indole scaffold - common in drugs
    substructure = "c1ccc2c(c1)cc[nH]2"

    print(f"\nSearching for compounds containing indole scaffold")
    print("Use case: Finding indole-based drug candidates")

    result = await adapter(
        substructure,
        mode="substructure",
        max_results=30
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} compounds with indole substructure")
        print(f"\n  Sample CIDs:")
        for i, cid in enumerate(data['cids'][:10], 1):
            print(f"    {i}. CID {cid}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_bioassay_search():
    """
    Example 4: Get bioactivity data from PubChem BioAssay
    Use case: Finding screening data for compounds
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Bioassay Data Retrieval")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    print(f"\nRetrieving bioassay data for aspirin")
    print("Use case: Finding experimental activity data")

    result = await adapter(smiles, mode="bioassays")

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_assays']} bioassay records")
        print(f"  PubChem CID: {data['cid']}")

        if data['bioassays']:
            print("\n  Sample Assays:")
            for i, assay in enumerate(data['bioassays'][:5], 1):
                print(f"\n    {i}. AID: {assay.get('aid', 'N/A')}")
                print(f"       Activity: {assay.get('activity_outcome', 'N/A')}")
                if assay.get('activity_value'):
                    print(f"       Value: {assay['activity_value']}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_vendor_lookup():
    """
    Example 5: Get vendor and purchasability information
    Use case: Sourcing compounds for testing
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: Vendor and Purchasability Lookup")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # A common research compound
    smiles = "c1ccccc1"  # Benzene

    print(f"\nLooking up vendors for: {smiles}")
    print("Use case: Finding commercial sources for compound procurement")

    result = await adapter(smiles, mode="vendors")

    if result.success:
        data = result.data
        print(f"\n✓ Vendor Information:")
        print(f"  PubChem CID: {data['cid']}")
        print(f"  Number of vendors: {data['num_vendors']}")
        print(f"  Purchasable: {'Yes' if data['purchasable'] else 'No'}")

        if data['vendors']:
            print("\n  Commercial Sources:")
            for i, vendor in enumerate(data['vendors'][:10], 1):
                print(f"    {i}. {vendor}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_patent_search():
    """
    Example 6: Find patents mentioning a compound
    Use case: IP analysis and literature review
    """
    print("\n" + "="*80)
    print("EXAMPLE 6: Patent Mentions")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # Atorvastatin (Lipitor)
    smiles = "CC(C)c1c(c(c(c(c1)F)c2ccccc2)c3ccc(cc3)O)C(=O)Nc4ccccc4C(=O)O"

    print(f"\nSearching for patent mentions of atorvastatin")
    print("Use case: Patent landscape analysis")

    result = await adapter(smiles, mode="patents")

    if result.success:
        data = result.data
        print(f"\n✓ Patent Information:")
        print(f"  PubChem CID: {data['cid']}")
        print(f"  Number of patents: {data['num_patents']}")

        if data['patents']:
            print("\n  Sample Patents:")
            for i, patent_id in enumerate(data['patents'][:10], 1):
                print(f"    {i}. {patent_id}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_comprehensive_profile():
    """
    Example 7: Build comprehensive compound profile
    Use case: Complete compound characterization for drug discovery
    """
    print("\n" + "="*80)
    print("EXAMPLE 7: Comprehensive Compound Profile")
    print("="*80)

    adapter = PubChemEnhancedAdapter()

    # Metformin
    smiles = "CN(C)C(=N)NC(=N)N"

    print(f"\nBuilding comprehensive profile for: {smiles} (Metformin)")
    print("Use case: Complete characterization for drug discovery decision-making")

    # Get multiple types of data
    tasks = [
        adapter(smiles, mode="properties"),
        adapter(smiles, mode="bioassays"),
        adapter(smiles, mode="vendors"),
        adapter(smiles, mode="patents")
    ]

    results = await asyncio.gather(*tasks)

    print("\n" + "="*60)
    print("COMPREHENSIVE COMPOUND PROFILE")
    print("="*60)

    # Properties
    if results[0].success:
        data = results[0].data
        print(f"\n1. MOLECULAR PROPERTIES")
        print(f"   CID: {data.get('cid')}")
        print(f"   MW: {data.get('molecular_weight')} Da")
        print(f"   LogP: {data.get('logp')}")
        print(f"   TPSA: {data.get('tpsa')} Ų")

    # Bioassays
    if results[1].success:
        data = results[1].data
        print(f"\n2. BIOACTIVITY DATA")
        print(f"   Bioassays: {data.get('num_assays', 0)} records")

    # Vendors
    if results[2].success:
        data = results[2].data
        print(f"\n3. COMMERCIAL AVAILABILITY")
        print(f"   Vendors: {data.get('num_vendors', 0)}")
        print(f"   Purchasable: {'Yes' if data.get('purchasable') else 'No'}")

    # Patents
    if results[3].success:
        data = results[3].data
        print(f"\n4. PATENT LANDSCAPE")
        print(f"   Patents: {data.get('num_patents', 0)} mentions")

    print("\n" + "="*60)


async def main():
    """Run all examples"""
    print("\n" + "#"*80)
    print("# Enhanced PubChem Adapter - Example Usage")
    print("# Advanced compound search and characterization")
    print("#"*80)

    await example_basic_properties()
    await example_similarity_search()
    await example_substructure_search()
    await example_bioassay_search()
    await example_vendor_lookup()
    await example_patent_search()
    await example_comprehensive_profile()

    print("\n" + "#"*80)
    print("# Examples Complete!")
    print("#"*80)
    print("\nKey Features Demonstrated:")
    print("  ✓ Enhanced property lookup")
    print("  ✓ Similarity and substructure search")
    print("  ✓ Bioassay data retrieval")
    print("  ✓ Vendor and purchasability information")
    print("  ✓ Patent mentions")
    print("  ✓ Comprehensive compound profiling")
    print("\nUse Cases:")
    print("  • Compound characterization")
    print("  • Analog searching and SAR")
    print("  • Bioactivity data mining")
    print("  • Commercial sourcing")
    print("  • IP analysis")


if __name__ == "__main__":
    asyncio.run(main())
