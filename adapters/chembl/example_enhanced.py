"""
Enhanced ChEMBL Adapter - Example Usage
Demonstrates advanced bioactivity queries, target data, drug info, and mechanism of action
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../..")))

from adapters.chembl.adapter_enhanced import ChEMBLEnhancedAdapter


async def example_bioactivity_search():
    """
    Example 1: Search for bioactivity data with filtering
    Use case: Finding potent compounds for a target
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Bioactivity Search with Filtering")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # A known kinase inhibitor
    smiles = "c1ccc2c(c1)ncnc2N"

    print(f"\nSearching bioactivity data for: {smiles}")
    print("Filter: IC50 ≤ 100 nM")

    result = await adapter(
        smiles,
        mode="bioactivity",
        activity_type="IC50",
        max_value=100.0  # Only highly potent compounds
    )

    if result.success:
        data = result.data
        print(f"\n✓ Bioactivity Summary:")
        print(f"  ChEMBL ID: {data.get('chembl_id')}")
        print(f"  Total activities: {data.get('num_activities')}")
        print(f"  Unique targets: {data.get('num_targets')}")
        print(f"  Best IC50: {data.get('best_ic50_nm')} nM")
        print(f"  Best Ki: {data.get('best_ki_nm')} nM")
        print(f"  Best Kd: {data.get('best_kd_nm')} nM")

        if data.get('assay_types'):
            print(f"\n  Assay types: {', '.join(data['assay_types'][:5])}")

        if data.get('targets'):
            print(f"\n  Top targets: {', '.join(data['targets'][:5])}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_target_query():
    """
    Example 2: Get detailed target information
    Use case: Target validation and characterization
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Target Information Query")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Example: EGFR target
    target_id = "TARGET:CHEMBL203"

    print(f"\nQuerying target information for: {target_id}")
    print("Use case: Target validation for drug discovery")

    result = await adapter(target_id, mode="target")

    if result.success:
        data = result.data
        print(f"\n✓ Target Information:")
        print(f"  Target ID: {data.get('target_id')}")
        print(f"  Name: {data.get('pref_name')}")
        print(f"  Type: {data.get('target_type')}")
        print(f"  Organism: {data.get('organism')}")

        if data.get('target_components'):
            print(f"\n  Components: {len(data['target_components'])}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_drug_lookup():
    """
    Example 3: Check if compound is a drug or clinical candidate
    Use case: Drug repurposing and competitive intelligence
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Drug Status Lookup")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    print(f"\nChecking drug status for: {smiles} (Aspirin)")
    print("Use case: Drug repurposing opportunity assessment")

    result = await adapter(smiles, mode="drug")

    if result.success:
        data = result.data
        print(f"\n✓ Drug Status:")
        print(f"  ChEMBL ID: {data.get('chembl_id')}")
        print(f"  Is Drug: {'Yes' if data.get('is_drug') else 'No'}")

        if data.get('drug_info'):
            drug = data['drug_info']
            print(f"\n  Drug Information:")
            print(f"    Max Phase: {drug.get('max_phase', 'N/A')}")
            print(f"    First Approval: {drug.get('first_approval', 'N/A')}")
            if drug.get('indication_class'):
                print(f"    Indication: {drug['indication_class']}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_mechanism_search():
    """
    Example 4: Get mechanism of action data
    Use case: Understanding drug mechanisms for repurposing
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Mechanism of Action Data")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Imatinib (Gleevec) - known mechanism
    chembl_id = "CHEMBL941"

    print(f"\nRetrieving mechanism of action for: {chembl_id}")
    print("Use case: Mechanism-based drug repurposing")

    result = await adapter(chembl_id, mode="mechanism")

    if result.success:
        data = result.data
        print(f"\n✓ Mechanism of Action:")
        print(f"  ChEMBL ID: {data.get('chembl_id')}")
        print(f"  Number of mechanisms: {data.get('num_mechanisms')}")

        if data.get('mechanisms'):
            print("\n  Mechanisms:")
            for i, mech in enumerate(data['mechanisms'][:5], 1):
                print(f"\n    {i}. Target: {mech.get('target_chembl_id', 'N/A')}")
                print(f"       Mechanism: {mech.get('mechanism_of_action', 'N/A')}")
                print(f"       Action Type: {mech.get('action_type', 'N/A')}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_similarity_search():
    """
    Example 5: Similarity search in ChEMBL
    Use case: Finding similar bioactive compounds
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: ChEMBL Similarity Search")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Caffeine
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    print(f"\nSearching for similar bioactive compounds")
    print("Use case: Finding bioactive analogs for SAR")

    result = await adapter(
        smiles,
        mode="similarity",
        similarity=80  # 80% similarity
    )

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} similar compounds")
        print(f"  Similarity threshold: {data['similarity_threshold']}%")

        if data['molecules']:
            print("\n  Top 5 Similar Compounds:")
            for i, mol in enumerate(data['molecules'][:5], 1):
                print(f"\n    {i}. {mol.get('molecule_chembl_id', 'N/A')}")
                print(f"       Similarity: {mol.get('similarity', 'N/A')}")
                if mol.get('pref_name'):
                    print(f"       Name: {mol['pref_name']}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_substructure_search():
    """
    Example 6: Substructure search for scaffold-based discovery
    Use case: Finding all compounds with a specific scaffold
    """
    print("\n" + "="*80)
    print("EXAMPLE 6: Substructure Search")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Quinoline scaffold
    substructure = "c1ccc2ncccc2c1"

    print(f"\nSearching for quinoline-containing bioactive compounds")
    print("Use case: Scaffold-based drug discovery")

    result = await adapter(substructure, mode="substructure")

    if result.success:
        data = result.data
        print(f"\n✓ Found {data['num_results']} compounds with quinoline scaffold")

        if data['molecules']:
            print("\n  Sample Compounds:")
            for i, mol in enumerate(data['molecules'][:5], 1):
                print(f"    {i}. {mol.get('molecule_chembl_id', 'N/A')}")
                if mol.get('pref_name'):
                    print(f"       Name: {mol['pref_name']}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_advanced_filtering():
    """
    Example 7: Advanced activity filtering
    Use case: Finding highly potent and selective compounds
    """
    print("\n" + "="*80)
    print("EXAMPLE 7: Advanced Activity Filtering")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Known EGFR inhibitor scaffold
    smiles = "c1cc2c(cc1N)ncnc2N"

    print(f"\nSearching for highly potent compounds (IC50 < 10 nM)")
    print("Use case: Lead optimization - finding potent hits")

    result = await adapter(
        smiles,
        mode="bioactivity",
        activity_type="IC50",
        max_value=10.0,  # Very potent
        target_type="PROTEIN"
    )

    if result.success:
        data = result.data
        print(f"\n✓ High-Potency Bioactivity Results:")
        print(f"  Total activities: {data.get('num_activities')}")
        print(f"  Unique targets: {data.get('num_targets')}")
        print(f"  Best IC50: {data.get('best_ic50_nm')} nM")

        if data.get('bioactivities'):
            print(f"\n  Sample Activities (IC50 < 10 nM):")
            for i, act in enumerate(data['bioactivities'][:5], 1):
                print(f"\n    {i}. Target: {act.get('target_chembl_id', 'N/A')}")
                print(f"       IC50: {act.get('standard_value', 'N/A')} {act.get('standard_units', 'N/A')}")
                print(f"       Assay Type: {act.get('assay_type', 'N/A')}")
    else:
        print(f"\n✗ Error: {result.error}")


async def example_comprehensive_profile():
    """
    Example 8: Build comprehensive compound profile from ChEMBL
    Use case: Complete bioactivity characterization
    """
    print("\n" + "="*80)
    print("EXAMPLE 8: Comprehensive Bioactivity Profile")
    print("="*80)

    adapter = ChEMBLEnhancedAdapter()

    # Imatinib
    chembl_id = "CHEMBL941"

    print(f"\nBuilding comprehensive profile for: {chembl_id} (Imatinib)")
    print("Use case: Complete bioactivity characterization for drug repurposing")

    # Get multiple types of data
    tasks = [
        adapter(chembl_id, mode="bioactivity"),
        adapter(chembl_id, mode="drug"),
        adapter(chembl_id, mode="mechanism")
    ]

    results = await asyncio.gather(*tasks)

    print("\n" + "="*60)
    print("COMPREHENSIVE BIOACTIVITY PROFILE")
    print("="*60)

    # Bioactivity
    if results[0].success:
        data = results[0].data
        print(f"\n1. BIOACTIVITY DATA")
        print(f"   Total activities: {data.get('num_activities')}")
        print(f"   Unique targets: {data.get('num_targets')}")
        print(f"   Best IC50: {data.get('best_ic50_nm')} nM")
        print(f"   Best Ki: {data.get('best_ki_nm')} nM")

    # Drug status
    if results[1].success:
        data = results[1].data
        print(f"\n2. DRUG STATUS")
        print(f"   Is Drug: {'Yes' if data.get('is_drug') else 'No'}")
        if data.get('drug_info'):
            print(f"   Max Phase: {data['drug_info'].get('max_phase', 'N/A')}")

    # Mechanism
    if results[2].success:
        data = results[2].data
        print(f"\n3. MECHANISM OF ACTION")
        print(f"   Known mechanisms: {data.get('num_mechanisms', 0)}")
        if data.get('mechanisms'):
            for mech in data['mechanisms'][:3]:
                print(f"   • {mech.get('mechanism_of_action', 'N/A')}")

    print("\n" + "="*60)


async def main():
    """Run all examples"""
    print("\n" + "#"*80)
    print("# Enhanced ChEMBL Adapter - Example Usage")
    print("# Advanced bioactivity queries and target analysis")
    print("#"*80)

    await example_bioactivity_search()
    await example_target_query()
    await example_drug_lookup()
    await example_mechanism_search()
    await example_similarity_search()
    await example_substructure_search()
    await example_advanced_filtering()
    await example_comprehensive_profile()

    print("\n" + "#"*80)
    print("# Examples Complete!")
    print("#"*80)
    print("\nKey Features Demonstrated:")
    print("  ✓ Advanced bioactivity filtering")
    print("  ✓ Target information queries")
    print("  ✓ Drug status and clinical data")
    print("  ✓ Mechanism of action")
    print("  ✓ Similarity and substructure search")
    print("  ✓ Multi-dimensional filtering")
    print("\nUse Cases:")
    print("  • Target validation")
    print("  • Lead optimization")
    print("  • Drug repurposing")
    print("  • SAR analysis")
    print("  • Competitive intelligence")


if __name__ == "__main__":
    asyncio.run(main())
