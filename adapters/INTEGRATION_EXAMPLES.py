"""
Database Adapters - Retrosynthesis Integration Examples
Demonstrates how to integrate chemical database adapters with retrosynthesis workflows
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter
from adapters.pubchem.adapter_enhanced import PubChemEnhancedAdapter
from adapters.chembl.adapter_enhanced import ChEMBLEnhancedAdapter
from adapters.surechembl.adapter import SureChEMBLAdapter


async def workflow_1_building_block_validation():
    """
    WORKFLOW 1: Validate Building Block Availability

    Use case: Before running retrosynthesis, check if proposed building blocks
    are commercially available using ZINC and PubChem
    """
    print("\n" + "="*80)
    print("WORKFLOW 1: Building Block Availability Validation")
    print("="*80)

    zinc = ZINCFragmentsAdapter()
    pubchem = PubChemEnhancedAdapter()

    # Proposed building blocks from retrosynthesis
    building_blocks = [
        "c1ccccc1Br",  # Bromobenzene
        "CC(=O)Cl",    # Acetyl chloride
        "c1ccc2[nH]ccc2c1"  # Indole
    ]

    print("\nValidating commercial availability of building blocks...")
    print("Use case: Pre-synthesis validation for retrosynthesis routes")

    for i, smiles in enumerate(building_blocks, 1):
        print(f"\n{i}. Building Block: {smiles}")

        # Check ZINC for fragments
        zinc_result = await zinc(
            smiles,
            search_type="exact",
            get_purchasability=True
        )

        # Check PubChem for vendors
        pubchem_result = await pubchem(smiles, mode="vendors")

        # Aggregate results
        zinc_available = False
        pubchem_vendors = 0

        if zinc_result.success and zinc_result.data['fragments']:
            purch = zinc_result.data['fragments'][0].get('purchasability', {})
            zinc_available = purch.get('purchasable', False)
            zinc_vendors = purch.get('vendor_count', 0)

        if pubchem_result.success:
            pubchem_vendors = pubchem_result.data.get('num_vendors', 0)

        # Decision
        total_sources = zinc_vendors + pubchem_vendors if zinc_available else pubchem_vendors

        print(f"   ZINC: {'Available' if zinc_available else 'Not found'}")
        if zinc_available:
            print(f"         {zinc_vendors} vendors")
        print(f"   PubChem: {pubchem_vendors} vendors")
        print(f"   Status: {'✓ COMMERCIALLY AVAILABLE' if total_sources > 0 else '✗ NOT AVAILABLE - NEEDS SYNTHESIS'}")


async def workflow_2_fragment_based_retro():
    """
    WORKFLOW 2: Fragment-Based Retrosynthesis Planning

    Use case: Use ZINC fragments to suggest disconnections for retrosynthesis
    """
    print("\n" + "="*80)
    print("WORKFLOW 2: Fragment-Based Retrosynthesis Planning")
    print("="*80)

    zinc = ZINCFragmentsAdapter()
    chembl = ChEMBLEnhancedAdapter()

    # Target molecule
    target_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen

    print(f"\nTarget molecule: {target_smiles} (Ibuprofen)")
    print("Use case: Fragment-based disconnection analysis")

    # Step 1: Find similar fragments in ZINC
    print("\nStep 1: Finding similar fragments...")
    zinc_result = await zinc(
        target_smiles,
        search_type="similarity",
        similarity_threshold=0.6,
        max_results=20,
        mw_max=200  # Smaller than target
    )

    if zinc_result.success:
        fragments = zinc_result.data['fragments']
        print(f"✓ Found {len(fragments)} potential building blocks")

        # Step 2: Check bioactivity of fragments (for privileged scaffolds)
        print("\nStep 2: Checking fragment bioactivity (privileged scaffolds)...")

        bioactive_fragments = []
        for frag in fragments[:5]:  # Check top 5
            frag_smiles = frag.get('smiles')
            if frag_smiles:
                chembl_result = await chembl(
                    frag_smiles,
                    mode="bioactivity"
                )

                if chembl_result.success:
                    num_activities = chembl_result.data.get('num_activities', 0)
                    if num_activities > 0:
                        bioactive_fragments.append({
                            'smiles': frag_smiles,
                            'zinc_id': frag.get('zinc_id'),
                            'activities': num_activities
                        })

        print(f"✓ Found {len(bioactive_fragments)} bioactive fragments")
        print("\nPrivileged Scaffolds (bioactive fragments):")
        for i, frag in enumerate(bioactive_fragments, 1):
            print(f"  {i}. {frag['zinc_id']}")
            print(f"     Activities: {frag['activities']}")
            print(f"     SMILES: {frag['smiles'][:50]}...")

        print("\nRetrosynthesis Strategy:")
        if bioactive_fragments:
            print("  → Use privileged scaffolds as starting points")
            print("  → Focus disconnections around bioactive cores")
        else:
            print("  → No privileged scaffolds found")
            print("  → Use standard retrosynthesis approach")


async def workflow_3_ip_aware_retrosynthesis():
    """
    WORKFLOW 3: IP-Aware Retrosynthesis Planning

    Use case: Check patent landscape before committing to synthesis route
    """
    print("\n" + "="*80)
    print("WORKFLOW 3: IP-Aware Retrosynthesis Planning")
    print("="*80)

    surechembl = SureChEMBLAdapter()

    # Proposed synthetic intermediates
    intermediates = [
        ("c1ccc(cc1)C(=O)Cl", "Key intermediate 1"),
        ("CC(C)Cc1ccccc1", "Key intermediate 2")
    ]

    print("\nChecking patent landscape for synthetic intermediates...")
    print("Use case: Freedom-to-operate analysis for synthesis routes")

    for smiles, name in intermediates:
        print(f"\n{name}: {smiles}")

        result = await surechembl(
            smiles,
            mode="structure_search",
            search_type="exact",
            max_results=10
        )

        if result.success:
            num_patents = result.data.get('num_patents', 0)

            if num_patents == 0:
                print(f"   ✓ NO EXACT MATCHES - Clear for use")
            elif num_patents < 5:
                print(f"   ⚠️  {num_patents} patents found - Review recommended")
            else:
                print(f"   ✗ {num_patents} patents found - HIGH RISK")
                print(f"   → Consider alternative route")

                # Check applicants
                if result.data.get('top_applicants'):
                    print(f"   Key patent holders:")
                    for app in result.data['top_applicants'][:3]:
                        print(f"     • {app['name']}")


async def workflow_4_bioactivity_guided_synthesis():
    """
    WORKFLOW 4: Bioactivity-Guided Synthesis Planning

    Use case: Prioritize synthesis routes based on bioactivity of intermediates
    """
    print("\n" + "="*80)
    print("WORKFLOW 4: Bioactivity-Guided Synthesis Planning")
    print("="*80)

    chembl = ChEMBLEnhancedAdapter()

    # Two alternative synthetic routes with different intermediates
    routes = {
        "Route A": ["c1ccc2c(c1)ncnc2N", "c1ccc(cc1)NC(=O)C"],
        "Route B": ["c1ccccc1N", "CC(=O)Oc1ccccc1"]
    }

    print("\nComparing alternative synthesis routes by intermediate bioactivity...")
    print("Use case: Select route with most drug-like intermediates")

    route_scores = {}

    for route_name, intermediates in routes.items():
        print(f"\n{route_name}:")
        total_activities = 0
        total_targets = 0

        for i, smiles in enumerate(intermediates, 1):
            result = await chembl(
                smiles,
                mode="bioactivity"
            )

            if result.success:
                num_activities = result.data.get('num_activities', 0)
                num_targets = result.data.get('num_targets', 0)
                best_ic50 = result.data.get('best_ic50_nm')

                total_activities += num_activities
                total_targets += num_targets

                print(f"  Intermediate {i}: {smiles[:30]}...")
                print(f"    Activities: {num_activities}, Targets: {num_targets}")
                if best_ic50:
                    print(f"    Best IC50: {best_ic50} nM")

        route_scores[route_name] = {
            'activities': total_activities,
            'targets': total_targets
        }

        print(f"  Route Score: {total_activities} activities, {total_targets} targets")

    # Recommendation
    print("\nRecommendation:")
    best_route = max(route_scores.items(), key=lambda x: x[1]['activities'])
    print(f"  → {best_route[0]} has more bioactive intermediates")
    print(f"  → Intermediates are more 'drug-like' and validated")
    print(f"  → Higher confidence in synthesis success")


async def workflow_5_purchasability_optimization():
    """
    WORKFLOW 5: Optimize Retrosynthesis for Purchasable Building Blocks

    Use case: Modify retrosynthesis to use only commercially available materials
    """
    print("\n" + "="*80)
    print("WORKFLOW 5: Purchasability-Optimized Retrosynthesis")
    print("="*80)

    zinc = ZINCFragmentsAdapter()
    pubchem = PubChemEnhancedAdapter()

    target = "CC(C)(C)c1ccc(cc1)O"  # 4-tert-Butylphenol

    print(f"\nTarget: {target}")
    print("Goal: Find purchasable fragments for rapid synthesis")

    # Find similar purchasable fragments
    print("\nSearching for purchasable alternatives...")

    zinc_result = await zinc(
        target,
        search_type="similarity",
        similarity_threshold=0.7,
        get_purchasability=True,
        max_results=20
    )

    if zinc_result.success:
        fragments = zinc_result.data['fragments']

        # Filter for purchasable only
        purchasable = [
            f for f in fragments
            if f.get('purchasability', {}).get('purchasable', False)
        ]

        print(f"\n✓ Found {len(purchasable)} purchasable alternatives")

        # Rank by delivery time
        print("\nTop Purchasable Building Blocks (by availability):")

        for i, frag in enumerate(purchasable[:5], 1):
            purch = frag.get('purchasability', {})
            print(f"\n  {i}. {frag.get('zinc_id')}")
            print(f"     MW: {frag.get('mw')} Da")
            print(f"     Vendors: {purch.get('vendor_count', 0)}")
            print(f"     Delivery: {purch.get('delivery_days', 'N/A')} days")

            if purch.get('suppliers'):
                print(f"     Available from: {', '.join(purch['suppliers'][:3])}")

        print("\nSynthesis Strategy:")
        print("  → Use most available building block")
        print("  → Minimize custom synthesis steps")
        print("  → Reduce time to first compound")


async def workflow_6_comprehensive_route_evaluation():
    """
    WORKFLOW 6: Comprehensive Synthesis Route Evaluation

    Use case: Score retrosynthesis routes by multiple criteria
    """
    print("\n" + "="*80)
    print("WORKFLOW 6: Comprehensive Route Evaluation")
    print("="*80)

    zinc = ZINCFragmentsAdapter()
    pubchem = PubChemEnhancedAdapter()
    chembl = ChEMBLEnhancedAdapter()
    surechembl = SureChEMBLAdapter()

    # Proposed route intermediates
    route_intermediates = [
        "c1ccccc1Br",
        "CC(=O)Cl"
    ]

    print("\nEvaluating synthesis route comprehensively...")
    print("Criteria: Availability, Cost, Bioactivity, Patent Risk")

    total_score = 0
    max_score = len(route_intermediates) * 4  # 4 criteria

    for i, smiles in enumerate(route_intermediates, 1):
        print(f"\nIntermediate {i}: {smiles}")
        step_score = 0

        # 1. Purchasability (ZINC + PubChem)
        zinc_result = await zinc(smiles, search_type="exact", get_purchasability=True)
        pubchem_result = await pubchem(smiles, mode="vendors")

        purchasable = False
        if zinc_result.success and zinc_result.data['fragments']:
            purchasable = zinc_result.data['fragments'][0].get('purchasability', {}).get('purchasable', False)

        vendors = pubchem_result.data.get('num_vendors', 0) if pubchem_result.success else 0

        if purchasable or vendors > 0:
            step_score += 1
            print(f"  ✓ Availability: {vendors} vendors")
        else:
            print(f"  ✗ Availability: Not commercially available")

        # 2. Patent risk (SureChEMBL)
        patent_result = await surechembl(smiles, mode="structure_search", search_type="exact")

        patent_risk = patent_result.data.get('num_patents', 0) if patent_result.success else 0

        if patent_risk < 5:
            step_score += 1
            print(f"  ✓ Patent Risk: Low ({patent_risk} patents)")
        else:
            print(f"  ⚠️  Patent Risk: High ({patent_risk} patents)")

        # 3. Bioactivity (ChEMBL) - want low activity for intermediates
        chembl_result = await chembl(smiles, mode="bioactivity")

        activities = chembl_result.data.get('num_activities', 0) if chembl_result.success else 0

        # For intermediates, we actually want LOW activity (less toxic)
        if activities < 10:
            step_score += 1
            print(f"  ✓ Safety: Low bioactivity ({activities} records)")
        else:
            print(f"  ⚠️  Safety: Significant bioactivity ({activities} records)")

        # 4. Structural complexity (simpler is better)
        # Simple heuristic: fewer rings and heavy atoms
        if len(smiles) < 20:  # Simple structure
            step_score += 1
            print(f"  ✓ Complexity: Simple structure")
        else:
            print(f"  ⚠️  Complexity: Complex structure")

        total_score += step_score
        print(f"  Step Score: {step_score}/4")

    # Final evaluation
    route_score = (total_score / max_score) * 100

    print("\n" + "="*60)
    print("ROUTE EVALUATION SUMMARY")
    print("="*60)
    print(f"Total Score: {total_score}/{max_score} ({route_score:.1f}%)")

    if route_score >= 75:
        print("✓ EXCELLENT ROUTE - Proceed with confidence")
    elif route_score >= 50:
        print("⚠️  GOOD ROUTE - Minor optimizations recommended")
    else:
        print("✗ POOR ROUTE - Consider alternatives")


async def main():
    """Run all integration examples"""
    print("\n" + "#"*80)
    print("# Chemical Database Adapters - Retrosynthesis Integration")
    print("# Comprehensive workflows for synthesis planning")
    print("#"*80)

    await workflow_1_building_block_validation()
    await workflow_2_fragment_based_retro()
    await workflow_3_ip_aware_retrosynthesis()
    await workflow_4_bioactivity_guided_synthesis()
    await workflow_5_purchasability_optimization()
    await workflow_6_comprehensive_route_evaluation()

    print("\n" + "#"*80)
    print("# Integration Examples Complete!")
    print("#"*80)
    print("\nWorkflows Demonstrated:")
    print("  1. Building Block Availability Validation")
    print("  2. Fragment-Based Retrosynthesis Planning")
    print("  3. IP-Aware Retrosynthesis Planning")
    print("  4. Bioactivity-Guided Synthesis Planning")
    print("  5. Purchasability-Optimized Retrosynthesis")
    print("  6. Comprehensive Route Evaluation")
    print("\nIntegrated Adapters:")
    print("  • ZINC Fragments - Building block sourcing")
    print("  • PubChem - Vendor information and properties")
    print("  • ChEMBL - Bioactivity and drug-likeness")
    print("  • SureChEMBL - Patent landscape analysis")
    print("\nBenefits:")
    print("  • Reduce synthesis time with purchasable building blocks")
    print("  • Mitigate IP risks early in planning")
    print("  • Prioritize drug-like synthetic routes")
    print("  • Make data-driven synthesis decisions")


if __name__ == "__main__":
    asyncio.run(main())
