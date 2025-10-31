"""
BRENDA Integration Example
Shows how to combine BRENDA enzyme data with other PharmForge adapters
"""
import asyncio
import sys
sys.path.insert(0, '.')

from adapters.brenda import BRENDAAdapter


async def example_enzyme_target_workflow():
    """
    Complete enzyme target workflow:
    1. Query enzyme information
    2. Analyze kinetic parameters
    3. Assess druggability
    4. Identify lead compounds
    """
    print("=" * 70)
    print("BRENDA Integration Example: Enzyme Target Workflow")
    print("=" * 70)

    brenda = BRENDAAdapter()

    # Target: Alcohol dehydrogenase (EC 1.1.1.1)
    ec_number = "1.1.1.1"
    print(f"\nAnalyzing target enzyme: EC {ec_number}\n")

    # Step 1: Get enzyme information
    print("Step 1: Retrieving enzyme information...")
    result = await brenda.execute({
        "query": ec_number,
        "parameters": ["km", "kcat", "ki"],
        "organisms": ["Homo sapiens"],
        "max_results": 50
    })

    if not result.success:
        print(f"Error: {result.error}")
        return

    enzyme = result.data["enzyme"]
    kinetics = result.data["kinetics"]
    inhibitors = result.data["inhibitors"]

    print(f"✓ Found: {enzyme['common_name']}")
    print(f"  Systematic name: {enzyme['systematic_name']}")
    print(f"  Reaction: {enzyme['reaction']}")
    print(f"  Substrates: {', '.join(enzyme['substrates'][:3])}")

    # Step 2: Analyze kinetic parameters
    print("\nStep 2: Analyzing kinetic parameters...")

    km_values = [k for k in kinetics if k['parameter'] == 'km']
    kcat_values = [k for k in kinetics if k['parameter'] == 'kcat']

    if km_values:
        print(f"✓ Found {len(km_values)} Km values")
        avg_km = sum(k['value'] for k in km_values) / len(km_values)
        print(f"  Average Km: {avg_km:.2f} mM")
        print(f"  Range: {min(k['value'] for k in km_values):.2f} - {max(k['value'] for k in km_values):.2f} mM")
    else:
        print("  No Km values available")

    if kcat_values:
        print(f"✓ Found {len(kcat_values)} Kcat values")
        avg_kcat = sum(k['value'] for k in kcat_values) / len(kcat_values)
        print(f"  Average Kcat: {avg_kcat:.1f} s⁻¹")
    else:
        print("  No Kcat values available")

    # Calculate catalytic efficiency
    if km_values and kcat_values:
        print("\n  Catalytic Efficiency (Kcat/Km):")
        for km in km_values[:3]:  # Show top 3
            kcat = next(
                (k for k in kcat_values if k['substrate'] == km['substrate']),
                None
            )
            if kcat:
                efficiency = kcat['value'] / km['value']
                print(f"    {km['substrate']}: {efficiency:.2f} mM⁻¹s⁻¹")

    # Step 3: Assess druggability
    print("\nStep 3: Druggability Assessment...")

    # Count potent inhibitors
    potent_inhibitors = [i for i in inhibitors if i['ki_value'] < 100]  # < 100 nM
    competitive = [i for i in inhibitors if i['inhibition_type'] == 'competitive']

    print(f"  Total inhibitors: {len(inhibitors)}")
    print(f"  Potent inhibitors (Ki < 100 nM): {len(potent_inhibitors)}")
    print(f"  Competitive inhibitors: {len(competitive)}")

    # Calculate druggability score
    score = 0
    if len(inhibitors) >= 5:
        score += 1
    if len(potent_inhibitors) >= 1:
        score += 1
    if len(competitive) >= 1:
        score += 1
    if len(kinetics) >= 5:
        score += 1

    print(f"\n  Druggability Score: {score}/4")

    if score >= 3:
        assessment = "✓ STRONG TARGET - Proceed with drug design"
    elif score >= 2:
        assessment = "~ MODERATE TARGET - Further validation recommended"
    else:
        assessment = "? WEAK TARGET - Consider alternatives"

    print(f"  Assessment: {assessment}")

    # Step 4: Identify lead compounds
    print("\nStep 4: Lead Compound Identification...")

    if potent_inhibitors:
        print(f"  Top {min(3, len(potent_inhibitors))} lead compounds:")
        for i, inh in enumerate(sorted(potent_inhibitors, key=lambda x: x['ki_value'])[:3], 1):
            print(f"    {i}. {inh['compound']}")
            print(f"       Ki: {inh['ki_value']} {inh['unit']}")
            print(f"       Type: {inh['inhibition_type']}")
            print(f"       Reference: {inh['reference']}")
    else:
        print("  No potent lead compounds identified")
        print("  Recommendation: Design de novo inhibitors based on enzyme structure")

    # Next steps
    print("\n" + "=" * 70)
    print("Recommended Next Steps:")
    print("=" * 70)

    if potent_inhibitors:
        print("1. Obtain 3D structure from PDB/AlphaFold adapter")
        print("2. Dock lead compounds using Vina/GNINA adapter")
        print("3. Calculate ADMET properties using TDC adapter")
        print("4. Optimize leads using retrosynthesis adapter")
    else:
        print("1. Get enzyme structure from AlphaFold/PDB adapter")
        print("2. Identify binding pocket with ProLIF adapter")
        print("3. Design de novo inhibitors with REINVENT adapter")
        print("4. Screen virtual library with docking adapter")

    print("\n✓ Workflow complete!")


async def example_compare_cytochrome_p450():
    """
    Compare different CYP450 enzymes for drug metabolism prediction
    """
    print("\n" + "=" * 70)
    print("CYP450 Enzyme Comparison for Drug Metabolism")
    print("=" * 70)

    brenda = BRENDAAdapter()

    cyp_enzymes = {
        "1.14.13.39": "CYP2D6",  # Metabolizes ~25% of drugs
        # Add more CYP enzymes as mock data is expanded
    }

    print("\nAnalyzing CYP450 enzymes involved in drug metabolism...\n")

    for ec, name in cyp_enzymes.items():
        result = await brenda.execute({
            "query": ec,
            "parameters": ["km", "ki"],
            "organisms": ["Homo sapiens"]
        })

        if result.success:
            enzyme = result.data["enzyme"]
            kinetics = result.data["kinetics"]
            inhibitors = result.data["inhibitors"]

            print(f"{name} (EC {ec}):")
            print(f"  Substrates: {', '.join(enzyme['substrates'][:5])}")

            if kinetics:
                km_values = [k['value'] for k in kinetics if k['parameter'] == 'km']
                if km_values:
                    print(f"  Average Km: {sum(km_values)/len(km_values):.2f} mM")

            if inhibitors:
                print(f"  Known inhibitors: {len(inhibitors)}")
                potent = [i for i in inhibitors if i['ki_value'] < 1000]  # < 1 µM
                if potent:
                    print(f"  Potent inhibitors: {len(potent)}")

            print()


async def example_metabolic_pathway():
    """
    Analyze enzymes in a metabolic pathway
    """
    print("\n" + "=" * 70)
    print("Metabolic Pathway Analysis")
    print("=" * 70)

    brenda = BRENDAAdapter()

    # Example: Alcohol metabolism pathway
    pathway = [
        ("1.1.1.1", "Alcohol dehydrogenase", "Ethanol → Acetaldehyde"),
        # Add more pathway steps as mock data is expanded
    ]

    print("\nAnalyzing alcohol metabolism pathway:\n")

    for ec, name, reaction in pathway:
        result = await brenda.execute(ec)

        if result.success:
            enzyme = result.data["enzyme"]
            kinetics = result.data["kinetics"]

            print(f"Step: {reaction}")
            print(f"  Enzyme: {name} (EC {ec})")

            if kinetics:
                km_values = [k for k in kinetics if k['parameter'] == 'km']
                if km_values:
                    avg_km = sum(k['value'] for k in km_values) / len(km_values)
                    print(f"  Km: {avg_km:.2f} mM (average)")

            print()

    print("✓ Pathway analysis complete")


async def main():
    """Run all integration examples"""
    examples = [
        example_enzyme_target_workflow,
        example_compare_cytochrome_p450,
        example_metabolic_pathway
    ]

    for example in examples:
        try:
            await example()
        except Exception as e:
            print(f"\nError in {example.__name__}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("All integration examples completed!")
    print("=" * 70)


if __name__ == "__main__":
    asyncio.run(main())
