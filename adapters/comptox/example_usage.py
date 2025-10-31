"""
CompTox Chemistry Dashboard Adapter - Example Usage
Demonstrates how to use the CompTox adapter for EPA chemistry and toxicity data
"""
import asyncio
import json
from adapters.comptox import CompToxAdapter


async def example_basic_search():
    """Example: Basic chemical search by name"""
    print("\n" + "="*80)
    print("EXAMPLE 1: Basic Chemical Search by Name")
    print("="*80)

    adapter = CompToxAdapter()

    # Search for aspirin
    result = await adapter.execute("aspirin")

    if result.success:
        data = result.data
        chemical = data["chemical"]

        print(f"\nChemical: {chemical['preferred_name']}")
        print(f"DTXSID: {chemical['dtxsid']}")
        print(f"CAS: {chemical['casrn']}")
        print(f"SMILES: {chemical['smiles']}")
        print(f"Molecular Formula: {chemical['molecular_formula']}")
        print(f"Molecular Weight: {chemical['molecular_weight']}")
    else:
        print(f"Error: {result.error}")


async def example_query_types():
    """Example: Search using different identifier types"""
    print("\n" + "="*80)
    print("EXAMPLE 2: Different Query Types")
    print("="*80)

    adapter = CompToxAdapter()

    # Search by CAS number
    print("\n--- Search by CAS Number ---")
    result = await adapter.execute({
        "query": "50-78-2",
        "query_type": "cas"
    })

    if result.success:
        print(f"Found: {result.data['chemical']['preferred_name']}")

    # Search by SMILES
    print("\n--- Search by SMILES ---")
    result = await adapter.execute({
        "query": "CC(=O)Oc1ccccc1C(=O)O",
        "query_type": "smiles"
    })

    if result.success:
        print(f"Found: {result.data['chemical']['preferred_name']}")

    # Search by DTXSID
    print("\n--- Search by DTXSID ---")
    result = await adapter.execute({
        "query": "DTXSID7020182",
        "query_type": "dtxsid"
    })

    if result.success:
        print(f"Found: {result.data['chemical']['preferred_name']}")


async def example_toxicity_profiling():
    """Example: Comprehensive toxicity profiling"""
    print("\n" + "="*80)
    print("EXAMPLE 3: Toxicity Profiling")
    print("="*80)

    adapter = CompToxAdapter()

    result = await adapter.execute(
        "benzene",
        include_toxicity=True,
        include_hazards=True
    )

    if result.success:
        data = result.data

        print(f"\nChemical: {data['chemical']['preferred_name']}")
        print(f"CAS: {data['chemical']['casrn']}")

        # Toxicity predictions
        print("\nToxicity Predictions (OPERA Models):")
        toxicity = data["toxicity"]
        for endpoint in toxicity["predicted_endpoints"][:10]:
            print(f"  {endpoint['endpoint']}: {endpoint['value']} ({endpoint['source']})")

        # QSAR predictions
        print("\nQSAR Predictions:")
        qsar = toxicity["qsar_predictions"]
        for key, value in qsar.items():
            print(f"  {key}: {value}")

        # Hazard classifications
        print("\nHazard Information:")
        hazard = data["hazard"]
        if hazard["ghs_classification"]:
            print(f"  GHS Codes: {', '.join(hazard['ghs_classification'])}")

        env_fate = hazard["environmental_fate"]
        if env_fate:
            print(f"  Biodegradation: {env_fate.get('biodegradation', 'Unknown')}")
            print(f"  Persistence: {env_fate.get('persistence', 'Unknown')}")


async def example_bioactivity():
    """Example: Bioactivity assay data"""
    print("\n" + "="*80)
    print("EXAMPLE 4: Bioactivity Data (ToxCast/Tox21)")
    print("="*80)

    adapter = CompToxAdapter()

    result = await adapter.execute(
        "caffeine",
        include_bioactivity=True
    )

    if result.success:
        data = result.data

        print(f"\nChemical: {data['chemical']['preferred_name']}")

        bioactivity = data["bioactivity"]
        print(f"\nTotal Assays: {bioactivity['num_assays']}")
        print(f"Active Assays: {bioactivity['active_assays']}")

        if bioactivity["gene_targets"]:
            print(f"\nGene Targets (top 10):")
            for target in bioactivity["gene_targets"][:10]:
                print(f"  - {target}")


async def example_exposure_analysis():
    """Example: Exposure pathways and use categories"""
    print("\n" + "="*80)
    print("EXAMPLE 5: Exposure Analysis")
    print("="*80)

    adapter = CompToxAdapter()

    result = await adapter.execute(
        "triclosan",
        include_exposure=True
    )

    if result.success:
        data = result.data

        print(f"\nChemical: {data['chemical']['preferred_name']}")
        print(f"CAS: {data['chemical']['casrn']}")

        exposure = data["exposure"]

        if exposure["consumer_products"]:
            print(f"\nConsumer Products:")
            for product in exposure["consumer_products"][:10]:
                print(f"  - {product}")

        if exposure["use_categories"]:
            print(f"\nUse Categories:")
            for category in exposure["use_categories"][:10]:
                print(f"  - {category}")

        if exposure["exposure_pathways"]:
            print(f"\nExposure Pathways:")
            for pathway in exposure["exposure_pathways"]:
                print(f"  - {pathway}")


async def example_properties():
    """Example: Physicochemical properties"""
    print("\n" + "="*80)
    print("EXAMPLE 6: Physicochemical Properties")
    print("="*80)

    adapter = CompToxAdapter()

    result = await adapter.execute(
        "ethanol",
        include_properties=True
    )

    if result.success:
        data = result.data

        print(f"\nChemical: {data['chemical']['preferred_name']}")

        props = data["properties"]
        print(f"\nPhysicochemical Properties:")

        if "logp" in props:
            print(f"  LogP: {props['logp']}")
        if "water_solubility" in props:
            print(f"  Water Solubility: {props['water_solubility']} g/L")
        if "vapor_pressure" in props:
            print(f"  Vapor Pressure: {props['vapor_pressure']} mmHg")
        if "melting_point" in props:
            print(f"  Melting Point: {props['melting_point']} °C")
        if "boiling_point" in props:
            print(f"  Boiling Point: {props['boiling_point']} °C")


async def example_green_chemistry_screening():
    """Example: Green chemistry screening"""
    print("\n" + "="*80)
    print("EXAMPLE 7: Green Chemistry Screening")
    print("="*80)

    adapter = CompToxAdapter()

    test_chemicals = ["benzene", "ethanol", "water"]

    for chemical in test_chemicals:
        result = await adapter.execute(
            chemical,
            include_toxicity=True,
            include_bioactivity=True,
            include_hazards=True
        )

        if result.success:
            data = result.data

            print(f"\n{data['chemical']['preferred_name']}:")

            # Assess greenness
            green_flags = []
            red_flags = []

            # Check bioactivity
            bioactivity = data.get("bioactivity", {})
            if bioactivity.get("active_assays", 0) < 10:
                green_flags.append("Low bioactivity promiscuity")
            else:
                red_flags.append(f"High bioactivity ({bioactivity['active_assays']} active assays)")

            # Check biodegradation
            hazard = data.get("hazard", {})
            env_fate = hazard.get("environmental_fate", {})
            if env_fate.get("biodegradation") == "readily biodegradable":
                green_flags.append("Readily biodegradable")
            elif env_fate.get("biodegradation"):
                red_flags.append(f"Biodegradation: {env_fate['biodegradation']}")

            # Check GHS codes
            ghs_codes = hazard.get("ghs_classification", [])
            if not ghs_codes:
                green_flags.append("No GHS hazard codes")
            else:
                red_flags.append(f"GHS hazards: {', '.join(ghs_codes[:3])}")

            # Check QSAR predictions
            toxicity = data.get("toxicity", {})
            qsar = toxicity.get("qsar_predictions", {})
            if qsar.get("mutagenicity") == "negative":
                green_flags.append("Non-mutagenic")
            elif qsar.get("mutagenicity") == "positive":
                red_flags.append("Mutagenic")

            print(f"  Green Flags: {', '.join(green_flags) if green_flags else 'None'}")
            print(f"  Red Flags: {', '.join(red_flags) if red_flags else 'None'}")

            # Overall assessment
            is_green = len(green_flags) > len(red_flags)
            print(f"  Assessment: {'GREEN CHEMISTRY CANDIDATE' if is_green else 'NOT GREEN'}")


async def example_drug_safety_screening():
    """Example: Drug candidate safety screening"""
    print("\n" + "="*80)
    print("EXAMPLE 8: Drug Candidate Safety Screening")
    print("="*80)

    adapter = CompToxAdapter()

    # Screen a drug candidate (using aspirin as example)
    result = await adapter.execute(
        "aspirin",
        include_toxicity=True,
        include_bioactivity=True,
        include_hazards=True
    )

    if result.success:
        data = result.data

        print(f"\nDrug Candidate: {data['chemical']['preferred_name']}")
        print(f"SMILES: {data['chemical']['smiles']}")

        # Safety assessment
        print("\n--- Safety Assessment ---")

        red_flags = []
        yellow_flags = []
        green_flags = []

        # Toxicity flags
        toxicity = data.get("toxicity", {})
        qsar = toxicity.get("qsar_predictions", {})

        if qsar.get("mutagenicity") == "positive":
            red_flags.append("Mutagenicity predicted")
        else:
            green_flags.append("Non-mutagenic")

        if qsar.get("developmental_toxicity") == "positive":
            red_flags.append("Developmental toxicity predicted")

        # Bioactivity promiscuity
        bioactivity = data.get("bioactivity", {})
        active_assays = bioactivity.get("active_assays", 0)
        total_assays = bioactivity.get("num_assays", 1)

        if total_assays > 0:
            promiscuity = active_assays / total_assays
            if promiscuity > 0.5:
                red_flags.append(f"High promiscuity ({active_assays}/{total_assays} active)")
            elif promiscuity > 0.3:
                yellow_flags.append(f"Moderate promiscuity ({active_assays}/{total_assays} active)")
            else:
                green_flags.append(f"Low promiscuity ({active_assays}/{total_assays} active)")

        # GHS hazards
        hazard = data.get("hazard", {})
        ghs_codes = hazard.get("ghs_classification", [])
        if ghs_codes:
            yellow_flags.append(f"GHS hazards: {', '.join(ghs_codes[:3])}")

        # Print assessment
        print(f"\nRED FLAGS (critical): {len(red_flags)}")
        for flag in red_flags:
            print(f"  - {flag}")

        print(f"\nYELLOW FLAGS (caution): {len(yellow_flags)}")
        for flag in yellow_flags:
            print(f"  - {flag}")

        print(f"\nGREEN FLAGS (positive): {len(green_flags)}")
        for flag in green_flags:
            print(f"  - {flag}")

        # Overall recommendation
        if len(red_flags) > 0:
            recommendation = "HIGH RISK - Further investigation required"
        elif len(yellow_flags) > 2:
            recommendation = "MODERATE RISK - Proceed with caution"
        else:
            recommendation = "LOW RISK - Acceptable safety profile"

        print(f"\nRecommendation: {recommendation}")


async def example_with_caching():
    """Example: Using caching for repeated queries"""
    print("\n" + "="*80)
    print("EXAMPLE 9: Caching")
    print("="*80)

    adapter = CompToxAdapter()

    # First call - fetches from API
    print("\nFirst call (from API):")
    result1 = await adapter(
        "caffeine",
        use_cache=True,
        include_toxicity=True
    )
    print(f"Cache hit: {result1.cache_hit}")
    print(f"DTXSID: {result1.data['chemical']['dtxsid']}")

    # Second call - retrieves from cache
    print("\nSecond call (from cache):")
    result2 = await adapter(
        "caffeine",
        use_cache=True,
        include_toxicity=True
    )
    print(f"Cache hit: {result2.cache_hit}")
    print(f"DTXSID: {result2.data['chemical']['dtxsid']}")


async def example_error_handling():
    """Example: Error handling"""
    print("\n" + "="*80)
    print("EXAMPLE 10: Error Handling")
    print("="*80)

    adapter = CompToxAdapter()

    # Try to query non-existent chemical
    print("\n--- Querying non-existent chemical ---")
    result = await adapter.execute("nonexistent_chemical_xyz123")

    if not result.success:
        print(f"Error: {result.error}")
        print(f"Metadata: {result.metadata}")
    else:
        # Check warnings
        if result.data.get("warnings"):
            print(f"Warnings: {', '.join(result.data['warnings'])}")


async def example_batch_processing():
    """Example: Batch processing multiple chemicals"""
    print("\n" + "="*80)
    print("EXAMPLE 11: Batch Processing")
    print("="*80)

    adapter = CompToxAdapter()

    chemicals = ["aspirin", "caffeine", "ibuprofen", "paracetamol"]

    print("\nProcessing multiple chemicals...")

    results = []
    for chemical in chemicals:
        result = await adapter.execute(
            chemical,
            include_toxicity=True
        )
        results.append((chemical, result))

    # Print summary
    print("\n--- Results Summary ---")
    print(f"{'Chemical':<20} {'DTXSID':<20} {'Toxicity Endpoints':<20} {'Status':<10}")
    print("-" * 70)

    for query, result in results:
        if result.success:
            dtxsid = result.data["chemical"]["dtxsid"]
            num_endpoints = len(result.data.get("toxicity", {}).get("predicted_endpoints", []))
            status = "SUCCESS"
        else:
            dtxsid = "N/A"
            num_endpoints = 0
            status = "FAILED"

        print(f"{query:<20} {dtxsid:<20} {num_endpoints:<20} {status:<10}")


async def main():
    """Run all examples"""
    print("\n" + "="*80)
    print("CompTox Chemistry Dashboard Adapter - Example Usage")
    print("="*80)

    examples = [
        example_basic_search,
        example_query_types,
        example_toxicity_profiling,
        example_bioactivity,
        example_exposure_analysis,
        example_properties,
        example_green_chemistry_screening,
        example_drug_safety_screening,
        example_with_caching,
        example_error_handling,
        example_batch_processing
    ]

    for example in examples:
        try:
            await example()
            await asyncio.sleep(1)  # Brief pause between examples
        except Exception as e:
            print(f"\nError running example: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "="*80)
    print("Examples complete!")
    print("="*80 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
