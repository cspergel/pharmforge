"""
CompTox Integration Example - Combining with other PharmForge adapters
Demonstrates how to use CompTox with PubChem, ChEMBL, and other adapters
"""
import asyncio
import json
from adapters.comptox import CompToxAdapter
from adapters.pubchem.adapter import PubChemAdapter
from adapters.chembl.adapter import ChEMBLAdapter


async def example_pubchem_comptox_integration():
    """
    Example: Get compound from PubChem, then analyze toxicity with CompTox
    Use case: Complete compound profiling
    """
    print("\n" + "="*80)
    print("INTEGRATION 1: PubChem + CompTox")
    print("Get compound properties from PubChem, then toxicity from CompTox")
    print("="*80)

    pubchem = PubChemAdapter()
    comptox = CompToxAdapter()

    compound = "caffeine"

    # Step 1: Get basic properties from PubChem
    print(f"\nStep 1: Fetching {compound} from PubChem...")
    pc_result = await pubchem.execute(compound)

    if not pc_result.success:
        print(f"PubChem query failed: {pc_result.error}")
        return

    # Step 2: Get toxicity data from CompTox using SMILES
    print(f"Step 2: Fetching toxicity data from CompTox...")
    smiles = pc_result.data["canonical_smiles"]

    ct_result = await comptox.execute({
        "query": smiles,
        "query_type": "smiles",
        "include_toxicity": True,
        "include_bioactivity": True,
        "include_hazards": True
    })

    # Step 3: Combine and analyze
    if ct_result.success:
        print("\n--- Combined Profile ---")
        print(f"Compound: {compound.title()}")
        print(f"SMILES: {smiles}")

        # PubChem data
        print("\nPhysicochemical Properties (PubChem):")
        print(f"  Molecular Weight: {pc_result.data['molecular_weight']} g/mol")
        print(f"  LogP: {pc_result.data['logp']}")
        print(f"  TPSA: {pc_result.data['tpsa']} Ų")
        print(f"  H-Bond Donors: {pc_result.data['h_bond_donors']}")
        print(f"  H-Bond Acceptors: {pc_result.data['h_bond_acceptors']}")

        # CompTox toxicity data
        print("\nToxicity Assessment (CompTox):")
        toxicity = ct_result.data.get("toxicity", {})
        qsar = toxicity.get("qsar_predictions", {})

        for key, value in qsar.items():
            print(f"  {key.replace('_', ' ').title()}: {value}")

        # CompTox bioactivity
        print("\nBioactivity (CompTox ToxCast/Tox21):")
        bioactivity = ct_result.data.get("bioactivity", {})
        print(f"  Total Assays: {bioactivity.get('num_assays', 0)}")
        print(f"  Active Assays: {bioactivity.get('active_assays', 0)}")

        if bioactivity.get("gene_targets"):
            print(f"  Gene Targets: {', '.join(bioactivity['gene_targets'][:5])}")

        # Hazard assessment
        print("\nHazard Assessment (CompTox):")
        hazard = ct_result.data.get("hazard", {})
        ghs_codes = hazard.get("ghs_classification", [])

        if ghs_codes:
            print(f"  GHS Codes: {', '.join(ghs_codes)}")
        else:
            print("  GHS Codes: None (low hazard)")

        env_fate = hazard.get("environmental_fate", {})
        if env_fate:
            print(f"  Biodegradation: {env_fate.get('biodegradation', 'Unknown')}")
    else:
        print(f"CompTox query failed: {ct_result.error}")


async def example_drug_safety_pipeline():
    """
    Example: Complete drug safety screening pipeline
    Use case: Pre-clinical candidate evaluation
    """
    print("\n" + "="*80)
    print("INTEGRATION 2: Complete Drug Safety Pipeline")
    print("Comprehensive safety assessment using multiple adapters")
    print("="*80)

    pubchem = PubChemAdapter()
    comptox = CompToxAdapter()

    # Test multiple drug candidates
    candidates = ["aspirin", "ibuprofen", "paracetamol"]

    results = []

    for candidate in candidates:
        print(f"\n--- Analyzing {candidate.title()} ---")

        # Get basic properties
        pc_result = await pubchem.execute(candidate)

        if not pc_result.success:
            print(f"  PubChem: Failed")
            continue

        smiles = pc_result.data["canonical_smiles"]

        # Get comprehensive toxicity profile
        ct_result = await comptox.execute({
            "query": smiles,
            "query_type": "smiles",
            "include_toxicity": True,
            "include_bioactivity": True,
            "include_hazards": True
        })

        if not ct_result.success:
            print(f"  CompTox: Failed")
            continue

        # Analyze safety profile
        toxicity = ct_result.data.get("toxicity", {})
        bioactivity = ct_result.data.get("bioactivity", {})
        hazard = ct_result.data.get("hazard", {})

        # Calculate safety score
        red_flags = 0
        yellow_flags = 0

        # Check QSAR predictions
        qsar = toxicity.get("qsar_predictions", {})
        if qsar.get("mutagenicity") == "positive":
            red_flags += 1
        if qsar.get("developmental_toxicity") == "positive":
            red_flags += 1

        # Check bioactivity promiscuity
        num_assays = bioactivity.get("num_assays", 1)
        active_assays = bioactivity.get("active_assays", 0)
        if num_assays > 0:
            promiscuity = active_assays / num_assays
            if promiscuity > 0.5:
                red_flags += 1
            elif promiscuity > 0.3:
                yellow_flags += 1

        # Check GHS hazards
        ghs_codes = hazard.get("ghs_classification", [])
        if len(ghs_codes) > 3:
            red_flags += 1
        elif len(ghs_codes) > 0:
            yellow_flags += 1

        # Assess
        if red_flags == 0 and yellow_flags <= 1:
            assessment = "LOW RISK"
        elif red_flags == 0:
            assessment = "MODERATE RISK"
        else:
            assessment = "HIGH RISK"

        print(f"  Assessment: {assessment}")
        print(f"  Red Flags: {red_flags}, Yellow Flags: {yellow_flags}")

        results.append({
            "compound": candidate,
            "smiles": smiles,
            "assessment": assessment,
            "red_flags": red_flags,
            "yellow_flags": yellow_flags
        })

    # Summary
    print("\n" + "="*80)
    print("SAFETY ASSESSMENT SUMMARY")
    print("="*80)
    print(f"{'Compound':<20} {'Assessment':<20} {'Red':<8} {'Yellow':<8}")
    print("-" * 60)
    for r in results:
        print(f"{r['compound'].title():<20} {r['assessment']:<20} {r['red_flags']:<8} {r['yellow_flags']:<8}")


async def example_environmental_screening():
    """
    Example: Environmental chemistry screening
    Use case: Green chemistry assessment
    """
    print("\n" + "="*80)
    print("INTEGRATION 3: Environmental Chemistry Screening")
    print("Assess environmental impact and green chemistry potential")
    print("="*80)

    comptox = CompToxAdapter()

    compounds = [
        ("benzene", "Industrial solvent"),
        ("ethanol", "Green solvent"),
        ("water", "Universal solvent")
    ]

    for compound, description in compounds:
        print(f"\n--- {compound.title()} ({description}) ---")

        result = await comptox.execute(
            compound,
            include_toxicity=True,
            include_hazards=True
        )

        if not result.success:
            print(f"  Query failed: {result.error}")
            continue

        # Environmental assessment
        hazard = result.data.get("hazard", {})
        env_fate = hazard.get("environmental_fate", {})

        print(f"  Biodegradation: {env_fate.get('biodegradation', 'Unknown')}")
        print(f"  Persistence: {env_fate.get('persistence', 'Unknown')}")

        ghs_codes = hazard.get("ghs_classification", [])
        print(f"  GHS Hazards: {len(ghs_codes)} codes")

        # Toxicity
        toxicity = result.data.get("toxicity", {})
        endpoints = toxicity.get("predicted_endpoints", [])

        if endpoints:
            print(f"  Toxicity Endpoints Available: {len(endpoints)}")

        # Green chemistry score
        is_biodegradable = env_fate.get("biodegradation") == "readily biodegradable"
        is_low_persistence = env_fate.get("persistence") == "low"
        is_low_hazard = len(ghs_codes) == 0

        green_score = sum([is_biodegradable, is_low_persistence, is_low_hazard])

        if green_score == 3:
            rating = "EXCELLENT"
        elif green_score == 2:
            rating = "GOOD"
        elif green_score == 1:
            rating = "MODERATE"
        else:
            rating = "POOR"

        print(f"  Green Chemistry Rating: {rating} ({green_score}/3)")


async def example_comparative_toxicity():
    """
    Example: Comparative toxicity analysis
    Use case: Compare toxicity profiles of similar compounds
    """
    print("\n" + "="*80)
    print("INTEGRATION 4: Comparative Toxicity Analysis")
    print("Compare toxicity profiles across compound series")
    print("="*80)

    comptox = CompToxAdapter()

    # NSAIDs comparison
    nsaids = ["aspirin", "ibuprofen", "naproxen"]

    print("\nComparing NSAIDs Toxicity Profiles...")

    profiles = []

    for drug in nsaids:
        result = await comptox.execute(
            drug,
            include_toxicity=True,
            include_bioactivity=True
        )

        if result.success:
            toxicity = result.data.get("toxicity", {})
            bioactivity = result.data.get("bioactivity", {})

            profile = {
                "drug": drug,
                "qsar": toxicity.get("qsar_predictions", {}),
                "endpoints": len(toxicity.get("predicted_endpoints", [])),
                "assays": bioactivity.get("num_assays", 0),
                "active_assays": bioactivity.get("active_assays", 0)
            }
            profiles.append(profile)

    # Display comparison
    print("\n" + "-" * 80)
    print(f"{'Drug':<15} {'Mutagenic':<12} {'Dev. Toxic':<12} {'Endpoints':<12} {'Assays':<10}")
    print("-" * 80)

    for p in profiles:
        qsar = p["qsar"]
        mut = qsar.get("mutagenicity", "unknown")
        dev = qsar.get("developmental_toxicity", "unknown")

        print(f"{p['drug'].title():<15} {mut:<12} {dev:<12} {p['endpoints']:<12} {p['assays']:<10}")


async def example_export_report():
    """
    Example: Generate comprehensive toxicity report
    Use case: Regulatory submission preparation
    """
    print("\n" + "="*80)
    print("INTEGRATION 5: Comprehensive Toxicity Report")
    print("Generate detailed report for regulatory purposes")
    print("="*80)

    comptox = CompToxAdapter()
    pubchem = PubChemAdapter()

    compound = "aspirin"

    # Get all data
    print(f"\nGenerating comprehensive report for {compound}...")

    pc_result = await pubchem.execute(compound)
    ct_result = await comptox.execute(
        compound,
        include_toxicity=True,
        include_bioactivity=True,
        include_properties=True,
        include_exposure=True,
        include_hazards=True
    )

    if pc_result.success and ct_result.success:
        # Build comprehensive report
        report = {
            "compound_name": compound,
            "identifiers": {
                "smiles": pc_result.data.get("canonical_smiles"),
                "inchi": pc_result.data.get("inchi"),
                "inchikey": pc_result.data.get("inchikey"),
                "dtxsid": ct_result.data["chemical"].get("dtxsid"),
                "casrn": ct_result.data["chemical"].get("casrn")
            },
            "physicochemical": {
                "molecular_weight": pc_result.data.get("molecular_weight"),
                "logp": pc_result.data.get("logp"),
                "tpsa": pc_result.data.get("tpsa"),
                "h_bond_donors": pc_result.data.get("h_bond_donors"),
                "h_bond_acceptors": pc_result.data.get("h_bond_acceptors")
            },
            "toxicity": ct_result.data.get("toxicity"),
            "bioactivity": ct_result.data.get("bioactivity"),
            "exposure": ct_result.data.get("exposure"),
            "hazard": ct_result.data.get("hazard"),
            "metadata": {
                "generated_by": "PharmForge",
                "data_sources": ["PubChem", "EPA CompTox"],
                "pubchem_version": pc_result.metadata.get("adapter_version"),
                "comptox_version": ct_result.metadata.get("adapter_version")
            }
        }

        # Save to file
        filename = f"comptox_report_{compound}.json"
        with open(filename, "w") as f:
            json.dump(report, f, indent=2)

        print(f"\nReport saved to: {filename}")
        print(f"\nReport Summary:")
        print(f"  Identifiers: {len(report['identifiers'])} fields")
        print(f"  Toxicity Endpoints: {len(report['toxicity']['predicted_endpoints'])} endpoints")
        print(f"  Bioactivity Assays: {report['bioactivity']['num_assays']} assays")
        print(f"  GHS Hazards: {len(report['hazard']['ghs_classification'])} codes")


async def main():
    """Run all integration examples"""
    print("\n" + "="*80)
    print("CompTox Integration Examples")
    print("Demonstrating multi-adapter workflows")
    print("="*80)

    examples = [
        example_pubchem_comptox_integration,
        example_drug_safety_pipeline,
        example_environmental_screening,
        example_comparative_toxicity,
        example_export_report
    ]

    for example in examples:
        try:
            await example()
            await asyncio.sleep(1)  # Brief pause
        except Exception as e:
            print(f"\nError in example: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "="*80)
    print("Integration examples complete!")
    print("="*80 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
