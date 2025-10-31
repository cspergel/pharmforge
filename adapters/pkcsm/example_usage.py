"""
Example usage of the pkCSM adapter

This script demonstrates how to use the pkCSM adapter for ADMET predictions.
"""

import asyncio
import sys
from pathlib import Path

# Add parent directories to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.pkcsm.adapter import PkCSMAdapter


async def example_single_prediction():
    """Example: Single molecule prediction"""
    print("="*80)
    print("Example 1: Single Molecule Prediction")
    print("="*80)

    # Initialize adapter
    adapter = PkCSMAdapter()

    # Predict ADMET for aspirin
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    print(f"\nPredicting ADMET properties for Aspirin")
    print(f"SMILES: {smiles}")

    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        predictions = data['predictions']

        print(f"\n[SUCCESS] Prediction successful!")
        print(f"Total properties predicted: {data['property_count']}")
        print(f"Canonical SMILES: {data['canonical_smiles']}")

        # Display absorption properties
        print("\nAbsorption Properties:")
        for prop, value in predictions['absorption'].items():
            print(f"  - {prop}: {value}")

        # Display toxicity properties
        print("\nToxicity Properties:")
        for prop, value in predictions['toxicity'].items():
            print(f"  - {prop}: {value}")
    else:
        print(f"[ERROR] Prediction failed: {result.error}")


async def example_batch_prediction():
    """Example: Batch prediction for multiple molecules"""
    print("\n" + "="*80)
    print("Example 2: Batch Prediction")
    print("="*80)

    adapter = PkCSMAdapter()

    # Common drug molecules
    molecules = {
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    print(f"\nPredicting ADMET for {len(molecules)} molecules...")

    for name, smiles in molecules.items():
        print(f"\n{name}:")
        result = await adapter.execute(smiles)

        if result.success:
            predictions = result.data['predictions']

            # Display key properties
            print(f"  Absorption: {predictions['absorption']['intestinal_absorption_human']:.1f}% absorbed")
            print(f"  BBB Permeability: {predictions['distribution']['bbb_permeability']:.2f}")
            print(f"  AMES toxicity: {'Yes' if predictions['toxicity']['ames_toxicity'] else 'No'}")
            print(f"  Hepatotoxicity: {'Yes' if predictions['toxicity']['hepatotoxicity'] else 'No'}")
        else:
            print(f"  Error: {result.error}")

        # Small delay between requests
        await asyncio.sleep(2)


async def example_category_access():
    """Example: Accessing properties by category"""
    print("\n" + "="*80)
    print("Example 3: Category-based Property Access")
    print("="*80)

    adapter = PkCSMAdapter()

    # Get properties by category
    print("\nAvailable property categories:")
    categories = ['absorption', 'distribution', 'metabolism', 'excretion', 'toxicity']

    for category in categories:
        props = adapter.get_properties_by_category(category)
        print(f"\n{category.upper()} ({len(props)} properties):")
        for prop in props:
            print(f"  - {prop}")


async def example_metadata():
    """Example: Accessing adapter metadata"""
    print("\n" + "="*80)
    print("Example 4: Adapter Metadata")
    print("="*80)

    adapter = PkCSMAdapter()
    metadata = adapter.get_metadata()

    print(f"\nAdapter Information:")
    print(f"  Name: {metadata['name']}")
    print(f"  Version: {metadata['version']}")
    print(f"  Type: {metadata['type']}")
    print(f"  Description: {metadata['description']}")

    print(f"\nConfiguration:")
    for key, value in metadata['config'].items():
        print(f"  {key}: {value}")

    print(f"\nNotes:")
    for note in metadata['notes']:
        print(f"  - {note}")

    print(f"\nCitation:")
    print(f"  {metadata['citation']}")


async def example_property_filtering():
    """Example: Filtering specific properties"""
    print("\n" + "="*80)
    print("Example 5: Property Filtering")
    print("="*80)

    adapter = PkCSMAdapter()

    # Predict for a molecule
    smiles = "CCO"  # Ethanol
    print(f"\nPredicting ADMET for Ethanol (CCO)...")

    result = await adapter.execute(smiles)

    if result.success:
        predictions = result.data['predictions']

        # Filter for drug-like properties
        print("\nDrug-likeness Indicators:")
        print(f"  Water Solubility: {predictions['absorption']['water_solubility']:.2f} log mol/L")
        print(f"  Intestinal Absorption: {predictions['absorption']['intestinal_absorption_human']:.1f}%")
        print(f"  BBB Permeability: {predictions['distribution']['bbb_permeability']:.2f}")

        # Check for safety red flags
        print("\nSafety Red Flags:")
        toxicity = predictions['toxicity']
        red_flags = []

        if toxicity['ames_toxicity']:
            red_flags.append("AMES positive (mutagenic)")
        if toxicity['hepatotoxicity']:
            red_flags.append("Hepatotoxic")
        if toxicity['herg_i_inhibitor'] or toxicity['herg_ii_inhibitor']:
            red_flags.append("hERG inhibitor (cardiac risk)")

        if red_flags:
            for flag in red_flags:
                print(f"  [WARNING] {flag}")
        else:
            print("  [OK] No major safety concerns detected")

        # Check CYP interactions
        print("\nCYP450 Interactions:")
        metabolism = predictions['metabolism']
        cyp_interactions = []

        if metabolism['cyp3a4_substrate']:
            cyp_interactions.append("CYP3A4 substrate")
        if metabolism['cyp3a4_inhibitor']:
            cyp_interactions.append("CYP3A4 inhibitor")

        if cyp_interactions:
            for interaction in cyp_interactions:
                print(f"  - {interaction}")
        else:
            print("  [OK] No major CYP interactions")


async def main():
    """Run all examples"""
    print("\npkCSM Adapter - Usage Examples")
    print("="*80)

    try:
        # Run examples
        await example_single_prediction()
        await example_batch_prediction()
        await example_category_access()
        await example_metadata()
        await example_property_filtering()

        print("\n" + "="*80)
        print("All examples completed successfully!")
        print("="*80)

    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(main())
