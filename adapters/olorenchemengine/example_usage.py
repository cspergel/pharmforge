"""
Example usage of the Oloren ChemEngine adapter for PharmForge

Demonstrates various use cases and features of the adapter.
"""

import asyncio
import json
from adapters.olorenchemengine.adapter import OlorenChemEngineAdapter


async def example_1_basic_prediction():
    """Example 1: Basic single molecule prediction"""
    print("\n" + "="*70)
    print("Example 1: Basic Single Molecule Prediction")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    # Predict properties for aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    result = await adapter.execute(smiles)

    if result.success:
        prediction = result.data["predictions"][0]
        print(f"\nSMILES: {prediction['smiles']}")
        print("\nPredicted Properties:")
        for prop, values in prediction["properties"].items():
            print(f"  {prop:20s}: {values['value']:7.3f} {values['unit']}")
            if 'uncertainty' in values:
                print(f"  {'':20s}  (uncertainty: ±{values['uncertainty']:.3f}, "
                      f"confidence: {values['confidence']:.3f})")
    else:
        print(f"Error: {result.error}")


async def example_2_batch_predictions():
    """Example 2: Batch predictions for multiple molecules"""
    print("\n" + "="*70)
    print("Example 2: Batch Predictions")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    # Common drug molecules
    molecules = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Ethanol": "CCO"
    }

    smiles_list = list(molecules.values())
    result = await adapter.execute(smiles_list)

    if result.success:
        print(f"\nPredicted properties for {result.data['num_molecules']} molecules\n")

        for i, (name, smiles) in enumerate(molecules.items()):
            prediction = result.data["predictions"][i]
            print(f"{name} ({smiles})")

            # Show selected properties
            for prop in ["solubility", "logp", "permeability"]:
                if prop in prediction["properties"]:
                    values = prediction["properties"][prop]
                    print(f"  {prop:20s}: {values['value']:7.3f} {values['unit']}")

            print()


async def example_3_specific_properties():
    """Example 3: Predict specific properties only"""
    print("\n" + "="*70)
    print("Example 3: Specific Properties")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

    # Predict only ADMET-relevant properties
    result = await adapter.execute(
        smiles,
        properties=["solubility", "permeability", "bioavailability"]
    )

    if result.success:
        prediction = result.data["predictions"][0]
        print(f"\nADMET Profile for: {smiles}")
        print(f"Properties predicted: {result.data['properties_predicted']}\n")

        for prop, values in prediction["properties"].items():
            print(f"{prop}:")
            print(f"  Value: {values['value']:.3f} {values['unit']}")
            if 'confidence' in values:
                print(f"  Confidence: {values['confidence']:.1%}")


async def example_4_toxicity_screening():
    """Example 4: Toxicity screening"""
    print("\n" + "="*70)
    print("Example 4: Toxicity Screening")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    compounds = {
        "Compound_A": "CC(=O)Oc1ccccc1C(=O)O",
        "Compound_B": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Compound_C": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    }

    # Screen for toxicity
    smiles_list = list(compounds.values())
    result = await adapter.execute(
        smiles_list,
        properties=["herg", "ames"]
    )

    if result.success:
        print("\nToxicity Screening Results\n")
        print(f"{'Compound':<15} {'hERG (pIC50)':<15} {'AMES (prob)':<15} {'Risk'}")
        print("-" * 60)

        for i, (name, smiles) in enumerate(compounds.items()):
            prediction = result.data["predictions"][i]

            herg = prediction["properties"]["herg"]["value"]
            ames = prediction["properties"]["ames"]["value"]

            # Risk assessment (simplified)
            risk = "High" if herg > 6.0 or ames > 0.5 else "Low"

            print(f"{name:<15} {herg:7.3f}         {ames:7.3f}         {risk}")


async def example_5_model_comparison():
    """Example 5: Compare different model architectures"""
    print("\n" + "="*70)
    print("Example 5: Model Architecture Comparison")
    print("="*70)

    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    property_name = "solubility"

    models = ["default", "gcn", "attentivefp", "mpnn"]
    results = {}

    for model in models:
        adapter = OlorenChemEngineAdapter(
            config={"model": model}
        )
        result = await adapter.execute(
            smiles,
            properties=[property_name]
        )

        if result.success:
            prediction = result.data["predictions"][0]
            results[model] = prediction["properties"][property_name]

    print(f"\nSolubility predictions for Aspirin using different models:\n")
    print(f"{'Model':<15} {'Value':<12} {'Uncertainty':<12} {'Confidence'}")
    print("-" * 60)

    for model, values in results.items():
        print(f"{model:<15} {values['value']:7.3f}     "
              f"±{values.get('uncertainty', 0.0):.3f}      "
              f"{values.get('confidence', 0.0):.1%}")


async def example_6_uncertainty_analysis():
    """Example 6: Analyze prediction uncertainty"""
    print("\n" + "="*70)
    print("Example 6: Uncertainty Analysis")
    print("="*70)

    adapter = OlorenChemEngineAdapter(
        config={"include_uncertainty": True}
    )

    compounds = [
        "CC(=O)Oc1ccccc1C(=O)O",  # Well-studied (Aspirin)
        "Cc1c(c(=O)n([nH]1)C)CCN(CC)CC"  # Less common structure
    ]

    result = await adapter.execute(
        compounds,
        properties=["solubility", "logp"]
    )

    if result.success:
        print("\nUncertainty Analysis\n")

        for i, prediction in enumerate(result.data["predictions"]):
            print(f"Compound {i+1}: {prediction['smiles']}")

            for prop, values in prediction["properties"].items():
                uncertainty = values.get('uncertainty', 0.0)
                confidence = values.get('confidence', 0.0)

                status = "High confidence" if confidence > 0.8 else \
                         "Medium confidence" if confidence > 0.6 else \
                         "Low confidence"

                print(f"  {prop:15s}: {values['value']:7.3f} ± {uncertainty:.3f} "
                      f"({status})")
            print()


async def example_7_property_categories():
    """Example 7: Explore property categories"""
    print("\n" + "="*70)
    print("Example 7: Property Categories")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    # Get properties by category
    categories = adapter.get_property_categories()

    print("\nAvailable Properties by Category:\n")

    for category, properties in sorted(categories.items()):
        print(f"{category.upper()}:")
        for prop in properties:
            prop_info = adapter.get_property_info(prop)
            print(f"  - {prop:20s}: {prop_info['description']} ({prop_info['unit']})")
        print()


async def example_8_metadata_inspection():
    """Example 8: Inspect adapter metadata"""
    print("\n" + "="*70)
    print("Example 8: Adapter Metadata")
    print("="*70)

    adapter = OlorenChemEngineAdapter()
    metadata = adapter.get_metadata()

    print("\nAdapter Information:")
    print(f"  Name: {metadata['name']}")
    print(f"  Type: {metadata['type']}")
    print(f"  Version: {metadata['version']}")
    print(f"  Description: {metadata['description']}")

    print("\nCapabilities:")
    for key, value in metadata['capabilities'].items():
        if isinstance(value, list):
            print(f"  {key}: {len(value)} available")
        else:
            print(f"  {key}: {value}")

    print("\nConfiguration:")
    for key, value in metadata['config'].items():
        print(f"  {key}: {value}")


async def example_9_error_handling():
    """Example 9: Error handling"""
    print("\n" + "="*70)
    print("Example 9: Error Handling")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    # Test 1: Invalid SMILES
    print("\nTest 1: Invalid SMILES")
    result = await adapter.execute("invalid_smiles_string")
    print(f"Success: {result.success}")
    print(f"Error: {result.error}")

    # Test 2: Unsupported property
    print("\nTest 2: Unsupported Property")
    result = await adapter.execute(
        "CCO",
        properties=["unknown_property"]
    )
    print(f"Success: {result.success}")
    print(f"Error: {result.error}")

    # Test 3: Empty input
    print("\nTest 3: Empty Input")
    result = await adapter.execute("")
    print(f"Success: {result.success}")
    print(f"Error: {result.error}")


async def example_10_cache_demonstration():
    """Example 10: Cache key generation"""
    print("\n" + "="*70)
    print("Example 10: Cache Key Generation")
    print("="*70)

    adapter = OlorenChemEngineAdapter()

    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    # Same input, different parameters = different cache keys
    params_list = [
        {"properties": ["solubility", "logp"]},
        {"properties": ["solubility"]},
        {"properties": ["solubility", "logp"], "model": "gcn"},
        {"properties": ["solubility", "logp"], "include_uncertainty": False}
    ]

    print("\nCache keys for different parameter combinations:\n")

    for i, params in enumerate(params_list, 1):
        cache_key = adapter.generate_cache_key(smiles, **params)
        print(f"Config {i}: {params}")
        print(f"Cache key: {cache_key[:32]}...")
        print()


async def main():
    """Run all examples"""
    examples = [
        ("Basic Prediction", example_1_basic_prediction),
        ("Batch Predictions", example_2_batch_predictions),
        ("Specific Properties", example_3_specific_properties),
        ("Toxicity Screening", example_4_toxicity_screening),
        ("Model Comparison", example_5_model_comparison),
        ("Uncertainty Analysis", example_6_uncertainty_analysis),
        ("Property Categories", example_7_property_categories),
        ("Metadata Inspection", example_8_metadata_inspection),
        ("Error Handling", example_9_error_handling),
        ("Cache Demonstration", example_10_cache_demonstration),
    ]

    print("\n" + "="*70)
    print("OLOREN CHEMENGINE ADAPTER - EXAMPLE USAGE")
    print("="*70)

    for name, example_func in examples:
        try:
            await example_func()
        except Exception as e:
            print(f"\nError in {name}: {e}")

    print("\n" + "="*70)
    print("All examples completed!")
    print("="*70)


if __name__ == "__main__":
    asyncio.run(main())
