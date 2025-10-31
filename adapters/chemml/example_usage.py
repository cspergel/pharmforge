"""
ChemML Adapter - Example Usage

Demonstrates various use cases for the ChemML adapter
"""
import asyncio
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.chemml import ChemMLAdapter, CHEMML_AVAILABLE
import numpy as np


async def example_1_basic_coulomb_matrix():
    """
    Example 1: Generate Coulomb matrix for a single molecule
    """
    print("\n" + "="*70)
    print("Example 1: Basic Coulomb Matrix Generation")
    print("="*70)

    adapter = ChemMLAdapter()

    # Simple molecule: ethanol
    input_data = {
        "smiles": "CCO",
        "operation": "coulomb_matrix",
        "parameters": {
            "max_n_atoms": 23,
            "sorting": "row_norm"
        }
    }

    result = await adapter(input_data)

    if result.success:
        print(f"✓ Success!")
        print(f"  Representation type: {result.data['representation_type']}")
        print(f"  Shape: {result.data['shape']}")
        print(f"  Feature dimension: {result.data['feature_dimension']}")
        print(f"  First 10 values: {result.data['representations'][0][:10]}")
    else:
        print(f"✗ Error: {result.error}")


async def example_2_batch_processing():
    """
    Example 2: Process multiple molecules at once
    """
    print("\n" + "="*70)
    print("Example 2: Batch Processing Multiple Molecules")
    print("="*70)

    adapter = ChemMLAdapter()

    molecules = [
        "CCO",           # Ethanol
        "CC(=O)O",       # Acetic acid
        "c1ccccc1",      # Benzene
        "CC(C)O",        # Isopropanol
        "CCC(=O)O"       # Propionic acid
    ]

    input_data = {
        "smiles": molecules,
        "operation": "morgan_fingerprint",
        "parameters": {
            "radius": 2,
            "n_bits": 2048
        }
    }

    result = await adapter(input_data)

    if result.success:
        print(f"✓ Success!")
        print(f"  Total molecules: {result.data['num_molecules']}")
        print(f"  Successful: {result.data['num_successful']}")
        print(f"  Failed: {len(result.data['failed_indices'])}")
        print(f"  Feature dimension: {result.data['feature_dimension']}")
        print(f"  Output shape: {result.data['shape']}")

        # Show sparsity of fingerprints
        fps = np.array(result.data['representations'])
        sparsity = 1.0 - (np.count_nonzero(fps) / fps.size)
        print(f"  Fingerprint sparsity: {sparsity:.3f}")
    else:
        print(f"✗ Error: {result.error}")


async def example_3_rdkit_descriptors():
    """
    Example 3: Calculate comprehensive RDKit descriptors
    """
    print("\n" + "="*70)
    print("Example 3: RDKit Descriptor Calculation")
    print("="*70)

    adapter = ChemMLAdapter()

    molecules = ["CCO", "CC(=O)O", "c1ccccc1"]

    input_data = {
        "smiles": molecules,
        "operation": "rdkit_descriptors"
    }

    result = await adapter(input_data)

    if result.success:
        print(f"✓ Success!")
        print(f"  Number of descriptors: {result.data['num_descriptors']}")
        print(f"  Shape: {result.data['shape']}")
        print(f"\n  First 10 descriptor names:")
        for i, name in enumerate(result.data['descriptor_names'][:10], 1):
            print(f"    {i}. {name}")

        # Show descriptor values for first molecule
        print(f"\n  First 5 descriptor values for '{molecules[0]}':")
        for i in range(5):
            name = result.data['descriptor_names'][i]
            value = result.data['descriptors'][0][i]
            print(f"    {name}: {value:.4f}")
    else:
        print(f"✗ Error: {result.error}")


async def example_4_bag_of_bonds():
    """
    Example 4: Generate bag of bonds representation
    """
    print("\n" + "="*70)
    print("Example 4: Bag of Bonds Representation")
    print("="*70)

    adapter = ChemMLAdapter()

    molecules = ["CCO", "CCCO"]

    input_data = {
        "smiles": molecules,
        "operation": "bag_of_bonds",
        "parameters": {
            "const": 1.0
        }
    }

    result = await adapter(input_data)

    if result.success:
        print(f"✓ Success!")
        print(f"  Representation type: {result.data['representation_type']}")
        print(f"  Shape: {result.data['shape']}")
        print(f"  Molecules processed: {result.data['num_successful']}/{result.data['num_molecules']}")
        print(f"  Feature dimension: {result.data['feature_dimension']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_5_simple_ml_pipeline():
    """
    Example 5: Simple ML pipeline with ChemML features
    """
    print("\n" + "="*70)
    print("Example 5: Simple ML Pipeline")
    print("="*70)

    adapter = ChemMLAdapter()

    # Training data
    training_molecules = [
        "CCO", "CCCO", "CC(C)O", "CC(C)(C)O",
        "CC(=O)O", "CCC(=O)O", "CCCC(=O)O",
        "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1"
    ]

    # Dummy property values (e.g., solubility)
    property_values = [2.1, 2.0, 2.3, 2.5, 1.5, 1.6, 1.7, 3.2, 3.0, 2.8]

    # Generate features
    print("\nGenerating Morgan fingerprints...")
    result = await adapter({
        "smiles": training_molecules,
        "operation": "morgan_fingerprint",
        "parameters": {"radius": 2, "n_bits": 1024}
    })

    if result.success:
        X = np.array(result.data['representations'])
        y = np.array(property_values)

        print(f"✓ Features generated: {X.shape}")

        # Simple model training
        try:
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.model_selection import cross_val_score

            model = RandomForestRegressor(n_estimators=50, random_state=42)

            # Cross-validation
            scores = cross_val_score(model, X, y, cv=3, scoring='r2')

            print(f"\n✓ Model trained successfully")
            print(f"  Cross-validation R² scores: {[f'{s:.3f}' for s in scores]}")
            print(f"  Mean R²: {scores.mean():.3f} (+/- {scores.std() * 2:.3f})")

            # Fit final model
            model.fit(X, y)

            # Predict on new molecules
            test_molecules = ["CCCC", "CCCCO"]

            print(f"\nPredicting properties for new molecules...")
            test_result = await adapter({
                "smiles": test_molecules,
                "operation": "morgan_fingerprint",
                "parameters": {"radius": 2, "n_bits": 1024}
            })

            if test_result.success:
                X_test = np.array(test_result.data['representations'])
                predictions = model.predict(X_test)

                print(f"\n  Predictions:")
                for smiles, pred in zip(test_molecules, predictions):
                    print(f"    {smiles}: {pred:.2f}")

        except ImportError:
            print("\n⚠ sklearn not available - skipping model training")
    else:
        print(f"✗ Error: {result.error}")


async def example_6_comparison_of_representations():
    """
    Example 6: Compare different molecular representations
    """
    print("\n" + "="*70)
    print("Example 6: Comparison of Representations")
    print("="*70)

    adapter = ChemMLAdapter()

    molecule = "CCO"  # Ethanol

    operations = [
        ("morgan_fingerprint", {"radius": 2, "n_bits": 2048}),
        ("rdkit_descriptors", {}),
        ("coulomb_matrix", {"max_n_atoms": 23, "sorting": "row_norm"})
    ]

    print(f"Testing molecule: {molecule}")
    print()

    for operation, params in operations:
        print(f"  Testing {operation}...")

        result = await adapter({
            "smiles": molecule,
            "operation": operation,
            "parameters": params
        })

        if result.success:
            if operation == "rdkit_descriptors":
                dim = result.data['num_descriptors']
            else:
                dim = result.data['feature_dimension']

            print(f"    ✓ Feature dimension: {dim}")
            print(f"    ✓ Computation time: {result.metadata.get('computation_time', 'N/A')}")
        else:
            print(f"    ✗ Failed: {result.error}")


async def example_7_error_handling():
    """
    Example 7: Proper error handling
    """
    print("\n" + "="*70)
    print("Example 7: Error Handling")
    print("="*70)

    adapter = ChemMLAdapter()

    # Test various error cases

    # Invalid SMILES
    print("\n1. Testing with invalid SMILES...")
    result = await adapter({
        "smiles": "INVALID_SMILES_123",
        "operation": "morgan_fingerprint"
    })
    print(f"   Result: {'Success' if result.success else f'Error - {result.error}'}")

    # Unsupported operation
    print("\n2. Testing with unsupported operation...")
    result = await adapter({
        "smiles": "CCO",
        "operation": "nonexistent_operation"
    })
    print(f"   Result: {'Success' if result.success else f'Error - {result.error}'}")

    # Mixed valid/invalid molecules
    print("\n3. Testing with mixed valid/invalid molecules...")
    result = await adapter({
        "smiles": ["CCO", "INVALID", "CC(=O)O", "ALSO_INVALID"],
        "operation": "morgan_fingerprint"
    })

    if result.success:
        print(f"   ✓ Partial success")
        print(f"     Total molecules: {result.data['num_molecules']}")
        print(f"     Successful: {result.data['num_successful']}")
        print(f"     Failed indices: {result.data['failed_indices']}")
    else:
        print(f"   ✗ Complete failure: {result.error}")


async def example_8_caching():
    """
    Example 8: Demonstrate caching behavior
    """
    print("\n" + "="*70)
    print("Example 8: Caching Demonstration")
    print("="*70)

    adapter = ChemMLAdapter()

    input_data = {
        "smiles": "CCO",
        "operation": "morgan_fingerprint",
        "parameters": {"radius": 2, "n_bits": 2048}
    }

    # First call
    print("\nFirst call (should compute)...")
    result1 = await adapter(input_data, use_cache=True)
    print(f"  Cache hit: {result1.cache_hit}")
    print(f"  Cache key: {result1.metadata.get('cache_key', 'N/A')[:16]}...")

    # Second call (should hit cache)
    print("\nSecond call (should use cache)...")
    result2 = await adapter(input_data, use_cache=True)
    print(f"  Cache hit: {result2.cache_hit}")

    # Third call (no cache)
    print("\nThird call (cache disabled)...")
    result3 = await adapter(input_data, use_cache=False)
    print(f"  Cache hit: {result3.cache_hit}")


async def main():
    """
    Run all examples
    """
    print("\n" + "="*70)
    print("ChemML Adapter - Example Usage")
    print("="*70)

    if not CHEMML_AVAILABLE:
        print("\n⚠ WARNING: ChemML is not installed!")
        print("Some examples will fail. Install with: pip install chemml")
        print("\nHowever, some operations (Morgan fingerprints, RDKit descriptors)")
        print("will still work as they only require RDKit.")

    # Run examples
    examples = [
        example_1_basic_coulomb_matrix,
        example_2_batch_processing,
        example_3_rdkit_descriptors,
        example_4_bag_of_bonds,
        example_5_simple_ml_pipeline,
        example_6_comparison_of_representations,
        example_7_error_handling,
        example_8_caching
    ]

    for i, example in enumerate(examples, 1):
        try:
            await example()
        except Exception as e:
            print(f"\n✗ Example {i} failed with exception: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "="*70)
    print("Examples completed!")
    print("="*70)


if __name__ == "__main__":
    asyncio.run(main())
