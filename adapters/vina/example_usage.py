"""
Example usage of the Vina Adapter

This script demonstrates how to use the VinaAdapter for molecular docking.

Requirements:
    1. AutoDock Vina installed and accessible in PATH
    2. A receptor PDBQT file
    3. Docking box coordinates (center and size)

Setup:
    # Download Vina
    # Linux/Mac:
    wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
    chmod +x vina_1.2.5_linux_x86_64
    sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

    # Windows:
    # Download from https://github.com/ccsb-scripps/AutoDock-Vina/releases
    # Add to PATH

Usage:
    python adapters/vina/example_usage.py
"""

import asyncio
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.vina.adapter import VinaAdapter


async def example_basic():
    """Basic docking example"""
    print("=" * 60)
    print("Example 1: Basic Vina Docking")
    print("=" * 60)

    # Configure adapter
    config = {
        'receptor_path': '/path/to/receptor.pdbqt',  # REPLACE WITH YOUR RECEPTOR
        'center_x': 10.0,  # REPLACE WITH YOUR DOCKING BOX CENTER
        'center_y': 20.0,
        'center_z': 15.0,
        'size_x': 25,
        'size_y': 25,
        'size_z': 25,
        'exhaustiveness': 8,
        'num_modes': 9
    }

    adapter = VinaAdapter(config=config)

    # Test molecules
    test_molecules = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    ]

    for name, smiles in test_molecules:
        print(f"\nDocking {name}...")
        print(f"SMILES: {smiles}")

        result = await adapter.execute(smiles)

        if result.success:
            data = result.data
            print(f"✓ Docking successful!")
            print(f"  Binding affinity: {data['binding_affinity']:.2f} kcal/mol")
            print(f"  Normalized score: {data['binding_score']:.3f}")
            print(f"  Number of poses: {data['num_poses']}")

            # Show all poses
            print(f"\n  All poses:")
            for pose in data['all_poses']:
                print(f"    Mode {pose['mode']}: {pose['affinity']:.2f} kcal/mol")

        else:
            print(f"✗ Docking failed: {result.error}")


async def example_with_custom_box():
    """Example with custom docking box per molecule"""
    print("\n" + "=" * 60)
    print("Example 2: Custom Docking Box per Molecule")
    print("=" * 60)

    # Create adapter with default receptor
    adapter = VinaAdapter(config={
        'receptor_path': '/path/to/receptor.pdbqt',
        'exhaustiveness': 16  # Higher exhaustiveness for better results
    })

    # Molecule with custom docking box
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

    print(f"\nDocking with custom box coordinates...")

    result = await adapter.execute(
        smiles,
        center_x=15.0,  # Override default
        center_y=22.0,
        center_z=18.0,
        size_x=20,
        size_y=20,
        size_z=20
    )

    if result.success:
        print(f"✓ Custom box docking successful!")
        print(f"  Affinity: {result.data['binding_affinity']:.2f} kcal/mol")
    else:
        print(f"✗ Docking failed: {result.error}")


async def example_batch_screening():
    """Example batch screening of multiple compounds"""
    print("\n" + "=" * 60)
    print("Example 3: Batch Virtual Screening")
    print("=" * 60)

    adapter = VinaAdapter(config={
        'receptor_path': '/path/to/receptor.pdbqt',
        'center_x': 10.0,
        'center_y': 20.0,
        'center_z': 15.0,
        'exhaustiveness': 8
    })

    # Library of compounds to screen
    compound_library = [
        ("CHEMBL1", "CCO"),
        ("CHEMBL2", "CC(=O)Oc1ccccc1C(=O)O"),
        ("CHEMBL3", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("CHEMBL4", "c1ccc2c(c1)ccc3c2cccc3"),
        ("CHEMBL5", "CC(C)Cc1ccc(cc1)C(C)C(=O)O")
    ]

    print(f"\nScreening {len(compound_library)} compounds...\n")

    results = []
    for compound_id, smiles in compound_library:
        result = await adapter.execute(smiles)

        if result.success:
            affinity = result.data['binding_affinity']
            score = result.data['binding_score']
            results.append((compound_id, smiles, affinity, score))
            print(f"✓ {compound_id}: {affinity:.2f} kcal/mol (score: {score:.3f})")
        else:
            print(f"✗ {compound_id}: Failed - {result.error}")

    # Rank by binding affinity
    if results:
        print("\n" + "-" * 60)
        print("Ranked Results (best to worst):")
        print("-" * 60)

        ranked = sorted(results, key=lambda x: x[2])  # Sort by affinity (lower=better)

        for i, (compound_id, smiles, affinity, score) in enumerate(ranked, 1):
            print(f"{i}. {compound_id}: {affinity:.2f} kcal/mol (score: {score:.3f})")


async def example_error_handling():
    """Example showing error handling"""
    print("\n" + "=" * 60)
    print("Example 4: Error Handling")
    print("=" * 60)

    adapter = VinaAdapter(config={
        'receptor_path': '/path/to/receptor.pdbqt',
        'center_x': 10.0,
        'center_y': 20.0,
        'center_z': 15.0
    })

    # Test invalid inputs
    test_cases = [
        ("Empty string", ""),
        ("Invalid SMILES", "INVALID_XYZ"),
        ("None", None),
        ("Number", 12345),
    ]

    for name, test_input in test_cases:
        print(f"\nTesting {name}: {test_input}")
        result = await adapter.execute(test_input)

        if not result.success:
            print(f"  ✓ Correctly rejected: {result.error}")
        else:
            print(f"  ✗ Unexpectedly succeeded!")


def example_metadata():
    """Example showing adapter metadata"""
    print("\n" + "=" * 60)
    print("Example 5: Adapter Metadata")
    print("=" * 60)

    adapter = VinaAdapter(config={
        'receptor_path': '/path/to/receptor.pdbqt',
        'center_x': 10.0,
        'center_y': 20.0,
        'center_z': 15.0
    })

    metadata = adapter.get_metadata()

    print(f"\nAdapter: {metadata['name']}")
    print(f"Type: {metadata['type']}")
    print(f"Version: {metadata['version']}")
    print(f"Description: {metadata['description']}")

    print(f"\nConfiguration:")
    for key, value in metadata['config'].items():
        print(f"  {key}: {value}")

    print(f"\nScoring:")
    for key, value in metadata['scoring'].items():
        print(f"  {key}: {value}")

    print(f"\nRequirements:")
    for req in metadata['requirements']:
        print(f"  - {req}")


async def main():
    """Run all examples"""
    print("\n" + "=" * 60)
    print("Vina Adapter - Example Usage")
    print("=" * 60)

    # Show metadata first
    example_metadata()

    # Note: The following examples require actual Vina setup
    print("\n" + "=" * 60)
    print("NOTE: The following examples require:")
    print("  1. AutoDock Vina installed")
    print("  2. A receptor PDBQT file")
    print("  3. Valid docking box coordinates")
    print("\nUpdate the configuration in each example before running.")
    print("=" * 60)

    # Uncomment to run actual docking examples:
    # await example_basic()
    # await example_with_custom_box()
    # await example_batch_screening()
    # await example_error_handling()


if __name__ == "__main__":
    asyncio.run(main())
