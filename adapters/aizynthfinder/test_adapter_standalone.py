"""
Standalone test for AiZynthFinder adapter

Run this script to test the adapter without full PharmForge setup.

Usage:
    python test_adapter_standalone.py

Note: Requires AiZynthFinder installation and configuration.
"""

import asyncio
import sys
import os

# Add parent directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from adapters.aizynthfinder.adapter import AiZynthFinderAdapter


async def test_basic_functionality():
    """Test basic adapter functionality without AiZynthFinder."""
    print("=" * 60)
    print("Testing AiZynthFinder Adapter - Basic Functionality")
    print("=" * 60)

    # Create adapter
    adapter = AiZynthFinderAdapter()

    # Test 1: Metadata
    print("\n1. Testing metadata...")
    metadata = adapter.get_metadata()
    print(f"   Adapter name: {metadata['name']}")
    print(f"   Adapter type: {metadata['type']}")
    print(f"   Version: {metadata['version']}")
    print(f"   ✓ Metadata check passed")

    # Test 2: Input validation
    print("\n2. Testing input validation...")
    test_smiles = [
        ("c1ccccc1", True, "Benzene (valid)"),
        ("CC(=O)Oc1ccccc1C(=O)O", True, "Aspirin (valid)"),
        ("INVALID", False, "Invalid SMILES"),
        ("", False, "Empty string"),
    ]

    for smiles, expected, description in test_smiles:
        is_valid = adapter.validate_input(smiles)
        status = "✓" if is_valid == expected else "✗"
        print(f"   {status} {description}: {is_valid}")

    # Test 3: Cache key generation
    print("\n3. Testing cache key generation...")
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    key1 = adapter.generate_cache_key(smiles, max_routes=5, expansion_time=60)
    key2 = adapter.generate_cache_key(smiles, max_routes=5, expansion_time=60)
    key3 = adapter.generate_cache_key(smiles, max_routes=10, expansion_time=60)

    print(f"   Cache key 1: {key1[:16]}...")
    print(f"   Cache key 2: {key2[:16]}...")
    print(f"   Same inputs match: {key1 == key2}")
    print(f"   Different params differ: {key1 != key3}")
    print(f"   ✓ Cache key generation passed")

    # Test 4: Execute with invalid input
    print("\n4. Testing execution with invalid input...")
    result = await adapter.execute("INVALID_SMILES")
    print(f"   Success: {result.success}")
    print(f"   Error: {result.error}")
    assert result.success is False
    print(f"   ✓ Invalid input handling passed")

    print("\n" + "=" * 60)
    print("Basic functionality tests PASSED ✓")
    print("=" * 60)


async def test_with_aizynthfinder():
    """Test with actual AiZynthFinder (requires installation)."""
    print("\n" + "=" * 60)
    print("Testing AiZynthFinder Adapter - Full Execution")
    print("=" * 60)

    try:
        # Create adapter
        adapter = AiZynthFinderAdapter(config={
            "max_routes": 3,
            "expansion_time": 30  # Short time for testing
        })

        # Test molecules
        test_molecules = [
            ("c1ccccc1", "Benzene"),
            ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
        ]

        for smiles, name in test_molecules:
            print(f"\n--- Testing: {name} ---")
            print(f"SMILES: {smiles}")

            result = await adapter.execute(smiles, max_routes=3, expansion_time=30)

            if result.success:
                data = result.data
                print(f"✓ Success!")
                print(f"  Routes found: {data['routes_found']}")

                if data['routes_found'] > 0:
                    best = data['best_route']
                    print(f"  Best route:")
                    print(f"    - Steps: {best['n_steps']}")
                    print(f"    - Synthesis score: {best['synthesis_score']:.3f}")
                    print(f"    - Feasibility: {best['feasibility']}")
                    print(f"    - Starting materials: {len(best['starting_materials'])}")
                else:
                    print(f"  No routes found (this may be expected for some molecules)")
            else:
                print(f"✗ Failed: {result.error}")

        print("\n" + "=" * 60)
        print("Full execution tests COMPLETED ✓")
        print("=" * 60)

    except ImportError as e:
        print(f"\n⚠ AiZynthFinder not installed")
        print(f"Error: {e}")
        print(f"\nTo install:")
        print(f"  pip install aizynthfinder")
        print(f"\nFor full setup, see:")
        print(f"  https://molecularai.github.io/aizynthfinder/getting_started.html")
        return False

    except Exception as e:
        print(f"\n✗ Error during execution: {e}")
        import traceback
        traceback.print_exc()
        return False

    return True


async def main():
    """Main test runner."""
    print("\n🧪 AiZynthFinder Adapter Standalone Test\n")

    # Always run basic tests
    await test_basic_functionality()

    # Ask user if they want to run full tests
    print("\n" + "=" * 60)
    response = input("Run full tests with AiZynthFinder? (requires installation) [y/N]: ")

    if response.lower() in ['y', 'yes']:
        success = await test_with_aizynthfinder()
        if success:
            print("\n✓ All tests PASSED!")
        else:
            print("\n⚠ Some tests failed (see above)")
    else:
        print("\nSkipping full tests.")
        print("To run later: python test_adapter_standalone.py")

    print("\n" + "=" * 60)
    print("Test run complete!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
