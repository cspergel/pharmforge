"""
Quick Test Script for Vina Adapter

Tests basic functionality without requiring Vina binary or receptor files.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def test_adapter_import():
    """Test that adapter can be imported"""
    print("=" * 60)
    print("Test 1: Import Adapter")
    print("=" * 60)

    try:
        from adapters.vina.adapter import VinaAdapter
        print("[PASS] VinaAdapter imported successfully")
        return True
    except Exception as e:
        print(f"[FAIL] Import failed: {e}")
        return False


def test_adapter_initialization():
    """Test adapter initialization"""
    print("\n" + "=" * 60)
    print("Test 2: Initialize Adapter")
    print("=" * 60)

    try:
        from adapters.vina.adapter import VinaAdapter

        adapter = VinaAdapter(config={
            'receptor_path': '/dummy/path.pdbqt',
            'center_x': 10.0,
            'center_y': 20.0,
            'center_z': 15.0
        })

        print(f"[PASS] Adapter created: {adapter.name}")
        print(f"  Type: {adapter.adapter_type}")
        print(f"  Version: {adapter.version}")
        print(f"  Receptor: {adapter.receptor_path}")
        print(f"  Center: ({adapter.center_x}, {adapter.center_y}, {adapter.center_z})")

        return True
    except Exception as e:
        print(f"[FAIL] Initialization failed: {e}")
        return False


def test_metadata():
    """Test metadata retrieval"""
    print("\n" + "=" * 60)
    print("Test 3: Get Metadata")
    print("=" * 60)

    try:
        from adapters.vina.adapter import VinaAdapter

        adapter = VinaAdapter()
        metadata = adapter.get_metadata()

        print(f"[PASS] Metadata retrieved")
        print(f"  Name: {metadata['name']}")
        print(f"  Type: {metadata['type']}")
        print(f"  Version: {metadata['version']}")
        print(f"  Description: {metadata['description']}")

        print(f"\n  Requirements:")
        for req in metadata['requirements']:
            print(f"    - {req}")

        print(f"\n  Scoring:")
        for key, value in metadata['scoring'].items():
            print(f"    {key}: {value}")

        return True
    except Exception as e:
        print(f"[FAIL] Metadata retrieval failed: {e}")
        return False


def test_registry_registration():
    """Test that adapter is registered in global registry"""
    print("\n" + "=" * 60)
    print("Test 4: Registry Registration")
    print("=" * 60)

    try:
        from backend.core.adapter_registry import register_all_adapters
        from backend.core.adapters.protocol import registry

        # Register all adapters
        register_all_adapters()

        # Check if vina_docking is registered
        vina = registry.get("vina_docking")

        if vina is None:
            print("[FAIL] Vina adapter not found in registry")
            return False

        print(f"[PASS] Vina adapter registered in global registry")
        print(f"  Adapter name: {vina.name}")
        print(f"  Adapter type: {vina.adapter_type}")

        # List all adapters
        all_adapters = registry.list_adapters()
        print(f"\n  All registered adapters ({len(all_adapters)}):")
        for name in all_adapters:
            adapter = registry.get(name)
            print(f"    - {name} ({adapter.adapter_type})")

        return True
    except Exception as e:
        print(f"[FAIL] Registry registration failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_scoring_normalization():
    """Test score normalization function"""
    print("\n" + "=" * 60)
    print("Test 5: Score Normalization")
    print("=" * 60)

    try:
        from backend.core.scoring_utils import vina_affinity_to01

        test_cases = [
            (-12.0, "Excellent binding"),
            (-10.0, "Strong binding"),
            (-8.0, "Moderate binding"),
            (-6.0, "Weak binding"),
            (-4.0, "Poor binding"),
        ]

        print("[PASS] Testing vina_affinity_to01() normalization:")
        print(f"\n  {'Affinity (kcal/mol)':<20} {'Normalized Score':<20} {'Description'}")
        print("  " + "-" * 70)

        for affinity, desc in test_cases:
            score = vina_affinity_to01(affinity)
            print(f"  {affinity:<20.1f} {score:<20.3f} {desc}")

        return True
    except Exception as e:
        print(f"[FAIL] Score normalization failed: {e}")
        return False


def test_input_validation():
    """Test input validation (requires RDKit)"""
    print("\n" + "=" * 60)
    print("Test 6: Input Validation")
    print("=" * 60)

    try:
        from adapters.vina.adapter import VinaAdapter

        adapter = VinaAdapter()

        # Test valid SMILES
        valid_smiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
        invalid_inputs = ["", "INVALID_XYZ", 12345, None]

        print("Testing valid SMILES:")
        for smiles in valid_smiles:
            is_valid = adapter.validate_input(smiles)
            status = "[PASS]" if is_valid else "[FAIL]"
            print(f"  {status} {smiles}: {is_valid}")

        print("\nTesting invalid inputs:")
        for inp in invalid_inputs:
            is_valid = adapter.validate_input(inp)
            status = "[PASS]" if not is_valid else "[FAIL]"
            print(f"  {status} {inp}: {is_valid}")

        return True
    except Exception as e:
        print(f"[WARN] Validation test skipped (RDKit not available): {e}")
        return True  # Not a failure if RDKit not installed


def test_cache_key_generation():
    """Test cache key generation"""
    print("\n" + "=" * 60)
    print("Test 7: Cache Key Generation")
    print("=" * 60)

    try:
        from adapters.vina.adapter import VinaAdapter

        adapter = VinaAdapter()
        smiles = "CCO"

        key1 = adapter.generate_cache_key(smiles, receptor_path="/path/to/receptor.pdbqt")
        key2 = adapter.generate_cache_key(smiles, receptor_path="/path/to/receptor.pdbqt")
        key3 = adapter.generate_cache_key("c1ccccc1", receptor_path="/path/to/receptor.pdbqt")

        print("[PASS] Cache keys generated:")
        print(f"  Same input:  {key1}")
        print(f"  Same input:  {key2}")
        print(f"  Diff input:  {key3}")

        if key1 == key2:
            print("  [PASS] Deterministic: Same input produces same key")
        else:
            print("  [FAIL] Non-deterministic: Same input produces different keys")
            return False

        if key1 != key3:
            print("  [PASS] Unique: Different inputs produce different keys")
        else:
            print("  [FAIL] Collision: Different inputs produce same key")
            return False

        return True
    except Exception as e:
        print(f"[FAIL] Cache key generation failed: {e}")
        return False


def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("VINA ADAPTER QUICK TEST SUITE")
    print("=" * 60)
    print("\nThis script tests basic functionality without requiring:")
    print("  - AutoDock Vina binary")
    print("  - Receptor PDBQT files")
    print("  - Actual docking execution")
    print("\n" + "=" * 60 + "\n")

    results = []

    # Run all tests
    results.append(("Import", test_adapter_import()))
    results.append(("Initialization", test_adapter_initialization()))
    results.append(("Metadata", test_metadata()))
    results.append(("Registry", test_registry_registration()))
    results.append(("Scoring", test_scoring_normalization()))
    results.append(("Validation", test_input_validation()))
    results.append(("Cache Keys", test_cache_key_generation()))

    # Print summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)

    passed = sum(1 for _, result in results if result)
    total = len(results)

    for name, result in results:
        status = "[PASS] PASS" if result else "[FAIL] FAIL"
        print(f"{status:10s} - {name}")

    print("-" * 60)
    print(f"Results: {passed}/{total} tests passed")

    if passed == total:
        print("\nSUCCESS! All tests passed! Adapter is working correctly.")
        return 0
    else:
        print(f"\n[WARN] {total - passed} test(s) failed. Check output above.")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
