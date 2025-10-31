#!/usr/bin/env python3
"""
Verify Batch 8 adapter integration into backend
"""
import sys

def main():
    print("\n" + "="*70)
    print("BATCH 8 BACKEND INTEGRATION VERIFICATION")
    print("="*70 + "\n")

    # Test 1: Import and register all adapters
    print("Test 1: Registering all adapters...")
    try:
        from backend.core.adapter_registry import register_all_adapters, registry
        register_all_adapters()
        print("[OK] Registration completed\n")
    except Exception as e:
        print(f"[FAIL] Registration failed: {e}\n")
        return 1

    # Test 2: Check total adapter count
    print("Test 2: Checking adapter count...")
    total = len(registry.list_adapters())
    print(f"Total adapters registered: {total}")
    if total >= 75:
        print(f"[OK] Expected 75+, got {total}\n")
    else:
        print(f"[WARN]  Expected 75+, got {total}\n")

    # Test 3: Verify Batch 8 adapters specifically
    print("Test 3: Verifying Batch 8 adapters...")
    batch8_adapters = {
        'rnacentral': 'RNAcentral (RNA database)',
        'ord': 'ORD (Open Reaction Database)',
        'intact': 'IntAct (Protein-Protein Interactions)',
        'xtb': 'xTB (Quantum Chemistry)'
    }

    all_present = True
    for adapter_name, description in batch8_adapters.items():
        adapter = registry.get(adapter_name)
        if adapter:
            print(f"  [OK] {adapter_name}: {description}")
        else:
            print(f"  [FAIL] {adapter_name}: NOT FOUND")
            all_present = False

    if all_present:
        print("\n[OK] All Batch 8 adapters successfully integrated!\n")
    else:
        print("\n[FAIL] Some Batch 8 adapters are missing!\n")
        return 1

    # Test 4: Verify earlier batches (5-7)
    print("Test 4: Verifying Batch 5-7 adapters...")
    batch57_adapters = {
        'coconut': 'COCONUT (Natural Products)',
        'hmdb': 'HMDB (Metabolomics)',
        'tox21': 'Tox21 (Toxicity)',
        'comptox': 'CompTox (EPA Chemistry)',
        'sabdab': 'SAbDab (Antibody Structures)',
        'immunebuilder': 'ImmuneBuilder (Antibody Prediction)'
    }

    batch57_present = True
    for adapter_name, description in batch57_adapters.items():
        adapter = registry.get(adapter_name)
        if adapter:
            print(f"  [OK] {adapter_name}: {description}")
        else:
            print(f"  [FAIL] {adapter_name}: NOT FOUND")
            batch57_present = False

    if batch57_present:
        print("\n[OK] All Batch 5-7 adapters confirmed!\n")
    else:
        print("\n[WARN]  Some Batch 5-7 adapters are missing!\n")

    # Test 5: Test adapter metadata
    print("Test 5: Testing adapter metadata...")
    test_adapter = registry.get('rnacentral')
    if test_adapter:
        metadata = test_adapter.get_metadata()
        print(f"  Name: {metadata.get('name')}")
        print(f"  Type: {metadata.get('adapter_type')}")
        print(f"  Version: {metadata.get('version')}")
        print("  [OK] Metadata accessible\n")
    else:
        print("  [FAIL] Could not access metadata\n")

    # Summary
    print("="*70)
    print("INTEGRATION SUMMARY")
    print("="*70)
    print(f"Total Adapters: {total}")
    print(f"Batch 8 Adapters: {'[OK] All 4 integrated' if all_present else '[FAIL] Some missing'}")
    print(f"Batch 5-7 Adapters: {'[OK] All 6 confirmed' if batch57_present else '[WARN] Some missing'}")
    print(f"Status: {'[OK] BACKEND INTEGRATION COMPLETE' if all_present else '[FAIL] INTEGRATION INCOMPLETE'}")
    print("="*70 + "\n")

    return 0 if all_present else 1

if __name__ == "__main__":
    sys.exit(main())
