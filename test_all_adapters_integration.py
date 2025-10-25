"""
Comprehensive integration test for ALL PharmForge adapters
Tests real compounds through the complete adapter stack
"""
import asyncio
import sys
from pathlib import Path

# Add backend to path
sys.path.insert(0, str(Path(__file__).parent / "backend"))

print("=" * 80)
print("PHARMFORGE ADAPTER INTEGRATION TEST")
print("=" * 80)

# Test compounds
TEST_COMPOUNDS = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
}

async def test_all_adapters():
    """Test all adapters with real compounds"""

    print("\n" + "=" * 80)
    print("STEP 1: Import and Initialize")
    print("=" * 80)

    try:
        from backend.core.adapter_registry import registry, register_all_adapters
        print("✓ Imported adapter registry")
    except Exception as e:
        print(f"✗ Failed to import: {e}")
        return

    # Register all adapters
    try:
        register_all_adapters()
        print("✓ Registered all adapters")
    except Exception as e:
        print(f"✗ Failed to register: {e}")
        return

    # List registered adapters
    print(f"\nRegistered adapters: {len(registry.adapters)}")
    for name, adapter in registry.adapters.items():
        meta = adapter.get_metadata()
        print(f"  - {name} ({meta.get('type', 'unknown')})")

    print("\n" + "=" * 80)
    print("STEP 2: Test Individual Adapters")
    print("=" * 80)

    test_smiles = TEST_COMPOUNDS["aspirin"]
    print(f"\nTest compound: Aspirin ({test_smiles})")

    results = {}

    # Test each adapter
    adapters_to_test = [
        ("rdkit_local", "RDKit Local"),
        ("pubchem", "PubChem API"),
        ("chembl", "ChEMBL API"),
        ("admet_ai", "ADMET-AI"),
        ("vina_docking", "Vina Docking"),
        ("aizynthfinder", "AiZynthFinder"),
    ]

    for adapter_name, display_name in adapters_to_test:
        print(f"\n{'─' * 80}")
        print(f"Testing: {display_name} ({adapter_name})")
        print('─' * 80)

        adapter = registry.get(adapter_name)
        if not adapter:
            print(f"✗ Adapter '{adapter_name}' not found in registry")
            results[adapter_name] = {"error": "not_found"}
            continue

        try:
            # Validate input
            validation = adapter.validate_input(test_smiles)
            if not validation.is_valid:
                print(f"✗ Input validation failed: {validation.error}")
                results[adapter_name] = {"error": "validation_failed", "message": validation.error}
                continue
            print(f"✓ Input validation passed")

            # Generate cache key
            cache_key = adapter.generate_cache_key(test_smiles)
            print(f"✓ Cache key: {cache_key[:16]}...")

            # Execute adapter
            print(f"  Executing {adapter_name}...")
            result = await adapter.execute(test_smiles)

            if result.success:
                print(f"✓ Execution successful")
                print(f"  Data keys: {list(result.data.keys()) if result.data else 'None'}")

                # Show key results
                if result.data:
                    if adapter_name == "rdkit_local":
                        print(f"  MW: {result.data.get('molecular_weight', 'N/A'):.2f}")
                        print(f"  LogP: {result.data.get('logp', 'N/A'):.2f}")
                    elif adapter_name == "admet_ai":
                        props = result.data.get('properties', {})
                        print(f"  Properties: {len(props)}")
                        if 'Caco2_Wang' in props:
                            print(f"  Caco2: {props['Caco2_Wang']:.3f}")
                    elif adapter_name == "vina_docking":
                        print(f"  Binding affinity: {result.data.get('binding_affinity', 'N/A')}")
                        print(f"  Binding score: {result.data.get('binding_score', 'N/A'):.3f}")
                    elif adapter_name == "aizynthfinder":
                        print(f"  Steps: {result.data.get('n_steps', 'N/A')}")
                        print(f"  Synthesis score: {result.data.get('synthesis_score', 'N/A'):.3f}")
                        print(f"  Feasibility: {result.data.get('feasibility', 'N/A')}")

                results[adapter_name] = {"success": True, "data": result.data}
            else:
                print(f"✗ Execution failed: {result.error}")
                results[adapter_name] = {"error": "execution_failed", "message": result.error}

        except Exception as e:
            print(f"✗ Exception: {type(e).__name__}: {str(e)}")
            results[adapter_name] = {"error": "exception", "message": str(e)}

    print("\n" + "=" * 80)
    print("STEP 3: Test Custom Filters")
    print("=" * 80)

    try:
        from adapters.custom.filters.toxicity import PAINSEvaluator, BrenkAlertEvaluator
        from adapters.custom.filters.druglikeness import VeberEvaluator
        from adapters.custom.synthesis import SynthesisAccessibilityEvaluator

        print("\nTesting custom filters with Aspirin...")

        # PAINS
        pains = PAINSEvaluator()
        pains_result = pains.evaluate(test_smiles)
        print(f"\n✓ PAINS: {pains_result['pass']} (score: {pains_result['score']})")

        # Brenk
        brenk = BrenkAlertEvaluator()
        brenk_result = brenk.evaluate(test_smiles)
        print(f"✓ Brenk: {brenk_result['pass']} (risk: {brenk_result['summary']['brenk_risk_level']})")

        # Veber
        veber = VeberEvaluator()
        veber_result = veber.evaluate(test_smiles)
        print(f"✓ Veber: {veber_result['pass']} (bioavailability: {veber_result['summary']['oral_bioavailability']})")

        # SAScore
        sascore = SynthesisAccessibilityEvaluator()
        sascore_result = sascore.evaluate(test_smiles)
        print(f"✓ SAScore: {sascore_result['score']:.2f} (difficulty: {sascore_result['summary']['synthesis_difficulty']})")

    except Exception as e:
        print(f"✗ Custom filters error: {e}")

    print("\n" + "=" * 80)
    print("STEP 4: Summary")
    print("=" * 80)

    success_count = sum(1 for r in results.values() if r.get("success"))
    error_count = len(results) - success_count

    print(f"\nTotal adapters tested: {len(results)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {error_count}")

    print("\nDetailed results:")
    for adapter_name, result in results.items():
        status = "✓ PASS" if result.get("success") else "✗ FAIL"
        error_msg = f" ({result.get('message', result.get('error'))})" if not result.get("success") else ""
        print(f"  {status} - {adapter_name}{error_msg}")

    print("\n" + "=" * 80)
    print("INTEGRATION TEST COMPLETE")
    print("=" * 80)

    return results

if __name__ == "__main__":
    results = asyncio.run(test_all_adapters())
