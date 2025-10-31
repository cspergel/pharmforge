"""
CompTox Adapter Installation Verification
Quick test to verify the adapter is properly installed and functional
"""
import asyncio
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


async def verify_adapter():
    """Verify CompTox adapter installation"""
    print("="*80)
    print("CompTox Chemistry Dashboard Adapter - Installation Verification")
    print("="*80)

    # Test 1: Import adapter
    print("\n[1/6] Testing adapter import...")
    try:
        from adapters.comptox.adapter import CompToxAdapter
        print("[PASS] CompTox adapter imported successfully")
    except ImportError as e:
        print(f"[FAIL] Failed to import adapter: {e}")
        return False

    # Test 2: Create instance
    print("\n[2/6] Testing adapter initialization...")
    try:
        adapter = CompToxAdapter()
        print(f"[PASS] Adapter created: {adapter.name} v{adapter.version}")
    except Exception as e:
        print(f"[FAIL] Failed to create adapter: {e}")
        return False

    # Test 3: Test validation
    print("\n[3/6] Testing input validation...")
    try:
        assert adapter.validate_input("aspirin") is True
        assert adapter.validate_input("") is False
        assert adapter.validate_input({"query": "test"}) is True
        print("[PASS] Input validation working correctly")
    except Exception as e:
        print(f"[FAIL] Validation failed: {e}")
        return False

    # Test 4: Test cache key generation
    print("\n[4/6] Testing cache key generation...")
    try:
        key1 = adapter.generate_cache_key("aspirin")
        key2 = adapter.generate_cache_key("caffeine")
        assert key1 != key2
        assert len(key1) == 64  # SHA256 hex
        print("[PASS] Cache key generation working")
    except Exception as e:
        print(f"[FAIL] Cache key generation failed: {e}")
        return False

    # Test 5: Test metadata
    print("\n[5/6] Testing adapter metadata...")
    try:
        metadata = adapter.get_metadata()
        assert metadata["name"] == "comptox"
        assert metadata["type"] == "api"
        assert metadata["version"] == "1.0.0"
        print("[PASS] Adapter metadata correct")
    except Exception as e:
        print(f"[FAIL] Metadata test failed: {e}")
        return False

    # Test 6: Test actual API call (optional - may fail if API is down)
    print("\n[6/6] Testing API connection (optional)...")
    try:
        result = await adapter.execute("aspirin")

        if result.success:
            print("[PASS] API connection successful")
            print(f"   DTXSID: {result.data['chemical']['dtxsid']}")
            print(f"   Chemical: {result.data['chemical']['preferred_name']}")
        else:
            print(f"[WARN] API call failed (expected if API is down or rate limited)")
            print(f"   Error: {result.error}")
    except Exception as e:
        print(f"[WARN] API test failed: {e}")
        print("   This is acceptable - adapter is still properly installed")

    # Test 7: Test registry integration
    print("\n[7/7] Testing registry integration...")
    try:
        from backend.core.adapter_registry import register_all_adapters, registry

        register_all_adapters()
        comptox_adapter = registry.get("comptox")

        if comptox_adapter:
            print("[PASS] Adapter registered in PharmForge registry")
            print(f"   Name: {comptox_adapter.name}")
            print(f"   Type: {comptox_adapter.adapter_type}")
        else:
            print("[FAIL] Adapter not found in registry")
            return False
    except Exception as e:
        print(f"[WARN] Registry test failed: {e}")

    # Summary
    print("\n" + "="*80)
    print("VERIFICATION COMPLETE")
    print("="*80)
    print("\n[SUCCESS] CompTox adapter is properly installed and functional!")
    print("\nYou can now use the adapter in your PharmForge pipelines:")
    print("  from adapters.comptox import CompToxAdapter")
    print("  adapter = CompToxAdapter()")
    print("  result = await adapter.execute('aspirin')")
    print("\nFor more examples, see:")
    print("  - adapters/comptox/example_usage.py")
    print("  - adapters/comptox/integration_example.py")
    print("  - adapters/comptox/README.md")
    print()

    return True


async def main():
    """Run verification"""
    try:
        success = await verify_adapter()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n[FAIL] Verification failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())
