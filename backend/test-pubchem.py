"""
Quick test script for PubChem adapter
Run this to test the adapter outside of the full API
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adapters.pubchem.adapter import PubChemAdapter


async def main():
    """Test the PubChem adapter with a few compounds"""

    print("=" * 60)
    print("Testing PubChem Adapter")
    print("=" * 60)
    print()

    # Create adapter instance
    adapter = PubChemAdapter()

    # Test compounds (SMILES)
    test_compounds = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
        ("Ethanol", "CCO"),
    ]

    for name, smiles in test_compounds:
        print(f"Testing {name}: {smiles}")
        print("-" * 60)

        try:
            result = await adapter.execute(smiles)

            if result.success:
                print(f"✓ Success!")
                print(f"  Molecular Weight: {result.data.get('molecular_weight')}")
                print(f"  LogP: {result.data.get('logp')}")
                print(f"  TPSA: {result.data.get('tpsa')}")
                print(f"  H-Bond Donors: {result.data.get('h_bond_donors')}")
                print(f"  H-Bond Acceptors: {result.data.get('h_bond_acceptors')}")
                print(f"  Canonical SMILES: {result.data.get('canonical_smiles')}")
            else:
                print(f"✗ Failed: {result.error}")

        except Exception as e:
            print(f"✗ Error: {e}")

        print()

    print("=" * 60)
    print("Testing complete!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
