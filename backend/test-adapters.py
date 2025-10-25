"""
Test script for both PubChem and ChEMBL adapters
Run this to test the adapters outside of the full API
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from adapters.pubchem.adapter import PubChemAdapter
from adapters.chembl.adapter import ChEMBLAdapter
from adapters.rdkit_local.adapter import RDKitAdapter


async def test_pubchem():
    """Test the PubChem adapter"""

    print("=" * 70)
    print("Testing PubChem Adapter (Molecular Properties)")
    print("=" * 70)
    print()

    adapter = PubChemAdapter()

    test_compounds = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ]

    for name, smiles in test_compounds:
        print(f"Testing {name}: {smiles}")
        print("-" * 70)

        try:
            result = await adapter.execute(smiles)

            if result.success:
                print(f"✓ Success!")
                print(f"  Molecular Weight: {result.data.get('molecular_weight')}")
                print(f"  LogP: {result.data.get('logp')}")
                print(f"  TPSA: {result.data.get('tpsa')}")
                print(f"  H-Bond Donors: {result.data.get('h_bond_donors')}")
                print(f"  H-Bond Acceptors: {result.data.get('h_bond_acceptors')}")
            else:
                print(f"✗ Failed: {result.error}")

        except Exception as e:
            print(f"✗ Error: {e}")

        print()


async def test_chembl():
    """Test the ChEMBL adapter"""

    print("=" * 70)
    print("Testing ChEMBL Adapter (Bioactivity Data)")
    print("=" * 70)
    print()

    adapter = ChEMBLAdapter()

    test_compounds = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Imatinib (Gleevec)", "CN1CCN(CC1)Cc2ccc(cc2)C(=O)Nc3ccc(c(c3)Nc4nccc(n4)c5cccnc5)C(F)(F)F"),
    ]

    for name, smiles in test_compounds:
        print(f"Testing {name}: {smiles}")
        print("-" * 70)

        try:
            result = await adapter.execute(smiles)

            if result.success:
                print(f"✓ Success!")
                print(f"  Number of Activities: {result.data.get('num_activities')}")
                print(f"  Number of Targets: {result.data.get('num_targets')}")
                print(f"  Assay Types: {result.data.get('assay_types')}")
                print(f"  Best IC50 (nM): {result.data.get('best_ic50_nm')}")
                print(f"  Best Ki (nM): {result.data.get('best_ki_nm')}")
                if result.data.get('targets'):
                    print(f"  Sample Targets: {result.data.get('targets')[:3]}")
            else:
                print(f"✗ Failed: {result.error}")

        except Exception as e:
            print(f"✗ Error: {e}")

        print()


async def test_rdkit():
    """Test the RDKit adapter"""

    print("=" * 70)
    print("Testing RDKit Adapter (Local Calculations)")
    print("=" * 70)
    print()

    adapter = RDKitAdapter()

    test_compounds = [
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ]

    for name, smiles in test_compounds:
        print(f"Testing {name}: {smiles}")
        print("-" * 70)

        try:
            result = await adapter.execute(smiles)

            if result.success:
                print(f"✓ Success!")
                print(f"  Molecular Weight: {result.data.get('molecular_weight'):.2f}")
                print(f"  LogP: {result.data.get('logp'):.2f}")
                print(f"  TPSA: {result.data.get('tpsa'):.2f}")
                print(f"  H-Bond Donors: {result.data.get('h_bond_donors')}")
                print(f"  H-Bond Acceptors: {result.data.get('h_bond_acceptors')}")
                print(f"  Rotatable Bonds: {result.data.get('rotatable_bonds')}")
                print(f"  Lipinski Violations: {result.data.get('lipinski_violations')}")
                print(f"  Num Atoms: {result.data.get('num_atoms')}")
            else:
                print(f"✗ Failed: {result.error}")

        except Exception as e:
            print(f"✗ Error: {e}")

        print()


async def main():
    """Run all adapter tests"""

    print()
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 20 + "PharmForge Adapter Tests" + " " * 24 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    # Test PubChem
    await test_pubchem()

    # Test ChEMBL (note: this may take longer due to API response times)
    await test_chembl()

    # Test RDKit (local, should be fast)
    await test_rdkit()

    print("=" * 70)
    print("All adapter tests complete!")
    print("=" * 70)
    print()


if __name__ == "__main__":
    asyncio.run(main())
