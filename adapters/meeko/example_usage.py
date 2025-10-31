"""
Example usage of the Meeko adapter for PharmForge

This file demonstrates various use cases for the Meeko ligand preparation adapter.
"""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

from adapters.meeko import MeekoAdapter


async def basic_example():
    """Basic ligand preparation from SMILES"""
    print("\n=== Basic Ligand Preparation ===")

    # Initialize adapter
    meeko = MeekoAdapter()

    # Prepare aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"Preparing: {smiles}")

    result = await meeko.execute(smiles)

    if result.success:
        print(f"✓ Success!")
        print(f"  - PDBQT length: {len(result.data['pdbqt_string'])} chars")
        print(f"  - Rotatable bonds: {result.data['conformers'][0]['metadata']['num_torsions']}")
        print(f"  - Heavy atoms: {result.data['conformers'][0]['metadata']['num_heavy_atoms']}")
        print(f"\nFirst 500 chars of PDBQT:")
        print(result.data['pdbqt_string'][:500])
    else:
        print(f"✗ Failed: {result.error}")


async def macrocycle_example():
    """Prepare a macrocyclic compound with flexible rings"""
    print("\n=== Macrocycle Preparation ===")

    # Cyclosporin A-like macrocycle (simplified)
    smiles = "CC1CCCCCCCCCC1"  # Simple 12-membered ring
    print(f"Preparing macrocycle: {smiles}")

    # Configure for flexible macrocycles
    meeko = MeekoAdapter(config={
        "rigid_macrocycles": False,  # Allow ring flexibility
        "min_ring_size": 7  # Rings >= 7 atoms can be flexible
    })

    result = await meeko.execute(smiles)

    if result.success:
        metadata = result.data['conformers'][0]['metadata']
        print(f"✓ Success!")
        print(f"  - Macrocycles detected: {metadata['num_macrocycles']}")
        print(f"  - Macrocycle info: {metadata['macrocycles_detected']}")
        print(f"  - Rotatable bonds: {metadata['num_torsions']}")
    else:
        print(f"✗ Failed: {result.error}")


async def multi_conformer_example():
    """Generate multiple conformers for ensemble docking"""
    print("\n=== Multiple Conformer Generation ===")

    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
    print(f"Preparing: {smiles}")
    print(f"Generating 5 conformers...")

    meeko = MeekoAdapter()
    result = await meeko.execute(smiles, num_conformers=5)

    if result.success:
        print(f"✓ Success! Generated {result.data['num_conformers']} conformers")
        for conf in result.data['conformers']:
            conf_id = conf['conformer_id']
            torsions = conf['metadata']['num_torsions']
            print(f"  - Conformer {conf_id}: {torsions} rotatable bonds")
    else:
        print(f"✗ Failed: {result.error}")


async def flexible_amides_example():
    """Prepare ligand with flexible amide bonds"""
    print("\n=== Flexible Amides Example ===")

    # Peptide-like molecule with amide bonds
    smiles = "CC(=O)NC(C)C(=O)NC(C)C(=O)O"
    print(f"Preparing: {smiles}")

    # Standard (rigid amides)
    meeko_rigid = MeekoAdapter(config={"flexible_amides": False})
    result_rigid = await meeko_rigid.execute(smiles)

    # Flexible amides
    meeko_flex = MeekoAdapter(config={"flexible_amides": True})
    result_flex = await meeko_flex.execute(smiles)

    if result_rigid.success and result_flex.success:
        rigid_torsions = result_rigid.data['conformers'][0]['metadata']['num_torsions']
        flex_torsions = result_flex.data['conformers'][0]['metadata']['num_torsions']

        print(f"✓ Success!")
        print(f"  - Rigid amides: {rigid_torsions} rotatable bonds")
        print(f"  - Flexible amides: {flex_torsions} rotatable bonds")
        print(f"  - Difference: {flex_torsions - rigid_torsions} additional torsions")
    else:
        print(f"✗ One or both preparations failed")


async def metadata_example():
    """Explore adapter metadata"""
    print("\n=== Adapter Metadata ===")

    meeko = MeekoAdapter()
    metadata = meeko.get_metadata()

    print(f"Adapter: {metadata['name']} v{metadata['version']}")
    print(f"Type: {metadata['type']}")
    print(f"Description: {metadata['description']}")
    print(f"\nFeatures:")
    for feature in metadata['features']:
        print(f"  - {feature}")
    print(f"\nRequirements:")
    for req in metadata['requirements']:
        print(f"  - {req}")


async def error_handling_example():
    """Demonstrate error handling"""
    print("\n=== Error Handling ===")

    meeko = MeekoAdapter()

    # Test invalid SMILES
    invalid_smiles = "INVALID_SMILES_123"
    print(f"Testing invalid SMILES: {invalid_smiles}")

    result = await meeko.execute(invalid_smiles)

    if not result.success:
        print(f"✓ Correctly handled error: {result.error}")
    else:
        print(f"✗ Should have failed but didn't")

    # Test empty input
    print(f"\nTesting empty input...")
    result = await meeko.execute("")

    if not result.success:
        print(f"✓ Correctly handled error: {result.error}")
    else:
        print(f"✗ Should have failed but didn't")


async def batch_processing_example():
    """Process multiple molecules in batch"""
    print("\n=== Batch Processing ===")

    smiles_list = [
        "CC(=O)Oc1ccccc1C(=O)O",        # Aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine
        "CC(C)NCC(COc1ccccc1)O"         # Propranolol
    ]

    meeko = MeekoAdapter()
    results = []

    for smiles in smiles_list:
        result = await meeko.execute(smiles)
        results.append({
            "smiles": smiles,
            "success": result.success,
            "torsions": result.data['conformers'][0]['metadata']['num_torsions'] if result.success else None
        })

    print(f"Processed {len(results)} molecules:")
    for r in results:
        status = "✓" if r['success'] else "✗"
        torsions = f"{r['torsions']} torsions" if r['success'] else "failed"
        print(f"  {status} {r['smiles'][:30]}... - {torsions}")


async def advanced_config_example():
    """Advanced configuration options"""
    print("\n=== Advanced Configuration ===")

    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    # Custom configuration
    config = {
        "rigid_macrocycles": False,
        "keep_nonpolar_hydrogens": True,   # Keep all hydrogens
        "hydrate": False,                  # No hydration sites
        "flexible_amides": False,
        "min_ring_size": 7,
        "num_conformers": 1,
        "double_bond_penalty": 50.0
    }

    meeko = MeekoAdapter(config=config)
    result = await meeko.execute(smiles)

    if result.success:
        print(f"✓ Success with custom config!")
        print(f"  Configuration used:")
        for key, value in result.metadata['meeko_config'].items():
            print(f"    - {key}: {value}")


async def main():
    """Run all examples"""
    print("=" * 60)
    print("Meeko Adapter - Example Usage")
    print("=" * 60)

    try:
        await basic_example()
        await macrocycle_example()
        await multi_conformer_example()
        await flexible_amides_example()
        await metadata_example()
        await error_handling_example()
        await batch_processing_example()
        await advanced_config_example()

    except ImportError as e:
        print(f"\n✗ Missing dependency: {e}")
        print("Install with: pip install meeko rdkit")
    except Exception as e:
        print(f"\n✗ Error: {e}")

    print("\n" + "=" * 60)
    print("Examples complete!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
