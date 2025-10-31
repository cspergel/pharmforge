"""
Example usage of the OpenBabel adapter
Demonstrates format conversion, 3D generation, and geometry optimization
"""
import asyncio
from adapters.openbabel import OpenBabelAdapter


async def example_basic_conversion():
    """Example: Basic SMILES to MOL2 conversion"""
    print("\n=== Example 1: Basic Format Conversion ===")

    adapter = OpenBabelAdapter()

    # Convert ethanol SMILES to MOL2
    result = await adapter.execute(
        "CCO",
        output_format="mol2",
        gen_3d=False,
        add_hydrogens=False
    )

    if result.success:
        print(f"Success! Converted to {result.data['format']}")
        print(f"Molecular formula: {result.data['properties']['molecular_formula']}")
        print(f"Molecular weight: {result.data['properties']['molecular_weight']:.2f}")
        print(f"Canonical SMILES: {result.data['canonical_smiles']}")
    else:
        print(f"Error: {result.error}")


async def example_3d_generation():
    """Example: Generate 3D coordinates with hydrogens"""
    print("\n=== Example 2: 3D Structure Generation ===")

    adapter = OpenBabelAdapter()

    # Generate 3D structure of benzene with hydrogens
    result = await adapter.execute(
        "c1ccccc1",  # Benzene
        output_format="pdb",
        gen_3d=True,
        add_hydrogens=True
    )

    if result.success:
        print(f"Success! Generated 3D structure")
        print(f"Number of atoms: {result.data['properties']['num_atoms']}")
        print(f"Has 3D coordinates: {result.data['properties']['has_3d']}")
        print(f"\nFirst 5 lines of PDB structure:")
        print('\n'.join(result.data['structure'].split('\n')[:5]))
    else:
        print(f"Error: {result.error}")


async def example_geometry_optimization():
    """Example: Optimize molecular geometry"""
    print("\n=== Example 3: Geometry Optimization ===")

    adapter = OpenBabelAdapter()

    # Optimize aspirin structure
    result = await adapter.execute(
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        output_format="mol2",
        gen_3d=True,
        add_hydrogens=True,
        optimize=True,
        force_field="mmff94"
    )

    if result.success:
        print(f"Success! Optimized geometry")
        print(f"Molecular formula: {result.data['properties']['molecular_formula']}")
        print(f"Number of atoms: {result.data['properties']['num_atoms']}")
        print(f"Number of rotatable bonds: {result.data['properties']['num_rotors']}")
        print(f"Energy: {result.data['properties']['energy']} kcal/mol")
        print(f"InChIKey: {result.data['inchikey']}")
    else:
        print(f"Error: {result.error}")


async def example_inchi_conversion():
    """Example: Convert from InChI to SMILES"""
    print("\n=== Example 4: InChI to SMILES Conversion ===")

    adapter = OpenBabelAdapter()

    # Convert InChI to canonical SMILES
    result = await adapter.execute(
        {
            "structure": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "format": "inchi"
        },
        output_format="can"
    )

    if result.success:
        print(f"Success! Converted InChI to SMILES")
        print(f"Canonical SMILES: {result.data['structure']}")
        print(f"Molecular formula: {result.data['properties']['molecular_formula']}")
    else:
        print(f"Error: {result.error}")


async def example_batch_conversion():
    """Example: Batch convert multiple molecules"""
    print("\n=== Example 5: Batch Conversion ===")

    adapter = OpenBabelAdapter()

    molecules = {
        "Ethanol": "CCO",
        "Benzene": "c1ccccc1",
        "Acetic acid": "CC(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    print("Converting molecules to SDF format with 3D coordinates...")

    for name, smiles in molecules.items():
        result = await adapter.execute(
            smiles,
            output_format="sdf",
            gen_3d=True,
            add_hydrogens=True
        )

        if result.success:
            props = result.data['properties']
            print(f"\n{name}:")
            print(f"  Formula: {props['molecular_formula']}")
            print(f"  MW: {props['molecular_weight']:.2f}")
            print(f"  Atoms: {props['num_atoms']}")
            print(f"  Rotors: {props['num_rotors']}")
        else:
            print(f"\n{name}: Failed - {result.error}")


async def example_ligand_preparation():
    """Example: Prepare ligand for docking"""
    print("\n=== Example 6: Ligand Preparation for Docking ===")

    adapter = OpenBabelAdapter()

    # Prepare a drug-like molecule for docking
    result = await adapter.execute(
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        output_format="mol2",
        gen_3d=True,
        add_hydrogens=True,
        optimize=True,
        force_field="mmff94"
    )

    if result.success:
        props = result.data['properties']
        print(f"Success! Ligand prepared for docking")
        print(f"\nMolecular properties:")
        print(f"  Formula: {props['molecular_formula']}")
        print(f"  Molecular weight: {props['molecular_weight']:.2f}")
        print(f"  Heavy atoms: {props['num_heavy_atoms']}")
        print(f"  Rotatable bonds: {props['num_rotors']}")
        print(f"  Energy: {props['energy']} kcal/mol")
        print(f"\nStructure ready for docking simulation!")
    else:
        print(f"Error: {result.error}")


async def main():
    """Run all examples"""
    print("=" * 60)
    print("OpenBabel Adapter - Example Usage")
    print("=" * 60)

    await example_basic_conversion()
    await example_3d_generation()
    await example_geometry_optimization()
    await example_inchi_conversion()
    await example_batch_conversion()
    await example_ligand_preparation()

    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
