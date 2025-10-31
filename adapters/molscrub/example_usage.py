"""
MolScrub Adapter - Example Usage

Demonstrates various use cases for the MolScrub conformer generation adapter.
"""

import asyncio
import logging
from adapters.molscrub import MolScrubAdapter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


async def example_basic():
    """Basic conformer generation"""
    print("\n" + "="*60)
    print("Example 1: Basic Conformer Generation")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"Generating conformers for aspirin: {smiles}\n")

    result = await adapter.execute(smiles)

    if result.success:
        data = result.data
        print(f"✓ Success!")
        print(f"  - Generated: {data['num_generated']} conformers")
        print(f"  - After filtering: {data['num_filtered']} conformers")
        print(f"  - Lowest energy: {data['lowest_energy']:.2f} kcal/mol")

        # Show first conformer
        print(f"\n  Best conformer (ID {data['conformers'][0]['conformer_id']}):")
        print(f"    Energy: {data['conformers'][0]['energy']:.2f} kcal/mol")
        print(f"    RMSD to lowest: {data['conformers'][0]['rmsd_to_lowest']:.3f} Å")
    else:
        print(f"✗ Failed: {result.error}")


async def example_custom_parameters():
    """Generate conformers with custom parameters"""
    print("\n" + "="*60)
    print("Example 2: Custom Parameters")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Caffeine
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    print(f"Generating 20 conformers for caffeine: {smiles}\n")

    result = await adapter.execute(
        smiles,
        num_conformers=20,
        energy_window=5.0,      # Tight energy window
        rms_threshold=1.0,      # Keep diverse conformers
        output_format="mol2"    # MOL2 format
    )

    if result.success:
        data = result.data
        print(f"✓ Success!")
        print(f"  - Generated: {data['num_generated']} conformers")
        print(f"  - After filtering: {data['num_filtered']} conformers")

        # Energy statistics
        energies = [c["energy"] for c in data["conformers"]]
        print(f"\n  Energy statistics:")
        print(f"    Lowest: {min(energies):.2f} kcal/mol")
        print(f"    Highest: {max(energies):.2f} kcal/mol")
        print(f"    Range: {max(energies) - min(energies):.2f} kcal/mol")
    else:
        print(f"✗ Failed: {result.error}")


async def example_dictionary_input():
    """Using dictionary input"""
    print("\n" + "="*60)
    print("Example 3: Dictionary Input")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Ibuprofen
    input_data = {
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "num_conformers": 15,
        "energy_window": 8.0,
        "rms_threshold": 0.75,
        "optimize": True,
        "output_format": "pdb"
    }

    print(f"Generating conformers for ibuprofen\n")
    print(f"Parameters:")
    for key, value in input_data.items():
        if key != "smiles":
            print(f"  - {key}: {value}")

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data
        print(f"\n✓ Success!")
        print(f"  - Molecular weight: {data['properties']['molecular_weight']:.2f}")
        print(f"  - Heavy atoms: {data['properties']['num_heavy_atoms']}")
        print(f"  - Rotatable bonds: {data['properties']['num_rotatable_bonds']}")
        print(f"  - Conformers: {data['num_filtered']}")
    else:
        print(f"✗ Failed: {result.error}")


async def example_multiple_molecules():
    """Process multiple molecules"""
    print("\n" + "="*60)
    print("Example 4: Multiple Molecules")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Common drugs
    molecules = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1"
    }

    print(f"Processing {len(molecules)} molecules...\n")

    results = {}
    for name, smiles in molecules.items():
        result = await adapter.execute(
            smiles,
            num_conformers=10,
            energy_window=10.0
        )

        if result.success:
            data = result.data
            results[name] = {
                "smiles": smiles,
                "conformers": data["num_filtered"],
                "lowest_energy": data["lowest_energy"],
                "molecular_weight": data["properties"]["molecular_weight"]
            }
            print(f"✓ {name}: {data['num_filtered']} conformers, "
                  f"MW={data['properties']['molecular_weight']:.1f}")
        else:
            print(f"✗ {name}: Failed - {result.error}")

    # Summary
    print(f"\n{'='*60}")
    print("Summary:")
    print(f"{'='*60}")
    for name, info in results.items():
        print(f"{name:15} | Conformers: {info['conformers']:2} | "
              f"Energy: {info['lowest_energy']:7.2f} kcal/mol | "
              f"MW: {info['molecular_weight']:6.1f}")


async def example_quality_control():
    """Quality control and validation"""
    print("\n" + "="*60)
    print("Example 5: Quality Control")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Complex molecule with multiple rotatable bonds
    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen
    print(f"Quality control for: {smiles}\n")

    result = await adapter.execute(
        smiles,
        num_conformers=30,
        energy_window=12.0,
        rms_threshold=0.5
    )

    if result.success:
        data = result.data

        print(f"Molecular Properties:")
        print(f"  - Formula: C13H18O2")
        print(f"  - MW: {data['properties']['molecular_weight']:.2f}")
        print(f"  - Heavy atoms: {data['properties']['num_heavy_atoms']}")
        print(f"  - Rotatable bonds: {data['properties']['num_rotatable_bonds']}")

        print(f"\nConformer Generation:")
        print(f"  - Requested: 30")
        print(f"  - Generated: {data['num_generated']}")
        print(f"  - After filtering: {data['num_filtered']}")

        # Energy analysis
        energies = [c["energy"] for c in data["conformers"]]
        print(f"\nEnergy Analysis:")
        print(f"  - Lowest: {min(energies):.2f} kcal/mol")
        print(f"  - Highest: {max(energies):.2f} kcal/mol")
        print(f"  - Mean: {sum(energies)/len(energies):.2f} kcal/mol")
        print(f"  - Range: {max(energies) - min(energies):.2f} kcal/mol")

        # RMSD analysis
        rmsds = [c["rmsd_to_lowest"] for c in data["conformers"][1:]]
        if rmsds:
            print(f"\nRMSD to Lowest Energy:")
            print(f"  - Min: {min(rmsds):.3f} Å")
            print(f"  - Max: {max(rmsds):.3f} Å")
            print(f"  - Mean: {sum(rmsds)/len(rmsds):.3f} Å")

        # Warnings
        if data["warnings"]:
            print(f"\nWarnings:")
            for warning in data["warnings"]:
                print(f"  ⚠ {warning}")
    else:
        print(f"✗ Failed: {result.error}")


async def example_export_formats():
    """Test different export formats"""
    print("\n" + "="*60)
    print("Example 6: Export Formats")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

    formats = ["pdb", "mol2", "sdf"]

    for fmt in formats:
        print(f"Exporting to {fmt.upper()} format...")

        result = await adapter.execute(
            smiles,
            num_conformers=3,
            output_format=fmt
        )

        if result.success:
            data = result.data
            structure_sample = data["conformers"][0]["structure"][:100]
            print(f"  ✓ Success! ({data['num_filtered']} conformers)")
            print(f"  Sample output: {structure_sample}...\n")
        else:
            print(f"  ✗ Failed: {result.error}\n")


async def example_docking_preparation():
    """Prepare conformers for docking"""
    print("\n" + "="*60)
    print("Example 7: Docking Preparation")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Generate high-quality conformers for docking
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
    print(f"Preparing ligand for docking: {smiles}\n")

    result = await adapter.execute(
        smiles,
        num_conformers=5,        # Few high-quality conformers
        energy_window=5.0,       # Tight energy window
        rms_threshold=0.5,       # Remove similar
        optimize=True,           # Full optimization
        output_format="pdb"
    )

    if result.success:
        data = result.data

        print(f"✓ Ligand prepared for docking!")
        print(f"  - {data['num_filtered']} diverse conformers generated")
        print(f"  - Energy range: {data['lowest_energy']:.2f} to "
              f"{max(c['energy'] for c in data['conformers']):.2f} kcal/mol")

        # Save best conformer for docking
        best_conformer = data["conformers"][0]
        print(f"\n  Best conformer (recommended for docking):")
        print(f"    Energy: {best_conformer['energy']:.2f} kcal/mol")
        print(f"    Ready for Vina/GNINA docking")

        # Example: Save to file
        output_file = "ligand_for_docking.pdb"
        with open(output_file, "w") as f:
            f.write(best_conformer["structure"])
        print(f"\n  Saved to: {output_file}")
    else:
        print(f"✗ Failed: {result.error}")


async def example_conformer_ensemble():
    """Generate diverse conformer ensemble"""
    print("\n" + "="*60)
    print("Example 8: Conformer Ensemble for MD")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()

    # Generate diverse ensemble
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    print(f"Generating conformer ensemble: {smiles}\n")

    result = await adapter.execute(
        smiles,
        num_conformers=50,       # Many conformers
        energy_window=15.0,      # Wide window
        rms_threshold=2.0,       # Very diverse
        optimize=True,
        output_format="pdb"
    )

    if result.success:
        data = result.data

        print(f"✓ Ensemble generated!")
        print(f"  - {data['num_filtered']} diverse conformers")
        print(f"  - Covers {max(c['rmsd_to_lowest'] for c in data['conformers'][1:]):.2f} Å conformational space")

        # Energy distribution
        energies = [c["energy"] for c in data["conformers"]]
        print(f"\n  Energy distribution:")
        print(f"    Range: {min(energies):.2f} to {max(energies):.2f} kcal/mol")
        print(f"    Span: {max(energies) - min(energies):.2f} kcal/mol")

        # Save ensemble
        for i, conf in enumerate(data["conformers"]):
            filename = f"ensemble_conformer_{i:02d}.pdb"
            with open(filename, "w") as f:
                f.write(conf["structure"])

        print(f"\n  Saved {len(data['conformers'])} PDB files")
    else:
        print(f"✗ Failed: {result.error}")


async def example_metadata():
    """Display adapter metadata"""
    print("\n" + "="*60)
    print("Example 9: Adapter Metadata")
    print("="*60 + "\n")

    adapter = MolScrubAdapter()
    metadata = adapter.get_metadata()

    print(f"Adapter: {metadata['name']}")
    print(f"Type: {metadata['type']}")
    print(f"Version: {metadata['version']}")
    print(f"\nDescription: {metadata['description']}")

    print(f"\nFeatures:")
    for feature in metadata['features']:
        print(f"  - {feature}")

    print(f"\nSupported formats: {', '.join(metadata['supported_formats'])}")

    print(f"\nUse cases:")
    for use_case in metadata['use_cases']:
        print(f"  - {use_case}")

    print(f"\nDefault configuration:")
    for key, value in metadata['config'].items():
        print(f"  - {key}: {value}")


async def main():
    """Run all examples"""
    examples = [
        ("Basic Usage", example_basic),
        ("Custom Parameters", example_custom_parameters),
        ("Dictionary Input", example_dictionary_input),
        ("Multiple Molecules", example_multiple_molecules),
        ("Quality Control", example_quality_control),
        ("Export Formats", example_export_formats),
        ("Docking Preparation", example_docking_preparation),
        ("Conformer Ensemble", example_conformer_ensemble),
        ("Adapter Metadata", example_metadata),
    ]

    print("\n" + "="*60)
    print("MolScrub Adapter - Example Usage")
    print("="*60)

    for i, (name, example_func) in enumerate(examples, 1):
        try:
            await example_func()
        except Exception as e:
            print(f"\n✗ Example failed: {e}")
            logger.error(f"Error in {name}", exc_info=True)

        # Pause between examples
        if i < len(examples):
            print("\n" + "-"*60)
            await asyncio.sleep(0.5)

    print("\n" + "="*60)
    print("All examples completed!")
    print("="*60 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
