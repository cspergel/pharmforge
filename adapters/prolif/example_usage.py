"""
Example usage of the ProLIF adapter for PharmForge

This script demonstrates various use cases for the ProLIF adapter:
1. Single structure analysis (docking pose)
2. Trajectory analysis (MD simulation)
3. Custom configuration
4. Integration with other adapters
"""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.prolif.adapter import ProLIFAdapter


async def example_single_structure():
    """Example: Analyze a single protein-ligand complex (docking pose)"""
    print("\n" + "="*60)
    print("Example 1: Single Structure Analysis")
    print("="*60)

    adapter = ProLIFAdapter(config={
        'protein_selection': 'protein',
        'ligand_selection': 'resname LIG',
        'compute_frequency': True
    })

    # Example input: complex PDB file
    input_data = {
        'complex_file': '/path/to/docked_complex.pdb',
        'ligand_selection': 'resname LIG'
    }

    # Note: This would fail without actual files
    # result = await adapter.execute(input_data)

    print("Input format: Complex file with ligand selection")
    print(f"Configuration: {adapter.get_metadata()['config']}")
    print("\nExpected output:")
    print("  - Total interactions count")
    print("  - Interaction type breakdown")
    print("  - Residue-interaction pairs")
    print("  - Most common interactions")


async def example_separate_files():
    """Example: Analyze with separate protein and ligand files"""
    print("\n" + "="*60)
    print("Example 2: Separate Protein and Ligand Files")
    print("="*60)

    adapter = ProLIFAdapter()

    input_data = {
        'protein_file': '/path/to/protein.pdb',
        'ligand_file': '/path/to/ligand.mol2'
    }

    print("Input format: Separate files")
    print(f"Protein: {input_data['protein_file']}")
    print(f"Ligand: {input_data['ligand_file']}")
    print("\nThis format is useful when:")
    print("  - Protein and ligand are in different files")
    print("  - Files come from different sources")
    print("  - You want to analyze multiple ligands against same protein")


async def example_trajectory_analysis():
    """Example: Analyze protein-ligand interactions across MD trajectory"""
    print("\n" + "="*60)
    print("Example 3: MD Trajectory Analysis")
    print("="*60)

    adapter = ProLIFAdapter(config={
        'protein_selection': 'protein and around 15 resname LIG',  # Focus on binding site
        'ligand_selection': 'resname LIG',
        'compute_frequency': True
    })

    # Would require MDAnalysis Universe
    """
    import MDAnalysis as mda
    universe = mda.Universe('topology.pdb', 'trajectory.dcd')

    input_data = {
        'universe': universe,
        'ligand_selection': 'resname LIG'
    }

    result = await adapter.execute(input_data)

    if result.success:
        print(f"Analyzed {result.data['n_frames']} frames")
        print(f"Average interactions per frame: {result.data['avg_interactions_per_frame']:.2f}")
        print(f"Most persistent interactions:")
        for interaction in result.data['most_common_interactions']:
            frequency = interaction['count'] / result.data['n_frames']
            print(f"  - {interaction['type']}: {frequency*100:.1f}% of frames")
    """

    print("Input format: MDAnalysis Universe")
    print("Configuration: Focused on binding site (15Å around ligand)")
    print("\nTrajectory analysis provides:")
    print("  - Interaction frequencies across frames")
    print("  - Most persistent interactions")
    print("  - Dynamic interaction patterns")


async def example_custom_interactions():
    """Example: Analyze specific interaction types only"""
    print("\n" + "="*60)
    print("Example 4: Custom Interaction Types")
    print("="*60)

    # Only analyze hydrogen bonds and hydrophobic contacts
    adapter = ProLIFAdapter(config={
        'interactions': ['HBDonor', 'HBAcceptor', 'Hydrophobic'],
        'compute_frequency': True
    })

    print("Configured interaction types:")
    for interaction in adapter.config['interactions']:
        print(f"  - {interaction}")

    print("\nBenefits of limiting interactions:")
    print("  - Faster computation")
    print("  - Focused analysis")
    print("  - Reduced output complexity")

    print("\nAll available interaction types:")
    for interaction in adapter.INTERACTION_TYPES:
        print(f"  - {interaction}")


async def example_post_docking_workflow():
    """Example: Integration with Vina docking adapter"""
    print("\n" + "="*60)
    print("Example 5: Post-Docking Analysis Workflow")
    print("="*60)

    print("Workflow:")
    print("1. Dock ligand with Vina adapter")
    print("2. Get best pose from docking results")
    print("3. Analyze interactions with ProLIF adapter")
    print("4. Filter/rank by interaction quality")

    print("\nPseudo-code:")
    print("""
    # Step 1: Docking
    vina_adapter = VinaAdapter(config={'receptor_path': 'protein.pdbqt', ...})
    docking_result = await vina_adapter.execute(smiles)

    # Step 2: Get best pose
    best_pose_file = docking_result.data['output_file']
    binding_affinity = docking_result.data['binding_affinity']

    # Step 3: Analyze interactions
    prolif_adapter = ProLIFAdapter()
    interaction_result = await prolif_adapter.execute({
        'complex_file': best_pose_file,
        'ligand_selection': 'resname LIG'
    })

    # Step 4: Evaluate
    if interaction_result.success:
        n_hbonds = (
            interaction_result.data['interaction_counts'].get('HBDonor', 0) +
            interaction_result.data['interaction_counts'].get('HBAcceptor', 0)
        )
        n_hydrophobic = interaction_result.data['interaction_counts'].get('Hydrophobic', 0)

        print(f"Binding affinity: {binding_affinity:.2f} kcal/mol")
        print(f"H-bonds: {n_hbonds}")
        print(f"Hydrophobic contacts: {n_hydrophobic}")
        print(f"Total interactions: {interaction_result.data['total_interactions']}")
    """)


async def example_error_handling():
    """Example: Proper error handling"""
    print("\n" + "="*60)
    print("Example 6: Error Handling")
    print("="*60)

    adapter = ProLIFAdapter()

    print("Common error scenarios and handling:")
    print()

    # Invalid input
    print("1. Invalid input format:")
    result = await adapter.execute("invalid_input")
    print(f"   Success: {result.success}")
    print(f"   Error: {result.error}")
    print()

    # Missing files
    print("2. Missing files:")
    input_data = {
        'protein_file': '/nonexistent/protein.pdb',
        'ligand_file': '/nonexistent/ligand.mol2'
    }
    result = await adapter.execute(input_data)
    print(f"   Success: {result.success}")
    print(f"   Error: {result.error}")
    print()

    print("Best practices:")
    print("  - Always check result.success before accessing result.data")
    print("  - Use try-except blocks for file operations")
    print("  - Validate file paths before execution")
    print("  - Check for required dependencies at startup")


async def example_metadata():
    """Example: Retrieve adapter metadata"""
    print("\n" + "="*60)
    print("Example 7: Adapter Metadata")
    print("="*60)

    adapter = ProLIFAdapter()
    metadata = adapter.get_metadata()

    print(f"Adapter: {metadata['name']}")
    print(f"Type: {metadata['type']}")
    print(f"Version: {metadata['version']}")
    print(f"Description: {metadata['description']}")
    print()

    print("Capabilities:")
    for capability, supported in metadata['capabilities'].items():
        status = "✓" if supported else "✗"
        print(f"  {status} {capability}")
    print()

    print("Supported Interactions:")
    for interaction in metadata['supported_interactions'][:5]:  # Show first 5
        print(f"  - {interaction}")
    print(f"  ... and {len(metadata['supported_interactions']) - 5} more")
    print()

    print("Requirements:")
    for req in metadata['requirements']:
        print(f"  - {req}")


async def main():
    """Run all examples"""
    print("\n" + "="*60)
    print("ProLIF Adapter Usage Examples")
    print("="*60)

    await example_single_structure()
    await example_separate_files()
    await example_trajectory_analysis()
    await example_custom_interactions()
    await example_post_docking_workflow()
    await example_error_handling()
    await example_metadata()

    print("\n" + "="*60)
    print("Examples Complete")
    print("="*60)
    print("\nNote: Most examples show pseudo-code as they require actual")
    print("structure files and ProLIF/MDAnalysis to be installed.")
    print("\nFor working examples, ensure:")
    print("  1. pip install prolif MDAnalysis pandas")
    print("  2. Provide valid protein-ligand structure files")
    print("  3. Update file paths in the examples")


if __name__ == "__main__":
    asyncio.run(main())
