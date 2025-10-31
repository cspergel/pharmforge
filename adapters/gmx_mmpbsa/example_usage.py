"""
Example Usage of gmx_MMPBSA Adapter for PharmForge

This script demonstrates various use cases for the gmx_MMPBSA adapter,
including basic MM/PBSA calculations, MM/GBSA for faster screening,
and per-residue decomposition analysis.
"""

import asyncio
import os
import logging
from adapters.gmx_mmpbsa import GmxMMPBSAAdapter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


async def example_basic_mmpbsa():
    """
    Example 1: Basic MM/PBSA calculation

    Calculate binding free energy using Poisson-Boltzmann solvation model.
    This is the most accurate method but slower.
    """
    logger.info("=" * 80)
    logger.info("Example 1: Basic MM/PBSA Calculation")
    logger.info("=" * 80)

    # Initialize adapter with MM/PBSA configuration
    adapter = GmxMMPBSAAdapter(
        config={
            "method": "pb",  # Poisson-Boltzmann
            "saltcon": 0.15,  # 150 mM salt concentration
            "temperature": 298.15,  # 298.15 K (25°C)
            "startframe": 1,  # Analyze from first frame
            "endframe": -1,  # Analyze all frames
            "interval": 1  # Analyze every frame
        }
    )

    # Prepare input data
    # NOTE: Replace these paths with your actual file paths
    input_data = {
        "topology_file": "data/md/complex.prmtop",  # AMBER topology
        "trajectory_file": "data/md/production.xtc",  # GROMACS trajectory
        "index_file": "data/md/index.ndx",  # GROMACS index file
        "receptor_group": "Protein",  # Receptor group name from index
        "ligand_group": "Ligand"  # Ligand group name from index
    }

    # Execute calculation
    result = await adapter.execute(input_data)

    # Process results
    if result.success:
        data = result.data

        logger.info("\n--- RESULTS ---")
        logger.info(f"Method: MM/{data['method']}SA")
        logger.info(f"Binding Free Energy (ΔG): {data['binding_free_energy']:.2f} kcal/mol")

        if data['delta_h'] is not None:
            logger.info(f"Enthalpy (ΔH): {data['delta_h']:.2f} kcal/mol")

        if data['delta_s'] is not None:
            logger.info(f"Entropy (-TΔS): {data['delta_s']:.2f} kcal/mol")

        logger.info("\n--- ENERGY COMPONENTS ---")
        components = data['components']
        logger.info(f"Van der Waals:      {components.get('van_der_waals', 0):.2f} kcal/mol")
        logger.info(f"Electrostatic:      {components.get('electrostatic', 0):.2f} kcal/mol")
        logger.info(f"Polar Solvation:    {components.get('polar_solvation', 0):.2f} kcal/mol")
        logger.info(f"SASA (Non-polar):   {components.get('sasa', 0):.2f} kcal/mol")

        # Interpret results
        logger.info("\n--- INTERPRETATION ---")
        dg = data['binding_free_energy']
        if dg < -10:
            logger.info("Strong binding (ΔG < -10 kcal/mol)")
        elif dg < -7:
            logger.info("Moderate binding (-10 < ΔG < -7 kcal/mol)")
        elif dg < -5:
            logger.info("Weak binding (-7 < ΔG < -5 kcal/mol)")
        else:
            logger.info("Very weak or no binding (ΔG > -5 kcal/mol)")

        if data['warnings']:
            logger.warning(f"Warnings: {data['warnings']}")

    else:
        logger.error(f"Calculation failed: {result.error}")


async def example_mmgbsa_screening():
    """
    Example 2: MM/GBSA for faster screening

    Use Generalized Born solvation model for quick initial validation.
    10x faster than MM/PBSA, suitable for screening many compounds.
    """
    logger.info("\n" + "=" * 80)
    logger.info("Example 2: MM/GBSA Screening (Fast Mode)")
    logger.info("=" * 80)

    # Initialize adapter with MM/GBSA configuration
    adapter = GmxMMPBSAAdapter(
        config={
            "method": "gb",  # Generalized Born
            "igb": 5,  # GB model (5 = mbondi2, recommended)
            "saltcon": 0.15,
            "temperature": 298.15,
            "startframe": 10,  # Skip first 10 frames (equilibration)
            "endframe": -1,  # Use all remaining frames
            "interval": 5  # Analyze every 5th frame for speed
        }
    )

    # Input data (same structure as MM/PBSA)
    input_data = {
        "topology_file": "data/md/complex.prmtop",
        "trajectory_file": "data/md/production.xtc",
        "index_file": "data/md/index.ndx",
        "receptor_group": "Protein",
        "ligand_group": "Ligand"
    }

    # Execute calculation
    result = await adapter.execute(input_data)

    if result.success:
        data = result.data

        logger.info("\n--- MM/GBSA RESULTS ---")
        logger.info(f"Binding Free Energy (ΔG): {data['binding_free_energy']:.2f} kcal/mol")
        logger.info(f"Frames analyzed: {data['frames_analyzed']}")
        logger.info(f"Method: MM/GBSA (igb={adapter.config['igb']})")

        logger.info("\nNote: MM/GBSA is faster but less accurate than MM/PBSA.")
        logger.info("Use for initial screening, then validate hits with MM/PBSA.")

    else:
        logger.error(f"Calculation failed: {result.error}")


async def example_residue_decomposition():
    """
    Example 3: Per-residue energy decomposition

    Identify which residues contribute most to binding.
    Useful for structure-based optimization and hotspot identification.
    """
    logger.info("\n" + "=" * 80)
    logger.info("Example 3: Per-Residue Decomposition")
    logger.info("=" * 80)

    # Initialize adapter with decomposition enabled
    adapter = GmxMMPBSAAdapter(
        config={
            "method": "gb",  # Use GB for faster decomposition
            "decomp": True,  # Enable per-residue decomposition
            "igb": 5,
            "saltcon": 0.15,
            "startframe": 50,  # Use equilibrated portion
            "endframe": 150,  # Limit frames for faster decomposition
            "interval": 5
        }
    )

    input_data = {
        "topology_file": "data/md/complex.prmtop",
        "trajectory_file": "data/md/production.xtc",
        "index_file": "data/md/index.ndx",
        "receptor_group": "Protein",
        "ligand_group": "Ligand"
    }

    result = await adapter.execute(input_data)

    if result.success:
        data = result.data

        logger.info("\n--- OVERALL BINDING ENERGY ---")
        logger.info(f"ΔG_bind: {data['binding_free_energy']:.2f} kcal/mol")

        if data['per_residue_decomposition']:
            logger.info("\n--- KEY RESIDUES CONTRIBUTING TO BINDING ---")
            logger.info("Top 10 residues (most favorable):")

            # Sort by energy (most negative = most favorable)
            sorted_residues = sorted(
                data['per_residue_decomposition'],
                key=lambda x: x.get('energy', 0)
            )

            for i, residue in enumerate(sorted_residues[:10], 1):
                logger.info(
                    f"  {i:2d}. {residue['name']:8s} "
                    f"{residue['energy']:8.2f} kcal/mol"
                )

            logger.info("\nThese residues are candidates for mutation studies or optimization.")

    else:
        logger.error(f"Calculation failed: {result.error}")


async def example_batch_screening():
    """
    Example 4: Batch screening multiple compounds

    Demonstrate how to screen multiple protein-ligand complexes
    from MD simulations efficiently.
    """
    logger.info("\n" + "=" * 80)
    logger.info("Example 4: Batch Screening Multiple Compounds")
    logger.info("=" * 80)

    # Initialize adapter for fast screening
    adapter = GmxMMPBSAAdapter(
        config={
            "method": "gb",  # Fast method for screening
            "igb": 5,
            "saltcon": 0.15,
            "startframe": 10,
            "endframe": -1,
            "interval": 10  # Every 10th frame for speed
        }
    )

    # List of compounds to screen
    compounds = [
        {
            "name": "Compound_A",
            "topology_file": "data/md/compound_a/complex.prmtop",
            "trajectory_file": "data/md/compound_a/production.xtc",
            "index_file": "data/md/compound_a/index.ndx",
            "receptor_group": "Protein",
            "ligand_group": "Ligand"
        },
        {
            "name": "Compound_B",
            "topology_file": "data/md/compound_b/complex.prmtop",
            "trajectory_file": "data/md/compound_b/production.xtc",
            "index_file": "data/md/compound_b/index.ndx",
            "receptor_group": "Protein",
            "ligand_group": "Ligand"
        },
        {
            "name": "Compound_C",
            "topology_file": "data/md/compound_c/complex.prmtop",
            "trajectory_file": "data/md/compound_c/production.xtc",
            "index_file": "data/md/compound_c/index.ndx",
            "receptor_group": "Protein",
            "ligand_group": "Ligand"
        }
    ]

    results = []

    # Screen all compounds
    for compound in compounds:
        logger.info(f"\nScreening: {compound['name']}")

        input_data = {
            "topology_file": compound["topology_file"],
            "trajectory_file": compound["trajectory_file"],
            "index_file": compound["index_file"],
            "receptor_group": compound["receptor_group"],
            "ligand_group": compound["ligand_group"]
        }

        result = await adapter.execute(input_data)

        if result.success:
            dg = result.data['binding_free_energy']
            results.append({
                "name": compound['name'],
                "dg": dg,
                "components": result.data['components']
            })
            logger.info(f"  ΔG = {dg:.2f} kcal/mol")
        else:
            logger.error(f"  Failed: {result.error}")
            results.append({
                "name": compound['name'],
                "dg": None,
                "error": result.error
            })

    # Rank compounds by binding affinity
    logger.info("\n--- RANKING ---")
    successful_results = [r for r in results if r['dg'] is not None]
    ranked = sorted(successful_results, key=lambda x: x['dg'])

    for i, compound in enumerate(ranked, 1):
        logger.info(f"{i}. {compound['name']:15s} ΔG = {compound['dg']:7.2f} kcal/mol")

    logger.info("\nRecommendation: Validate top hits with MM/PBSA for higher accuracy.")


async def example_comparison_pb_vs_gb():
    """
    Example 5: Compare MM/PBSA vs MM/GBSA

    Demonstrate accuracy vs speed tradeoff between methods.
    """
    logger.info("\n" + "=" * 80)
    logger.info("Example 5: Comparison of MM/PBSA vs MM/GBSA")
    logger.info("=" * 80)

    input_data = {
        "topology_file": "data/md/complex.prmtop",
        "trajectory_file": "data/md/production.xtc",
        "index_file": "data/md/index.ndx",
        "receptor_group": "Protein",
        "ligand_group": "Ligand"
    }

    # Common configuration
    common_config = {
        "saltcon": 0.15,
        "temperature": 298.15,
        "startframe": 10,
        "endframe": 100,  # Use same frame range
        "interval": 5
    }

    # Test MM/PBSA
    logger.info("\nRunning MM/PBSA...")
    adapter_pb = GmxMMPBSAAdapter(
        config={**common_config, "method": "pb"}
    )
    result_pb = await adapter_pb.execute(input_data)

    # Test MM/GBSA
    logger.info("\nRunning MM/GBSA...")
    adapter_gb = GmxMMPBSAAdapter(
        config={**common_config, "method": "gb", "igb": 5}
    )
    result_gb = await adapter_gb.execute(input_data)

    # Compare results
    logger.info("\n--- COMPARISON ---")

    if result_pb.success:
        logger.info(f"MM/PBSA:  ΔG = {result_pb.data['binding_free_energy']:7.2f} kcal/mol")
    else:
        logger.error(f"MM/PBSA failed: {result_pb.error}")

    if result_gb.success:
        logger.info(f"MM/GBSA:  ΔG = {result_gb.data['binding_free_energy']:7.2f} kcal/mol")
    else:
        logger.error(f"MM/GBSA failed: {result_gb.error}")

    if result_pb.success and result_gb.success:
        diff = abs(result_pb.data['binding_free_energy'] - result_gb.data['binding_free_energy'])
        logger.info(f"\nDifference: {diff:.2f} kcal/mol")

        logger.info("\nInterpretation:")
        if diff < 2.0:
            logger.info("  - Good agreement between methods")
            logger.info("  - MM/GBSA suitable for screening")
        elif diff < 5.0:
            logger.info("  - Moderate agreement")
            logger.info("  - Use MM/PBSA for final validation")
        else:
            logger.info("  - Large difference detected")
            logger.info("  - System may have complex solvation")
            logger.info("  - Recommend MM/PBSA for accuracy")


async def main():
    """
    Run all examples

    Note: Examples will fail if input files don't exist.
    This is expected - examples are for demonstration purposes.
    """
    logger.info("gmx_MMPBSA Adapter - Example Usage")
    logger.info("=" * 80)
    logger.info("\nNOTE: These examples require MD simulation files.")
    logger.info("If files are not available, examples will fail with file errors.")
    logger.info("Modify file paths to match your actual data location.\n")

    try:
        # Example 1: Basic MM/PBSA
        await example_basic_mmpbsa()

        # Example 2: MM/GBSA screening
        await example_mmgbsa_screening()

        # Example 3: Per-residue decomposition
        await example_residue_decomposition()

        # Example 4: Batch screening
        await example_batch_screening()

        # Example 5: Method comparison
        await example_comparison_pb_vs_gb()

    except Exception as e:
        logger.error(f"Error running examples: {e}")
        logger.info("\nThis is expected if input files don't exist.")
        logger.info("Create test data or modify file paths in the examples.")

    logger.info("\n" + "=" * 80)
    logger.info("Examples complete!")
    logger.info("=" * 80)


if __name__ == "__main__":
    asyncio.run(main())
