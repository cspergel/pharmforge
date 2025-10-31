"""
Integration Example: gmx_MMPBSA with PharmForge Pipeline

This example demonstrates how to integrate gmx_MMPBSA into a complete
PharmForge drug discovery workflow:

1. Docking with Vina
2. MD simulation (external step)
3. Binding free energy calculation with gmx_MMPBSA
4. Result comparison and ranking
"""

import asyncio
import logging
from typing import List, Dict

from adapters.vina import VinaAdapter
from adapters.gmx_mmpbsa import GmxMMPBSAAdapter
from adapters.rdkit_local import RDKitAdapter

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


async def complete_workflow_example():
    """
    Complete workflow: Docking → MD → MM/PBSA

    This example shows how to combine multiple PharmForge adapters
    for a comprehensive binding affinity assessment.
    """
    logger.info("=" * 80)
    logger.info("COMPLETE WORKFLOW: Docking → MD → MM/PBSA")
    logger.info("=" * 80)

    # Test compound (aspirin)
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    compound_name = "Aspirin"

    logger.info(f"\nCompound: {compound_name}")
    logger.info(f"SMILES: {smiles}")

    # Step 1: Calculate molecular properties
    logger.info("\n" + "-" * 80)
    logger.info("STEP 1: Calculate Molecular Properties")
    logger.info("-" * 80)

    rdkit = RDKitAdapter()
    properties = await rdkit.execute(smiles)

    if properties.success:
        logger.info("Molecular Properties:")
        logger.info(f"  MW: {properties.data['molecular_weight']:.2f}")
        logger.info(f"  LogP: {properties.data['logP']:.2f}")
        logger.info(f"  HBD: {properties.data['hbd']}")
        logger.info(f"  HBA: {properties.data['hba']}")
        logger.info(f"  TPSA: {properties.data['tpsa']:.2f}")
    else:
        logger.error(f"Property calculation failed: {properties.error}")
        return

    # Step 2: Molecular docking
    logger.info("\n" + "-" * 80)
    logger.info("STEP 2: Molecular Docking (Vina)")
    logger.info("-" * 80)

    # Note: In production, you would have a real receptor file
    vina = VinaAdapter(
        config={
            "receptor_path": "data/receptors/protein.pdbqt",
            "center_x": 15.0,
            "center_y": 20.0,
            "center_z": 10.0,
            "size_x": 25,
            "size_y": 25,
            "size_z": 25
        }
    )

    # This would normally run with real receptor
    # For this example, we'll simulate the expected output
    logger.info("Note: This example requires receptor file for actual docking")
    logger.info("Expected Vina output:")
    logger.info("  Binding Affinity: -7.5 kcal/mol")
    logger.info("  Best Pose: Mode 1")

    docking_affinity = -7.5  # Simulated result

    # Step 3: MD Simulation Setup
    logger.info("\n" + "-" * 80)
    logger.info("STEP 3: MD Simulation (External)")
    logger.info("-" * 80)

    logger.info("The following steps would be performed externally:")
    logger.info("  1. Extract docked pose from Vina output")
    logger.info("  2. Prepare MD system (solvate, add ions)")
    logger.info("  3. Energy minimization")
    logger.info("  4. Equilibration (NVT, NPT)")
    logger.info("  5. Production MD (10-100 ns)")
    logger.info("  6. Create GROMACS index file")
    logger.info("\nFor this example, we assume MD files are ready:")
    logger.info("  - complex.prmtop")
    logger.info("  - production.xtc")
    logger.info("  - index.ndx")

    # Step 4: MM/PBSA Calculation
    logger.info("\n" + "-" * 80)
    logger.info("STEP 4: Binding Free Energy (MM/PBSA)")
    logger.info("-" * 80)

    mmpbsa = GmxMMPBSAAdapter(
        config={
            "method": "pb",  # Accurate MM/PBSA
            "startframe": 10,  # Skip equilibration
            "endframe": -1,  # Use all remaining frames
            "interval": 5,  # Every 5th frame for speed
            "saltcon": 0.15,
            "temperature": 298.15
        }
    )

    # Prepare input (would use real files in production)
    mmpbsa_input = {
        "topology_file": "data/md/complex.prmtop",
        "trajectory_file": "data/md/production.xtc",
        "index_file": "data/md/index.ndx",
        "receptor_group": "Protein",
        "ligand_group": "Ligand"
    }

    logger.info("Note: This example requires MD simulation files")
    logger.info("Expected MM/PBSA output:")
    logger.info("  ΔG_bind: -12.3 kcal/mol")
    logger.info("  ΔH: -18.5 kcal/mol")
    logger.info("  -TΔS: 6.2 kcal/mol")
    logger.info("\nEnergy Components:")
    logger.info("  Van der Waals: -35.2 kcal/mol")
    logger.info("  Electrostatic: -28.3 kcal/mol")
    logger.info("  Polar Solvation: 48.5 kcal/mol")
    logger.info("  SASA: -3.5 kcal/mol")

    mmpbsa_affinity = -12.3  # Simulated result

    # Step 5: Result Comparison
    logger.info("\n" + "-" * 80)
    logger.info("STEP 5: Results Summary")
    logger.info("-" * 80)

    logger.info(f"\nCompound: {compound_name}")
    logger.info(f"Molecular Weight: {properties.data['molecular_weight']:.2f}")
    logger.info(f"\nBinding Affinity Estimates:")
    logger.info(f"  Vina Docking:  {docking_affinity:7.2f} kcal/mol  (Fast, approximate)")
    logger.info(f"  MM/PBSA:       {mmpbsa_affinity:7.2f} kcal/mol  (Accurate, MD-based)")

    logger.info(f"\nDifference: {abs(docking_affinity - mmpbsa_affinity):.2f} kcal/mol")

    if abs(docking_affinity - mmpbsa_affinity) < 3.0:
        logger.info("Interpretation: Good agreement between methods")
    else:
        logger.info("Interpretation: Significant difference - MD captures additional effects")

    logger.info("\nConclusion:")
    if mmpbsa_affinity < -10:
        logger.info("  ✓ Strong binder (ΔG < -10 kcal/mol)")
        logger.info("  ✓ Proceed to experimental validation")
    elif mmpbsa_affinity < -7:
        logger.info("  ~ Moderate binder (-10 < ΔG < -7 kcal/mol)")
        logger.info("  ~ Consider optimization")
    else:
        logger.info("  ✗ Weak binder (ΔG > -7 kcal/mol)")
        logger.info("  ✗ Not suitable for this target")


async def batch_screening_example():
    """
    Batch screening with quick MM/GBSA validation

    Demonstrates how to screen multiple compounds efficiently
    using fast MM/GBSA calculations.
    """
    logger.info("\n" + "=" * 80)
    logger.info("BATCH SCREENING: Quick MM/GBSA Validation")
    logger.info("=" * 80)

    # Test compounds
    compounds = [
        {"name": "Compound_A", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        {"name": "Compound_B", "smiles": "CCO"},
        {"name": "Compound_C", "smiles": "c1ccccc1"},
    ]

    # Initialize adapters
    mmpbsa_screening = GmxMMPBSAAdapter(
        config={
            "method": "gb",  # Fast Generalized Born
            "igb": 5,
            "startframe": 10,
            "endframe": 100,  # Limited frames for speed
            "interval": 10  # Every 10th frame
        }
    )

    results = []

    for compound in compounds:
        logger.info(f"\nScreening: {compound['name']}")
        logger.info(f"  SMILES: {compound['smiles']}")

        # Simulated MM/GBSA calculation
        # In production, this would use real MD files
        logger.info("  Status: Would run MM/GBSA calculation")
        logger.info("  Expected runtime: 1-2 minutes per compound")

        # Simulated results
        simulated_dg = -8.0 + (hash(compound['name']) % 10) * 0.5
        results.append({
            "name": compound['name'],
            "dg": simulated_dg,
            "method": "MM/GBSA (simulated)"
        })

        logger.info(f"  ΔG: {simulated_dg:.2f} kcal/mol")

    # Rank compounds
    logger.info("\n" + "-" * 80)
    logger.info("RANKING (by ΔG)")
    logger.info("-" * 80)

    ranked = sorted(results, key=lambda x: x['dg'])
    for i, compound in enumerate(ranked, 1):
        logger.info(f"{i}. {compound['name']:15s} ΔG = {compound['dg']:7.2f} kcal/mol")

    logger.info("\nRecommendation:")
    logger.info(f"  → Validate top 3 hits with accurate MM/PBSA")
    logger.info(f"  → Consider experimental testing for ΔG < -8 kcal/mol")


async def hotspot_analysis_example():
    """
    Per-residue decomposition for hotspot analysis

    Identifies key residues contributing to binding for
    structure-based optimization.
    """
    logger.info("\n" + "=" * 80)
    logger.info("HOTSPOT ANALYSIS: Per-Residue Decomposition")
    logger.info("=" * 80)

    mmpbsa_decomp = GmxMMPBSAAdapter(
        config={
            "method": "gb",  # Use GB for faster decomposition
            "decomp": True,  # Enable per-residue decomposition
            "igb": 5,
            "startframe": 50,
            "endframe": 150,
            "interval": 5
        }
    )

    logger.info("\nConfiguration:")
    logger.info("  Method: MM/GBSA with per-residue decomposition")
    logger.info("  Frames: 50-150, interval=5")

    logger.info("\nNote: This example requires MD simulation files")
    logger.info("Expected output:")

    logger.info("\nOverall Binding Energy:")
    logger.info("  ΔG_bind: -11.5 kcal/mol")

    logger.info("\nTop Contributing Residues:")
    logger.info("  1. TYR 105   -3.2 kcal/mol  (π-π stacking)")
    logger.info("  2. ASP 189   -2.8 kcal/mol  (salt bridge)")
    logger.info("  3. PHE 234   -2.1 kcal/mol  (hydrophobic)")
    logger.info("  4. ARG 156   -1.9 kcal/mol  (H-bond donor)")
    logger.info("  5. TRP 290   -1.6 kcal/mol  (π-π stacking)")

    logger.info("\nInterpretation:")
    logger.info("  → TYR 105 is a key hotspot residue")
    logger.info("  → ASP 189 provides electrostatic stabilization")
    logger.info("  → PHE 234 contributes hydrophobic interactions")

    logger.info("\nStructure-Based Design Recommendations:")
    logger.info("  1. Optimize π-π interactions with TYR 105")
    logger.info("  2. Strengthen salt bridge with ASP 189")
    logger.info("  3. Consider mutations at TYR 105 for resistance studies")


async def workflow_comparison():
    """
    Compare different workflow strategies

    Shows when to use different approaches based on project needs.
    """
    logger.info("\n" + "=" * 80)
    logger.info("WORKFLOW COMPARISON: Choosing the Right Approach")
    logger.info("=" * 80)

    workflows = [
        {
            "name": "Fast Screening",
            "steps": ["Docking (Vina)", "Short MD (10 ns)", "MM/GBSA"],
            "time": "~4 hours per compound",
            "accuracy": "Moderate",
            "use_case": "Initial screening of 100+ compounds"
        },
        {
            "name": "Balanced",
            "steps": ["Docking (Vina)", "Medium MD (50 ns)", "MM/PBSA"],
            "time": "~12 hours per compound",
            "accuracy": "Good",
            "use_case": "Lead optimization (10-20 compounds)"
        },
        {
            "name": "High Accuracy",
            "steps": ["Docking", "Long MD (100+ ns)", "MM/PBSA", "Decomposition"],
            "time": "~24 hours per compound",
            "accuracy": "Excellent",
            "use_case": "Final validation (top 3-5 candidates)"
        }
    ]

    for i, workflow in enumerate(workflows, 1):
        logger.info(f"\n{i}. {workflow['name']}")
        logger.info(f"   Steps: {' → '.join(workflow['steps'])}")
        logger.info(f"   Time: {workflow['time']}")
        logger.info(f"   Accuracy: {workflow['accuracy']}")
        logger.info(f"   Use case: {workflow['use_case']}")

    logger.info("\n" + "-" * 80)
    logger.info("RECOMMENDATIONS")
    logger.info("-" * 80)

    logger.info("\n1. Early Discovery (>100 compounds)")
    logger.info("   → Use Fast Screening workflow")
    logger.info("   → MM/GBSA is sufficient for ranking")

    logger.info("\n2. Lead Optimization (10-20 compounds)")
    logger.info("   → Use Balanced workflow")
    logger.info("   → MM/PBSA for accurate ranking")

    logger.info("\n3. Final Validation (top 3-5)")
    logger.info("   → Use High Accuracy workflow")
    logger.info("   → Add per-residue decomposition")
    logger.info("   → Consider experimental validation")


async def main():
    """Run all integration examples"""
    try:
        # Example 1: Complete workflow
        await complete_workflow_example()

        # Example 2: Batch screening
        await batch_screening_example()

        # Example 3: Hotspot analysis
        await hotspot_analysis_example()

        # Example 4: Workflow comparison
        await workflow_comparison()

        logger.info("\n" + "=" * 80)
        logger.info("All integration examples complete!")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Error in integration examples: {e}")


if __name__ == "__main__":
    asyncio.run(main())
