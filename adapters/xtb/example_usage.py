"""
Example usage of the xTB adapter for PharmForge

This script demonstrates various use cases for quantum chemistry calculations
using the xTB adapter.
"""

import asyncio
import logging
from adapters.xtb.adapter import xTBAdapter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


async def example_basic():
    """Basic xTB calculation for a simple molecule."""
    print("\n" + "="*70)
    print("Example 1: Basic Calculation")
    print("="*70)

    xtb = xTBAdapter()

    # Ethanol
    smiles = "CCO"
    print(f"\nCalculating properties for ethanol: {smiles}")

    result = await xtb.execute(smiles)

    if result.success:
        data = result.data
        print("\nResults:")
        print(f"  Energy: {data['energy_kcal_mol']:.2f} kcal/mol")
        print(f"  HOMO energy: {data['homo_energy_eV']:.2f} eV")
        print(f"  LUMO energy: {data['lumo_energy_eV']:.2f} eV")
        print(f"  HOMO-LUMO gap: {data['homo_lumo_gap_eV']:.2f} eV")
        print(f"  Dipole moment: {data['dipole_moment_debye']:.2f} Debye")
        print(f"  Converged: {data['converged']}")
        print(f"  Cache hit: {result.cache_hit}")
    else:
        print(f"Error: {result.error}")


async def example_method_comparison():
    """Compare different GFN methods."""
    print("\n" + "="*70)
    print("Example 2: GFN Method Comparison")
    print("="*70)

    # Aspirin
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    print(f"\nComparing GFN methods for aspirin: {smiles}")

    methods = ["GFN0-xTB", "GFN1-xTB", "GFN2-xTB"]
    results = {}

    for method in methods:
        xtb = xTBAdapter(config={"method": method})
        result = await xtb.execute(smiles)

        if result.success:
            results[method] = result.data
            print(f"\n{method}:")
            print(f"  Energy: {result.data['energy_kcal_mol']:.2f} kcal/mol")
            print(f"  HOMO-LUMO gap: {result.data['homo_lumo_gap_eV']:.2f} eV")

    # Compare results
    if len(results) == 3:
        print("\nComparison:")
        energies = [results[m]['energy_kcal_mol'] for m in methods]
        gaps = [results[m]['homo_lumo_gap_eV'] for m in methods]
        print(f"  Energy range: {max(energies) - min(energies):.2f} kcal/mol")
        print(f"  Gap range: {max(gaps) - min(gaps):.2f} eV")


async def example_solvent_effects():
    """Calculate solvent effects using GBSA."""
    print("\n" + "="*70)
    print("Example 3: Solvent Effects")
    print("="*70)

    # Acetic acid
    smiles = "CC(=O)O"
    print(f"\nCalculating solvent effects for acetic acid: {smiles}")

    xtb = xTBAdapter()

    # Gas phase
    result_gas = await xtb.execute(smiles, solvent=None)

    # Water
    result_water = await xtb.execute(smiles, solvent="water")

    # DMSO
    result_dmso = await xtb.execute(smiles, solvent="dmso")

    if result_gas.success and result_water.success and result_dmso.success:
        E_gas = result_gas.data['energy_kcal_mol']
        E_water = result_water.data['energy_kcal_mol']
        E_dmso = result_dmso.data['energy_kcal_mol']

        print("\nEnergies:")
        print(f"  Gas phase: {E_gas:.2f} kcal/mol")
        print(f"  Water: {E_water:.2f} kcal/mol")
        print(f"  DMSO: {E_dmso:.2f} kcal/mol")

        print("\nSolvation free energies:")
        print(f"  ΔG_solv (water): {E_water - E_gas:.2f} kcal/mol")
        print(f"  ΔG_solv (DMSO): {E_dmso - E_gas:.2f} kcal/mol")


async def example_charged_species():
    """Calculate properties of charged molecules."""
    print("\n" + "="*70)
    print("Example 4: Charged Species")
    print("="*70)

    xtb = xTBAdapter()

    # Acetate anion
    smiles_anion = "CC(=O)[O-]"
    print(f"\nAcetate anion: {smiles_anion}")
    result_anion = await xtb.execute(smiles_anion, charge=-1, solvent="water")

    if result_anion.success:
        print("  Results:")
        print(f"    Energy: {result_anion.data['energy_kcal_mol']:.2f} kcal/mol")
        print(f"    HOMO-LUMO gap: {result_anion.data['homo_lumo_gap_eV']:.2f} eV")

    # Ammonium cation
    smiles_cation = "C[NH3+]"
    print(f"\nMethylammonium cation: {smiles_cation}")
    result_cation = await xtb.execute(smiles_cation, charge=+1, solvent="water")

    if result_cation.success:
        print("  Results:")
        print(f"    Energy: {result_cation.data['energy_kcal_mol']:.2f} kcal/mol")
        print(f"    HOMO-LUMO gap: {result_cation.data['homo_lumo_gap_eV']:.2f} eV")


async def example_batch_screening():
    """Screen multiple compounds for HOMO-LUMO gaps."""
    print("\n" + "="*70)
    print("Example 5: Batch Screening - HOMO-LUMO Gaps")
    print("="*70)

    # Common drug molecules
    molecules = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(cc1)O",
        "Benzene": "c1ccccc1",
    }

    xtb = xTBAdapter(config={"method": "GFN1-xTB"})  # Use faster method

    print("\nScreening molecules for HOMO-LUMO gaps...")

    results = []
    for name, smiles in molecules.items():
        result = await xtb.execute(smiles)
        if result.success:
            gap = result.data['homo_lumo_gap_eV']
            results.append((name, gap))
            print(f"  {name}: {gap:.2f} eV")

    # Sort by gap
    results.sort(key=lambda x: x[1])

    print("\nRanking (smallest to largest gap):")
    for i, (name, gap) in enumerate(results, 1):
        print(f"  {i}. {name}: {gap:.2f} eV")


async def example_optimization_comparison():
    """Compare single-point vs. optimized calculations."""
    print("\n" + "="*70)
    print("Example 6: Optimization Comparison")
    print("="*70)

    # Propanol (flexible molecule)
    smiles = "CCCO"
    print(f"\nComparing single-point vs. optimization for: {smiles}")

    xtb = xTBAdapter()

    # Single-point calculation
    print("\nSingle-point calculation...")
    result_sp = await xtb.execute(smiles, optimize=False)

    # Optimized calculation
    print("Geometry optimization...")
    result_opt = await xtb.execute(smiles, optimize=True)

    if result_sp.success and result_opt.success:
        E_sp = result_sp.data['energy_kcal_mol']
        E_opt = result_opt.data['energy_kcal_mol']

        print("\nResults:")
        print(f"  Single-point energy: {E_sp:.2f} kcal/mol")
        print(f"  Optimized energy: {E_opt:.2f} kcal/mol")
        print(f"  Energy difference: {E_sp - E_opt:.2f} kcal/mol")
        print(f"  Optimization converged: {result_opt.data['converged']}")


async def example_drug_molecule():
    """Comprehensive analysis of a drug molecule."""
    print("\n" + "="*70)
    print("Example 7: Comprehensive Drug Analysis")
    print("="*70)

    # Imatinib (Gleevec - cancer drug)
    smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
    print(f"\nAnalyzing Imatinib (Gleevec)")
    print(f"SMILES: {smiles}")

    xtb = xTBAdapter(config={"method": "GFN2-xTB"})

    # Calculate in water (physiological conditions)
    result = await xtb.execute(smiles, solvent="water", optimize=True)

    if result.success:
        data = result.data
        print("\nQuantum Chemical Properties:")
        print(f"  Total energy: {data['energy_kcal_mol']:.2f} kcal/mol")
        print(f"  HOMO energy: {data['homo_energy_eV']:.2f} eV")
        print(f"  LUMO energy: {data['lumo_energy_eV']:.2f} eV")
        print(f"  HOMO-LUMO gap: {data['homo_lumo_gap_eV']:.2f} eV")
        print(f"  Dipole moment: {data['dipole_moment_debye']:.2f} Debye")
        print(f"  Number of atoms: {data['n_atoms']}")
        print(f"  Optimization converged: {data['converged']}")

        print("\nInterpretation:")
        gap = data['homo_lumo_gap_eV']
        if gap > 5.0:
            print("  - Large HOMO-LUMO gap suggests chemical stability")
        elif gap > 3.0:
            print("  - Moderate HOMO-LUMO gap (typical for drug molecules)")
        else:
            print("  - Small HOMO-LUMO gap suggests high reactivity")

        dipole = data['dipole_moment_debye']
        if dipole > 5.0:
            print("  - High dipole moment suggests good water solubility")
        elif dipole > 2.0:
            print("  - Moderate polarity")
        else:
            print("  - Low polarity (lipophilic)")


async def main():
    """Run all examples."""
    print("\n" + "="*70)
    print("xTB Adapter - Example Usage")
    print("="*70)

    examples = [
        ("Basic Calculation", example_basic),
        ("GFN Method Comparison", example_method_comparison),
        ("Solvent Effects", example_solvent_effects),
        ("Charged Species", example_charged_species),
        ("Batch Screening", example_batch_screening),
        ("Optimization Comparison", example_optimization_comparison),
        ("Drug Molecule Analysis", example_drug_molecule),
    ]

    print("\nAvailable examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")

    print("\nRunning all examples...\n")

    for name, example_func in examples:
        try:
            await example_func()
        except Exception as e:
            logger.error(f"Error in {name}: {e}", exc_info=True)

    print("\n" + "="*70)
    print("All examples completed!")
    print("="*70)


if __name__ == "__main__":
    # Check if xtb is available
    try:
        import xtb
        print("✓ xtb-python is installed")
    except ImportError:
        print("\n" + "!"*70)
        print("WARNING: xtb-python is not installed!")
        print("Install with: conda install -c conda-forge xtb-python")
        print("!"*70 + "\n")
        exit(1)

    # Run examples
    asyncio.run(main())
