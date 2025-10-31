"""
Example Workflows for De Novo Drug Design Adapters

Demonstrates common use cases:
1. Basic de novo generation with REINVENT
2. Property-optimized generation with MolGAN
3. Target-guided design with TargetNet
4. Multi-objective optimization
5. Integration with retrosynthesis and docking
"""

import asyncio
import logging
from typing import List, Dict, Any

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# ============================================================================
# Example 1: Basic De Novo Generation with REINVENT
# ============================================================================

async def example_reinvent_basic_generation():
    """
    Generate novel molecules with REINVENT targeting specific properties.
    """
    from adapters.reinvent import REINVENTAdapter

    print("\n" + "=" * 60)
    print("Example 1: Basic REINVENT Generation")
    print("=" * 60 + "\n")

    # Initialize adapter
    adapter = REINVENTAdapter(config={
        'mode': 'generate',
        'use_local': True
    })

    # Define target properties
    target_properties = {
        'molecular_weight': 400.0,  # Target MW around 400
        'logp': 3.0,                # Target LogP around 3
        'qed': 0.7                  # Target QED > 0.7
    }

    # Define constraints
    constraints = {
        'mw_range': (300, 500),
        'logp_range': (1, 5),
        'max_atoms': 40
    }

    # Generate molecules
    result = await adapter.execute(
        input_data={},
        num_molecules=50,
        target_properties=target_properties,
        constraints=constraints,
        optimization_steps=100,
        diversity_filter=True
    )

    if result.success:
        molecules = result.data['molecules']
        stats = result.data['statistics']

        print(f"Generated {len(molecules)} molecules")
        print(f"Validity: {stats['valid_ratio']:.1%}")
        print(f"Mean score: {stats['mean_score']:.3f}")

        # Show top 5
        print("\nTop 5 molecules:")
        for i, mol in enumerate(molecules[:5], 1):
            print(f"{i}. {mol['smiles'][:60]}... (score: {mol['score']:.3f})")
    else:
        print(f"Generation failed: {result.error}")

    return result


# ============================================================================
# Example 2: Property-Optimized Generation with MolGAN
# ============================================================================

async def example_molgan_property_optimization():
    """
    Generate drug-like molecules optimized for specific properties.
    """
    from adapters.molgan import MolGANAdapter

    print("\n" + "=" * 60)
    print("Example 2: MolGAN Property Optimization")
    print("=" * 60 + "\n")

    # Initialize adapter
    adapter = MolGANAdapter(config={
        'temperature': 0.8,
        'max_atoms': 35
    })

    # Target CNS drug-like properties
    target_properties = {
        'molecular_weight': 350.0,
        'logp': 2.5,
        'tpsa': 60.0,  # CNS penetration
        'qed': 0.75
    }

    # Generate molecules
    result = await adapter.execute(
        input_data=None,
        num_molecules=100,
        target_properties=target_properties,
        temperature=0.8,
        validate=True
    )

    if result.success:
        molecules = result.data['molecules']
        stats = result.data['statistics']

        print(f"Generated {len(molecules)} molecules")
        print(f"Validity: {stats['validity_ratio']:.1%}")
        print(f"Uniqueness: {stats['uniqueness_ratio']:.1%}")
        print(f"Mean QED: {stats['mean_qed']:.3f}")

        # Analyze drug-likeness
        druglike = sum(1 for m in molecules if m.get('properties', {}).get('qed', 0) > 0.5)
        print(f"Drug-like (QED > 0.5): {druglike}/{len(molecules)} ({druglike/len(molecules):.1%})")

    else:
        print(f"Generation failed: {result.error}")

    return result


# ============================================================================
# Example 3: Target-Guided Design with TargetNet
# ============================================================================

async def example_targetnet_guided_design():
    """
    Generate molecules and predict their protein targets.
    """
    from adapters.targetnet import TargetNetAdapter
    from adapters.denovo import DeNovoAdapter

    print("\n" + "=" * 60)
    print("Example 3: Target-Guided Design")
    print("=" * 60 + "\n")

    # Initialize adapters
    denovo_adapter = DeNovoAdapter(config={'model_type': 'rnn'})
    target_adapter = TargetNetAdapter()

    # Generate initial molecules
    print("Generating molecules...")
    gen_result = await denovo_adapter.execute(
        input_data=None,
        num_molecules=20,
        target_profile={'qed': 0.7}
    )

    if not gen_result.success:
        print(f"Generation failed: {gen_result.error}")
        return

    molecules = gen_result.data['molecules']
    print(f"Generated {len(molecules)} molecules")

    # Predict targets for each molecule
    print("\nPredicting targets...")
    for i, mol_data in enumerate(molecules[:5], 1):  # Top 5
        smiles = mol_data['smiles']

        target_result = await target_adapter.execute(
            input_data=smiles,
            confidence_threshold=0.6,
            max_targets=5
        )

        if target_result.success:
            targets = target_result.data['primary_targets']
            print(f"\n{i}. {smiles[:50]}...")
            print(f"   Targets: {len(targets)}")
            for target in targets[:3]:
                print(f"   - {target['target_name']} ({target['confidence']:.2f})")
        else:
            print(f"\n{i}. Target prediction failed")

    return gen_result, target_result


# ============================================================================
# Example 4: Multi-Objective Optimization
# ============================================================================

async def example_multi_objective_optimization():
    """
    Optimize molecules for multiple objectives simultaneously.
    """
    from adapters.denovo import DeNovoAdapter

    print("\n" + "=" * 60)
    print("Example 4: Multi-Objective Optimization")
    print("=" * 60 + "\n")

    # Initialize adapter
    adapter = DeNovoAdapter(config={
        'model_type': 'vae',
        'optimize_properties': True
    })

    # Define multiple objectives
    target_profile = {
        'molecular_weight': 400.0,
        'logp': 2.5,
        'tpsa': 70.0,
        'qed': 0.75
    }

    objectives = ['qed', 'molecular_weight', 'logp', 'tpsa']

    # Generate with optimization
    result = await adapter.execute(
        input_data=None,
        num_molecules=50,
        target_profile=target_profile,
        optimization_objectives=objectives,
        diversity_threshold=0.7
    )

    if result.success:
        molecules = result.data['molecules']
        metrics = result.data['metrics']

        print(f"Generated {len(molecules)} molecules")
        print(f"Validity: {metrics['validity']:.1%}")
        print(f"Diversity: {metrics['diversity']:.1%}")

        # Analyze Pareto front
        print("\nTop 5 molecules (multi-objective score):")
        for i, mol in enumerate(molecules[:5], 1):
            props = mol.get('properties', {})
            print(f"{i}. Score: {mol.get('score', 0):.3f}")
            print(f"   MW: {props.get('molecular_weight', 0):.1f}, "
                  f"LogP: {props.get('logp', 0):.2f}, "
                  f"QED: {props.get('qed', 0):.3f}")

    else:
        print(f"Optimization failed: {result.error}")

    return result


# ============================================================================
# Example 5: Integration with Retrosynthesis
# ============================================================================

async def example_integration_retrosynthesis():
    """
    Generate molecules and evaluate their synthetic accessibility.
    """
    from adapters.denovo import DeNovoAdapter
    from adapters.llm_retrosynthesis import LLMRetrosynthesisAdapter

    print("\n" + "=" * 60)
    print("Example 5: Integration with Retrosynthesis")
    print("=" * 60 + "\n")

    # Initialize adapters
    denovo = DeNovoAdapter(config={'model_type': 'fragment'})

    try:
        retro = LLMRetrosynthesisAdapter()
    except Exception as e:
        print(f"Retrosynthesis adapter not available: {e}")
        print("Using SA score estimation instead")
        retro = None

    # Generate molecules
    print("Generating molecules...")
    gen_result = await denovo.execute(
        input_data=None,
        num_molecules=10,
        target_profile={'qed': 0.7, 'molecular_weight': 350}
    )

    if not gen_result.success:
        print(f"Generation failed: {gen_result.error}")
        return

    molecules = gen_result.data['molecules']
    print(f"Generated {len(molecules)} molecules")

    # Evaluate synthetic accessibility
    print("\nEvaluating synthetic accessibility...")

    if retro:
        # Use retrosynthesis analysis
        for i, mol_data in enumerate(molecules[:5], 1):
            smiles = mol_data['smiles']

            retro_result = await retro.execute(
                input_data=smiles,
                max_iterations=3
            )

            if retro_result.success:
                analysis = retro_result.data
                print(f"\n{i}. {smiles[:50]}...")
                print(f"   Synthetic routes found: {len(analysis.get('routes', []))}")
                print(f"   Convergent synthesis: {analysis.get('convergent', False)}")
            else:
                print(f"\n{i}. Retrosynthesis analysis failed")
    else:
        # Use SA score
        from rdkit import Chem
        try:
            from adapters.custom.synthesis.sascore import calculateScore

            for i, mol_data in enumerate(molecules[:5], 1):
                smiles = mol_data['smiles']
                mol = Chem.MolFromSmiles(smiles)

                if mol:
                    sa_score = calculateScore(mol)
                    print(f"\n{i}. {smiles[:50]}...")
                    print(f"   SA Score: {sa_score:.2f} "
                          f"({'Easy' if sa_score <= 3 else 'Moderate' if sa_score <= 6 else 'Difficult'})")
        except ImportError:
            print("SA score calculation not available")

    return gen_result


# ============================================================================
# Example 6: Scaffold Hopping
# ============================================================================

async def example_scaffold_hopping():
    """
    Generate molecules with modified scaffolds.
    """
    from adapters.reinvent import REINVENTAdapter

    print("\n" + "=" * 60)
    print("Example 6: Scaffold Hopping")
    print("=" * 60 + "\n")

    # Initialize adapter
    adapter = REINVENTAdapter(config={'mode': 'optimize'})

    # Starting molecule (e.g., aspirin scaffold)
    seed_smiles = "CC(=O)Oc1ccccc1C(=O)O"

    print(f"Seed molecule: {seed_smiles}")

    # Optimize/modify
    result = await adapter.execute(
        input_data=seed_smiles,
        num_molecules=20,
        target_properties={
            'molecular_weight': 300.0,
            'logp': 2.5,
            'qed': 0.7
        }
    )

    if result.success:
        molecules = result.data['molecules']

        print(f"\nGenerated {len(molecules)} scaffold variants")

        # Show top variants
        print("\nTop 5 variants:")
        for i, mol in enumerate(molecules[:5], 1):
            is_original = mol.get('is_original', False)
            marker = " [ORIGINAL]" if is_original else ""
            print(f"{i}. {mol['smiles'][:60]}...{marker}")
            print(f"   Score: {mol['score']:.3f}")

    else:
        print(f"Scaffold hopping failed: {result.error}")

    return result


# ============================================================================
# Example 7: Comprehensive Pipeline
# ============================================================================

async def example_comprehensive_pipeline():
    """
    Complete pipeline: Generate -> Filter -> Predict Targets -> Evaluate Synthesis
    """
    from adapters.denovo import DeNovoAdapter, GenerationMetrics
    from adapters.targetnet import TargetNetAdapter

    print("\n" + "=" * 60)
    print("Example 7: Comprehensive Drug Design Pipeline")
    print("=" * 60 + "\n")

    # Step 1: Generate molecules
    print("Step 1: Generating molecules...")
    denovo = DeNovoAdapter(config={'model_type': 'rnn'})

    gen_result = await denovo.execute(
        input_data=None,
        num_molecules=50,
        target_profile={
            'molecular_weight': 400.0,
            'logp': 2.5,
            'qed': 0.7
        },
        diversity_threshold=0.7
    )

    if not gen_result.success:
        print(f"Generation failed: {gen_result.error}")
        return

    molecules = gen_result.data['molecules']
    print(f"Generated {len(molecules)} molecules")

    # Step 2: Calculate comprehensive metrics
    print("\nStep 2: Calculating generation metrics...")
    metrics_calc = GenerationMetrics(calculate_synthetic_accessibility=True)

    smiles_list = [m['smiles'] for m in molecules]
    metrics = metrics_calc.calculate_all_metrics(smiles_list, verbose=True)

    # Step 3: Filter by quality
    print("\nStep 3: Filtering by quality...")
    filtered_molecules = [
        m for m in molecules
        if m.get('properties', {}).get('qed', 0) > 0.6
    ]
    print(f"Filtered to {len(filtered_molecules)} high-quality molecules")

    # Step 4: Predict targets
    print("\nStep 4: Predicting protein targets...")
    target_adapter = TargetNetAdapter()

    molecules_with_targets = []
    for mol_data in filtered_molecules[:10]:  # Top 10
        target_result = await target_adapter.execute(
            input_data=mol_data['smiles'],
            confidence_threshold=0.6,
            max_targets=5
        )

        if target_result.success:
            mol_data['targets'] = target_result.data['primary_targets']
            molecules_with_targets.append(mol_data)

    print(f"Target prediction completed for {len(molecules_with_targets)} molecules")

    # Step 5: Summary
    print("\n" + "=" * 60)
    print("PIPELINE SUMMARY")
    print("=" * 60)
    print(f"Generated: {len(molecules)} molecules")
    print(f"High quality: {len(filtered_molecules)} molecules")
    print(f"With targets: {len(molecules_with_targets)} molecules")

    # Top candidates
    print("\nTop 3 candidates:")
    for i, mol in enumerate(molecules_with_targets[:3], 1):
        print(f"\n{i}. {mol['smiles'][:50]}...")
        print(f"   Score: {mol.get('score', 0):.3f}")
        print(f"   QED: {mol.get('properties', {}).get('qed', 0):.3f}")
        print(f"   Targets: {len(mol.get('targets', []))}")

    return molecules_with_targets


# ============================================================================
# Main execution
# ============================================================================

async def run_all_examples():
    """Run all example workflows."""
    examples = [
        ("REINVENT Basic Generation", example_reinvent_basic_generation),
        ("MolGAN Property Optimization", example_molgan_property_optimization),
        ("Target-Guided Design", example_targetnet_guided_design),
        ("Multi-Objective Optimization", example_multi_objective_optimization),
        ("Integration with Retrosynthesis", example_integration_retrosynthesis),
        ("Scaffold Hopping", example_scaffold_hopping),
        ("Comprehensive Pipeline", example_comprehensive_pipeline)
    ]

    results = {}

    for name, example_func in examples:
        print(f"\n{'=' * 70}")
        print(f"Running: {name}")
        print(f"{'=' * 70}")

        try:
            result = await example_func()
            results[name] = result
            print(f"\n{'=' * 70}")
            print(f"Completed: {name}")
            print(f"{'=' * 70}\n")
        except Exception as e:
            print(f"\nExample failed: {e}")
            import traceback
            traceback.print_exc()
            results[name] = None

    return results


def main():
    """Main entry point."""
    print("\n" + "=" * 70)
    print("DE NOVO DRUG DESIGN ADAPTER EXAMPLES")
    print("=" * 70 + "\n")

    print("This script demonstrates various de novo drug design workflows.")
    print("Choose an example to run:\n")

    print("1. REINVENT Basic Generation")
    print("2. MolGAN Property Optimization")
    print("3. Target-Guided Design")
    print("4. Multi-Objective Optimization")
    print("5. Integration with Retrosynthesis")
    print("6. Scaffold Hopping")
    print("7. Comprehensive Pipeline")
    print("8. Run all examples")
    print("0. Exit\n")

    choice = input("Enter choice (0-8): ").strip()

    if choice == "0":
        print("Exiting...")
        return

    examples_map = {
        "1": example_reinvent_basic_generation,
        "2": example_molgan_property_optimization,
        "3": example_targetnet_guided_design,
        "4": example_multi_objective_optimization,
        "5": example_integration_retrosynthesis,
        "6": example_scaffold_hopping,
        "7": example_comprehensive_pipeline,
        "8": run_all_examples
    }

    if choice in examples_map:
        asyncio.run(examples_map[choice]())
    else:
        print("Invalid choice")


if __name__ == "__main__":
    main()
