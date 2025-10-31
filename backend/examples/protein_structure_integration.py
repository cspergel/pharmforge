"""
Integration Examples: Protein Structure Databases + OpenMM

Demonstrates how to use protein structure adapters with molecular dynamics
simulations for drug discovery workflows.

Examples:
1. Retrieve protein structure and run MD simulation
2. Compare predicted vs experimental structures
3. Prepare protein-ligand complex for docking
4. Quality assessment workflow
"""

import asyncio
from pathlib import Path

# Import protein structure adapters
from adapters.alphafold.adapter import AlphaFoldAdapter
from adapters.rcsb_pdb.adapter import RCSBPDBAdapter
from adapters.pdb_redo.adapter import PDBRedoAdapter
from adapters.swissmodel.adapter import SwissModelAdapter

# Import OpenMM adapter
from adapters.openmm.adapter import OpenMMAdapter


async def example1_alphafold_md_workflow():
    """
    Example 1: Retrieve AlphaFold structure and run MD simulation

    Workflow:
    1. Download AlphaFold predicted structure
    2. Check confidence scores
    3. Run energy minimization with OpenMM
    4. Assess stability
    """
    print("=" * 80)
    print("Example 1: AlphaFold Structure + OpenMM MD Workflow")
    print("=" * 80)

    # Step 1: Get AlphaFold structure
    alphafold = AlphaFoldAdapter(config={
        "download_pdb": True,
        "download_pae": True
    })

    uniprot_id = "P12345"  # Example: Aspartate aminotransferase
    print(f"\n1. Fetching AlphaFold structure for {uniprot_id}...")

    af_result = await alphafold.execute(uniprot_id)

    if not af_result.success:
        print(f"   ERROR: {af_result.error}")
        return

    print(f"   ✓ Structure retrieved")

    # Step 2: Check confidence scores
    plddt = af_result.data.get("plddt_scores", {})
    print(f"\n2. Confidence Assessment:")
    print(f"   Mean pLDDT: {plddt.get('mean', 'N/A')}")
    print(f"   High confidence regions: {plddt.get('high_confidence_pct', 'N/A')}%")
    print(f"   Low confidence regions: {plddt.get('low_confidence_pct', 'N/A')}%")

    # Check if structure is suitable for MD
    mean_plddt = plddt.get('mean', 0)
    if mean_plddt < 70:
        print(f"   ⚠ WARNING: Low average confidence (pLDDT={mean_plddt})")
        print(f"   MD results may be unreliable for low-confidence regions")
    else:
        print(f"   ✓ Good confidence, suitable for MD simulation")

    # Step 3: Get PDB structure for OpenMM
    pdb_file = af_result.data.get("file_paths", {}).get("pdb")

    if not pdb_file:
        print("\n   ERROR: PDB file not available")
        return

    print(f"\n3. Structure file: {pdb_file}")

    # Note: OpenMM adapter works with SMILES, not PDB files
    # For protein MD, you would use OpenMM directly or create a protein-specific adapter
    print("\n4. OpenMM Protein Simulation:")
    print("   Note: OpenMM adapter currently supports small molecules.")
    print("   For protein MD, use OpenMM directly with the PDB file:")
    print(f"   - Load structure: openmm.app.PDBFile('{pdb_file}')")
    print(f"   - Add force field: amber14/protein.ff14SB.xml")
    print(f"   - Run minimization and MD")

    print("\n" + "=" * 80)


async def example2_compare_structure_sources():
    """
    Example 2: Compare structures from different sources

    Compare AlphaFold prediction vs experimental structure
    """
    print("=" * 80)
    print("Example 2: Compare Structure Sources")
    print("=" * 80)

    # Insulin: Good example with both experimental and predicted structures
    pdb_id = "1ABC"  # Experimental structure
    uniprot_id = "P01308"  # Insulin UniProt ID

    # Initialize adapters
    rcsb = RCSBPDBAdapter(config={"download_pdb": True})
    alphafold = AlphaFoldAdapter(config={"download_pdb": True})

    print(f"\n1. Fetching experimental structure (PDB {pdb_id})...")
    exp_result = await rcsb.execute(pdb_id)

    print(f"\n2. Fetching predicted structure (AlphaFold {uniprot_id})...")
    pred_result = await alphafold.execute(uniprot_id)

    # Compare results
    print("\n3. Comparison:")
    print("\n   Experimental Structure (RCSB PDB):")
    if exp_result.success:
        exp_metrics = exp_result.data.get("quality_metrics", {})
        print(f"   - Method: {exp_result.data.get('experimental_method', 'N/A')}")
        print(f"   - Resolution: {exp_metrics.get('resolution', 'N/A')} Å")
        print(f"   - R-value: {exp_metrics.get('r_value', 'N/A')}")
        print(f"   - Ligands: {len(exp_result.data.get('ligands', []))}")
    else:
        print(f"   ERROR: {exp_result.error}")

    print("\n   Predicted Structure (AlphaFold):")
    if pred_result.success:
        plddt = pred_result.data.get("plddt_scores", {})
        print(f"   - Model version: {pred_result.data.get('model_version', 'N/A')}")
        print(f"   - Mean pLDDT: {plddt.get('mean', 'N/A')}")
        print(f"   - High confidence: {plddt.get('high_confidence_pct', 'N/A')}%")
    else:
        print(f"   ERROR: {pred_result.error}")

    print("\n   Use Cases:")
    print("   - Experimental: Best for drug docking (crystallized with ligands)")
    print("   - Predicted: Good when no experimental structure available")
    print("   - Combine: Use AlphaFold for missing regions, experimental for binding site")

    print("\n" + "=" * 80)


async def example3_structure_refinement_workflow():
    """
    Example 3: Structure refinement workflow using PDB-REDO

    Compare original PDB with re-refined structure
    """
    print("=" * 80)
    print("Example 3: Structure Refinement with PDB-REDO")
    print("=" * 80)

    pdb_id = "1ABC"

    # Get original structure
    rcsb = RCSBPDBAdapter(config={"download_pdb": True})
    print(f"\n1. Fetching original PDB structure ({pdb_id})...")
    orig_result = await rcsb.execute(pdb_id)

    # Get re-refined structure
    pdb_redo = PDBRedoAdapter(config={
        "download_pdb": True,
        "include_validation": True
    })
    print(f"\n2. Fetching PDB-REDO re-refined structure...")
    redo_result = await pdb_redo.execute(pdb_id)

    # Compare quality metrics
    print("\n3. Quality Comparison:")

    if orig_result.success:
        orig_metrics = orig_result.data.get("quality_metrics", {})
        print("\n   Original PDB:")
        print(f"   - R-work: {orig_metrics.get('r_value', 'N/A')}")
        print(f"   - R-free: {orig_metrics.get('r_free', 'N/A')}")
        print(f"   - Resolution: {orig_metrics.get('resolution', 'N/A')} Å")

    if redo_result.success:
        redo_metrics = redo_result.data.get("quality_metrics", {})
        improvements = redo_result.data.get("improvements", {})

        print("\n   PDB-REDO Refined:")
        if "refined" in redo_metrics:
            print(f"   - R-work: {redo_metrics['refined'].get('r_work', 'N/A')}")
            print(f"   - R-free: {redo_metrics['refined'].get('r_free', 'N/A')}")

        if improvements:
            print("\n   Improvements:")
            print(f"   - R-work: {improvements.get('r_work_improvement_pct', 'N/A')}% better")
            print(f"   - R-free: {improvements.get('r_free_improvement_pct', 'N/A')}% better")

        print("\n   Recommendation:")
        print("   Use PDB-REDO structure for better geometry and reduced errors")
    else:
        print(f"\n   PDB-REDO not available: {redo_result.error}")

    print("\n" + "=" * 80)


async def example4_homology_modeling_workflow():
    """
    Example 4: Homology modeling with SWISS-MODEL

    Get homology models when experimental structure unavailable
    """
    print("=" * 80)
    print("Example 4: Homology Modeling with SWISS-MODEL")
    print("=" * 80)

    uniprot_id = "P12345"

    # Try to get experimental structure first
    rcsb = RCSBPDBAdapter()
    print(f"\n1. Searching for experimental structures...")
    search_result = await rcsb.search_structures(
        query=f"uniprot:{uniprot_id}",
        max_results=1
    )

    has_experimental = search_result.success and len(search_result.data.get("pdb_ids", [])) > 0

    if has_experimental:
        print(f"   ✓ Experimental structure available: {search_result.data['pdb_ids'][0]}")
    else:
        print(f"   ✗ No experimental structure found")

    # Get SWISS-MODEL homology models
    swissmodel = SwissModelAdapter(config={
        "download_pdb": True,
        "min_qmean": -3.0,
        "min_gmqe": 0.4
    })

    print(f"\n2. Fetching SWISS-MODEL homology models...")
    sm_result = await swissmodel.execute(uniprot_id)

    if sm_result.success:
        models = sm_result.data.get("models", [])
        print(f"   ✓ Found {len(models)} models")

        best = sm_result.data.get("best_model", {})
        quality = sm_result.data.get("quality_metrics", {})

        print(f"\n3. Best Model Quality:")
        print(f"   - QMEAN: {quality.get('qmean', 'N/A')} (higher is better, max 0)")
        print(f"   - GMQE: {quality.get('gmqe', 'N/A')} (0-1, higher is better)")
        print(f"   - Sequence identity: {quality.get('sequence_identity', 'N/A')}%")
        print(f"   - Coverage: {quality.get('coverage', 'N/A')}%")

        assessment = quality.get("quality_assessment", [])
        if assessment:
            print(f"\n   Assessment: {', '.join(assessment)}")

        template = sm_result.data.get("template_info", {})
        if template:
            print(f"\n4. Template Information:")
            print(f"   - Template PDB: {template.get('pdb_id', 'N/A')}")
            print(f"   - Method: {template.get('method', 'N/A')}")

        print("\n5. Recommendation:")
        qmean = quality.get("qmean", -10)
        if qmean >= -1.5:
            print("   ✓ Very high quality model - suitable for structure-based design")
        elif qmean >= -2.5:
            print("   ✓ Good quality model - suitable for most applications")
        elif qmean >= -3.5:
            print("   ⚠ Medium quality - verify key regions before use")
        else:
            print("   ⚠ Low quality - use with caution, consider AlphaFold instead")
    else:
        print(f"   ERROR: {sm_result.error}")

    print("\n" + "=" * 80)


async def example5_quality_assessment_workflow():
    """
    Example 5: Comprehensive quality assessment workflow

    Compare quality metrics across all sources
    """
    print("=" * 80)
    print("Example 5: Comprehensive Quality Assessment")
    print("=" * 80)

    uniprot_id = "P01308"  # Insulin
    pdb_id = "1ABC"

    print(f"\nAssessing structure quality for {uniprot_id} ({pdb_id})...")

    # Initialize all adapters
    adapters = {
        "AlphaFold": AlphaFoldAdapter(config={"download_pdb": True}),
        "RCSB PDB": RCSBPDBAdapter(config={"download_pdb": True}),
        "PDB-REDO": PDBRedoAdapter(config={"download_pdb": True}),
        "SWISS-MODEL": SwissModelAdapter(config={"download_pdb": True})
    }

    # Fetch structures
    results = {}

    print("\n1. Fetching structures from all sources...")
    results["AlphaFold"] = await adapters["AlphaFold"].execute(uniprot_id)
    results["RCSB PDB"] = await adapters["RCSB PDB"].execute(pdb_id)
    results["PDB-REDO"] = await adapters["PDB-REDO"].execute(pdb_id)
    results["SWISS-MODEL"] = await adapters["SWISS-MODEL"].execute(uniprot_id)

    # Compare quality metrics
    print("\n2. Quality Metrics Summary:")
    print("\n" + "-" * 80)

    for source, result in results.items():
        print(f"\n{source}:")
        if result.success:
            data = result.data

            if source == "AlphaFold":
                plddt = data.get("plddt_scores", {})
                print(f"  ✓ Available")
                print(f"  - Mean pLDDT: {plddt.get('mean', 'N/A')}")
                print(f"  - High confidence: {plddt.get('high_confidence_pct', 'N/A')}%")

            elif source == "RCSB PDB":
                metrics = data.get("quality_metrics", {})
                print(f"  ✓ Available")
                print(f"  - Method: {data.get('experimental_method', 'N/A')}")
                print(f"  - Resolution: {metrics.get('resolution', 'N/A')} Å")
                print(f"  - R-free: {metrics.get('r_free', 'N/A')}")

            elif source == "PDB-REDO":
                improvements = data.get("improvements", {})
                print(f"  ✓ Available")
                if improvements:
                    print(f"  - R-work improvement: {improvements.get('r_work_improvement_pct', 'N/A')}%")
                    print(f"  - R-free improvement: {improvements.get('r_free_improvement_pct', 'N/A')}%")

            elif source == "SWISS-MODEL":
                quality = data.get("quality_metrics", {})
                print(f"  ✓ Available")
                print(f"  - QMEAN: {quality.get('qmean', 'N/A')}")
                print(f"  - GMQE: {quality.get('gmqe', 'N/A')}")
                print(f"  - Models: {data.get('model_count', 'N/A')}")
        else:
            print(f"  ✗ Not available: {result.error}")

    print("\n" + "-" * 80)

    print("\n3. Recommendations:")
    print("\n   For Drug Docking:")
    if results["RCSB PDB"].success:
        print("   ✓ Use RCSB PDB experimental structure (best for binding site accuracy)")
        if results["PDB-REDO"].success:
            print("   ✓ Consider PDB-REDO for improved geometry")
    elif results["AlphaFold"].success:
        print("   ⚠ Use AlphaFold if no experimental structure (check binding site confidence)")

    print("\n   For Homology Modeling:")
    if results["SWISS-MODEL"].success:
        print("   ✓ SWISS-MODEL provides template-based models with quality scores")

    print("\n   For Molecular Dynamics:")
    if results["PDB-REDO"].success:
        print("   ✓ PDB-REDO structure (best geometry for MD)")
    elif results["AlphaFold"].success:
        print("   ✓ AlphaFold structure (good for high-confidence regions)")

    print("\n" + "=" * 80)


async def main():
    """Run all examples"""
    print("\n" + "=" * 80)
    print("PROTEIN STRUCTURE DATABASE ADAPTERS - INTEGRATION EXAMPLES")
    print("=" * 80)

    # Run examples
    examples = [
        ("AlphaFold + OpenMM Workflow", example1_alphafold_md_workflow),
        ("Compare Structure Sources", example2_compare_structure_sources),
        ("Structure Refinement (PDB-REDO)", example3_structure_refinement_workflow),
        ("Homology Modeling (SWISS-MODEL)", example4_homology_modeling_workflow),
        ("Comprehensive Quality Assessment", example5_quality_assessment_workflow)
    ]

    for i, (title, func) in enumerate(examples, 1):
        print(f"\n\nRunning Example {i}: {title}")
        print("=" * 80)
        try:
            await func()
        except Exception as e:
            print(f"\nExample failed: {e}")

        if i < len(examples):
            print("\n\nPress Enter to continue to next example...")
            input()

    print("\n\n" + "=" * 80)
    print("ALL EXAMPLES COMPLETED")
    print("=" * 80)


if __name__ == "__main__":
    asyncio.run(main())
