"""
Protein Structure Adapters - Integration Example

This example demonstrates how to use all three protein structure adapters
together for comprehensive structure-based drug discovery workflows.
"""

import asyncio
import logging
import sys
from pathlib import Path
from typing import Dict, List, Any

# Add parent directories to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from adapters.alphafold import AlphaFoldAdapter
from adapters.rcsb_pdb import RCSBPDBAdapter
from adapters.swissmodel import SwissModelAdapter

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Example 1: Complete Structure Retrieval Workflow
# ============================================================================

async def example1_complete_workflow():
    """
    Retrieve structure from all three sources and compare quality.
    """
    print("\n" + "="*80)
    print("EXAMPLE 1: Complete Structure Retrieval Workflow")
    print("="*80)

    # Target: P04637 (TP53_HUMAN - Tumor protein p53)
    uniprot_id = "P04637"
    pdb_id = "1TUP"  # Known experimental structure

    print(f"\nTarget: {uniprot_id} (TP53 Tumor Suppressor)")
    print(f"Known PDB structure: {pdb_id}")

    # Initialize all adapters
    alphafold = AlphaFoldAdapter(config={"download_pdb": True})
    rcsb = RCSBPDBAdapter(config={"download_pdb": True, "include_ligands": True})
    swissmodel = SwissModelAdapter(config={"download_pdb": True})

    # Fetch from all sources
    print("\nFetching structures from all sources...")

    af_result = await alphafold.execute(uniprot_id)
    pdb_result = await rcsb.execute(pdb_id)
    sm_result = await swissmodel.execute(uniprot_id)

    # Compare results
    print("\n" + "-"*80)
    print("STRUCTURE COMPARISON")
    print("-"*80)

    if af_result.success:
        print(f"\n✓ AlphaFold (AI-predicted):")
        print(f"    Mean pLDDT: {af_result.data['plddt_scores']['mean']}")
        print(f"    High confidence: {af_result.data['plddt_scores']['high_confidence_pct']}%")
        print(f"    Model version: {af_result.data['model_version']}")
    else:
        print(f"\n✗ AlphaFold failed: {af_result.error}")

    if pdb_result.success:
        print(f"\n✓ RCSB PDB (Experimental):")
        print(f"    Resolution: {pdb_result.data['resolution']} Å")
        print(f"    Method: {pdb_result.data['experimental_method']}")
        print(f"    R-value: {pdb_result.data['quality_metrics'].get('r_value', 'N/A')}")
        print(f"    Ligands: {len(pdb_result.data['ligands'])}")
    else:
        print(f"\n✗ RCSB PDB failed: {pdb_result.error}")

    if sm_result.success:
        print(f"\n✓ SWISS-MODEL (Homology):")
        print(f"    Models found: {sm_result.data['model_count']}")
        print(f"    Best QMEAN: {sm_result.data['quality_metrics']['qmean']}")
        print(f"    Best GMQE: {sm_result.data['quality_metrics']['gmqe']}")
    else:
        print(f"\n✗ SWISS-MODEL failed: {sm_result.error}")

    # Recommendation
    print("\n" + "-"*80)
    print("RECOMMENDATION")
    print("-"*80)

    if pdb_result.success and pdb_result.data['resolution'] and pdb_result.data['resolution'] < 2.5:
        print("\n✓ Use RCSB PDB structure (high-quality experimental data available)")
    elif af_result.success and af_result.data['plddt_scores']['mean'] > 80:
        print("\n✓ Use AlphaFold structure (high-confidence prediction)")
    elif sm_result.success:
        print("\n✓ Use SWISS-MODEL structure (homology model available)")
    else:
        print("\n⚠ No high-quality structures available")

    return af_result, pdb_result, sm_result


# ============================================================================
# Example 2: Drug Target Analysis
# ============================================================================

async def example2_drug_target_analysis():
    """
    Analyze a drug target: structure quality, ligands, binding sites.
    """
    print("\n" + "="*80)
    print("EXAMPLE 2: Drug Target Analysis")
    print("="*80)

    # Target: SARS-CoV-2 Main Protease
    pdb_id = "6LU7"
    uniprot_id = "P0DTD1"

    print(f"\nTarget: SARS-CoV-2 Main Protease (Mpro)")
    print(f"PDB: {pdb_id}, UniProt: {uniprot_id}")

    # Get experimental structure with ligands
    rcsb = RCSBPDBAdapter(config={
        "download_pdb": True,
        "include_ligands": True,
        "include_binding_sites": True
    })

    result = await rcsb.execute(pdb_id)

    if result.success:
        data = result.data

        print("\n" + "-"*80)
        print("STRUCTURE QUALITY")
        print("-"*80)
        print(f"  Resolution: {data['resolution']} Å")
        print(f"  Method: {data['experimental_method']}")
        print(f"  R-value: {data['quality_metrics'].get('r_value', 'N/A')}")
        print(f"  R-free: {data['quality_metrics'].get('r_free', 'N/A')}")

        print("\n" + "-"*80)
        print("BOUND LIGANDS")
        print("-"*80)
        if data['ligands']:
            for i, ligand in enumerate(data['ligands'], 1):
                print(f"\n  Ligand {i}:")
                print(f"    Name: {ligand.get('name', 'Unknown')}")
                print(f"    Type: {ligand.get('type', 'Unknown')}")
                if 'smiles' in ligand:
                    print(f"    SMILES: {ligand['smiles'][:80]}...")
        else:
            print("  No ligands found")

        print("\n" + "-"*80)
        print("BINDING SITES")
        print("-"*80)
        if data['binding_sites']:
            for site in data['binding_sites']:
                print(f"\n  Site {site['id']}:")
                print(f"    Residues: {len(site['residues'])}")
                print(f"    Key residues: {', '.join(site['residues'][:5])}...")
        else:
            print("  No binding sites identified")

        # Get AlphaFold prediction for comparison
        alphafold = AlphaFoldAdapter()
        af_result = await alphafold.execute(uniprot_id)

        if af_result.success:
            print("\n" + "-"*80)
            print("ALPHAFOLD COMPARISON")
            print("-"*80)
            print(f"  Mean pLDDT: {af_result.data['plddt_scores']['mean']}")
            print(f"  Regions with high confidence (>90): {af_result.data['plddt_scores']['high_confidence_pct']}%")
            print(f"  Regions with low confidence (<70): {af_result.data['plddt_scores']['low_confidence_pct']}%")

    return result


# ============================================================================
# Example 3: Batch Protein Analysis
# ============================================================================

async def example3_batch_analysis():
    """
    Analyze multiple proteins to find best structures for drug discovery.
    """
    print("\n" + "="*80)
    print("EXAMPLE 3: Batch Protein Analysis")
    print("="*80)

    # List of targets
    targets = [
        ("P04637", "TP53 - Tumor suppressor"),
        ("P01308", "Insulin"),
        ("P69905", "Hemoglobin alpha"),
        ("P00533", "EGFR - Drug target")
    ]

    print(f"\nAnalyzing {len(targets)} protein targets...")

    alphafold = AlphaFoldAdapter(config={"download_pdb": False})  # Faster without downloads

    results = []
    for uniprot_id, description in targets:
        print(f"\n  Analyzing {uniprot_id} ({description})...")
        result = await alphafold.execute(uniprot_id)
        results.append((uniprot_id, description, result))

    # Summarize results
    print("\n" + "-"*80)
    print("BATCH RESULTS SUMMARY")
    print("-"*80)

    high_quality = []
    medium_quality = []
    low_quality = []

    for uniprot_id, desc, result in results:
        if result.success:
            mean_plddt = result.data['plddt_scores']['mean']

            if mean_plddt > 80:
                high_quality.append((uniprot_id, desc, mean_plddt))
            elif mean_plddt > 60:
                medium_quality.append((uniprot_id, desc, mean_plddt))
            else:
                low_quality.append((uniprot_id, desc, mean_plddt))

    print(f"\nHigh Quality (pLDDT > 80): {len(high_quality)} proteins")
    for uid, desc, score in high_quality:
        print(f"  ✓ {uid} ({desc}): pLDDT = {score:.1f}")

    print(f"\nMedium Quality (pLDDT 60-80): {len(medium_quality)} proteins")
    for uid, desc, score in medium_quality:
        print(f"  ⚠ {uid} ({desc}): pLDDT = {score:.1f}")

    print(f"\nLow Quality (pLDDT < 60): {len(low_quality)} proteins")
    for uid, desc, score in low_quality:
        print(f"  ✗ {uid} ({desc}): pLDDT = {score:.1f}")

    return results


# ============================================================================
# Example 4: Structure Search and Discovery
# ============================================================================

async def example4_structure_search():
    """
    Search for structures related to a drug target.
    """
    print("\n" + "="*80)
    print("EXAMPLE 4: Structure Search and Discovery")
    print("="*80)

    rcsb = RCSBPDBAdapter()

    # Search for COVID-19 drug targets
    query = "SARS-CoV-2 main protease"

    print(f"\nSearching PDB for: '{query}'")
    result = await rcsb.search_structures(query, max_results=10)

    if result.success:
        pdb_ids = result.data['pdb_ids']

        print(f"\n✓ Found {len(pdb_ids)} structures")
        print("\nTop results:")
        for i, pdb_id in enumerate(pdb_ids, 1):
            print(f"  {i}. {pdb_id}")

        # Get details for first result
        if pdb_ids:
            print(f"\nRetrieving detailed information for {pdb_ids[0]}...")
            detail_result = await rcsb.execute(pdb_ids[0])

            if detail_result.success:
                data = detail_result.data
                print(f"\n  Title: {data['quality_metrics'].get('title', 'N/A')}")
                print(f"  Resolution: {data['resolution']} Å")
                print(f"  Method: {data['experimental_method']}")
                print(f"  Ligands: {len(data['ligands'])}")

    return result


# ============================================================================
# Example 5: Quality-Based Model Selection
# ============================================================================

async def example5_quality_selection():
    """
    Select best structure based on quality criteria.
    """
    print("\n" + "="*80)
    print("EXAMPLE 5: Quality-Based Model Selection")
    print("="*80)

    uniprot_id = "P01308"  # Insulin
    print(f"\nTarget: {uniprot_id} (Insulin)")

    # Get SWISS-MODEL models with quality filtering
    swissmodel = SwissModelAdapter(config={
        "download_pdb": False,
        "min_qmean": -2.5,  # Medium quality and above
        "min_gmqe": 0.4
    })

    result = await swissmodel.execute(uniprot_id)

    if result.success:
        data = result.data

        print(f"\n✓ Found {data['model_count']} models meeting quality criteria")
        print(f"  Filters: QMEAN >= -2.5, GMQE >= 0.4")

        print("\n" + "-"*80)
        print("ALL MODELS (sorted by GMQE)")
        print("-"*80)

        # Sort models by GMQE
        models = sorted(data['models'], key=lambda m: m.get('gmqe', 0), reverse=True)

        for i, model in enumerate(models[:5], 1):  # Show top 5
            print(f"\n  Model {i}:")
            print(f"    QMEAN: {model.get('qmean', 'N/A')}")
            print(f"    GMQE: {model.get('gmqe', 'N/A')}")
            print(f"    Sequence Identity: {model.get('identity', 'N/A')}%")
            print(f"    Coverage: {model.get('coverage', 'N/A')}%")

        print("\n" + "-"*80)
        print("SELECTED MODEL (Best GMQE)")
        print("-"*80)
        best = data['quality_metrics']
        print(f"  QMEAN: {best['qmean']}")
        print(f"  GMQE: {best['gmqe']}")
        print(f"  Quality: {', '.join(best['quality_assessment'])}")

    return result


# ============================================================================
# Example 6: Comprehensive Drug Target Report
# ============================================================================

async def example6_drug_target_report():
    """
    Generate a comprehensive report for a drug target combining all sources.
    """
    print("\n" + "="*80)
    print("EXAMPLE 6: Comprehensive Drug Target Report")
    print("="*80)

    # Target: EGFR (Epidermal Growth Factor Receptor)
    uniprot_id = "P00533"
    search_query = "EGFR kinase domain"

    print(f"\nTarget: EGFR (Epidermal Growth Factor Receptor)")
    print(f"UniProt ID: {uniprot_id}")

    # Initialize adapters
    alphafold = AlphaFoldAdapter(config={"download_pdb": False})
    rcsb = RCSBPDBAdapter()
    swissmodel = SwissModelAdapter(config={"download_pdb": False})

    # Gather data
    print("\nGathering structure data from all sources...")

    af_result = await alphafold.execute(uniprot_id)
    search_result = await rcsb.search_structures(search_query, max_results=5)
    sm_result = await swissmodel.execute(uniprot_id)

    # Generate report
    print("\n" + "="*80)
    print("DRUG TARGET REPORT: EGFR")
    print("="*80)

    # AlphaFold section
    if af_result.success:
        print("\n1. AI-PREDICTED STRUCTURE (AlphaFold)")
        print("-"*80)
        data = af_result.data
        print(f"  Model Version: {data['model_version']}")
        print(f"  Mean pLDDT: {data['plddt_scores']['mean']:.1f}")
        print(f"  High confidence regions: {data['plddt_scores']['high_confidence_pct']:.1f}%")
        print(f"  Medium confidence regions: {data['plddt_scores']['medium_confidence_pct']:.1f}%")
        print(f"  Low confidence regions: {data['plddt_scores']['low_confidence_pct']:.1f}%")

        if data['plddt_scores']['mean'] > 80:
            print("\n  ✓ ASSESSMENT: High-quality prediction suitable for structure-based design")
        elif data['plddt_scores']['mean'] > 60:
            print("\n  ⚠ ASSESSMENT: Medium-quality prediction, validate with experimental data")
        else:
            print("\n  ✗ ASSESSMENT: Low-quality prediction, use experimental structures")

    # RCSB PDB section
    if search_result.success:
        print("\n2. EXPERIMENTAL STRUCTURES (RCSB PDB)")
        print("-"*80)
        pdb_ids = search_result.data['pdb_ids']
        print(f"  Structures Found: {len(pdb_ids)}")
        print(f"  Top PDB IDs: {', '.join(pdb_ids)}")
        print("\n  ✓ ASSESSMENT: Experimental structures available for validation")

    # SWISS-MODEL section
    if sm_result.success:
        print("\n3. HOMOLOGY MODELS (SWISS-MODEL)")
        print("-"*80)
        data = sm_result.data
        print(f"  Models Available: {data['model_count']}")
        print(f"  Best QMEAN: {data['quality_metrics']['qmean']}")
        print(f"  Best GMQE: {data['quality_metrics']['gmqe']}")
        print(f"  Sequence Identity: {data['quality_metrics'].get('sequence_identity', 'N/A')}%")

    # Overall recommendation
    print("\n4. RECOMMENDATIONS")
    print("-"*80)
    print("\n  For structure-based drug design:")
    if search_result.success and len(search_result.data['pdb_ids']) > 0:
        print("  1. Use experimental structures from RCSB PDB (highest accuracy)")
        print("  2. Validate with AlphaFold prediction")
        print("  3. Consider multiple conformations from different PDB structures")
    elif af_result.success and af_result.data['plddt_scores']['mean'] > 70:
        print("  1. Use AlphaFold prediction (no experimental structures available)")
        print("  2. Focus on high-confidence regions")
        print("  3. Validate binding site predictions")
    else:
        print("  ⚠ Limited high-quality structural data available")
        print("  Consider homology modeling or experimental structure determination")

    print("\n" + "="*80)


# ============================================================================
# Main Menu
# ============================================================================

async def main():
    """
    Run integration examples.
    """
    print("\n" + "="*80)
    print("PROTEIN STRUCTURE ADAPTERS - INTEGRATION EXAMPLES")
    print("="*80)
    print("\nThese examples demonstrate how to use AlphaFold, RCSB PDB,")
    print("and SWISS-MODEL adapters for structure-based drug discovery.")

    examples = [
        ("Complete Structure Retrieval Workflow", example1_complete_workflow),
        ("Drug Target Analysis", example2_drug_target_analysis),
        ("Batch Protein Analysis", example3_batch_analysis),
        ("Structure Search and Discovery", example4_structure_search),
        ("Quality-Based Model Selection", example5_quality_selection),
        ("Comprehensive Drug Target Report", example6_drug_target_report)
    ]

    print("\n" + "-"*80)
    print("Available Examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")
    print("  0. Run all examples")
    print("-"*80)

    try:
        choice = input("\nSelect example (0-6): ").strip()

        if choice == "0":
            # Run all examples
            for name, func in examples:
                await func()
                print("\n" + "="*80)
                input("Press Enter to continue to next example...")
        elif choice.isdigit() and 1 <= int(choice) <= len(examples):
            # Run selected example
            _, func = examples[int(choice) - 1]
            await func()
        else:
            print("Invalid choice. Running first example...")
            await examples[0][1]()

        print("\n" + "="*80)
        print("EXAMPLES COMPLETED")
        print("="*80)

    except KeyboardInterrupt:
        print("\n\nExamples interrupted by user.")
    except Exception as e:
        logger.error(f"Example failed: {e}", exc_info=True)
        print(f"\n✗ EXAMPLE FAILED: {e}")


if __name__ == "__main__":
    asyncio.run(main())
