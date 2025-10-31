"""
HMDB Adapter - Example Usage Scripts

Demonstrates various use cases for the HMDB adapter in drug discovery workflows.
"""

import asyncio
import logging
from typing import List, Dict, Any

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import adapter
try:
    from adapters.hmdb import HMDBAdapter
except ImportError:
    # Fallback for running from adapters/hmdb directory
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent.parent))
    from adapters.hmdb import HMDBAdapter


# ============================================================================
# Example 1: Basic Metabolite Lookup
# ============================================================================

async def example_basic_lookup():
    """
    Look up a single metabolite by HMDB ID
    """
    print("\n" + "="*70)
    print("EXAMPLE 1: Basic Metabolite Lookup")
    print("="*70)

    adapter = HMDBAdapter()

    # Query 1-Methylhistidine
    result = await adapter.execute("HMDB0000001")

    if result.success:
        metabolite = result.data["metabolite"]
        summary = result.data["summary"]

        print(f"\nMetabolite: {metabolite['name']}")
        print(f"HMDB ID: {metabolite['hmdb_id']}")
        print(f"Formula: {metabolite['chemical_formula']}")
        print(f"Molecular Weight: {metabolite['average_molecular_weight']}")
        print(f"SMILES: {metabolite['smiles']}")
        print(f"\nSummary Statistics:")
        print(f"  - Found in {summary['num_biofluids']} biofluids")
        print(f"  - Found in {summary['num_tissues']} tissues")
        print(f"  - {summary['num_diseases']} disease associations")
        print(f"  - {summary['num_pathways']} metabolic pathways")
        print(f"  - {summary['num_proteins']} protein associations")

        # Show some synonyms
        if metabolite.get("synonyms"):
            print(f"\nAlternative Names:")
            for syn in metabolite["synonyms"][:5]:
                print(f"  - {syn}")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 2: Biofluid Concentration Analysis
# ============================================================================

async def example_biofluid_concentrations():
    """
    Analyze metabolite concentrations in different biofluids
    """
    print("\n" + "="*70)
    print("EXAMPLE 2: Biofluid Concentration Analysis")
    print("="*70)

    adapter = HMDBAdapter()

    # Query with specific biofluid filtering
    result = await adapter.execute(
        "HMDB0000001",
        biofluids=["blood", "urine", "cerebrospinal_fluid"],
        include_concentrations=True
    )

    if result.success:
        metabolite = result.data["metabolite"]
        print(f"\nMetabolite: {metabolite['name']}")
        print(f"\nNormal Concentration Ranges:")

        concentrations = metabolite.get("concentrations", {})
        if concentrations:
            for biofluid, data in concentrations.items():
                biofluid_name = biofluid.replace("_", " ").title()
                print(f"\n  {biofluid_name}:")
                print(f"    Range: {data.get('value', 'N/A')} {data.get('unit', '')}")

                if data.get('subject_age'):
                    print(f"    Age: {data['subject_age']}")
                if data.get('subject_condition'):
                    print(f"    Condition: {data['subject_condition']}")
                if data.get('references'):
                    refs = data['references'][:3]
                    print(f"    References: {', '.join(refs)}")
        else:
            print("  No concentration data available")

        # Show where metabolite is found
        print(f"\nBiofluid Locations:")
        for location in metabolite.get("biofluid_locations", []):
            print(f"  - {location}")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 3: Disease Association Discovery
# ============================================================================

async def example_disease_associations():
    """
    Find disease associations for a metabolite
    """
    print("\n" + "="*70)
    print("EXAMPLE 3: Disease Association Discovery")
    print("="*70)

    adapter = HMDBAdapter()

    result = await adapter.execute(
        "HMDB0000001",
        include_diseases=True
    )

    if result.success:
        metabolite = result.data["metabolite"]
        print(f"\nMetabolite: {metabolite['name']}")
        print(f"\nAssociated Diseases:")

        diseases = metabolite.get("diseases", [])
        if diseases:
            for i, disease in enumerate(diseases[:10], 1):
                print(f"\n{i}. {disease['name']}")
                if disease.get('omim_id'):
                    print(f"   OMIM ID: {disease['omim_id']}")
                if disease.get('references'):
                    print(f"   References: {len(disease['references'])} PubMed articles")
        else:
            print("  No disease associations found")

        # Biomarker potential
        ontology = metabolite.get("ontology", {})
        if ontology.get("applications"):
            print(f"\nApplications:")
            for app in ontology["applications"]:
                print(f"  - {app}")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 4: Metabolic Pathway Analysis
# ============================================================================

async def example_pathway_analysis():
    """
    Analyze metabolic pathways and protein interactions
    """
    print("\n" + "="*70)
    print("EXAMPLE 4: Metabolic Pathway Analysis")
    print("="*70)

    adapter = HMDBAdapter()

    result = await adapter.execute(
        "HMDB0000001",
        include_pathways=True,
        include_proteins=True
    )

    if result.success:
        metabolite = result.data["metabolite"]
        print(f"\nMetabolite: {metabolite['name']}")

        # Pathways
        print(f"\nMetabolic Pathways:")
        pathways = metabolite.get("pathways", [])
        if pathways:
            for pathway in pathways:
                print(f"\n  - {pathway['name']}")
                if pathway.get('smpdb_id'):
                    print(f"    SMPDB ID: {pathway['smpdb_id']}")
                if pathway.get('kegg_map_id'):
                    print(f"    KEGG Map: {pathway['kegg_map_id']}")
        else:
            print("  No pathway data available")

        # Proteins
        print(f"\nAssociated Proteins/Enzymes:")
        proteins = metabolite.get("proteins", [])
        if proteins:
            for i, protein in enumerate(proteins[:10], 1):
                gene_name = protein.get('gene_name', 'Unknown')
                uniprot_id = protein.get('uniprot_id', 'N/A')
                protein_type = protein.get('protein_type', 'N/A')
                name = protein.get('name', '')

                print(f"\n  {i}. {gene_name} ({uniprot_id})")
                print(f"     Type: {protein_type}")
                if name:
                    print(f"     Name: {name}")
        else:
            print("  No protein associations found")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 5: Cross-Reference Integration
# ============================================================================

async def example_cross_references():
    """
    Get cross-references to other databases
    """
    print("\n" + "="*70)
    print("EXAMPLE 5: Cross-Reference Integration")
    print("="*70)

    adapter = HMDBAdapter()

    result = await adapter.execute("HMDB0000001")

    if result.success:
        metabolite = result.data["metabolite"]
        ext_ids = metabolite.get("external_ids", {})

        print(f"\nMetabolite: {metabolite['name']}")
        print(f"\nExternal Database Cross-References:")
        print(f"  PubChem CID: {ext_ids.get('pubchem_compound_id', 'N/A')}")
        print(f"  KEGG ID: {ext_ids.get('kegg_id', 'N/A')}")
        print(f"  ChEBI ID: {ext_ids.get('chebi_id', 'N/A')}")
        print(f"  ChemSpider ID: {ext_ids.get('chemspider_id', 'N/A')}")
        print(f"  DrugBank ID: {ext_ids.get('drugbank_id', 'N/A')}")
        print(f"  FooDB ID: {ext_ids.get('foodb_id', 'N/A')}")

        print(f"\nStructure Identifiers:")
        print(f"  SMILES: {metabolite.get('smiles', 'N/A')}")
        print(f"  InChI: {metabolite.get('inchi', 'N/A')[:80]}...")
        print(f"  InChIKey: {metabolite.get('inchikey', 'N/A')}")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 6: Metabolite Search
# ============================================================================

async def example_search():
    """
    Search for metabolites by name
    """
    print("\n" + "="*70)
    print("EXAMPLE 6: Metabolite Search")
    print("="*70)

    adapter = HMDBAdapter()

    # Search for glucose-related metabolites
    result = await adapter.execute(
        {"query": "glucose", "mode": "search"},
        include_concentrations=False,
        include_diseases=False
    )

    if result.success:
        print(f"\nSearch query: '{result.data['query']}'")
        print(f"Total results found: {result.data['num_results']}")
        print(f"Details fetched for: {result.data['num_fetched']}")

        if result.data.get("warning"):
            print(f"\nNote: {result.data['warning']}")

        print(f"\nAll HMDB IDs found:")
        for hmdb_id in result.data["all_hmdb_ids"][:10]:
            print(f"  - {hmdb_id}")

        if len(result.data["all_hmdb_ids"]) > 10:
            print(f"  ... and {len(result.data['all_hmdb_ids']) - 10} more")

        print(f"\nDetailed Results:")
        for i, metabolite in enumerate(result.data["metabolites"], 1):
            print(f"\n{i}. {metabolite['name']} ({metabolite['hmdb_id']})")
            print(f"   Formula: {metabolite.get('chemical_formula', 'N/A')}")
            print(f"   MW: {metabolite.get('average_molecular_weight', 'N/A')}")

            locations = metabolite.get('biofluid_locations', [])
            if locations:
                print(f"   Found in: {', '.join(locations[:5])}")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 7: Drug Metabolism Analysis
# ============================================================================

async def example_drug_metabolism():
    """
    Analyze drug metabolism using metabolite data
    """
    print("\n" + "="*70)
    print("EXAMPLE 7: Drug Metabolism Analysis")
    print("="*70)

    adapter = HMDBAdapter()

    # Example: Analyze creatinine (common in kidney function tests)
    result = await adapter.execute(
        {"query": "creatinine", "mode": "search"},
        include_concentrations=True,
        include_diseases=True,
        include_proteins=True,
        biofluids=["blood", "urine"]
    )

    if result.success and result.data["metabolites"]:
        metabolite = result.data["metabolites"][0]

        print(f"\nMetabolite: {metabolite['name']}")
        print(f"Clinical Significance:")

        # Disease associations
        diseases = metabolite.get("diseases", [])
        if diseases:
            print(f"\n  Associated with {len(diseases)} disease conditions:")
            for disease in diseases[:5]:
                print(f"    - {disease['name']}")

        # Concentration data (important for clinical interpretation)
        concentrations = metabolite.get("concentrations", {})
        if concentrations:
            print(f"\n  Normal Ranges:")
            for biofluid, data in concentrations.items():
                print(f"    {biofluid.replace('_', ' ').title()}: {data.get('value', 'N/A')} {data.get('unit', '')}")

        # Protein interactions (enzymes involved in metabolism)
        proteins = metabolite.get("proteins", [])
        if proteins:
            print(f"\n  Metabolizing Enzymes:")
            for protein in proteins[:5]:
                print(f"    - {protein.get('gene_name', 'Unknown')} ({protein.get('protein_type', 'N/A')})")

    else:
        print(f"Error or no results: {result.error if not result.success else 'No metabolites found'}")


# ============================================================================
# Example 8: Batch Processing
# ============================================================================

async def example_batch_processing():
    """
    Process multiple metabolites efficiently
    """
    print("\n" + "="*70)
    print("EXAMPLE 8: Batch Processing Multiple Metabolites")
    print("="*70)

    adapter = HMDBAdapter()

    # List of HMDB IDs to process
    hmdb_ids = [
        "HMDB0000001",  # 1-Methylhistidine
        "HMDB0000122",  # D-Glucose
        "HMDB0000064",  # Creatinine
        "HMDB0000094",  # Citric acid
        "HMDB0000148",  # L-Glutamic acid
    ]

    print(f"\nProcessing {len(hmdb_ids)} metabolites...\n")

    results = []
    for hmdb_id in hmdb_ids:
        result = await adapter.execute(
            hmdb_id,
            include_concentrations=False,
            include_diseases=True,
            include_pathways=True
        )

        if result.success:
            metabolite = result.data["metabolite"]
            summary = result.data["summary"]

            results.append({
                "hmdb_id": metabolite["hmdb_id"],
                "name": metabolite["name"],
                "formula": metabolite["chemical_formula"],
                "mw": metabolite["average_molecular_weight"],
                "diseases": summary["num_diseases"],
                "pathways": summary["num_pathways"],
                "proteins": summary["num_proteins"]
            })

            print(f"✓ {metabolite['name']}")

    # Summary table
    print(f"\n{'Metabolite':<25} {'Formula':<15} {'MW':<10} {'Diseases':<10} {'Pathways':<10}")
    print("-" * 70)
    for r in results:
        print(f"{r['name']:<25} {r['formula']:<15} {r['mw']:<10} {r['diseases']:<10} {r['pathways']:<10}")


# ============================================================================
# Example 9: Clinical Biomarker Analysis
# ============================================================================

async def example_biomarker_analysis():
    """
    Analyze metabolites as clinical biomarkers
    """
    print("\n" + "="*70)
    print("EXAMPLE 9: Clinical Biomarker Analysis")
    print("="*70)

    adapter = HMDBAdapter()

    result = await adapter.execute(
        "HMDB0000001",
        include_concentrations=True,
        include_diseases=True,
        biofluids=["blood", "urine"]
    )

    if result.success:
        metabolite = result.data["metabolite"]

        print(f"\nBiomarker: {metabolite['name']}")
        print(f"HMDB ID: {metabolite['hmdb_id']}")

        # Biomarker potential
        ontology = metabolite.get("ontology", {})
        if ontology.get("applications"):
            print(f"\nApplications:")
            for app in ontology["applications"]:
                print(f"  - {app}")

        # Normal ranges
        print(f"\nNormal Reference Ranges:")
        concentrations = metabolite.get("concentrations", {})
        for biofluid, data in concentrations.items():
            print(f"  {biofluid.replace('_', ' ').title()}: {data.get('value', 'N/A')} {data.get('unit', '')}")

        # Disease associations
        print(f"\nClinical Significance:")
        diseases = metabolite.get("diseases", [])
        if diseases:
            for disease in diseases[:5]:
                print(f"  - {disease['name']}")
        else:
            print("  No specific disease associations")

        # Description
        if metabolite.get("description"):
            print(f"\nDescription:")
            desc = metabolite["description"][:300]
            print(f"  {desc}...")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Example 10: Integration with Other Adapters
# ============================================================================

async def example_multi_adapter_integration():
    """
    Demonstrate integration with other PharmForge adapters
    """
    print("\n" + "="*70)
    print("EXAMPLE 10: Multi-Adapter Integration")
    print("="*70)

    adapter = HMDBAdapter()

    # Step 1: Get metabolite from HMDB
    result = await adapter.execute("HMDB0000122")  # D-Glucose

    if result.success:
        metabolite = result.data["metabolite"]
        print(f"\nStep 1: Retrieved {metabolite['name']} from HMDB")
        print(f"  HMDB ID: {metabolite['hmdb_id']}")

        # Step 2: Get external IDs for cross-referencing
        ext_ids = metabolite.get("external_ids", {})

        print(f"\nStep 2: Found cross-references:")
        if ext_ids.get("pubchem_compound_id"):
            print(f"  ✓ PubChem CID: {ext_ids['pubchem_compound_id']}")
            print(f"    → Could query PubChem adapter for additional properties")

        if ext_ids.get("kegg_id"):
            print(f"  ✓ KEGG ID: {ext_ids['kegg_id']}")
            print(f"    → Could query KEGG adapter for pathway details")

        if ext_ids.get("chebi_id"):
            print(f"  ✓ ChEBI ID: {ext_ids['chebi_id']}")
            print(f"    → Could query ChEBI adapter for ontology")

        if ext_ids.get("drugbank_id"):
            print(f"  ✓ DrugBank ID: {ext_ids['drugbank_id']}")
            print(f"    → Could query DrugBank adapter for drug information")

        # Step 3: Show protein associations for UniProt integration
        proteins = metabolite.get("proteins", [])
        if proteins:
            print(f"\nStep 3: Found {len(proteins)} protein associations:")
            for protein in proteins[:3]:
                uniprot_id = protein.get("uniprot_id")
                if uniprot_id:
                    print(f"  ✓ UniProt ID: {uniprot_id} ({protein.get('gene_name', 'Unknown')})")
                    print(f"    → Could query UniProt adapter for protein details")

        # Step 4: Show pathway integration
        pathways = metabolite.get("pathways", [])
        if pathways:
            print(f"\nStep 4: Found {len(pathways)} metabolic pathways:")
            for pathway in pathways[:3]:
                if pathway.get("kegg_map_id"):
                    print(f"  ✓ KEGG Map: {pathway['kegg_map_id']} ({pathway['name']})")
                    print(f"    → Could query KEGG adapter for pathway visualization")

        print(f"\n{'─'*70}")
        print(f"Complete workflow demonstrates HMDB as a hub for multi-database integration!")

    else:
        print(f"Error: {result.error}")


# ============================================================================
# Main Execution
# ============================================================================

async def main():
    """
    Run all examples
    """
    print("\n")
    print("╔" + "═"*68 + "╗")
    print("║" + " "*15 + "HMDB Adapter - Example Usage" + " "*25 + "║")
    print("║" + " "*10 + "PharmForge Human Metabolome Database Integration" + " "*7 + "║")
    print("╚" + "═"*68 + "╝")

    examples = [
        ("Basic Metabolite Lookup", example_basic_lookup),
        ("Biofluid Concentration Analysis", example_biofluid_concentrations),
        ("Disease Association Discovery", example_disease_associations),
        ("Metabolic Pathway Analysis", example_pathway_analysis),
        ("Cross-Reference Integration", example_cross_references),
        ("Metabolite Search", example_search),
        ("Drug Metabolism Analysis", example_drug_metabolism),
        ("Batch Processing", example_batch_processing),
        ("Clinical Biomarker Analysis", example_biomarker_analysis),
        ("Multi-Adapter Integration", example_multi_adapter_integration),
    ]

    print("\nAvailable Examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")

    print("\n" + "="*70)
    print("Running all examples...")
    print("="*70)

    for name, example_func in examples:
        try:
            await example_func()
            await asyncio.sleep(1)  # Rate limiting between examples
        except Exception as e:
            logger.error(f"Error in {name}: {e}", exc_info=True)

    print("\n" + "="*70)
    print("All examples completed!")
    print("="*70 + "\n")


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
