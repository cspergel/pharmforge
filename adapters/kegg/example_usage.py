"""
KEGG Adapter - Quick Usage Examples
Demonstrates common use cases for the KEGG adapter
"""
import asyncio
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.kegg.adapter import KEGGAdapter


async def example_1_simple_queries():
    """Example 1: Simple queries using auto-detection"""
    print("\n" + "="*80)
    print("Example 1: Simple Queries with Auto-Detection")
    print("="*80)

    adapter = KEGGAdapter()

    # Query a compound by ID
    print("\n1. Query Glucose (C00031):")
    result = await adapter.execute("C00031")
    if result.success:
        print(f"   Name: {result.data['data']['NAME']}")
        print(f"   Formula: {result.data['data']['FORMULA']}")

    # Query a gene by ID
    print("\n2. Query TP53 gene (hsa:7157):")
    result = await adapter.execute("hsa:7157")
    if result.success:
        print(f"   Name: {result.data['data']['NAME']}")
        print(f"   Associated pathways: {len(result.data['pathways'])}")

    # Query a disease by ID
    print("\n3. Query Alzheimer's disease (H00056):")
    result = await adapter.execute("H00056")
    if result.success:
        print(f"   Name: {result.data['data']['NAME']}")

    # Query a pathway by ID
    print("\n4. Query Glycolysis pathway (hsa00010):")
    result = await adapter.execute("hsa00010")
    if result.success:
        print(f"   Name: {result.data['data']['NAME']}")


async def example_2_search_operations():
    """Example 2: Searching databases"""
    print("\n" + "="*80)
    print("Example 2: Search Operations")
    print("="*80)

    adapter = KEGGAdapter()

    # Search for compounds
    print("\n1. Search compounds for 'caffeine':")
    result = await adapter.execute({
        "operation": "query_compound",
        "search_term": "caffeine"
    })
    if result.success:
        for compound in result.data['results'][:3]:
            print(f"   {compound['id']}: {compound['description'][:60]}...")

    await asyncio.sleep(0.5)

    # Search for pathways
    print("\n2. Search pathways for 'metabolism':")
    result = await adapter.execute({
        "operation": "query_pathway",
        "search_term": "metabolism"
    })
    if result.success:
        print(f"   Found {result.data['count']} pathways")
        for pathway in result.data['results'][:3]:
            print(f"   {pathway['id']}: {pathway['description'][:60]}...")

    await asyncio.sleep(0.5)

    # Search for diseases
    print("\n3. Search diseases for 'diabetes':")
    result = await adapter.execute({
        "operation": "query_disease",
        "search_term": "diabetes"
    })
    if result.success:
        for disease in result.data['results'][:3]:
            print(f"   {disease['id']}: {disease['description'][:60]}...")


async def example_3_pathway_analysis():
    """Example 3: Pathway analysis workflow"""
    print("\n" + "="*80)
    print("Example 3: Pathway Analysis Workflow")
    print("="*80)

    adapter = KEGGAdapter()

    # Step 1: Search for a compound
    print("\n1. Search for 'glucose':")
    result = await adapter.execute({
        "operation": "query_compound",
        "search_term": "glucose"
    })

    if result.success and result.data['results']:
        compound_id = result.data['results'][0]['id'].replace('cpd:', '')
        print(f"   Found: {compound_id}")

        await asyncio.sleep(0.5)

        # Step 2: Get compound details with pathways
        print(f"\n2. Get pathways for {compound_id}:")
        result = await adapter.execute({
            "operation": "query_compound",
            "compound_id": compound_id
        })

        if result.success:
            pathways = result.data.get('pathways', [])[:5]
            print(f"   Found in {len(result.data.get('pathways', []))} pathways")

            # Step 3: Get details for first pathway
            if pathways:
                pathway_id = pathways[0]['target'].replace('path:', '')
                print(f"\n3. Get details for pathway {pathway_id}:")

                await asyncio.sleep(0.5)

                pathway_result = await adapter.execute({
                    "operation": "query_pathway",
                    "pathway_id": pathway_id
                })

                if pathway_result.success:
                    print(f"   Name: {pathway_result.data['data'].get('NAME', 'N/A')}")
                    print(f"   Class: {pathway_result.data['data'].get('CLASS', 'N/A')}")


async def example_4_gene_to_pathway():
    """Example 4: Gene to pathway mapping"""
    print("\n" + "="*80)
    print("Example 4: Gene to Pathway Mapping")
    print("="*80)

    adapter = KEGGAdapter()

    gene_names = ["BRCA1", "TP53"]

    for gene_name in gene_names:
        print(f"\n1. Search for gene '{gene_name}':")
        result = await adapter.execute({
            "operation": "query_gene",
            "search_term": gene_name
        })

        if result.success and result.data['results']:
            gene_id = result.data['results'][0]['id']
            print(f"   Found: {gene_id}")

            await asyncio.sleep(0.5)

            # Get gene details
            print(f"\n2. Get pathways for {gene_id}:")
            gene_result = await adapter.execute({
                "operation": "query_gene",
                "gene_id": gene_id
            })

            if gene_result.success:
                pathways = gene_result.data.get('pathways', [])
                print(f"   Gene name: {gene_result.data['data'].get('NAME', 'N/A')}")
                print(f"   Found in {len(pathways)} pathways")

                # Show first 3 pathways
                for i, p in enumerate(pathways[:3], 1):
                    print(f"      {i}. {p['target']}")

        await asyncio.sleep(0.5)


async def example_5_disease_analysis():
    """Example 5: Disease to genes and pathways"""
    print("\n" + "="*80)
    print("Example 5: Disease Analysis - Genes and Pathways")
    print("="*80)

    adapter = KEGGAdapter()

    # Search for cancer-related disease
    print("\n1. Search for 'breast cancer':")
    result = await adapter.execute({
        "operation": "query_disease",
        "search_term": "breast cancer"
    })

    if result.success and result.data['results']:
        disease_id = result.data['results'][0]['id'].replace('ds:', '')
        print(f"   Found: {disease_id}")
        print(f"   Description: {result.data['results'][0]['description'][:70]}...")

        await asyncio.sleep(0.5)

        # Get disease details
        print(f"\n2. Get genes and pathways for {disease_id}:")
        disease_result = await adapter.execute({
            "operation": "query_disease",
            "disease_id": disease_id
        })

        if disease_result.success:
            genes = disease_result.data.get('genes', [])
            pathways = disease_result.data.get('pathways', [])

            print(f"   Disease: {disease_result.data['data'].get('NAME', 'N/A')}")
            print(f"   Associated genes: {len(genes)}")
            print(f"   Associated pathways: {len(pathways)}")

            if pathways:
                print("\n   Sample pathways:")
                for i, p in enumerate(pathways[:3], 1):
                    print(f"      {i}. {p['target']}")


async def example_6_chemical_formula_search():
    """Example 6: Search by chemical formula"""
    print("\n" + "="*80)
    print("Example 6: Chemical Formula Search")
    print("="*80)

    adapter = KEGGAdapter()

    # Search by formula
    formula = "C6H12O6"
    print(f"\n1. Search compounds with formula {formula}:")
    result = await adapter.execute({
        "operation": "query_compound",
        "formula": formula
    })

    if result.success:
        print(f"   Found {result.data['count']} compounds")
        for i, compound in enumerate(result.data['results'][:5], 1):
            print(f"   {i}. {compound['id']}: {compound['description'][:60]}...")


async def example_7_drug_target_analysis():
    """Example 7: Drug and target analysis"""
    print("\n" + "="*80)
    print("Example 7: Drug and Target Analysis")
    print("="*80)

    adapter = KEGGAdapter()

    # Search for a drug
    print("\n1. Search for 'imatinib' (cancer drug):")
    result = await adapter.execute({
        "operation": "query_drug",
        "search_term": "imatinib"
    })

    if result.success and result.data['results']:
        drug_id = result.data['results'][0]['id'].replace('dr:', '')
        print(f"   Found: {drug_id}")

        await asyncio.sleep(0.5)

        # Get drug details
        print(f"\n2. Get targets for {drug_id}:")
        drug_result = await adapter.execute({
            "operation": "query_drug",
            "drug_id": drug_id
        })

        if drug_result.success:
            print(f"   Drug: {drug_result.data['data'].get('NAME', 'N/A')}")
            targets = drug_result.data.get('targets', [])
            print(f"   Number of targets: {len(targets)}")

            if targets:
                print("\n   Sample targets:")
                for i, t in enumerate(targets[:5], 1):
                    print(f"      {i}. {t['target']}")


async def run_all_examples():
    """Run all example workflows"""
    print("\n" + "#"*80)
    print("# KEGG ADAPTER - USAGE EXAMPLES")
    print("#"*80)

    try:
        await example_1_simple_queries()
        await example_2_search_operations()
        await example_3_pathway_analysis()
        await example_4_gene_to_pathway()
        await example_5_disease_analysis()
        await example_6_chemical_formula_search()
        await example_7_drug_target_analysis()

        print("\n" + "="*80)
        print("All examples completed!")
        print("="*80)

    except Exception as e:
        print(f"\n!!! ERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(run_all_examples())
