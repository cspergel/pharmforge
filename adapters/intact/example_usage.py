"""
IntAct Adapter - Example Usage Scripts

Demonstrates various ways to query the IntAct molecular interaction database
for protein-protein interaction data with experimental evidence.
"""
import asyncio
import json
from adapters.intact import IntActAdapter


async def example_1_basic_query():
    """Example 1: Basic query for single protein interactions"""
    print("\n" + "="*80)
    print("Example 1: Basic Query - TP53 Interactions")
    print("="*80)

    adapter = IntActAdapter()

    # Query TP53 (UniProt ID: P04637)
    result = await adapter.execute(
        "P04637",
        organism="human",
        min_mi_score=0.4
    )

    if result.success:
        data = result.data
        print(f"\nFound {data['num_interactions']} interactions for TP53")
        print(f"Involving {data['num_proteins']} unique proteins")

        print("\nTop 10 interactions by MI score:")
        for i, interaction in enumerate(data['interactions'][:10], 1):
            print(f"{i}. {interaction['protein_a_name']} <-> {interaction['protein_b_name']}")
            print(f"   MI Score: {interaction['mi_score']:.3f}")
            print(f"   Method: {interaction['detection_method']}")
            print(f"   Type: {interaction['interaction_type']}")
            print(f"   Publication: {interaction['publication']}")
    else:
        print(f"Error: {result.error}")


async def example_2_multiple_proteins():
    """Example 2: Query interactions for multiple proteins"""
    print("\n" + "="*80)
    print("Example 2: Multiple Proteins - TP53, MDM2, BRCA1 Network")
    print("="*80)

    adapter = IntActAdapter()

    # Query TP53, MDM2, and BRCA1
    proteins = ["P04637", "Q00987", "P38398"]  # TP53, MDM2, BRCA1
    gene_names = ["TP53", "MDM2", "BRCA1"]

    result = await adapter.execute(
        proteins,
        organism="human",
        min_mi_score=0.5
    )

    if result.success:
        data = result.data
        print(f"\nAnalyzing network for: {', '.join(gene_names)}")
        print(f"Total interactions found: {data['num_interactions']}")
        print(f"Total proteins in network: {data['num_proteins']}")

        # Find interactions between query proteins
        print("\nDirect interactions between query proteins:")
        for interaction in data['interactions']:
            prot_a = interaction['protein_a']
            prot_b = interaction['protein_b']
            if prot_a in proteins and prot_b in proteins:
                print(f"  {interaction['protein_a_name']} <-> {interaction['protein_b_name']}")
                print(f"    MI Score: {interaction['mi_score']:.3f}")
                print(f"    Method: {interaction['detection_method']}")

        # Detection methods used
        print(f"\nDetection methods used ({len(data['detection_methods'])}):")
        for method in data['detection_methods'][:10]:
            print(f"  - {method}")
    else:
        print(f"Error: {result.error}")


async def example_3_high_confidence_only():
    """Example 3: Filter for high-confidence interactions only"""
    print("\n" + "="*80)
    print("Example 3: High-Confidence Interactions (MI score >= 0.7)")
    print("="*80)

    adapter = IntActAdapter()

    result = await adapter.execute(
        "P04637",  # TP53
        organism="human",
        min_mi_score=0.7,  # High confidence only
        max_results=50
    )

    if result.success:
        data = result.data
        print(f"\nFound {data['num_interactions']} high-confidence interactions")

        # Group by interaction type
        by_type = {}
        for interaction in data['interactions']:
            itype = interaction['interaction_type']
            if itype not in by_type:
                by_type[itype] = []
            by_type[itype].append(interaction)

        print("\nInteractions grouped by type:")
        for itype, interactions in sorted(by_type.items()):
            print(f"\n{itype} ({len(interactions)} interactions):")
            for interaction in interactions[:3]:
                print(f"  - {interaction['protein_a_name']} <-> {interaction['protein_b_name']}")
                print(f"    MI Score: {interaction['mi_score']:.3f}")
    else:
        print(f"Error: {result.error}")


async def example_4_phosphorylation_events():
    """Example 4: Find phosphorylation interactions"""
    print("\n" + "="*80)
    print("Example 4: Phosphorylation Events")
    print("="*80)

    adapter = IntActAdapter()

    result = await adapter.execute(
        ["P04637", "Q00987"],  # TP53, MDM2
        organism="human",
        interaction_type="phosphorylation",
        min_mi_score=0.5
    )

    if result.success:
        data = result.data
        print(f"\nFound {data['num_interactions']} phosphorylation events")

        for interaction in data['interactions'][:10]:
            print(f"\n{interaction['protein_a_name']} -> {interaction['protein_b_name']}")
            print(f"  MI Score: {interaction['mi_score']:.3f}")
            print(f"  Evidence: {interaction['detection_method']}")
            print(f"  Source: {interaction['publication']}")
    else:
        print(f"Error: {result.error}")


async def example_5_custom_psicquic_query():
    """Example 5: Advanced PSICQUIC query"""
    print("\n" + "="*80)
    print("Example 5: Custom PSICQUIC Query")
    print("="*80)

    adapter = IntActAdapter()

    # Complex query: TP53 interactions in humans with experimental evidence
    query = "id:P04637 AND taxid:9606 AND detmethod:\"anti tag coimmunoprecipitation\""

    result = await adapter.execute(
        {"query": query},
        min_mi_score=0.5
    )

    if result.success:
        data = result.data
        print(f"\nQuery: {query}")
        print(f"Results: {data['num_interactions']} interactions")

        print("\nTop 5 results:")
        for i, interaction in enumerate(data['interactions'][:5], 1):
            print(f"\n{i}. {interaction['protein_a_name']} <-> {interaction['protein_b_name']}")
            print(f"   MI Score: {interaction['mi_score']:.3f}")
            print(f"   Method: {interaction['detection_method']}")
    else:
        print(f"Error: {result.error}")


async def example_6_gene_name_search():
    """Example 6: Search by gene names instead of UniProt IDs"""
    print("\n" + "="*80)
    print("Example 6: Search by Gene Names")
    print("="*80)

    adapter = IntActAdapter()

    result = await adapter.execute(
        {"identifiers": ["TP53", "MDM2"], "organism": "human"},
        min_mi_score=0.6
    )

    if result.success:
        data = result.data
        print(f"\nSearching for: TP53, MDM2")
        print(f"Found {data['num_interactions']} interactions")
        print(f"Unique proteins: {data['num_proteins']}")

        print("\nAll proteins in network:")
        for protein in sorted(data['proteins'])[:20]:
            print(f"  - {protein}")
    else:
        print(f"Error: {result.error}")


async def example_7_compare_organisms():
    """Example 7: Compare same protein across organisms"""
    print("\n" + "="*80)
    print("Example 7: Cross-Organism Comparison - TP53 homologs")
    print("="*80)

    adapter = IntActAdapter()

    organisms = {
        "Human": ("P04637", "9606"),
        "Mouse": ("P02340", "10090"),
    }

    for org_name, (uniprot_id, taxid) in organisms.items():
        result = await adapter.execute(
            uniprot_id,
            organism=taxid,
            min_mi_score=0.5,
            max_results=50
        )

        if result.success:
            data = result.data
            print(f"\n{org_name} (Tax ID: {taxid}):")
            print(f"  UniProt ID: {uniprot_id}")
            print(f"  Interactions: {data['num_interactions']}")
            print(f"  Unique partners: {data['num_proteins']}")
        else:
            print(f"\n{org_name}: Error - {result.error}")


async def example_8_export_network():
    """Example 8: Export network data for visualization"""
    print("\n" + "="*80)
    print("Example 8: Export Network for Visualization")
    print("="*80)

    adapter = IntActAdapter()

    result = await adapter.execute(
        ["P04637", "Q00987", "P38398"],  # TP53, MDM2, BRCA1
        organism="human",
        min_mi_score=0.6
    )

    if result.success:
        data = result.data

        # Create node list
        nodes = []
        for protein in data['proteins']:
            nodes.append({
                "id": protein,
                "label": protein
            })

        # Create edge list
        edges = []
        for i, interaction in enumerate(data['interactions']):
            edges.append({
                "id": i,
                "source": interaction['protein_a'],
                "target": interaction['protein_b'],
                "weight": interaction['mi_score'],
                "type": interaction['interaction_type'],
                "method": interaction['detection_method']
            })

        network = {
            "nodes": nodes,
            "edges": edges,
            "metadata": {
                "query_proteins": ["TP53", "MDM2", "BRCA1"],
                "organism": "human",
                "min_mi_score": 0.6
            }
        }

        # Save to JSON
        with open("intact_network.json", "w") as f:
            json.dump(network, f, indent=2)

        print(f"\nNetwork exported to intact_network.json")
        print(f"Nodes: {len(nodes)}")
        print(f"Edges: {len(edges)}")
        print("\nThis file can be imported into Cytoscape, Gephi, or other network visualization tools")
    else:
        print(f"Error: {result.error}")


async def example_9_batch_analysis():
    """Example 9: Batch analysis of cancer-related proteins"""
    print("\n" + "="*80)
    print("Example 9: Batch Analysis - Cancer Pathway Proteins")
    print("="*80)

    adapter = IntActAdapter()

    cancer_proteins = {
        "TP53": "P04637",
        "BRCA1": "P38398",
        "PTEN": "P60484",
        "RB1": "P06400",
        "APC": "P25054"
    }

    summary = []

    for gene, uniprot in cancer_proteins.items():
        result = await adapter.execute(
            uniprot,
            organism="human",
            min_mi_score=0.5,
            max_results=100
        )

        if result.success:
            data = result.data
            summary.append({
                "gene": gene,
                "uniprot": uniprot,
                "interactions": data['num_interactions'],
                "partners": data['num_proteins'],
                "top_partner": data['interactions'][0]['protein_b_name'] if data['interactions'] else "None",
                "top_score": data['interactions'][0]['mi_score'] if data['interactions'] else 0.0
            })

    print("\nCancer Protein Interaction Summary:")
    print(f"{'Gene':<10} {'UniProt':<10} {'Interactions':<15} {'Partners':<10} {'Top Partner':<20} {'Score':<10}")
    print("-" * 85)
    for item in summary:
        print(f"{item['gene']:<10} {item['uniprot']:<10} {item['interactions']:<15} {item['partners']:<10} {item['top_partner']:<20} {item['top_score']:<10.3f}")


async def example_10_caching_demo():
    """Example 10: Demonstrate caching functionality"""
    print("\n" + "="*80)
    print("Example 10: Caching Performance")
    print("="*80)

    adapter = IntActAdapter()
    protein_id = "P04637"

    import time

    # First query (no cache)
    start = time.time()
    result1 = await adapter.execute(protein_id, organism="human", use_cache=True)
    time1 = time.time() - start

    # Second query (from cache)
    start = time.time()
    result2 = await adapter.execute(protein_id, organism="human", use_cache=True)
    time2 = time.time() - start

    print(f"\nFirst query (no cache): {time1:.3f} seconds")
    print(f"Second query (cached): {time2:.3f} seconds")
    print(f"Speedup: {time1/time2:.1f}x faster")

    print(f"\nCache hit on second query: {result2.cache_hit}")
    print(f"Results identical: {result1.data['num_interactions'] == result2.data['num_interactions']}")


async def main():
    """Run all examples"""
    examples = [
        ("Basic Query", example_1_basic_query),
        ("Multiple Proteins", example_2_multiple_proteins),
        ("High Confidence Only", example_3_high_confidence_only),
        ("Phosphorylation Events", example_4_phosphorylation_events),
        ("Custom PSICQUIC Query", example_5_custom_psicquic_query),
        ("Gene Name Search", example_6_gene_name_search),
        ("Cross-Organism Comparison", example_7_compare_organisms),
        ("Export Network", example_8_export_network),
        ("Batch Analysis", example_9_batch_analysis),
        ("Caching Demo", example_10_caching_demo)
    ]

    print("\n" + "#"*80)
    print("#" + " "*78 + "#")
    print("#" + "  IntAct Adapter - Example Usage Demonstrations".center(78) + "#")
    print("#" + " "*78 + "#")
    print("#"*80)

    for i, (name, func) in enumerate(examples, 1):
        print(f"\n[{i}/{len(examples)}] Running: {name}")
        try:
            await func()
        except Exception as e:
            print(f"Error in {name}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "#"*80)
    print("#" + "  All examples completed!".center(78) + "#")
    print("#"*80 + "\n")


if __name__ == "__main__":
    # Run all examples
    asyncio.run(main())

    # Or run individual examples:
    # asyncio.run(example_1_basic_query())
    # asyncio.run(example_2_multiple_proteins())
    # asyncio.run(example_8_export_network())
