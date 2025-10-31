"""
BioGRID Adapter - Simple Example Usage
Quick start guide for using the BioGRID adapter
"""
import asyncio
import os
from adapters.biogrid import BioGRIDAdapter


async def simple_example():
    """
    Simple example: Query TP53 interactions
    """
    # Get your access key from: https://webservice.thebiogrid.org/
    access_key = os.environ.get("BIOGRID_ACCESS_KEY", "YOUR_ACCESS_KEY_HERE")

    # Initialize adapter
    adapter = BioGRIDAdapter(access_key=access_key)

    # Query TP53 protein interactions in humans
    result = await adapter.execute(
        "TP53",
        organism="human",
        max_results=50
    )

    # Process results
    if result.success:
        data = result.data
        print(f"\nTP53 Interactions Summary:")
        print(f"  Total interactions: {data['num_interactions']}")
        print(f"  Unique interacting genes: {data['summary']['num_unique_genes']}")
        print(f"  Publications: {data['summary']['num_publications']}")

        print(f"\nTop 10 interacting partners:")
        for i, gene in enumerate(data['summary']['unique_genes'][:10], 1):
            print(f"  {i}. {gene}")

        print(f"\nFirst 3 interactions:")
        for interaction in data['interactions'][:3]:
            print(f"  - {interaction['gene_a']} <-> {interaction['gene_b']}")
            print(f"    Method: {interaction['experimental_system']}")
            print(f"    PubMed: {interaction['pubmed_id']}")
            print()
    else:
        print(f"Error: {result.error}")


async def multi_gene_example():
    """
    Query multiple genes at once
    """
    access_key = os.environ.get("BIOGRID_ACCESS_KEY", "YOUR_ACCESS_KEY_HERE")
    adapter = BioGRIDAdapter(access_key=access_key)

    # Query multiple cancer-related genes
    genes = ["BRCA1", "BRCA2", "TP53"]

    result = await adapter.execute(
        genes,
        organism="human",
        max_results=100
    )

    if result.success:
        data = result.data
        print(f"\nQueried genes: {', '.join(genes)}")
        print(f"Total interactions found: {data['num_interactions']}")
        print(f"Network size: {data['summary']['num_unique_genes']} genes")


async def physical_interactions_example():
    """
    Query only physical (direct) interactions
    """
    access_key = os.environ.get("BIOGRID_ACCESS_KEY", "YOUR_ACCESS_KEY_HERE")
    adapter = BioGRIDAdapter(access_key=access_key)

    # Get only high-confidence physical interactions
    result = await adapter.execute(
        "EGFR",
        organism="human",
        evidence_types=[
            "Affinity Capture-MS",
            "Affinity Capture-Western",
            "Co-crystal Structure"
        ],
        max_results=50
    )

    if result.success:
        data = result.data
        print(f"\nEGFR Physical Interactions:")
        print(f"  Total: {data['num_interactions']}")
        print(f"  Methods used: {', '.join(data['summary']['experimental_systems'])}")


async def get_available_organisms():
    """
    List all available organisms in BioGRID
    """
    access_key = os.environ.get("BIOGRID_ACCESS_KEY", "YOUR_ACCESS_KEY_HERE")
    adapter = BioGRIDAdapter(access_key=access_key)

    result = await adapter.execute("", query_type="organisms")

    if result.success:
        organisms = result.data["organisms"]
        print(f"\nBioGRID contains data for {len(organisms)} organisms")
        print("\nSample organisms:")
        for org in organisms[:10]:
            if isinstance(org, dict):
                name = org.get("ORGANISM_OFFICIAL_NAME", "Unknown")
                tax_id = org.get("ORGANISM_ID", "N/A")
                print(f"  - {name} (Tax ID: {tax_id})")


if __name__ == "__main__":
    print("BioGRID Adapter - Example Usage")
    print("="*50)
    print("\nGet your free access key at:")
    print("https://webservice.thebiogrid.org/")
    print("="*50)

    # Run examples
    asyncio.run(simple_example())
    # asyncio.run(multi_gene_example())
    # asyncio.run(physical_interactions_example())
    # asyncio.run(get_available_organisms())
