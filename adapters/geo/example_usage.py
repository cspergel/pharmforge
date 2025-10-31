"""
Simple examples demonstrating GEO Adapter usage
Run these examples to see the adapter in action
"""
import asyncio
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from adapters.geo.adapter import GEOAdapter


async def example_1_gene_search():
    """Example 1: Search for datasets related to a gene"""
    print("\n" + "="*60)
    print("EXAMPLE 1: Search for TP53 gene expression datasets")
    print("="*60 + "\n")

    adapter = GEOAdapter(email="test@example.com")
    result = await adapter.execute({
        "query": "TP53",
        "max_results": 5
    })

    if result.success:
        data = result.data
        print(f"Found {data['total_count']} total datasets")
        print(f"Showing first {data['returned_count']}:\n")

        for i, ds in enumerate(data['datasets'], 1):
            print(f"{i}. {ds['accession']}")
            print(f"   Title: {ds['title'][:70]}...")
            print(f"   Organism: {ds['taxon']}")
            print(f"   Samples: {ds['n_samples']}\n")
    else:
        print(f"Error: {result.error}")


async def example_2_disease_search():
    """Example 2: Search for disease-related datasets"""
    print("\n" + "="*60)
    print("EXAMPLE 2: Search for Alzheimer's disease datasets")
    print("="*60 + "\n")

    adapter = GEOAdapter(email="test@example.com")
    result = await adapter.execute({
        "query": "Alzheimer disease",
        "max_results": 5,
        "filters": {
            "organism": "Homo sapiens"
        }
    })

    if result.success:
        data = result.data
        print(f"Found {data['total_count']} total datasets in humans")
        print(f"\nShowing first {data['returned_count']}:\n")

        for i, ds in enumerate(data['datasets'], 1):
            print(f"{i}. {ds['accession']} - {ds['gpl']}")
            print(f"   {ds['title'][:70]}...")
            print(f"   Samples: {ds['n_samples']}")

            # Check if there are PubMed links
            if ds.get('pubmedids'):
                pmids = ds['pubmedids'][:2]
                print(f"   PubMed: {', '.join(map(str, pmids))}")
            print()
    else:
        print(f"Error: {result.error}")


async def example_3_accession_lookup():
    """Example 3: Get specific dataset by accession"""
    print("\n" + "="*60)
    print("EXAMPLE 3: Lookup specific dataset GSE1234")
    print("="*60 + "\n")

    adapter = GEOAdapter(email="test@example.com")
    result = await adapter.execute("GSE1234")

    if result.success:
        dataset = result.data['dataset']
        print(f"Accession: {dataset['accession']}")
        print(f"Title: {dataset['title']}")
        print(f"Organism: {dataset['taxon']}")
        print(f"Platform: {dataset['gpl']}")
        print(f"Samples: {dataset['n_samples']}")
        print(f"\nDataset Type: {dataset['gdstype']}")
        print(f"Platform Tech: {dataset['ptechtype']}")

        print(f"\nDownload URLs:")
        for file_type, url in dataset['download_urls'].items():
            if url:
                print(f"  {file_type}: {url}")
    else:
        print(f"Error: {result.error}")


async def example_4_drug_response():
    """Example 4: Search for drug response studies"""
    print("\n" + "="*60)
    print("EXAMPLE 4: Search for drug response studies")
    print("="*60 + "\n")

    adapter = GEOAdapter(email="test@example.com")
    result = await adapter.execute({
        "query": "doxorubicin AND response",
        "max_results": 5,
        "filters": {
            "organism": "Homo sapiens"
        }
    })

    if result.success:
        data = result.data
        print(f"Found {data['total_count']} total datasets")

        summary = data['summary']
        print(f"\nPlatform types used:")
        for ptype in summary['platform_types'][:3]:
            print(f"  - {ptype['type']}: {ptype['count']} datasets")

        print(f"\nDatasets:\n")
        for i, ds in enumerate(data['datasets'], 1):
            print(f"{i}. {ds['accession']}")
            print(f"   {ds['title'][:70]}...")
            print()
    else:
        print(f"Error: {result.error}")


async def example_5_complex_query():
    """Example 5: Complex boolean query with filters"""
    print("\n" + "="*60)
    print("EXAMPLE 5: Complex query - Cancer AND Gene")
    print("="*60 + "\n")

    adapter = GEOAdapter(email="test@example.com")
    result = await adapter.execute({
        "query": "breast cancer AND BRCA1",
        "max_results": 5,
        "filters": {
            "organism": "Homo sapiens",
            "entry_type": "GSE"
        }
    })

    if result.success:
        data = result.data
        print(f"Query: {data['query']}")
        print(f"Found {data['total_count']} total datasets\n")

        for i, ds in enumerate(data['datasets'], 1):
            print(f"{i}. {ds['accession']}")
            print(f"   {ds['title'][:70]}...")
            print(f"   Type: {ds['gdstype']}")
            print(f"   Samples: {ds['n_samples']}")
            print()
    else:
        print(f"Error: {result.error}")


async def main():
    """Run all examples"""
    print("\n" + "="*60)
    print("GEO ADAPTER USAGE EXAMPLES")
    print("="*60)

    examples = [
        example_1_gene_search,
        example_2_disease_search,
        example_3_accession_lookup,
        example_4_drug_response,
        example_5_complex_query
    ]

    for example in examples:
        try:
            await example()
            await asyncio.sleep(0.5)  # Rate limiting
        except Exception as e:
            print(f"\nError in {example.__name__}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "="*60)
    print("Examples completed!")
    print("="*60 + "\n")


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\n\nExamples interrupted by user")
    except Exception as e:
        print(f"\n\nFatal error: {e}")
        import traceback
        traceback.print_exc()
