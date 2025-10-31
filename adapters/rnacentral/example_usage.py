"""
Example usage of the RNAcentral Adapter
Demonstrates various ways to query RNA data for drug discovery applications
"""
import asyncio
import logging
from adapters.rnacentral.adapter import RNAcentralAdapter

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


async def example_1_rna_id_lookup():
    """
    Example 1: Look up a specific RNA by RNAcentral ID
    """
    print("\n" + "="*70)
    print("Example 1: RNA ID Lookup")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Look up a specific human microRNA (let-7a)
    # Note: Use a real RNAcentral ID in practice
    result = await adapter.execute("URS0000527F89_9606")  # let-7a in human

    if result.success:
        rna = result.data
        print(f"\nRNA ID: {rna['rnacentral_id']}")
        print(f"Type: {rna['rna_type']}")
        print(f"Description: {rna['description']}")
        print(f"Length: {rna['length']} nucleotides")
        print(f"Sequence: {rna['sequence'][:50]}..." if len(rna['sequence']) > 50 else rna['sequence'])
        print(f"URL: {rna['url']}")

        if rna['species']:
            print(f"\nSpecies:")
            for species in rna['species'][:3]:
                print(f"  - {species['name']} (taxid: {species['taxid']})")
    else:
        print(f"Error: {result.error}")


async def example_2_keyword_search():
    """
    Example 2: Search for RNAs by keyword
    """
    print("\n" + "="*70)
    print("Example 2: Keyword Search")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Search for microRNAs related to cancer
    result = await adapter.execute("cancer")

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} results for 'cancer'")

        print("\nRNA Types Distribution:")
        for rna_type in data['rna_types'][:5]:
            print(f"  {rna_type['type']}: {rna_type['count']} entries")

        print("\nTop 5 Results:")
        for i, entry in enumerate(data['entries'][:5], 1):
            print(f"\n{i}. {entry['rnacentral_id']}")
            print(f"   Type: {entry['rna_type']}, Length: {entry['length']} nt")
            print(f"   Description: {entry['description'][:80]}...")
            print(f"   URL: {entry['url']}")
    else:
        print(f"Error: {result.error}")


async def example_3_filter_by_rna_type():
    """
    Example 3: Search with RNA type filter
    """
    print("\n" + "="*70)
    print("Example 3: Filter by RNA Type")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Search specifically for long non-coding RNAs
    result = await adapter.execute(
        "gene regulation",
        rna_type="lncRNA"
    )

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} lncRNAs related to gene regulation")

        print("\nTop 5 lncRNAs:")
        for entry in data['entries'][:5]:
            print(f"\n- {entry['rnacentral_id']}")
            print(f"  Length: {entry['length']} nucleotides")
            print(f"  {entry['description'][:100]}...")
    else:
        print(f"Error: {result.error}")


async def example_4_organism_filter():
    """
    Example 4: Filter by organism
    """
    print("\n" + "="*70)
    print("Example 4: Filter by Organism")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Search for human microRNAs
    result = await adapter.execute(
        "regulatory",
        rna_type="miRNA",
        organism="Homo sapiens"
    )

    if result.success:
        data = result.data
        print(f"\nFound {data['num_results']} human microRNAs")

        print("\nSample Results:")
        for entry in data['entries'][:5]:
            print(f"- {entry['rnacentral_id']}: {entry['description'][:60]}...")
    else:
        print(f"Error: {result.error}")


async def example_5_get_cross_references():
    """
    Example 5: Get cross-references to other databases
    """
    print("\n" + "="*70)
    print("Example 5: Get Cross-References")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Get cross-references for a specific RNA
    result = await adapter.execute(
        "URS0000527F89_9606",  # let-7a
        get_xrefs=True
    )

    if result.success and "cross_references" in result.data:
        print(f"\nRNA: {result.data['rnacentral_id']}")
        print(f"Description: {result.data['description']}")

        print("\nCross-references to other databases:")
        for xref in result.data['cross_references'][:10]:
            db = xref.get('database', 'Unknown')
            acc = xref.get('accession', 'N/A')
            print(f"  {db}: {acc}")
    else:
        print(f"Error or no cross-references: {result.error}")


async def example_6_antisense_target_identification():
    """
    Example 6: Real-world use case - Identify antisense drug targets
    """
    print("\n" + "="*70)
    print("Example 6: Antisense Drug Target Identification")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Search for disease-associated microRNAs
    print("\nStep 1: Find disease-associated miRNAs...")
    result = await adapter.execute(
        "Alzheimer",
        rna_type="miRNA",
        organism="Homo sapiens"
    )

    if result.success:
        candidates = result.data['entries'][:5]
        print(f"Found {len(candidates)} candidate miRNAs")

        print("\nStep 2: Evaluate candidates for antisense targeting...")
        for candidate in candidates:
            rna_id = candidate['rnacentral_id']
            length = candidate['length']

            # Get detailed information
            detail_result = await adapter.execute(rna_id)

            if detail_result.success:
                sequence = detail_result.data.get('sequence', '')

                # Simple scoring (in practice, use more sophisticated analysis)
                score = 0
                if 18 <= length <= 25:  # Ideal miRNA length
                    score += 10
                if sequence:  # Has sequence data
                    score += 5

                print(f"\n{rna_id}")
                print(f"  Length: {length} nt (ideal: 18-25)")
                print(f"  Sequence available: {'Yes' if sequence else 'No'}")
                print(f"  Target Score: {score}/15")

                if score >= 10:
                    print(f"  ✓ Good candidate for ASO design")
                else:
                    print(f"  ✗ May need further evaluation")


async def example_7_rna_type_survey():
    """
    Example 7: Survey RNA types in a biological process
    """
    print("\n" + "="*70)
    print("Example 7: RNA Type Survey")
    print("="*70)

    adapter = RNAcentralAdapter()

    # Query different RNA types involved in immune response
    rna_types = ["miRNA", "lncRNA", "snoRNA", "piRNA"]

    print("\nRNA involvement in 'immune response':\n")

    results_summary = {}

    for rna_type in rna_types:
        result = await adapter.execute(
            "immune response",
            rna_type=rna_type,
            organism="Homo sapiens"
        )

        if result.success:
            count = result.data['num_results']
            results_summary[rna_type] = count
            print(f"{rna_type:10s}: {count:4d} entries")
        else:
            print(f"{rna_type:10s}: Error - {result.error}")

    # Show which RNA type is most prevalent
    if results_summary:
        most_common = max(results_summary.items(), key=lambda x: x[1])
        print(f"\nMost common: {most_common[0]} with {most_common[1]} entries")


async def example_8_batch_processing():
    """
    Example 8: Batch processing of multiple RNA IDs
    """
    print("\n" + "="*70)
    print("Example 8: Batch Processing")
    print("="*70)

    adapter = RNAcentralAdapter()

    # List of RNA IDs to process (example IDs)
    rna_ids = [
        "URS0000527F89_9606",  # let-7a
        "URS0000598F88_9606",  # example ID
        "URS00005A37DD_9606",  # example ID
    ]

    print(f"\nProcessing {len(rna_ids)} RNA IDs...\n")

    results = []
    for rna_id in rna_ids:
        result = await adapter.execute(rna_id)

        if result.success:
            rna = result.data
            print(f"✓ {rna_id}")
            print(f"  Type: {rna['rna_type']}, Length: {rna['length']} nt")
            results.append(rna)
        else:
            print(f"✗ {rna_id}: {result.error}")

    print(f"\nSuccessfully processed {len(results)}/{len(rna_ids)} RNAs")

    # Calculate statistics
    if results:
        avg_length = sum(r['length'] for r in results) / len(results)
        print(f"Average length: {avg_length:.1f} nucleotides")


async def example_9_caching_demonstration():
    """
    Example 9: Demonstrate caching behavior
    """
    print("\n" + "="*70)
    print("Example 9: Caching Demonstration")
    print("="*70)

    adapter = RNAcentralAdapter()

    query = "cancer"

    print("\nFirst query (will fetch from API)...")
    import time
    start = time.time()
    result1 = await adapter.execute(query)
    time1 = time.time() - start
    print(f"Time: {time1:.2f}s, Cache hit: {result1.cache_hit}")

    print("\nSecond query (should be cached)...")
    start = time.time()
    result2 = await adapter.execute(query)
    time2 = time.time() - start
    print(f"Time: {time2:.2f}s, Cache hit: {result2.cache_hit}")

    print(f"\nSpeedup: {time1/time2:.1f}x faster with cache")


async def main():
    """
    Run all examples
    """
    print("\n" + "#"*70)
    print("# RNAcentral Adapter - Example Usage")
    print("#"*70)

    examples = [
        ("RNA ID Lookup", example_1_rna_id_lookup),
        ("Keyword Search", example_2_keyword_search),
        ("Filter by RNA Type", example_3_filter_by_rna_type),
        ("Filter by Organism", example_4_organism_filter),
        ("Get Cross-References", example_5_get_cross_references),
        ("Antisense Target Identification", example_6_antisense_target_identification),
        ("RNA Type Survey", example_7_rna_type_survey),
        ("Batch Processing", example_8_batch_processing),
        ("Caching Demonstration", example_9_caching_demonstration),
    ]

    print("\nAvailable examples:")
    for i, (name, _) in enumerate(examples, 1):
        print(f"  {i}. {name}")

    print("\nRunning all examples...\n")

    for name, example_func in examples:
        try:
            await example_func()
            await asyncio.sleep(0.5)  # Brief pause between examples
        except Exception as e:
            print(f"\nError in {name}: {e}")
            logger.exception(f"Exception in {name}")

    print("\n" + "#"*70)
    print("# All examples complete!")
    print("#"*70 + "\n")


if __name__ == "__main__":
    # Run the examples
    asyncio.run(main())
