"""
SAbDab Adapter - Example Usage

This script demonstrates various ways to use the SAbDab adapter
for antibody structure and sequence retrieval.
"""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path to import adapter
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.sabdab import SAbDabAdapter


async def example_pdb_lookup():
    """Example 1: Look up antibody by PDB ID"""
    print("\n" + "="*70)
    print("Example 1: PDB ID Lookup")
    print("="*70)

    adapter = SAbDabAdapter()

    # Look up a specific antibody structure
    pdb_id = "7BWJ"  # COVID-19 spike antibody
    result = await adapter(pdb_id, include_cdrs=True)

    if result.success:
        data = result.data
        ab_data = data['antibody_data']

        print(f"\n✓ Successfully retrieved: {data['pdb_id']}")
        print(f"\nAntibody Information:")
        print(f"  Resolution:      {ab_data['resolution']} Å")
        print(f"  Method:          {ab_data['method']}")
        print(f"  R-factor:        {ab_data['rfactor']}")
        print(f"  Antibody Type:   {ab_data['antibody_type']}")
        print(f"  Antigen:         {ab_data['antigen']}")
        print(f"  Species:         {ab_data['species']}")
        print(f"  Heavy Chain:     {ab_data['heavy_chain']}")
        print(f"  Light Chain:     {ab_data['light_chain']}")
        print(f"  Engineered:      {ab_data['engineered']}")

        # CDR information
        if 'cdrs' in ab_data:
            print(f"\nCDR Information:")
            cdrs = ab_data['cdrs']
            print(f"  Heavy Chain CDRs: {list(cdrs.get('heavy_chain', {}).keys())}")
            print(f"  Light Chain CDRs: {list(cdrs.get('light_chain', {}).keys())}")
    else:
        print(f"✗ Error: {result.error}")


async def example_antigen_search():
    """Example 2: Search antibodies by antigen name"""
    print("\n" + "="*70)
    print("Example 2: Search by Antigen")
    print("="*70)

    adapter = SAbDabAdapter()

    # Search for COVID-19 antibodies
    search_params = {
        "antigen": "spike",
        "species": "human",
        "resolution_max": 3.0
    }

    result = await adapter(search_params, max_results=10)

    if result.success:
        data = result.data
        antibodies = data['antibodies']

        print(f"\n✓ Found {data['total_results']} antibodies")
        print(f"\nSearch Parameters:")
        for key, value in data['search_params'].items():
            print(f"  {key}: {value}")

        print(f"\nTop Results:")
        for i, ab in enumerate(antibodies[:5], 1):
            print(f"\n  {i}. {ab['pdb_id']}")
            print(f"     Resolution: {ab['resolution']} Å")
            print(f"     Type: {ab['antibody_type']}")
            print(f"     Antigen: {ab['antigen']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_nanobody_search():
    """Example 3: Search for nanobodies (VHH)"""
    print("\n" + "="*70)
    print("Example 3: Nanobody Search")
    print("="*70)

    adapter = SAbDabAdapter()

    # Search for high-quality nanobodies
    search_params = {
        "ab_type": "VHH",  # Single-domain antibodies
        "resolution_max": 2.5,
        "rfactor_max": 0.25
    }

    result = await adapter(search_params, max_results=20)

    if result.success:
        data = result.data
        nanobodies = data['antibodies']

        print(f"\n✓ Found {len(nanobodies)} high-quality nanobodies")

        print(f"\nNanobody Structures:")
        for i, nb in enumerate(nanobodies[:5], 1):
            print(f"\n  {i}. {nb['pdb_id']}")
            print(f"     Resolution: {nb['resolution']} Å")
            print(f"     Species: {nb['species']}")
            print(f"     Antigen: {nb['antigen']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_quality_filtering():
    """Example 4: Find high-quality structures"""
    print("\n" + "="*70)
    print("Example 4: Quality Filtering")
    print("="*70)

    adapter = SAbDabAdapter()

    # Search for very high-quality structures
    search_params = {
        "antigen": "Her2",
        "species": "human",
        "resolution_max": 2.0,  # Very high resolution
        "rfactor_max": 0.20     # Excellent R-factor
    }

    result = await adapter(search_params)

    if result.success:
        data = result.data
        antibodies = data['antibodies']

        print(f"\n✓ Found {len(antibodies)} high-quality structures")

        if antibodies:
            print(f"\nQuality Metrics:")
            for ab in antibodies[:3]:
                print(f"\n  {ab['pdb_id']}:")
                print(f"    Resolution: {ab['resolution']} Å")
                print(f"    Type: {ab['antibody_type']}")
                print(f"    Antigen: {ab['antigen']}")
        else:
            print("\n  No structures meet these strict quality criteria")
    else:
        print(f"✗ Error: {result.error}")


async def example_antibody_type_comparison():
    """Example 5: Compare different antibody formats"""
    print("\n" + "="*70)
    print("Example 5: Antibody Format Comparison")
    print("="*70)

    adapter = SAbDabAdapter()

    antibody_types = ["IgG", "Fab", "scFv", "VHH"]

    print(f"\nComparing antibody formats against same target:")

    for ab_type in antibody_types:
        search_params = {
            "ab_type": ab_type,
            "antigen": "SARS-CoV-2",
            "resolution_max": 3.5
        }

        result = await adapter(search_params, max_results=50)

        if result.success:
            count = result.data['total_results']
            antibodies = result.data['antibodies']

            # Calculate average resolution
            resolutions = [ab['resolution'] for ab in antibodies if ab['resolution']]
            avg_res = sum(resolutions) / len(resolutions) if resolutions else 0

            print(f"\n  {ab_type:6s}: {count:3d} structures, avg resolution: {avg_res:.2f} Å")
        else:
            print(f"\n  {ab_type:6s}: Error - {result.error}")


async def example_batch_lookup():
    """Example 6: Batch lookup of multiple PDB IDs"""
    print("\n" + "="*70)
    print("Example 6: Batch PDB Lookup")
    print("="*70)

    adapter = SAbDabAdapter()

    # Multiple antibody structures to look up
    pdb_ids = ["7BWJ", "6XDG", "7KMG", "7C01", "6M0J"]

    print(f"\nLooking up {len(pdb_ids)} antibody structures...")

    results = []
    for pdb_id in pdb_ids:
        result = await adapter(pdb_id)
        results.append((pdb_id, result))

    print(f"\nBatch Results:")
    print(f"{'PDB ID':<8} {'Status':<10} {'Resolution':<12} {'Antigen':<30}")
    print("-" * 70)

    for pdb_id, result in results:
        if result.success:
            ab_data = result.data['antibody_data']
            status = "✓ Success"
            resolution = f"{ab_data['resolution']} Å" if ab_data['resolution'] else "N/A"
            antigen = ab_data['antigen'][:28]
        else:
            status = "✗ Failed"
            resolution = "N/A"
            antigen = result.error[:28]

        print(f"{pdb_id:<8} {status:<10} {resolution:<12} {antigen:<30}")


async def example_metadata():
    """Example 7: Get adapter metadata"""
    print("\n" + "="*70)
    print("Example 7: Adapter Metadata")
    print("="*70)

    adapter = SAbDabAdapter()
    metadata = adapter.get_metadata()

    print(f"\nAdapter Information:")
    print(f"  Name:        {metadata['name']}")
    print(f"  Version:     {metadata['version']}")
    print(f"  Type:        {metadata['type']}")
    print(f"  Description: {metadata['description']}")

    print(f"\nCapabilities:")
    for cap, enabled in metadata['capabilities'].items():
        status = "✓" if enabled else "✗"
        print(f"  {status} {cap}")

    print(f"\nDatabase Information:")
    db_info = metadata['database']
    print(f"  Name:       {db_info['name']}")
    print(f"  Maintainer: {db_info['maintainer']}")
    print(f"  URL:        {db_info['url']}")
    print(f"  Updates:    {db_info['update_frequency']}")

    print(f"\nSupported Antibody Types:")
    for ab_type in metadata['antibody_types']:
        print(f"  • {ab_type}")


async def example_caching():
    """Example 8: Demonstrate caching"""
    print("\n" + "="*70)
    print("Example 8: Caching Demonstration")
    print("="*70)

    adapter = SAbDabAdapter()
    pdb_id = "7BWJ"

    print(f"\nFirst request (cache miss):")
    import time
    start = time.time()
    result1 = await adapter(pdb_id, use_cache=True)
    time1 = time.time() - start

    if result1.success:
        print(f"  ✓ Retrieved in {time1:.2f}s")
        print(f"  Cache hit: {result1.cache_hit}")

    print(f"\nSecond request (cache hit):")
    start = time.time()
    result2 = await adapter(pdb_id, use_cache=True)
    time2 = time.time() - start

    if result2.success:
        print(f"  ✓ Retrieved in {time2:.2f}s")
        print(f"  Cache hit: {result2.cache_hit}")
        print(f"  Speedup: {time1/time2:.1f}x faster")


async def main():
    """Run all examples"""
    print("\n" + "="*70)
    print("SAbDab Adapter - Usage Examples")
    print("="*70)

    examples = [
        example_pdb_lookup,
        example_antigen_search,
        example_nanobody_search,
        example_quality_filtering,
        example_antibody_type_comparison,
        example_batch_lookup,
        example_metadata,
        example_caching
    ]

    for example_func in examples:
        try:
            await example_func()
        except Exception as e:
            print(f"\n✗ Example failed with error: {e}")

    print("\n" + "="*70)
    print("All examples completed!")
    print("="*70 + "\n")


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
