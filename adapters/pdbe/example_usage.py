"""
PDBe Adapter - Example Usage

This script demonstrates various use cases for the PDBe adapter.
"""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from adapters.pdbe.adapter import PDBEAdapter


async def example_1_basic_retrieval():
    """Example 1: Basic structure retrieval by PDB ID"""
    print("\n" + "=" * 70)
    print("Example 1: Basic Structure Retrieval")
    print("=" * 70)

    adapter = PDBEAdapter()

    # Fetch COVID-19 main protease structure
    pdb_id = "6LU7"
    print(f"\nFetching structure: {pdb_id}")

    result = await adapter.execute(pdb_id)

    if result.success:
        data = result.data
        print(f"\n[OK] Successfully retrieved: {data['pdb_id']}")
        print(f"  Title: {data['quality_metrics'].get('title', 'N/A')[:80]}...")
        print(f"  Method: {data['experimental_method']}")
        print(f"  Resolution: {data['resolution']} Å")
        print(f"  Organism: {data['organism']}")
        print(f"  Number of ligands: {len(data['ligands'])}")
        print(f"  Number of binding sites: {len(data['binding_sites'])}")

        if 'pdb' in data['structure_files']:
            print(f"  PDB file size: {len(data['structure_files']['pdb'])} bytes")
    else:
        print(f"[ERROR] Error: {result.error}")


async def example_2_search_structures():
    """Example 2: Search for structures by keyword"""
    print("\n" + "=" * 70)
    print("Example 2: Search Structures by Keyword")
    print("=" * 70)

    adapter = PDBEAdapter()

    # Search for kinase structures
    query = {
        "query_type": "keyword",
        "query": "SARS-CoV-2",
        "max_results": 5
    }

    print(f"\nSearching for: {query['query']}")

    result = await adapter.execute(query)

    if result.success:
        data = result.data
        print(f"\n[OK] Found {data['num_results']} results (out of {data['total_available']} total)")

        for i, entry in enumerate(data['entries'], 1):
            print(f"\n{i}. {entry['pdb_id']}: {entry['title'][:60]}...")
            print(f"   Method: {entry['method']}")
            if entry['resolution']:
                print(f"   Resolution: {entry['resolution']} Å")
            print(f"   Organism: {entry['organism']}")
    else:
        print(f"[ERROR] Error: {result.error}")


async def example_3_ligand_analysis():
    """Example 3: Analyze ligands and binding sites"""
    print("\n" + "=" * 70)
    print("Example 3: Ligand and Binding Site Analysis")
    print("=" * 70)

    adapter = PDBEAdapter(config={
        "include_ligands": True,
        "include_binding_sites": True
    })

    # HIV-1 protease with inhibitor
    pdb_id = "1HSG"
    print(f"\nAnalyzing: {pdb_id} (HIV-1 protease with inhibitor)")

    result = await adapter.execute(pdb_id)

    if result.success:
        data = result.data

        print(f"\n[OK] Retrieved: {data['pdb_id']}")
        print(f"  Title: {data['quality_metrics'].get('title', 'N/A')[:80]}...")

        print(f"\n  Ligands ({len(data['ligands'])}):")
        for ligand in data['ligands']:
            print(f"    - {ligand['id']}: {ligand['name'][:60]}")
            print(f"      Formula: {ligand['formula']}")
            if ligand.get('weight'):
                print(f"      Molecular Weight: {ligand['weight']:.2f}")

        print(f"\n  Binding Sites ({len(data['binding_sites'])}):")
        for site in data['binding_sites']:
            print(f"    - Site {site['site_id']}")
            print(f"      Ligand: {site['ligand_id']}")
            print(f"      Residues: {site['num_residues']}")
    else:
        print(f"[ERROR] Error: {result.error}")


async def example_4_download_files():
    """Example 4: Download and cache structure files"""
    print("\n" + "=" * 70)
    print("Example 4: Download and Cache Structure Files")
    print("=" * 70)

    adapter = PDBEAdapter(config={
        "download_pdb": True,
        "download_cif": True,
        "cache_dir": "./test_cache/pdbe"
    })

    pdb_id = "3CL0"
    print(f"\nDownloading structure files for: {pdb_id}")

    result = await adapter.execute(pdb_id)

    if result.success:
        data = result.data

        print(f"\n[OK] Files downloaded and cached:")

        if 'pdb' in data['file_paths']:
            pdb_path = Path(data['file_paths']['pdb'])
            print(f"  PDB file: {pdb_path}")
            print(f"    Size: {pdb_path.stat().st_size} bytes")
            print(f"    Exists: {pdb_path.exists()}")

        if 'cif' in data['file_paths']:
            cif_path = Path(data['file_paths']['cif'])
            print(f"  mmCIF file: {cif_path}")
            print(f"    Size: {cif_path.stat().st_size} bytes")
            print(f"    Exists: {cif_path.exists()}")
    else:
        print(f"[ERROR] Error: {result.error}")


async def example_5_quality_metrics():
    """Example 5: Extract and analyze quality metrics"""
    print("\n" + "=" * 70)
    print("Example 5: Quality Metrics Analysis")
    print("=" * 70)

    adapter = PDBEAdapter()

    # Compare multiple structures
    pdb_ids = ["6LU7", "1HSG", "3CL0"]

    print("\nComparing quality metrics:")

    for pdb_id in pdb_ids:
        result = await adapter.execute(pdb_id)

        if result.success:
            data = result.data
            metrics = data['quality_metrics']

            print(f"\n{pdb_id}:")
            print(f"  Method: {metrics.get('experimental_method', 'N/A')}")
            print(f"  Resolution: {metrics.get('resolution', 'N/A')} Å")
            print(f"  R-value: {metrics.get('r_value', 'N/A')}")
            print(f"  R-free: {metrics.get('r_free', 'N/A')}")
            print(f"  Deposit date: {metrics.get('deposit_date', 'N/A')}")

            # Quality assessment
            resolution = metrics.get('resolution')
            if resolution:
                if resolution < 2.0:
                    quality = "HIGH"
                elif resolution < 3.0:
                    quality = "MEDIUM"
                else:
                    quality = "LOW"
                print(f"  Quality: {quality}")


async def example_6_caching_demo():
    """Example 6: Demonstrate caching behavior"""
    print("\n" + "=" * 70)
    print("Example 6: Caching Demonstration")
    print("=" * 70)

    adapter = PDBEAdapter()

    pdb_id = "6LU7"

    print(f"\nFirst request for {pdb_id} (will fetch from API):")
    result1 = await adapter.execute(pdb_id)
    cache_hit_1 = result1.cache_hit

    print(f"  Cache hit: {cache_hit_1}")
    print(f"  Success: {result1.success}")

    print(f"\nSecond request for {pdb_id} (should use cache):")
    result2 = await adapter.execute(pdb_id)
    cache_hit_2 = result2.cache_hit

    print(f"  Cache hit: {cache_hit_2}")
    print(f"  Success: {result2.success}")

    if cache_hit_2:
        print("\n[OK] Caching is working correctly!")


async def example_7_error_handling():
    """Example 7: Error handling"""
    print("\n" + "=" * 70)
    print("Example 7: Error Handling")
    print("=" * 70)

    adapter = PDBEAdapter()

    # Try invalid PDB ID
    invalid_id = "XXXX"
    print(f"\nTrying invalid PDB ID: {invalid_id}")

    result = await adapter.execute(invalid_id)

    if result.success:
        print("[OK] Success (unexpected)")
    else:
        print(f"[ERROR] Expected error: {result.error}")

    # Try non-existent PDB ID
    nonexistent_id = "9ZZZ"
    print(f"\nTrying non-existent PDB ID: {nonexistent_id}")

    result = await adapter.execute(nonexistent_id)

    if result.success:
        print("[OK] Success (unexpected)")
    else:
        print(f"[ERROR] Expected error: {result.error}")


async def example_8_metadata():
    """Example 8: Get adapter metadata"""
    print("\n" + "=" * 70)
    print("Example 8: Adapter Metadata")
    print("=" * 70)

    adapter = PDBEAdapter()

    metadata = adapter.get_metadata()

    print(f"\nAdapter: {metadata['name']}")
    print(f"Type: {metadata['type']}")
    print(f"Version: {metadata['version']}")
    print(f"Description: {metadata['description']}")

    print("\nCapabilities:")
    for capability, enabled in metadata['capabilities'].items():
        print(f"  - {capability}: {enabled}")

    print("\nDatabase Information:")
    db_info = metadata['database']
    print(f"  Name: {db_info['name']}")
    print(f"  Full Name: {db_info['full_name']}")
    print(f"  Location: {db_info['location']}")
    print(f"  Structures: {db_info['structures']}")

    print("\nAdvantages:")
    for advantage in db_info['advantages']:
        print(f"  - {advantage}")


async def main():
    """Run all examples"""
    print("\n" + "=" * 70)
    print("PDBe Adapter - Example Usage")
    print("=" * 70)

    examples = [
        ("Basic Retrieval", example_1_basic_retrieval),
        ("Search Structures", example_2_search_structures),
        ("Ligand Analysis", example_3_ligand_analysis),
        ("Download Files", example_4_download_files),
        ("Quality Metrics", example_5_quality_metrics),
        ("Caching Demo", example_6_caching_demo),
        ("Error Handling", example_7_error_handling),
        ("Adapter Metadata", example_8_metadata),
    ]

    for name, example_func in examples:
        try:
            await example_func()
        except Exception as e:
            print(f"\n[ERROR] Example '{name}' failed: {e}")

    print("\n" + "=" * 70)
    print("All examples completed!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
