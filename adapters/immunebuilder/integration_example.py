"""
ImmuneBuilder Integration Examples

Demonstrates how to integrate ImmuneBuilder with other PharmForge adapters
for complete antibody discovery and design workflows.
"""

import asyncio
from adapters.immunebuilder import ImmuneBuilderAdapter
from adapters.sabdab import SAbDabAdapter
from adapters.alphafold import AlphaFoldAdapter


async def example_1_sabdab_to_immunebuilder():
    """
    Example 1: Get antibody from SAbDab and predict structure with ImmuneBuilder.

    Workflow:
    1. Query SAbDab for known antibody
    2. Extract sequences
    3. Predict structure with ImmuneBuilder
    4. Compare to experimental structure
    """
    print("\n" + "="*80)
    print("Example 1: SAbDab → ImmuneBuilder Workflow")
    print("="*80)

    # Initialize adapters
    sabdab = SAbDabAdapter()
    immunebuilder = ImmuneBuilderAdapter()

    # Query SAbDab for specific antibody
    pdb_id = "7BWJ"  # COVID-19 antibody example
    print(f"\nStep 1: Querying SAbDab for PDB {pdb_id}...")

    sabdab_result = await sabdab.execute(pdb_id)

    if not sabdab_result.success:
        print(f"✗ SAbDab query failed: {sabdab_result.error}")
        print("Note: This is expected if adapter needs configuration")
        return

    print(f"✓ Found antibody in SAbDab")

    # Extract sequences (in real implementation, these would come from SAbDab)
    # For now, using example sequences
    heavy_chain = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS"
    light_chain = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK"

    print(f"\nStep 2: Predicting structure with ImmuneBuilder...")
    pred_result = await immunebuilder.execute({
        "heavy_chain": heavy_chain,
        "light_chain": light_chain,
        "name": f"{pdb_id}_predicted"
    })

    if pred_result.success:
        print(f"✓ Structure prediction successful!")
        print(f"  Prediction time: {pred_result.data['prediction_time_seconds']}s")
        print(f"  Mean confidence: {pred_result.data['confidence_scores'].get('mean', 'N/A')}")
        print(f"  Output file: {pred_result.data['file_path']}")
    else:
        print(f"✗ Prediction failed: {pred_result.error}")


async def example_2_batch_antibody_library():
    """
    Example 2: Predict structures for an antibody library.

    Use case: High-throughput screening of antibody variants.
    """
    print("\n" + "="*80)
    print("Example 2: High-Throughput Antibody Library Screening")
    print("="*80)

    immunebuilder = ImmuneBuilderAdapter()

    # Simulate antibody library (in practice, this would come from lab)
    antibody_library = [
        {
            "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
            "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
            "name": "variant_1"
        },
        {
            "heavy_chain": "QVQLQESGPGLVKPSETLSLTCTVSGGSISSSSYYWGWIRQPPGKGLEWIGSIYYSGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCARDLGPYYYYYGMDVWGQGTTVTVSS",
            "light_chain": "DIVMTQSPLSLPVTLGQPASISCRSSQSLVHSNGNTYLHWYLQKPGQSPQLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCFQGSHVPFTFGQGTKLEIK",
            "name": "variant_2"
        },
        {
            "heavy_chain": "EVQLVESGGGLVQAGGSLRLSCAASGRTFSTYAMGWFRQAPGKEREFVARITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPFDYWGQGTQVTVSS",
            "antibody_type": "vhh",
            "name": "nanobody_1"
        }
    ]

    print(f"\nProcessing library of {len(antibody_library)} antibodies...")

    results = await immunebuilder.batch_predict(antibody_library)

    # Analyze results
    successful = [r for r in results if r.success]
    failed = [r for r in results if not r.success]

    print(f"\n✓ Successfully predicted: {len(successful)}/{len(results)}")
    if failed:
        print(f"✗ Failed: {len(failed)}")

    # Rank by confidence
    if successful:
        print("\nRanking by confidence score:")
        sorted_results = sorted(
            successful,
            key=lambda r: r.data['confidence_scores'].get('mean', 0),
            reverse=True
        )

        for i, result in enumerate(sorted_results[:5], 1):
            name = result.data['name']
            conf = result.data['confidence_scores'].get('mean', 'N/A')
            ab_type = result.data['antibody_type']
            print(f"  {i}. {name} ({ab_type}): confidence={conf}")


async def example_3_compare_immunebuilder_alphafold():
    """
    Example 3: Compare ImmuneBuilder vs AlphaFold predictions.

    Demonstrates speed/accuracy tradeoffs.
    """
    print("\n" + "="*80)
    print("Example 3: ImmuneBuilder vs AlphaFold Comparison")
    print("="*80)

    immunebuilder = ImmuneBuilderAdapter()
    alphafold = AlphaFoldAdapter()

    # Example antibody
    antibody_data = {
        "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
        "name": "comparison_ab"
    }

    print("\n1. ImmuneBuilder Prediction (Fast):")
    import time
    start = time.time()
    ib_result = await immunebuilder.execute(antibody_data)
    ib_time = time.time() - start

    if ib_result.success:
        print(f"   ✓ Time: {ib_time:.2f}s")
        print(f"   ✓ Confidence: {ib_result.data['confidence_scores'].get('mean', 'N/A')}")
    else:
        print(f"   ✗ Failed: {ib_result.error}")

    print("\n2. AlphaFold Prediction (Accurate):")
    print("   Note: AlphaFold requires UniProt ID, not raw sequences")
    print("   For antibody sequences, use ImmuneBuilder or upload to AlphaFold server")
    print("   Typical AlphaFold time: 5-30 minutes per structure")

    print("\nRecommendations:")
    print("  • Use ImmuneBuilder for:")
    print("    - High-throughput screening")
    print("    - Rapid prototyping")
    print("    - Initial structure generation")
    print("  • Use AlphaFold for:")
    print("    - Final structure validation")
    print("    - Complex antibody-antigen modeling")
    print("    - Publication-quality structures")


async def example_4_structure_based_design():
    """
    Example 4: Structure-based antibody design workflow.

    Complete workflow for antibody optimization:
    1. Predict structure
    2. Analyze binding interface
    3. Design variants
    4. Predict variant structures
    5. Rank by predicted properties
    """
    print("\n" + "="*80)
    print("Example 4: Structure-Based Antibody Design Workflow")
    print("="*80)

    immunebuilder = ImmuneBuilderAdapter()

    # Step 1: Predict parent antibody structure
    print("\nStep 1: Predict parent antibody structure")
    parent_ab = {
        "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
        "name": "parent"
    }

    parent_result = await immunebuilder.execute(parent_ab)

    if not parent_result.success:
        print(f"✗ Parent prediction failed: {parent_result.error}")
        return

    print(f"✓ Parent structure predicted")
    print(f"  Confidence: {parent_result.data['confidence_scores'].get('mean', 'N/A')}")

    # Step 2: Design variants (simulated - in practice use ML or rational design)
    print("\nStep 2: Design CDR variants")
    variants = [
        {
            **parent_ab,
            "heavy_chain": parent_ab["heavy_chain"][:100] + "MUTATED" + parent_ab["heavy_chain"][107:],
            "name": "variant_CDR_H3_1"
        },
        {
            **parent_ab,
            "heavy_chain": parent_ab["heavy_chain"][:100] + "CHANGED" + parent_ab["heavy_chain"][107:],
            "name": "variant_CDR_H3_2"
        }
    ]
    print(f"  Generated {len(variants)} CDR-H3 variants")

    # Step 3: Predict variant structures
    print("\nStep 3: Predict variant structures")
    variant_results = await immunebuilder.batch_predict(variants)

    successful_variants = [r for r in variant_results if r.success]
    print(f"  ✓ Successfully predicted {len(successful_variants)} variants")

    # Step 4: Rank variants
    print("\nStep 4: Rank variants by confidence")
    all_results = [parent_result] + successful_variants

    ranked = sorted(
        all_results,
        key=lambda r: r.data['confidence_scores'].get('mean', 0),
        reverse=True
    )

    for i, result in enumerate(ranked, 1):
        name = result.data['name']
        conf = result.data['confidence_scores'].get('mean', 'N/A')
        print(f"  {i}. {name}: confidence={conf}")

    print("\nNext steps would include:")
    print("  • Docking to target antigen")
    print("  • Interface analysis")
    print("  • ADMET prediction")
    print("  • Experimental validation")


async def example_5_nanobody_discovery():
    """
    Example 5: VHH nanobody discovery pipeline.

    Demonstrates nanobody-specific workflow.
    """
    print("\n" + "="*80)
    print("Example 5: VHH Nanobody Discovery Pipeline")
    print("="*80)

    immunebuilder = ImmuneBuilderAdapter()

    # Nanobody library (from immunization campaign)
    nanobody_sequences = [
        "EVQLVESGGGLVQAGGSLRLSCAASGRTFSTYAMGWFRQAPGKEREFVARITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPFDYWGQGTQVTVSS",
        "QVQLQESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKGREFVAAITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWFGPFDYWGQGTQVTVSS",
        "EVQLVESGGGLVQAGGSLRLSCAASGRTFSDYAMGWFRQAPGKEREFVARITWSGGGTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPYDYWGQGTQVTVSS"
    ]

    print(f"\nProcessing {len(nanobody_sequences)} nanobodies...")

    nanobodies = [
        {
            "heavy_chain": seq,
            "antibody_type": "vhh",
            "name": f"nanobody_{i+1}"
        }
        for i, seq in enumerate(nanobody_sequences)
    ]

    results = await immunebuilder.batch_predict(nanobodies)

    successful = [r for r in results if r.success]

    if successful:
        print(f"\n✓ Successfully predicted {len(successful)} nanobody structures")

        # Analyze nanobody-specific features
        print("\nNanobody Analysis:")
        for result in successful:
            name = result.data['name']
            conf = result.data['confidence_scores'].get('mean', 'N/A')
            length = result.data['sequences']['heavy_chain_length']

            print(f"  {name}:")
            print(f"    Length: {length} aa")
            print(f"    Confidence: {conf}")
            print(f"    Single-domain: Yes (VHH)")

        print("\nNanobody advantages:")
        print("  • Smaller size (~15 kDa)")
        print("  • Higher stability")
        print("  • Easier production")
        print("  • Better tissue penetration")


async def main():
    """Run all integration examples."""
    print("\n" + "="*80)
    print("ImmuneBuilder Integration Examples")
    print("="*80)

    try:
        await example_1_sabdab_to_immunebuilder()
        await example_2_batch_antibody_library()
        await example_3_compare_immunebuilder_alphafold()
        await example_4_structure_based_design()
        await example_5_nanobody_discovery()

        print("\n" + "="*80)
        print("All integration examples completed!")
        print("="*80)

    except ImportError as e:
        print(f"\n✗ Import Error: {e}")
        print("\nSome adapters may not be installed.")
        print("Install with: pip install immunebuilder")

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    asyncio.run(main())
