"""
ImmuneBuilder Adapter - Example Usage

Demonstrates various ways to use the ImmuneBuilder adapter for antibody structure prediction.
"""

import asyncio
from adapters.immunebuilder import ImmuneBuilderAdapter


async def example_1_basic_fab_prediction():
    """
    Example 1: Basic Fab structure prediction from heavy and light chains.
    """
    print("\n" + "="*80)
    print("Example 1: Basic Fab Structure Prediction")
    print("="*80)

    adapter = ImmuneBuilderAdapter()

    # Example antibody sequences (truncated for brevity)
    # In practice, use full variable domain sequences
    input_data = {
        "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
        "name": "example_fab"
    }

    result = await adapter.execute(input_data)

    if result.success:
        print(f"✓ Structure prediction successful!")
        print(f"  Antibody name: {result.data['name']}")
        print(f"  Antibody type: {result.data['antibody_type']}")
        print(f"  Prediction time: {result.data['prediction_time_seconds']}s")
        print(f"  Output file: {result.data['file_path']}")
        print(f"\nConfidence Scores:")
        conf = result.data['confidence_scores']
        print(f"  Mean: {conf.get('mean', 'N/A')}")
        print(f"  Median: {conf.get('median', 'N/A')}")
        print(f"  Range: {conf.get('min', 'N/A')} - {conf.get('max', 'N/A')}")
        print(f"\nStructure Stats:")
        stats = result.data['structure_stats']
        print(f"  Residues: {stats['num_residues']}")
        print(f"  Atoms: {stats['num_atoms']}")
        print(f"  Chains: {', '.join(stats['chain_ids'])}")
    else:
        print(f"✗ Prediction failed: {result.error}")


async def example_2_vhh_nanobody():
    """
    Example 2: VHH nanobody prediction (heavy chain only).
    """
    print("\n" + "="*80)
    print("Example 2: VHH Nanobody Structure Prediction")
    print("="*80)

    adapter = ImmuneBuilderAdapter()

    # VHH nanobody sequence (no light chain)
    input_data = {
        "heavy_chain": "EVQLVESGGGLVQAGGSLRLSCAASGRTFSTYAMGWFRQAPGKEREFVARITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPFDYWGQGTQVTVSS",
        "antibody_type": "vhh",  # Explicitly specify VHH
        "name": "nanobody_1"
    }

    result = await adapter.execute(input_data)

    if result.success:
        print(f"✓ Nanobody structure prediction successful!")
        print(f"  Name: {result.data['name']}")
        print(f"  Type: {result.data['antibody_type']}")
        print(f"  Heavy chain length: {result.data['sequences']['heavy_chain_length']} aa")
        print(f"  Light chain: {result.data['sequences']['light_chain'] or 'None (VHH)'}")
        print(f"  File: {result.data['file_path']}")
    else:
        print(f"✗ Prediction failed: {result.error}")


async def example_3_batch_prediction():
    """
    Example 3: Batch prediction for multiple antibodies.
    """
    print("\n" + "="*80)
    print("Example 3: Batch Structure Prediction")
    print("="*80)

    adapter = ImmuneBuilderAdapter()

    # Multiple antibodies to predict
    antibodies = [
        {
            "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
            "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
            "name": "antibody_1"
        },
        {
            "heavy_chain": "QVQLQESGPGLVKPSETLSLTCTVSGGSISSSSYYWGWIRQPPGKGLEWIGSIYYSGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCARDLGPYYYYYGMDVWGQGTTVTVSS",
            "light_chain": "DIVMTQSPLSLPVTLGQPASISCRSSQSLVHSNGNTYLHWYLQKPGQSPQLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCFQGSHVPFTFGQGTKLEIK",
            "name": "antibody_2"
        },
        {
            "heavy_chain": "EVQLVESGGGLVQAGGSLRLSCAASGRTFSTYAMGWFRQAPGKEREFVARITWSGGSTYFADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAAGGWLGPFDYWGQGTQVTVSS",
            "antibody_type": "vhh",
            "name": "nanobody_1"
        }
    ]

    results = await adapter.batch_predict(antibodies)

    print(f"\nBatch Prediction Results:")
    print(f"Total antibodies: {len(results)}")
    successful = sum(1 for r in results if r.success)
    print(f"Successful: {successful}")
    print(f"Failed: {len(results) - successful}")

    print("\nDetailed Results:")
    for i, result in enumerate(results, 1):
        if result.success:
            print(f"\n  {i}. {result.data['name']}")
            print(f"     Type: {result.data['antibody_type']}")
            print(f"     Time: {result.data['prediction_time_seconds']}s")
            print(f"     Confidence: {result.data['confidence_scores'].get('mean', 'N/A')}")
        else:
            print(f"\n  {i}. Failed - {result.error}")


async def example_4_custom_config():
    """
    Example 4: Using custom configuration.
    """
    print("\n" + "="*80)
    print("Example 4: Custom Configuration")
    print("="*80)

    # Initialize with custom configuration
    adapter = ImmuneBuilderAdapter(config={
        "output_dir": "./my_antibodies",
        "num_models": 1,
        "save_structures": True,
        "use_gpu": True
    })

    input_data = {
        "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
        "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
        "name": "custom_config_example"
    }

    result = await adapter.execute(input_data, save_structures=True)

    if result.success:
        print(f"✓ Structure saved to: {result.data['file_path']}")
        print(f"  Configuration used:")
        metadata = adapter.get_metadata()
        for key, value in metadata['config'].items():
            print(f"    {key}: {value}")
    else:
        print(f"✗ Prediction failed: {result.error}")


async def example_5_validation():
    """
    Example 5: Input validation and error handling.
    """
    print("\n" + "="*80)
    print("Example 5: Input Validation")
    print("="*80)

    adapter = ImmuneBuilderAdapter()

    # Test various invalid inputs
    test_cases = [
        {
            "name": "Missing heavy chain",
            "input": {"light_chain": "DIQMTQS..."},
            "expected": "Should fail - no heavy chain"
        },
        {
            "name": "Invalid amino acids",
            "input": {
                "heavy_chain": "EVQLVESGGGXYZ123",
                "name": "invalid_aa"
            },
            "expected": "Should fail - invalid characters"
        },
        {
            "name": "Empty sequence",
            "input": {
                "heavy_chain": "",
                "name": "empty"
            },
            "expected": "Should fail - empty sequence"
        },
        {
            "name": "Valid input",
            "input": {
                "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
                "name": "valid"
            },
            "expected": "Should succeed"
        }
    ]

    for test in test_cases:
        print(f"\n  Testing: {test['name']}")
        print(f"  Expected: {test['expected']}")

        is_valid = adapter.validate_input(test['input'])
        print(f"  Validation result: {'✓ Valid' if is_valid else '✗ Invalid'}")


async def example_6_metadata():
    """
    Example 6: Accessing adapter metadata.
    """
    print("\n" + "="*80)
    print("Example 6: Adapter Metadata")
    print("="*80)

    adapter = ImmuneBuilderAdapter()
    metadata = adapter.get_metadata()

    print(f"Adapter Name: {metadata['name']}")
    print(f"Version: {metadata['version']}")
    print(f"Type: {metadata['type']}")
    print(f"Description: {metadata['description']}")

    print(f"\nCapabilities:")
    for cap, enabled in metadata['capabilities'].items():
        print(f"  {cap}: {enabled}")

    print(f"\nSupported Antibody Types:")
    for ab_type, description in metadata['antibody_types'].items():
        print(f"  {ab_type}: {description}")

    print(f"\nMethod Information:")
    for key, value in metadata['method'].items():
        print(f"  {key}: {value}")

    print(f"\nReference:")
    print(f"  Paper: {metadata['reference']['paper']}")
    print(f"  DOI: {metadata['reference']['doi']}")
    print(f"  GitHub: {metadata['reference']['github']}")


async def example_7_integration_with_sabdab():
    """
    Example 7: Integration with SAbDab adapter.
    """
    print("\n" + "="*80)
    print("Example 7: Integration with SAbDab")
    print("="*80)

    print("This example shows how ImmuneBuilder could integrate with SAbDab:")
    print("\n1. Query SAbDab for antibody sequences")
    print("2. Extract heavy/light chain sequences")
    print("3. Predict structure with ImmuneBuilder")
    print("4. Compare predicted vs experimental structure")
    print("\nNote: This requires the SAbDab adapter to be available")

    # Pseudo-code for integration
    print("\nExample workflow:")
    print("""
    from adapters.sabdab import SAbDabAdapter
    from adapters.immunebuilder import ImmuneBuilderAdapter

    # Get sequences from SAbDab
    sabdab = SAbDabAdapter()
    sabdab_result = await sabdab.execute("7BWJ")

    # Extract sequences
    heavy = sabdab_result.data["antibody_data"]["heavy_sequence"]
    light = sabdab_result.data["antibody_data"]["light_sequence"]

    # Predict structure
    immunebuilder = ImmuneBuilderAdapter()
    pred_result = await immunebuilder.execute({
        "heavy_chain": heavy,
        "light_chain": light,
        "name": "7BWJ_predicted"
    })
    """)


async def main():
    """
    Run all examples.
    """
    print("\n" + "="*80)
    print("ImmuneBuilder Adapter - Example Usage")
    print("="*80)
    print("\nRunning examples...")

    try:
        # Run examples
        await example_1_basic_fab_prediction()
        await example_2_vhh_nanobody()
        await example_3_batch_prediction()
        await example_4_custom_config()
        await example_5_validation()
        await example_6_metadata()
        await example_7_integration_with_sabdab()

        print("\n" + "="*80)
        print("All examples completed!")
        print("="*80)

    except ImportError as e:
        print(f"\n✗ Error: {e}")
        print("\nImmuneBuilder may not be installed.")
        print("Install with: pip install immunebuilder")
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
