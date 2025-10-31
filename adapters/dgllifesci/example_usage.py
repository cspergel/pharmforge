"""
Example usage of the DGL-LifeSci adapter
This file demonstrates how to use the adapter - DO NOT RUN without installing dgllife first
"""
import asyncio
from adapters.dgllifesci.adapter import DGLLifeSciAdapter


async def example_basic_usage():
    """
    Example 1: Basic molecular graph generation
    """
    print("=" * 80)
    print("Example 1: Basic Molecular Graph Generation")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()

    # Generate molecular graph for ethanol
    smiles = "CCO"
    result = await adapter.execute(smiles, featurizer_type="canonical")

    if result.success:
        data = result.data
        print(f"\nSMILES: {smiles}")
        print(f"Featurizer: {data['featurizer_type']}")
        print(f"Number of nodes: {data['graph_features']['num_nodes']}")
        print(f"Number of edges: {data['graph_features']['num_edges']}")
        print(f"Node feature dimension: {data['graph_features']['node_features']['feature_dim']}")
        print(f"Edge feature dimension: {data['graph_features']['edge_features']['feature_dim']}")
    else:
        print(f"Error: {result.error}")


async def example_multiple_featurizers():
    """
    Example 2: Compare different featurizers
    """
    print("\n" + "=" * 80)
    print("Example 2: Compare Different Featurizers")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()
    smiles = "c1ccccc1"  # Benzene

    featurizers = ["canonical", "attentivefp", "weave", "minimal"]

    for featurizer_type in featurizers:
        result = await adapter.execute(smiles, featurizer_type=featurizer_type)

        if result.success:
            data = result.data
            node_dim = data['graph_features']['node_features']['feature_dim']
            edge_dim = data['graph_features']['edge_features']['feature_dim']
            print(f"\n{featurizer_type:15s} -> Node features: {node_dim:3d}, Edge features: {edge_dim:3d}")
        else:
            print(f"\n{featurizer_type:15s} -> Error: {result.error}")


async def example_batch_processing():
    """
    Example 3: Batch processing of multiple molecules
    """
    print("\n" + "=" * 80)
    print("Example 3: Batch Processing")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()

    # Process multiple drug-like molecules
    smiles_list = [
        "CCO",                              # Ethanol
        "c1ccccc1",                         # Benzene
        "CC(=O)O",                          # Acetic acid
        "CC(C)CC1=CC=C(C=C1)C(C)C",        # Ibuprofen
        "CC(C)NCC(COC1=CC=CC=C1)O",        # Propranolol
    ]

    result = await adapter.execute(smiles_list, featurizer_type="canonical")

    if result.success:
        data = result.data
        print(f"\nTotal molecules: {data['n_molecules']}")
        print(f"Successfully processed: {data['n_successful']}")
        print(f"Failed: {data['n_failed']}")

        print("\nGraph Statistics:")
        print(f"{'SMILES':<40s} {'Nodes':<8s} {'Edges':<8s}")
        print("-" * 60)

        for features in data['graph_features']:
            smiles = features['smiles']
            nodes = features['num_nodes']
            edges = features['num_edges']
            print(f"{smiles:<40s} {nodes:<8d} {edges:<8d}")
    else:
        print(f"Error: {result.error}")


async def example_with_raw_graphs():
    """
    Example 4: Get raw DGL graphs for further processing
    """
    print("\n" + "=" * 80)
    print("Example 4: Access Raw DGL Graphs")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()
    smiles = "CCO"

    result = await adapter.execute(
        smiles,
        featurizer_type="canonical",
        return_graphs=True
    )

    if result.success:
        graph = result.data['graphs']
        print(f"\nSMILES: {smiles}")
        print(f"Graph object type: {type(graph)}")
        print(f"Number of nodes: {graph.num_nodes()}")
        print(f"Number of edges: {graph.num_edges()}")

        # Access node features
        if 'h' in graph.ndata:
            node_features = graph.ndata['h']
            print(f"Node features shape: {node_features.shape}")
            print(f"Node features dtype: {node_features.dtype}")

        # Access edge features
        if 'e' in graph.edata:
            edge_features = graph.edata['e']
            print(f"Edge features shape: {edge_features.shape}")
            print(f"Edge features dtype: {edge_features.dtype}")
    else:
        print(f"Error: {result.error}")


async def example_model_predictions():
    """
    Example 5: Load model and make predictions
    """
    print("\n" + "=" * 80)
    print("Example 5: Model Predictions (Requires Pre-trained Model)")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()

    # Note: This example shows the API, but requires an actual model file
    model_path = "path/to/gcn_model.pt"
    model_type = "gcn"

    print(f"\nAttempting to load model: {model_path}")
    print(f"Model type: {model_type}")

    # Load model (will fail without actual model file)
    success = adapter.load_model(model_path, model_type=model_type)

    if success:
        print("Model loaded successfully!")

        # Make predictions
        smiles_list = ["CCO", "c1ccccc1", "CC(=O)O"]

        result = await adapter.execute(
            smiles_list,
            featurizer_type="canonical",
            include_predictions=True
        )

        if result.success:
            predictions = result.data['predictions']
            print(f"\nPredictions available: {predictions['predictions_available']}")
            print(f"Model type: {predictions['model_type']}")
            print(f"Number of molecules: {predictions['n_molecules']}")
        else:
            print(f"Prediction error: {result.error}")
    else:
        print("Model loading failed (expected - no model file provided)")


async def example_error_handling():
    """
    Example 6: Error handling and validation
    """
    print("\n" + "=" * 80)
    print("Example 6: Error Handling")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()

    # Test with invalid SMILES
    invalid_smiles = [
        "CCO",           # Valid
        "INVALID",       # Invalid
        "c1ccccc1",      # Valid
        "XYZ123",        # Invalid
    ]

    print("\nProcessing mixture of valid and invalid SMILES:")
    result = await adapter.execute(invalid_smiles, featurizer_type="canonical")

    if result.success:
        data = result.data
        print(f"Successfully processed: {data['n_successful']}/{data['n_molecules']}")

        if data['n_failed'] > 0:
            print(f"\nFailed SMILES:")
            for smiles in data.get('failed_smiles', []):
                print(f"  - {smiles}")

        print(f"\nSuccessfully processed:")
        for features in data['graph_features']:
            print(f"  - {features['smiles']}")
    else:
        print(f"Complete failure: {result.error}")


async def example_metadata():
    """
    Example 7: Accessing adapter metadata
    """
    print("\n" + "=" * 80)
    print("Example 7: Adapter Metadata")
    print("=" * 80)

    adapter = DGLLifeSciAdapter()

    # Get adapter metadata
    metadata = adapter.get_metadata()

    print("\nAdapter Information:")
    print(f"Name: {metadata['name']}")
    print(f"Type: {metadata['type']}")
    print(f"Version: {metadata['version']}")
    print(f"Enabled: {metadata['enabled']}")

    print("\nConfiguration:")
    for key, value in metadata['config'].items():
        print(f"  {key}: {value}")

    # Available featurizers
    print("\nSupported Featurizers:")
    for name, config in adapter.FEATURIZER_TYPES.items():
        print(f"  - {name}: {config['description']}")

    # Available models
    print("\nSupported Pre-trained Models:")
    for name, info in adapter.PRETRAINED_MODELS.items():
        print(f"  - {name}: {info['description']}")
        print(f"    Tasks: {', '.join(info['tasks'])}")


async def main():
    """
    Run all examples
    """
    print("\n")
    print("=" * 80)
    print(" DGL-LifeSci Adapter - Usage Examples")
    print("=" * 80)
    print("\nNOTE: These examples require dgllife to be installed:")
    print("      pip install dgllife")
    print("\n")

    try:
        # Run examples
        await example_basic_usage()
        await example_multiple_featurizers()
        await example_batch_processing()
        await example_with_raw_graphs()
        await example_model_predictions()
        await example_error_handling()
        await example_metadata()

        print("\n" + "=" * 80)
        print("All examples completed!")
        print("=" * 80 + "\n")

    except ImportError as e:
        print(f"\nImportError: {e}")
        print("\nPlease install required packages:")
        print("  pip install dgllife")
        print("  pip install rdkit")
    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Run examples
    asyncio.run(main())
