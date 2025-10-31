"""
TorchDrug Adapter - Usage Examples

This file demonstrates various ways to use the TorchDrug adapter
for molecular feature extraction and embedding generation.
"""

import asyncio
import numpy as np
from adapters.torchdrug import TorchDrugAdapter


async def example_basic_usage():
    """Basic usage: Generate embeddings for a single molecule"""
    print("\n=== Example 1: Basic Usage ===")

    adapter = TorchDrugAdapter()

    # Generate embeddings for ethanol
    result = await adapter.execute("CCO")

    if result.success:
        print(f"✓ Success!")
        print(f"  SMILES: {result.data['smiles']}")
        print(f"  Embedding dimension: {result.data['embedding_dim']}")
        print(f"  Number of atoms: {result.data['graph_features']['num_atoms']}")
        print(f"  Number of bonds: {result.data['graph_features']['num_bonds']}")
        print(f"  Model used: {result.metadata['model_type']}")
        print(f"  Device: {result.metadata['device']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_batch_processing():
    """Batch processing: Generate embeddings for multiple molecules"""
    print("\n=== Example 2: Batch Processing ===")

    adapter = TorchDrugAdapter()

    # List of common organic molecules
    molecules = [
        "CCO",           # Ethanol
        "CC(=O)O",       # Acetic acid
        "c1ccccc1",      # Benzene
        "CC(C)O",        # Isopropanol
        "C1CCCCC1"       # Cyclohexane
    ]

    result = await adapter.execute(molecules)

    if result.success:
        print(f"✓ Success!")
        print(f"  Processed {result.data['n_molecules']} molecules")
        print(f"  Failed: {result.data['n_failed']} molecules")
        print(f"  Embedding shape: {result.data['embedding_shape']}")

        # Display first molecule's features
        first_features = result.data['graph_features'][0]
        print(f"\n  First molecule (CCO):")
        print(f"    Atoms: {first_features['num_atoms']}")
        print(f"    Bonds: {first_features['num_bonds']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_different_models():
    """Using different GNN architectures"""
    print("\n=== Example 3: Different GNN Models ===")

    adapter = TorchDrugAdapter()
    smiles = "c1ccccc1"  # Benzene

    models = ["gin", "gcn", "gat", "graphsage"]

    for model_type in models:
        result = await adapter.execute(smiles, model_type=model_type)

        if result.success:
            print(f"✓ {model_type.upper()}: {result.metadata['model_description']}")
            print(f"  Embedding dim: {result.data['embedding_dim']}")
        else:
            print(f"✗ {model_type.upper()}: {result.error}")


async def example_graph_features_only():
    """Extract graph features without generating embeddings"""
    print("\n=== Example 4: Graph Features Only ===")

    adapter = TorchDrugAdapter()

    result = await adapter.execute(
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        generate_embeddings=False,
        include_graph_features=True
    )

    if result.success:
        features = result.data['graph_features']
        print(f"✓ Aspirin molecular graph features:")
        print(f"  Atoms: {features['num_atoms']}")
        print(f"  Bonds: {features['num_bonds']}")
        print(f"  Atom feature dimension: {features['atom_feature_dim']}")
        print(f"  Bond feature dimension: {features['bond_feature_dim']}")
        print(f"  Unique atom types: {features['unique_atom_types']}")
        print(f"  Unique bond types: {features['unique_bond_types']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_custom_configuration():
    """Using custom model configuration"""
    print("\n=== Example 5: Custom Configuration ===")

    adapter = TorchDrugAdapter()

    result = await adapter.execute(
        "CCO",
        model_type="gin",
        hidden_dims=[128, 128, 128],  # Custom hidden layer dimensions
        batch_size=16,
        generate_embeddings=True,
        include_graph_features=True
    )

    if result.success:
        print(f"✓ Custom GIN configuration:")
        print(f"  Embedding dimension: {result.data['embedding_dim']}")
        print(f"  Embedding stats:")
        stats = result.data['embedding_stats']
        print(f"    Mean: {stats['mean']:.4f}")
        print(f"    Std: {stats['std']:.4f}")
        print(f"    Min: {stats['min']:.4f}")
        print(f"    Max: {stats['max']:.4f}")
    else:
        print(f"✗ Error: {result.error}")


async def example_with_caching():
    """Using the adapter with caching enabled"""
    print("\n=== Example 6: Caching ===")

    adapter = TorchDrugAdapter()
    smiles = "CCO"

    # First call - will compute and cache
    print("First call (computing)...")
    result1 = await adapter(smiles, use_cache=True, model_type="gin")

    if result1.success:
        print(f"  Cache hit: {result1.cache_hit}")
        print(f"  Embedding dim: {result1.data['embedding_dim']}")

    # Second call - should use cache
    print("\nSecond call (from cache)...")
    result2 = await adapter(smiles, use_cache=True, model_type="gin")

    if result2.success:
        print(f"  Cache hit: {result2.cache_hit}")
        print(f"  Embedding dim: {result2.data['embedding_dim']}")


async def example_error_handling():
    """Error handling for invalid input"""
    print("\n=== Example 7: Error Handling ===")

    adapter = TorchDrugAdapter()

    # Invalid SMILES
    result = await adapter.execute("INVALID_SMILES_123")

    if not result.success:
        print(f"✓ Properly handled invalid SMILES")
        print(f"  Error: {result.error}")

    # Mix of valid and invalid SMILES
    result = await adapter.execute(["CCO", "INVALID", "c1ccccc1"])

    if result.success:
        print(f"\n✓ Partial success with mixed input:")
        print(f"  Successful: {result.data['n_molecules']}")
        print(f"  Failed: {result.data['n_failed']}")
        if result.data['n_failed'] > 0:
            print(f"  Failed SMILES: {result.data['failed_smiles']}")


async def example_molecular_similarity():
    """Computing molecular similarity using embeddings"""
    print("\n=== Example 8: Molecular Similarity ===")

    adapter = TorchDrugAdapter()

    # Query molecule
    query = "CCO"  # Ethanol

    # Database molecules
    database = [
        "CC(C)O",      # Isopropanol (similar alcohol)
        "c1ccccc1",    # Benzene (different)
        "CCCO",        # Propanol (similar alcohol)
        "CC(=O)O",     # Acetic acid (somewhat similar)
    ]

    # Generate embeddings
    query_result = await adapter.execute(query, model_type="gin")
    db_result = await adapter.execute(database, model_type="gin")

    if query_result.success and db_result.success:
        # Convert to numpy arrays
        query_emb = np.array(query_result.data['embeddings'])
        db_embs = np.array(db_result.data['embeddings'])

        # Compute cosine similarities
        from sklearn.metrics.pairwise import cosine_similarity
        similarities = cosine_similarity([query_emb], db_embs)[0]

        print(f"✓ Similarity to {query} (Ethanol):")
        for i, (smiles, sim) in enumerate(zip(database, similarities)):
            print(f"  {smiles:15s}: {sim:.4f}")


async def main():
    """Run all examples"""
    print("=" * 60)
    print("TorchDrug Adapter - Usage Examples")
    print("=" * 60)

    try:
        await example_basic_usage()
        await example_batch_processing()
        await example_different_models()
        await example_graph_features_only()
        await example_custom_configuration()
        await example_with_caching()
        await example_error_handling()
        await example_molecular_similarity()

    except Exception as e:
        print(f"\n✗ Error running examples: {e}")
        print("\nNote: Make sure TorchDrug is installed:")
        print("  pip install torch torchdrug")

    print("\n" + "=" * 60)
    print("Examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
