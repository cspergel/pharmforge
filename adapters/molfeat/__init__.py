"""
MolFeat Adapter - Unified molecular featurization library

MolFeat is a comprehensive library providing unified access to 100+ molecular
featurizers including 2D/3D fingerprints, physicochemical descriptors, and
pre-trained deep learning model embeddings.

Installation:
    pip install molfeat

Dependencies:
    - rdkit
    - molfeat

Features:
    - 100+ featurizers (Morgan, MACCS, RDKit, Avalon, etc.)
    - 2D and 3D descriptors
    - Pre-trained model embeddings (ChemBERTa, GIN, MolT5, ChemGPT)
    - Unified API for all featurizers
    - Batch processing support
    - Fast local computation
    - Access to MolFeat's model zoo

Supported Featurizer Types:
    2D Fingerprints:
        - morgan: Morgan (circular) fingerprints
        - ecfp: Extended Connectivity Fingerprints
        - fcfp: Functional Connectivity Fingerprints
        - maccs: MACCS keys (166-bit structural keys)
        - rdkit: RDKit fingerprints
        - atompair: Atom pair fingerprints
        - topological: Topological torsion fingerprints
        - avalon: Avalon fingerprints
        - erg: ErG fingerprints
        - estate: EState fingerprints

    Descriptors:
        - desc2d: 2D physicochemical descriptors
        - desc3d: 3D physicochemical descriptors (requires 3D coordinates)

    3D Shape-based:
        - usr: Ultrafast Shape Recognition
        - usrcat: USR with CREDO Atom Types
        - electroshape: Electroshape descriptors

    Pre-trained Models:
        - pretrained: Use with pretrained_model parameter
          Examples: "ChemBERTa-77M-MLM", "gin_supervised_contextpred",
                   "MolT5", "ChemGPT-1.2B"

Example usage:
    from adapters.molfeat import MolFeatAdapter

    # Initialize adapter
    adapter = MolFeatAdapter()

    # Example 1: Morgan fingerprints (default)
    result = await adapter.execute("CCO")  # Ethanol
    if result.success:
        features = result.data["features"]
        print(f"Generated {len(features)} features")

    # Example 2: MACCS keys
    result = await adapter.execute(
        "CCO",
        featurizer_type="maccs"
    )

    # Example 3: Custom Morgan parameters
    result = await adapter.execute(
        "CCO",
        featurizer_type="morgan",
        radius=3,
        nBits=4096
    )

    # Example 4: Batch processing
    smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
    result = await adapter.execute(
        smiles_list,
        featurizer_type="ecfp"
    )
    if result.success:
        features = result.data["features"]  # List of feature arrays
        print(f"Generated features for {len(features)} molecules")

    # Example 5: Pre-trained model embeddings
    result = await adapter.execute(
        "CCO",
        featurizer_type="pretrained",
        pretrained_model="ChemBERTa-77M-MLM"
    )
    if result.success:
        embeddings = result.data["features"]
        print(f"Generated embeddings: {embeddings.shape}")

    # Example 6: 2D descriptors
    result = await adapter.execute(
        "CCO",
        featurizer_type="desc2d"
    )
    if result.success:
        descriptors = result.data["features"]
        metadata = result.data
        print(f"Calculated {metadata['n_features']} descriptors")

Parameters:
    - featurizer_type: Type of featurizer to use (default: "morgan")
    - pretrained_model: Name of pre-trained model (when featurizer_type="pretrained")
    - radius: Fingerprint radius for circular fingerprints (default: 2)
    - nBits: Fingerprint size in bits (default: 2048)
    - return_raw: If True, return raw MolFeat features; if False, return
                  JSON-serializable format (default: False)

Returns:
    AdapterResult with:
        - features: Feature array or list of feature arrays
        - smiles: Input SMILES string(s)
        - featurizer_type: Type of featurizer used
        - n_molecules: Number of molecules processed
        - n_features: Number of features per molecule
        - is_binary: Whether features are binary (for fingerprints)
        - sparsity: Feature sparsity (for binary features)
        - pretrained_model: Name of pre-trained model (if used)
"""

from .adapter import MolFeatAdapter

__all__ = ["MolFeatAdapter"]
