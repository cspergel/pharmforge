"""
ProLIF Adapter - Protein-Ligand Interaction Fingerprints

This adapter provides protein-ligand interaction fingerprint analysis using ProLIF.
ProLIF analyzes protein-ligand complexes to identify and quantify various types
of non-covalent interactions.

Key Features:
-------------
- Comprehensive interaction detection (H-bonds, hydrophobic, pi-stacking, etc.)
- Support for single structures (docking poses) and trajectories (MD simulations)
- Residue-level interaction mapping
- Interaction frequency calculations
- Integration with MDAnalysis Universe objects

Supported Interaction Types:
----------------------------
- Hydrogen bonds (donor/acceptor)
- Hydrophobic contacts
- Pi-stacking interactions
- Pi-cation interactions
- Salt bridges (anionic/cationic)
- Halogen bonds
- Metal coordination
- Water bridges (if configured)

Input Formats:
-------------
The adapter accepts three input formats:

1. Separate protein and ligand files:
   {
       'protein_file': '/path/to/protein.pdb',
       'ligand_file': '/path/to/ligand.mol2'
   }

2. Complex file with ligand selection:
   {
       'complex_file': '/path/to/complex.pdb',
       'ligand_selection': 'resname LIG'
   }

3. MDAnalysis Universe object:
   {
       'universe': mda.Universe(...),
       'ligand_selection': 'resname LIG'
   }

Output Format:
-------------
The adapter returns a dictionary with:
- n_frames: Number of frames analyzed
- interaction_counts: Total count of each interaction type
- interaction_frequencies: Average occurrence per frame
- residue_pairs: List of interacting residue-interaction type pairs
- n_unique_interactions: Number of unique residue-interaction pairs
- total_interactions: Total interaction count across all frames
- avg_interactions_per_frame: Average interactions per frame
- most_common_interactions: Top 5 most frequent interaction types

Usage Example:
-------------
```python
from adapters.prolif.adapter import ProLIFAdapter

# Initialize adapter
adapter = ProLIFAdapter(config={
    'protein_selection': 'protein',
    'ligand_selection': 'resname LIG',
    'compute_frequency': True
})

# Analyze a docking pose
input_data = {
    'complex_file': '/path/to/docked_complex.pdb',
    'ligand_selection': 'resname LIG'
}

result = await adapter.execute(input_data)

if result.success:
    print(f"Found {result.data['total_interactions']} interactions")
    print(f"Interaction types: {result.data['interaction_counts']}")
    print(f"Unique residue pairs: {result.data['n_unique_interactions']}")
else:
    print(f"Error: {result.error}")
```

Requirements:
------------
- ProLIF: pip install prolif
- MDAnalysis: pip install MDAnalysis
- Pandas (optional): pip install pandas

Reference:
---------
ProLIF documentation: https://prolif.readthedocs.io/
"""

from .adapter import ProLIFAdapter

__all__ = ["ProLIFAdapter"]
__version__ = "1.0.0"
