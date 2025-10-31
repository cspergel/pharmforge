"""
PDB-REDO Adapter for PharmForge

Provides access to re-refined and optimized PDB structures from PDB-REDO.

Features:
- Download re-refined structures by PDB ID
- Improved structure quality metrics
- Comparison with original PDB structure
- Quality validation reports
"""

from .adapter import PDBRedoAdapter

__all__ = ["PDBRedoAdapter"]
