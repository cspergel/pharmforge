"""
RCSB PDB Adapter for PharmForge

Provides access to experimental protein structures from the
Research Collaboratory for Structural Bioinformatics Protein Data Bank.

Features:
- Download experimental structures by PDB ID
- Search by protein name or keywords
- Access ligand information and binding sites
- Structure quality metrics (resolution, R-factors)
"""

from .adapter import RCSBPDBAdapter

__all__ = ["RCSBPDBAdapter"]
