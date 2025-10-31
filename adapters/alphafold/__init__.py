"""
AlphaFold DB Adapter for PharmForge

Provides access to AlphaFold protein structure predictions from the
AlphaFold Protein Structure Database (https://alphafold.ebi.ac.uk/).

Features:
- Download predicted structures by UniProt ID
- Access confidence scores (pLDDT)
- Quality metrics and model information
- Structure file caching
"""

from .adapter import AlphaFoldAdapter

__all__ = ["AlphaFoldAdapter"]
