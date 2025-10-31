"""
SWISS-MODEL Adapter for PharmForge

Provides access to homology modeling and structure quality assessment
from the SWISS-MODEL Repository.

Features:
- Homology model retrieval by UniProt ID
- Structure quality assessment (QMEAN, GMQE)
- Template information
- Model confidence scores
"""

from .adapter import SwissModelAdapter

__all__ = ["SwissModelAdapter"]
