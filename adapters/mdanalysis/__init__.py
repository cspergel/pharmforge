"""
MDAnalysis Adapter for PharmForge

Provides molecular dynamics trajectory analysis using the MDAnalysis toolkit.

Supports analysis of:
- RMSD (Root Mean Square Deviation)
- RMSF (Root Mean Square Fluctuation)
- Radius of gyration
- Hydrogen bonds
- Multiple trajectory formats (DCD, XTC, TRR, NetCDF, etc.)
"""

from .adapter import MDAnalysisAdapter

__all__ = ["MDAnalysisAdapter"]
