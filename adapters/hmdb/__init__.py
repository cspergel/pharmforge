"""
HMDB Adapter - Human Metabolome Database integration for PharmForge

Provides access to comprehensive metabolomics data including:
- Metabolite structures and properties
- Biofluid concentrations (blood, urine, CSF)
- Disease associations and biomarkers
- Metabolic pathway information
- Protein/enzyme interactions
"""

from .adapter import HMDBAdapter

__all__ = ["HMDBAdapter"]
