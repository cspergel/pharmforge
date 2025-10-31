"""
COCONUT Natural Products Database Adapter

COlleCtion of Open Natural prodUcTs - the largest open natural products database

Features:
- 400,000+ natural product structures
- Biological source information (organism, taxonomy)
- Molecular properties and descriptors
- Natural product annotations
- CC0 public domain license

API: https://coconut.naturalproducts.net/api/
"""

from adapters.coconut.adapter import COCONUTAdapter

__all__ = ["COCONUTAdapter"]
__version__ = "1.0.0"
