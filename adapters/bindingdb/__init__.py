"""
BindingDB Adapter - Experimental binding affinity database

Queries BindingDB for experimental binding data including Ki, Kd, IC50 values.
"""

from .adapter import BindingDBAdapter

__all__ = ["BindingDBAdapter"]
