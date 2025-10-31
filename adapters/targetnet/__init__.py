"""
TargetNet Adapter for PharmForge

Deep learning-based target prediction for small molecules.
Predicts protein targets and off-target interactions.

Reference: Target prediction networks for drug discovery
"""

from .adapter import TargetNetAdapter

__all__ = ["TargetNetAdapter"]
