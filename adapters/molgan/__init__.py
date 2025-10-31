"""
MolGAN Adapter for PharmForge

Generative Adversarial Network for molecular graph generation.
Generates drug-like molecules with desired properties.

Reference: https://arxiv.org/abs/1805.11973
"""

from .adapter import MolGANAdapter

__all__ = ["MolGANAdapter"]
