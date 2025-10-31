"""
Meeko Ligand Preparation Adapter

Prepares small molecules for docking with AutoDock Vina/GPU using Meeko.
Handles PDBQT generation, macrocycles, flexible residues, and conformer generation.
"""

from adapters.meeko.adapter import MeekoAdapter

__all__ = ['MeekoAdapter']
