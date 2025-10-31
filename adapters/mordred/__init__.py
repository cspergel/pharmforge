"""
Mordred Adapter - Comprehensive molecular descriptor calculator

Calculates 1800+ molecular descriptors using the Mordred library.
Mordred is a molecular descriptor calculator that provides a comprehensive
set of 2D and 3D descriptors for molecular characterization.

Installation:
    pip install mordred

Dependencies:
    - rdkit
    - mordred

Features:
    - 1800+ molecular descriptors
    - 2D and 3D descriptors (3D optional)
    - Fast local computation
    - No API calls required
    - Comprehensive molecular characterization

Example usage:
    from adapters.mordred import MordredAdapter

    adapter = MordredAdapter()
    result = await adapter.execute("CCO")  # Ethanol

    if result.success:
        descriptors = result.data["descriptors"]
        metadata = result.data["metadata"]
        print(f"Calculated {metadata['valid_descriptors']} descriptors")
"""

from .adapter import MordredAdapter

__all__ = ["MordredAdapter"]
