"""
xTB (extended Tight-Binding) Adapter for PharmForge

Provides fast semiempirical quantum mechanical calculations for molecular properties.
Uses the xtb-python interface to the xtb program.

Features:
- Geometry optimization
- Energy calculations
- Molecular orbital analysis (HOMO/LUMO)
- Molecular properties (dipole, polarizability, etc.)
- Multiple GFN methods (GFN0-xTB, GFN1-xTB, GFN2-xTB)

Reference: https://github.com/grimme-lab/xtb-python
"""

import hashlib
import logging
import tempfile
import os
from typing import Dict, Any, Optional, List
from pathlib import Path

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class xTBAdapter(AdapterProtocol):
    """
    xTB adapter for quantum chemistry calculations.

    Provides fast semiempirical QM calculations using GFN methods.
    Optimized for drug-like molecules (organic chemistry).

    Typical calculation time: 1-10 seconds per molecule (depending on size)

    Available GFN methods:
    - GFN2-xTB: Most accurate, recommended for most applications (default)
    - GFN1-xTB: Faster, good for large molecules
    - GFN0-xTB: Fastest, suitable for quick screening
    """

    # GFN method codes for xtb
    GFN_METHODS = {
        "GFN2-xTB": 2,
        "GFN1-xTB": 1,
        "GFN0-xTB": 0,
    }

    def __init__(self, name: str = "xtb", adapter_type: str = "local", config: Optional[Dict[str, Any]] = None):
        """
        Initialize xTB adapter.

        Args:
            name: Adapter name (default: "xtb")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - method: GFN method to use ("GFN2-xTB", "GFN1-xTB", "GFN0-xTB")
                   - optimize: Whether to optimize geometry (default: True)
                   - charge: Molecular charge (default: 0)
                   - uhf: Number of unpaired electrons (default: 0)
                   - solvent: Solvent for GBSA (default: None, options: "water", "methanol", etc.)
                   - accuracy: Numerical accuracy (default: 1.0, lower = more accurate)
        """
        super().__init__(name, adapter_type, config)
        self.version = "1.0.0"

        # Default configuration
        self.method = self.config.get("method", "GFN2-xTB")
        self.optimize = self.config.get("optimize", True)
        self.charge = self.config.get("charge", 0)
        self.uhf = self.config.get("uhf", 0)
        self.solvent = self.config.get("solvent", None)
        self.accuracy = self.config.get("accuracy", 1.0)

        # Validate method
        if self.method not in self.GFN_METHODS:
            logger.warning(f"Invalid GFN method '{self.method}', defaulting to GFN2-xTB")
            self.method = "GFN2-xTB"

        # Lazy-load xtb (heavy import)
        self._xtb_available = None

    def _check_xtb_available(self) -> bool:
        """Check if xtb-python is available."""
        if self._xtb_available is None:
            try:
                import xtb
                import xtb.interface
                self._xtb_available = True
                logger.info("✓ xtb-python is available")
            except ImportError:
                self._xtb_available = False
                logger.warning(
                    "xtb-python not installed. Install with: "
                    "conda install -c conda-forge xtb-python"
                )
        return self._xtb_available

    def validate_input(self, smiles: str) -> bool:
        """
        Validate SMILES string input.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not smiles or not isinstance(smiles, str):
            return False

        # Basic SMILES validation with RDKit
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.warning(f"SMILES validation error: {e}")
            return False

    def _smiles_to_coords(self, smiles: str) -> tuple:
        """
        Convert SMILES to 3D coordinates using RDKit.

        Args:
            smiles: SMILES string

        Returns:
            Tuple of (atomic_numbers, coordinates) where:
                - atomic_numbers: List of atomic numbers
                - coordinates: Nx3 numpy array of coordinates in Angstroms
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        import numpy as np

        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            # Fallback: try without random seed
            result = AllChem.EmbedMolecule(mol)
            if result != 0:
                raise ValueError(f"Failed to generate 3D coordinates for: {smiles}")

        # Optimize geometry with UFF
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)

        # Extract atomic numbers and coordinates
        conf = mol.GetConformer()
        atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        coordinates = np.array([
            [conf.GetAtomPosition(i).x,
             conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ])

        return atomic_numbers, coordinates

    async def execute(self, smiles: str, **params) -> AdapterResult:
        """
        Execute xTB calculation for a molecule.

        Args:
            smiles: SMILES string of the molecule
            **params: Additional parameters:
                - method: Override default GFN method
                - optimize: Override default optimization setting
                - charge: Override default molecular charge
                - uhf: Override default unpaired electrons
                - solvent: Override default solvent

        Returns:
            AdapterResult containing quantum chemical properties
        """
        try:
            # Check if xtb is available
            if not self._check_xtb_available():
                return AdapterResult(
                    success=False,
                    data={},
                    error="xtb-python not installed. Install with: conda install -c conda-forge xtb-python"
                )

            # Validate input
            if not self.validate_input(smiles):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid SMILES string"
                )

            # Get parameters (allow override via params)
            method = params.get("method", self.method)
            optimize = params.get("optimize", self.optimize)
            charge = params.get("charge", self.charge)
            uhf = params.get("uhf", self.uhf)
            solvent = params.get("solvent", self.solvent)

            # Generate cache key
            cache_key = self.generate_cache_key(
                smiles,
                method=method,
                optimize=optimize,
                charge=charge,
                uhf=uhf,
                solvent=solvent
            )

            logger.info(f"Running xTB calculation for: {smiles[:50]}...")
            logger.info(f"  Method: {method}, Optimize: {optimize}, Charge: {charge}, UHF: {uhf}")

            # Import xtb modules
            import xtb
            from xtb.interface import Calculator, Param
            import numpy as np

            # Convert SMILES to 3D coordinates
            logger.info("  Generating 3D structure...")
            atomic_numbers, coordinates = self._smiles_to_coords(smiles)

            # Create xtb calculator
            calc = Calculator(
                Param.GFN2xTB if method == "GFN2-xTB" else
                Param.GFN1xTB if method == "GFN1-xTB" else
                Param.GFN0xTB,
                np.array(atomic_numbers, dtype=np.int32),
                coordinates,
                charge=float(charge),
                uhf=int(uhf)
            )

            # Set accuracy
            calc.set_accuracy(self.accuracy)

            # Apply solvent if specified
            if solvent:
                logger.info(f"  Applying GBSA solvent model: {solvent}")
                calc.set_solvent(solvent)

            # Run calculation
            if optimize:
                logger.info("  Optimizing geometry...")
                calc.singlepoint()  # Initial single point
                res = calc.optimize()  # Geometry optimization
                logger.info(f"  Optimization converged: {res.get('converged', False)}")
            else:
                logger.info("  Running single point calculation...")
                res = calc.singlepoint()

            # Get optimized coordinates
            opt_coordinates = calc.get_positions()

            # Get results
            results = calc.get_results()

            # Extract key properties
            energy = results.get("energy", None)  # Total energy in Hartree
            gradient = results.get("gradient", None)  # Energy gradient

            # Get orbital energies
            orbital_energies = results.get("orbital energies", None)
            homo_lumo_gap = None
            homo_energy = None
            lumo_energy = None

            if orbital_energies is not None and len(orbital_energies) > 0:
                # Find HOMO and LUMO
                n_electrons = sum(atomic_numbers) - charge
                homo_idx = n_electrons // 2 - 1  # 0-indexed
                lumo_idx = homo_idx + 1

                if homo_idx >= 0 and homo_idx < len(orbital_energies):
                    homo_energy = float(orbital_energies[homo_idx])
                if lumo_idx < len(orbital_energies):
                    lumo_energy = float(orbital_energies[lumo_idx])
                if homo_energy is not None and lumo_energy is not None:
                    homo_lumo_gap = lumo_energy - homo_energy

            # Get dipole moment
            dipole = results.get("dipole", None)
            dipole_magnitude = None
            if dipole is not None:
                dipole_magnitude = float(np.linalg.norm(dipole))

            # Build result data
            result_data = {
                "smiles": smiles,
                "method": method,
                "optimized": optimize,
                "charge": charge,
                "uhf": uhf,
                "solvent": solvent,
                "energy_hartree": float(energy) if energy is not None else None,
                "energy_kcal_mol": float(energy * 627.509) if energy is not None else None,  # Convert to kcal/mol
                "homo_energy_eV": float(homo_energy * 27.2114) if homo_energy is not None else None,  # Convert to eV
                "lumo_energy_eV": float(lumo_energy * 27.2114) if lumo_energy is not None else None,
                "homo_lumo_gap_eV": float(homo_lumo_gap * 27.2114) if homo_lumo_gap is not None else None,
                "dipole_moment_debye": float(dipole_magnitude * 2.54177) if dipole_magnitude is not None else None,  # Convert to Debye
                "n_atoms": len(atomic_numbers),
                "converged": res.get("converged", True) if optimize else True,
            }

            # Add optimized structure
            if optimize:
                result_data["optimized_coordinates"] = opt_coordinates.tolist()

            logger.info(f"✓ xTB calculation complete")
            logger.info(f"  Energy: {result_data['energy_kcal_mol']:.2f} kcal/mol")
            if homo_lumo_gap is not None:
                logger.info(f"  HOMO-LUMO gap: {result_data['homo_lumo_gap_eV']:.2f} eV")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "cache_key": cache_key,
                    "version": self.version,
                    "xtb_version": method
                }
            )

        except ImportError as e:
            logger.error(f"Import error: {e}")
            return AdapterResult(
                success=False,
                data={},
                error=f"Missing dependency: {str(e)}. Install with: conda install -c conda-forge xtb-python rdkit"
            )
        except Exception as e:
            logger.error(f"xTB calculation failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name
                }
            )

    def generate_cache_key(
        self,
        smiles: str,
        method: str = "GFN2-xTB",
        optimize: bool = True,
        charge: int = 0,
        uhf: int = 0,
        solvent: Optional[str] = None
    ) -> str:
        """
        Generate deterministic cache key for calculations.

        Args:
            smiles: SMILES string
            method: GFN method
            optimize: Whether geometry was optimized
            charge: Molecular charge
            uhf: Unpaired electrons
            solvent: Solvent model

        Returns:
            SHA256 hash as cache key
        """
        # Canonicalize SMILES for consistent caching
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            pass  # Use original SMILES if canonicalization fails

        # Create cache key string
        key_string = (
            f"xtb_v{self.version}:{smiles}:{method}:"
            f"opt={optimize}:charge={charge}:uhf={uhf}:solvent={solvent}"
        )

        return hashlib.sha256(key_string.encode()).hexdigest()

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata.

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "Fast semiempirical quantum chemistry calculations using xTB",
            "methods": list(self.GFN_METHODS.keys()),
            "capabilities": [
                "Geometry optimization",
                "Energy calculations",
                "HOMO/LUMO analysis",
                "Molecular properties",
                "Solvent models (GBSA)"
            ],
            "typical_runtime": "1-10 seconds per molecule",
            "reference": "https://github.com/grimme-lab/xtb",
            "citation": "Grimme, S. et al. J. Chem. Theory Comput. 2017, 13, 1989-2009",
            "config": {
                "method": self.method,
                "optimize": self.optimize,
                "charge": self.charge,
                "uhf": self.uhf,
                "solvent": self.solvent,
                "accuracy": self.accuracy
            }
        }
