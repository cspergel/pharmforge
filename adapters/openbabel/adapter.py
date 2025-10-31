"""
OpenBabel Adapter - Format conversion and ligand preparation
Provides molecular format conversion, 3D coordinate generation, and geometry optimization
"""
from typing import Any, Dict, Optional
import logging

try:
    from openbabel import openbabel as ob
    from openbabel import pybel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False
    logging.warning("OpenBabel not available - install with: pip install openbabel-wheel (or conda install openbabel)")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class OpenBabelAdapter(AdapterProtocol):
    """
    Adapter for OpenBabel format conversion and ligand preparation

    Capabilities:
    - Format conversion (SMILES, InChI, MOL, PDB, SDF, MOL2, etc.)
    - 3D coordinate generation
    - Add/remove hydrogens
    - Geometry optimization
    - Energy minimization
    """

    # Supported input formats
    SUPPORTED_INPUT_FORMATS = [
        'smi', 'smiles', 'inchi', 'mol', 'mol2', 'sdf', 'pdb',
        'xyz', 'cml', 'can', 'inchikey'
    ]

    # Supported output formats
    SUPPORTED_OUTPUT_FORMATS = [
        'smi', 'smiles', 'inchi', 'inchikey', 'mol', 'mol2', 'sdf',
        'pdb', 'xyz', 'cml', 'can'
    ]

    def __init__(self):
        super().__init__(
            name="openbabel",
            adapter_type="local",
            config={
                "timeout": 30,  # Local calculations
                "force_field": "mmff94",  # Default force field for optimization
                "gen_3d": True,  # Generate 3D coordinates by default
                "add_hydrogens": True,  # Add hydrogens by default
                "optimize_steps": 500  # Geometry optimization steps
            }
        )
        self.version = "1.0.0"

        if not OPENBABEL_AVAILABLE:
            logger.error("OpenBabel is not installed! Install with: pip install openbabel-wheel")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid molecular structure string or dict

        Args:
            input_data: SMILES string, InChI string, or dict with 'structure' and 'format'

        Returns:
            True if valid, False otherwise
        """
        if not OPENBABEL_AVAILABLE:
            return False

        # Handle dictionary input with format specification
        if isinstance(input_data, dict):
            structure = input_data.get('structure')
            input_format = input_data.get('format', 'smi').lower()

            if not structure or not isinstance(structure, str):
                return False
            if input_format not in self.SUPPORTED_INPUT_FORMATS:
                logger.warning(f"Unsupported input format: {input_format}")
                return False

            # Try to parse with OpenBabel
            try:
                mol = pybel.readstring(input_format, structure)
                return mol is not None
            except Exception as e:
                logger.warning(f"OpenBabel: Could not parse structure: {e}")
                return False

        # Handle string input (assume SMILES)
        elif isinstance(input_data, str):
            if len(input_data) == 0:
                return False

            # Try to parse as SMILES
            try:
                mol = pybel.readstring('smi', input_data)
                return mol is not None
            except Exception:
                # Try as InChI
                try:
                    mol = pybel.readstring('inchi', input_data)
                    return mol is not None
                except Exception:
                    return False

        return False

    def _convert_format(
        self,
        structure: str,
        input_format: str,
        output_format: str,
        gen_3d: bool = False,
        add_hydrogens: bool = False,
        optimize: bool = False,
        force_field: str = "mmff94"
    ) -> Optional[Dict[str, Any]]:
        """
        Convert molecular structure between formats

        Args:
            structure: Input molecular structure string
            input_format: Format of input (e.g., 'smi', 'inchi', 'mol')
            output_format: Desired output format
            gen_3d: Generate 3D coordinates
            add_hydrogens: Add hydrogen atoms
            optimize: Perform geometry optimization
            force_field: Force field for optimization (mmff94, uff, gaff, ghemical)

        Returns:
            Dictionary with converted structure and metadata
        """
        if not OPENBABEL_AVAILABLE:
            raise Exception("OpenBabel is not installed")

        try:
            # Parse input structure
            mol = pybel.readstring(input_format, structure)

            if mol is None:
                logger.warning(f"OpenBabel: Could not parse structure in format {input_format}")
                return None

            # Add hydrogens if requested
            if add_hydrogens:
                mol.addh()
                logger.debug("Added hydrogens to molecule")

            # Generate 3D coordinates if requested
            if gen_3d:
                mol.make3D(forcefield=force_field, steps=50)
                logger.debug("Generated 3D coordinates")

            # Optimize geometry if requested
            if optimize:
                ff = ob.OBForceField.FindForceField(force_field)
                if ff:
                    success = mol.localopt(forcefield=force_field, steps=self.config.get("optimize_steps", 500))
                    if success:
                        logger.debug(f"Optimized geometry with {force_field}")
                    else:
                        logger.warning(f"Geometry optimization with {force_field} failed")
                else:
                    logger.warning(f"Force field {force_field} not available")

            # Convert to output format
            output_str = mol.write(output_format).strip()

            # Calculate some basic properties
            properties = {
                "molecular_formula": mol.formula,
                "molecular_weight": mol.molwt,
                "num_atoms": len(mol.atoms),
                "num_heavy_atoms": len([a for a in mol.atoms if a.atomicnum > 1]),
                "num_bonds": len(mol.bonds),
                "num_rotors": mol.OBMol.NumRotors(),
                "exact_mass": mol.exactmass,
                "has_3d": gen_3d or optimize,
                "energy": mol.energy if (gen_3d or optimize) else None
            }

            # Generate multiple output formats if converting to/from SMILES
            result = {
                "structure": output_str,
                "format": output_format,
                "properties": properties
            }

            # Add canonical SMILES if not already the output format
            if output_format not in ['smi', 'smiles', 'can']:
                try:
                    result["canonical_smiles"] = mol.write('can').strip()
                except Exception:
                    pass

            # Add InChI if available and not the output format
            if output_format not in ['inchi', 'inchikey']:
                try:
                    result["inchi"] = mol.write('inchi').strip()
                    result["inchikey"] = mol.write('inchikey').strip()
                except Exception:
                    pass

            return result

        except Exception as e:
            logger.error(f"OpenBabel: Error converting structure: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute OpenBabel format conversion and/or ligand preparation

        Args:
            input_data: SMILES/InChI string or dict with 'structure' and 'format'
            **kwargs: Additional parameters:
                - output_format: Target format (default: 'mol2')
                - gen_3d: Generate 3D coordinates (default: True)
                - add_hydrogens: Add hydrogens (default: True)
                - optimize: Optimize geometry (default: False)
                - force_field: Force field for optimization (default: 'mmff94')

        Returns:
            AdapterResult containing converted structure and properties
        """
        # Check if OpenBabel is available
        if not OPENBABEL_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="OpenBabel is not installed. Install with: pip install openbabel-wheel"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid molecular structure or format"
            )

        # Parse input data
        if isinstance(input_data, dict):
            structure = input_data.get('structure')
            input_format = input_data.get('format', 'smi').lower()
        else:
            structure = input_data
            # Try to detect format (SMILES vs InChI)
            if structure.startswith('InChI='):
                input_format = 'inchi'
            else:
                input_format = 'smi'

        # Get conversion parameters
        output_format = kwargs.get('output_format', 'mol2').lower()
        gen_3d = kwargs.get('gen_3d', self.config.get('gen_3d', True))
        add_hydrogens = kwargs.get('add_hydrogens', self.config.get('add_hydrogens', True))
        optimize = kwargs.get('optimize', False)
        force_field = kwargs.get('force_field', self.config.get('force_field', 'mmff94'))

        # Validate output format
        if output_format not in self.SUPPORTED_OUTPUT_FORMATS:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported output format: {output_format}. Supported: {', '.join(self.SUPPORTED_OUTPUT_FORMATS)}"
            )

        # Perform conversion
        result = self._convert_format(
            structure=structure,
            input_format=input_format,
            output_format=output_format,
            gen_3d=gen_3d,
            add_hydrogens=add_hydrogens,
            optimize=optimize,
            force_field=force_field
        )

        if result is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to convert molecular structure with OpenBabel",
                metadata={
                    "source": "openbabel",
                    "input_format": input_format,
                    "output_format": output_format
                }
            )

        return AdapterResult(
            success=True,
            data=result,
            cache_hit=False,
            metadata={
                "source": "openbabel",
                "input_format": input_format,
                "output_format": output_format,
                "adapter_version": self.version,
                "computation_type": "local",
                "operations": {
                    "gen_3d": gen_3d,
                    "add_hydrogens": add_hydrogens,
                    "optimize": optimize,
                    "force_field": force_field if (gen_3d or optimize) else None
                }
            }
        )
