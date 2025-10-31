"""
PyMOL Adapter - Publication-quality molecular visualization
Generates high-quality images and ray-traced renders of molecular structures
Supports SMILES, PDB, and MOL formats with multiple visualization styles
"""
from typing import Any, Dict, Optional, List, Tuple
import logging
import os
import tempfile
import uuid
import asyncio
from pathlib import Path

try:
    import pymol
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    logging.warning("PyMOL not available - install with: pip install pymol-open-source")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - some features require it: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class PyMOLAdapter(AdapterProtocol):
    """
    Adapter for publication-quality molecular visualization using PyMOL
    Generates high-quality images with customizable styles and rendering options
    Supports SMILES (via RDKit conversion), PDB files, and MOL files
    Designed for headless/server deployment with no GUI
    """

    # Supported visualization styles
    SUPPORTED_STYLES = ["sticks", "spheres", "cartoon", "surface", "lines", "ribbon", "mesh", "dots"]

    # Supported color schemes
    SUPPORTED_COLOR_SCHEMES = ["by_element", "by_residue", "rainbow", "spectrum", "white", "gray", "cyan"]

    # Supported output formats
    SUPPORTED_FORMATS = ["png", "jpg", "jpeg"]

    def __init__(self):
        super().__init__(
            name="pymol",
            adapter_type="local",
            config={
                "timeout": 60,
                "default_style": "sticks",
                "default_color_scheme": "by_element",
                "default_background": "white",
                "default_width": 800,
                "default_height": 600,
                "default_ray_trace": True,
                "temp_dir": None,  # Will use system temp if None
                "cleanup_temp_files": True
            }
        )
        self.version = "1.0.0"

        if not PYMOL_AVAILABLE:
            logger.error("PyMOL is not installed! Install with: pip install pymol-open-source")
        else:
            # Initialize PyMOL in headless mode
            try:
                pymol.finish_launching(['pymol', '-c'])  # -c flag for command-line mode
                logger.info("PyMOL initialized in headless mode")
            except Exception as e:
                logger.warning(f"PyMOL headless initialization warning: {e}")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data for PyMOL visualization

        Args:
            input_data: Expected to be a dict with 'structure' key

        Returns:
            True if valid, False otherwise
        """
        if not PYMOL_AVAILABLE:
            logger.warning("PyMOL: PyMOL is not installed")
            return False

        # Check if input is a dictionary
        if not isinstance(input_data, dict):
            logger.warning("PyMOL: Input must be a dictionary")
            return False

        # Check if structure key exists
        if "structure" not in input_data:
            logger.warning("PyMOL: Input must contain 'structure' key")
            return False

        structure = input_data["structure"]

        # Check if structure is a string
        if not isinstance(structure, str):
            logger.warning("PyMOL: 'structure' must be a string")
            return False

        # Check if structure is not empty
        if len(structure.strip()) == 0:
            logger.warning("PyMOL: 'structure' cannot be empty")
            return False

        return True

    def _detect_structure_format(self, structure: str) -> str:
        """
        Detect the format of the input structure

        Args:
            structure: Structure string or file path

        Returns:
            Format string: 'pdb', 'mol', 'sdf', or 'smiles'
        """
        # Check if it's a file path
        if os.path.isfile(structure):
            ext = os.path.splitext(structure)[1].lower()
            if ext in ['.pdb']:
                return 'pdb'
            elif ext in ['.mol', '.sdf']:
                return 'mol'
            elif ext in ['.mol2']:
                return 'mol2'

        # Check if it looks like PDB format (has ATOM or HETATM lines)
        if 'ATOM' in structure or 'HETATM' in structure:
            return 'pdb'

        # Check if it looks like SDF/MOL format
        if structure.count('\n') > 3 and 'M  END' in structure:
            return 'mol'

        # Otherwise assume it's SMILES
        return 'smiles'

    def _smiles_to_pdb(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to PDB format using RDKit

        Args:
            smiles: SMILES string

        Returns:
            PDB string or None on error
        """
        if not RDKIT_AVAILABLE:
            logger.error("PyMOL: RDKit required for SMILES conversion")
            return None

        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"PyMOL: Could not parse SMILES: {smiles}")
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                logger.warning("PyMOL: Could not generate 3D coordinates, trying with random coords")
                AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)

            # Optimize geometry
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

            # Convert to PDB
            pdb_string = Chem.MolToPDBBlock(mol)

            return pdb_string

        except Exception as e:
            logger.error(f"PyMOL: Error converting SMILES to PDB: {e}")
            return None

    def _prepare_structure(self, structure: str, format_hint: Optional[str] = None) -> Tuple[Optional[str], str, str]:
        """
        Prepare structure data for visualization

        Args:
            structure: Structure string or file path
            format_hint: Optional format hint ('pdb', 'mol', 'smiles')

        Returns:
            Tuple of (structure_data, format, file_path)
        """
        # Check if it's a file path
        if os.path.isfile(structure):
            detected_format = self._detect_structure_format(structure)
            return structure, detected_format, structure

        # Use format hint or detect format
        detected_format = format_hint if format_hint else self._detect_structure_format(structure)

        # Handle SMILES separately (needs conversion)
        if detected_format == 'smiles':
            pdb_data = self._smiles_to_pdb(structure)
            if pdb_data is None:
                return None, 'unknown', ''

            # Write to temp file
            temp_dir = self.config.get("temp_dir") or tempfile.gettempdir()
            temp_path = os.path.join(temp_dir, f"pymol_temp_{uuid.uuid4().hex}.pdb")
            with open(temp_path, 'w') as f:
                f.write(pdb_data)

            return pdb_data, 'pdb', temp_path

        # For other formats, write to temp file
        temp_dir = self.config.get("temp_dir") or tempfile.gettempdir()
        ext = 'pdb' if detected_format == 'pdb' else 'mol'
        temp_path = os.path.join(temp_dir, f"pymol_temp_{uuid.uuid4().hex}.{ext}")

        with open(temp_path, 'w') as f:
            f.write(structure)

        return structure, detected_format, temp_path

    def _apply_style(self, obj_name: str, style: str, color_scheme: str) -> None:
        """
        Apply visualization style to PyMOL object

        Args:
            obj_name: Name of PyMOL object
            style: Visualization style
            color_scheme: Color scheme to apply
        """
        # Hide everything first
        cmd.hide("everything", obj_name)

        # Apply the requested style
        if style == "sticks":
            cmd.show("sticks", obj_name)
        elif style == "spheres":
            cmd.show("spheres", obj_name)
        elif style == "cartoon":
            cmd.show("cartoon", obj_name)
        elif style == "surface":
            cmd.show("surface", obj_name)
        elif style == "lines":
            cmd.show("lines", obj_name)
        elif style == "ribbon":
            cmd.show("ribbon", obj_name)
        elif style == "mesh":
            cmd.show("mesh", obj_name)
        elif style == "dots":
            cmd.show("dots", obj_name)
        else:
            logger.warning(f"PyMOL: Unknown style '{style}', using sticks")
            cmd.show("sticks", obj_name)

        # Apply color scheme
        if color_scheme == "by_element":
            cmd.color("atomic", obj_name)
        elif color_scheme == "by_residue":
            cmd.spectrum("resi", "rainbow", obj_name)
        elif color_scheme == "rainbow":
            cmd.spectrum("count", "rainbow", obj_name)
        elif color_scheme == "spectrum":
            cmd.spectrum("b", "rainbow", obj_name)
        elif color_scheme == "white":
            cmd.color("white", obj_name)
        elif color_scheme == "gray":
            cmd.color("gray", obj_name)
        elif color_scheme == "cyan":
            cmd.color("cyan", obj_name)

    def _render_image(
        self,
        structure_file: str,
        output_path: str,
        style: str,
        color_scheme: str,
        background: str,
        width: int,
        height: int,
        ray_trace: bool
    ) -> Dict[str, Any]:
        """
        Render molecular image using PyMOL

        Args:
            structure_file: Path to structure file
            output_path: Path for output image
            style: Visualization style
            color_scheme: Color scheme
            background: Background color
            width: Image width
            height: Image height
            ray_trace: Whether to use ray tracing

        Returns:
            Dictionary with rendering metadata
        """
        import time
        start_time = time.time()

        # Clear PyMOL session
        cmd.reinitialize()

        # Set background color
        cmd.bg_color(background)

        # Load structure
        obj_name = "molecule"
        cmd.load(structure_file, obj_name)

        # Apply style and colors
        self._apply_style(obj_name, style, color_scheme)

        # Center and zoom
        cmd.center(obj_name)
        cmd.zoom(obj_name, buffer=2)

        # Set image dimensions
        cmd.viewport(width, height)

        # Render and save
        if ray_trace:
            # Ray-traced high-quality render
            cmd.ray(width, height)

        # Save image
        cmd.png(output_path, width=width, height=height, dpi=300, ray=1 if ray_trace else 0)

        render_time = time.time() - start_time

        # Get some metadata
        metadata = {
            "num_atoms": cmd.count_atoms(obj_name),
            "num_residues": cmd.count_atoms(f"{obj_name} and name ca"),
            "render_time": round(render_time, 2)
        }

        return metadata

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute PyMOL visualization and generate image

        Args:
            input_data: Dictionary with 'structure' key (SMILES, PDB string, or file path)
            **kwargs: Additional parameters:
                - format: Format hint ('smiles', 'pdb', 'mol')
                - style: Visualization style (default: 'sticks')
                - color_scheme: Color scheme (default: 'by_element')
                - background: Background color (default: 'white')
                - width: Image width in pixels (default: 800)
                - height: Image height in pixels (default: 600)
                - ray_trace: Use ray tracing for high quality (default: True)
                - output_path: Custom output path (optional)

        Returns:
            AdapterResult containing image path and metadata
        """
        # Check if PyMOL is available
        if not PYMOL_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="PyMOL is not installed. Install with: pip install pymol-open-source"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input format. Expected: {'structure': 'SMILES/PDB/file_path'}"
            )

        structure = input_data["structure"]
        format_hint = kwargs.get("format", None)
        style = kwargs.get("style", self.config["default_style"]).lower()
        color_scheme = kwargs.get("color_scheme", self.config["default_color_scheme"]).lower()
        background = kwargs.get("background", self.config["default_background"]).lower()
        width = kwargs.get("width", self.config["default_width"])
        height = kwargs.get("height", self.config["default_height"])
        ray_trace = kwargs.get("ray_trace", self.config["default_ray_trace"])
        output_path = kwargs.get("output_path", None)

        # Validate parameters
        if style not in self.SUPPORTED_STYLES:
            logger.warning(f"PyMOL: Style '{style}' not supported, using sticks")
            style = "sticks"

        if color_scheme not in self.SUPPORTED_COLOR_SCHEMES:
            logger.warning(f"PyMOL: Color scheme '{color_scheme}' not supported, using by_element")
            color_scheme = "by_element"

        # Prepare output path if not provided
        if output_path is None:
            temp_dir = self.config.get("temp_dir") or tempfile.gettempdir()
            output_path = os.path.join(temp_dir, f"pymol_render_{uuid.uuid4().hex}.png")

        # Prepare structure data
        structure_data, structure_format, temp_file_path = self._prepare_structure(structure, format_hint)

        if structure_data is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to prepare structure data. Check if file exists or SMILES is valid.",
                metadata={
                    "source": "pymol",
                    "structure_format": structure_format
                }
            )

        # Render in executor to avoid blocking
        temp_files_to_cleanup = [temp_file_path]

        try:
            # Run rendering in thread pool since PyMOL is synchronous
            loop = asyncio.get_event_loop()
            render_metadata = await loop.run_in_executor(
                None,
                self._render_image,
                temp_file_path,
                output_path,
                style,
                color_scheme,
                background,
                width,
                height,
                ray_trace
            )

            # Check if output file was created
            if not os.path.exists(output_path):
                raise Exception("Output file was not created")

            # Get file size
            file_size = os.path.getsize(output_path)

            result_data = {
                "image_path": output_path,
                "width": width,
                "height": height,
                "format": "png",
                "file_size": file_size,
                "render_time": render_metadata.get("render_time", 0),
                "num_atoms": render_metadata.get("num_atoms", 0),
                "style": style,
                "color_scheme": color_scheme,
                "background": background,
                "ray_traced": ray_trace,
                "warnings": []
            }

            # Cleanup temp files if configured
            if self.config.get("cleanup_temp_files", True):
                for temp_file in temp_files_to_cleanup:
                    try:
                        if os.path.exists(temp_file) and temp_file != output_path:
                            os.remove(temp_file)
                    except Exception as e:
                        logger.warning(f"PyMOL: Could not cleanup temp file {temp_file}: {e}")

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "pymol",
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "structure_format": structure_format,
                    "style": style,
                    "color_scheme": color_scheme,
                    "ray_traced": ray_trace
                }
            )

        except Exception as e:
            logger.error(f"PyMOL: Error rendering image: {e}")

            # Cleanup temp files on error
            if self.config.get("cleanup_temp_files", True):
                for temp_file in temp_files_to_cleanup:
                    try:
                        if os.path.exists(temp_file):
                            os.remove(temp_file)
                    except Exception as cleanup_error:
                        logger.warning(f"PyMOL: Could not cleanup temp file {temp_file}: {cleanup_error}")

            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to render image: {str(e)}",
                metadata={
                    "source": "pymol",
                    "structure_format": structure_format
                }
            )

    def render_batch(self, structures: List[Dict[str, Any]], **common_kwargs) -> List[AdapterResult]:
        """
        Render multiple structures in batch (convenience method)

        Args:
            structures: List of structure dictionaries
            **common_kwargs: Common parameters for all structures

        Returns:
            List of AdapterResults
        """
        results = []
        for structure_input in structures:
            result = asyncio.run(self.execute(structure_input, **common_kwargs))
            results.append(result)
        return results
