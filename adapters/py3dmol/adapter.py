"""
py3Dmol Adapter - Interactive 3D molecular visualization
Generates HTML/JavaScript for embedding 3D molecular visualizations in the frontend
Supports PDB, mol2, SMILES inputs with various visual styles
"""
from typing import Any, Dict, Optional, Union, Tuple
import logging
import os

try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError:
    PY3DMOL_AVAILABLE = False
    logging.warning("py3Dmol not available - install with: pip install py3Dmol")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - some features require it: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class Py3DmolAdapter(AdapterProtocol):
    """
    Adapter for 3D molecular visualization using py3Dmol
    Generates interactive HTML/JavaScript visualizations for molecular structures
    Supports PDB files, PDB strings, mol2 files, and SMILES strings
    Returns embeddable HTML/JS for frontend integration
    """

    def __init__(self):
        super().__init__(
            name="py3dmol",
            adapter_type="local",
            config={
                "timeout": 30,
                "default_style": "stick",
                "default_width": 400,
                "default_height": 400,
                "supported_styles": ["stick", "sphere", "cartoon", "surface", "line", "cross"],
                "supported_formats": ["pdb", "mol2", "sdf", "smiles"]
            }
        )
        self.version = "1.0.0"

        if not PY3DMOL_AVAILABLE:
            logger.error("py3Dmol is not installed! Install with: pip install py3Dmol")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a dictionary with required fields

        Args:
            input_data: Expected to be a dict with 'structure' key

        Returns:
            True if valid, False otherwise
        """
        if not PY3DMOL_AVAILABLE:
            return False

        # Check if input is a dictionary
        if not isinstance(input_data, dict):
            logger.warning("Py3Dmol: Input must be a dictionary")
            return False

        # Check if structure key exists
        if "structure" not in input_data:
            logger.warning("Py3Dmol: Input must contain 'structure' key")
            return False

        structure = input_data["structure"]

        # Check if structure is a string
        if not isinstance(structure, str):
            logger.warning("Py3Dmol: 'structure' must be a string")
            return False

        # Check if structure is not empty
        if len(structure.strip()) == 0:
            logger.warning("Py3Dmol: 'structure' cannot be empty")
            return False

        return True

    def _detect_structure_format(self, structure: str) -> str:
        """
        Detect the format of the input structure

        Args:
            structure: Structure string or file path

        Returns:
            Format string: 'pdb', 'mol2', 'sdf', or 'smiles'
        """
        # Check if it's a file path
        if os.path.isfile(structure):
            ext = os.path.splitext(structure)[1].lower()
            if ext in ['.pdb']:
                return 'pdb'
            elif ext in ['.mol2']:
                return 'mol2'
            elif ext in ['.sdf', '.mol']:
                return 'sdf'

        # Check if it looks like PDB format (has ATOM or HETATM lines)
        if 'ATOM' in structure or 'HETATM' in structure:
            return 'pdb'

        # Check if it looks like mol2 format
        if '@<TRIPOS>MOLECULE' in structure:
            return 'mol2'

        # Check if it looks like SDF format
        if structure.count('\n') > 3 and 'M  END' in structure:
            return 'sdf'

        # Otherwise assume it's SMILES
        return 'smiles'

    def _read_structure_file(self, file_path: str) -> Optional[str]:
        """
        Read structure from file

        Args:
            file_path: Path to structure file

        Returns:
            File contents as string or None on error
        """
        try:
            with open(file_path, 'r') as f:
                return f.read()
        except Exception as e:
            logger.error(f"Py3Dmol: Error reading file {file_path}: {e}")
            return None

    def _smiles_to_pdb(self, smiles: str) -> Optional[str]:
        """
        Convert SMILES to PDB format using RDKit

        Args:
            smiles: SMILES string

        Returns:
            PDB string or None on error
        """
        if not RDKIT_AVAILABLE:
            logger.error("Py3Dmol: RDKit required for SMILES conversion")
            return None

        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Py3Dmol: Could not parse SMILES: {smiles}")
                return None

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                logger.warning("Py3Dmol: Could not generate 3D coordinates, trying with random coords")
                AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)

            # Optimize geometry
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)

            # Convert to PDB
            pdb_string = Chem.MolToPDBBlock(mol)

            return pdb_string

        except Exception as e:
            logger.error(f"Py3Dmol: Error converting SMILES to PDB: {e}")
            return None

    def _prepare_structure(self, structure: str, format_hint: Optional[str] = None) -> Tuple[Optional[str], str]:
        """
        Prepare structure data for visualization

        Args:
            structure: Structure string or file path
            format_hint: Optional format hint ('pdb', 'mol2', 'sdf', 'smiles')

        Returns:
            Tuple of (structure_data, format)
        """
        # Check if it's a file path
        if os.path.isfile(structure):
            structure_data = self._read_structure_file(structure)
            if structure_data is None:
                return None, 'unknown'
            detected_format = self._detect_structure_format(structure_data)
            return structure_data, detected_format

        # Use format hint or detect format
        detected_format = format_hint if format_hint else self._detect_structure_format(structure)

        # Handle SMILES separately (needs conversion)
        if detected_format == 'smiles':
            pdb_data = self._smiles_to_pdb(structure)
            return pdb_data, 'pdb'

        # For other formats, return as-is
        return structure, detected_format

    def _determine_structure_type(self, structure_data: str, format: str) -> str:
        """
        Determine if structure is protein, ligand, or complex

        Args:
            structure_data: Structure data string
            format: Structure format

        Returns:
            'protein', 'ligand', or 'complex'
        """
        if format != 'pdb':
            return 'ligand'

        # Count atoms
        lines = structure_data.split('\n')
        atom_lines = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]

        if len(atom_lines) == 0:
            return 'ligand'

        # Check for protein markers (amino acid residues)
        protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}

        has_protein = False
        has_ligand = False

        for line in atom_lines:
            if len(line) < 20:
                continue

            # Extract residue name (columns 17-20 in PDB format)
            try:
                residue = line[17:20].strip()
                if residue in protein_residues:
                    has_protein = True
                elif line.startswith('HETATM'):
                    has_ligand = True
            except:
                continue

        if has_protein and has_ligand:
            return 'complex'
        elif has_protein:
            return 'protein'
        else:
            return 'ligand'

    def _generate_visualization_html(
        self,
        structure_data: str,
        structure_format: str,
        structure_type: str,
        style: str = "stick",
        width: int = 400,
        height: int = 400,
        background_color: str = "white",
        show_labels: bool = False
    ) -> Dict[str, str]:
        """
        Generate HTML/JavaScript for 3D visualization

        Args:
            structure_data: Structure data string
            structure_format: Format of structure ('pdb', 'mol2', 'sdf')
            structure_type: Type of structure ('protein', 'ligand', 'complex')
            style: Visual style
            width: Viewer width in pixels
            height: Viewer height in pixels
            background_color: Background color
            show_labels: Whether to show atom labels

        Returns:
            Dictionary with 'html' and 'javascript' keys
        """
        # Generate unique viewer ID
        import uuid
        viewer_id = f"viewer_{uuid.uuid4().hex[:8]}"

        # Build style configuration based on structure type and requested style
        style_config = self._get_style_config(style, structure_type)

        # Escape structure data for JavaScript
        import json
        structure_data_escaped = json.dumps(structure_data)

        # Generate JavaScript code
        javascript = f"""
// py3Dmol visualization
var element_{viewer_id} = $('#{viewer_id}');
var config_{viewer_id} = {{ backgroundColor: '{background_color}' }};
var viewer_{viewer_id} = $3Dmol.createViewer(element_{viewer_id}, config_{viewer_id});

// Add structure
viewer_{viewer_id}.addModel({structure_data_escaped}, '{structure_format}');

// Apply style
{style_config}

// Zoom to fit
viewer_{viewer_id}.zoomTo();

// Render
viewer_{viewer_id}.render();

// Make viewer interactive
viewer_{viewer_id}.setClickable({{}}, true, function(atom, viewer) {{
    if (atom) {{
        console.log('Clicked atom:', atom);
    }}
}});
"""

        # Generate HTML
        html = f"""
<div id="{viewer_id}" style="width: {width}px; height: {height}px; position: relative;"></div>
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<script>
{javascript}
</script>
"""

        return {
            "html": html,
            "javascript": javascript,
            "viewer_id": viewer_id
        }

    def _get_style_config(self, style: str, structure_type: str) -> str:
        """
        Get style configuration JavaScript code

        Args:
            style: Visual style name
            structure_type: Type of structure

        Returns:
            JavaScript code for applying style
        """
        if structure_type == "protein":
            # For proteins, use cartoon for backbone and stick for ligands
            return f"""
viewer_{"{viewer_id}"}.setStyle({{cartoon: {{color: 'spectrum'}}}});
"""
        elif structure_type == "complex":
            # For complexes, cartoon for protein, stick for ligand
            return f"""
// Style protein as cartoon
viewer_{"{viewer_id}"}.setStyle({{atom: ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']}}, {{cartoon: {{color: 'spectrum'}}}});
// Style ligand as stick
viewer_{"{viewer_id}"}.setStyle({{hetflag: true}}, {{{self._get_style_spec(style)}}});
"""
        else:
            # For ligands, use requested style
            return f"""
viewer_{"{viewer_id}"}.setStyle({{{self._get_style_spec(style)}}});
"""

    def _get_style_spec(self, style: str) -> str:
        """
        Get py3Dmol style specification

        Args:
            style: Style name

        Returns:
            JavaScript style specification
        """
        style_specs = {
            "stick": "stick: {colorscheme: 'default', radius: 0.15}",
            "sphere": "sphere: {colorscheme: 'default', scale: 0.3}",
            "cartoon": "cartoon: {color: 'spectrum'}",
            "surface": "surface: {opacity: 0.8, color: 'white'}",
            "line": "line: {colorscheme: 'default'}",
            "cross": "cross: {colorscheme: 'default', radius: 0.1}"
        }

        return style_specs.get(style.lower(), style_specs["stick"])

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute 3D molecular visualization

        Args:
            input_data: Dictionary with 'structure' key (PDB string, file path, or SMILES)
            **kwargs: Additional parameters:
                - style: Visual style ('stick', 'sphere', 'cartoon', 'surface', 'line', 'cross')
                - width: Viewer width in pixels (default: 400)
                - height: Viewer height in pixels (default: 400)
                - background_color: Background color (default: 'white')
                - show_labels: Show atom labels (default: False)
                - format: Format hint ('pdb', 'mol2', 'sdf', 'smiles')

        Returns:
            AdapterResult containing HTML/JavaScript for embedding
        """
        # Check if py3Dmol is available
        if not PY3DMOL_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="py3Dmol is not installed. Install with: pip install py3Dmol"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input format. Expected: {'structure': 'PDB/SMILES/file_path'}"
            )

        structure = input_data["structure"]
        style = kwargs.get("style", self.config["default_style"]).lower()
        width = kwargs.get("width", self.config["default_width"])
        height = kwargs.get("height", self.config["default_height"])
        background_color = kwargs.get("background_color", "white")
        show_labels = kwargs.get("show_labels", False)
        format_hint = kwargs.get("format", None)

        # Validate style
        if style not in self.config["supported_styles"]:
            logger.warning(f"Py3Dmol: Style '{style}' not supported, using stick")
            style = "stick"

        # Prepare structure data
        structure_data, structure_format = self._prepare_structure(structure, format_hint)

        if structure_data is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to prepare structure data. Check if file exists or SMILES is valid.",
                metadata={
                    "source": "py3dmol",
                    "structure_format": structure_format
                }
            )

        # Determine structure type
        structure_type = self._determine_structure_type(structure_data, structure_format)

        # Generate visualization
        try:
            viz_data = self._generate_visualization_html(
                structure_data=structure_data,
                structure_format=structure_format,
                structure_type=structure_type,
                style=style,
                width=width,
                height=height,
                background_color=background_color,
                show_labels=show_labels
            )

            result_data = {
                "html": viz_data["html"],
                "javascript": viz_data["javascript"],
                "viewer_id": viz_data["viewer_id"],
                "structure_type": structure_type,
                "structure_format": structure_format,
                "style": style,
                "width": width,
                "height": height
            }

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "py3dmol",
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "structure_type": structure_type,
                    "structure_format": structure_format,
                    "style": style
                }
            )

        except Exception as e:
            logger.error(f"Py3Dmol: Error generating visualization: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to generate visualization: {str(e)}",
                metadata={
                    "source": "py3dmol",
                    "structure_type": structure_type,
                    "structure_format": structure_format
                }
            )
