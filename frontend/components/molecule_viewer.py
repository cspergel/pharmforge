"""
Molecule Viewer Component
Provides 2D molecule visualization using RDKit
"""
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64
from typing import Optional
import logging

logger = logging.getLogger(__name__)


def render_molecule_2d(
    smiles: str,
    width: int = 300,
    height: int = 300,
    show_smiles: bool = True
) -> Optional[str]:
    """
    Render 2D molecule structure from SMILES

    Args:
        smiles: SMILES string
        width: Image width in pixels
        height: Image height in pixels
        show_smiles: Whether to display SMILES below image

    Returns:
        Base64 encoded PNG image or None on error
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None

        # Generate 2D coordinates
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)

        # Draw molecule
        img = Draw.MolToImage(mol, size=(width, height))

        # Convert to base64
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()

        return img_str

    except Exception as e:
        logger.error(f"Error rendering molecule {smiles}: {e}")
        return None


def molecule_viewer_component(smiles: str, width: int = 300, height: int = 300):
    """
    Streamlit component for displaying molecule structure

    Args:
        smiles: SMILES string
        width: Image width
        height: Image height
    """
    img_str = render_molecule_2d(smiles, width, height)

    if img_str:
        # Display image
        st.markdown(
            f'<img src="data:image/png;base64,{img_str}" style="border-radius: 8px; border: 1px solid #e0e0e0;">',
            unsafe_allow_html=True
        )

        # Display SMILES below
        st.markdown(
            f'<div style="text-align: center; margin-top: 8px; font-family: monospace; font-size: 12px; color: #616161;">{smiles}</div>',
            unsafe_allow_html=True
        )
    else:
        # Error placeholder
        st.error(f"Unable to render molecule: {smiles}")
        st.markdown(
            f'<div style="text-align: center; padding: 40px; background: #f5f5f5; border-radius: 8px; border: 1px solid #e0e0e0;">'
            f'<p style="color: #9e9e9e;">Molecule rendering failed</p>'
            f'<code>{smiles}</code>'
            f'</div>',
            unsafe_allow_html=True
        )


def molecule_grid(smiles_list: list, columns: int = 4, width: int = 200, height: int = 200):
    """
    Display a grid of molecules

    Args:
        smiles_list: List of SMILES strings
        columns: Number of columns in grid
        width: Molecule image width
        height: Molecule image height
    """
    cols = st.columns(columns)

    for idx, smiles in enumerate(smiles_list):
        with cols[idx % columns]:
            molecule_viewer_component(smiles, width, height)
            st.caption(f"Molecule {idx + 1}")


def molecule_viewer_placeholder(smiles: str = "CC(=O)O"):
    """
    Placeholder for future 3D molecule viewer (Next.js version)

    Args:
        smiles: SMILES string to display
    """
    st.info("🧪 3D Molecule Viewer (Coming in Next.js Version)")

    st.markdown("""
    **Features in Next.js version:**
    - Interactive 3D rotation
    - Zoom and pan controls
    - Atom/bond highlighting
    - Property annotations
    - Export to PDB/SDF
    - Protein binding site overlay

    **For now, view 2D structure:**
    """)

    molecule_viewer_component(smiles)
