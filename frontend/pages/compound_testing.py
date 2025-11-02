"""
Compound Testing Page for PharmForge - Enhanced Modern UI
Simple interface: SMILES input -> adapter selection -> results display
"""
import streamlit as st
from typing import Dict, List, Optional, Any
import time
import json
from datetime import datetime
import pandas as pd
from rdkit import Chem
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import PharmForge components
import sys
from pathlib import Path
parent_dir = Path(__file__).parent.parent
sys.path.insert(0, str(parent_dir))

from components.api_client import get_api_client
from components.molecule_viewer import molecule_viewer_component

# ============================================================================
# ADAPTER CATEGORIES
# ============================================================================

DOCKING_ADAPTERS = {
    'vina', 'gnina', 'diffdock', 'smina', 'vina_docking', 'gnina_docking'
}

MD_ADAPTERS = {
    'openmm'
}

STRUCTURE_ADAPTERS = {
    'alphafold', 'swissmodel', 'rcsb_pdb', 'pdb_redo'
}

# ============================================================================
# FAVORITES & HISTORY MANAGEMENT
# ============================================================================

def get_favorite_adapters():
    """Get favorite adapters from session state (localStorage simulation)"""
    if 'favorite_adapters' not in st.session_state:
        st.session_state.favorite_adapters = set()
    return st.session_state.favorite_adapters

def toggle_favorite(adapter_name):
    """Toggle favorite status for an adapter"""
    favorites = get_favorite_adapters()
    if adapter_name in favorites:
        favorites.remove(adapter_name)
    else:
        favorites.add(adapter_name)
    st.session_state.favorite_adapters = favorites

def save_to_history(results, smiles, protein_data=None):
    """Save test results to history"""
    if 'test_history' not in st.session_state:
        st.session_state.test_history = []

    history_entry = {
        'timestamp': datetime.now().isoformat(),
        'smiles': smiles,
        'protein_data': protein_data,
        'results': results,
        'summary': {
            'total': len(results),
            'successful': sum(1 for r in results if r['status'] == 'success'),
            'failed': sum(1 for r in results if r['status'] == 'error')
        }
    }

    st.session_state.test_history.insert(0, history_entry)
    st.session_state.test_history = st.session_state.test_history[:10]  # Keep last 10

def get_test_history():
    """Get test history"""
    if 'test_history' not in st.session_state:
        st.session_state.test_history = []
    return st.session_state.test_history

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string using RDKit"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception:
        return False

def group_adapters_by_category(adapters: List[Dict]) -> Dict[str, List[Dict]]:
    """Group adapters by category"""
    categories = {}
    for adapter in adapters:
        category = adapter.get('category', 'Other')
        if category not in categories:
            categories[category] = []
        categories[category].append(adapter)
    return categories

def get_status_icon(status: str) -> str:
    """Get emoji icon for adapter status"""
    status_lower = status.lower()
    if status_lower == 'healthy' or status_lower == 'online':
        return '✅'
    elif status_lower == 'degraded' or status_lower == 'warning':
        return '⚠️'
    elif status_lower == 'offline' or status_lower == 'error':
        return '❌'
    else:
        return '❓'

def format_execution_time(seconds: float) -> str:
    """Format execution time in human-readable format"""
    if seconds < 1:
        return f"{seconds*1000:.0f}ms"
    elif seconds < 60:
        return f"{seconds:.1f}s"
    else:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m {secs:.0f}s"

def requires_protein_input(selected_adapters: set) -> bool:
    """Check if any selected adapter requires protein input"""
    return bool(selected_adapters.intersection(DOCKING_ADAPTERS.union(MD_ADAPTERS)))

def export_results_json(results: List[Dict]) -> str:
    """Export results as JSON string"""
    export_data = []
    for r in results:
        export_data.append({
            'adapter': r['adapter'],
            'status': r['status'],
            'execution_time_seconds': r['time'],
            'cache_hit': r.get('cache_hit', False),
            'result': r['result']
        })
    return json.dumps(export_data, indent=2)

def export_results_csv(results: List[Dict]) -> str:
    """Export results as CSV string"""
    import io
    output = io.StringIO()

    # Create DataFrame
    rows = []
    for r in results:
        row = {
            'Adapter': r['adapter'],
            'Status': r['status'],
            'Time (s)': r['time'],
            'Cache Hit': r.get('cache_hit', False)
        }

        # Add key result fields if success
        if r['status'] == 'success' and isinstance(r['result'], dict):
            for key, value in r['result'].items():
                if isinstance(value, (int, float, str, bool)):
                    row[key] = value

        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output, index=False)
    return output.getvalue()

# ============================================================================
# UI COMPONENTS
# ============================================================================

def render_favorites_section(available_adapters):
    """Render favorites section at the top"""
    favorites = get_favorite_adapters()

    if favorites:
        st.markdown("### ⭐ Favorite Adapters")
        st.caption("Quick access to your favorite adapters - check to select")

        cols = st.columns(min(len(favorites), 4))
        for idx, adapter_name in enumerate(list(favorites)[:4]):
            # Find adapter info
            adapter_info = next((a for a in available_adapters if a['name'] == adapter_name), None)
            if adapter_info:
                with cols[idx % 4]:
                    status_icon = get_status_icon(adapter_info.get('status', 'unknown'))
                    is_selected = adapter_name in st.session_state.selected_adapters

                    # Show checkbox for selection
                    is_checked = st.checkbox(
                        f"{status_icon} {adapter_info['display_name']}",
                        value=is_selected,
                        key=f"fav_{adapter_name}"
                    )

                    # Update selection state
                    if is_checked and adapter_name not in st.session_state.selected_adapters:
                        st.session_state.selected_adapters.add(adapter_name)
                    elif not is_checked and adapter_name in st.session_state.selected_adapters:
                        st.session_state.selected_adapters.remove(adapter_name)

        st.divider()

def render_results_cards(results: List[Dict]):
    """Render results as modern cards instead of table"""

    # Summary metrics with better styling
    col1, col2, col3, col4 = st.columns(4)

    successful = sum(1 for r in results if r['status'] == 'success')
    failed = sum(1 for r in results if r['status'] == 'error')
    total_time = sum(r['time'] for r in results)
    cached = sum(1 for r in results if r.get('cache_hit', False))

    with col1:
        st.metric("✅ Successful", f"{successful}/{len(results)}")
    with col2:
        st.metric("❌ Failed", failed)
    with col3:
        st.metric("⏱️ Total Time", format_execution_time(total_time))
    with col4:
        st.metric("💾 Cached", cached)

    st.divider()

    # Results cards in 2-column layout
    for i in range(0, len(results), 2):
        cols = st.columns(2)

        for j, col in enumerate(cols):
            if i + j < len(results):
                result = results[i + j]

                with col:
                    # Card container with status-based styling
                    card_style = "success" if result['status'] == 'success' else "error"

                    with st.container():
                        # Header
                        header_col1, header_col2, header_col3 = st.columns([3, 1, 1])

                        with header_col1:
                            st.markdown(f"**{result['icon']} {result['adapter']}**")

                        with header_col2:
                            cache_badge = "💾" if result.get('cache_hit', False) else ""
                            st.caption(f"{cache_badge} {format_execution_time(result['time'])}")

                        with header_col3:
                            status_badge = "✅" if result['status'] == 'success' else "❌"
                            st.caption(status_badge)

                        # Expandable result details
                        with st.expander("View Details", expanded=False):
                            if result['status'] == 'success':
                                result_data = result['result']

                                # Try to extract key metrics for better display
                                if isinstance(result_data, dict):
                                    # Look for common result fields
                                    metrics_to_show = {}

                                    for key in ['binding_affinity', 'score', 'confidence', 'similarity', 'qed', 'sa_score', 'logp']:
                                        if key in result_data:
                                            metrics_to_show[key.replace('_', ' ').title()] = result_data[key]

                                    if metrics_to_show:
                                        # Show metrics in mini-columns
                                        metric_cols = st.columns(len(metrics_to_show))
                                        for idx, (metric_name, value) in enumerate(metrics_to_show.items()):
                                            with metric_cols[idx]:
                                                if isinstance(value, (int, float)):
                                                    st.metric(metric_name, f"{value:.2f}" if isinstance(value, float) else value)
                                                else:
                                                    st.metric(metric_name, value)

                                        st.divider()

                                    # Show full JSON
                                    st.json(result_data)
                                else:
                                    st.write(result_data)
                            else:
                                st.error(result['result'].get('error', 'Unknown error'))

                        st.markdown("---")

def load_hiv_protease_example():
    """Load the HIV Protease docking example"""
    st.session_state.smiles = "CC(C)CC1NC(=O)C(Cc2ccccc2)NC(=O)C(CC(N)=O)NC(=O)C(CC(C)C)NC(=O)C(Cc2c[nH]c3ccccc23)NC(=O)C(Cc2c[nH]c3ccccc23)NC(=O)C(CC(=O)O)NC(=O)C(CC(C)C)NC1=O"  # Saquinavir
    st.session_state.smiles_valid = True
    st.session_state.protein_input_mode = "PDB ID"
    st.session_state.protein_data = {
        "type": "pdb_id",
        "value": "1HSG"
    }
    # Auto-select docking adapters
    st.session_state.selected_adapters = {'vina', 'gnina'}

# ============================================================================
# MAIN PAGE
# ============================================================================

def render_compound_testing_page():
    """Main compound testing page - Modern, sleek design"""

    # Page config
    st.title("🧪 Compound Testing")
    st.caption("Test compounds across molecular property and docking adapters")

    # Initialize session state
    if 'smiles' not in st.session_state:
        st.session_state.smiles = ""
    if 'selected_adapters' not in st.session_state:
        st.session_state.selected_adapters = set()
    if 'results' not in st.session_state:
        st.session_state.results = None
    if 'smiles_valid' not in st.session_state:
        st.session_state.smiles_valid = False
    if 'protein_data' not in st.session_state:
        st.session_state.protein_data = None
    if 'protein_input_mode' not in st.session_state:
        st.session_state.protein_input_mode = "PDB ID"

    # Get API client
    api_client = get_api_client()

    # ========================================================================
    # SIDEBAR - RECENT RESULTS
    # ========================================================================

    with st.sidebar:
        st.markdown("### 📊 Recent Tests")

        history = get_test_history()

        if history:
            for idx, entry in enumerate(history[:5]):  # Show last 5
                with st.container():
                    timestamp = datetime.fromisoformat(entry['timestamp'])
                    st.caption(timestamp.strftime("%b %d, %H:%M"))
                    st.markdown(f"**{entry['summary']['successful']}/{entry['summary']['total']} passed**")
                    st.caption(f"SMILES: {entry['smiles'][:20]}...")

                    if st.button(f"Load", key=f"history_{idx}", use_container_width=True):
                        st.session_state.smiles = entry['smiles']
                        st.session_state.smiles_valid = True
                        if entry.get('protein_data'):
                            st.session_state.protein_data = entry['protein_data']
                        # Streamlit will naturally update on next interaction

                    st.divider()
        else:
            st.info("No recent tests")

    # ========================================================================
    # SMILES INPUT
    # ========================================================================

    col1, col2 = st.columns([4, 1])

    with col1:
        smiles_input = st.text_input(
            "SMILES String",
            value=st.session_state.smiles,
            placeholder="CC(=O)Oc1ccccc1C(=O)O (Aspirin)",
            label_visibility="collapsed"
        )

        if smiles_input != st.session_state.smiles:
            st.session_state.smiles = smiles_input
            st.session_state.smiles_valid = validate_smiles(smiles_input) if smiles_input else False

    with col2:
        if st.session_state.smiles:
            if st.session_state.smiles_valid:
                st.success("Valid SMILES")
            else:
                st.error("Invalid")

    # Quick examples with HIV Protease
    example_col1, example_col2, example_col3, example_col4 = st.columns(4)

    with example_col1:
        if st.button("💊 Aspirin", key="example_aspirin", use_container_width=True):
            st.session_state.smiles = "CC(=O)Oc1ccccc1C(=O)O"
            st.session_state.smiles_valid = True

    with example_col2:
        if st.button("☕ Caffeine", key="example_caffeine", use_container_width=True):
            st.session_state.smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
            st.session_state.smiles_valid = True

    with example_col3:
        if st.button("🏥 Ibuprofen", key="example_ibuprofen", use_container_width=True):
            st.session_state.smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
            st.session_state.smiles_valid = True

    with example_col4:
        if st.button("🎯 HIV Protease Example", key="example_hiv", use_container_width=True, type="primary"):
            load_hiv_protease_example()

    # ========================================================================
    # PROTEIN INPUT
    # ========================================================================

    if requires_protein_input(st.session_state.selected_adapters):
        st.divider()
        st.markdown("### 🧬 Protein Target")
        st.caption("Required for molecular docking")

        protein_input_mode = st.radio(
            "How would you like to provide the protein?",
            ["PDB ID", "Upload PDB File", "UniProt ID"],
            horizontal=True,
            key="protein_input_mode_radio"
        )

        st.session_state.protein_input_mode = protein_input_mode

        if protein_input_mode == "PDB ID":
            pdb_id = st.text_input(
                "Enter PDB ID",
                value=st.session_state.protein_data.get('value', '') if st.session_state.protein_data and st.session_state.protein_data.get('type') == 'pdb_id' else '',
                placeholder="e.g., 1HSG (HIV protease)",
                help="4-character PDB identifier"
            )
            if pdb_id and len(pdb_id) == 4:
                st.session_state.protein_data = {
                    "type": "pdb_id",
                    "value": pdb_id.upper()
                }
                st.success(f"✓ Will fetch {pdb_id.upper()} from RCSB PDB")
            elif pdb_id:
                st.warning("PDB IDs are 4 characters (e.g., 1HSG)")

        elif protein_input_mode == "Upload PDB File":
            protein_file = st.file_uploader(
                "Upload PDB file",
                type=['pdb'],
                help="Upload a protein structure in PDB format"
            )
            if protein_file:
                protein_content = protein_file.read().decode('utf-8')
                st.session_state.protein_data = {
                    "type": "pdb_file",
                    "value": protein_content,
                    "filename": protein_file.name
                }
                st.success(f"✓ Loaded {protein_file.name}")

        elif protein_input_mode == "UniProt ID":
            uniprot_id = st.text_input(
                "Enter UniProt ID",
                placeholder="e.g., P03366 (HIV protease)",
                help="Will fetch from AlphaFold database"
            )
            if uniprot_id:
                st.session_state.protein_data = {
                    "type": "uniprot_id",
                    "value": uniprot_id.upper()
                }
                st.success(f"✓ Will fetch {uniprot_id.upper()} from AlphaFold")

    st.divider()

    # ========================================================================
    # ADAPTER SELECTION
    # ========================================================================

    st.markdown("### 🔧 Select Adapters")

    # Fetch adapters
    response = api_client.list_adapters()

    if response.success and response.data:
        available_adapters = response.data.get('adapters', [])

        # Show favorites first
        render_favorites_section(available_adapters)

        # Group by category
        categories = group_adapters_by_category(available_adapters)

        # Category tabs
        tabs = st.tabs(list(categories.keys()))

        for tab, (category, adapters) in zip(tabs, categories.items()):
            with tab:
                # Display adapters in grid (3 columns)
                for i in range(0, len(adapters), 3):
                    cols = st.columns(3)

                    for j in range(3):
                        if i + j < len(adapters):
                            adapter = adapters[i + j]
                            adapter_name = adapter['name']
                            is_selected = adapter_name in st.session_state.selected_adapters
                            is_favorite = adapter_name in get_favorite_adapters()

                            with cols[j]:
                                # Adapter card
                                with st.container():
                                    # Header with favorite star
                                    header_col1, header_col2 = st.columns([4, 1])

                                    with header_col1:
                                        status_icon = get_status_icon(adapter.get('status', 'unknown'))
                                        st.markdown(f"{status_icon} **{adapter['display_name']}**")

                                    with header_col2:
                                        if st.button(
                                            "⭐" if is_favorite else "☆",
                                            key=f"star_{adapter_name}",
                                            help="Add to favorites"
                                        ):
                                            toggle_favorite(adapter_name)

                                    # Description
                                    st.caption(adapter.get('description', 'No description')[:60] + "...")

                                    # Checkbox for selection
                                    is_checked = st.checkbox(
                                        "Select this adapter",
                                        value=is_selected,
                                        key=f"select_{adapter_name}",
                                        label_visibility="collapsed"
                                    )

                                    # Update selection state
                                    if is_checked and adapter_name not in st.session_state.selected_adapters:
                                        st.session_state.selected_adapters.add(adapter_name)
                                    elif not is_checked and adapter_name in st.session_state.selected_adapters:
                                        st.session_state.selected_adapters.remove(adapter_name)
    else:
        st.error("Failed to load adapters")

    st.divider()

    # ========================================================================
    # RUN BUTTON
    # ========================================================================

    if st.session_state.selected_adapters:
        st.markdown(f"**Selected: {len(st.session_state.selected_adapters)} adapters**")

        # Check if protein data is needed but missing
        needs_protein = requires_protein_input(st.session_state.selected_adapters)
        has_protein = st.session_state.protein_data is not None

        can_run = st.session_state.smiles_valid and (not needs_protein or has_protein)

        if not can_run:
            if not st.session_state.smiles_valid:
                st.warning("⚠️ Please enter a valid SMILES string")
            elif needs_protein and not has_protein:
                st.warning("⚠️ Selected adapters require protein target data")

        if st.button(
            f"🚀 Run {len(st.session_state.selected_adapters)} Adapters",
            type="primary",
            use_container_width=True,
            disabled=not can_run
        ):
            # Create progress display
            progress_bar = st.progress(0)
            status_text = st.empty()

            results = []
            total_adapters = len(st.session_state.selected_adapters)

            status_text.text(f"Running {total_adapters} adapters...")

            # Use ThreadPoolExecutor for parallel execution
            with ThreadPoolExecutor(max_workers=5) as executor:
                futures = {}

                for adapter_name in st.session_state.selected_adapters:
                    future = executor.submit(
                        api_client.test_adapter,
                        adapter_name,
                        st.session_state.smiles,
                        st.session_state.protein_data
                    )
                    futures[future] = adapter_name

                # Collect results as they complete
                completed = 0
                for future in as_completed(futures):
                    adapter_name = futures[future]
                    try:
                        response = future.result()

                        if response.success:
                            results.append({
                                'adapter': adapter_name,
                                'status': 'success',
                                'result': response.data,
                                'time': response.data.get('execution_time', 0),
                                'cache_hit': response.data.get('cache_hit', False),
                                'icon': '✅'
                            })
                        else:
                            results.append({
                                'adapter': adapter_name,
                                'status': 'error',
                                'result': {'error': response.error},
                                'time': 0,
                                'icon': '❌'
                            })
                    except Exception as e:
                        results.append({
                            'adapter': adapter_name,
                            'status': 'error',
                            'result': {'error': str(e)},
                            'time': 0,
                            'icon': '❌'
                        })

                    completed += 1
                    progress_bar.progress(completed / total_adapters)
                    status_text.text(f"Completed {completed}/{total_adapters} adapters...")

            # Clear progress indicators
            progress_bar.empty()
            status_text.empty()

            # Save results
            st.session_state.results = results
            save_to_history(results, st.session_state.smiles, st.session_state.protein_data)

            st.success(f"✅ Completed! {len(results)} adapters executed.")
            # Results will display below automatically

    # ========================================================================
    # RESULTS DISPLAY
    # ========================================================================

    if st.session_state.results:
        st.divider()
        st.markdown("### 📊 Results")

        render_results_cards(st.session_state.results)

        # Export options
        st.divider()
        st.markdown("### 💾 Export Results")

        col1, col2, col3 = st.columns(3)

        with col1:
            csv_data = export_results_csv(st.session_state.results)
            st.download_button(
                label="📄 Download CSV",
                data=csv_data,
                file_name=f"pharmforge_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv",
                use_container_width=True
            )

        with col2:
            json_data = export_results_json(st.session_state.results)
            st.download_button(
                label="📋 Download JSON",
                data=json_data,
                file_name=f"pharmforge_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                mime="application/json",
                use_container_width=True
            )

        with col3:
            if st.button("🔄 Clear Results", use_container_width=True):
                st.session_state.results = None
                # Streamlit will naturally update on next interaction

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    render_compound_testing_page()
