"""
PharmForge Adapter Browser
Complete visualization and testing interface for all 75 PharmForge adapters
"""
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from typing import Dict, List, Optional, Any
import json

# Import API client
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from components.api_client import get_api_client

# ============================================================================
# ADAPTER DATA - Complete 75 Adapter Catalog
# ============================================================================

def get_all_adapters() -> List[Dict[str, Any]]:
    """
    Get comprehensive list of all PharmForge adapters from the API.
    """
    try:
        api_client = get_api_client()
        response = api_client._make_request("GET", "/api/adapters/list")

        if response.success and response.data:
            # Convert API format to expected format
            adapters = []
            for adapter in response.data.get('adapters', []):
                adapters.append({
                    'name': adapter.get('name', ''),
                    'display_name': adapter.get('display_name', ''),
                    'category': adapter.get('category', 'Other'),
                    'status': adapter.get('status', 'unknown'),
                    'version': adapter.get('version', '1.0.0'),
                    'description': adapter.get('description', ''),
                    'inputs': [],  # API doesn't provide this yet
                    'outputs': [],  # API doesn't provide this yet
                    'benchmarks': adapter.get('benchmarks', {}),
                    'required_params': [],
                    'optional_params': [],
                    'example_usage': f'adapter.execute({{"smiles": "CCO"}})',
                    'api_endpoint': 'API',
                    'adapter_type': adapter.get('adapter_type', 'api'),
                    'enabled': adapter.get('enabled', True)
                })
            return adapters
        else:
            st.error("Failed to load adapters from API")
            return []
    except Exception as e:
        st.error(f"Error fetching adapters: {str(e)}")
        return []


# ============================================================================
# FAVORITES MANAGEMENT
# ============================================================================

def get_favorite_adapters():
    """Get favorite adapters from session state"""
    if 'favorite_adapters' not in st.session_state:
        st.session_state.favorite_adapters = set()
    return st.session_state.favorite_adapters


def toggle_favorite(adapter_name: str):
    """Toggle favorite status for an adapter"""
    favorites = get_favorite_adapters()
    if adapter_name in favorites:
        favorites.remove(adapter_name)
    else:
        favorites.add(adapter_name)
    st.session_state.favorite_adapters = favorites


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def filter_adapters(
    adapters: List[Dict],
    search: str = "",
    category: str = "All",
    status: str = "All",
    sort_by: str = "name"
) -> List[Dict]:
    """Filter and sort adapters based on criteria."""
    filtered = adapters

    # Search filter
    if search:
        search_lower = search.lower()
        filtered = [
            a for a in filtered
            if search_lower in a['name'].lower()
            or search_lower in a['display_name'].lower()
            or search_lower in a['description'].lower()
        ]

    # Category filter
    if category != "All":
        filtered = [a for a in filtered if a['category'] == category]

    # Status filter
    if status != "All":
        filtered = [a for a in filtered if a['status'] == status.lower()]

    # Sort
    if sort_by == "name":
        filtered = sorted(filtered, key=lambda x: x['display_name'])
    elif sort_by == "category":
        filtered = sorted(filtered, key=lambda x: (x['category'], x['display_name']))
    elif sort_by == "status":
        status_order = {"healthy": 0, "degraded": 1, "offline": 2}
        filtered = sorted(filtered, key=lambda x: (status_order.get(x['status'], 3), x['display_name']))

    return filtered


def get_adapter_stats(adapters: List[Dict]) -> Dict[str, Any]:
    """Calculate adapter statistics."""
    total = len(adapters)

    # Count by status
    healthy = len([a for a in adapters if a['status'] == 'healthy'])
    degraded = len([a for a in adapters if a['status'] == 'degraded'])
    offline = len([a for a in adapters if a['status'] == 'offline'])

    # Count by category
    categories = {}
    for adapter in adapters:
        cat = adapter['category']
        categories[cat] = categories.get(cat, 0) + 1

    return {
        "total": total,
        "healthy": healthy,
        "degraded": degraded,
        "offline": offline,
        "categories": categories,
        "avg_response_time": "1.2s"  # Mock
    }


def show_test_modal_simple(adapter: Dict):
    """Display simplified test interface."""
    st.markdown(f"### Test: {adapter['display_name']}")

    with st.form(f"test_form_{adapter['name']}"):
        # Common SMILES input
        if "smiles" in adapter.get('inputs', []):
            smiles = st.text_input(
                "SMILES",
                value="CC(=O)Oc1ccccc1C(=O)O",
                help="Enter a SMILES string to test"
            )

        submitted = st.form_submit_button("Run Test", type="primary")

        if submitted:
            with st.spinner(f"Testing {adapter['display_name']}..."):
                # Mock test execution
                import time
                time.sleep(1)

                st.success("Test completed!")
                st.json({
                    "status": "success",
                    "execution_time": "1.23s",
                    "result": "Mock test result - API integration pending"
                })


def show_adapter_detail_page():
    """Display detailed adapter information on a separate page."""
    adapters = get_all_adapters()
    adapter_name = st.session_state.selected_adapter_detail

    # Find the adapter
    adapter = next((a for a in adapters if a['name'] == adapter_name), None)

    if not adapter:
        st.error("Adapter not found")
        if st.button("← Back to Browser"):
            del st.session_state.selected_adapter_detail
            # Streamlit will naturally return to main browser view
        return

    # Back button
    if st.button("← Back to Browser"):
        del st.session_state.selected_adapter_detail
        # Streamlit will naturally return to main browser view

    st.divider()

    # Header
    status_emoji = {
        "healthy": "✅",
        "degraded": "⚠️",
        "offline": "❌"
    }.get(adapter['status'], "❓")

    st.title(f"{adapter['display_name']}")
    st.caption(f"{status_emoji} {adapter['status'].title()} | {adapter['category']} | v{adapter['version']}")

    st.divider()

    # Description
    st.markdown("### Description")
    st.write(adapter['description'])

    st.divider()

    # Inputs & Outputs
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### Inputs")
        for inp in adapter.get('inputs', []):
            st.markdown(f"- `{inp}`")

    with col2:
        st.markdown("### Outputs")
        for out in adapter.get('outputs', []):
            st.markdown(f"- `{out}`")

    st.divider()

    # Benchmarks
    if 'benchmarks' in adapter:
        st.markdown("### Benchmarks")
        col1, col2, col3 = st.columns(3)

        benchmarks = adapter['benchmarks']
        metrics_list = list(benchmarks.items())

        for i, (metric, value) in enumerate(metrics_list):
            with [col1, col2, col3][i % 3]:
                st.metric(metric.replace('_', ' ').title(), value)

        st.divider()

    # Usage Example
    st.markdown("### Usage Example")
    st.code(adapter.get('example_usage', 'No example available'), language='python')

    st.divider()

    # Test Section
    st.markdown("### Test Adapter")
    show_test_modal_simple(adapter)


# ============================================================================
# MAIN PAGE
# ============================================================================

def adapter_browser_page():
    """Main adapter browser page - Clean, minimal design."""

    # Check if we're showing a detail page
    if 'selected_adapter_detail' in st.session_state:
        show_adapter_detail_page()
        return

    # Page header
    st.title("Adapter Browser")
    st.caption("Browse and test PharmForge adapters")

    # Get all adapters
    adapters = get_all_adapters()
    stats = get_adapter_stats(adapters)
    favorites = get_favorite_adapters()

    # ==================== FAVORITES SECTION ====================
    if favorites:
        st.markdown("### ⭐ Favorite Adapters")
        st.caption("Quick access to your favorite adapters")

        fav_cols = st.columns(min(len(favorites), 4))
        for idx, adapter_name in enumerate(list(favorites)[:4]):
            adapter_info = next((a for a in adapters if a['name'] == adapter_name), None)
            if adapter_info:
                with fav_cols[idx]:
                    status_emoji = {
                        "healthy": "✅",
                        "degraded": "⚠️",
                        "offline": "❌"
                    }.get(adapter_info.get('status', 'unknown'), "❓")

                    st.markdown(f"**{status_emoji} {adapter_info['display_name']}**")

                    # Quick action buttons
                    if st.button("🧪 Test", key=f"fav_test_{adapter_name}", use_container_width=True):
                        st.session_state.test_adapter = adapter_info
                    if st.button("📖 Details", key=f"fav_details_{adapter_name}", use_container_width=True):
                        st.session_state.selected_adapter_detail = adapter_name

        st.divider()

    # ==================== MODERN STATS ROW ====================
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("📦 Total Adapters", stats['total'])

    with col2:
        healthy_pct = f"{int(stats['healthy'] / stats['total'] * 100)}%" if stats['total'] > 0 else "0%"
        st.metric("✅ Healthy", f"{stats['healthy']} ({healthy_pct})")

    with col3:
        if stats['degraded'] > 0:
            st.metric("⚠️ Degraded", stats['degraded'])
        else:
            st.metric("⚠️ Degraded", "0")

    with col4:
        if stats['offline'] > 0:
            st.metric("❌ Offline", stats['offline'])
        else:
            st.metric("❌ Offline", "0")

    st.divider()

    # ==================== MODERN FILTERS ====================
    st.markdown("### 🔍 Filter & Search")

    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        search = st.text_input(
            "Search",
            placeholder="Search by name or description...",
            label_visibility="collapsed"
        )

    with col2:
        categories = ["All Categories"] + sorted(list(stats['categories'].keys()))
        category_filter = st.selectbox(
            "Category",
            categories,
            label_visibility="collapsed"
        )

    with col3:
        status_options = ["All Status", "Healthy", "Degraded", "Offline"]
        status_filter = st.selectbox(
            "Status",
            status_options,
            label_visibility="collapsed"
        )

    # Filter adapters
    category = category_filter if category_filter != "All Categories" else "All"
    status = status_filter if status_filter != "All Status" else "All"

    filtered_adapters = filter_adapters(
        adapters,
        search,
        category,
        status,
        "name"
    )

    # ==================== RESULTS HEADER ====================
    st.divider()

    col1, col2 = st.columns([3, 1])
    with col1:
        st.markdown(f"### 📦 Adapters")
        st.caption(f"Showing {len(filtered_adapters)} of {len(adapters)} adapters")
    with col2:
        # View toggle (could add compact/expanded view in future)
        pass

    st.divider()

    # ==================== CLEAN ADAPTER GRID ====================
    # Display 2 adapters per row for clean layout
    for i in range(0, len(filtered_adapters), 2):
        cols = st.columns(2)

        for col_idx, adapter_idx in enumerate([i, i + 1]):
            if adapter_idx >= len(filtered_adapters):
                break

            adapter = filtered_adapters[adapter_idx]
            adapter_name = adapter['name']
            favorites = get_favorite_adapters()
            is_favorite = adapter_name in favorites

            with cols[col_idx]:
                # Status badge
                status_emoji = {
                    "healthy": "✅",
                    "degraded": "⚠️",
                    "offline": "❌"
                }.get(adapter['status'], "❓")

                # Modern card container with subtle border
                with st.container():
                    # Header: Name + Favorite + Status
                    col_a, col_b, col_c = st.columns([3, 0.5, 1])
                    with col_a:
                        st.markdown(f"**{adapter['display_name']}**")
                    with col_b:
                        # Favorite star (no rerun needed)
                        star_button = st.button(
                            "⭐" if is_favorite else "☆",
                            key=f"fav_star_{adapter_name}",
                            help="Add to favorites",
                            use_container_width=True
                        )
                        if star_button:
                            toggle_favorite(adapter_name)
                            # Update will happen naturally
                    with col_c:
                        st.markdown(f"{status_emoji} {adapter['status'].title()}")

                    # One-line description
                    description = adapter.get('description', '')
                    if len(description) > 80:
                        description = description[:80] + "..."
                    st.caption(description)

                    # Category badge with icon
                    category_icons = {
                        "Docking": "🎯",
                        "Filtering": "🔍",
                        "Synthesis": "⚗️",
                        "Property": "📊",
                        "Generation": "✨",
                        "Optimization": "⚙️"
                    }
                    cat_icon = category_icons.get(adapter['category'], "📦")
                    st.caption(f"{cat_icon} {adapter['category']}")

                    # Action buttons
                    btn_col1, btn_col2 = st.columns(2)
                    with btn_col1:
                        if st.button("🧪 Test", key=f"test_{adapter['name']}", use_container_width=True):
                            st.session_state.test_adapter = adapter
                    with btn_col2:
                        if st.button("📖 Details", key=f"details_{adapter['name']}", use_container_width=True):
                            st.session_state.selected_adapter_detail = adapter['name']

                    st.divider()

    # ==================== TEST MODAL ====================
    if 'test_adapter' in st.session_state:
        show_test_modal_simple(st.session_state.test_adapter)
        if st.button("Close Test"):
            del st.session_state.test_adapter
            # Streamlit will naturally close modal on next interaction


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    adapter_browser_page()
