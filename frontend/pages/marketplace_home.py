"""
PharmForge Marketplace - Home Page
Mock implementation showing UX design
"""
import streamlit as st
from typing import List, Dict

# This is a MOCKUP - shows design, not connected to real backend yet

# Note: Page config is set in main streamlit_app.py, not here


def mock_trending_items() -> List[Dict]:
    """Mock data for trending items"""
    return [
        {
            "name": "EGFR-Docking-v2",
            "author": "pharma-lab-stanford",
            "type": "model",
            "category": "Docking",
            "description": "AutoDock Vina model trained on 12,000 EGFR-ligand complexes",
            "rating": 4.8,
            "reviews": 234,
            "downloads": 12300,
            "metrics": {"AUC": 0.89, "Runtime": "3.2s"},
            "license": "MIT",
            "updated": "2 weeks ago"
        },
        {
            "name": "EGFR-Inhibitor-Dataset",
            "author": "medchem-ai",
            "type": "dataset",
            "category": "Binding Affinity",
            "description": "45K EGFR inhibitors with IC50 values from ChEMBL + patents",
            "rating": 4.9,
            "reviews": 156,
            "downloads": 8900,
            "metrics": {"Compounds": 45000, "Size": "120 MB"},
            "license": "CC-BY-4.0",
            "updated": "1 month ago"
        },
        {
            "name": "Tox21-Predictor",
            "author": "deepchem-org",
            "type": "model",
            "category": "ADMET",
            "description": "Multi-task toxicity prediction across 12 endpoints",
            "rating": 4.7,
            "reviews": 189,
            "downloads": 15600,
            "metrics": {"AUC": 0.83, "MAE": 0.12},
            "license": "Apache-2.0",
            "updated": "3 days ago"
        },
        {
            "name": "USPTO-Retrosynthesis",
            "author": "rxn-ai",
            "type": "adapter",
            "category": "Retrosynthesis",
            "description": "IBM RXN API adapter with caching and rate limiting",
            "rating": 4.6,
            "reviews": 92,
            "downloads": 4200,
            "metrics": {"Success Rate": "78%", "Avg Routes": 5},
            "license": "MIT",
            "updated": "1 week ago"
        }
    ]


def mock_featured_collections() -> List[Dict]:
    """Mock data for featured collections"""
    return [
        {
            "name": "Complete EGFR Pipeline",
            "author": "pharma-lab-stanford",
            "description": "Full workflow for EGFR inhibitor discovery",
            "items": 4,
            "stars": 89,
            "updated": "2 weeks ago"
        },
        {
            "name": "Kinase Inhibitor Toolkit",
            "author": "medchem-ai",
            "description": "Models and datasets for kinase drug discovery",
            "items": 7,
            "stars": 156,
            "updated": "1 month ago"
        },
        {
            "name": "ADMET Prediction Suite",
            "author": "deepchem-org",
            "description": "Comprehensive toxicity and property prediction",
            "items": 5,
            "stars": 203,
            "updated": "3 days ago"
        }
    ]


def render_hero_section():
    """Hero section with search"""
    st.markdown("""
    <style>
    .hero {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 60px 40px;
        border-radius: 16px;
        color: white;
        margin-bottom: 40px;
    }
    .hero h1 {
        font-size: 48px;
        font-weight: 700;
        margin-bottom: 20px;
    }
    .hero p {
        font-size: 20px;
        opacity: 0.95;
        margin-bottom: 30px;
    }
    </style>

    <div class="hero">
        <h1>🏪 PharmForge Marketplace</h1>
        <p>Discover, share, and collaborate on drug discovery models, datasets, and tools</p>
    </div>
    """, unsafe_allow_html=True)

    # Search bar
    col1, col2 = st.columns([3, 1])
    with col1:
        search_query = st.text_input(
            "Search models, datasets, adapters, pipelines...",
            placeholder="e.g., EGFR docking model",
            label_visibility="collapsed"
        )
    with col2:
        search_type = st.selectbox(
            "Type",
            ["All", "Models", "Datasets", "Adapters", "Pipelines"],
            label_visibility="collapsed"
        )

    if search_query:
        st.info(f"🔍 Searching for '{search_query}' in {search_type}...")
        st.markdown("*Search functionality coming soon - this is a mockup*")


def render_stats_bar():
    """Stats overview"""
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Models", "1,234", "+89 this week")
    with col2:
        st.metric("Datasets", "567", "+23 this week")
    with col3:
        st.metric("Adapters", "189", "+12 this week")
    with col4:
        st.metric("Contributors", "3,456", "+145 this week")


def render_trending_section():
    """Trending items"""
    st.markdown("### 🔥 Trending This Week")

    items = mock_trending_items()

    # Display as cards in grid
    cols = st.columns(2)
    for idx, item in enumerate(items):
        with cols[idx % 2]:
            render_item_card(item)


def render_item_card(item: Dict):
    """Render a marketplace item card"""
    icon_map = {
        "model": "🤖",
        "dataset": "📊",
        "adapter": "🔌",
        "pipeline": "🔄"
    }

    # Build metrics HTML separately to avoid f-string issues
    metrics_html = ""
    for k, v in item['metrics'].items():
        metrics_html += f'<span style="background: #f3f4f6; padding: 4px 8px; border-radius: 4px; font-size: 12px; margin-right: 8px;">{k}: {v}</span>'

    st.markdown(f"""
    <div style="
        border: 1px solid #e0e0e0;
        border-radius: 12px;
        padding: 20px;
        margin-bottom: 16px;
        background: white;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        transition: box-shadow 0.2s;
    ">
        <div style="display: flex; justify-content: space-between; align-items: start;">
            <div style="flex: 1;">
                <h4 style="margin: 0 0 8px 0; color: #1f2937;">
                    {icon_map[item['type']]} {item['name']}
                </h4>
                <p style="margin: 0 0 4px 0; color: #6b7280; font-size: 14px;">
                    by @{item['author']} • {item['category']}
                </p>
                <p style="margin: 0 0 12px 0; color: #374151; font-size: 14px;">
                    {item['description']}
                </p>
            </div>
        </div>

        <div style="display: flex; gap: 20px; margin-bottom: 12px;">
            <span style="color: #6b7280; font-size: 13px;">
                ⭐ {item['rating']}/5 ({item['reviews']} reviews)
            </span>
            <span style="color: #6b7280; font-size: 13px;">
                ⬇️ {item['downloads']:,} downloads
            </span>
        </div>

        <div style="display: flex; gap: 12px; margin-bottom: 12px;">
            {metrics_html}
        </div>

        <div style="display: flex; justify-content: space-between; align-items: center;">
            <span style="color: #9ca3af; font-size: 12px;">
                🔒 {item['license']} • Updated {item['updated']}
            </span>
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Action buttons
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("View Details", key=f"view_{item['name']}", use_container_width=True):
            st.info(f"Opening {item['name']} details...")
    with col2:
        if st.button("Download", key=f"download_{item['name']}", use_container_width=True):
            st.success(f"Downloading {item['name']}...")
    with col3:
        if st.button("⭐ Star", key=f"star_{item['name']}", use_container_width=True):
            st.success("Starred!")


def render_featured_collections():
    """Featured collections section"""
    st.markdown("### 📚 Featured Collections")

    collections = mock_featured_collections()

    cols = st.columns(3)
    for idx, collection in enumerate(collections):
        with cols[idx]:
            st.markdown(f"""
            <div style="
                border: 1px solid #e0e0e0;
                border-radius: 12px;
                padding: 20px;
                background: white;
                height: 180px;
            ">
                <h5 style="margin: 0 0 8px 0;">{collection['name']}</h5>
                <p style="color: #6b7280; font-size: 13px; margin: 0 0 8px 0;">
                    by @{collection['author']}
                </p>
                <p style="color: #374151; font-size: 13px; margin: 0 0 16px 0;">
                    {collection['description']}
                </p>
                <div style="color: #9ca3af; font-size: 12px;">
                    📦 {collection['items']} items • ⭐ {collection['stars']} stars
                </div>
            </div>
            """, unsafe_allow_html=True)

            if st.button("View Collection", key=f"col_{idx}", use_container_width=True):
                st.info(f"Opening {collection['name']}...")


def render_categories():
    """Category quick links"""
    st.markdown("### 🗂️ Browse by Category")

    categories = [
        ("🤖 Docking Models", "234 items"),
        ("💊 ADMET Predictors", "189 items"),
        ("🔬 Retrosynthesis", "156 items"),
        ("📊 Binding Datasets", "342 items"),
        ("🧬 Protein Structure", "98 items"),
        ("🔌 Custom Adapters", "189 items")
    ]

    cols = st.columns(3)
    for idx, (category, count) in enumerate(categories):
        with cols[idx % 3]:
            if st.button(f"{category}\n{count}", key=f"cat_{idx}", use_container_width=True):
                st.info(f"Browsing {category}...")


def render_cta_section():
    """Call to action for contributors"""
    st.markdown("---")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 40px;
            border-radius: 12px;
            color: white;
        ">
            <h3>🚀 Share Your Work</h3>
            <p>Upload your models, datasets, and adapters to help the community</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Upload to Marketplace", use_container_width=True):
            st.info("Upload functionality coming soon!")

    with col2:
        st.markdown("""
        <div style="
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            padding: 40px;
            border-radius: 12px;
            color: white;
        ">
            <h3>🏆 Submit to Leaderboard</h3>
            <p>Test your models on standard benchmarks and compete</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("View Leaderboards", use_container_width=True):
            st.info("Leaderboards coming soon!")


def main():
    """Main marketplace home page"""

    # Notice banner
    st.info("📋 **This is a preview mockup.** The PharmForge Marketplace is planned for a future release. This page demonstrates the design and user experience.")

    # Hero section
    render_hero_section()

    # Stats bar
    render_stats_bar()

    st.markdown("---")

    # Trending section
    render_trending_section()

    st.markdown("---")

    # Featured collections
    render_featured_collections()

    st.markdown("---")

    # Categories
    render_categories()

    st.markdown("---")

    # Call to action
    render_cta_section()

    # Footer
    st.markdown("---")
    st.markdown("""
    <div style="text-align: center; color: #9ca3af; padding: 20px;">
        <p>PharmForge Marketplace • The Hugging Face of Drug Discovery</p>
        <p style="font-size: 12px;">
            This is a mockup preview. Real marketplace coming post-MVP.
        </p>
    </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
