import streamlit as st
import requests
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import time
import os
import numpy as np

# Import custom components
from components.api_client import get_api_client
from components.molecule_viewer import molecule_viewer_component, molecule_grid
from components.progress_tracker import render_progress_view, render_run_metrics, render_stage_progress

# Import pages
from pages.compound_testing import render_compound_testing_page
from pages.adapter_browser import adapter_browser_page as show_adapter_browser_page
from pages.marketplace_home import main as show_marketplace_page

# ============================================================================
# DESIGN CONFIGURATION (from PharmForge_Frontend_Design_Spec.md)
# ============================================================================

# Configure page with branding
st.set_page_config(
    page_title="PharmForge - Drug Discovery Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS - Implementing design system fundamentals
st.markdown("""
<style>
    /* Import Inter font (Design System) */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

    /* Design System Variables */
    :root {
        --color-primary-500: #00bcd4;
        --color-primary-600: #00acc1;
        --color-secondary-500: #9c27b0;
        --color-success-500: #4caf50;
        --color-warning-500: #ffc107;
        --color-error-500: #f44336;
        --color-gray-50: #fafafa;
        --color-gray-100: #f5f5f5;
        --color-gray-300: #e0e0e0;
        --color-gray-700: #616161;
        --color-gray-900: #212121;
    }

    /* Global Typography */
    html, body, [class*="css"] {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    }

    /* Headers (from Design System) */
    h1 {
        font-size: 2.25rem;
        font-weight: 700;
        color: var(--color-gray-900);
        margin-bottom: 1.5rem;
    }

    h2 {
        font-size: 1.875rem;
        font-weight: 600;
        color: var(--color-gray-900);
        margin-bottom: 1rem;
    }

    h3 {
        font-size: 1.5rem;
        font-weight: 600;
        color: var(--color-primary-500);
        margin-bottom: 0.75rem;
    }

    /* Primary Button (from Component Library) */
    .stButton > button {
        background-color: var(--color-primary-500);
        color: white;
        border-radius: 6px;
        padding: 0.75rem 1.5rem;
        font-weight: 500;
        border: none;
        transition: all 200ms ease-out;
        cursor: pointer !important;
    }

    .stButton > button:hover {
        background-color: var(--color-primary-600);
        transform: translateY(-1px);
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        cursor: pointer !important;
    }

    /* Ensure normal cursor everywhere */
    *, *::before, *::after {
        cursor: default !important;
    }

    /* Buttons and clickable elements should show pointer */
    button, a, .stButton > button, [role="button"] {
        cursor: pointer !important;
    }

    /* Input fields should show text cursor */
    input, textarea, .stTextInput input, .stTextArea textarea {
        cursor: text !important;
    }

    /* Checkboxes and radio buttons */
    input[type="checkbox"], input[type="radio"] {
        cursor: pointer !important;
    }

    /* Input Fields (from Component Library) */
    .stTextInput > div > div > input,
    .stTextArea > div > div > textarea {
        background-color: var(--color-gray-50);
        border: 1px solid var(--color-gray-300);
        border-radius: 6px;
        padding: 0.75rem 1rem;
        font-size: 1rem;
        transition: all 200ms;
    }

    .stTextInput > div > div > input:focus,
    .stTextArea > div > div > textarea:focus {
        background-color: white;
        border-color: var(--color-primary-500);
        box-shadow: 0 0 0 3px rgba(0, 188, 212, 0.5);
    }

    /* Cards */
    .metric-card {
        background: white;
        padding: 1.5rem;
        border-radius: 8px;
        box-shadow: 0 1px 3px 0 rgba(0, 0, 0, 0.1);
        transition: all 200ms;
    }

    .metric-card:hover {
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        transform: translateY(-2px);
    }

    /* Progress indicators */
    .stProgress > div > div > div {
        background-color: var(--color-primary-500);
    }

    /* Sidebar */
    [data-testid="stSidebar"] {
        background-color: var(--color-gray-50);
    }

    /* Status badges */
    .status-badge {
        display: inline-block;
        padding: 0.25rem 0.75rem;
        border-radius: 9999px;
        font-size: 0.875rem;
        font-weight: 500;
    }

    .status-running {
        background-color: var(--color-warning-500);
        color: white;
    }

    .status-completed {
        background-color: var(--color-success-500);
        color: white;
    }

    .status-failed {
        background-color: var(--color-error-500);
        color: white;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# API CONFIGURATION
# ============================================================================

API_URL = os.getenv("PHARMFORGE_API_URL", "http://backend:8000")

# Initialize API client
api_client = get_api_client()

# ============================================================================
# MAIN APP
# ============================================================================

def main():
    """Main application entry point."""

    # Header with branding
    st.markdown("""
        <div style='text-align: center; padding: 2rem 0;'>
            <h1>🧬 PharmForge</h1>
            <p style='font-size: 1.125rem; color: var(--color-gray-700);'>
                Open-Source Drug Discovery Workflow Orchestrator
            </p>
        </div>
    """, unsafe_allow_html=True)

    # Initialize navigation state
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "🏠 Home"

    # Sidebar navigation
    with st.sidebar:
        st.markdown("### Navigation")
        st.caption("✅ = Functional | 🚧 = In Development | 📋 = Preview")
        page = st.radio(
            "Go to",
            [
                "🏠 Home",
                "🧪 Compound Testing ✅",
                "🔍 Adapter Browser ✅",
                "🚀 New Run 🚧",
                "📊 My Runs 🚧",
                "🔌 Adapters ✅",
                "🏪 Marketplace 📋",
                "📈 Analytics 🚧"
            ],
            label_visibility="collapsed",
            index=[
                "🏠 Home",
                "🧪 Compound Testing ✅",
                "🔍 Adapter Browser ✅",
                "🚀 New Run 🚧",
                "📊 My Runs 🚧",
                "🔌 Adapters ✅",
                "🏪 Marketplace 📋",
                "📈 Analytics 🚧"
            ].index(st.session_state.current_page) if st.session_state.current_page in [
                "🏠 Home",
                "🧪 Compound Testing ✅",
                "🔍 Adapter Browser ✅",
                "🚀 New Run 🚧",
                "📊 My Runs 🚧",
                "🔌 Adapters ✅",
                "🏪 Marketplace 📋",
                "📈 Analytics 🚧"
            ] else 0
        )

        # Update session state when navigation changes
        if page != st.session_state.current_page:
            st.session_state.current_page = page

        st.markdown("---")

        # System health indicator
        st.markdown("### System Health")
        health_response = api_client.get_health()
        if health_response.success:
            st.success("✅ System Online")
            health_data = health_response.data
            if health_data:
                db_status = health_data.get("db", "unknown")
                redis_status = health_data.get("redis", "unknown")
                st.caption(f"DB: {db_status} | Redis: {redis_status}")
        else:
            st.error("❌ System Offline")
            st.caption(health_response.error)

        st.markdown("---")

        # Quick stats (mock data for MVP)
        st.markdown("### Quick Stats")
        st.metric("Total Runs", "47")
        st.metric("Compounds Screened", "12,389")

    # Route to pages (clean up page names)
    page_clean = page.split(" ")[0] + " " + page.split(" ")[1] if len(page.split(" ")) > 1 else page

    if "Home" in page:
        show_home_page()
    elif "New Run" in page:
        show_new_run_page()
    elif "Compound Testing" in page:
        render_compound_testing_page()
    elif "My Runs" in page:
        show_runs_page()
    elif "Adapters" in page and "Browser" not in page:
        show_adapters_page()
    elif "Adapter Browser" in page:
        show_adapter_browser_page()
    elif "Marketplace" in page:
        show_marketplace_page()
    elif "Analytics" in page:
        show_analytics_page()

# ============================================================================
# PAGE: HOME
# ============================================================================

def show_home_page():
    """
    Home page - Overview and quick actions.

    Design: Follows "Dashboard Layout" from design spec
    """
    st.markdown("## Welcome to PharmForge")
    st.markdown("### AI-Powered Drug Discovery Workflow Orchestrator")

    # Status message
    st.info("✅ **Active Features**: Compound Testing, Adapter Browser, System Health Monitoring | 🚧 **In Development**: Pipeline Orchestration, Batch Processing")

    st.markdown("---")

    # Quick action cards - Row 1
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("""
            <div class='metric-card'>
                <h3>🧪 Compound Testing ✅</h3>
                <p>Test a single compound across multiple adapters. Quick analysis and validation.</p>
                <span style="color: #4caf50; font-size: 12px;">FULLY FUNCTIONAL</span>
            </div>
        """, unsafe_allow_html=True)
        if st.button("Test Compound", key="home_compound", type="primary"):
            st.rerun()

    with col2:
        st.markdown("""
            <div class='metric-card'>
                <h3>🔍 Adapter Browser ✅</h3>
                <p>Browse all 39+ drug discovery adapters. View details, test functionality, and explore capabilities.</p>
                <span style="color: #4caf50; font-size: 12px;">FULLY FUNCTIONAL</span>
            </div>
        """, unsafe_allow_html=True)
        if st.button("Browse Adapters", key="home_adapters", type="primary"):
            st.rerun()

    with col3:
        st.markdown("""
            <div class='metric-card'>
                <h3>💬 Natural Language 🚧</h3>
                <p>Describe your goal in plain English. PharmForge plans the pipeline automatically.</p>
                <span style="color: #ffc107; font-size: 12px;">COMING SOON</span>
            </div>
        """, unsafe_allow_html=True)
        if st.button("View Demo", key="home_nl", disabled=False):
            st.info("Natural language pipeline orchestration is in development. Use Compound Testing or Adapter Browser for now.")

    st.markdown("---")

    # Recent runs
    st.markdown("## Recent Runs")

    # Fetch recent runs (mock for MVP)
    recent_runs = [
        {"id": 1, "name": "EGFR Inhibitors", "status": "completed", "progress": 1.0},
        {"id": 2, "name": "BRAF Kinase", "status": "running", "progress": 0.6},
        {"id": 3, "name": "CNS Penetration", "status": "completed", "progress": 1.0},
    ]

    for run in recent_runs:
        with st.expander(f"Run #{run['id']}: {run['name']}"):
            col1, col2 = st.columns([3, 1])

            with col1:
                if run['status'] == 'running':
                    st.progress(run['progress'])
                else:
                    st.success("✅ Completed")

            with col2:
                if st.button("View", key=f"view_{run['id']}"):
                    st.info(f"Viewing run #{run['id']}")

# ============================================================================
# PAGE: NEW RUN
# ============================================================================

def show_new_run_page():
    """
    New run creation page.

    Design: Follows "Input Modes" and "Natural Language Input" from design spec
    """
    st.markdown("## Create New Run")

    st.warning("🚧 **This feature is in development.** Full pipeline orchestration will be available soon. For now, use **Compound Testing** to test individual compounds across adapters.")

    st.markdown("---")

    # Run name
    run_name = st.text_input(
        "Run Name",
        value="My Drug Discovery Run",
        help="Give your run a descriptive name"
    )

    # Input mode tabs (Progressive Disclosure principle)
    tab1, tab2, tab3 = st.tabs([
        "💬 Natural Language",
        "📦 Batch Upload",
        "🧬 Evolution"
    ])

    with tab1:
        st.markdown("""
            ### Natural Language Query
            Describe your drug discovery goal in plain English.
            PharmForge will automatically plan and execute the pipeline.
        """)

        # NL input (from design spec)
        nl_query = st.text_area(
            "What would you like to discover?",
            placeholder="Example: Design CNS-penetrant EGFR inhibitors with good oral bioavailability and minimal hERG liability",
            height=150,
            help="Be specific about target, properties, and constraints"
        )

        # Quick templates (from design spec)
        st.markdown("**💡 Templates:**")
        template_cols = st.columns(3)

        with template_cols[0]:
            if st.button("🎯 Kinase Inhibitor", key="template1"):
                nl_query = "Design selective kinase inhibitors for EGFR with good BBB penetration"
                st.rerun()

        with template_cols[1]:
            if st.button("🧠 CNS Drug", key="template2"):
                nl_query = "Find CNS-penetrant compounds targeting NMDA receptors with low toxicity"
                st.rerun()

        with template_cols[2]:
            if st.button("💊 Oral Drug", key="template3"):
                nl_query = "Design orally bioavailable GPCR ligands with good PK properties"
                st.rerun()

        st.markdown("---")

        # Advanced options (collapsible - Progressive Disclosure)
        with st.expander("⚙️ Advanced Options"):
            col1, col2 = st.columns(2)

            with col1:
                target_protein = st.text_input("Target Protein", value="EGFR")
                n_candidates = st.slider("Top Candidates", 5, 50, 10)

            with col2:
                use_caching = st.checkbox("Use Caching", value=True)
                timeout_mins = st.slider("Timeout (minutes)", 15, 120, 60)

        # Launch button (Primary CTA)
        st.markdown("<br>", unsafe_allow_html=True)
        col1, col2, col3 = st.columns([2, 1, 2])

        with col2:
            if st.button("🚀 Launch Pipeline", type="primary", use_container_width=True):
                if nl_query:
                    with st.spinner("Starting pipeline..."):
                        # Create run via API client
                        response = api_client.create_run(
                            name=run_name,
                            input_type="nl",
                            nl_query=nl_query,
                            target_protein=target_protein,
                            n_candidates=n_candidates
                        )

                        if response.success:
                            run = response.data
                            st.success(f"✅ Run created! ID: {run['run_id']}")
                            st.info("Pipeline is running in the background. Check 'My Runs' for progress.")
                            time.sleep(2)
                            st.rerun()
                        else:
                            st.error(f"❌ Error: {response.error}")
                else:
                    st.warning("⚠️ Please enter a query")

    with tab2:
        st.markdown("""
            ### Batch Upload
            Upload a CSV file with SMILES strings to process multiple compounds.
        """)

        uploaded_file = st.file_uploader(
            "Upload SMILES file",
            type=["csv"],
            help="CSV should have a 'smiles' column"
        )

        if uploaded_file:
            df = pd.read_csv(uploaded_file)
            st.success(f"✅ Loaded {len(df)} compounds")
            st.dataframe(df.head(10), use_container_width=True)

            if st.button("🚀 Process Batch", type="primary"):
                st.info("Batch processing will be implemented in Phase 3")

    with tab3:
        st.markdown("""
            ### Evolution Mode
            Use genetic algorithms to evolve better drug candidates.
        """)

        st.info("🧬 Evolution mode coming soon in Phase 3")

        # Placeholder for evolution config
        col1, col2 = st.columns(2)

        with col1:
            st.number_input("Population Size", value=100, disabled=True)
            st.number_input("Generations", value=50, disabled=True)

        with col2:
            st.number_input("Mutation Rate", value=0.1, disabled=True)
            st.number_input("Crossover Rate", value=0.7, disabled=True)

# ============================================================================
# PAGE: MY RUNS
# ============================================================================

def show_runs_page():
    """
    Runs dashboard page.

    Design: Follows "Progress Tracking View" from design spec
    """
    st.markdown("## My Runs")

    st.warning("🚧 **This feature is in development.** Run management and tracking will be available once pipeline orchestration is complete.")

    # Fetch runs from API
    response = api_client.list_runs()

    if response.success:
        runs_data = response.data
        if runs_data and isinstance(runs_data, dict):
            runs = runs_data.get("runs", [])

            if runs and len(runs) > 0:
                # Filter options
                col1, col2, col3 = st.columns([2, 1, 1])

                with col1:
                    search = st.text_input("🔍 Search runs", placeholder="Search by name...")

                with col2:
                    status_filter = st.selectbox("Status", ["All", "Running", "Completed", "Failed"])

                with col3:
                    sort_by = st.selectbox("Sort by", ["Recent", "Name", "Status"])

                st.markdown("---")

                # Display runs
                for run in runs:
                    # Status badge
                    run_id = run.get('id', run.get('run_id', 'unknown'))
                    run_name = run.get('name', 'Unnamed Run')
                    run_status = run.get('status', 'unknown')

                    status_emoji = {
                        "running": "⏳",
                        "completed": "✅",
                        "failed": "❌",
                        "queued": "📝"
                    }.get(run_status, "❓")

                    with st.expander(f"Run #{run_id}: {run_name} {status_emoji}"):
                        # Metrics row
                        col1, col2, col3, col4 = st.columns(4)

                        with col1:
                            st.metric("Status", run_status.upper())

                        with col2:
                            progress_pct = run.get('progress', 0) * 100
                            st.metric("Progress", f"{progress_pct:.0f}%")

                        with col3:
                            n_compounds = run.get('n_compounds_processed', 0)
                            st.metric("Compounds", n_compounds)

                        with col4:
                            duration = run.get('duration_seconds', 0)
                            st.metric("Duration", f"{duration//60}m {duration%60}s")

                        # Progress bar for running jobs
                        if run_status == 'running':
                            st.progress(run.get('progress', 0))

                            # Current step indicator
                            current_step = run.get('current_step', 'Initializing')
                            st.info(f"🔄 Current step: {current_step}")

                        # Action buttons
                        button_cols = st.columns(4)

                        with button_cols[0]:
                            if run_status == 'completed':
                                if st.button("📊 View Results", key=f"results_{run_id}"):
                                    show_results_page(run_id)

                        with button_cols[1]:
                            if st.button("📄 Details", key=f"details_{run_id}"):
                                st.json(run)

                        with button_cols[2]:
                            if st.button("🔄 Refresh", key=f"refresh_{run_id}"):
                                st.rerun()

                        with button_cols[3]:
                            if run_status == 'running':
                                if st.button("⏹️ Cancel", key=f"cancel_{run_id}"):
                                    st.warning("Cancel functionality coming soon")
            else:
                st.info("📝 No runs yet. Create one in the 'New Run' page!")
        else:
            st.info("📝 No runs yet. Create one in the 'New Run' page!")
    else:
        st.error(f"❌ Failed to fetch runs: {response.error}")

# ============================================================================
# PAGE: RESULTS
# ============================================================================

def show_results_page(run_id: int):
    """
    Results display page.

    Design: Follows "Results Table" and "Pareto Plot" from design spec
    """
    st.markdown(f"## Results for Run #{run_id}")

    # Fetch results from API
    response = api_client.get_run_results(run_id, top_n=20)

    if response.success:
        results = response.data

        # Summary metrics
        st.markdown("### Summary")
        col1, col2, col3, col4 = st.columns(4)

        with col1:
            st.metric("Top Compounds", len(results.get('top_compounds', [])))

        with col2:
            st.metric("Compounds Screened", results.get('n_input', 0))

        with col3:
            avg_score = 0.87  # Mock
            st.metric("Avg Score", f"{avg_score:.2f}")

        with col4:
            st.metric("Runtime", "24m 13s")  # Mock

        st.markdown("---")

        # Tabs for different views
        tab1, tab2, tab3 = st.tabs([
            "📊 Ranked List",
            "📈 Pareto Plot",
            "🔬 Compound Details"
        ])

        with tab1:
            # Results table
            st.markdown("### Top Candidates")

            if 'top_compounds' in results:
                df = pd.DataFrame(results['top_compounds'])

                # Format scores
                for col in ['binding_score', 'admet_score', 'synthesis_score', 'novelty_score', 'composite_score']:
                    if col in df.columns:
                        df[col] = df[col].apply(lambda x: f"{x:.3f}")

                # Display table
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=False
                )

                # Export button
                csv = df.to_csv(index=False)
                st.download_button(
                    label="📥 Download CSV",
                    data=csv,
                    file_name=f"pharmforge_run_{run_id}_results.csv",
                    mime="text/csv"
                )
            else:
                st.info("No results available yet")

        with tab2:
            # Pareto plot (2D: Binding vs ADMET)
            st.markdown("### Pareto Frontier")

            if 'top_compounds' in results:
                df = pd.DataFrame(results['top_compounds'])

                # Create scatter plot
                fig = px.scatter(
                    df,
                    x='binding_score',
                    y='admet_score',
                    color='pareto_rank',
                    size='composite_score',
                    hover_data=['smiles', 'synthesis_score', 'novelty_score'],
                    title="Multi-Objective Optimization Space",
                    labels={
                        'binding_score': 'Binding Affinity',
                        'admet_score': 'ADMET Score',
                        'pareto_rank': 'Pareto Rank'
                    },
                    color_continuous_scale='Viridis'
                )

                fig.update_layout(
                    height=600,
                    hovermode='closest'
                )

                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No data available for plotting")

        with tab3:
            # Detailed compound viewer
            st.markdown("### Compound Details")

            if 'top_compounds' in results:
                df = pd.DataFrame(results['top_compounds'])

                # Compound selector
                selected_idx = st.selectbox(
                    "Select compound",
                    range(len(df)),
                    format_func=lambda i: f"Rank {i+1}: {df.iloc[i]['smiles']}"
                )

                if selected_idx is not None:
                    compound = df.iloc[selected_idx]

                    # Display details
                    col1, col2 = st.columns([1, 2])

                    with col1:
                        st.markdown(f"**SMILES:** `{compound['smiles']}`")
                        st.markdown(f"**Rank:** {selected_idx + 1}")
                        st.markdown(f"**Pareto Rank:** {compound.get('pareto_rank', 'N/A')}")

                    with col2:
                        # Score breakdown
                        scores = {
                            'Binding': float(compound.get('binding_score', 0)),
                            'ADMET': float(compound.get('admet_score', 0)),
                            'Synthesis': float(compound.get('synthesis_score', 0)),
                            'Novelty': float(compound.get('novelty_score', 0))
                        }

                        fig = go.Figure(data=[
                            go.Bar(
                                x=list(scores.keys()),
                                y=list(scores.values()),
                                marker_color='#00bcd4'
                            )
                        ])

                        fig.update_layout(
                            title="Score Breakdown",
                            yaxis_title="Score",
                            height=300
                        )

                        st.plotly_chart(fig, use_container_width=True)

                    # Molecule viewer using our component
                    st.markdown("---")
                    st.markdown("#### 2D Structure")
                    molecule_viewer_component(compound['smiles'], width=400, height=400)
            else:
                st.info("No compound data available")
    else:
        st.error(f"❌ Failed to fetch results: {response.error}")

# ============================================================================
# PAGE: ADAPTERS
# ============================================================================

def show_adapters_page():
    """
    Adapters dashboard page.

    Design: Follows "Adapter Health" monitoring from design spec
    """
    st.markdown("## Available Adapters")

    # System health
    try:
        response = requests.get(f"{API_URL}/health", timeout=2)

        if response.status_code == 200:
            health = response.json()

            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric("System Status", health.get('status', 'unknown').upper())

            with col2:
                st.metric("Database", health.get('database', 'unknown'))

            with col3:
                st.metric("Redis", health.get('redis', 'unknown'))
        else:
            st.error("⚠️ System health check failed")

    except Exception as e:
        st.error(f"❌ Connection error: {str(e)}")

    st.markdown("---")

    # Adapter categories
    st.markdown("### Adapter Categories")

    # Mock adapter data (will come from API in production)
    adapters = [
        {
            "name": "vina_docking",
            "category": "Docking",
            "status": "healthy",
            "version": "1.2.3",
            "description": "AutoDock Vina molecular docking"
        },
        {
            "name": "tdc_admet",
            "category": "ADMET",
            "status": "healthy",
            "version": "0.4.1",
            "description": "Therapeutics Data Commons ADMET predictions"
        },
        {
            "name": "aizynthfinder",
            "category": "Retrosynthesis",
            "status": "degraded",
            "version": "4.2.0",
            "description": "Retrosynthetic route planning"
        },
        {
            "name": "chembl_similarity",
            "category": "Novelty",
            "status": "healthy",
            "version": "1.0.0",
            "description": "ChEMBL similarity search for novelty scoring"
        }
    ]

    # Group by category
    categories = {}
    for adapter in adapters:
        cat = adapter['category']
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(adapter)

    # Display adapters by category
    for category, adapters_list in categories.items():
        st.markdown(f"#### {category}")

        for adapter in adapters_list:
            with st.expander(f"{adapter['name']} v{adapter['version']}"):
                col1, col2 = st.columns([2, 1])

                with col1:
                    st.markdown(f"**Description:** {adapter['description']}")
                    st.markdown(f"**Category:** {adapter['category']}")

                with col2:
                    status_emoji = {
                        "healthy": "✅",
                        "degraded": "⚠️",
                        "offline": "❌"
                    }.get(adapter['status'], "❓")

                    st.metric("Status", f"{status_emoji} {adapter['status'].upper()}")

                # Action buttons
                button_cols = st.columns(3)

                with button_cols[0]:
                    if st.button("🔍 Test", key=f"test_{adapter['name']}"):
                        st.info("Running health check...")

                with button_cols[1]:
                    if st.button("📊 Stats", key=f"stats_{adapter['name']}"):
                        st.info("Adapter statistics coming soon")

                with button_cols[2]:
                    if st.button("⚙️ Config", key=f"config_{adapter['name']}"):
                        st.info("Adapter configuration coming soon")

# ============================================================================
# PAGE: ANALYTICS
# ============================================================================

def show_analytics_page():
    """Analytics dashboard (placeholder for MVP)."""
    st.markdown("## Analytics")

    st.warning("🚧 **This feature is in development.** Usage analytics and insights will be available in a future release.")

    st.markdown("---")

    # Placeholder charts
    st.markdown("### Usage Over Time")

    # Mock data
    dates = pd.date_range(start='2025-01-01', end='2025-01-31', freq='D')
    runs_per_day = np.random.randint(5, 25, size=len(dates))

    df = pd.DataFrame({
        'Date': dates,
        'Runs': runs_per_day
    })

    fig = px.line(df, x='Date', y='Runs', title='Daily Runs')
    st.plotly_chart(fig, use_container_width=True)

    st.markdown("### Popular Targets")

    targets = ['EGFR', 'BRAF', 'ALK', 'RET', 'MET']
    counts = [45, 32, 28, 19, 15]

    fig = px.bar(
        x=targets,
        y=counts,
        labels={'x': 'Target', 'y': 'Number of Runs'},
        title='Most Searched Targets'
    )
    st.plotly_chart(fig, use_container_width=True)

# ============================================================================
# RUN APP
# ============================================================================

if __name__ == "__main__":
    main()
