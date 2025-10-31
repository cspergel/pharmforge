"""
Progress Tracker Component
Displays pipeline execution progress with stage indicators
"""
import streamlit as st
from typing import List, Dict, Optional
from datetime import datetime


def format_duration(seconds: Optional[float]) -> str:
    """
    Format duration in seconds to human-readable string

    Args:
        seconds: Duration in seconds

    Returns:
        Formatted string (e.g., "2m 34s")
    """
    if seconds is None or seconds == 0:
        return "-"

    minutes = int(seconds // 60)
    secs = int(seconds % 60)

    if minutes > 0:
        return f"{minutes}m {secs}s"
    else:
        return f"{secs}s"


def render_stage_progress(
    stages: List[Dict],
    current_stage_idx: Optional[int] = None
):
    """
    Render pipeline stage progress indicators

    Args:
        stages: List of stage dictionaries with keys:
            - name: Stage name
            - status: "completed", "running", "queued", "failed"
            - time: Optional completion time
            - eta: Optional ETA for queued stages
            - progress: Optional progress fraction (0-1) for running stages
        current_stage_idx: Index of currently running stage
    """
    for idx, stage in enumerate(stages):
        col1, col2, col3 = st.columns([1, 8, 2])

        status = stage.get("status", "queued")
        name = stage.get("name", f"Stage {idx+1}")

        # Status icon
        with col1:
            if status == "completed":
                st.success("✅")
            elif status == "running":
                st.info("🔄")
            elif status == "failed":
                st.error("❌")
            else:  # queued
                st.warning("⏳")

        # Stage name and progress
        with col2:
            # Highlight current stage
            if idx == current_stage_idx:
                st.markdown(f"**{idx+1}. {name}** (current)")
            else:
                st.markdown(f"{idx+1}. {name}")

            # Show progress bar for running stages
            if status == "running" and "progress" in stage:
                st.progress(stage["progress"])

        # Time/ETA
        with col3:
            if "time" in stage:
                st.caption(format_duration(stage["time"]))
            elif "eta" in stage:
                st.caption(f"ETA {stage['eta']}")
            else:
                st.caption("")


def render_run_metrics(
    n_compounds_processed: int = 0,
    n_compounds_total: int = 0,
    top_score: Optional[float] = None,
    elapsed_time: Optional[float] = None
):
    """
    Render quick stats for a running pipeline

    Args:
        n_compounds_processed: Number of compounds processed
        n_compounds_total: Total number of compounds
        top_score: Best score so far
        elapsed_time: Time elapsed in seconds
    """
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        progress_pct = (n_compounds_processed / n_compounds_total * 100) if n_compounds_total > 0 else 0
        st.metric(
            "Progress",
            f"{n_compounds_processed} / {n_compounds_total}",
            f"{progress_pct:.0f}%"
        )

    with col2:
        if top_score is not None:
            st.metric("Top Score", f"{top_score:.3f}")
        else:
            st.metric("Top Score", "-")

    with col3:
        if elapsed_time is not None:
            st.metric("Elapsed", format_duration(elapsed_time))
        else:
            st.metric("Elapsed", "-")

    with col4:
        # Placeholder for future metrics
        st.metric("Status", "Running")


def render_progress_view(run_data: Dict):
    """
    Render complete progress view for a pipeline run

    Args:
        run_data: Dictionary containing run information:
            - name: Run name
            - run_id: Unique identifier
            - status: Current status
            - stages: List of stage dictionaries
            - metrics: Optional metrics dictionary
    """
    # Header
    st.markdown(f"### 🏃 Running: {run_data['name']}")
    st.caption(f"Run ID: {run_data['run_id']} • Status: {run_data['status']}")

    st.markdown("---")

    # Metrics (if available)
    if "metrics" in run_data:
        metrics = run_data["metrics"]
        render_run_metrics(
            n_compounds_processed=metrics.get("n_processed", 0),
            n_compounds_total=metrics.get("n_total", 0),
            top_score=metrics.get("top_score"),
            elapsed_time=metrics.get("elapsed_time")
        )

        st.markdown("---")

    # Stage progress
    st.markdown("### Pipeline Stages")

    stages = run_data.get("stages", [])
    current_stage_idx = run_data.get("current_stage_idx")

    if stages:
        render_stage_progress(stages, current_stage_idx)
    else:
        st.info("No stage information available")

    # Auto-refresh button
    st.markdown("---")
    if st.button("🔄 Refresh Status", use_container_width=True):
        st.rerun()


# Example usage for testing
def _example_progress_view():
    """Example progress view for testing"""
    example_data = {
        "name": "EGFR Kinase Inhibitors",
        "run_id": "run_abc123",
        "status": "running",
        "current_stage_idx": 2,
        "stages": [
            {"name": "Target Validation", "status": "completed", "time": 2.5},
            {"name": "Scaffold Generation", "status": "completed", "time": 45.2},
            {"name": "ADMET Filtering", "status": "running", "progress": 0.67},
            {"name": "Molecular Docking", "status": "queued", "eta": "2h 15m"},
            {"name": "Ranking", "status": "queued"},
        ],
        "metrics": {
            "n_processed": 247,
            "n_total": 500,
            "top_score": 0.89,
            "elapsed_time": 123.5
        }
    }

    render_progress_view(example_data)


if __name__ == "__main__":
    # Run example when testing
    _example_progress_view()
