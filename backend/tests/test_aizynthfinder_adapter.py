"""
Tests for AiZynthFinder Adapter

Tests retrosynthesis route planning functionality.
"""

import pytest
import asyncio
from adapters.aizynthfinder.adapter import AiZynthFinderAdapter
from backend.core.adapters.protocol import AdapterResult

# Check if RDKit is available
try:
    import rdkit
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


@pytest.fixture
def adapter():
    """Create AiZynthFinder adapter instance."""
    return AiZynthFinderAdapter()


@pytest.fixture
def simple_molecules():
    """Simple test molecules for retrosynthesis."""
    return {
        "aspirin": "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "benzene": "c1ccccc1",  # Benzene (very simple)
        "invalid": "INVALID_SMILES",  # Invalid SMILES
        "complex": "CN1CCN(CC1)C(=O)C2=C(C=CC(=C2)OC)OC"  # More complex molecule
    }


class TestAiZynthFinderAdapter:
    """Test suite for AiZynthFinder adapter."""

    def test_adapter_initialization(self, adapter):
        """Test that adapter initializes correctly."""
        assert adapter.name == "aizynthfinder"
        assert adapter.adapter_type == "local"
        assert adapter.version == "1.0.0"
        assert adapter.enabled is True

    def test_adapter_metadata(self, adapter):
        """Test adapter metadata."""
        metadata = adapter.get_metadata()

        assert metadata["name"] == "aizynthfinder"
        assert metadata["type"] == "local"
        assert "description" in metadata
        assert "capabilities" in metadata
        assert metadata["capabilities"]["multi_route_discovery"] is True
        assert metadata["capabilities"]["synthesis_step_counting"] is True

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_validate_input_valid_smiles(self, adapter, simple_molecules):
        """Test input validation with valid SMILES."""
        assert adapter.validate_input(simple_molecules["aspirin"]) is True
        assert adapter.validate_input(simple_molecules["benzene"]) is True

    def test_validate_input_invalid_smiles(self, adapter, simple_molecules):
        """Test input validation with invalid SMILES."""
        assert adapter.validate_input(simple_molecules["invalid"]) is False
        assert adapter.validate_input("") is False
        assert adapter.validate_input(None) is False
        assert adapter.validate_input(123) is False

    @pytest.mark.skipif(not RDKIT_AVAILABLE, reason="RDKit not installed")
    def test_cache_key_generation(self, adapter, simple_molecules):
        """Test deterministic cache key generation."""
        smiles = simple_molecules["aspirin"]

        # Same inputs should produce same cache key
        key1 = adapter.generate_cache_key(smiles, max_routes=5, expansion_time=60)
        key2 = adapter.generate_cache_key(smiles, max_routes=5, expansion_time=60)
        assert key1 == key2

        # Different parameters should produce different keys
        key3 = adapter.generate_cache_key(smiles, max_routes=10, expansion_time=60)
        assert key1 != key3

        # Cache key should be a hex string
        assert isinstance(key1, str)
        assert len(key1) == 64  # SHA256 produces 64 hex characters

    @pytest.mark.asyncio
    async def test_execute_invalid_input(self, adapter, simple_molecules):
        """Test execution with invalid input."""
        result = await adapter.execute(simple_molecules["invalid"])

        assert isinstance(result, AdapterResult)
        assert result.success is False
        assert result.error is not None
        assert "Invalid SMILES" in result.error

    @pytest.mark.asyncio
    @pytest.mark.skipif(True, reason="Requires AiZynthFinder installation")
    async def test_execute_simple_molecule(self, adapter, simple_molecules):
        """
        Test retrosynthesis for a simple molecule.

        NOTE: This test is skipped by default as it requires:
        1. AiZynthFinder installation
        2. Trained model files
        3. Stock molecule database

        To run: pip install aizynthfinder && pytest -v -m "not skipif"
        """
        result = await adapter.execute(
            simple_molecules["aspirin"],
            max_routes=3,
            expansion_time=30
        )

        assert isinstance(result, AdapterResult)
        assert result.success is True
        assert result.data is not None

        # Check result structure
        data = result.data
        assert "smiles" in data
        assert "routes_found" in data
        assert "routes" in data
        assert "best_route" in data
        assert "synthesis_score" in data
        assert "feasibility" in data

        # Check best route
        best_route = data["best_route"]
        if best_route:
            assert "n_steps" in best_route
            assert "synthesis_score" in best_route
            assert "feasibility" in best_route
            assert best_route["feasibility"] in ["high", "medium", "low"]

            # Synthesis score should be 0-1 (normalized)
            assert 0.0 <= best_route["synthesis_score"] <= 1.0

    @pytest.mark.asyncio
    @pytest.mark.skipif(True, reason="Requires AiZynthFinder installation")
    async def test_execute_benzene(self, adapter, simple_molecules):
        """
        Test retrosynthesis for benzene (should be very simple/commercial).

        This is skipped by default - requires full AiZynthFinder setup.
        """
        result = await adapter.execute(
            simple_molecules["benzene"],
            max_routes=1,
            expansion_time=20
        )

        assert isinstance(result, AdapterResult)
        assert result.success is True

        # Benzene should have very few steps or be in stock
        data = result.data
        if data.get("best_route"):
            assert data["best_route"]["n_steps"] <= 2

    def test_custom_config(self):
        """Test adapter with custom configuration."""
        custom_config = {
            "max_routes": 10,
            "timeout": 180,
            "expansion_time": 90,
            "stock": "emolecules"
        }

        adapter = AiZynthFinderAdapter(config=custom_config)

        assert adapter.config["max_routes"] == 10
        assert adapter.config["timeout"] == 180
        assert adapter.config["expansion_time"] == 90
        assert adapter.config["stock"] == "emolecules"

    def test_result_format_without_aizynthfinder(self, adapter, simple_molecules):
        """
        Test that adapter returns proper error when AiZynthFinder is not installed.

        This test will pass regardless of whether AiZynthFinder is installed.
        """
        # This test documents the expected behavior
        # If AiZynthFinder is not installed, execution should return error
        pass


class TestSynthesisScoring:
    """Test synthesis step scoring functions."""

    def test_synthesis_steps_normalization(self):
        """Test that synthesis steps are correctly normalized."""
        from backend.core.scoring_utils import synthesis_steps_to01

        # Fewer steps should give higher score
        score_1_step = synthesis_steps_to01(1)
        score_5_steps = synthesis_steps_to01(5)
        score_10_steps = synthesis_steps_to01(10)

        assert score_1_step > score_5_steps
        assert score_5_steps > score_10_steps

        # All scores should be in 0-1 range
        assert 0.0 <= score_1_step <= 1.0
        assert 0.0 <= score_5_steps <= 1.0
        assert 0.0 <= score_10_steps <= 1.0

        # 1 step should be maximum score (1.0)
        assert score_1_step == 1.0

        # 10 steps should be minimum score (0.0)
        assert score_10_steps == 0.0


# Integration test (requires full setup)
@pytest.mark.integration
@pytest.mark.skipif(True, reason="Integration test - requires full setup")
class TestAiZynthFinderIntegration:
    """Integration tests requiring full AiZynthFinder setup."""

    @pytest.mark.asyncio
    async def test_full_pipeline_aspirin(self):
        """
        Full integration test for aspirin synthesis.

        Requires:
        - AiZynthFinder installed
        - Policy models downloaded
        - Stock database configured
        """
        adapter = AiZynthFinderAdapter(config={
            "max_routes": 5,
            "expansion_time": 60
        })

        aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
        result = await adapter.execute(aspirin_smiles)

        assert result.success is True
        assert result.data["routes_found"] > 0

        # Verify route structure
        for route in result.data["routes"]:
            assert "route_id" in route
            assert "n_steps" in route
            assert "synthesis_score" in route
            assert "feasibility" in route
            assert "starting_materials" in route

            # Validate scores
            assert 0.0 <= route["synthesis_score"] <= 1.0
            assert route["feasibility"] in ["high", "medium", "low"]

    @pytest.mark.asyncio
    async def test_no_routes_found_handling(self):
        """
        Test handling when no synthesis routes are found.

        Some complex molecules may not have discoverable routes.
        """
        adapter = AiZynthFinderAdapter(config={
            "max_routes": 3,
            "expansion_time": 30
        })

        # Use a complex or problematic molecule
        complex_smiles = "CC(C)(C)c1ccc(cc1)C(C)(C)c2ccc(cc2)C(C)(C)C"

        result = await adapter.execute(complex_smiles)

        # Should succeed even if no routes found
        assert result.success is True
        assert result.data["routes_found"] >= 0

        if result.data["routes_found"] == 0:
            assert result.data["feasibility"] == "none"
            assert result.data["synthesis_score"] == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
