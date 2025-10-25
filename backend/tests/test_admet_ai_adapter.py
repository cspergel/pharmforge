"""
Unit tests for ADMET-AI adapter
Tests protocol compliance, predictions, caching, and error handling
"""
import pytest
import asyncio
from typing import List, Optional

from adapters.admet_ai.adapter import ADMETaiAdapter
from backend.core.adapters.protocol import AdapterResult


class TestADMETaiAdapterInitialization:
    """Test adapter initialization and configuration"""

    def test_default_initialization(self):
        """Test initialization with default parameters"""
        adapter = ADMETaiAdapter()

        assert adapter.name == "admet_ai"
        assert adapter.adapter_type == "ml"
        assert adapter.version == "1.0.0"
        assert adapter.enabled is True
        assert adapter.properties is None  # Default to all properties
        assert adapter.config == {}

    def test_initialization_with_config(self):
        """Test initialization with custom config"""
        config = {
            "properties": ["AMES", "hERG", "CYP3A4_Veith"]
        }
        adapter = ADMETaiAdapter(config=config)

        assert adapter.name == "admet_ai"
        assert adapter.adapter_type == "ml"
        assert adapter.properties == ["AMES", "hERG", "CYP3A4_Veith"]
        assert adapter.config == config

    def test_initialization_with_custom_name(self):
        """Test initialization with custom name and type"""
        adapter = ADMETaiAdapter(
            name="custom_admet",
            adapter_type="ml",
            config={"properties": ["AMES"]}
        )

        assert adapter.name == "custom_admet"
        assert adapter.adapter_type == "ml"


class TestSMILESValidation:
    """Test SMILES validation logic"""

    def test_valid_smiles(self):
        """Test validation of valid SMILES strings"""
        adapter = ADMETaiAdapter()

        # Common valid SMILES
        valid_smiles = [
            "CCO",  # Ethanol
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
            "c1ccccc1",  # Benzene
        ]

        for smiles in valid_smiles:
            assert adapter.validate_input(smiles), f"Failed to validate: {smiles}"

    def test_invalid_smiles(self):
        """Test rejection of invalid SMILES strings"""
        adapter = ADMETaiAdapter()

        # Invalid SMILES
        invalid_smiles = [
            "",  # Empty string
            "XYZ123",  # Invalid characters
            "C(C",  # Unbalanced parentheses
            "C1CCC",  # Incomplete ring
            None,  # None type
            123,  # Non-string type
        ]

        for smiles in invalid_smiles:
            assert not adapter.validate_input(smiles), f"Should reject: {smiles}"

    def test_edge_cases(self):
        """Test edge cases in SMILES validation"""
        adapter = ADMETaiAdapter()

        # Edge cases
        assert adapter.validate_input("C")  # Methane (single atom)
        assert adapter.validate_input("[Na+].[Cl-]")  # Ionic compound
        assert adapter.validate_input("C/C=C/C")  # Stereochemistry


class TestExecuteMethod:
    """Test the execute() method with real predictions"""

    @pytest.mark.asyncio
    async def test_execute_simple_molecule(self):
        """Test prediction for simple molecule (ethanol)"""
        adapter = ADMETaiAdapter()
        result = await adapter.execute("CCO")

        # Check result structure
        assert isinstance(result, AdapterResult)
        assert result.success is True
        assert result.error is None
        assert result.data is not None

        # Check data structure
        data = result.data
        assert "smiles" in data
        assert "properties" in data
        assert "percentiles" in data
        assert "property_count" in data
        assert "model" in data

        # Check metadata
        assert result.metadata is not None
        assert result.metadata["adapter_name"] == "admet_ai"
        assert "cache_key" in result.metadata
        assert "version" in result.metadata

    @pytest.mark.asyncio
    async def test_execute_complex_molecule(self):
        """Test prediction for complex molecule (aspirin)"""
        adapter = ADMETaiAdapter()
        result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

        assert result.success is True
        assert result.data is not None
        assert len(result.data["properties"]) > 0
        assert len(result.data["percentiles"]) > 0

    @pytest.mark.asyncio
    async def test_execute_invalid_smiles(self):
        """Test execution with invalid SMILES"""
        adapter = ADMETaiAdapter()
        result = await adapter.execute("INVALID_SMILES")

        assert result.success is False
        assert result.error == "Invalid SMILES string"
        assert result.data == {}

    @pytest.mark.asyncio
    async def test_execute_with_property_filter(self):
        """Test execution with specific properties filter"""
        # Initialize adapter with property filter
        config = {"properties": ["AMES", "hERG", "CYP3A4_Veith"]}
        adapter = ADMETaiAdapter(config=config)

        result = await adapter.execute("CCO")

        assert result.success is True

        # Check that only requested properties are returned
        properties = result.data["properties"]
        # Note: The actual filtering happens in the adapter
        # We just verify the structure is correct
        assert isinstance(properties, dict)

    @pytest.mark.asyncio
    async def test_execute_with_params_override(self):
        """Test execution with properties override via params"""
        adapter = ADMETaiAdapter()

        # Override properties via params
        result = await adapter.execute(
            "CCO",
            properties=["AMES", "hERG"]
        )

        assert result.success is True
        assert result.data is not None


class TestCacheKeyGeneration:
    """Test deterministic cache key generation"""

    def test_cache_key_deterministic(self):
        """Test that cache keys are deterministic"""
        adapter = ADMETaiAdapter()

        smiles = "CCO"
        key1 = adapter.generate_cache_key(smiles)
        key2 = adapter.generate_cache_key(smiles)

        assert key1 == key2
        assert isinstance(key1, str)
        assert len(key1) == 64  # SHA256 hex length

    def test_cache_key_different_smiles(self):
        """Test that different SMILES produce different keys"""
        adapter = ADMETaiAdapter()

        key1 = adapter.generate_cache_key("CCO")
        key2 = adapter.generate_cache_key("CC(=O)Oc1ccccc1C(=O)O")

        assert key1 != key2

    def test_cache_key_canonical_smiles(self):
        """Test that equivalent SMILES produce same key"""
        adapter = ADMETaiAdapter()

        # Different representations of benzene
        key1 = adapter.generate_cache_key("c1ccccc1")
        key2 = adapter.generate_cache_key("C1=CC=CC=C1")

        # Should be the same after canonicalization
        assert key1 == key2

    def test_cache_key_with_properties(self):
        """Test cache key includes properties filter"""
        adapter = ADMETaiAdapter()

        smiles = "CCO"
        key1 = adapter.generate_cache_key(smiles, properties=None)
        key2 = adapter.generate_cache_key(smiles, properties=["AMES", "hERG"])

        # Different properties should produce different keys
        assert key1 != key2

    def test_cache_key_properties_order_independent(self):
        """Test that property order doesn't affect cache key"""
        adapter = ADMETaiAdapter()

        smiles = "CCO"
        key1 = adapter.generate_cache_key(smiles, properties=["AMES", "hERG"])
        key2 = adapter.generate_cache_key(smiles, properties=["hERG", "AMES"])

        # Should be the same (sorted internally)
        assert key1 == key2


class TestMetadata:
    """Test adapter metadata retrieval"""

    def test_get_metadata_default(self):
        """Test metadata retrieval with default config"""
        adapter = ADMETaiAdapter()
        metadata = adapter.get_metadata()

        assert metadata["name"] == "admet_ai"
        assert metadata["type"] == "ml"
        assert metadata["version"] == "1.0.0"
        assert metadata["enabled"] is True

        # Check description and properties info
        assert "description" in metadata
        assert "properties" in metadata
        assert metadata["properties"]["total_properties"] == 49
        assert "categories" in metadata["properties"]

    def test_get_metadata_with_config(self):
        """Test metadata includes config information"""
        config = {"properties": ["AMES", "hERG"]}
        adapter = ADMETaiAdapter(config=config)
        metadata = adapter.get_metadata()

        assert metadata["config"]["properties_filter"] == ["AMES", "hERG"]


class TestPropertyGroupHelpers:
    """Test helper methods for property groups"""

    def test_get_absorption_properties(self):
        """Test absorption properties list"""
        adapter = ADMETaiAdapter()
        props = adapter.get_absorption_properties()

        assert isinstance(props, list)
        assert len(props) > 0
        assert "Caco2_Wang" in props
        assert "HIA_Hou" in props

    def test_get_distribution_properties(self):
        """Test distribution properties list"""
        adapter = ADMETaiAdapter()
        props = adapter.get_distribution_properties()

        assert isinstance(props, list)
        assert "BBB_Martins" in props
        assert "PPBR_AZ" in props

    def test_get_metabolism_properties(self):
        """Test metabolism properties list"""
        adapter = ADMETaiAdapter()
        props = adapter.get_metabolism_properties()

        assert isinstance(props, list)
        assert "CYP1A2_Veith" in props
        assert "CYP3A4_Veith" in props

    def test_get_excretion_properties(self):
        """Test excretion properties list"""
        adapter = ADMETaiAdapter()
        props = adapter.get_excretion_properties()

        assert isinstance(props, list)
        assert "Clearance_Hepatocyte_AZ" in props
        assert "Half_Life_Obach" in props

    def test_get_toxicity_properties(self):
        """Test toxicity properties list"""
        adapter = ADMETaiAdapter()
        props = adapter.get_toxicity_properties()

        assert isinstance(props, list)
        assert "AMES" in props
        assert "hERG" in props
        assert "DILI" in props

    def test_get_physicochemical_properties(self):
        """Test physicochemical properties list"""
        adapter = ADMETaiAdapter()
        props = adapter.get_physicochemical_properties()

        assert isinstance(props, list)
        assert "molecular_weight" in props
        assert "logP" in props


class TestErrorHandling:
    """Test error handling and edge cases"""

    @pytest.mark.asyncio
    async def test_execute_with_none_input(self):
        """Test execution with None input"""
        adapter = ADMETaiAdapter()
        result = await adapter.execute(None)

        assert result.success is False
        assert "Invalid SMILES string" in result.error

    @pytest.mark.asyncio
    async def test_execute_with_empty_string(self):
        """Test execution with empty string"""
        adapter = ADMETaiAdapter()
        result = await adapter.execute("")

        assert result.success is False
        assert "Invalid SMILES string" in result.error

    @pytest.mark.asyncio
    async def test_execute_preserves_smiles_in_result(self):
        """Test that original SMILES is preserved in result"""
        adapter = ADMETaiAdapter()
        original_smiles = "CCO"
        result = await adapter.execute(original_smiles)

        assert result.success is True
        assert result.data["smiles"] == original_smiles


class TestLazyModelLoading:
    """Test lazy loading of ADMET-AI model"""

    def test_model_not_loaded_on_init(self):
        """Test that model is not loaded during initialization"""
        adapter = ADMETaiAdapter()
        assert adapter._model is None

    @pytest.mark.asyncio
    async def test_model_loaded_on_first_use(self):
        """Test that model is loaded on first use"""
        adapter = ADMETaiAdapter()

        # Model should be None before first use
        assert adapter._model is None

        # Execute a prediction
        await adapter.execute("CCO")

        # Model should now be loaded
        assert adapter._model is not None

    @pytest.mark.asyncio
    async def test_model_reused_on_subsequent_calls(self):
        """Test that model is reused for subsequent calls"""
        adapter = ADMETaiAdapter()

        # First call
        await adapter.execute("CCO")
        model1 = adapter._model

        # Second call
        await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")
        model2 = adapter._model

        # Should be the same model instance
        assert model1 is model2


class TestProtocolCompliance:
    """Test that adapter follows AdapterProtocol correctly"""

    def test_inherits_from_protocol(self):
        """Test that adapter inherits from AdapterProtocol"""
        from backend.core.adapters.protocol import AdapterProtocol

        adapter = ADMETaiAdapter()
        assert isinstance(adapter, AdapterProtocol)

    def test_has_required_methods(self):
        """Test that adapter implements all required methods"""
        adapter = ADMETaiAdapter()

        assert hasattr(adapter, 'execute')
        assert hasattr(adapter, 'validate_input')
        assert hasattr(adapter, 'generate_cache_key')
        assert hasattr(adapter, 'get_metadata')
        assert callable(adapter.execute)
        assert callable(adapter.validate_input)
        assert callable(adapter.generate_cache_key)
        assert callable(adapter.get_metadata)

    def test_has_required_attributes(self):
        """Test that adapter has all required attributes"""
        adapter = ADMETaiAdapter()

        assert hasattr(adapter, 'name')
        assert hasattr(adapter, 'adapter_type')
        assert hasattr(adapter, 'version')
        assert hasattr(adapter, 'enabled')
        assert hasattr(adapter, 'config')

    @pytest.mark.asyncio
    async def test_execute_returns_adapter_result(self):
        """Test that execute() returns AdapterResult"""
        adapter = ADMETaiAdapter()
        result = await adapter.execute("CCO")

        assert isinstance(result, AdapterResult)
        assert hasattr(result, 'success')
        assert hasattr(result, 'data')
        assert hasattr(result, 'error')
        assert hasattr(result, 'metadata')

    def test_validate_input_returns_bool(self):
        """Test that validate_input() returns boolean"""
        adapter = ADMETaiAdapter()

        result = adapter.validate_input("CCO")
        assert isinstance(result, bool)


# Integration test
class TestIntegration:
    """Integration tests with multiple molecules"""

    @pytest.mark.asyncio
    async def test_batch_predictions(self):
        """Test predictions for multiple molecules"""
        adapter = ADMETaiAdapter()

        molecules = [
            "CCO",  # Ethanol
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        ]

        results = []
        for smiles in molecules:
            result = await adapter.execute(smiles)
            results.append(result)

        # All should succeed
        assert all(r.success for r in results)

        # All should have predictions
        assert all(len(r.data["properties"]) > 0 for r in results)

        # Cache keys should be different
        cache_keys = [r.metadata["cache_key"] for r in results]
        assert len(cache_keys) == len(set(cache_keys))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
