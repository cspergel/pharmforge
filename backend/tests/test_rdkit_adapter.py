"""
Unit tests for RDKit adapter
"""
import pytest
from adapters.rdkit_local.adapter import RDKitAdapter
from backend.core.adapters.protocol import AdapterResult

class TestRDKitAdapter:
    def test_initialization(self):
        adapter = RDKitAdapter()
        assert adapter.name == "rdkit_local"
        assert adapter.adapter_type == "local"
        assert adapter.version == "1.0.0"
    
    def test_validation(self):
        adapter = RDKitAdapter()
        assert adapter.validate_input("CCO") is True
        assert adapter.validate_input("") is False
        assert adapter.validate_input(None) is False
    
    @pytest.mark.asyncio
    async def test_execute(self):
        adapter = RDKitAdapter()
        result = await adapter.execute("CCO")
        assert isinstance(result, AdapterResult)
        assert result.success is True
        assert "molecular_weight" in result.data

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
