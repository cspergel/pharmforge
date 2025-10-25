"""
Unit tests for ChEMBL adapter
"""
import pytest
from adapters.chembl.adapter import ChEMBLAdapter
from backend.core.adapters.protocol import AdapterResult

class TestChEMBLAdapter:
    def test_initialization(self):
        adapter = ChEMBLAdapter()
        assert adapter.name == "chembl"
        assert adapter.adapter_type == "api"
        assert adapter.version == "1.0.0"
    
    def test_validation(self):
        adapter = ChEMBLAdapter()
        assert adapter.validate_input("CCO") is True
        assert adapter.validate_input("") is False
        assert adapter.validate_input(None) is False
    
    @pytest.mark.asyncio
    async def test_execute(self):
        adapter = ChEMBLAdapter()
        result = await adapter.execute("CCO")
        assert isinstance(result, AdapterResult)
        assert result.success is True

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
