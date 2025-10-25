"""
Tests for PubChem Adapter
"""
import pytest
import asyncio
from adapters.pubchem.adapter import PubChemAdapter


@pytest.mark.asyncio
async def test_pubchem_adapter_aspirin():
    """Test PubChem adapter with aspirin SMILES"""
    adapter = PubChemAdapter()

    # Aspirin SMILES
    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    result = await adapter.execute(smiles)

    assert result.success is True
    assert result.data is not None
    assert "molecular_weight" in result.data
    assert result.data["molecular_weight"] is not None
    # Aspirin MW should be around 180 g/mol
    assert 175 < result.data["molecular_weight"] < 185


@pytest.mark.asyncio
async def test_pubchem_adapter_caffeine():
    """Test PubChem adapter with caffeine SMILES"""
    adapter = PubChemAdapter()

    # Caffeine SMILES
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    result = await adapter.execute(smiles)

    assert result.success is True
    assert result.data is not None
    assert "canonical_smiles" in result.data
    assert "inchikey" in result.data


@pytest.mark.asyncio
async def test_pubchem_adapter_invalid_smiles():
    """Test PubChem adapter with invalid SMILES"""
    adapter = PubChemAdapter()

    # Invalid SMILES
    smiles = "INVALID_SMILES_$%^"

    result = await adapter.execute(smiles)

    assert result.success is False
    assert result.error is not None


def test_pubchem_adapter_validation():
    """Test SMILES validation"""
    adapter = PubChemAdapter()

    # Valid SMILES
    assert adapter.validate_input("CCO") is True
    assert adapter.validate_input("c1ccccc1") is True

    # Invalid inputs
    assert adapter.validate_input("") is False
    assert adapter.validate_input(123) is False
    assert adapter.validate_input(None) is False


def test_pubchem_adapter_metadata():
    """Test adapter metadata"""
    adapter = PubChemAdapter()

    metadata = adapter.get_metadata()

    assert metadata["name"] == "pubchem"
    assert metadata["type"] == "api"
    assert metadata["version"] == "1.0.0"
    assert metadata["enabled"] is True
