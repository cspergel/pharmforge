"""
Tests for Vina Molecular Docking Adapter

Tests the VinaAdapter implementation including:
- Input validation
- 3D structure generation
- PDBQT preparation
- Docking execution (mocked)
- Result parsing and normalization
"""

import pytest
import asyncio
import os
import tempfile
from unittest.mock import Mock, patch, AsyncMock

from adapters.vina.adapter import VinaAdapter
from backend.core.adapters.protocol import AdapterResult
from backend.core.scoring_utils import vina_affinity_to01


class TestVinaAdapter:
    """Test suite for VinaAdapter"""

    @pytest.fixture
    def adapter(self):
        """Create a VinaAdapter instance with test configuration"""
        config = {
            'receptor_path': '/path/to/receptor.pdbqt',
            'center_x': 10.0,
            'center_y': 20.0,
            'center_z': 15.0,
            'size_x': 25,
            'size_y': 25,
            'size_z': 25,
            'exhaustiveness': 8,
            'num_modes': 9
        }
        return VinaAdapter(config=config)

    def test_adapter_initialization(self, adapter):
        """Test adapter initializes with correct configuration"""
        assert adapter.name == "vina_docking"
        assert adapter.adapter_type == "local"
        assert adapter.version == "1.0.0"
        assert adapter.receptor_path == '/path/to/receptor.pdbqt'
        assert adapter.center_x == 10.0
        assert adapter.exhaustiveness == 8

    def test_validate_input_valid_smiles(self, adapter):
        """Test validation accepts valid SMILES"""
        valid_smiles = [
            "CCO",  # Ethanol
            "c1ccccc1",  # Benzene
            "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
        ]

        for smiles in valid_smiles:
            assert adapter.validate_input(smiles), f"Should accept valid SMILES: {smiles}"

    def test_validate_input_invalid(self, adapter):
        """Test validation rejects invalid input"""
        invalid_inputs = [
            "",  # Empty string
            None,  # None
            123,  # Not a string
            "INVALID_SMILES_XYZ",  # Invalid SMILES
            "C(C)(C",  # Unclosed parenthesis
        ]

        for invalid in invalid_inputs:
            assert not adapter.validate_input(invalid), f"Should reject invalid input: {invalid}"

    def test_smiles_to_3d(self, adapter):
        """Test conversion of SMILES to 3D structure"""
        smiles = "CCO"  # Ethanol

        mol_3d = adapter._smiles_to_3d(smiles)

        assert mol_3d is not None, "Should generate 3D molecule"
        assert mol_3d.GetNumAtoms() > 0, "Molecule should have atoms"
        assert mol_3d.GetConformer() is not None, "Should have 3D coordinates"

    def test_smiles_to_3d_complex(self, adapter):
        """Test 3D generation for complex molecule"""
        smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

        mol_3d = adapter._smiles_to_3d(smiles)

        assert mol_3d is not None
        assert mol_3d.GetNumAtoms() >= 21, "Aspirin should have 21+ atoms (with H)"

    def test_write_pdb(self, adapter):
        """Test writing molecule to PDB file"""
        smiles = "CCO"
        mol_3d = adapter._smiles_to_3d(smiles)

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
            pdb_path = f.name

        try:
            success = adapter._write_pdb(mol_3d, pdb_path)
            assert success, "Should write PDB successfully"
            assert os.path.exists(pdb_path), "PDB file should exist"
            assert os.path.getsize(pdb_path) > 0, "PDB file should not be empty"

            # Verify PDB content
            with open(pdb_path, 'r') as f:
                content = f.read()
                assert "ATOM" in content or "HETATM" in content, "Should contain atom records"

        finally:
            if os.path.exists(pdb_path):
                os.remove(pdb_path)

    def test_prepare_ligand_pdbqt(self, adapter):
        """Test PDBQT preparation from PDB"""
        smiles = "CCO"
        mol_3d = adapter._smiles_to_3d(smiles)

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
            pdb_path = f.name

        with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False) as f:
            pdbqt_path = f.name

        try:
            adapter._write_pdb(mol_3d, pdb_path)
            success = adapter._prepare_ligand_pdbqt(pdb_path, pdbqt_path)

            assert success, "Should prepare PDBQT successfully"
            assert os.path.exists(pdbqt_path), "PDBQT file should exist"
            assert os.path.getsize(pdbqt_path) > 0, "PDBQT file should not be empty"

            # Verify PDBQT content
            with open(pdbqt_path, 'r') as f:
                content = f.read()
                assert "ATOM" in content, "Should contain ATOM records"

        finally:
            if os.path.exists(pdb_path):
                os.remove(pdb_path)
            if os.path.exists(pdbqt_path):
                os.remove(pdbqt_path)

    def test_parse_vina_log(self, adapter):
        """Test parsing of Vina log file"""
        # Create mock Vina log content
        log_content = """
Detected 8 CPUs
Reading input ... done.
Setting up the scoring function ... done.
Analyzing the binding site ... done.
Using random seed: 0
Performing search ... done.
Refining results ... done.

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1        -7.5      0.000      0.000
   2        -7.2      1.234      2.567
   3        -6.8      2.345      3.456
"""

        with tempfile.NamedTemporaryFile(mode='w', suffix=".log", delete=False) as f:
            f.write(log_content)
            log_path = f.name

        try:
            results = adapter._parse_vina_log(log_path)

            assert len(results) == 3, "Should parse 3 modes"
            assert results[0]['mode'] == 1
            assert results[0]['affinity'] == -7.5
            assert results[0]['rmsd_lb'] == 0.0
            assert results[1]['affinity'] == -7.2
            assert results[2]['affinity'] == -6.8

        finally:
            if os.path.exists(log_path):
                os.remove(log_path)

    def test_scoring_normalization(self):
        """Test that score normalization follows expected behavior"""
        # Test cases: (affinity_kcal, expected_normalized)
        test_cases = [
            (-12.0, 1.0),  # Very strong binding -> 1.0
            (-8.0, 0.5),   # Mid-range -> 0.5
            (-4.0, 0.0),   # Weak binding -> 0.0
            (-10.0, 0.75), # Strong binding -> 0.75
        ]

        for affinity, expected in test_cases:
            normalized = vina_affinity_to01(affinity)
            assert abs(normalized - expected) < 0.01, \
                f"Affinity {affinity} should normalize to ~{expected}, got {normalized}"

    @pytest.mark.asyncio
    async def test_execute_no_receptor(self):
        """Test execution fails gracefully without receptor"""
        adapter = VinaAdapter()  # No receptor configured
        smiles = "CCO"

        result = await adapter.execute(smiles)

        assert not result.success
        assert "receptor" in result.error.lower()

    @pytest.mark.asyncio
    async def test_execute_invalid_smiles(self, adapter):
        """Test execution fails with invalid SMILES"""
        invalid_smiles = "INVALID_XYZ"

        result = await adapter.execute(invalid_smiles)

        assert not result.success
        assert "invalid" in result.error.lower()

    @pytest.mark.asyncio
    async def test_execute_success_mocked(self, adapter):
        """Test successful execution with mocked Vina"""
        smiles = "CCO"

        # Mock the Vina execution
        mock_log_content = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1        -8.5      0.000      0.000
   2        -8.2      1.234      2.567
"""

        # Create temporary receptor file
        with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False) as f:
            f.write(b"ATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00\n")
            receptor_path = f.name

        try:
            # Update adapter config with real receptor
            adapter.receptor_path = receptor_path

            # Mock _run_vina to avoid actual Vina execution
            async def mock_run_vina(receptor, ligand, output, log):
                # Create mock log file
                with open(log, 'w') as f:
                    f.write(mock_log_content)
                # Create empty output file
                with open(output, 'w') as f:
                    f.write("ENDMDL\n")
                return True, ""

            with patch.object(adapter, '_run_vina', new=mock_run_vina):
                result = await adapter.execute(smiles)

            # Verify result
            assert result.success, f"Execution should succeed: {result.error}"
            assert result.data is not None
            assert result.data['binding_affinity'] == -8.5
            assert 'binding_score' in result.data
            assert 0 <= result.data['binding_score'] <= 1
            assert result.data['num_poses'] == 2

        finally:
            if os.path.exists(receptor_path):
                os.remove(receptor_path)

    def test_get_metadata(self, adapter):
        """Test adapter metadata retrieval"""
        metadata = adapter.get_metadata()

        assert metadata['name'] == "vina_docking"
        assert metadata['type'] == "local"
        assert metadata['version'] == "1.0.0"
        assert 'description' in metadata
        assert 'requirements' in metadata
        assert 'AutoDock Vina' in str(metadata['requirements'])
        assert 'scoring' in metadata
        assert metadata['scoring']['raw_unit'] == "kcal/mol"

    def test_cache_key_generation(self, adapter):
        """Test deterministic cache key generation"""
        smiles = "CCO"

        key1 = adapter.generate_cache_key(smiles)
        key2 = adapter.generate_cache_key(smiles)

        assert key1 == key2, "Cache keys should be deterministic"
        assert len(key1) == 64, "Should be SHA256 hash (64 hex chars)"

        # Different SMILES should have different keys
        key3 = adapter.generate_cache_key("c1ccccc1")
        assert key1 != key3, "Different inputs should have different keys"


# Integration test (requires actual Vina binary - skipped in CI)
@pytest.mark.skip(reason="Requires Vina binary and receptor file")
@pytest.mark.asyncio
async def test_vina_integration():
    """Integration test with real Vina binary (manual testing only)"""
    # Configuration for real test
    config = {
        'receptor_path': '/path/to/real/receptor.pdbqt',
        'center_x': 10.0,
        'center_y': 20.0,
        'center_z': 15.0,
        'vina_binary': 'vina'  # or full path
    }

    adapter = VinaAdapter(config=config)
    smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

    result = await adapter.execute(smiles)

    assert result.success
    assert result.data['binding_affinity'] < 0  # Should be negative
    assert 0 <= result.data['binding_score'] <= 1
    print(f"Binding affinity: {result.data['binding_affinity']} kcal/mol")
    print(f"Normalized score: {result.data['binding_score']}")


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "-s"])
