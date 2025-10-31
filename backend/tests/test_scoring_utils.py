
import pytest
from backend.core.scoring_utils import vina_affinity_to01, synthesis_steps_to01, _validate_score

def test_vina_normalization():
    assert vina_affinity_to01(-12.0) == 1.0
    assert vina_affinity_to01(-4.0) == 0.0
    assert 0 < vina_affinity_to01(-8.0) < 1

def test_synthesis_normalization():
    assert synthesis_steps_to01(1) == 1.0
    assert synthesis_steps_to01(10) == 0.0
    assert 0 < synthesis_steps_to01(5) < 1

def test_validate_score_valid():
    """Test that valid scores pass validation."""
    # Should not raise for valid scores
    _validate_score(0.0, "test_score")
    _validate_score(0.5, "test_score")
    _validate_score(1.0, "test_score")

def test_validate_score_invalid_low():
    """Test that scores below 0 raise ValueError."""
    with pytest.raises(ValueError, match="must be in \\[0, 1\\] range"):
        _validate_score(-0.1, "test_score")

def test_validate_score_invalid_high():
    """Test that scores above 1 raise ValueError."""
    with pytest.raises(ValueError, match="must be in \\[0, 1\\] range"):
        _validate_score(1.5, "test_score")

def test_validate_score_custom_name():
    """Test that custom score names appear in error messages."""
    with pytest.raises(ValueError, match="binding_score must be in"):
        _validate_score(2.0, "binding_score")
