
from backend.core.scoring_utils import vina_affinity_to01, synthesis_steps_to01

def test_vina_normalization():
    assert vina_affinity_to01(-12.0) == 1.0
    assert vina_affinity_to01(-4.0) == 0.0
    assert 0 < vina_affinity_to01(-8.0) < 1

def test_synthesis_normalization():
    assert synthesis_steps_to01(1) == 1.0
    assert synthesis_steps_to01(10) == 0.0
    assert 0 < synthesis_steps_to01(5) < 1
