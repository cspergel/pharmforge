"""
Score normalization utilities for PharmForge

All scoring functions normalize values to [0, 1] range where higher is better.
This ensures consistent multi-objective optimization across different metrics.
"""


def to01(val: float, lo: float, hi: float) -> float:
    """
    Normalize a value to [0, 1] range using linear interpolation.

    Args:
        val: Value to normalize
        lo: Lower bound (maps to 0.0)
        hi: Upper bound (maps to 1.0)

    Returns:
        Normalized value in [0, 1] range

    Examples:
        >>> to01(5.0, 0.0, 10.0)
        0.5
        >>> to01(-8.0, -12.0, -4.0)
        0.5
    """
    if lo == hi:
        return 0.0
    v = max(min(val, hi), lo)
    return (v - lo) / (hi - lo)


def vina_affinity_to01(kcal: float) -> float:
    """
    Normalize AutoDock Vina binding affinity to [0, 1] scale.

    Vina reports binding affinity in kcal/mol where more negative = better binding.
    This function converts to 0-1 scale where higher = better.

    Maps [-12 kcal/mol, -4 kcal/mol] -> [1.0, 0.0]

    Args:
        kcal: Vina binding affinity in kcal/mol (typically -12 to -4)

    Returns:
        Normalized score where 1.0 = best binding (-12 kcal/mol),
                              0.0 = worst binding (-4 kcal/mol)

    Examples:
        >>> vina_affinity_to01(-12.0)  # Strong binder
        1.0
        >>> vina_affinity_to01(-4.0)   # Weak binder
        0.0
        >>> vina_affinity_to01(-8.0)   # Medium binder
        0.5
    """
    return to01(-kcal, 4.0, 12.0)


def synthesis_steps_to01(steps: int) -> float:
    """
    Normalize retrosynthesis steps to [0, 1] scale.

    Fewer synthesis steps = better (easier to synthesize).
    This function converts step count to 0-1 scale where higher = better.

    Maps [1 step, 10 steps] -> [1.0, 0.0]

    Args:
        steps: Number of synthesis steps (clamped to 1-10 range)

    Returns:
        Normalized score where 1.0 = 1 step (easiest to synthesize),
                              0.0 = 10 steps (hardest to synthesize)

    Examples:
        >>> synthesis_steps_to01(1)   # Very easy
        1.0
        >>> synthesis_steps_to01(10)  # Very difficult
        0.0
        >>> synthesis_steps_to01(5)   # Medium difficulty
        0.5555555555555556
    """
    steps = max(1, min(int(steps), 10))
    return to01(11 - steps, 1.0, 10.0)

def _validate_score(score: float, score_name: str = "score") -> None:
    """
    Validate that a score is in the [0, 1] range.

    Args:
        score: The score value to validate
        score_name: Name of the score for error messages

    Raises:
        ValueError: If score is not in [0, 1] range
    """
    if not (0.0 <= score <= 1.0):
        raise ValueError(f"{score_name} must be in [0, 1] range, got {score}")
