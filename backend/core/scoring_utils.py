
def to01(val: float, lo: float, hi: float) -> float:
    if lo == hi: 
        return 0.0
    v = max(min(val, hi), lo)
    return (v - lo) / (hi - lo)

def vina_affinity_to01(kcal: float) -> float:
    return to01(-kcal, 4.0, 12.0)

def synthesis_steps_to01(steps: int) -> float:
    steps = max(1, min(int(steps), 10))
    return to01(11 - steps, 1.0, 10.0)
