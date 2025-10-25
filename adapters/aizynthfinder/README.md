# AiZynthFinder Adapter

**Retrosynthesis route planning using AiZynthFinder**

This adapter integrates AiZynthFinder into PharmForge to provide AI-powered retrosynthetic analysis for drug candidates.

## Overview

**AiZynthFinder** is a tool for retrosynthetic planning that uses:
- Monte Carlo Tree Search (MCTS) for route exploration
- Neural network policies trained on reaction databases
- Stock molecule filtering (ZINC, eMolecules, etc.)

**Reference:** Genheden et al., *J. Chem. Inf. Model.* 2020, 60, 12, 5910-5919

## Features

- **Multi-route discovery**: Find multiple synthesis pathways
- **Step counting**: Count synthetic steps from commercial starting materials
- **Feasibility scoring**: 0-1 normalized scores (higher = easier synthesis)
- **Reaction templates**: Extract reaction SMILES for each step
- **Starting material identification**: List commercially available precursors

## Installation

### Basic Installation

```bash
pip install aizynthfinder
```

### Full Setup (Required for Production)

AiZynthFinder requires additional setup:

1. **Install the package:**
   ```bash
   pip install aizynthfinder
   ```

2. **Download trained models:**
   ```bash
   # Follow instructions at:
   # https://molecularai.github.io/aizynthfinder/getting_started.html
   ```

3. **Configure stock databases:**
   Create `aizynthfinder_config.yml`:
   ```yaml
   properties:
     policy:
       files:
         - path/to/uspto_model.hdf5
     stock:
       files:
         - path/to/zinc_stock.hdf5
   ```

4. **Set environment variables (optional):**
   ```bash
   export AIZYNTHFINDER_CONFIG_PATH=/path/to/aizynthfinder_config.yml
   ```

## Usage

### Basic Usage

```python
from adapters.aizynthfinder.adapter import AiZynthFinderAdapter

# Initialize adapter
adapter = AiZynthFinderAdapter()

# Run retrosynthesis
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin

# Check results
if result.success:
    data = result.data
    print(f"Found {data['routes_found']} synthesis routes")
    print(f"Best route: {data['best_route']['n_steps']} steps")
    print(f"Feasibility: {data['best_route']['feasibility']}")
    print(f"Synthesis score: {data['best_route']['synthesis_score']:.3f}")
```

### Custom Configuration

```python
adapter = AiZynthFinderAdapter(config={
    "max_routes": 10,         # Find up to 10 routes
    "expansion_time": 90,     # Search for 90 seconds
    "stock": "emolecules",    # Use eMolecules stock
})

result = await adapter.execute(
    "CN1CCN(CC1)C(=O)C2=C(C=CC(=C2)OC)OC",
    max_routes=5,             # Override config
    expansion_time=60
)
```

### Integration with PharmForge Pipeline

```python
from backend.core.adapter_registry import registry

# Get registered adapter
retro_adapter = registry.get("aizynthfinder")

# Use in pipeline
result = await retro_adapter.execute(
    candidate_smiles,
    max_routes=3,
    expansion_time=60
)

# Extract synthesis metrics
if result.success and result.data["routes_found"] > 0:
    best_route = result.data["best_route"]

    # Use normalized score for ranking
    synthesis_score = best_route["synthesis_score"]  # 0-1, higher is better
    n_steps = best_route["n_steps"]
    feasibility = best_route["feasibility"]  # "high", "medium", "low"
```

## Output Format

```python
{
    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    "routes_found": 3,
    "routes": [
        {
            "route_id": 1,
            "n_steps": 2,
            "n_reactions": 2,
            "synthesis_score": 0.853,          # 0-1, normalized (higher=better)
            "feasibility": "high",             # "high", "medium", "low"
            "tree_score": 0.921,               # MCTS confidence
            "reaction_smiles": [
                "CC(=O)Cl.Oc1ccccc1C(=O)O>>CC(=O)Oc1ccccc1C(=O)O",
                "..."
            ],
            "starting_materials": [
                "CC(=O)Cl",
                "Oc1ccccc1C(=O)O"
            ],
            "is_complete": true
        },
        # ... more routes
    ],
    "best_route": { /* route with highest synthesis_score */ },
    "n_steps": 2,                              # From best route
    "synthesis_score": 0.853,                  # From best route
    "feasibility": "high",                     # From best route
    "model": "AiZynthFinder MCTS",
    "reference": "Genheden et al., J. Chem. Inf. Model. 2020"
}
```

### If No Routes Found

```python
{
    "smiles": "...",
    "routes_found": 0,
    "routes": [],
    "best_route": null,
    "n_steps": null,
    "synthesis_score": 0.0,
    "feasibility": "none",
    "model": "AiZynthFinder MCTS",
    "reference": "Genheden et al., J. Chem. Inf. Model. 2020"
}
```

## Scoring System

### Synthesis Score (0-1 scale)

The adapter uses `synthesis_steps_to01()` from `backend/core/scoring_utils.py`:

- **Fewer steps = Higher score**
- 1 step → 1.0 (best)
- 5 steps → 0.6
- 10 steps → 0.0 (worst)

### Combined Feasibility Score

The adapter combines two metrics:
```python
combined_score = 0.6 * tree_search_score + 0.4 * step_efficiency_score
```

- **tree_search_score**: MCTS confidence (from AiZynthFinder)
- **step_efficiency_score**: Normalized step count (fewer is better)

### Feasibility Classification

- **High** (≥0.7): Easy synthesis, likely feasible
- **Medium** (0.4-0.7): Moderate difficulty
- **Low** (<0.4): Difficult synthesis, may be impractical

## Configuration Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_routes` | int | 5 | Maximum number of routes to return |
| `timeout` | int | 120 | Overall timeout in seconds |
| `expansion_time` | int | 60 | MCTS expansion time in seconds |
| `stock` | str | "zinc" | Stock database ("zinc", "emolecules", "molport") |
| `return_first` | bool | false | Return immediately after first route found |

## Caching

The adapter implements deterministic caching:

```python
cache_key = SHA256(
    "aizynthfinder_v1.0.0" +
    canonical_smiles +
    max_routes +
    expansion_time
)
```

This ensures:
- Identical inputs produce identical cache keys
- Version changes invalidate cache
- Parameter changes create new cache entries

## Error Handling

### Import Errors

If AiZynthFinder is not installed:

```python
AdapterResult(
    success=False,
    error="AiZynthFinder library not installed. Install with: pip install aizynthfinder",
    metadata={"installation_help": "pip install aizynthfinder"}
)
```

### Invalid SMILES

```python
AdapterResult(
    success=False,
    error="Invalid SMILES string",
    data={}
)
```

### No Routes Found

```python
AdapterResult(
    success=True,  # Not an error!
    data={
        "routes_found": 0,
        "synthesis_score": 0.0,
        "feasibility": "none"
    },
    metadata={"message": "No synthesis routes found"}
)
```

## Testing

### Run Basic Tests

```bash
cd backend
pytest tests/test_aizynthfinder_adapter.py -v
```

### Run with AiZynthFinder Installed

```bash
pytest tests/test_aizynthfinder_adapter.py -v -m "not skipif"
```

### Run Integration Tests

```bash
pytest tests/test_aizynthfinder_adapter.py -v -m integration
```

## Performance Notes

- **Search time**: 30-120 seconds typical for MCTS expansion
- **Memory**: Requires 2-4 GB for model files
- **CPU**: Benefits from multi-core CPUs
- **Caching**: Essential for repeated queries (saves minutes per molecule)

## Troubleshooting

### "AiZynthFinder library not installed"

**Solution:**
```bash
pip install aizynthfinder
```

### "No module named 'aizynthfinder.context.config'"

**Solution:** Update to latest version:
```bash
pip install --upgrade aizynthfinder
```

### "No policy model found"

**Solution:** Download and configure models:
1. Visit: https://molecularai.github.io/aizynthfinder/getting_started.html
2. Download USPTO or other policy models
3. Create `aizynthfinder_config.yml` pointing to model files

### No routes found for simple molecules

**Possible causes:**
- Stock database not configured
- Model files missing or outdated
- Search time too short (increase `expansion_time`)

**Solution:**
```python
adapter = AiZynthFinderAdapter(config={
    "expansion_time": 120,  # Longer search
    "max_routes": 10        # More routes
})
```

## References

- **Paper:** Genheden, S. et al. "AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning." *J. Chem. Inf. Model.* 2020, 60, 12, 5910-5919
- **GitHub:** https://github.com/MolecularAI/aizynthfinder
- **Documentation:** https://molecularai.github.io/aizynthfinder/

## License

This adapter is part of PharmForge (MIT License).

AiZynthFinder itself is licensed under Apache 2.0.
