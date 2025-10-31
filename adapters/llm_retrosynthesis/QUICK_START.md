# LLM Retrosynthesis Adapter - Quick Start Guide

## 30-Second Setup

```bash
# 1. Install dependencies
pip install rdkit anthropic openai

# 2. Set API key (choose one)
export ANTHROPIC_API_KEY=your-key  # For Claude
export OPENAI_API_KEY=your-key     # For OpenAI

# 3. Run example
cd adapters/llm_retrosynthesis
python example_usage.py
```

## Minimal Example

```python
import asyncio
from adapters.llm_retrosynthesis import LLMRetrosynthesisAdapter

async def main():
    # Initialize
    adapter = LLMRetrosynthesisAdapter()

    # Run retrosynthesis for aspirin
    result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

    # Print results
    if result.success:
        print(f"Routes: {result.data['routes_found']}")
        print(f"Best route: {result.data['n_steps']} steps")
        print(f"Materials: {result.data['requirements']['starting_materials']}")
    else:
        print(f"Error: {result.error}")

asyncio.run(main())
```

## Common Use Cases

### 1. Basic Usage (Claude)
```python
adapter = LLMRetrosynthesisAdapter()
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")
```

### 2. Use OpenAI Instead
```python
adapter = LLMRetrosynthesisAdapter(config={
    "provider": "openai",
    "model": "gpt-4o"
})
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")
```

### 3. Request More Routes
```python
adapter = LLMRetrosynthesisAdapter(config={"num_routes": 5})
result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")
```

### 4. With Caching
```python
# First call - hits API
result1 = await adapter("CC(=O)Oc1ccccc1C(=O)O", use_cache=True)

# Second call - instant from cache
result2 = await adapter("CC(=O)Oc1ccccc1C(=O)O", use_cache=True)
```

### 5. Fallback Pattern
```python
# Try AiZynthFinder first
aizsynth = await aizynthfinder.execute(smiles)

if aizsynth.data['routes_found'] == 0:
    # Fallback to LLM
    llm_result = await llm_adapter.execute(smiles)
```

## Quick Reference

### Configuration Options
```python
config = {
    "provider": "claude",              # or "openai"
    "model": "claude-3-5-sonnet-...",  # or "gpt-4o"
    "num_routes": 3,                   # number of routes
    "temperature": 0.7,                # LLM temperature
    "max_tokens": 4000,                # max response length
    "api_key": "your-key"              # optional
}
```

### Response Structure
```python
result.data = {
    "routes_found": 3,
    "n_steps": 1,
    "synthesis_score": 0.9,
    "feasibility": "high",
    "routes": [...],
    "best_route": {...},
    "requirements": {
        "starting_materials": [...],
        "n_steps": 1
    }
}
```

### Testing Commands

```bash
# Unit tests (no API)
python test_standalone.py

# Mock tests (no API)
python test_adapter.py

# Real API tests
python test_adapter_real.py --provider claude

# Examples
python example_usage.py
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| No API key | `export ANTHROPIC_API_KEY=your-key` |
| Import error | `pip install rdkit anthropic openai` |
| Invalid SMILES | Check SMILES validity with RDKit |
| No routes | Check logs, try different temperature |
| Slow | Enable caching: `use_cache=True` |

## File Locations

```
adapters/llm_retrosynthesis/
├── adapter.py           # Main implementation
├── example_usage.py     # 7 usage examples
├── test_standalone.py   # Unit tests
├── test_adapter_real.py # API tests
└── README.md           # Full documentation
```

## Key Features

✅ Works with ANY molecule
✅ Claude or OpenAI support
✅ Human-readable explanations
✅ Multiple alternative routes
✅ Automatic caching
✅ Comprehensive testing
✅ Full documentation

## Next Steps

1. **Install**: `pip install -r requirements.txt`
2. **Test**: `python test_standalone.py`
3. **Try**: `python example_usage.py`
4. **Read**: See `README.md` for full docs

---

**Version**: 1.0.0
**Status**: Production Ready ✅
**Location**: `C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2\adapters\llm_retrosynthesis\`
