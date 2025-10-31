# LLM-Based Retrosynthesis Adapter

## Overview

The LLM Retrosynthesis Adapter uses Large Language Models (Claude or GPT-4) to generate synthesis routes for any molecule. This adapter provides an alternative to template-based methods like AiZynthFinder, offering flexibility for novel or uncommon molecules.

## Features

- **Universal Coverage**: Works for ANY molecule without template database limitations
- **Human-Readable Explanations**: Provides natural language descriptions of synthetic strategies
- **Detailed Routes**: Suggests reagents, conditions, and step-by-step procedures
- **Multiple Routes**: Generates alternative synthesis pathways
- **Dual Provider Support**: Compatible with both Anthropic Claude and OpenAI GPT models
- **Fallback Option**: Can be used when template-based methods return no routes

## Advantages vs Template-Based Methods

| Feature | LLM-Based | Template-Based (AiZynthFinder) |
|---------|-----------|-------------------------------|
| Coverage | Any molecule | Limited to template database |
| Novel molecules | ✓ Excellent | ✗ Poor |
| Explanations | ✓ Human-readable | ✗ Technical only |
| Starting materials | ✓ Flexible suggestions | ✓ Stock-validated |
| Precision | ⚠ May hallucinate | ✓ Computationally precise |
| Cost | $ API costs | Free (local compute) |

## Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Required packages:
# - rdkit >= 2023.9.5
# - anthropic >= 0.39.0  (for Claude)
# - openai >= 1.50.0     (for GPT)
```

## API Keys

Set your API key as an environment variable:

```bash
# For Claude
export ANTHROPIC_API_KEY=your-anthropic-key

# For OpenAI
export OPENAI_API_KEY=your-openai-key
```

## Usage

### Basic Usage

```python
import asyncio
from adapters.llm_retrosynthesis import LLMRetrosynthesisAdapter

async def main():
    # Initialize adapter (uses Claude by default)
    adapter = LLMRetrosynthesisAdapter()

    # Run retrosynthesis for aspirin
    result = await adapter.execute("CC(=O)Oc1ccccc1C(=O)O")

    if result.success:
        print(f"Routes found: {result.data['routes_found']}")
        print(f"Best route: {result.data['n_steps']} steps")
        print(f"Score: {result.data['synthesis_score']}")
    else:
        print(f"Error: {result.error}")

asyncio.run(main())
```

### Custom Configuration

```python
# Use OpenAI GPT-4
adapter = LLMRetrosynthesisAdapter(config={
    "provider": "openai",
    "model": "gpt-4o",
    "num_routes": 5,
    "temperature": 0.7,
    "api_key": "your-key"  # Optional, reads from env by default
})

# Generate routes
result = await adapter.execute(
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
    num_routes=3  # Override default
)
```

### Providers and Models

**Claude (Anthropic):**
- `claude-3-5-sonnet-20241022` (default, recommended)
- `claude-3-opus-20240229`
- `claude-3-sonnet-20240229`

**OpenAI:**
- `gpt-4o` (recommended)
- `gpt-4-turbo`
- `gpt-4`

## Response Format

The adapter returns an `AdapterResult` with the following structure:

```python
{
    "success": True,
    "data": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "routes_found": 3,
        "routes": [
            {
                "route_id": 1,
                "n_steps": 1,
                "synthesis_score": 0.9,  # Normalized 0-1
                "feasibility": "high",   # high/medium/low
                "feasibility_score": 9,  # Original 1-10
                "strategy": "Esterification of salicylic acid...",
                "starting_materials": ["salicylic acid", "acetic anhydride"],
                "steps": [
                    {
                        "step_number": 1,
                        "description": "Acetylation of salicylic acid",
                        "reagents": ["acetic anhydride", "sulfuric acid"],
                        "conditions": "Heat at 80-90°C for 30 minutes",
                        "product_smiles": "CC(=O)Oc1ccccc1C(=O)O"
                    }
                ],
                "challenges": "Requires careful temperature control...",
                "advantages": "Classical synthesis, high yielding..."
            },
            # ... more routes
        ],
        "best_route": {  # Highest scoring route
            "route_id": 1,
            "n_steps": 1,
            "synthesis_score": 0.9,
            # ... full route details
        },
        "requirements": {
            "starting_materials": ["salicylic acid", "acetic anhydride"],
            "n_steps": 1,
            "all_possible_starting_materials": ["salicylic acid", "acetic anhydride", ...]
        }
    },
    "metadata": {
        "adapter_name": "llm_retrosynthesis",
        "provider": "claude",
        "model": "claude-3-5-sonnet-20241022",
        "version": "1.0.0"
    }
}
```

## Feasibility Scoring

Routes are scored based on the LLM's feasibility assessment:

- **1-10 Score**: LLM rates each route from 1 (infeasible) to 10 (highly feasible)
- **Normalized Score**: Converted to 0.0-1.0 range for consistency with other adapters
- **Classification**:
  - `high`: Score > 7 (feasible in most labs)
  - `medium`: Score 4-7 (requires specialized equipment/expertise)
  - `low`: Score < 4 (very difficult or impractical)

## Testing

### Comprehensive Test Suite (Recommended)

The main test script tests the adapter with multiple molecules, providers, and scenarios:

```bash
# Run all tests in mock mode (no API calls, no cost)
python test_llm_retrosynthesis.py --mock

# Test with Claude API (requires ANTHROPIC_API_KEY)
python test_llm_retrosynthesis.py --provider claude

# Test with OpenAI API (requires OPENAI_API_KEY)
python test_llm_retrosynthesis.py --provider openai

# Test both providers with all molecules
python test_llm_retrosynthesis.py --provider both

# Test specific molecules only
python test_llm_retrosynthesis.py --molecules aspirin,caffeine

# Verbose output
python test_llm_retrosynthesis.py --verbose
```

The comprehensive test suite includes:
- Adapter initialization testing
- Input validation (valid/invalid SMILES)
- Output format validation
- Mock API testing (no costs)
- Real API testing (if keys available)
- Error handling scenarios
- Provider comparison
- Multiple molecule testing (aspirin, ibuprofen, caffeine, paracetamol, penicillin)

### Unit Tests (No API Required)

```bash
# Standalone tests - verify core logic
python adapters/llm_retrosynthesis/test_standalone.py
```

### Mock Tests (No API Required)

```bash
# Full adapter tests with mock responses
python adapters/llm_retrosynthesis/test_adapter.py
```

### Real API Tests

```bash
# Test with Claude
export ANTHROPIC_API_KEY=your-key
python adapters/llm_retrosynthesis/test_adapter_real.py --provider claude

# Test with OpenAI
export OPENAI_API_KEY=your-key
python adapters/llm_retrosynthesis/test_adapter_real.py --provider openai

# Test both providers and molecules
python adapters/llm_retrosynthesis/test_adapter_real.py --provider both --molecule both
```

## Integration with PharmForge

The adapter follows the standard PharmForge `AdapterProtocol`:

```python
from backend.core.adapter_registry import registry
from adapters.llm_retrosynthesis import LLMRetrosynthesisAdapter

# Register the adapter
adapter = LLMRetrosynthesisAdapter()
registry.register(adapter)

# Use in pipeline
result = await adapter(
    "CC(=O)Oc1ccccc1C(=O)O",
    use_cache=True  # Enable caching
)
```

### Caching

The adapter supports automatic caching:

- Cache keys include: SMILES, provider, model, and parameters
- Different providers/models have separate cache entries
- Cache improves response time and reduces API costs

### As a Fallback

Use LLM retrosynthesis when template-based methods fail:

```python
# Try AiZynthFinder first
aizsynth_result = await aizynthfinder.execute(smiles)

if aizsynth_result.data['routes_found'] == 0:
    # Fall back to LLM
    llm_result = await llm_adapter.execute(smiles)
```

### Integration Example

For a complete integration example that demonstrates the fallback pattern, see:
- `examples/retrosynthesis_integration.py` - Full integration with AiZynthFinder fallback logic

This example shows:
- How to try AiZynthFinder first (fast, precise)
- When to fall back to LLM (no routes or low confidence)
- How to combine results from both methods
- Quality thresholds and cost optimization

## Error Handling

The adapter handles various error conditions gracefully:

```python
# Missing API key
result = await adapter.execute(smiles)
# Returns: success=False, error="No API key configured"

# Invalid SMILES
result = await adapter.execute("invalid")
# Returns: success=False, error="Invalid SMILES string"

# API errors
# Returns: success=False, error="Claude API error: ..."

# Malformed response
# Returns: success=True, routes_found=0 (empty routes list)
```

## Implementation Details

### Files

```
adapters/llm_retrosynthesis/
├── __init__.py              # Package exports
├── adapter.py               # Main adapter implementation
├── requirements.txt         # Python dependencies
├── README.md               # This file
├── test_standalone.py      # Unit tests (no deps)
├── test_adapter.py         # Mock tests
└── test_adapter_real.py    # Real API tests
```

### Methods

- `__init__()`: Initialize adapter with config
- `validate_input()`: Validate SMILES strings
- `execute()`: Main execution method (async)
- `_build_retrosynthesis_prompt()`: Generate LLM prompt
- `_call_claude_api()`: Call Anthropic Claude API
- `_call_openai_api()`: Call OpenAI GPT API
- `_parse_llm_response()`: Parse JSON from LLM response
- `generate_cache_key()`: Generate deterministic cache keys
- `get_metadata()`: Return adapter metadata

### Prompt Engineering

The adapter uses a carefully crafted prompt that:
1. Provides molecule information (SMILES, formula, weight)
2. Requests a specific number of routes
3. Defines required output format (JSON schema)
4. Emphasizes practical, realistic routes
5. Requests feasibility scores and explanations

## Limitations

1. **May Hallucinate**: LLMs can generate plausible but incorrect routes
2. **Less Precise**: Not as algorithmically rigorous as template-based methods
3. **API Costs**: Each request costs money (Claude/OpenAI pricing)
4. **Variable Quality**: Route quality depends on LLM training data
5. **No Stock Validation**: Doesn't check commercial availability automatically

## Cost Considerations

### API Pricing (as of 2025)

**Claude (Anthropic):**
- Claude 3.5 Sonnet: ~$3 per million input tokens, ~$15 per million output tokens
- Typical request: ~1,500 input tokens, ~3,000 output tokens
- Estimated cost per molecule: $0.005 - $0.05 (0.5 to 5 cents)

**OpenAI:**
- GPT-4o: ~$2.50 per million input tokens, ~$10 per million output tokens
- Similar token usage as Claude
- Estimated cost per molecule: $0.004 - $0.04

### Cost Optimization Strategies

1. **Use as Fallback Only**: Try AiZynthFinder first (free) before LLM
2. **Enable Caching**: Cache results for repeated queries
3. **Batch Processing**: Process multiple molecules in a single session
4. **Quality Thresholds**: Set thresholds to minimize unnecessary LLM calls
5. **Monitor Usage**: Track API costs and set budgets/alerts

### Example Cost Calculation

Processing 1,000 molecules:
- If 80% have template routes (use AiZynthFinder): 800 molecules = $0
- If 20% need LLM fallback: 200 molecules × $0.02 = $4.00
- Total cost: ~$4 for 1,000 molecules

## Production Deployment Recommendations

### 1. Infrastructure Setup

```python
# Use environment variables for API keys (never hardcode)
export ANTHROPIC_API_KEY=your-production-key
export OPENAI_API_KEY=your-production-key

# Set up Redis for distributed caching
export REDIS_URL=redis://localhost:6379
```

### 2. Rate Limiting

```python
from adapters.llm_retrosynthesis import LLMRetrosynthesisAdapter
import asyncio
from asyncio import Semaphore

# Limit concurrent API calls
semaphore = Semaphore(5)  # Max 5 concurrent requests

async def rate_limited_execute(adapter, smiles):
    async with semaphore:
        return await adapter.execute(smiles)
```

### 3. Error Handling and Retries

```python
from tenacity import retry, stop_after_attempt, wait_exponential

@retry(stop=stop_after_attempt(3), wait=wait_exponential(min=1, max=10))
async def execute_with_retry(adapter, smiles):
    return await adapter.execute(smiles)
```

### 4. Monitoring and Logging

```python
import logging
from datetime import datetime

# Track API usage
api_usage_log = {
    "timestamp": datetime.now(),
    "molecule": smiles,
    "provider": "claude",
    "cost_estimate": 0.02,
    "success": True,
    "routes_found": 3
}

# Log to monitoring system
logger.info("LLM API call", extra=api_usage_log)
```

### 5. Quality Validation

```python
def validate_llm_routes(routes):
    """Validate LLM-generated routes for basic chemical sanity."""
    validated = []

    for route in routes:
        # Check for valid SMILES in steps
        if all(is_valid_smiles(step.get("product_smiles", ""))
               for step in route.get("steps", [])):
            validated.append(route)
        else:
            logger.warning(f"Route {route['route_id']} failed validation")

    return validated
```

## Best Practices

1. **Use as Complement**: Combine with template-based methods for best results
2. **Validate Routes**: Review suggested routes with chemical expertise
3. **Cache Results**: Enable caching to reduce costs and improve speed
4. **Choose Appropriate Model**: Use Claude 3.5 Sonnet or GPT-4o for best results
5. **Temperature Setting**: Use 0.7 for balance between creativity and reliability
6. **Monitor Costs**: Set up cost tracking and alerts for API usage
7. **Implement Retries**: Handle transient API failures gracefully
8. **Rate Limit**: Prevent overwhelming the API with too many concurrent requests
9. **Validate Output**: Always validate LLM responses before using in production
10. **Fallback Strategy**: Have a clear fallback pattern (AiZynthFinder → LLM → Manual)

## Version History

- **1.0.0** (2025-10-25): Initial implementation
  - Claude and OpenAI support
  - Complete route generation
  - Caching and error handling
  - Comprehensive testing suite

## License

Part of the PharmForge project.

## Support

For issues or questions:
1. Check test scripts for examples
2. Review error messages in logs
3. Verify API keys are set correctly
4. Ensure dependencies are installed

## References

- Anthropic Claude API: https://docs.anthropic.com/
- OpenAI API: https://platform.openai.com/docs/
- PharmForge Documentation: ../README.md
