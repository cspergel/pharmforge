"""
LLM-based Retrosynthesis Adapter for PharmForge

Uses Large Language Models (Claude, GPT-4) to generate synthesis routes
when template-based methods fail or as an alternative approach.

Features:
- Works for ANY molecule (no template database limitations)
- Provides human-readable explanations
- Suggests reagents, conditions, and procedures
- Can generate multiple alternative routes
- Fallback when AiZynthFinder returns 0 routes

Reference: Leveraging LLM knowledge of chemistry literature
"""

import hashlib
import logging
import json
import re
import os
from typing import Dict, Any, List, Optional
from rdkit import Chem

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class LLMRetrosynthesisAdapter(AdapterProtocol):
    """
    LLM-based retrosynthesis adapter using Claude or GPT-4 APIs.

    This adapter uses the knowledge embedded in large language models
    to suggest synthesis routes for any molecule, without the limitations
    of template-based approaches.

    Advantages over template-based:
    - No template database limitations
    - Works for novel/uncommon molecules
    - Provides human-readable explanations
    - Can suggest alternative synthetic strategies

    Disadvantages:
    - Less precise than computational methods
    - May hallucinate or provide infeasible routes
    - Requires API access (costs money)
    """

    def __init__(self, name: str = "llm_retrosynthesis", adapter_type: str = "api", config: Optional[Dict[str, Any]] = None):
        """
        Initialize LLM retrosynthesis adapter.

        Args:
            name: Adapter name
            adapter_type: Type (default: "api")
            config: Configuration dictionary with keys:
                   - provider: "claude" or "openai" (default: "claude")
                   - model: Model name (default: "claude-3-5-sonnet-20241022")
                   - api_key: API key (optional, reads from env)
                   - num_routes: Number of routes to generate (default: 3)
                   - temperature: LLM temperature (default: 0.7)
        """
        default_config = {
            "provider": "claude",  # or "openai"
            "model": "claude-3-5-sonnet-20241022",  # or "gpt-4o"
            "api_key": None,  # Read from environment if not provided
            "num_routes": 3,
            "temperature": 0.7,
            "max_tokens": 4000,
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Get API key from config or environment
        self.api_key = self.config.get("api_key") or self._get_api_key_from_env()

    def _get_api_key_from_env(self) -> Optional[str]:
        """Get API key from environment variables."""
        provider = self.config.get("provider", "claude")

        if provider == "claude":
            return os.getenv("ANTHROPIC_API_KEY") or os.getenv("CLAUDE_API_KEY")
        elif provider == "openai":
            return os.getenv("OPENAI_API_KEY")

        return None

    def validate_input(self, smiles: str) -> bool:
        """Validate SMILES string input."""
        if not smiles or not isinstance(smiles, str):
            return False

        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.warning(f"SMILES validation error: {e}")
            return False

    def _build_retrosynthesis_prompt(self, smiles: str, num_routes: int = 3) -> str:
        """
        Build prompt for LLM to generate retrosynthesis routes.

        Args:
            smiles: Target molecule SMILES
            num_routes: Number of routes to generate

        Returns:
            Prompt string
        """
        # Get molecule name if possible
        try:
            mol = Chem.MolFromSmiles(smiles)
            mol_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
        except:
            mol_formula = "Unknown"
            mol_weight = 0

        prompt = f"""You are an expert synthetic organic chemist. Please suggest {num_routes} feasible retrosynthetic routes to synthesize the following target molecule.

Target Molecule:
- SMILES: {smiles}
- Molecular Formula: {mol_formula}
- Molecular Weight: {mol_weight:.2f} g/mol

For each route, provide:
1. A brief overview of the synthetic strategy
2. Step-by-step synthesis with reagents and conditions
3. Starting materials needed (commercially available or simple compounds)
4. Estimated number of steps
5. Key challenges or considerations
6. Approximate feasibility rating (1-10, where 10 is most feasible)

Format your response as a structured JSON object with this schema:
{{
  "routes": [
    {{
      "route_id": 1,
      "strategy": "Brief description of the synthetic strategy",
      "steps": [
        {{
          "step_number": 1,
          "description": "Reaction description",
          "reagents": ["reagent1", "reagent2"],
          "conditions": "Temperature, solvent, etc.",
          "product_smiles": "SMILES of the intermediate or final product"
        }}
      ],
      "starting_materials": ["compound1", "compound2"],
      "n_steps": 3,
      "feasibility_score": 8,
      "challenges": "Key challenges for this route",
      "advantages": "Why this route might be preferred"
    }}
  ]
}}

Provide practical, chemically sound routes that could realistically be performed in a well-equipped synthetic laboratory. Focus on established reactions and commercially available starting materials.
"""
        return prompt

    async def execute(self, smiles: str, **params) -> AdapterResult:
        """
        Execute LLM-based retrosynthesis route finding.

        Args:
            smiles: SMILES string of the target molecule
            **params: Additional parameters:
                     - num_routes: Override default num_routes
                     - temperature: Override default temperature

        Returns:
            AdapterResult containing synthesis routes and metrics
        """
        try:
            # Validate input
            if not self.validate_input(smiles):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid SMILES string"
                )

            # Check API key
            if not self.api_key:
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"No API key configured for {self.config.get('provider', 'claude')}"
                )

            # Get parameters
            num_routes = params.get('num_routes', self.config.get('num_routes', 3))
            temperature = params.get('temperature', self.config.get('temperature', 0.7))

            # Generate cache key
            cache_key = self.generate_cache_key(smiles, num_routes=num_routes)

            logger.info(f"Running LLM retrosynthesis for: {smiles[:50]}...")
            logger.info(f"Provider: {self.config.get('provider')}, Model: {self.config.get('model')}")

            # Build prompt
            prompt = self._build_retrosynthesis_prompt(smiles, num_routes)

            # Call appropriate LLM API
            provider = self.config.get("provider", "claude")
            if provider == "claude":
                response_text = await self._call_claude_api(prompt, temperature)
            elif provider == "openai":
                response_text = await self._call_openai_api(prompt, temperature)
            else:
                return AdapterResult(
                    success=False,
                    data={},
                    error=f"Unsupported provider: {provider}"
                )

            # Parse LLM response
            routes_data = self._parse_llm_response(response_text)

            if not routes_data or "routes" not in routes_data:
                logger.warning(f"No valid routes parsed from LLM response for {smiles}")
                return AdapterResult(
                    success=True,
                    data={
                        "smiles": smiles,
                        "routes_found": 0,
                        "routes": [],
                        "best_route": None,
                        "synthesis_score": 0.0,
                        "n_steps": None,
                        "feasibility": "none"
                    },
                    metadata={
                        "adapter_name": self.name,
                        "cache_key": cache_key,
                        "version": self.version,
                        "provider": provider,
                        "message": "No synthesis routes found"
                    }
                )

            # Process routes
            route_list = []
            for route in routes_data["routes"]:
                # Normalize feasibility score to 0-1 range
                feasibility_score = route.get("feasibility_score", 5)
                synthesis_score = feasibility_score / 10.0

                # Classify feasibility
                if feasibility_score > 7:
                    feasibility = "high"
                elif feasibility_score >= 4:
                    feasibility = "medium"
                else:
                    feasibility = "low"

                route_list.append({
                    "route_id": route.get("route_id", len(route_list) + 1),
                    "n_steps": route.get("n_steps", len(route.get("steps", []))),
                    "synthesis_score": round(synthesis_score, 3),
                    "feasibility": feasibility,
                    "feasibility_score": feasibility_score,
                    "strategy": route.get("strategy", ""),
                    "starting_materials": route.get("starting_materials", []),
                    "steps": route.get("steps", []),
                    "challenges": route.get("challenges", ""),
                    "advantages": route.get("advantages", "")
                })

            # Find best route (highest synthesis score)
            best_route = max(route_list, key=lambda x: x["synthesis_score"])

            # Extract unique starting materials
            all_starting_materials = set()
            for route in route_list:
                all_starting_materials.update(route["starting_materials"])

            # Create result data structure
            result_data = {
                "smiles": smiles,
                "routes_found": len(route_list),
                "routes": route_list,
                "best_route": best_route,
                "n_steps": best_route["n_steps"],
                "synthesis_score": best_route["synthesis_score"],
                "feasibility": best_route["feasibility"],
                "model": f"{provider}/{self.config.get('model')}",
                "reference": "LLM-based retrosynthesis",
                "requirements": {
                    "starting_materials": best_route["starting_materials"],
                    "n_steps": best_route["n_steps"],
                    "all_possible_starting_materials": list(all_starting_materials)
                }
            }

            logger.info(f"✓ Found {len(route_list)} synthesis routes from LLM")
            logger.info(f"  Best route: {best_route['n_steps']} steps, score={best_route['synthesis_score']:.3f}")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "cache_key": cache_key,
                    "version": self.version,
                    "provider": provider,
                    "model": self.config.get('model')
                }
            )

        except Exception as e:
            logger.error(f"LLM retrosynthesis failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name,
                    "provider": self.config.get("provider")
                }
            )

    async def _call_claude_api(self, prompt: str, temperature: float = 0.7) -> str:
        """
        Call Claude API to generate synthesis routes.

        Args:
            prompt: Retrosynthesis prompt
            temperature: Sampling temperature

        Returns:
            LLM response text
        """
        try:
            import anthropic
        except ImportError:
            raise ImportError(
                "anthropic package not installed. Install with: pip install anthropic"
            )

        try:
            client = anthropic.Anthropic(api_key=self.api_key)

            message = client.messages.create(
                model=self.config.get("model", "claude-3-5-sonnet-20241022"),
                max_tokens=self.config.get("max_tokens", 4000),
                temperature=temperature,
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )

            # Extract text from response
            response_text = message.content[0].text
            logger.debug(f"Claude API response length: {len(response_text)} chars")
            return response_text

        except anthropic.APIError as e:
            logger.error(f"Claude API error: {e}")
            raise Exception(f"Claude API error: {str(e)}")
        except Exception as e:
            logger.error(f"Error calling Claude API: {e}")
            raise

    async def _call_openai_api(self, prompt: str, temperature: float = 0.7) -> str:
        """
        Call OpenAI API to generate synthesis routes.

        Args:
            prompt: Retrosynthesis prompt
            temperature: Sampling temperature

        Returns:
            LLM response text
        """
        try:
            import openai
        except ImportError:
            raise ImportError(
                "openai package not installed. Install with: pip install openai"
            )

        try:
            client = openai.OpenAI(api_key=self.api_key)

            response = client.chat.completions.create(
                model=self.config.get("model", "gpt-4o"),
                messages=[
                    {
                        "role": "system",
                        "content": "You are an expert synthetic organic chemist."
                    },
                    {
                        "role": "user",
                        "content": prompt
                    }
                ],
                temperature=temperature,
                max_tokens=self.config.get("max_tokens", 4000)
            )

            response_text = response.choices[0].message.content
            logger.debug(f"OpenAI API response length: {len(response_text)} chars")
            return response_text

        except openai.APIError as e:
            logger.error(f"OpenAI API error: {e}")
            raise Exception(f"OpenAI API error: {str(e)}")
        except Exception as e:
            logger.error(f"Error calling OpenAI API: {e}")
            raise

    def _parse_llm_response(self, response_text: str) -> Dict[str, Any]:
        """
        Parse JSON response from LLM.

        Args:
            response_text: Raw LLM response

        Returns:
            Parsed routes dictionary
        """
        try:
            # Try to extract JSON from markdown code blocks
            json_match = re.search(r'```json\s*(.*?)\s*```', response_text, re.DOTALL)
            if json_match:
                json_text = json_match.group(1)
            else:
                # Try to find JSON object directly
                json_match = re.search(r'\{.*\}', response_text, re.DOTALL)
                if json_match:
                    json_text = json_match.group(0)
                else:
                    json_text = response_text

            # Parse JSON
            routes_data = json.loads(json_text)
            logger.debug(f"Successfully parsed {len(routes_data.get('routes', []))} routes from LLM")
            return routes_data

        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM response as JSON: {e}")
            logger.debug(f"Response text: {response_text[:500]}...")
            return {"routes": []}
        except Exception as e:
            logger.error(f"Error parsing LLM response: {e}")
            return {"routes": []}

    def generate_cache_key(self, smiles: str, **kwargs) -> str:
        """
        Generate cache key for LLM retrosynthesis.

        Args:
            smiles: Target molecule SMILES
            **kwargs: Additional parameters

        Returns:
            Cache key hash
        """
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "smiles": smiles,
            "provider": self.config.get("provider"),
            "model": self.config.get("model"),
            "params": kwargs
        }

        cache_string = json.dumps(cache_dict, sort_keys=True)
        cache_key = hashlib.sha256(cache_string.encode()).hexdigest()
        return cache_key

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata.

        Returns:
            Metadata dictionary
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "provider": self.config.get("provider", "claude"),
            "model": self.config.get("model"),
            "description": "LLM-based retrosynthesis using Claude or GPT-4",
            "capabilities": [
                "retrosynthesis",
                "route_generation",
                "synthesis_planning"
            ]
        }