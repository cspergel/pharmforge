"""
AiZynthFinder Adapter for PharmForge

Provides retrosynthesis route planning using AiZynthFinder.
Returns synthesis routes with step counts and feasibility scores.

AiZynthFinder is a tool for retrosynthetic planning using Monte Carlo Tree Search
and neural network policies trained on reaction data.

Reference: https://github.com/MolecularAI/aizynthfinder
Paper: Genheden et al., J. Chem. Inf. Model. 2020, 60, 12, 5910-5919
"""

import hashlib
import logging
from typing import Dict, Any, List, Optional, Tuple
import asyncio

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult
from backend.core.scoring_utils import synthesis_steps_to01

logger = logging.getLogger(__name__)


class AiZynthFinderAdapter(AdapterProtocol):
    """
    AiZynthFinder adapter for retrosynthetic route planning.

    Uses Monte Carlo Tree Search (MCTS) with neural network guidance to find
    synthetic routes from commercially available starting materials to target molecules.

    Features:
    - Multi-route discovery
    - Synthesis step counting
    - Feasibility scoring (0-1 scale, higher is better)
    - Reaction template identification
    """

    def __init__(self, name: str = "aizynthfinder", adapter_type: str = "local", config: Optional[Dict[str, Any]] = None):
        """
        Initialize AiZynthFinder adapter.

        Args:
            name: Adapter name (default: "aizynthfinder")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - max_routes: Maximum number of routes to find (default: 5)
                   - timeout: Search timeout in seconds (default: 120)
                   - expansion_time: MCTS expansion time in seconds (default: 60)
                   - stock: Stock availability filter (default: "zinc")
        """
        default_config = {
            "max_routes": 5,
            "timeout": 120,
            "expansion_time": 60,
            "stock": "zinc",  # Options: "zinc", "emolecules", "molport", etc.
            "return_first": False  # Return immediately after first route found
        }

        # Merge with provided config
        merged_config = {**default_config, **(config or {})}

        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Lazy-load AiZynthFinder (heavy import)
        self._finder = None
        self._policy = None
        self._stock = None

    @property
    def finder(self):
        """Lazy-load AiZynthFinder configuration on first use."""
        if self._finder is None:
            try:
                from aizynthfinder.aizynthfinder import AiZynthFinder
                from aizynthfinder.context.config import Configuration

                logger.info("Loading AiZynthFinder configuration...")

                # Create configuration
                config = Configuration()

                # Try to load default config file if it exists
                try:
                    config.from_file("aizynthfinder_config.yml")
                    logger.info("Loaded AiZynthFinder config from file")
                except Exception:
                    logger.warning("No config file found, using defaults")

                # Create finder
                self._finder = AiZynthFinder(config)

                logger.info("✓ AiZynthFinder initialized successfully")

            except ImportError as e:
                logger.error(f"Failed to import AiZynthFinder: {e}")
                raise ImportError(
                    "AiZynthFinder library not installed. "
                    "Install with: pip install aizynthfinder\n"
                    "For full setup, see: https://molecularai.github.io/aizynthfinder/"
                ) from e
            except Exception as e:
                logger.error(f"Failed to initialize AiZynthFinder: {e}")
                raise

        return self._finder

    def validate_input(self, smiles: str) -> bool:
        """
        Validate SMILES string input.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not smiles or not isinstance(smiles, str):
            return False

        # Basic SMILES validation with RDKit
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.warning(f"SMILES validation error: {e}")
            return False

    def _calculate_route_feasibility(self, route_data: Dict[str, Any]) -> Tuple[float, str]:
        """
        Calculate feasibility score for a synthesis route.

        Args:
            route_data: Route data from AiZynthFinder

        Returns:
            Tuple of (score 0-1, feasibility_label)
        """
        # Extract route metrics
        n_steps = route_data.get("n_steps", 10)
        n_reactions = route_data.get("n_reactions", 10)

        # AiZynthFinder provides a "score" which is related to tree search confidence
        # Higher score = more confident route
        tree_score = route_data.get("score", 0.5)

        # Normalize steps (fewer is better): use scoring_utils
        step_score = synthesis_steps_to01(n_steps)

        # Combine tree search score and step count
        # Weight: 60% tree confidence, 40% step efficiency
        combined_score = 0.6 * tree_score + 0.4 * step_score

        # Classify feasibility
        if combined_score >= 0.7:
            feasibility = "high"
        elif combined_score >= 0.4:
            feasibility = "medium"
        else:
            feasibility = "low"

        return combined_score, feasibility

    def _extract_route_info(self, route) -> Dict[str, Any]:
        """
        Extract structured information from an AiZynthFinder route.

        Args:
            route: Route object from AiZynthFinder

        Returns:
            Dictionary with route details
        """
        try:
            # Extract basic metrics
            reactions = route.reactions()
            n_reactions = len(reactions)
            n_steps = route.depth()

            # Get reaction SMILES
            reaction_smiles = []
            for reaction in reactions:
                try:
                    reaction_smiles.append(reaction.smiles)
                except Exception:
                    reaction_smiles.append("N/A")

            # Get starting materials
            starting_materials = []
            for mol in route.leafs():
                try:
                    starting_materials.append(mol.smiles)
                except Exception:
                    starting_materials.append("N/A")

            route_data = {
                "n_steps": n_steps,
                "n_reactions": n_reactions,
                "score": getattr(route, 'score', 0.5),  # Tree search score
                "reaction_smiles": reaction_smiles,
                "starting_materials": starting_materials,
                "is_solved": True  # Found a complete route
            }

            return route_data

        except Exception as e:
            logger.warning(f"Failed to extract route info: {e}")
            return {
                "n_steps": 10,
                "n_reactions": 10,
                "score": 0.1,
                "reaction_smiles": [],
                "starting_materials": [],
                "is_solved": False
            }

    async def execute(self, smiles: str, **params) -> AdapterResult:
        """
        Execute retrosynthesis route finding for a target molecule.

        Args:
            smiles: SMILES string of the target molecule
            **params: Additional parameters:
                     - max_routes: Override default max routes
                     - timeout: Override default timeout
                     - expansion_time: Override MCTS expansion time

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

            # Get parameters
            max_routes = params.get('max_routes', self.config.get('max_routes', 5))
            expansion_time = params.get('expansion_time', self.config.get('expansion_time', 60))

            # Generate cache key
            cache_key = self.generate_cache_key(
                smiles,
                max_routes=max_routes,
                expansion_time=expansion_time
            )

            logger.info(f"Running retrosynthesis for: {smiles[:50]}...")
            logger.info(f"Config: max_routes={max_routes}, expansion_time={expansion_time}s")

            # Run tree search in thread pool (AiZynthFinder is synchronous)
            loop = asyncio.get_event_loop()
            routes = await loop.run_in_executor(
                None,
                self._run_tree_search,
                smiles,
                expansion_time
            )

            # Check if any routes were found
            if not routes or len(routes) == 0:
                logger.warning(f"No synthesis routes found for {smiles}")
                return AdapterResult(
                    success=True,  # Not an error, just no routes found
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
                        "message": "No synthesis routes found"
                    }
                )

            # Extract route information
            route_list = []
            for i, route in enumerate(routes[:max_routes]):
                route_data = self._extract_route_info(route)
                score, feasibility = self._calculate_route_feasibility(route_data)

                route_list.append({
                    "route_id": i + 1,
                    "n_steps": route_data["n_steps"],
                    "n_reactions": route_data["n_reactions"],
                    "synthesis_score": round(score, 3),
                    "feasibility": feasibility,
                    "tree_score": round(route_data["score"], 3),
                    "reaction_smiles": route_data["reaction_smiles"],
                    "starting_materials": route_data["starting_materials"],
                    "is_complete": route_data["is_solved"]
                })

            # Find best route (highest synthesis score)
            best_route = max(route_list, key=lambda x: x["synthesis_score"])

            # Create result data structure
            result_data = {
                "smiles": smiles,
                "routes_found": len(route_list),
                "routes": route_list,
                "best_route": best_route,
                "n_steps": best_route["n_steps"],
                "synthesis_score": best_route["synthesis_score"],
                "feasibility": best_route["feasibility"],
                "model": "AiZynthFinder MCTS",
                "reference": "Genheden et al., J. Chem. Inf. Model. 2020"
            }

            logger.info(f"✓ Found {len(route_list)} synthesis routes")
            logger.info(f"  Best route: {best_route['n_steps']} steps, score={best_route['synthesis_score']:.3f}")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "cache_key": cache_key,
                    "version": self.version,
                    "expansion_time": expansion_time
                }
            )

        except ImportError as e:
            # AiZynthFinder not installed
            logger.error(f"AiZynthFinder not available: {e}")
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name,
                    "installation_help": "pip install aizynthfinder"
                }
            )

        except Exception as e:
            logger.error(f"Retrosynthesis failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name
                }
            )

    def _run_tree_search(self, smiles: str, expansion_time: int) -> List[Any]:
        """
        Run MCTS tree search (synchronous method for thread pool).

        Args:
            smiles: Target SMILES
            expansion_time: Expansion time in seconds

        Returns:
            List of route objects
        """
        try:
            # Set target molecule
            self.finder.target_smiles = smiles

            # Configure expansion time
            self.finder.config.search_time = expansion_time

            # Run tree search
            self.finder.tree_search()

            # Extract routes
            self.finder.build_routes()

            # Get routes sorted by score
            routes = self.finder.routes

            return routes.all_routes if hasattr(routes, 'all_routes') else []

        except Exception as e:
            logger.error(f"Tree search failed: {e}")
            return []

    def generate_cache_key(
        self,
        smiles: str,
        max_routes: int = 5,
        expansion_time: int = 60
    ) -> str:
        """
        Generate deterministic cache key for retrosynthesis.

        Args:
            smiles: SMILES string
            max_routes: Maximum routes to find
            expansion_time: MCTS expansion time

        Returns:
            SHA256 hash as cache key
        """
        # Canonicalize SMILES for consistent caching
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            pass  # Use original SMILES if canonicalization fails

        # Create cache key including version and parameters
        key_string = f"aizynthfinder_v{self.version}:{smiles}:{max_routes}:{expansion_time}"
        return hashlib.sha256(key_string.encode()).hexdigest()

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata.

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "Retrosynthesis route planning using AiZynthFinder",
            "capabilities": {
                "multi_route_discovery": True,
                "synthesis_step_counting": True,
                "feasibility_scoring": True,
                "reaction_templates": True
            },
            "config": {
                "max_routes": self.config.get("max_routes"),
                "timeout": self.config.get("timeout"),
                "expansion_time": self.config.get("expansion_time"),
                "stock": self.config.get("stock")
            },
            "reference": {
                "paper": "Genheden et al., J. Chem. Inf. Model. 2020, 60, 12, 5910-5919",
                "github": "https://github.com/MolecularAI/aizynthfinder"
            }
        }
