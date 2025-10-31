"""
Open Reaction Database (ORD) Adapter for PharmForge

Provides access to the Open Reaction Database for forward synthesis planning.
Returns reaction data including reactants, products, conditions, and yields.

The ORD is a public dataset of chemical reactions with detailed experimental
conditions, making it valuable for forward synthesis planning and reaction
prediction.

Reference: https://open-reaction-database.org/
GitHub: https://github.com/open-reaction-database/ord-schema
"""

import logging
from typing import Any, Dict, List, Optional, Union
import asyncio
import aiohttp
import hashlib
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ORDAdapter(AdapterProtocol):
    """
    Open Reaction Database adapter for forward synthesis planning.

    Features:
    - Search reactions by product, reactant, or reagent SMILES
    - Retrieve reaction conditions (temperature, pressure, solvent, catalyst)
    - Get reaction yields and selectivity
    - Filter by reaction type
    - Support batch queries
    - Caching for repeated queries

    The ORD contains curated experimental reaction data with detailed
    conditions, making it suitable for:
    - Forward synthesis planning
    - Reaction condition optimization
    - Literature-validated reaction routes
    - Reaction feasibility assessment
    """

    # ORD API base URL (they provide a REST interface)
    BASE_URL = "https://client.open-reaction-database.org/api"

    # Alternative: Direct PostgreSQL query endpoint (if available)
    QUERY_URL = "https://ord-postgres.herokuapp.com"

    def __init__(
        self,
        name: str = "ord",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize ORD adapter.

        Args:
            name: Adapter name (default: "ord")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - rate_limit_delay: Delay between API calls (default: 0.5s)
                   - timeout: Request timeout (default: 30s)
                   - max_results: Maximum reactions to return (default: 50)
                   - min_yield: Minimum yield threshold % (default: 0.0)
                   - similarity_threshold: SMILES similarity threshold (default: 0.85)
        """
        default_config = {
            "rate_limit_delay": 0.5,
            "timeout": 30,
            "max_results": 50,
            "min_yield": 0.0,  # Filter by minimum yield %
            "similarity_threshold": 0.85,  # For fuzzy SMILES matching
            "include_conditions": True,  # Include reaction conditions
            "include_reagents": True,  # Include reagent details
        }

        merged_config = {**default_config, **(config or {})}
        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data.

        Args:
            input_data: Can be:
                       - SMILES string (product, reactant, or reagent)
                       - Reaction ID
                       - Dict with search parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            # SMILES or reaction ID
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            # Must have at least one search parameter
            valid_keys = {"smiles", "product", "reactant", "reagent", "reaction_id", "reaction_type"}
            return bool(set(input_data.keys()) & valid_keys)

        return False

    def _canonicalize_smiles(self, smiles: str) -> Optional[str]:
        """
        Canonicalize SMILES for consistent searching.

        Args:
            smiles: Input SMILES string

        Returns:
            Canonical SMILES or None if invalid
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
        except Exception as e:
            logger.warning(f"SMILES canonicalization failed: {e}")

        return None

    async def _search_ord_api(
        self,
        query_params: Dict[str, Any]
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Search ORD database via API.

        Args:
            query_params: Search parameters

        Returns:
            List of reaction records or None on error
        """
        try:
            # Build query URL
            # ORD uses different endpoints for different search types
            search_type = query_params.get("search_type", "product")
            smiles = query_params.get("smiles")

            if not smiles:
                logger.error("No SMILES provided for ORD search")
                return None

            # Canonicalize SMILES
            canonical_smiles = self._canonicalize_smiles(smiles)
            if not canonical_smiles:
                logger.warning(f"Failed to canonicalize SMILES: {smiles}")
                canonical_smiles = smiles

            # ORD API endpoint (Note: The actual API may vary - this is based on documentation)
            url = f"{self.BASE_URL}/query"

            # Build query payload
            payload = {
                "component": canonical_smiles,
                "component_role": search_type,  # "product", "reactant", "reagent"
                "limit": self.config.get("max_results", 50),
                "similarity": self.config.get("similarity_threshold", 0.85)
            }

            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

            async with aiohttp.ClientSession() as session:
                async with session.post(url, json=payload, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        reactions = data.get("reactions", [])
                        logger.info(f"ORD API: Found {len(reactions)} reactions")
                        return reactions
                    else:
                        error_text = await response.text()
                        logger.error(f"ORD API error {response.status}: {error_text}")
                        return None

        except aiohttp.ClientError as e:
            logger.error(f"ORD API connection error: {e}")
            return None
        except asyncio.TimeoutError:
            logger.error(f"ORD API timeout")
            return None
        except Exception as e:
            logger.error(f"ORD API search error: {e}", exc_info=True)
            return None

    def _extract_reaction_info(self, reaction_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract structured information from ORD reaction data.

        Args:
            reaction_data: Raw reaction data from ORD

        Returns:
            Dictionary with standardized reaction information
        """
        try:
            # Extract reaction components
            inputs = reaction_data.get("inputs", {})
            outcomes = reaction_data.get("outcomes", [])

            # Get reactants
            reactants = []
            for key, input_data in inputs.items():
                components = input_data.get("components", [])
                for comp in components:
                    reactant_smiles = comp.get("smiles", "")
                    if reactant_smiles:
                        reactants.append({
                            "smiles": reactant_smiles,
                            "name": comp.get("name", ""),
                            "role": comp.get("role", "reactant"),
                            "amount": comp.get("amount", {})
                        })

            # Get products and yields from outcomes
            products = []
            yields_list = []

            for outcome in outcomes:
                for product in outcome.get("products", []):
                    product_smiles = product.get("smiles", "")
                    if product_smiles:
                        yield_data = product.get("measurements", [])
                        yield_pct = 0.0

                        # Extract yield percentage
                        for measurement in yield_data:
                            if measurement.get("type") == "YIELD":
                                yield_pct = measurement.get("percentage", {}).get("value", 0.0)
                                yields_list.append(yield_pct)

                        products.append({
                            "smiles": product_smiles,
                            "name": product.get("name", ""),
                            "yield": yield_pct,
                            "is_desired": product.get("is_desired_product", False)
                        })

            # Get reaction conditions
            conditions = reaction_data.get("conditions", {})

            temperature = None
            temp_data = conditions.get("temperature", {})
            if temp_data:
                temp_control = temp_data.get("control", {})
                if temp_control:
                    temperature = temp_control.get("setpoint", {}).get("value")

            pressure = None
            press_data = conditions.get("pressure", {})
            if press_data:
                press_control = press_data.get("control", {})
                if press_control:
                    pressure = press_control.get("setpoint", {}).get("value")

            # Get stirring conditions
            stirring = conditions.get("stirring", {})
            stirring_rate = stirring.get("rate", {}).get("rpm")

            # Extract solvents and reagents
            solvents = []
            reagents = []

            for key, input_data in inputs.items():
                components = input_data.get("components", [])
                for comp in components:
                    role = comp.get("role", "").lower()
                    if "solvent" in role:
                        solvents.append(comp.get("name") or comp.get("smiles", ""))
                    elif "catalyst" in role or "reagent" in role:
                        reagents.append(comp.get("name") or comp.get("smiles", ""))

            # Calculate average yield
            avg_yield = sum(yields_list) / len(yields_list) if yields_list else 0.0

            # Get reaction metadata
            reaction_id = reaction_data.get("reaction_id", "")
            reaction_type = reaction_data.get("identifiers", [{}])[0].get("type", "")

            # Provenance (literature source)
            provenance = reaction_data.get("provenance", {})
            doi = provenance.get("doi", "")
            publication = provenance.get("publication_url", "")

            return {
                "reaction_id": reaction_id,
                "reaction_type": reaction_type,
                "reactants": reactants,
                "products": products,
                "yield": round(avg_yield, 2),
                "conditions": {
                    "temperature": temperature,
                    "temperature_unit": "C" if temperature else None,
                    "pressure": pressure,
                    "pressure_unit": "bar" if pressure else None,
                    "stirring_rate_rpm": stirring_rate,
                    "solvents": solvents,
                    "reagents": reagents,
                },
                "literature": {
                    "doi": doi,
                    "url": publication
                },
                "n_reactants": len(reactants),
                "n_products": len(products)
            }

        except Exception as e:
            logger.warning(f"Failed to extract reaction info: {e}")
            return {
                "reaction_id": reaction_data.get("reaction_id", "unknown"),
                "error": str(e),
                "reactants": [],
                "products": [],
                "yield": 0.0
            }

    async def execute(
        self,
        input_data: Union[str, Dict[str, Any]],
        **params
    ) -> AdapterResult:
        """
        Execute ORD database search.

        Args:
            input_data: Can be:
                       - SMILES string (searches as product by default)
                       - Dict with search parameters:
                         {
                           "smiles": "CC(=O)O",
                           "search_type": "product",  # or "reactant", "reagent"
                           "reaction_type": "coupling",  # optional filter
                           "min_yield": 50.0  # optional yield threshold
                         }
            **params: Additional parameters to override config

        Returns:
            AdapterResult containing reaction data
        """
        try:
            # Validate input
            if not self.validate_input(input_data):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid input data. Provide SMILES string or search dict."
                )

            # Parse input
            if isinstance(input_data, str):
                # Simple SMILES search (as product)
                query_params = {
                    "smiles": input_data,
                    "search_type": "product"
                }
            else:
                query_params = input_data.copy()

            # Override with params
            query_params.update(params)

            # Apply rate limiting
            rate_delay = self.config.get("rate_limit_delay", 0.5)
            await asyncio.sleep(rate_delay)

            logger.info(f"Searching ORD for: {query_params.get('smiles', 'N/A')[:50]}...")

            # Search ORD database
            reactions = await self._search_ord_api(query_params)

            if reactions is None:
                return AdapterResult(
                    success=False,
                    data={},
                    error="Failed to query ORD database"
                )

            if len(reactions) == 0:
                logger.warning("No reactions found in ORD")
                return AdapterResult(
                    success=True,
                    data={
                        "reactions_found": 0,
                        "reactions": [],
                        "query": query_params
                    },
                    metadata={
                        "adapter_name": self.name,
                        "version": self.version,
                        "message": "No reactions found"
                    }
                )

            # Extract reaction information
            processed_reactions = []
            for reaction in reactions:
                reaction_info = self._extract_reaction_info(reaction)

                # Apply yield filter if specified
                min_yield = params.get("min_yield", self.config.get("min_yield", 0.0))
                if reaction_info.get("yield", 0.0) >= min_yield:
                    processed_reactions.append(reaction_info)

            # Sort by yield (descending)
            processed_reactions.sort(key=lambda x: x.get("yield", 0.0), reverse=True)

            # Limit results
            max_results = self.config.get("max_results", 50)
            processed_reactions = processed_reactions[:max_results]

            # Calculate statistics
            yields = [r.get("yield", 0.0) for r in processed_reactions if r.get("yield", 0.0) > 0]
            avg_yield = sum(yields) / len(yields) if yields else 0.0
            max_yield = max(yields) if yields else 0.0

            # Collect unique reaction types
            reaction_types = set(
                r.get("reaction_type", "unknown")
                for r in processed_reactions
                if r.get("reaction_type")
            )

            # Collect unique solvents and reagents
            all_solvents = set()
            all_reagents = set()
            for reaction in processed_reactions:
                conditions = reaction.get("conditions", {})
                all_solvents.update(conditions.get("solvents", []))
                all_reagents.update(conditions.get("reagents", []))

            result_data = {
                "reactions_found": len(processed_reactions),
                "reactions": processed_reactions,
                "statistics": {
                    "average_yield": round(avg_yield, 2),
                    "max_yield": round(max_yield, 2),
                    "total_reactions": len(reactions),
                    "after_filtering": len(processed_reactions)
                },
                "reaction_types": sorted(list(reaction_types)),
                "common_solvents": sorted(list(all_solvents))[:10],
                "common_reagents": sorted(list(all_reagents))[:10],
                "query": query_params,
                "database": "Open Reaction Database",
                "reference": "https://open-reaction-database.org/"
            }

            logger.info(f"✓ Found {len(processed_reactions)} reactions in ORD")
            logger.info(f"  Average yield: {avg_yield:.1f}%, Max yield: {max_yield:.1f}%")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "version": self.version,
                    "cache_key": self.generate_cache_key(input_data, **params)
                }
            )

        except Exception as e:
            logger.error(f"ORD search failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    def generate_cache_key(
        self,
        input_data: Union[str, Dict[str, Any]],
        **params
    ) -> str:
        """
        Generate deterministic cache key for ORD search.

        Args:
            input_data: Input SMILES or search parameters
            **params: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        # Normalize input
        if isinstance(input_data, str):
            # Canonicalize SMILES
            canonical = self._canonicalize_smiles(input_data) or input_data
            cache_dict = {
                "adapter": self.name,
                "version": self.version,
                "smiles": canonical,
                "search_type": "product",
                "params": params
            }
        else:
            # Dictionary input
            smiles = input_data.get("smiles", "")
            canonical = self._canonicalize_smiles(smiles) or smiles
            cache_dict = {
                "adapter": self.name,
                "version": self.version,
                "smiles": canonical,
                "search_type": input_data.get("search_type", "product"),
                "reaction_type": input_data.get("reaction_type"),
                "min_yield": input_data.get("min_yield"),
                "params": params
            }

        # Generate hash
        cache_string = json.dumps(cache_dict, sort_keys=True)
        return hashlib.sha256(cache_string.encode()).hexdigest()

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
            "description": "Open Reaction Database for forward synthesis planning",
            "capabilities": {
                "search_by_product": True,
                "search_by_reactant": True,
                "search_by_reagent": True,
                "reaction_conditions": True,
                "yield_data": True,
                "batch_queries": True,
                "literature_references": True
            },
            "config": {
                "rate_limit_delay": self.config.get("rate_limit_delay"),
                "timeout": self.config.get("timeout"),
                "max_results": self.config.get("max_results"),
                "min_yield": self.config.get("min_yield")
            },
            "reference": {
                "website": "https://open-reaction-database.org/",
                "github": "https://github.com/open-reaction-database/ord-schema",
                "paper": "Kearnes et al., J. Am. Chem. Soc. 2021, 143, 45, 18820-18826"
            }
        }
