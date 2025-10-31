"""
STRING-DB Adapter - Protein-protein interaction network queries
API Documentation: https://string-db.org/help/api/

Provides access to:
- Protein-protein interaction networks
- Interaction confidence scores
- Multiple evidence channels (experimental, database, textmining, etc.)
- Network enrichment analysis
- Functional annotations
"""
from typing import Any, Dict, Optional, List, Union
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class StringDBAdapter(AdapterProtocol):
    """
    Adapter for STRING-DB Protein-Protein Interaction Network Database

    Free API access to protein interaction networks with confidence scores
    across multiple evidence channels.
    """

    BASE_URL = "https://string-db.org/api"

    # Common organism NCBI taxonomy IDs
    ORGANISMS = {
        "human": "9606",
        "homo_sapiens": "9606",
        "mouse": "10090",
        "mus_musculus": "10090",
        "rat": "10116",
        "rattus_norvegicus": "10116",
        "zebrafish": "7955",
        "danio_rerio": "7955",
        "fruit_fly": "7227",
        "drosophila_melanogaster": "7227",
        "c_elegans": "6239",
        "caenorhabditis_elegans": "6239",
        "yeast": "4932",
        "saccharomyces_cerevisiae": "4932",
        "e_coli": "511145",
        "escherichia_coli": "511145",
        "arabidopsis": "3702",
        "arabidopsis_thaliana": "3702"
    }

    def __init__(self):
        super().__init__(
            name="string_db",
            adapter_type="api",
            config={
                "rate_limit_delay": 1.0,  # STRING requires 1 second between requests
                "timeout": 60,
                "format": "json",
                "default_organism": "9606",  # Human
                "default_score": 400,  # Medium confidence (0-1000)
                "network_type": "functional"  # functional or physical
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data for STRING-DB queries

        Args:
            input_data: Protein identifier(s), gene name(s), or dict with parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, list):
            return len(input_data) > 0 and all(isinstance(x, str) for x in input_data)
        elif isinstance(input_data, dict):
            return ("identifiers" in input_data or
                    "proteins" in input_data or
                    "string_ids" in input_data)
        return False

    def _get_organism_id(self, organism: Union[str, int]) -> str:
        """
        Convert organism name or ID to NCBI taxonomy ID

        Args:
            organism: Organism name (e.g., "human", "mouse") or NCBI taxonomy ID

        Returns:
            NCBI taxonomy ID as string
        """
        if isinstance(organism, int):
            return str(organism)

        # Check if it's already a numeric ID
        if isinstance(organism, str) and organism.isdigit():
            return organism

        # Look up common name
        organism_lower = organism.lower().replace(" ", "_")
        if organism_lower in self.ORGANISMS:
            return self.ORGANISMS[organism_lower]

        # Default to human if not found
        logger.warning(f"Unknown organism '{organism}', defaulting to human (9606)")
        return self.config.get("default_organism", "9606")

    async def _map_identifiers(
        self,
        identifiers: List[str],
        species: str,
        limit: int = 1
    ) -> Optional[Dict[str, Any]]:
        """
        Map gene names or identifiers to STRING IDs

        Args:
            identifiers: List of protein/gene identifiers
            species: NCBI taxonomy ID
            limit: Maximum number of matches per identifier

        Returns:
            Mapping results or None on error
        """
        url = f"{self.BASE_URL}/json/get_string_ids"

        # Join identifiers with newline
        identifiers_str = "\r".join(identifiers)

        params = {
            "identifiers": identifiers_str,
            "species": species,
            "limit": limit,
            "echo_query": 1
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(url, data=params, timeout=timeout) as response:
                    if response.status == 200:
                        # STRING API returns 'text/json' which aiohttp doesn't recognize
                        data = await response.json(content_type=None)
                        return data
                    else:
                        error_text = await response.text()
                        logger.error(f"STRING-DB mapping error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"STRING-DB: Timeout mapping identifiers")
            return None
        except Exception as e:
            logger.error(f"STRING-DB: Error mapping identifiers: {e}")
            return None

    async def _get_network(
        self,
        identifiers: List[str],
        species: str,
        required_score: int = 400,
        network_type: str = "functional",
        add_nodes: int = 0
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get protein-protein interaction network

        Args:
            identifiers: List of protein identifiers (preferably STRING IDs)
            species: NCBI taxonomy ID
            required_score: Minimum confidence score (0-1000)
            network_type: "functional" or "physical"
            add_nodes: Number of additional nodes to add to network

        Returns:
            List of interactions or None on error
        """
        url = f"{self.BASE_URL}/json/network"

        # Join identifiers with newline
        identifiers_str = "\r".join(identifiers)

        params = {
            "identifiers": identifiers_str,
            "species": species,
            "required_score": required_score,
            "network_type": network_type
        }

        if add_nodes > 0:
            params["add_nodes"] = add_nodes

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(url, data=params, timeout=timeout) as response:
                    if response.status == 200:
                        # STRING API returns 'text/json' which aiohttp doesn't recognize
                        data = await response.json(content_type=None)
                        return data
                    else:
                        error_text = await response.text()
                        logger.error(f"STRING-DB network error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"STRING-DB: Timeout fetching network")
            return None
        except Exception as e:
            logger.error(f"STRING-DB: Error fetching network: {e}")
            return None

    async def _get_interaction_partners(
        self,
        identifiers: List[str],
        species: str,
        required_score: int = 400,
        limit: int = 10
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get interaction partners for given proteins

        Args:
            identifiers: List of protein identifiers
            species: NCBI taxonomy ID
            required_score: Minimum confidence score (0-1000)
            limit: Maximum number of partners per protein

        Returns:
            List of interaction partners or None on error
        """
        url = f"{self.BASE_URL}/json/interaction_partners"

        identifiers_str = "\r".join(identifiers)

        params = {
            "identifiers": identifiers_str,
            "species": species,
            "required_score": required_score,
            "limit": limit
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(url, data=params, timeout=timeout) as response:
                    if response.status == 200:
                        # STRING API returns 'text/json' which aiohttp doesn't recognize
                        data = await response.json(content_type=None)
                        return data
                    else:
                        error_text = await response.text()
                        logger.error(f"STRING-DB partners error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"STRING-DB: Timeout fetching partners")
            return None
        except Exception as e:
            logger.error(f"STRING-DB: Error fetching partners: {e}")
            return None

    async def _get_enrichment(
        self,
        identifiers: List[str],
        species: str
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get functional enrichment for protein list

        Args:
            identifiers: List of protein identifiers
            species: NCBI taxonomy ID

        Returns:
            Enrichment results or None on error
        """
        url = f"{self.BASE_URL}/json/enrichment"

        identifiers_str = "\r".join(identifiers)

        params = {
            "identifiers": identifiers_str,
            "species": species
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(url, data=params, timeout=timeout) as response:
                    if response.status == 200:
                        # STRING API returns 'text/json' which aiohttp doesn't recognize
                        data = await response.json(content_type=None)
                        return data
                    else:
                        error_text = await response.text()
                        logger.error(f"STRING-DB enrichment error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"STRING-DB: Timeout fetching enrichment")
            return None
        except Exception as e:
            logger.error(f"STRING-DB: Error fetching enrichment: {e}")
            return None

    async def _get_ppi_enrichment(
        self,
        identifiers: List[str],
        species: str,
        required_score: int = 400
    ) -> Optional[Dict[str, Any]]:
        """
        Test if protein network has more interactions than expected

        Args:
            identifiers: List of protein identifiers
            species: NCBI taxonomy ID
            required_score: Minimum confidence score (0-1000)

        Returns:
            PPI enrichment p-value or None on error
        """
        url = f"{self.BASE_URL}/json/ppi_enrichment"

        identifiers_str = "\r".join(identifiers)

        params = {
            "identifiers": identifiers_str,
            "species": species,
            "required_score": required_score
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(url, data=params, timeout=timeout) as response:
                    if response.status == 200:
                        # STRING API returns 'text/json' which aiohttp doesn't recognize
                        data = await response.json(content_type=None)
                        return data
                    else:
                        error_text = await response.text()
                        logger.error(f"STRING-DB PPI enrichment error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"STRING-DB: Timeout fetching PPI enrichment")
            return None
        except Exception as e:
            logger.error(f"STRING-DB: Error fetching PPI enrichment: {e}")
            return None

    def _format_interaction_data(
        self,
        interactions: List[Dict[str, Any]],
        include_scores: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Format interaction data into standardized structure

        Args:
            interactions: Raw interaction data from STRING
            include_scores: Include individual evidence scores

        Returns:
            Formatted interaction list
        """
        formatted = []

        for interaction in interactions:
            # Basic interaction info
            interaction_data = {
                "protein_a": interaction.get("preferredName_A") or interaction.get("stringId_A"),
                "protein_b": interaction.get("preferredName_B") or interaction.get("stringId_B"),
                "combined_score": interaction.get("score", 0) / 1000.0,  # Normalize to 0-1
                "string_id_a": interaction.get("stringId_A"),
                "string_id_b": interaction.get("stringId_B")
            }

            # Add evidence scores if requested
            if include_scores:
                interaction_data["evidence_scores"] = {
                    "neighborhood": interaction.get("nscore", 0) / 1000.0,
                    "fusion": interaction.get("fscore", 0) / 1000.0,
                    "cooccurrence": interaction.get("pscore", 0) / 1000.0,
                    "coexpression": interaction.get("ascore", 0) / 1000.0,
                    "experimental": interaction.get("escore", 0) / 1000.0,
                    "database": interaction.get("dscore", 0) / 1000.0,
                    "textmining": interaction.get("tscore", 0) / 1000.0
                }

            formatted.append(interaction_data)

        # Sort by combined score (highest first)
        formatted.sort(key=lambda x: x["combined_score"], reverse=True)

        return formatted

    def _format_enrichment_data(
        self,
        enrichment: List[Dict[str, Any]],
        max_results: int = 20
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Format enrichment data by category

        Args:
            enrichment: Raw enrichment data from STRING
            max_results: Maximum results per category

        Returns:
            Enrichment results grouped by category
        """
        categories = {}

        # Handle case where enrichment might not be a list
        if not isinstance(enrichment, list):
            logger.warning(f"Enrichment data is not a list: {type(enrichment)}")
            return categories

        for item in enrichment:
            # Skip if item is not a dict
            if not isinstance(item, dict):
                logger.warning(f"Enrichment item is not a dict: {type(item)}")
                continue

            category = item.get("category", "Unknown")

            if category not in categories:
                categories[category] = []

            # Only add if not too many already
            if len(categories[category]) < max_results:
                # Handle genes - can be string or list
                input_genes = item.get("inputGenes", [])
                if isinstance(input_genes, str):
                    genes = input_genes.split(",") if input_genes else []
                elif isinstance(input_genes, list):
                    genes = input_genes
                else:
                    genes = []

                enrichment_item = {
                    "term": item.get("term"),
                    "description": item.get("description"),
                    "number_of_genes": item.get("number_of_genes"),
                    "fdr": item.get("fdr"),  # False discovery rate
                    "p_value": item.get("p_value"),
                    "genes": genes
                }
                categories[category].append(enrichment_item)

        # Sort each category by p-value (most significant first)
        for category in categories:
            categories[category].sort(key=lambda x: x.get("p_value", 1.0))

        return categories

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute STRING-DB query for protein interaction networks

        Args:
            input_data: Protein identifier(s), gene name(s), or dict with parameters
                       Examples:
                       - "TP53" (single protein)
                       - ["TP53", "MDM2", "BRCA1"] (protein list)
                       - {"identifiers": ["TP53", "MDM2"], "species": "human"}
                       - {"string_ids": ["9606.ENSP00000269305"], "species": 9606}
            **kwargs: Additional parameters:
                     - species: str/int - Organism (name or NCBI taxonomy ID) [default: human/9606]
                     - required_score: int - Minimum confidence (0-1000) [default: 400]
                     - network_type: str - "functional" or "physical" [default: functional]
                     - add_nodes: int - Additional nodes to add to network [default: 0]
                     - include_partners: bool - Include interaction partners [default: False]
                     - include_enrichment: bool - Include functional enrichment [default: False]
                     - include_ppi_enrichment: bool - Include PPI enrichment p-value [default: False]
                     - limit_partners: int - Max partners per protein [default: 10]
                     - include_evidence_scores: bool - Include individual evidence scores [default: True]

        Returns:
            AdapterResult containing interaction network data with confidence scores
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be protein identifier(s) or dict with 'identifiers'"
            )

        # Rate limiting - STRING requires 1 second between requests
        rate_delay = self.config.get("rate_limit_delay", 1.0)
        await asyncio.sleep(rate_delay)

        # Parse input
        identifiers = []
        species = kwargs.get("species", self.config.get("default_organism", "9606"))

        if isinstance(input_data, dict):
            identifiers = (input_data.get("identifiers") or
                          input_data.get("proteins") or
                          input_data.get("string_ids"))
            species = input_data.get("species", species)
        elif isinstance(input_data, list):
            identifiers = input_data
        elif isinstance(input_data, str):
            identifiers = [input_data]

        # Convert species to taxonomy ID
        species_id = self._get_organism_id(species)

        # Get parameters
        required_score = kwargs.get("required_score", self.config.get("default_score", 400))
        network_type = kwargs.get("network_type", self.config.get("network_type", "functional"))
        add_nodes = kwargs.get("add_nodes", 0)
        include_partners = kwargs.get("include_partners", False)
        include_enrichment = kwargs.get("include_enrichment", False)
        include_ppi_enrichment = kwargs.get("include_ppi_enrichment", False)
        limit_partners = kwargs.get("limit_partners", 10)
        include_evidence_scores = kwargs.get("include_evidence_scores", True)

        try:
            # Step 1: Map identifiers to STRING IDs (if not already STRING IDs)
            if not any(id.startswith(f"{species_id}.") for id in identifiers):
                logger.info(f"STRING-DB: Mapping {len(identifiers)} identifiers for species {species_id}")
                mapping_result = await self._map_identifiers(identifiers, species_id)

                if not mapping_result:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="Failed to map protein identifiers",
                        metadata={"identifiers": identifiers, "species": species_id}
                    )

                # Extract STRING IDs
                string_ids = [item.get("stringId") for item in mapping_result if item.get("stringId")]

                if not string_ids:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"No proteins found for identifiers: {identifiers}",
                        metadata={"identifiers": identifiers, "species": species_id}
                    )

                # Create mapping info
                protein_mapping = {
                    item.get("queryItem"): {
                        "string_id": item.get("stringId"),
                        "preferred_name": item.get("preferredName"),
                        "annotation": item.get("annotation")
                    }
                    for item in mapping_result if item.get("stringId")
                }
            else:
                string_ids = identifiers
                protein_mapping = {}

            await asyncio.sleep(rate_delay)  # Rate limit

            # Step 2: Get interaction network
            logger.info(f"STRING-DB: Fetching network for {len(string_ids)} proteins")
            network_data = await self._get_network(
                string_ids,
                species_id,
                required_score,
                network_type,
                add_nodes
            )

            if network_data is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to fetch interaction network",
                    metadata={"string_ids": string_ids, "species": species_id}
                )

            # Format network data
            interactions = self._format_interaction_data(network_data, include_evidence_scores)

            result_data = {
                "interactions": interactions,
                "num_interactions": len(interactions),
                "protein_mapping": protein_mapping,
                "query_proteins": identifiers,
                "string_ids": string_ids,
                "parameters": {
                    "species": species_id,
                    "required_score": required_score,
                    "network_type": network_type
                }
            }

            # Step 3: Get interaction partners if requested
            if include_partners and string_ids:
                await asyncio.sleep(rate_delay)
                logger.info(f"STRING-DB: Fetching interaction partners")
                partners_data = await self._get_interaction_partners(
                    string_ids,
                    species_id,
                    required_score,
                    limit_partners
                )

                if partners_data:
                    result_data["interaction_partners"] = self._format_interaction_data(
                        partners_data,
                        include_evidence_scores
                    )

            # Step 4: Get functional enrichment if requested
            if include_enrichment and string_ids:
                await asyncio.sleep(rate_delay)
                logger.info(f"STRING-DB: Fetching functional enrichment")
                enrichment_data = await self._get_enrichment(string_ids, species_id)

                if enrichment_data:
                    result_data["enrichment"] = self._format_enrichment_data(enrichment_data)

            # Step 5: Get PPI enrichment if requested
            if include_ppi_enrichment and string_ids:
                await asyncio.sleep(rate_delay)
                logger.info(f"STRING-DB: Fetching PPI enrichment")
                ppi_data = await self._get_ppi_enrichment(string_ids, species_id, required_score)

                if ppi_data:
                    # PPI enrichment might return a list with one dict
                    if isinstance(ppi_data, list) and len(ppi_data) > 0:
                        ppi_data = ppi_data[0]

                    if isinstance(ppi_data, dict):
                        result_data["ppi_enrichment"] = {
                            "p_value": ppi_data.get("p_value"),
                            "number_of_nodes": ppi_data.get("number_of_nodes"),
                            "number_of_edges": ppi_data.get("number_of_edges"),
                            "expected_edges": ppi_data.get("expected_number_of_edges")
                        }
                    else:
                        logger.warning(f"Unexpected PPI enrichment data type: {type(ppi_data)}")

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "string_db",
                    "species": species_id,
                    "num_proteins": len(string_ids),
                    "num_interactions": len(interactions),
                    "adapter_version": self.version
                }
            )

        except Exception as e:
            logger.error(f"STRING-DB: Error processing request: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error processing STRING-DB query: {str(e)}",
                metadata={"identifiers": identifiers, "species": species_id}
            )
