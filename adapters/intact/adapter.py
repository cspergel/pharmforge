"""
IntAct Adapter - Protein-protein interaction database queries
API Documentation: https://www.ebi.ac.uk/intact/ws/

Provides access to:
- Curated protein-protein interactions from published literature
- Experimental evidence and detection methods
- Interaction confidence scores (MI scores)
- Interaction types (binary, complex, enzymatic reaction)
- Cross-references to other databases
- Publication metadata
"""
from typing import Any, Dict, Optional, List, Union
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class IntActAdapter(AdapterProtocol):
    """
    Adapter for IntAct Molecular Interaction Database

    Free API access to manually curated protein-protein interactions
    from published literature with detailed experimental evidence.
    """

    BASE_URL = "https://www.ebi.ac.uk/intact/ws"

    # Interaction types supported by IntAct
    INTERACTION_TYPES = [
        "physical association",
        "direct interaction",
        "association",
        "colocalization",
        "enzymatic reaction",
        "phosphorylation",
        "ubiquitination",
        "methylation"
    ]

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
        "e_coli": "562",
        "escherichia_coli": "562"
    }

    def __init__(self):
        super().__init__(
            name="intact",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,  # Be polite to EBI servers
                "timeout": 60,
                "format": "json",
                "default_organism": "9606",  # Human
                "min_mi_score": 0.4,  # Default MI score threshold (0-1)
                "max_results": 200,
                "negative_interactions": False  # Include negative results by default
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data for IntAct queries

        Args:
            input_data: Protein identifier(s), UniProt ID(s), or dict with parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, list):
            return len(input_data) > 0 and all(isinstance(x, str) for x in input_data)
        elif isinstance(input_data, dict):
            return ("identifiers" in input_data or
                    "uniprot_ids" in input_data or
                    "gene_names" in input_data or
                    "query" in input_data)
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

    async def _search_interactions(
        self,
        query: str,
        format: str = "json",
        first_result: int = 0,
        max_results: int = 200,
        negative: bool = False
    ) -> Optional[Dict[str, Any]]:
        """
        Search for interactions using PSICQUIC query syntax

        Args:
            query: PSICQUIC query (e.g., "P12345", "gene:TP53", "taxid:9606")
            format: Output format ("json", "xml", "tab25", "tab27")
            first_result: Index of first result (for pagination)
            max_results: Maximum number of results
            negative: Include negative interactions

        Returns:
            Interaction results or None on error
        """
        # Use PSICQUIC endpoint for flexible queries
        url = f"{self.BASE_URL}/psicquic/current/search/query/{quote(query)}"

        params = {
            "format": format,
            "firstResult": first_result,
            "maxResults": max_results
        }

        if negative:
            params["negative"] = "true"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        if format == "json":
                            data = await response.json()
                            return data
                        else:
                            text = await response.text()
                            return {"raw": text}
                    elif response.status == 404:
                        # No results found
                        logger.info(f"IntAct: No interactions found for query: {query}")
                        return {"data": []}
                    else:
                        error_text = await response.text()
                        logger.error(f"IntAct search error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"IntAct: Timeout searching interactions")
            return None
        except Exception as e:
            logger.error(f"IntAct: Error searching interactions: {e}")
            return None

    async def _get_interaction_details(
        self,
        interaction_ac: str
    ) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a specific interaction

        Args:
            interaction_ac: IntAct interaction accession (e.g., "EBI-123456")

        Returns:
            Detailed interaction data or None on error
        """
        url = f"{self.BASE_URL}/interaction/{interaction_ac}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    else:
                        error_text = await response.text()
                        logger.error(f"IntAct details error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"IntAct: Timeout fetching interaction details")
            return None
        except Exception as e:
            logger.error(f"IntAct: Error fetching interaction details: {e}")
            return None

    async def _search_by_interactor(
        self,
        interactor_ac: str,
        format: str = "json"
    ) -> Optional[Dict[str, Any]]:
        """
        Search interactions by interactor accession

        Args:
            interactor_ac: Interactor accession (UniProt ID or IntAct AC)
            format: Output format

        Returns:
            Interaction results or None on error
        """
        url = f"{self.BASE_URL}/interaction/findInteractionDetails/{interactor_ac}"

        params = {"format": format}

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    elif response.status == 404:
                        logger.info(f"IntAct: No interactions found for interactor: {interactor_ac}")
                        return {"data": []}
                    else:
                        error_text = await response.text()
                        logger.error(f"IntAct interactor error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"IntAct: Timeout searching by interactor")
            return None
        except Exception as e:
            logger.error(f"IntAct: Error searching by interactor: {e}")
            return None

    def _parse_mitab_line(self, line: str) -> Optional[Dict[str, Any]]:
        """
        Parse MITAB format interaction line

        Args:
            line: MITAB format line (tab-separated values)

        Returns:
            Parsed interaction data or None if invalid
        """
        try:
            fields = line.strip().split('\t')

            if len(fields) < 15:
                return None

            # Parse basic fields (MITAB 2.5 format)
            interaction = {
                "interactor_a_id": fields[0] if len(fields) > 0 else "",
                "interactor_b_id": fields[1] if len(fields) > 1 else "",
                "alt_id_a": fields[2] if len(fields) > 2 else "",
                "alt_id_b": fields[3] if len(fields) > 3 else "",
                "alias_a": fields[4] if len(fields) > 4 else "",
                "alias_b": fields[5] if len(fields) > 5 else "",
                "detection_method": fields[6] if len(fields) > 6 else "",
                "author": fields[7] if len(fields) > 7 else "",
                "publication": fields[8] if len(fields) > 8 else "",
                "organism_a": fields[9] if len(fields) > 9 else "",
                "organism_b": fields[10] if len(fields) > 10 else "",
                "interaction_type": fields[11] if len(fields) > 11 else "",
                "source_database": fields[12] if len(fields) > 12 else "",
                "interaction_ac": fields[13] if len(fields) > 13 else "",
                "confidence_score": fields[14] if len(fields) > 14 else ""
            }

            # Parse extended fields (MITAB 2.7+ format)
            if len(fields) > 15:
                interaction["expansion_method"] = fields[15] if len(fields) > 15 else ""
                interaction["biological_role_a"] = fields[16] if len(fields) > 16 else ""
                interaction["biological_role_b"] = fields[17] if len(fields) > 17 else ""
                interaction["experimental_role_a"] = fields[18] if len(fields) > 18 else ""
                interaction["experimental_role_b"] = fields[19] if len(fields) > 19 else ""
                interaction["interactor_type_a"] = fields[20] if len(fields) > 20 else ""
                interaction["interactor_type_b"] = fields[21] if len(fields) > 21 else ""
                interaction["xrefs_a"] = fields[22] if len(fields) > 22 else ""
                interaction["xrefs_b"] = fields[23] if len(fields) > 23 else ""

            return interaction

        except Exception as e:
            logger.error(f"Error parsing MITAB line: {e}")
            return None

    def _extract_mi_score(self, confidence_str: str) -> float:
        """
        Extract MI score from confidence string

        Args:
            confidence_str: Confidence score string (e.g., "intact-miscore:0.45")

        Returns:
            MI score as float (0-1), or 0.0 if not found
        """
        try:
            if not confidence_str:
                return 0.0

            # Split by pipe for multiple scores
            scores = confidence_str.split('|')

            for score in scores:
                if 'intact-miscore' in score.lower():
                    # Extract numeric value
                    parts = score.split(':')
                    if len(parts) > 1:
                        return float(parts[1])

            return 0.0

        except Exception as e:
            logger.warning(f"Error extracting MI score: {e}")
            return 0.0

    def _format_interaction_data(
        self,
        interactions: List[Dict[str, Any]],
        min_mi_score: float = 0.0
    ) -> List[Dict[str, Any]]:
        """
        Format interaction data into standardized structure

        Args:
            interactions: Raw interaction data from IntAct
            min_mi_score: Minimum MI score threshold (0-1)

        Returns:
            Formatted interaction list
        """
        formatted = []

        for interaction in interactions:
            # Extract MI score
            mi_score = self._extract_mi_score(interaction.get("confidence_score", ""))

            # Filter by MI score if specified
            if mi_score < min_mi_score:
                continue

            # Extract identifiers
            protein_a = interaction.get("interactor_a_id", "").split(':')[-1]
            protein_b = interaction.get("interactor_b_id", "").split(':')[-1]

            # Extract names from aliases
            alias_a = interaction.get("alias_a", "").split(':')[-1] if interaction.get("alias_a") else protein_a
            alias_b = interaction.get("alias_b", "").split(':')[-1] if interaction.get("alias_b") else protein_b

            # Extract method and type
            detection_method = interaction.get("detection_method", "").split('(')[-1].split(')')[0] if interaction.get("detection_method") else "unknown"
            interaction_type = interaction.get("interaction_type", "").split('(')[-1].split(')')[0] if interaction.get("interaction_type") else "unknown"

            # Extract publication
            publication = interaction.get("publication", "").split(':')[-1] if interaction.get("publication") else ""

            interaction_data = {
                "protein_a": protein_a,
                "protein_b": protein_b,
                "protein_a_name": alias_a,
                "protein_b_name": alias_b,
                "mi_score": mi_score,
                "detection_method": detection_method,
                "interaction_type": interaction_type,
                "publication": publication,
                "interaction_ac": interaction.get("interaction_ac", "").split(':')[-1],
                "source": interaction.get("source_database", "").split(':')[-1] if interaction.get("source_database") else "intact",
                "organism_a": interaction.get("organism_a", "").split('(')[-1].split(')')[0] if interaction.get("organism_a") else "",
                "organism_b": interaction.get("organism_b", "").split('(')[-1].split(')')[0] if interaction.get("organism_b") else ""
            }

            formatted.append(interaction_data)

        # Sort by MI score (highest first)
        formatted.sort(key=lambda x: x["mi_score"], reverse=True)

        return formatted

    def _build_query(
        self,
        identifiers: List[str],
        organism: Optional[str] = None,
        interaction_type: Optional[str] = None,
        detection_method: Optional[str] = None
    ) -> str:
        """
        Build PSICQUIC query string

        Args:
            identifiers: List of protein identifiers
            organism: NCBI taxonomy ID
            interaction_type: Interaction type to filter
            detection_method: Detection method to filter

        Returns:
            PSICQUIC query string
        """
        query_parts = []

        # Add identifiers (OR logic)
        if identifiers:
            id_queries = []
            for identifier in identifiers:
                # Try different identifier formats
                if identifier.startswith('P') or identifier.startswith('Q'):
                    # Likely UniProt ID
                    id_queries.append(f"id:{identifier}")
                elif identifier.startswith('EBI-'):
                    # IntAct AC
                    id_queries.append(f"idA:{identifier} OR idB:{identifier}")
                else:
                    # Gene name or other
                    id_queries.append(f"alias:{identifier}")

            if id_queries:
                query_parts.append(f"({' OR '.join(id_queries)})")

        # Add organism filter
        if organism:
            query_parts.append(f"taxid:{organism}")

        # Add interaction type filter
        if interaction_type:
            query_parts.append(f'type:"{interaction_type}"')

        # Add detection method filter
        if detection_method:
            query_parts.append(f'detmethod:"{detection_method}"')

        # Combine with AND logic
        return ' AND '.join(query_parts) if query_parts else "*"

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute IntAct query for protein interaction data

        Args:
            input_data: Protein identifier(s), UniProt ID(s), or dict with parameters
                       Examples:
                       - "P04637" (UniProt ID for TP53)
                       - ["P04637", "Q00987"] (UniProt IDs)
                       - {"identifiers": ["TP53", "MDM2"], "organism": "human"}
                       - {"uniprot_ids": ["P04637"], "min_mi_score": 0.5}
                       - {"query": "id:P04637 AND taxid:9606"}
            **kwargs: Additional parameters:
                     - organism: str/int - Organism (name or NCBI taxonomy ID) [default: None]
                     - min_mi_score: float - Minimum MI score (0-1) [default: 0.4]
                     - max_results: int - Maximum results [default: 200]
                     - interaction_type: str - Filter by interaction type
                     - detection_method: str - Filter by detection method
                     - negative_interactions: bool - Include negative results [default: False]
                     - format: str - Output format ("json", "tab25", "tab27") [default: "json"]

        Returns:
            AdapterResult containing interaction data with evidence and scores
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be protein identifier(s) or dict with 'identifiers'"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Parse input
        identifiers = []
        organism = kwargs.get("organism", None)
        query = None

        if isinstance(input_data, dict):
            identifiers = (input_data.get("identifiers") or
                          input_data.get("uniprot_ids") or
                          input_data.get("gene_names") or
                          [])
            query = input_data.get("query")
            organism = input_data.get("organism", organism)
        elif isinstance(input_data, list):
            identifiers = input_data
        elif isinstance(input_data, str):
            # Check if it's a raw query or identifier
            if ":" in input_data and any(op in input_data for op in ["AND", "OR"]):
                query = input_data
            else:
                identifiers = [input_data]

        # Get parameters
        min_mi_score = kwargs.get("min_mi_score", self.config.get("min_mi_score", 0.4))
        max_results = kwargs.get("max_results", self.config.get("max_results", 200))
        interaction_type = kwargs.get("interaction_type", None)
        detection_method = kwargs.get("detection_method", None)
        negative = kwargs.get("negative_interactions", self.config.get("negative_interactions", False))
        format_type = kwargs.get("format", "json")

        # Convert organism if provided
        organism_id = self._get_organism_id(organism) if organism else None

        try:
            # Build query if not provided
            if not query and identifiers:
                query = self._build_query(
                    identifiers,
                    organism_id,
                    interaction_type,
                    detection_method
                )
            elif not query:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="No identifiers or query provided"
                )

            logger.info(f"IntAct: Executing query: {query}")

            # Search for interactions
            search_result = await self._search_interactions(
                query,
                format=format_type,
                max_results=max_results,
                negative=negative
            )

            if search_result is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to search IntAct database"
                )

            # Handle JSON format
            if format_type == "json":
                interactions_raw = search_result.get("data", [])

                # Format interactions
                interactions = self._format_interaction_data(
                    interactions_raw,
                    min_mi_score
                )

                # Get unique proteins
                proteins = set()
                for interaction in interactions:
                    proteins.add(interaction["protein_a"])
                    proteins.add(interaction["protein_b"])

                # Get unique detection methods
                methods = list(set(i["detection_method"] for i in interactions))

                # Get unique interaction types
                types = list(set(i["interaction_type"] for i in interactions))

                result_data = {
                    "interactions": interactions,
                    "num_interactions": len(interactions),
                    "num_proteins": len(proteins),
                    "proteins": sorted(list(proteins)),
                    "detection_methods": sorted(methods),
                    "interaction_types": sorted(types),
                    "query": query,
                    "parameters": {
                        "organism": organism_id,
                        "min_mi_score": min_mi_score,
                        "max_results": max_results,
                        "interaction_type": interaction_type,
                        "detection_method": detection_method,
                        "negative_interactions": negative
                    }
                }

                return AdapterResult(
                    success=True,
                    data=result_data,
                    cache_hit=False,
                    metadata={
                        "source": "intact",
                        "query": query,
                        "num_interactions": len(interactions),
                        "adapter_version": self.version
                    }
                )

            # Handle TAB format (MITAB)
            else:
                raw_data = search_result.get("raw", "")
                lines = raw_data.strip().split('\n')

                # Parse MITAB lines
                interactions_raw = []
                for line in lines:
                    if line and not line.startswith('#'):
                        parsed = self._parse_mitab_line(line)
                        if parsed:
                            interactions_raw.append(parsed)

                # Format interactions
                interactions = self._format_interaction_data(
                    interactions_raw,
                    min_mi_score
                )

                result_data = {
                    "interactions": interactions,
                    "num_interactions": len(interactions),
                    "format": format_type,
                    "query": query
                }

                return AdapterResult(
                    success=True,
                    data=result_data,
                    cache_hit=False,
                    metadata={
                        "source": "intact",
                        "format": format_type,
                        "adapter_version": self.version
                    }
                )

        except Exception as e:
            logger.error(f"IntAct: Error processing request: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error processing IntAct query: {str(e)}",
                metadata={"query": query if 'query' in locals() else None}
            )
