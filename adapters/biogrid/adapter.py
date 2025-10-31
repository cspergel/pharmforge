"""
BioGRID Adapter - Fetches protein-protein interactions from BioGRID REST API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class BioGRIDAdapter(AdapterProtocol):
    """
    Adapter for BioGRID REST API
    Fetches protein-protein interactions, genetic interactions, and related data

    BioGRID (Biological General Repository for Interaction Datasets) is a
    comprehensive database of protein, genetic, and chemical interactions.
    """

    BASE_URL = "https://webservice.thebiogrid.org"

    # Common organism taxonomy IDs
    ORGANISMS = {
        "human": "9606",
        "mouse": "10090",
        "rat": "10116",
        "yeast": "559292",
        "fly": "7227",
        "worm": "6239",
        "zebrafish": "7955",
        "arabidopsis": "3702",
        "ecoli": "83333"
    }

    def __init__(self, access_key: Optional[str] = None):
        """
        Initialize the BioGRID adapter

        Args:
            access_key: BioGRID API access key (required for API access)
                       Get one at: https://webservice.thebiogrid.org/
        """
        super().__init__(
            name="biogrid",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.2,  # 5 requests/second (be nice to BioGRID)
                "timeout": 120,
                "max_results": 10000,  # BioGRID max per query
                "access_key": access_key
            }
        )
        self.version = "1.0.0"
        self.access_key = access_key

        if not self.access_key:
            logger.warning("BioGRID: No access key provided. Get one at https://webservice.thebiogrid.org/")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid gene name, gene ID, or list of genes

        Args:
            input_data: String or list of gene identifiers

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, list):
            return len(input_data) > 0 and all(isinstance(x, str) for x in input_data)
        return False

    def _get_organism_id(self, organism: str) -> str:
        """
        Convert organism name to taxonomy ID

        Args:
            organism: Organism name (e.g., "human", "mouse") or taxonomy ID

        Returns:
            Taxonomy ID string
        """
        # If already a numeric ID, return as is
        if organism.isdigit():
            return organism

        # Look up in common organisms
        organism_lower = organism.lower()
        return self.ORGANISMS.get(organism_lower, organism)

    async def _fetch_interactions(
        self,
        gene_list: List[str],
        organism: Optional[str] = None,
        evidence_list: Optional[List[str]] = None,
        include_interactors: bool = False,
        exclude_interspecies: bool = True,
        max_results: Optional[int] = None,
        format: str = "json"
    ) -> Optional[Dict[str, Any]]:
        """
        Fetch interactions from BioGRID API

        Args:
            gene_list: List of gene names or IDs to query
            organism: Organism name or taxonomy ID (e.g., "human", "9606")
            evidence_list: List of evidence types to include (e.g., ["Affinity Capture-MS"])
            include_interactors: Include first-order interactors
            exclude_interspecies: Exclude cross-species interactions
            max_results: Maximum number of results (default: 10000)
            format: Response format ("json", "tab1", "tab2", "count")

        Returns:
            Dictionary containing interaction data or None on error
        """
        if not self.access_key:
            logger.error("BioGRID: Access key required. Get one at https://webservice.thebiogrid.org/")
            return None

        url = f"{self.BASE_URL}/interactions/"

        # Prepare parameters
        params = {
            "accesskey": self.access_key,
            "format": format,
            "geneList": "|".join(gene_list),  # Pipe-separated gene list
            "searchNames": "true",  # Search official gene symbols
            "includeInteractors": "true" if include_interactors else "false",
            "interSpeciesExcluded": "true" if exclude_interspecies else "false",
            "selfInteractionsExcluded": "false",
            "max": max_results or self.config.get("max_results", 10000)
        }

        # Add organism filter if specified
        if organism:
            tax_id = self._get_organism_id(organism)
            params["taxId"] = tax_id

        # Add evidence filter if specified
        if evidence_list:
            params["evidenceList"] = "|".join(evidence_list)
            params["includeEvidence"] = "true"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 120))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        if format == "json":
                            data = await response.json()
                            return data
                        else:
                            text = await response.text()
                            return {"data": text}
                    else:
                        error_text = await response.text()
                        logger.error(f"BioGRID API error {response.status}: {error_text}")
                        return None

        except asyncio.TimeoutError:
            logger.error(f"BioGRID: Timeout fetching interactions for {gene_list}")
            return None
        except Exception as e:
            logger.error(f"BioGRID: Error fetching interactions: {e}")
            return None

    async def _get_organisms(self) -> Optional[List[Dict[str, Any]]]:
        """
        Get list of available organisms in BioGRID

        Returns:
            List of organism dictionaries or None on error
        """
        if not self.access_key:
            logger.error("BioGRID: Access key required")
            return None

        url = f"{self.BASE_URL}/organisms/"
        params = {"accesskey": self.access_key, "format": "json"}

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 120))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        return await response.json()
                    else:
                        logger.error(f"BioGRID organisms error {response.status}")
                        return None
        except Exception as e:
            logger.error(f"BioGRID: Error fetching organisms: {e}")
            return None

    async def _get_evidence_types(self) -> Optional[List[str]]:
        """
        Get list of experimental evidence types available in BioGRID

        Returns:
            List of evidence type strings or None on error
        """
        if not self.access_key:
            logger.error("BioGRID: Access key required")
            return None

        url = f"{self.BASE_URL}/evidence/"
        params = {"accesskey": self.access_key, "format": "json"}

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 120))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        return await response.json()
                    else:
                        logger.error(f"BioGRID evidence types error {response.status}")
                        return None
        except Exception as e:
            logger.error(f"BioGRID: Error fetching evidence types: {e}")
            return None

    def _parse_interactions(self, raw_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse and summarize interaction data from BioGRID response

        Args:
            raw_data: Raw JSON response from BioGRID API

        Returns:
            Summarized interaction data
        """
        if not raw_data or "error" in raw_data:
            return {
                "num_interactions": 0,
                "interactions": [],
                "summary": {}
            }

        interactions = []
        experimental_systems = set()
        interaction_types = set()
        organisms = set()
        interactor_genes = set()
        pubmed_ids = set()

        # Parse each interaction record
        for interaction_id, interaction in raw_data.items():
            if not isinstance(interaction, dict):
                continue

            # Extract key information
            gene_a = interaction.get("OFFICIAL_SYMBOL_A", "")
            gene_b = interaction.get("OFFICIAL_SYMBOL_B", "")
            entrez_a = interaction.get("ENTREZ_GENE_A", "")
            entrez_b = interaction.get("ENTREZ_GENE_B", "")
            organism_a = interaction.get("ORGANISM_A_ID", "")
            organism_b = interaction.get("ORGANISM_B_ID", "")
            experimental_system = interaction.get("EXPERIMENTAL_SYSTEM", "")
            experimental_system_type = interaction.get("EXPERIMENTAL_SYSTEM_TYPE", "")
            pubmed_id = interaction.get("PUBMED_ID", "")
            interaction_score = interaction.get("SCORE", "")

            # Collect summary statistics
            if experimental_system:
                experimental_systems.add(experimental_system)
            if experimental_system_type:
                interaction_types.add(experimental_system_type)
            if organism_a:
                organisms.add(organism_a)
            if organism_b:
                organisms.add(organism_b)
            if gene_a:
                interactor_genes.add(gene_a)
            if gene_b:
                interactor_genes.add(gene_b)
            if pubmed_id:
                pubmed_ids.add(pubmed_id)

            # Create simplified interaction record
            interactions.append({
                "interaction_id": interaction_id,
                "gene_a": gene_a,
                "gene_b": gene_b,
                "entrez_gene_a": entrez_a,
                "entrez_gene_b": entrez_b,
                "organism_a": organism_a,
                "organism_b": organism_b,
                "experimental_system": experimental_system,
                "experimental_system_type": experimental_system_type,
                "pubmed_id": pubmed_id,
                "score": interaction_score,
                "uniprot_a": interaction.get("SWISS_PROT_ACCESSIONS_A", ""),
                "uniprot_b": interaction.get("SWISS_PROT_ACCESSIONS_B", ""),
            })

        # Create summary
        summary = {
            "num_interactions": len(interactions),
            "num_unique_genes": len(interactor_genes),
            "num_publications": len(pubmed_ids),
            "experimental_systems": sorted(list(experimental_systems)),
            "interaction_types": sorted(list(interaction_types)),
            "organisms": sorted(list(organisms)),
            "unique_genes": sorted(list(interactor_genes))[:20]  # Top 20 genes
        }

        return {
            "num_interactions": len(interactions),
            "interactions": interactions,
            "summary": summary
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute BioGRID interaction query

        Args:
            input_data: Gene name(s) as string or list of strings
            **kwargs: Additional parameters
                - organism: Organism name or taxonomy ID (e.g., "human", "9606")
                - evidence_types: List of evidence types to filter
                - include_interactors: Include first-order interactors (default: False)
                - exclude_interspecies: Exclude cross-species interactions (default: True)
                - max_results: Maximum results to return (default: 10000)
                - query_type: Type of query ("interactions", "organisms", "evidence")

        Returns:
            AdapterResult containing interaction data
        """
        # Special queries
        query_type = kwargs.get("query_type", "interactions")

        if query_type == "organisms":
            organisms = await self._get_organisms()
            if organisms is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to fetch organism list from BioGRID",
                    metadata={"source": "biogrid", "query_type": "organisms"}
                )
            return AdapterResult(
                success=True,
                data={"organisms": organisms},
                cache_hit=False,
                metadata={"source": "biogrid", "query_type": "organisms"}
            )

        if query_type == "evidence":
            evidence_types = await self._get_evidence_types()
            if evidence_types is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to fetch evidence types from BioGRID",
                    metadata={"source": "biogrid", "query_type": "evidence"}
                )
            return AdapterResult(
                success=True,
                data={"evidence_types": evidence_types},
                cache_hit=False,
                metadata={"source": "biogrid", "query_type": "evidence"}
            )

        # Validate input for interaction queries
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a non-empty string or list of gene identifiers"
            )

        # Convert input to list if single string
        if isinstance(input_data, str):
            gene_list = [input_data]
        else:
            gene_list = input_data

        # Extract parameters
        organism = kwargs.get("organism", None)
        evidence_types = kwargs.get("evidence_types", None)
        include_interactors = kwargs.get("include_interactors", False)
        exclude_interspecies = kwargs.get("exclude_interspecies", True)
        max_results = kwargs.get("max_results", None)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.2)
        await asyncio.sleep(rate_delay)

        # Fetch interactions
        raw_data = await self._fetch_interactions(
            gene_list=gene_list,
            organism=organism,
            evidence_list=evidence_types,
            include_interactors=include_interactors,
            exclude_interspecies=exclude_interspecies,
            max_results=max_results
        )

        if raw_data is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to fetch interactions from BioGRID API",
                metadata={
                    "source": "biogrid",
                    "genes": gene_list,
                    "organism": organism
                }
            )

        # Parse and summarize results
        parsed_data = self._parse_interactions(raw_data)

        return AdapterResult(
            success=True,
            data=parsed_data,
            cache_hit=False,
            metadata={
                "source": "biogrid",
                "adapter_version": self.version,
                "genes_queried": gene_list,
                "organism": organism,
                "num_interactions_found": parsed_data["num_interactions"]
            }
        )
