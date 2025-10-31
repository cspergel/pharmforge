"""
Reactome Adapter - Pathway analysis and biological process enrichment
API Documentation: https://reactome.org/ContentService/
"""
from typing import Any, Dict, Optional, List, Union
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ReactomeAdapter(AdapterProtocol):
    """
    Adapter for Reactome Pathway Database
    Provides pathway analysis, enrichment, and biological process information
    """

    BASE_URL = "https://reactome.org/ContentService"

    def __init__(self):
        super().__init__(
            name="reactome",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,
                "timeout": 60,
                "species": "Homo sapiens"  # Default to human
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is valid gene list or protein IDs

        Args:
            input_data: Gene/protein identifier(s) or dict with search parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, list):
            return len(input_data) > 0 and all(isinstance(x, str) for x in input_data)
        elif isinstance(input_data, dict):
            return "genes" in input_data or "proteins" in input_data or "pathway_id" in input_data
        return False

    async def _search_pathways(self, query: str, species: str = "Homo sapiens") -> Optional[List[Dict[str, Any]]]:
        """
        Search for pathways by name or keyword

        Args:
            query: Search term
            species: Species name (default: Homo sapiens)

        Returns:
            List of matching pathways or None on error
        """
        url = f"{self.BASE_URL}/search/query"
        params = {
            "query": query,
            "species": species,
            "types": "Pathway"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.warning(f"Reactome search failed: {response.status}")
                        text = await response.text()
                        logger.debug(f"Response: {text[:500]}")
                        return None

                    data = await response.json()
                    return data.get("results", [])

        except asyncio.TimeoutError:
            logger.error(f"Reactome: Timeout searching for {query}")
            return None
        except Exception as e:
            logger.error(f"Reactome: Error searching for {query}: {e}")
            return None

    async def _get_pathway_details(self, pathway_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information about a specific pathway

        Args:
            pathway_id: Reactome pathway stable identifier (e.g., R-HSA-xxxxx)

        Returns:
            Pathway details or None on error
        """
        url = f"{self.BASE_URL}/data/query/{pathway_id}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 404:
                        logger.info(f"Pathway {pathway_id} not found")
                        return None
                    elif response.status != 200:
                        logger.warning(f"Reactome pathway detail failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        except Exception as e:
            logger.error(f"Reactome: Error fetching pathway {pathway_id}: {e}")
            return None

    async def _analyze_gene_list(self, gene_list: List[str], species: str = "Homo sapiens") -> Optional[Dict[str, Any]]:
        """
        Perform pathway enrichment analysis on a gene list

        Args:
            gene_list: List of gene symbols or identifiers
            species: Species name

        Returns:
            Enrichment results or None on error
        """
        url = f"{self.BASE_URL}/data/pathways/low/diagram/entity/list"

        # Format gene list for API
        gene_data = "\n".join(gene_list)

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    url,
                    data=gene_data,
                    headers={"Content-Type": "text/plain"},
                    timeout=timeout
                ) as response:
                    if response.status != 200:
                        logger.warning(f"Reactome analysis failed: {response.status}")
                        text = await response.text()
                        logger.debug(f"Response: {text[:500]}")
                        return None

                    data = await response.json()
                    return data

        except Exception as e:
            logger.error(f"Reactome: Error analyzing gene list: {e}")
            return None

    async def _get_pathway_hierarchy(self, pathway_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get the hierarchical structure of a pathway

        Args:
            pathway_id: Reactome pathway stable identifier

        Returns:
            List of parent/child pathway relationships or None on error
        """
        url = f"{self.BASE_URL}/data/pathway/{pathway_id}/containedEvents"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 404:
                        return []
                    elif response.status != 200:
                        logger.warning(f"Reactome hierarchy request failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        except Exception as e:
            logger.error(f"Reactome: Error fetching hierarchy for {pathway_id}: {e}")
            return None

    async def _get_pathway_participants(self, pathway_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get the participating molecules (proteins, complexes, etc.) in a pathway

        Args:
            pathway_id: Reactome pathway stable identifier

        Returns:
            List of participants or None on error
        """
        url = f"{self.BASE_URL}/data/pathway/{pathway_id}/participatingPhysicalEntities"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 404:
                        return []
                    elif response.status != 200:
                        logger.warning(f"Reactome participants request failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        except Exception as e:
            logger.error(f"Reactome: Error fetching participants for {pathway_id}: {e}")
            return None

    def _extract_pathway_info(
        self,
        pathway: Dict[str, Any],
        details: Optional[Dict[str, Any]] = None,
        hierarchy: Optional[List[Dict[str, Any]]] = None,
        participants: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, Any]:
        """
        Extract and compile pathway information

        Args:
            pathway: Basic pathway data
            details: Detailed pathway information
            hierarchy: Hierarchical structure
            participants: Participating molecules

        Returns:
            Compiled pathway information
        """
        info = {
            "pathway_id": pathway.get("stId") or pathway.get("dbId"),
            "name": pathway.get("displayName") or pathway.get("name"),
            "species": pathway.get("species", [{}])[0].get("displayName") if isinstance(pathway.get("species"), list) else None,
            "type": pathway.get("schemaClass")
        }

        # Add details if available
        if details:
            info["description"] = details.get("summation", [{}])[0].get("text") if details.get("summation") else None
            info["disease_related"] = details.get("disease")
            info["inferred"] = details.get("isInferred")

            # Extract literature references
            if details.get("literatureReference"):
                info["references"] = [
                    {
                        "pubmed_id": ref.get("pubMedIdentifier"),
                        "title": ref.get("title")
                    }
                    for ref in details.get("literatureReference", [])[:5]  # Top 5 refs
                ]

        # Add hierarchy information
        if hierarchy:
            info["sub_pathways"] = len(hierarchy)
            info["sub_pathway_list"] = [
                {
                    "id": h.get("stId"),
                    "name": h.get("displayName"),
                    "type": h.get("schemaClass")
                }
                for h in hierarchy[:10]  # Top 10 sub-pathways
            ]

        # Add participant information
        if participants:
            info["num_participants"] = len(participants)
            participant_types = {}
            for p in participants:
                p_type = p.get("schemaClass", "Unknown")
                participant_types[p_type] = participant_types.get(p_type, 0) + 1
            info["participant_types"] = participant_types

        return info

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Reactome pathway analysis

        Args:
            input_data: Gene/protein list, pathway search term, or dict with parameters
                       Examples:
                       - "TP53" (single gene)
                       - ["TP53", "BRCA1", "PTEN"] (gene list)
                       - {"genes": ["TP53", "BRCA1"], "species": "Homo sapiens"}
                       - {"pathway_id": "R-HSA-69278"} (specific pathway lookup)
            **kwargs: Additional parameters:
                     - species: str - Species name (default: Homo sapiens)
                     - enrichment: bool - Perform enrichment analysis (default: True for lists)
                     - include_hierarchy: bool - Include pathway hierarchy (default: True)
                     - include_participants: bool - Include participants (default: True)

        Returns:
            AdapterResult containing pathway analysis results
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be gene/protein identifier(s) or search parameters"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Parse input
        species = kwargs.get("species", self.config.get("species", "Homo sapiens"))

        # Handle different input types
        if isinstance(input_data, dict):
            if "pathway_id" in input_data:
                # Direct pathway lookup
                pathway_id = input_data["pathway_id"]
                pathway_details = await self._get_pathway_details(pathway_id)

                if pathway_details is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"Pathway {pathway_id} not found",
                        metadata={"pathway_id": pathway_id}
                    )

                # Get additional info
                hierarchy = None
                participants = None

                if kwargs.get("include_hierarchy", True):
                    hierarchy = await self._get_pathway_hierarchy(pathway_id)
                    await asyncio.sleep(0.3)

                if kwargs.get("include_participants", True):
                    participants = await self._get_pathway_participants(pathway_id)
                    await asyncio.sleep(0.3)

                pathway_info = self._extract_pathway_info(
                    pathway_details,
                    pathway_details,
                    hierarchy,
                    participants
                )

                return AdapterResult(
                    success=True,
                    data={"pathway": pathway_info},
                    cache_hit=False,
                    metadata={
                        "pathway_id": pathway_id,
                        "source": "reactome",
                        "adapter_version": self.version
                    }
                )

            elif "genes" in input_data or "proteins" in input_data:
                gene_list = input_data.get("genes") or input_data.get("proteins")
                species = input_data.get("species", species)
            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Invalid dict input: must contain 'pathway_id', 'genes', or 'proteins'"
                )

        elif isinstance(input_data, str):
            # Single gene or search term
            # Try as search first
            pathways = await self._search_pathways(input_data, species)

            if pathways is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to search Reactome",
                    metadata={"query": input_data}
                )

            if not pathways:
                return AdapterResult(
                    success=True,
                    data={
                        "pathways": [],
                        "total_found": 0
                    },
                    metadata={"query": input_data, "source": "reactome"}
                )

            # Return search results
            pathway_list = []
            for pathway in pathways[:10]:  # Top 10 results
                pathway_info = self._extract_pathway_info(pathway)
                pathway_list.append(pathway_info)

            return AdapterResult(
                success=True,
                data={
                    "pathways": pathway_list,
                    "total_found": len(pathways)
                },
                cache_hit=False,
                metadata={
                    "query": input_data,
                    "source": "reactome",
                    "adapter_version": self.version
                }
            )

        elif isinstance(input_data, list):
            gene_list = input_data
        else:
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input type"
            )

        # Perform enrichment analysis on gene list
        if kwargs.get("enrichment", True):
            enrichment_results = await self._analyze_gene_list(gene_list, species)

            if enrichment_results is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to perform enrichment analysis",
                    metadata={"genes": gene_list, "species": species}
                )

            # Parse enrichment results
            enriched_pathways = []
            for pathway_id in enrichment_results.get("pathways", [])[:20]:  # Top 20
                pathway_details = await self._get_pathway_details(pathway_id)
                if pathway_details:
                    pathway_info = self._extract_pathway_info(pathway_details, pathway_details)
                    enriched_pathways.append(pathway_info)
                await asyncio.sleep(0.3)  # Rate limit

            return AdapterResult(
                success=True,
                data={
                    "enriched_pathways": enriched_pathways,
                    "input_genes": gene_list,
                    "species": species,
                    "total_pathways": len(enrichment_results.get("pathways", []))
                },
                cache_hit=False,
                metadata={
                    "analysis_type": "enrichment",
                    "source": "reactome",
                    "adapter_version": self.version
                }
            )

        # If no enrichment, just search for each gene
        all_pathways = []
        for gene in gene_list[:5]:  # Limit to 5 genes for search
            pathways = await self._search_pathways(gene, species)
            if pathways:
                for pathway in pathways[:3]:  # Top 3 per gene
                    pathway_info = self._extract_pathway_info(pathway)
                    pathway_info["matched_gene"] = gene
                    all_pathways.append(pathway_info)
            await asyncio.sleep(0.3)

        return AdapterResult(
            success=True,
            data={
                "pathways": all_pathways,
                "input_genes": gene_list,
                "species": species
            },
            cache_hit=False,
            metadata={
                "analysis_type": "search",
                "source": "reactome",
                "adapter_version": self.version
            }
        )
