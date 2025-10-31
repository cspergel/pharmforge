"""
UniProt Adapter - Fetches protein sequences and annotations from UniProt REST API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class UniProtAdapter(AdapterProtocol):
    """
    Adapter for UniProt REST API
    Fetches protein sequences, annotations, functions, structures, and post-translational modifications
    """

    BASE_URL = "https://rest.uniprot.org"

    def __init__(self):
        super().__init__(
            name="uniprot",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.1,  # 10 requests/second
                "timeout": 60,
                "format": "json"
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid protein ID or gene name

        Args:
            input_data: String containing UniProt ID, gene name, or accession

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        return True

    async def _search_protein(self, query: str, max_results: int = 10) -> Optional[List[str]]:
        """
        Search for protein IDs matching a query

        Args:
            query: Search query (gene name, protein name, etc.)
            max_results: Maximum number of results

        Returns:
            List of UniProt accession IDs or None
        """
        url = f"{self.BASE_URL}/uniprotkb/search"

        params = {
            "query": query,
            "format": "json",
            "size": max_results,
            "fields": "accession,id,gene_names,protein_name"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        results = data.get("results", [])
                        return [r.get("primaryAccession") for r in results if r.get("primaryAccession")]
                    else:
                        error_text = await response.text()
                        logger.error(f"UniProt search error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"UniProt: Timeout searching for {query}")
            return None
        except Exception as e:
            logger.error(f"UniProt: Error searching for {query}: {e}")
            return None

    async def _get_protein_info(self, accession: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed protein information

        Args:
            accession: UniProt accession ID

        Returns:
            Dictionary of protein information or None
        """
        url = f"{self.BASE_URL}/uniprotkb/{accession}.json"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        return await response.json()
                    elif response.status == 404:
                        logger.warning(f"UniProt: Protein not found: {accession}")
                        return None
                    else:
                        error_text = await response.text()
                        logger.error(f"UniProt API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"UniProt: Timeout fetching {accession}")
            return None
        except Exception as e:
            logger.error(f"UniProt: Error fetching {accession}: {e}")
            return None

    async def _get_protein_sequence(self, accession: str, format: str = "fasta") -> Optional[str]:
        """
        Get protein sequence in specified format

        Args:
            accession: UniProt accession ID
            format: Sequence format ("fasta", "json", etc.)

        Returns:
            Sequence string or None
        """
        url = f"{self.BASE_URL}/uniprotkb/{accession}.{format}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        return await response.text()
                    else:
                        logger.error(f"UniProt sequence error {response.status}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"UniProt: Timeout fetching sequence for {accession}")
            return None
        except Exception as e:
            logger.error(f"UniProt: Error fetching sequence for {accession}: {e}")
            return None

    def _extract_protein_info(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract and format protein information from UniProt JSON

        Args:
            data: Raw UniProt JSON response

        Returns:
            Formatted protein information
        """
        # Basic identifiers
        accession = data.get("primaryAccession")
        entry_name = data.get("uniProtkbId")

        # Gene names
        genes = data.get("genes", [])
        gene_names = []
        if genes:
            for gene in genes:
                gene_name = gene.get("geneName", {}).get("value")
                if gene_name:
                    gene_names.append(gene_name)

        # Protein names
        protein_name = None
        protein_description = data.get("proteinDescription", {})
        recommended_name = protein_description.get("recommendedName")
        if recommended_name:
            full_name = recommended_name.get("fullName", {}).get("value")
            if full_name:
                protein_name = full_name

        # Organism
        organism = None
        organism_data = data.get("organism", {})
        scientific_name = organism_data.get("scientificName")
        if scientific_name:
            organism = scientific_name

        # Sequence
        sequence_data = data.get("sequence", {})
        sequence = sequence_data.get("value", "")
        sequence_length = sequence_data.get("length", 0)
        molecular_weight = sequence_data.get("molWeight", 0)

        # Function/Comments
        comments = data.get("comments", [])
        functions = []
        subcellular_locations = []
        disease_associations = []
        ptms = []  # Post-translational modifications

        for comment in comments:
            comment_type = comment.get("commentType")

            if comment_type == "FUNCTION":
                texts = comment.get("texts", [])
                for text in texts:
                    functions.append(text.get("value", ""))

            elif comment_type == "SUBCELLULAR LOCATION":
                locations = comment.get("subcellularLocations", [])
                for loc in locations:
                    location_value = loc.get("location", {}).get("value")
                    if location_value:
                        subcellular_locations.append(location_value)

            elif comment_type == "DISEASE":
                disease = comment.get("disease", {})
                disease_name = disease.get("diseaseId")
                if disease_name:
                    disease_associations.append(disease_name)

            elif comment_type == "PTM":
                texts = comment.get("texts", [])
                for text in texts:
                    ptms.append(text.get("value", ""))

        # Features (domains, modifications, etc.)
        features = data.get("features", [])
        domains = []
        active_sites = []
        binding_sites = []
        modifications = []

        for feature in features:
            feature_type = feature.get("type")
            description = feature.get("description")

            if feature_type == "Domain":
                if description:
                    domains.append(description)
            elif feature_type == "Active site":
                if description:
                    active_sites.append(description)
            elif feature_type == "Binding site":
                if description:
                    binding_sites.append(description)
            elif feature_type in ["Modified residue", "Glycosylation", "Disulfide bond"]:
                if description:
                    modifications.append(f"{feature_type}: {description}")

        # Cross-references
        references = data.get("uniProtKBCrossReferences", [])
        pdb_ids = []
        alphafold_ids = []

        for ref in references:
            database = ref.get("database")
            ref_id = ref.get("id")

            if database == "PDB" and ref_id:
                pdb_ids.append(ref_id)
            elif database == "AlphaFoldDB" and ref_id:
                alphafold_ids.append(ref_id)

        # Keywords
        keywords_data = data.get("keywords", [])
        keywords = [kw.get("name") for kw in keywords_data if kw.get("name")]

        return {
            "accession": accession,
            "entry_name": entry_name,
            "gene_names": gene_names,
            "protein_name": protein_name,
            "organism": organism,
            "sequence": sequence,
            "sequence_length": sequence_length,
            "molecular_weight": molecular_weight,
            "function": functions,
            "subcellular_location": subcellular_locations,
            "disease_associations": disease_associations,
            "post_translational_modifications": ptms,
            "domains": domains,
            "active_sites": active_sites,
            "binding_sites": binding_sites,
            "modifications": modifications,
            "pdb_structures": pdb_ids,
            "alphafold_structures": alphafold_ids,
            "keywords": keywords
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the UniProt query

        Args:
            input_data: UniProt accession, gene name, or protein name
            **kwargs: Additional parameters
                - include_sequence: Include full sequence (default: True)
                - search_first: If True and input is not accession, search first (default: True)

        Returns:
            AdapterResult containing UniProt protein information
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a non-empty string"
            )

        include_sequence = kwargs.get("include_sequence", True)
        search_first = kwargs.get("search_first", True)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.1)
        await asyncio.sleep(rate_delay)

        accession = input_data

        # If input doesn't look like a UniProt accession, search for it
        if search_first and not (len(input_data) == 6 or len(input_data) == 10):
            logger.info(f"UniProt: Searching for protein: {input_data}")
            search_results = await self._search_protein(input_data, max_results=1)

            if not search_results:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Protein not found: {input_data}",
                    metadata={"source": "uniprot", "query": input_data}
                )

            accession = search_results[0]
            logger.info(f"UniProt: Found accession: {accession}")

        # Get protein information
        try:
            protein_data = await self._get_protein_info(accession)

            if not protein_data:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Failed to fetch protein data for: {accession}",
                    metadata={"source": "uniprot", "accession": accession}
                )

            # Extract and format information
            result_data = self._extract_protein_info(protein_data)

            # Optionally remove sequence to reduce size
            if not include_sequence:
                result_data["sequence"] = None

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "uniprot",
                    "accession": accession,
                    "original_query": input_data,
                    "adapter_version": self.version
                }
            )

        except Exception as e:
            logger.error(f"UniProt: Error processing results: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error processing UniProt data: {str(e)}",
                metadata={"source": "uniprot", "accession": accession}
            )
