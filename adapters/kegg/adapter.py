"""
KEGG Adapter - Kyoto Encyclopedia of Genes and Genomes
API Documentation: https://www.kegg.jp/kegg/rest/keggapi.html

Provides access to KEGG pathway, gene, compound, and disease databases
using the FREE KEGG REST API.
"""
from typing import Any, Dict, Optional, List, Union
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class KEGGAdapter(AdapterProtocol):
    """
    Adapter for KEGG (Kyoto Encyclopedia of Genes and Genomes) Database

    Supports queries for:
    - Pathways: Metabolic and signaling pathway maps
    - Genes: Gene information across organisms
    - Compounds: Chemical structures and properties
    - Diseases: Disease information and associated genes/pathways
    - Drugs: Pharmaceutical compounds

    All queries use the FREE KEGG REST API.
    """

    BASE_URL = "https://rest.kegg.jp"

    def __init__(self):
        super().__init__(
            name="kegg",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.3,  # Be respectful to free API
                "timeout": 30,
                "default_organism": "hsa",  # Human (Homo sapiens)
                "max_entries": 10  # KEGG API limit for some operations
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid KEGG query

        Args:
            input_data: KEGG identifier, search term, or dict with query parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data.strip()) > 0
        elif isinstance(input_data, list):
            return len(input_data) > 0 and all(isinstance(x, str) for x in input_data)
        elif isinstance(input_data, dict):
            required_keys = ["operation"]
            return any(key in input_data for key in required_keys) or \
                   any(key in input_data for key in ["pathway_id", "gene_id", "compound_id", "disease_id", "query"])
        return False

    async def _make_request(self, endpoint: str) -> Optional[str]:
        """
        Make a request to KEGG REST API

        Args:
            endpoint: API endpoint path (e.g., "/list/pathway/hsa")

        Returns:
            Response text or None on error
        """
        url = f"{self.BASE_URL}{endpoint}"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 404:
                        logger.info(f"KEGG: Resource not found: {endpoint}")
                        return None
                    elif response.status != 200:
                        logger.warning(f"KEGG request failed: {response.status} for {endpoint}")
                        return None

                    text = await response.text()
                    return text

        except asyncio.TimeoutError:
            logger.error(f"KEGG: Timeout for {endpoint}")
            return None
        except Exception as e:
            logger.error(f"KEGG: Error fetching {endpoint}: {e}")
            return None

    async def _get_info(self, database: str) -> Optional[Dict[str, Any]]:
        """
        Get database information and statistics

        Args:
            database: Database name (e.g., "pathway", "hsa", "compound")

        Returns:
            Database info or None on error
        """
        endpoint = f"/info/{database}"
        text = await self._make_request(endpoint)

        if text is None:
            return None

        # Parse info response
        info = {"database": database, "raw_info": text}
        lines = text.strip().split('\n')

        for line in lines:
            if ':' in line:
                key, value = line.split(':', 1)
                info[key.strip().lower().replace(' ', '_')] = value.strip()

        return info

    async def _list_entries(self, database: str, organism: Optional[str] = None) -> Optional[List[Dict[str, str]]]:
        """
        List entries in a database

        Args:
            database: Database name (e.g., "pathway", "compound", "disease")
            organism: Optional organism code (e.g., "hsa" for human)

        Returns:
            List of entries or None on error
        """
        if organism and database == "pathway":
            endpoint = f"/list/pathway/{organism}"
        else:
            endpoint = f"/list/{database}"

        text = await self._make_request(endpoint)

        if text is None:
            return None

        # Parse tab-delimited results
        entries = []
        lines = text.strip().split('\n')

        for line in lines:
            if '\t' in line:
                entry_id, name = line.split('\t', 1)
                entries.append({
                    "id": entry_id.strip(),
                    "name": name.strip()
                })

        return entries

    async def _find_entries(self, database: str, query: str, option: Optional[str] = None) -> Optional[List[Dict[str, str]]]:
        """
        Search for entries by keyword or property

        Args:
            database: Database to search (e.g., "genes", "compound", "disease")
            query: Search query
            option: Optional search option (e.g., "formula", "exact_mass")

        Returns:
            List of matching entries or None on error
        """
        query_encoded = quote(query)

        if option:
            endpoint = f"/find/{database}/{query_encoded}/{option}"
        else:
            endpoint = f"/find/{database}/{query_encoded}"

        text = await self._make_request(endpoint)

        if text is None:
            return None

        # Parse results
        entries = []
        lines = text.strip().split('\n')

        for line in lines:
            if '\t' in line:
                entry_id, description = line.split('\t', 1)
                entries.append({
                    "id": entry_id.strip(),
                    "description": description.strip()
                })

        return entries

    async def _get_entry(self, entry_ids: Union[str, List[str]], option: Optional[str] = None) -> Optional[str]:
        """
        Get full entry information

        Args:
            entry_ids: Single entry ID or list of entry IDs (max 10)
            option: Optional format (e.g., "json", "image", "kgml")

        Returns:
            Entry data or None on error
        """
        if isinstance(entry_ids, list):
            entry_ids = "+".join(entry_ids[:self.config.get("max_entries", 10)])

        if option:
            endpoint = f"/get/{entry_ids}/{option}"
        else:
            endpoint = f"/get/{entry_ids}"

        return await self._make_request(endpoint)

    async def _link_databases(self, target_db: str, source_db: Optional[str] = None,
                             entry_ids: Optional[Union[str, List[str]]] = None) -> Optional[List[Dict[str, str]]]:
        """
        Find links between databases

        Args:
            target_db: Target database
            source_db: Source database (if querying all links)
            entry_ids: Specific entry IDs to query

        Returns:
            List of links or None on error
        """
        if entry_ids:
            if isinstance(entry_ids, list):
                entry_ids = "+".join(entry_ids[:self.config.get("max_entries", 10)])
            endpoint = f"/link/{target_db}/{entry_ids}"
        elif source_db:
            endpoint = f"/link/{target_db}/{source_db}"
        else:
            logger.error("KEGG link: Must provide either source_db or entry_ids")
            return None

        text = await self._make_request(endpoint)

        if text is None:
            return None

        # Parse tab-delimited results
        links = []
        lines = text.strip().split('\n')

        for line in lines:
            if '\t' in line:
                source, target = line.split('\t', 1)
                links.append({
                    "source": source.strip(),
                    "target": target.strip()
                })

        return links

    async def _convert_ids(self, target_db: str, source_db: Optional[str] = None,
                          entry_ids: Optional[Union[str, List[str]]] = None) -> Optional[Dict[str, str]]:
        """
        Convert between KEGG and external database IDs

        Args:
            target_db: Target database (e.g., "ncbi-geneid", "uniprot")
            source_db: Source database
            entry_ids: Specific entry IDs to convert

        Returns:
            Mapping of source to target IDs or None on error
        """
        if entry_ids:
            if isinstance(entry_ids, list):
                entry_ids = "+".join(entry_ids[:self.config.get("max_entries", 10)])
            endpoint = f"/conv/{target_db}/{entry_ids}"
        elif source_db:
            endpoint = f"/conv/{target_db}/{source_db}"
        else:
            logger.error("KEGG conv: Must provide either source_db or entry_ids")
            return None

        text = await self._make_request(endpoint)

        if text is None:
            return None

        # Parse tab-delimited results
        conversions = {}
        lines = text.strip().split('\n')

        for line in lines:
            if '\t' in line:
                source, target = line.split('\t', 1)
                conversions[source.strip()] = target.strip()

        return conversions

    def _parse_flat_file(self, text: str) -> Dict[str, Any]:
        """
        Parse KEGG flat file format

        Args:
            text: Flat file text

        Returns:
            Parsed data as dictionary
        """
        result = {}
        current_key = None
        current_value = []

        lines = text.split('\n')
        for line in lines:
            if not line.strip():
                continue

            # Check if line starts with a key (not indented)
            if line and not line[0].isspace():
                # Save previous key-value
                if current_key:
                    result[current_key] = '\n'.join(current_value).strip()

                # Parse new key
                parts = line.split(None, 1)
                if len(parts) >= 2:
                    current_key, value = parts
                    current_value = [value]
                elif len(parts) == 1:
                    # Key with no immediate value
                    current_key = parts[0].strip()
                    current_value = []
                else:
                    # Empty line or just whitespace
                    continue
            else:
                # Continuation of previous value
                if current_key and line.strip():
                    current_value.append(line.strip())

        # Save last key-value
        if current_key:
            result[current_key] = '\n'.join(current_value).strip()

        return result

    async def query_pathway(self, pathway_id: Optional[str] = None, organism: Optional[str] = None,
                           search_term: Optional[str] = None) -> Dict[str, Any]:
        """
        Query pathway information

        Args:
            pathway_id: Specific pathway ID (e.g., "hsa04010", "map00010")
            organism: Organism code (e.g., "hsa" for human)
            search_term: Search term to find pathways

        Returns:
            Pathway information
        """
        if pathway_id:
            # Get specific pathway
            text = await self._get_entry(pathway_id)
            if text:
                parsed = self._parse_flat_file(text)
                return {
                    "pathway_id": pathway_id,
                    "data": parsed,
                    "raw": text
                }
            return {"error": f"Pathway {pathway_id} not found"}

        elif search_term:
            # Search pathways
            results = await self._find_entries("pathway", search_term)
            return {
                "query": search_term,
                "results": results or [],
                "count": len(results) if results else 0
            }

        elif organism:
            # List pathways for organism
            results = await self._list_entries("pathway", organism)
            return {
                "organism": organism,
                "pathways": results or [],
                "count": len(results) if results else 0
            }

        else:
            # List all pathway maps
            results = await self._list_entries("pathway")
            return {
                "pathways": results[:50] if results else [],  # Limit to 50
                "count": len(results) if results else 0
            }

    async def query_gene(self, gene_id: Optional[str] = None, organism: Optional[str] = None,
                        search_term: Optional[str] = None) -> Dict[str, Any]:
        """
        Query gene information

        Args:
            gene_id: Specific gene ID (e.g., "hsa:10458")
            organism: Organism code (e.g., "hsa")
            search_term: Search term to find genes

        Returns:
            Gene information
        """
        if gene_id:
            # Get specific gene
            text = await self._get_entry(gene_id)
            if text:
                parsed = self._parse_flat_file(text)

                # Get associated pathways
                pathways = await self._link_databases("pathway", entry_ids=gene_id)

                return {
                    "gene_id": gene_id,
                    "data": parsed,
                    "pathways": pathways or [],
                    "raw": text
                }
            return {"error": f"Gene {gene_id} not found"}

        elif search_term:
            # Search genes
            results = await self._find_entries("genes", search_term)
            return {
                "query": search_term,
                "results": results or [],
                "count": len(results) if results else 0
            }

        elif organism:
            # List genes for organism (limited to first batch)
            results = await self._list_entries(organism)
            return {
                "organism": organism,
                "genes": results[:100] if results else [],  # Limit to 100
                "count": len(results) if results else 0
            }

        else:
            return {"error": "Must provide gene_id, organism, or search_term"}

    async def query_compound(self, compound_id: Optional[str] = None, search_term: Optional[str] = None,
                            formula: Optional[str] = None) -> Dict[str, Any]:
        """
        Query compound information

        Args:
            compound_id: Specific compound ID (e.g., "C00031")
            search_term: Search term to find compounds
            formula: Chemical formula to search

        Returns:
            Compound information
        """
        if compound_id:
            # Get specific compound
            text = await self._get_entry(compound_id)
            if text:
                parsed = self._parse_flat_file(text)

                # Get associated pathways
                pathways = await self._link_databases("pathway", entry_ids=compound_id)

                return {
                    "compound_id": compound_id,
                    "data": parsed,
                    "pathways": pathways or [],
                    "raw": text
                }
            return {"error": f"Compound {compound_id} not found"}

        elif formula:
            # Search by formula
            results = await self._find_entries("compound", formula, "formula")
            return {
                "formula": formula,
                "results": results or [],
                "count": len(results) if results else 0
            }

        elif search_term:
            # Search compounds
            results = await self._find_entries("compound", search_term)
            return {
                "query": search_term,
                "results": results or [],
                "count": len(results) if results else 0
            }

        else:
            # List compounds (limited)
            results = await self._list_entries("compound")
            return {
                "compounds": results[:100] if results else [],  # Limit to 100
                "count": len(results) if results else 0
            }

    async def query_disease(self, disease_id: Optional[str] = None, search_term: Optional[str] = None) -> Dict[str, Any]:
        """
        Query disease information

        Args:
            disease_id: Specific disease ID (e.g., "H00001")
            search_term: Search term to find diseases

        Returns:
            Disease information
        """
        if disease_id:
            # Get specific disease
            text = await self._get_entry(disease_id)
            if text:
                parsed = self._parse_flat_file(text)

                # Get associated genes
                genes = await self._link_databases("genes", entry_ids=disease_id)

                # Get associated pathways
                pathways = await self._link_databases("pathway", entry_ids=disease_id)

                return {
                    "disease_id": disease_id,
                    "data": parsed,
                    "genes": genes or [],
                    "pathways": pathways or [],
                    "raw": text
                }
            return {"error": f"Disease {disease_id} not found"}

        elif search_term:
            # Search diseases
            results = await self._find_entries("disease", search_term)
            return {
                "query": search_term,
                "results": results or [],
                "count": len(results) if results else 0
            }

        else:
            # List diseases
            results = await self._list_entries("disease")
            return {
                "diseases": results[:100] if results else [],  # Limit to 100
                "count": len(results) if results else 0
            }

    async def query_drug(self, drug_id: Optional[str] = None, search_term: Optional[str] = None) -> Dict[str, Any]:
        """
        Query drug information

        Args:
            drug_id: Specific drug ID (e.g., "D00001")
            search_term: Search term to find drugs

        Returns:
            Drug information
        """
        if drug_id:
            # Get specific drug
            text = await self._get_entry(drug_id)
            if text:
                parsed = self._parse_flat_file(text)

                # Get associated targets
                targets = await self._link_databases("genes", entry_ids=drug_id)

                return {
                    "drug_id": drug_id,
                    "data": parsed,
                    "targets": targets or [],
                    "raw": text
                }
            return {"error": f"Drug {drug_id} not found"}

        elif search_term:
            # Search drugs
            results = await self._find_entries("drug", search_term)
            return {
                "query": search_term,
                "results": results or [],
                "count": len(results) if results else 0
            }

        else:
            # List drugs (limited)
            results = await self._list_entries("drug")
            return {
                "drugs": results[:100] if results else [],  # Limit to 100
                "count": len(results) if results else 0
            }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute KEGG query

        Args:
            input_data: Query string, ID, or dict with query parameters
                       Examples:
                       - "glucose" (search term)
                       - "C00031" (compound ID)
                       - "hsa:10458" (gene ID)
                       - {"operation": "query_pathway", "pathway_id": "hsa04010"}
                       - {"operation": "query_compound", "search_term": "aspirin"}
                       - {"operation": "query_disease", "disease_id": "H00001"}
            **kwargs: Additional parameters:
                     - operation: str - Operation to perform (query_pathway, query_gene, query_compound, query_disease)
                     - organism: str - Organism code (default: hsa)

        Returns:
            AdapterResult containing query results
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be search term, ID, or dict with query parameters"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.3)
        await asyncio.sleep(rate_delay)

        try:
            # Parse input and determine operation
            if isinstance(input_data, dict):
                operation = input_data.get("operation", kwargs.get("operation"))

                if operation == "query_pathway":
                    result = await self.query_pathway(
                        pathway_id=input_data.get("pathway_id"),
                        organism=input_data.get("organism", self.config.get("default_organism")),
                        search_term=input_data.get("search_term")
                    )
                elif operation == "query_gene":
                    result = await self.query_gene(
                        gene_id=input_data.get("gene_id"),
                        organism=input_data.get("organism", self.config.get("default_organism")),
                        search_term=input_data.get("search_term")
                    )
                elif operation == "query_compound":
                    result = await self.query_compound(
                        compound_id=input_data.get("compound_id"),
                        search_term=input_data.get("search_term"),
                        formula=input_data.get("formula")
                    )
                elif operation == "query_disease":
                    result = await self.query_disease(
                        disease_id=input_data.get("disease_id"),
                        search_term=input_data.get("search_term")
                    )
                elif operation == "query_drug":
                    result = await self.query_drug(
                        drug_id=input_data.get("drug_id"),
                        search_term=input_data.get("search_term")
                    )
                else:
                    # Try to infer from dict keys
                    if "pathway_id" in input_data:
                        result = await self.query_pathway(pathway_id=input_data["pathway_id"])
                    elif "gene_id" in input_data:
                        result = await self.query_gene(gene_id=input_data["gene_id"])
                    elif "compound_id" in input_data:
                        result = await self.query_compound(compound_id=input_data["compound_id"])
                    elif "disease_id" in input_data:
                        result = await self.query_disease(disease_id=input_data["disease_id"])
                    elif "drug_id" in input_data:
                        result = await self.query_drug(drug_id=input_data["drug_id"])
                    else:
                        return AdapterResult(
                            success=False,
                            data=None,
                            error="Invalid dict input: must specify operation or provide specific ID field"
                        )

            elif isinstance(input_data, str):
                # Try to infer query type from input string
                input_str = input_data.strip()

                # Check if it looks like a KEGG ID
                if ':' in input_str:
                    # Gene ID (e.g., "hsa:10458")
                    result = await self.query_gene(gene_id=input_str)
                elif input_str.startswith('C') and len(input_str) == 6 and input_str[1:].isdigit():
                    # Compound ID (e.g., "C00031")
                    result = await self.query_compound(compound_id=input_str)
                elif input_str.startswith('D') and len(input_str) == 6 and input_str[1:].isdigit():
                    # Drug ID (e.g., "D00001")
                    result = await self.query_drug(drug_id=input_str)
                elif input_str.startswith('H') and len(input_str) == 6 and input_str[1:].isdigit():
                    # Disease ID (e.g., "H00001")
                    result = await self.query_disease(disease_id=input_str)
                elif input_str.startswith('hsa') or input_str.startswith('map'):
                    # Pathway ID (e.g., "hsa04010", "map00010")
                    result = await self.query_pathway(pathway_id=input_str)
                else:
                    # Generic search - try compound first, then pathway
                    operation = kwargs.get("operation", "query_compound")
                    if operation == "query_pathway":
                        result = await self.query_pathway(search_term=input_str)
                    elif operation == "query_gene":
                        result = await self.query_gene(search_term=input_str)
                    elif operation == "query_disease":
                        result = await self.query_disease(search_term=input_str)
                    elif operation == "query_drug":
                        result = await self.query_drug(search_term=input_str)
                    else:
                        result = await self.query_compound(search_term=input_str)

            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Invalid input type"
                )

            # Check if result contains error
            if isinstance(result, dict) and "error" in result:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=result["error"]
                )

            return AdapterResult(
                success=True,
                data=result,
                cache_hit=False,
                metadata={
                    "source": "kegg",
                    "adapter_version": self.version,
                    "query_type": type(input_data).__name__
                }
            )

        except Exception as e:
            logger.error(f"KEGG: Unexpected error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"KEGG query failed: {str(e)}"
            )
