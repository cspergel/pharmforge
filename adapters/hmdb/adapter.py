"""
HMDB Adapter for PharmForge

Access Human Metabolome Database for metabolite data, biofluid concentrations,
disease associations, and metabolic pathway information.

HMDB contains 220,000+ human metabolites with rich clinical data including:
- Normal concentration ranges in biofluids (blood, urine, CSF, saliva)
- Disease associations and biomarker relationships
- Metabolic pathway information
- Protein/enzyme associations
- NMR/MS spectral data references

Reference: https://hmdb.ca/
"""

from typing import Any, Dict, Optional, List, Union
import aiohttp
import asyncio
import logging
import xml.etree.ElementTree as ET
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class HMDBAdapter(AdapterProtocol):
    """
    HMDB (Human Metabolome Database) adapter for metabolomics data

    Provides access to:
    - Metabolite structures and properties
    - Biofluid concentration data (normal ranges)
    - Disease associations and biomarker information
    - Metabolic pathway data
    - Protein/enzyme interactions
    - Cross-references (KEGG, ChEBI, PubChem)

    The adapter queries HMDB's web API and parses XML responses.
    """

    BASE_URL = "https://hmdb.ca"
    API_URL = "https://hmdb.ca/metabolites"

    def __init__(self):
        super().__init__(
            name="hmdb",
            adapter_type="api",
            config={
                "rate_limit_delay": 1.0,  # Be respectful to HMDB servers
                "timeout": 60,
                "max_results": 20,
                "default_biofluids": ["blood", "urine", "cerebrospinal_fluid"]
            }
        )
        self.version = "1.0.0"

        # HMDB namespace for XML parsing
        self.ns = {"hmdb": "http://www.hmdb.ca"}

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data

        Args:
            input_data: HMDB ID, metabolite name, SMILES, or dict

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return any(k in input_data for k in [
                "hmdb_id", "name", "smiles", "formula",
                "pathway", "disease", "query"
            ])
        return False

    async def _fetch_metabolite_xml_async(self, hmdb_id: str) -> Optional[str]:
        """
        Fetch metabolite XML data from HMDB

        Args:
            hmdb_id: HMDB ID (e.g., HMDB0000001)

        Returns:
            XML string or None
        """
        # Ensure proper HMDB ID format
        if not hmdb_id.startswith("HMDB"):
            hmdb_id = f"HMDB{hmdb_id.zfill(7)}"

        url = f"{self.API_URL}/{hmdb_id}.xml"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        xml_data = await response.text()
                        logger.info(f"HMDB: Successfully fetched data for {hmdb_id}")
                        return xml_data
                    elif response.status == 404:
                        logger.warning(f"HMDB: Metabolite not found: {hmdb_id}")
                        return None
                    else:
                        logger.error(f"HMDB API error {response.status}")
                        return None
        except Exception as e:
            logger.error(f"HMDB: Error fetching metabolite {hmdb_id}: {e}")
            return None

    async def _search_metabolite_async(
        self,
        query: str,
        search_type: str = "name"
    ) -> Optional[List[str]]:
        """
        Search for metabolites by name, formula, or other criteria

        Args:
            query: Search query
            search_type: Type of search ("name", "formula", "mass")

        Returns:
            List of HMDB IDs or None
        """
        # Note: HMDB's search API is limited, this is a simplified implementation
        # In production, you might want to use the HMDB text search or download
        # the full database and implement local search

        url = f"{self.BASE_URL}/unearth/q"
        params = {
            "utf8": "✓",
            "query": query,
            "searcher": "metabolites"
        }
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        # Parse HTML response to extract HMDB IDs
                        # This is a simplified version - production code would need
                        # better HTML parsing
                        text = await response.text()
                        # Extract HMDB IDs from response
                        import re
                        hmdb_ids = re.findall(r'HMDB\d{7}', text)
                        unique_ids = list(set(hmdb_ids))[:self.config.get("max_results", 20)]
                        logger.info(f"HMDB: Found {len(unique_ids)} results for query: {query}")
                        return unique_ids
                    return None
        except Exception as e:
            logger.error(f"HMDB: Error searching for '{query}': {e}")
            return None

    def _parse_metabolite_xml(self, xml_data: str) -> Optional[Dict[str, Any]]:
        """
        Parse metabolite XML data into structured format

        Args:
            xml_data: XML string from HMDB

        Returns:
            Parsed metabolite data or None
        """
        try:
            root = ET.fromstring(xml_data)

            # Helper function to safely get text from XML element
            def get_text(element, path, default=None):
                elem = element.find(path, self.ns)
                return elem.text if elem is not None else default

            # Helper function to get list of texts
            def get_texts(element, path):
                return [e.text for e in element.findall(path, self.ns) if e.text]

            # Extract basic information
            metabolite = {
                "hmdb_id": get_text(root, "hmdb:accession"),
                "name": get_text(root, "hmdb:name"),
                "systematic_name": get_text(root, "hmdb:systematic_name"),
                "chemical_formula": get_text(root, "hmdb:chemical_formula"),
                "average_molecular_weight": get_text(root, "hmdb:average_molecular_weight"),
                "monoisotopic_molecular_weight": get_text(root, "hmdb:monoisotopic_molecular_weight"),
                "iupac_name": get_text(root, "hmdb:iupac_name"),
                "smiles": get_text(root, "hmdb:smiles"),
                "inchi": get_text(root, "hmdb:inchi"),
                "inchikey": get_text(root, "hmdb:inchikey"),
                "state": get_text(root, "hmdb:state"),
                "description": get_text(root, "hmdb:description"),
            }

            # Parse synonyms
            synonyms = get_texts(root, ".//hmdb:synonym")
            metabolite["synonyms"] = synonyms[:10] if synonyms else []

            # Parse biofluid locations
            biofluid_locations = get_texts(root, ".//hmdb:biofluid_locations/hmdb:biofluid")
            metabolite["biofluid_locations"] = biofluid_locations

            # Parse tissue locations
            tissue_locations = get_texts(root, ".//hmdb:tissue_locations/hmdb:tissue")
            metabolite["tissue_locations"] = tissue_locations

            # Parse normal concentrations
            concentrations = {}
            normal_conc = root.find("hmdb:normal_concentrations", self.ns)
            if normal_conc is not None:
                for conc in normal_conc.findall("hmdb:concentration", self.ns):
                    biofluid = get_text(conc, "hmdb:biofluid")
                    if biofluid:
                        concentrations[biofluid.lower().replace(" ", "_")] = {
                            "value": get_text(conc, "hmdb:concentration_value"),
                            "unit": get_text(conc, "hmdb:concentration_units"),
                            "subject_age": get_text(conc, "hmdb:subject_age"),
                            "subject_sex": get_text(conc, "hmdb:subject_sex"),
                            "subject_condition": get_text(conc, "hmdb:subject_condition"),
                            "references": get_texts(conc, ".//hmdb:reference/hmdb:pubmed_id")
                        }
            metabolite["concentrations"] = concentrations

            # Parse diseases
            diseases = []
            diseases_elem = root.find("hmdb:diseases", self.ns)
            if diseases_elem is not None:
                for disease in diseases_elem.findall("hmdb:disease", self.ns):
                    diseases.append({
                        "name": get_text(disease, "hmdb:name"),
                        "omim_id": get_text(disease, "hmdb:omim_id"),
                        "references": get_texts(disease, ".//hmdb:reference/hmdb:pubmed_id")[:5]
                    })
            metabolite["diseases"] = diseases[:20]  # Limit to 20 diseases

            # Parse pathways
            pathways = []
            pathways_elem = root.find("hmdb:pathways", self.ns)
            if pathways_elem is not None:
                for pathway in pathways_elem.findall("hmdb:pathway", self.ns):
                    pathways.append({
                        "name": get_text(pathway, "hmdb:name"),
                        "smpdb_id": get_text(pathway, "hmdb:smpdb_id"),
                        "kegg_map_id": get_text(pathway, "hmdb:kegg_map_id")
                    })
            metabolite["pathways"] = pathways

            # Parse protein associations
            proteins = []
            proteins_elem = root.find("hmdb:protein_associations", self.ns)
            if proteins_elem is not None:
                for protein in proteins_elem.findall("hmdb:protein", self.ns):
                    uniprot_id = get_text(protein, "hmdb:uniprot_id")
                    gene_name = get_text(protein, "hmdb:gene_name")
                    protein_type = get_text(protein, "hmdb:protein_type")
                    if uniprot_id or gene_name:
                        proteins.append({
                            "uniprot_id": uniprot_id,
                            "gene_name": gene_name,
                            "protein_type": protein_type,
                            "name": get_text(protein, "hmdb:name")
                        })
            metabolite["proteins"] = proteins[:20]  # Limit to 20 proteins

            # Parse ontology terms
            ontology = {
                "description": get_text(root, "hmdb:ontology/hmdb:description"),
                "origins": get_texts(root, ".//hmdb:ontology/hmdb:origins/hmdb:origin"),
                "biofunctions": get_texts(root, ".//hmdb:ontology/hmdb:biofunctions/hmdb:biofunction"),
                "applications": get_texts(root, ".//hmdb:ontology/hmdb:applications/hmdb:application")
            }
            metabolite["ontology"] = ontology

            # Parse external IDs
            external_ids = {
                "pubchem_compound_id": get_text(root, "hmdb:pubchem_compound_id"),
                "chemspider_id": get_text(root, "hmdb:chemspider_id"),
                "kegg_id": get_text(root, "hmdb:kegg_id"),
                "chebi_id": get_text(root, "hmdb:chebi_id"),
                "drugbank_id": get_text(root, "hmdb:drugbank_id"),
                "foodb_id": get_text(root, "hmdb:foodb_id")
            }
            metabolite["external_ids"] = external_ids

            logger.info(f"HMDB: Successfully parsed data for {metabolite['hmdb_id']}")
            return metabolite

        except Exception as e:
            logger.error(f"HMDB: Error parsing XML: {e}")
            return None

    def _filter_by_biofluids(
        self,
        metabolite: Dict[str, Any],
        biofluids: List[str]
    ) -> Dict[str, Any]:
        """
        Filter concentration data by specific biofluids

        Args:
            metabolite: Parsed metabolite data
            biofluids: List of biofluids to include

        Returns:
            Filtered metabolite data
        """
        if not biofluids or not metabolite.get("concentrations"):
            return metabolite

        # Normalize biofluid names
        normalized = [b.lower().replace(" ", "_") for b in biofluids]

        # Filter concentrations
        filtered_conc = {
            k: v for k, v in metabolite["concentrations"].items()
            if k in normalized
        }

        metabolite["concentrations"] = filtered_conc
        return metabolite

    def _summarize_metabolite(self, metabolite: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create summary statistics for a metabolite

        Args:
            metabolite: Parsed metabolite data

        Returns:
            Summary statistics
        """
        return {
            "hmdb_id": metabolite.get("hmdb_id"),
            "name": metabolite.get("name"),
            "formula": metabolite.get("chemical_formula"),
            "molecular_weight": metabolite.get("average_molecular_weight"),
            "num_biofluids": len(metabolite.get("biofluid_locations", [])),
            "num_tissues": len(metabolite.get("tissue_locations", [])),
            "num_concentrations": len(metabolite.get("concentrations", {})),
            "num_diseases": len(metabolite.get("diseases", [])),
            "num_pathways": len(metabolite.get("pathways", [])),
            "num_proteins": len(metabolite.get("proteins", [])),
            "has_pubchem": bool(metabolite.get("external_ids", {}).get("pubchem_compound_id")),
            "has_kegg": bool(metabolite.get("external_ids", {}).get("kegg_id")),
            "state": metabolite.get("state")
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute HMDB query

        Args:
            input_data: HMDB ID, name, or search dict
            **kwargs: Additional parameters:
                - mode: "metabolite", "search", "pathway", "disease"
                - biofluids: List of biofluids to include
                - include_concentrations: bool
                - include_diseases: bool
                - include_pathways: bool
                - include_proteins: bool
                - search_type: "name", "formula", "mass"

        Returns:
            AdapterResult containing metabolite data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input"
            )

        # Parse input
        if isinstance(input_data, str):
            # Check if it's an HMDB ID
            if input_data.startswith("HMDB") or input_data.isdigit():
                hmdb_id = input_data
                query = None
                mode = "metabolite"
            else:
                # It's a search query
                query = input_data
                hmdb_id = None
                mode = "search"
        else:
            hmdb_id = input_data.get("hmdb_id")
            query = input_data.get("query") or input_data.get("name")
            mode = input_data.get("mode", "metabolite" if hmdb_id else "search")

        # Get parameters
        mode = kwargs.get("mode", mode)
        biofluids = kwargs.get("biofluids", self.config.get("default_biofluids"))
        include_concentrations = kwargs.get("include_concentrations", True)
        include_diseases = kwargs.get("include_diseases", True)
        include_pathways = kwargs.get("include_pathways", True)
        include_proteins = kwargs.get("include_proteins", True)
        search_type = kwargs.get("search_type", "name")

        # Rate limiting
        await asyncio.sleep(self.config.get("rate_limit_delay", 1.0))

        try:
            if mode == "metabolite" and hmdb_id:
                # Fetch single metabolite
                xml_data = await self._fetch_metabolite_xml_async(hmdb_id)
                if not xml_data:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"Metabolite {hmdb_id} not found"
                    )

                metabolite = self._parse_metabolite_xml(xml_data)
                if not metabolite:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="Failed to parse metabolite data"
                    )

                # Filter by biofluids if specified
                if biofluids:
                    metabolite = self._filter_by_biofluids(metabolite, biofluids)

                # Remove fields if not requested
                if not include_concentrations:
                    metabolite.pop("concentrations", None)
                if not include_diseases:
                    metabolite.pop("diseases", None)
                if not include_pathways:
                    metabolite.pop("pathways", None)
                if not include_proteins:
                    metabolite.pop("proteins", None)

                summary = self._summarize_metabolite(metabolite)

                result_data = {
                    "metabolite": metabolite,
                    "summary": summary
                }

            elif mode == "search" and query:
                # Search for metabolites
                hmdb_ids = await self._search_metabolite_async(query, search_type)
                if hmdb_ids is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="Search failed"
                    )

                if not hmdb_ids:
                    return AdapterResult(
                        success=True,
                        data={
                            "query": query,
                            "num_results": 0,
                            "metabolites": [],
                            "warning": "No metabolites found"
                        }
                    )

                # Fetch details for each metabolite (limit to avoid overload)
                metabolites = []
                max_fetch = min(len(hmdb_ids), 5)  # Only fetch first 5 in detail

                for hmdb_id in hmdb_ids[:max_fetch]:
                    xml_data = await self._fetch_metabolite_xml_async(hmdb_id)
                    if xml_data:
                        metabolite = self._parse_metabolite_xml(xml_data)
                        if metabolite:
                            # Apply filters
                            if biofluids:
                                metabolite = self._filter_by_biofluids(metabolite, biofluids)
                            if not include_concentrations:
                                metabolite.pop("concentrations", None)
                            if not include_diseases:
                                metabolite.pop("diseases", None)
                            if not include_pathways:
                                metabolite.pop("pathways", None)
                            if not include_proteins:
                                metabolite.pop("proteins", None)

                            metabolites.append(metabolite)

                    # Rate limiting between requests
                    await asyncio.sleep(self.config.get("rate_limit_delay", 1.0))

                result_data = {
                    "query": query,
                    "search_type": search_type,
                    "num_results": len(hmdb_ids),
                    "num_fetched": len(metabolites),
                    "all_hmdb_ids": hmdb_ids,
                    "metabolites": metabolites,
                    "warning": f"Only fetched details for first {max_fetch} results" if len(hmdb_ids) > max_fetch else None
                }

            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Invalid mode '{mode}' or missing required input"
                )

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "hmdb",
                    "adapter_version": self.version,
                    "mode": mode,
                    "database_url": self.BASE_URL
                }
            )

        except Exception as e:
            logger.error(f"HMDB: Unexpected error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unexpected error: {str(e)}"
            )
