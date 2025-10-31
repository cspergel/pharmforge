"""
DrugCentral Adapter - Access drug information, targets, indications, and pharmacology
API Documentation: https://www.ebi.ac.uk/chembl/api/data/docs
Note: Migrated to ChEMBL API due to Pharos REST API deprecation
Updated: 2025-10 - Now using ChEMBL Web Services for drug data access
ChEMBL provides comprehensive drug information including targets, mechanisms, and indications
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DrugCentralAdapter(AdapterProtocol):
    """
    Adapter for DrugCentral database via ChEMBL API
    Provides drug information including targets, MOA, indications, and pharmacokinetics

    Note: Migrated from deprecated Pharos REST API to ChEMBL Web Services
    ChEMBL provides comprehensive drug discovery data with better API reliability
    """

    # ChEMBL API endpoint (provides comprehensive drug data)
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self):
        super().__init__(
            name="drugcentral",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,
                "timeout": 60,
                "max_retries": 3
            }
        )
        self.version = "3.0.0"  # Updated to reflect ChEMBL API migration

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid drug name or structure

        Args:
            input_data: Drug name string or dict with search parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return "drug" in input_data or "structure" in input_data or "inchikey" in input_data
        return False

    async def _retry_request(self, func, *args, **kwargs):
        """
        Retry a request with exponential backoff

        Args:
            func: Async function to retry
            *args, **kwargs: Arguments to pass to func

        Returns:
            Result of func or None on failure
        """
        max_retries = self.config.get("max_retries", 3)

        for attempt in range(max_retries):
            try:
                result = await func(*args, **kwargs)
                return result
            except (asyncio.TimeoutError, aiohttp.ClientError) as e:
                if attempt == max_retries - 1:
                    logger.error(f"DrugCentral: All {max_retries} retry attempts failed: {e}")
                    return None

                wait_time = 2 ** attempt  # Exponential backoff: 1s, 2s, 4s
                logger.warning(f"DrugCentral: Retry {attempt + 1}/{max_retries} after {wait_time}s delay")
                await asyncio.sleep(wait_time)

        return None

    async def _search_drug(self, query: str) -> Optional[List[Dict[str, Any]]]:
        """
        Search for drugs by name using ChEMBL API

        Args:
            query: Drug name to search

        Returns:
            List of matching drugs or None on error
        """
        # Use ChEMBL molecule search endpoint
        url = f"{self.BASE_URL}/molecule/search.json"
        params = {
            "q": query,
            "limit": 10  # Limit results
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        headers = {
            'User-Agent': 'PharmForge/3.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        logger.info(f"DrugCentral: No results found for {query}")
                        return []
                    elif response.status != 200:
                        logger.warning(f"DrugCentral search failed: {response.status}")
                        text = await response.text()
                        logger.debug(f"Response: {text[:500]}")
                        return None

                    data = await response.json()

                    # Parse ChEMBL response format
                    if isinstance(data, dict):
                        results = data.get("molecules", [])
                    elif isinstance(data, list):
                        results = data
                    else:
                        results = []

                    return results

        return await self._retry_request(make_request)

    async def _get_drug_targets(self, drug_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get target information for a drug from ChEMBL API

        Args:
            drug_id: ChEMBL molecule ID

        Returns:
            List of drug targets or None on error
        """
        url = f"{self.BASE_URL}/activity.json"
        params = {
            "molecule_chembl_id": drug_id,
            "limit": 100
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        headers = {
            'User-Agent': 'PharmForge/3.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        return []
                    elif response.status != 200:
                        logger.warning(f"DrugCentral targets request failed: {response.status}")
                        return None

                    data = await response.json()

                    # Parse response and extract unique targets
                    if isinstance(data, dict):
                        activities = data.get("activities", [])
                        # Extract unique targets from activities
                        targets_dict = {}
                        for activity in activities:
                            target_chembl_id = activity.get("target_chembl_id")
                            if target_chembl_id and target_chembl_id not in targets_dict:
                                targets_dict[target_chembl_id] = {
                                    "target_chembl_id": target_chembl_id,
                                    "target_pref_name": activity.get("target_pref_name"),
                                    "target_organism": activity.get("target_organism"),
                                    "target_type": activity.get("target_type")
                                }
                        return list(targets_dict.values())
                    return []

        return await self._retry_request(make_request)

    async def _get_drug_indications(self, drug_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get approved indications for a drug from ChEMBL API

        Args:
            drug_id: ChEMBL molecule ID

        Returns:
            List of indications or None on error
        """
        url = f"{self.BASE_URL}/drug_indication.json"
        params = {
            "molecule_chembl_id": drug_id,
            "limit": 50
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        headers = {
            'User-Agent': 'PharmForge/3.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        return []
                    elif response.status != 200:
                        logger.warning(f"DrugCentral indications request failed: {response.status}")
                        return None

                    data = await response.json()

                    # Extract indications from ChEMBL response
                    if isinstance(data, dict):
                        indications = data.get("drug_indications", [])
                        return indications if isinstance(indications, list) else []
                    return []

        return await self._retry_request(make_request)

    async def _get_drug_pharmacology(self, drug_id: str) -> Optional[Dict[str, Any]]:
        """
        Get pharmacology information for a drug from ChEMBL API

        Args:
            drug_id: ChEMBL molecule ID

        Returns:
            Pharmacology data or None on error
        """
        url = f"{self.BASE_URL}/mechanism.json"
        params = {
            "molecule_chembl_id": drug_id,
            "limit": 50
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        headers = {
            'User-Agent': 'PharmForge/3.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        return {}
                    elif response.status != 200:
                        logger.warning(f"DrugCentral pharmacology request failed: {response.status}")
                        return None

                    data = await response.json()

                    # Extract mechanisms from ChEMBL response
                    if isinstance(data, dict):
                        mechanisms = data.get("mechanisms", [])
                        pharmacology = {
                            "mechanisms": mechanisms,
                            "mechanism_of_action": mechanisms[0].get("mechanism_of_action") if mechanisms else None,
                            "action_type": mechanisms[0].get("action_type") if mechanisms else None
                        }
                        return pharmacology
                    return {}

        return await self._retry_request(make_request)

    async def _get_drug_properties(self, drug_id: str) -> Optional[Dict[str, Any]]:
        """
        Get chemical properties for a drug from ChEMBL API
        Properties are already included in molecule data, so this returns None
        and properties are extracted from the main molecule object

        Args:
            drug_id: ChEMBL molecule ID

        Returns:
            Properties data (None as they come with molecule data)
        """
        # Properties are included in the molecule search results
        # No separate API call needed
        return None

    def _extract_drug_info(
        self,
        drug: Dict[str, Any],
        targets: Optional[List[Dict[str, Any]]] = None,
        indications: Optional[List[Dict[str, Any]]] = None,
        pharmacology: Optional[Dict[str, Any]] = None,
        properties: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Compile comprehensive drug information from ChEMBL API response

        Args:
            drug: Basic drug data from search
            targets: Target information
            indications: Indication information
            pharmacology: Pharmacology data
            properties: Chemical properties (unused, included in drug object)

        Returns:
            Compiled drug information
        """
        # ChEMBL API field mapping
        mol_structures = drug.get("molecule_structures", {})
        mol_properties = drug.get("molecule_properties", {})

        info = {
            "drug_id": drug.get("molecule_chembl_id"),
            "name": drug.get("pref_name"),
            "generic_name": drug.get("pref_name"),
            "synonyms": [syn.get("molecule_synonym") for syn in drug.get("molecule_synonyms", [])[:10]],
            "smiles": mol_structures.get("canonical_smiles"),
            "inchikey": mol_structures.get("standard_inchi_key"),
            "max_phase": drug.get("max_phase"),
            "first_approval": drug.get("first_approval"),
            "therapeutic_flag": drug.get("therapeutic_flag"),
            "source": "chembl/drugcentral"
        }

        # Add targets (ChEMBL format)
        if targets:
            info["targets"] = [
                {
                    "target_chembl_id": t.get("target_chembl_id"),
                    "target_name": t.get("target_pref_name"),
                    "organism": t.get("target_organism"),
                    "target_type": t.get("target_type")
                }
                for t in targets
            ]
            info["num_targets"] = len(targets)
        else:
            info["targets"] = []
            info["num_targets"] = 0

        # Add indications (ChEMBL format)
        if indications:
            info["indications"] = [
                {
                    "disease": ind.get("mesh_heading"),
                    "mesh_id": ind.get("mesh_id"),
                    "efo_id": ind.get("efo_id"),
                    "max_phase": ind.get("max_phase_for_ind")
                }
                for ind in indications
            ]
            info["num_indications"] = len(indications)
        else:
            info["indications"] = []
            info["num_indications"] = 0

        # Add pharmacology (ChEMBL format)
        if pharmacology:
            info["pharmacology"] = {
                "mechanism_of_action": pharmacology.get("mechanism_of_action"),
                "action_type": pharmacology.get("action_type"),
                "mechanisms": pharmacology.get("mechanisms", [])
            }

        # Add chemical properties (ChEMBL format - from molecule_properties)
        info["properties"] = {
            "molecular_weight": mol_properties.get("full_mwt"),
            "alogp": mol_properties.get("alogp"),
            "psa": mol_properties.get("psa"),
            "hbd": mol_properties.get("hbd"),
            "hba": mol_properties.get("hba"),
            "rotatable_bonds": mol_properties.get("rtb"),
            "aromatic_rings": mol_properties.get("aromatic_rings"),
            "qed_weighted": mol_properties.get("qed_weighted"),
            "ro5_violations": mol_properties.get("num_ro5_violations")
        }

        return info

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute DrugCentral drug information lookup via ChEMBL API

        Args:
            input_data: Drug name string or dict with search parameters
                       Examples:
                       - "aspirin"
                       - {"drug": "aspirin"}
                       - {"inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}
            **kwargs: Additional parameters:
                     - include_targets: bool - Fetch target information (default: True)
                     - include_indications: bool - Fetch indications (default: True)
                     - include_pharmacology: bool - Fetch pharmacology (default: True)
                     - include_properties: bool - Always included in molecule data (default: True)

        Returns:
            AdapterResult containing comprehensive drug information from ChEMBL
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a drug name or search parameters dict"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Extract search query
        if isinstance(input_data, str):
            query = input_data
        else:
            query = input_data.get("drug") or input_data.get("structure") or input_data.get("inchikey")

        # Search for drug
        search_results = await self._search_drug(query)

        if search_results is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search DrugCentral",
                metadata={"query": query}
            )

        if not search_results or len(search_results) == 0:
            return AdapterResult(
                success=True,
                data={
                    "drugs": [],
                    "total_found": 0
                },
                metadata={"query": query, "source": "drugcentral"}
            )

        # Get detailed information for each drug found
        drugs_detailed = []

        for drug in search_results[:5]:  # Limit to top 5 results
            drug_id = drug.get("molecule_chembl_id")

            if not drug_id:
                continue

            # Fetch additional data if requested
            targets = None
            indications = None
            pharmacology = None
            properties = None

            if kwargs.get("include_targets", True):
                targets = await self._get_drug_targets(drug_id)
                await asyncio.sleep(0.2)  # Rate limit

            if kwargs.get("include_indications", True):
                indications = await self._get_drug_indications(drug_id)
                await asyncio.sleep(0.2)

            if kwargs.get("include_pharmacology", True):
                pharmacology = await self._get_drug_pharmacology(drug_id)
                await asyncio.sleep(0.2)

            if kwargs.get("include_properties", True):
                properties = await self._get_drug_properties(drug_id)
                await asyncio.sleep(0.2)

            # Compile drug info
            drug_info = self._extract_drug_info(
                drug,
                targets,
                indications,
                pharmacology,
                properties
            )
            drugs_detailed.append(drug_info)

        result_data = {
            "drugs": drugs_detailed,
            "total_found": len(search_results)
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "query": query,
                "source": "drugcentral",
                "adapter_version": self.version
            }
        )
