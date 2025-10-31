"""
Enhanced PubChem Adapter - Complete PubChem API integration
Adds similarity search, substructure search, bioassays, vendors, and patents
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class PubChemEnhancedAdapter(AdapterProtocol):
    """
    Enhanced adapter for PubChem REST API

    Capabilities:
    - Compound property lookup (existing)
    - Similarity search (new)
    - Substructure search (new)
    - Bioassay data (new)
    - Vendor/purchasability information (new)
    - Patent mentions (new)
    - Compound synonyms and names (new)
    """

    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self):
        super().__init__(
            name="pubchem_enhanced",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.2,  # 5 requests/second
                "timeout": 60,
                "default_similarity": 90,  # 90% similarity
                "max_results": 100
            }
        )
        self.version = "2.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input (SMILES, CID, name, or InChI)

        Args:
            input_data: Input to validate

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, int):
            return input_data > 0
        elif isinstance(input_data, dict):
            return any(k in input_data for k in ["smiles", "cid", "name", "inchi"])
        return False

    async def _fetch_properties_async(
        self,
        smiles: str,
        properties: Optional[List[str]] = None
    ) -> Optional[Dict[str, Any]]:
        """
        Fetch molecular properties from PubChem API

        Args:
            smiles: SMILES string
            properties: List of properties to fetch (or None for default)

        Returns:
            Dictionary of molecular properties or None
        """
        encoded_smiles = quote(smiles)

        # Default properties if not specified
        if properties is None:
            properties = [
                "MolecularWeight",
                "XLogP",
                "TPSA",
                "HBondDonorCount",
                "HBondAcceptorCount",
                "RotatableBondCount",
                "Complexity",
                "CanonicalSMILES",
                "InChI",
                "InChIKey",
                "IUPACName"
            ]

        url = f"{self.BASE_URL}/compound/smiles/{encoded_smiles}/property/{','.join(properties)}/JSON"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
                        logger.info(f"PubChem: Successfully fetched properties for {smiles}")
                        return props
                    elif response.status == 404:
                        logger.warning(f"PubChem: Compound not found for SMILES: {smiles}")
                        return None
                    else:
                        error_text = await response.text()
                        logger.error(f"PubChem API error {response.status}: {error_text[:200]}")
                        return None
        except Exception as e:
            logger.error(f"PubChem: Error fetching properties for {smiles}: {e}")
            return None

    async def _similarity_search_async(
        self,
        smiles: str,
        similarity_threshold: int = 90,
        max_results: int = 100
    ) -> Optional[List[int]]:
        """
        Perform similarity search and return list of CIDs

        Args:
            smiles: SMILES string
            similarity_threshold: Similarity % (0-100)
            max_results: Maximum number of results

        Returns:
            List of PubChem CIDs or None
        """
        encoded_smiles = quote(smiles)
        url = f"{self.BASE_URL}/compound/fastsimilarity_2d/smiles/{encoded_smiles}/cids/JSON"
        params = {
            "Threshold": similarity_threshold,
            "MaxRecords": max_results
        }
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        cids = data.get("IdentifierList", {}).get("CID", [])
                        logger.info(f"PubChem: Found {len(cids)} similar compounds")
                        return cids
                    elif response.status == 404:
                        logger.info(f"PubChem: No similar compounds found")
                        return []
                    else:
                        logger.error(f"PubChem similarity search error {response.status}")
                        return None
        except Exception as e:
            logger.error(f"PubChem: Error in similarity search: {e}")
            return None

    async def _substructure_search_async(
        self,
        smiles: str,
        max_results: int = 100
    ) -> Optional[List[int]]:
        """
        Perform substructure search and return list of CIDs

        Args:
            smiles: SMILES string (substructure query)
            max_results: Maximum number of results

        Returns:
            List of PubChem CIDs or None
        """
        encoded_smiles = quote(smiles)
        url = f"{self.BASE_URL}/compound/fastsubstructure/smiles/{encoded_smiles}/cids/JSON"
        params = {"MaxRecords": max_results}
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        cids = data.get("IdentifierList", {}).get("CID", [])
                        logger.info(f"PubChem: Found {len(cids)} compounds with substructure")
                        return cids
                    elif response.status == 404:
                        logger.info(f"PubChem: No compounds found with substructure")
                        return []
                    else:
                        logger.error(f"PubChem substructure search error {response.status}")
                        return None
        except Exception as e:
            logger.error(f"PubChem: Error in substructure search: {e}")
            return None

    async def _get_bioassays_async(self, cid: int) -> Optional[List[Dict[str, Any]]]:
        """
        Get bioassay data for a compound

        Args:
            cid: PubChem Compound ID

        Returns:
            List of bioassay records or None
        """
        url = f"{self.BASE_URL}/compound/cid/{cid}/assaysummary/JSON"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        assays = data.get("Table", {}).get("Row", [])
                        logger.info(f"PubChem: Found {len(assays)} bioassays for CID {cid}")
                        return assays
                    elif response.status == 404:
                        logger.info(f"PubChem: No bioassays found for CID {cid}")
                        return []
                    else:
                        logger.warning(f"PubChem bioassay query error {response.status}")
                        return []
        except Exception as e:
            logger.error(f"PubChem: Error fetching bioassays: {e}")
            return None

    async def _get_vendors_async(self, cid: int) -> Optional[List[str]]:
        """
        Get vendor/supplier information for a compound

        Args:
            cid: PubChem Compound ID

        Returns:
            List of vendor names or None
        """
        url = f"{self.BASE_URL}/compound/cid/{cid}/synonyms/JSON"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        info = data.get("InformationList", {}).get("Information", [{}])[0]
                        synonyms = info.get("Synonym", [])

                        # Filter for vendor-related synonyms
                        vendors = [s for s in synonyms if any(
                            v in s.upper() for v in ["SIGMA", "ALDRICH", "TCI", "ACROS", "ALFA"]
                        )]

                        logger.info(f"PubChem: Found {len(vendors)} vendor entries for CID {cid}")
                        return vendors
                    elif response.status == 404:
                        return []
                    else:
                        return []
        except Exception as e:
            logger.error(f"PubChem: Error fetching vendors: {e}")
            return None

    async def _get_patents_async(self, cid: int) -> Optional[List[str]]:
        """
        Get patents mentioning a compound

        Args:
            cid: PubChem Compound ID

        Returns:
            List of patent IDs or None
        """
        url = f"{self.BASE_URL}/compound/cid/{cid}/xrefs/PatentID/JSON"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        patents = data.get("InformationList", {}).get("Information", [{}])[0].get("PatentID", [])
                        logger.info(f"PubChem: Found {len(patents)} patents for CID {cid}")
                        return patents
                    elif response.status == 404:
                        return []
                    else:
                        return []
        except Exception as e:
            logger.error(f"PubChem: Error fetching patents: {e}")
            return None

    async def _get_cid_from_smiles_async(self, smiles: str) -> Optional[int]:
        """
        Get PubChem CID from SMILES

        Args:
            smiles: SMILES string

        Returns:
            CID or None
        """
        props = await self._fetch_properties_async(smiles, ["CID"])
        if props:
            return props.get("CID")
        return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute enhanced PubChem query

        Args:
            input_data: SMILES, CID, or search dict
            **kwargs: Additional parameters:
                - mode: "properties", "similarity", "substructure", "bioassays", "vendors", "patents"
                - similarity_threshold: int (0-100)
                - max_results: int
                - include_bioassays: bool
                - include_vendors: bool
                - include_patents: bool

        Returns:
            AdapterResult containing query results
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
            smiles = input_data
            cid = None
        elif isinstance(input_data, int):
            cid = input_data
            smiles = None
        else:
            smiles = input_data.get("smiles")
            cid = input_data.get("cid")

        # Get mode
        mode = kwargs.get("mode", "properties")
        similarity_threshold = kwargs.get(
            "similarity_threshold",
            self.config.get("default_similarity", 90)
        )
        max_results = kwargs.get("max_results", self.config.get("max_results", 100))

        # Rate limiting
        await asyncio.sleep(self.config.get("rate_limit_delay", 0.2))

        # Execute based on mode
        try:
            if mode == "properties" and smiles:
                props = await self._fetch_properties_async(smiles)
                if props is None:
                    return AdapterResult(success=False, data=None, error="Failed to fetch properties")

                # Convert to standard format
                result_data = {
                    "cid": props.get("CID"),
                    "molecular_weight": float(props.get("MolecularWeight")) if props.get("MolecularWeight") else None,
                    "logp": float(props.get("XLogP")) if props.get("XLogP") else None,
                    "tpsa": float(props.get("TPSA")) if props.get("TPSA") else None,
                    "h_bond_donors": int(props.get("HBondDonorCount")) if props.get("HBondDonorCount") else None,
                    "h_bond_acceptors": int(props.get("HBondAcceptorCount")) if props.get("HBondAcceptorCount") else None,
                    "rotatable_bonds": int(props.get("RotatableBondCount")) if props.get("RotatableBondCount") else None,
                    "complexity": float(props.get("Complexity")) if props.get("Complexity") else None,
                    "canonical_smiles": props.get("CanonicalSMILES"),
                    "inchi": props.get("InChI"),
                    "inchikey": props.get("InChIKey"),
                    "iupac_name": props.get("IUPACName")
                }

            elif mode == "similarity" and smiles:
                cids = await self._similarity_search_async(smiles, similarity_threshold, max_results)
                if cids is None:
                    return AdapterResult(success=False, data=None, error="Similarity search failed")

                result_data = {
                    "query_smiles": smiles,
                    "similarity_threshold": similarity_threshold,
                    "num_results": len(cids),
                    "cids": cids
                }

            elif mode == "substructure" and smiles:
                cids = await self._substructure_search_async(smiles, max_results)
                if cids is None:
                    return AdapterResult(success=False, data=None, error="Substructure search failed")

                result_data = {
                    "query_smiles": smiles,
                    "num_results": len(cids),
                    "cids": cids
                }

            elif mode == "bioassays":
                if not cid and smiles:
                    cid = await self._get_cid_from_smiles_async(smiles)
                if not cid:
                    return AdapterResult(success=False, data=None, error="CID required for bioassays")

                assays = await self._get_bioassays_async(cid)
                if assays is None:
                    return AdapterResult(success=False, data=None, error="Failed to fetch bioassays")

                result_data = {
                    "cid": cid,
                    "num_assays": len(assays),
                    "bioassays": assays
                }

            elif mode == "vendors":
                if not cid and smiles:
                    cid = await self._get_cid_from_smiles_async(smiles)
                if not cid:
                    return AdapterResult(success=False, data=None, error="CID required for vendors")

                vendors = await self._get_vendors_async(cid)
                if vendors is None:
                    return AdapterResult(success=False, data=None, error="Failed to fetch vendors")

                result_data = {
                    "cid": cid,
                    "num_vendors": len(vendors),
                    "vendors": vendors,
                    "purchasable": len(vendors) > 0
                }

            elif mode == "patents":
                if not cid and smiles:
                    cid = await self._get_cid_from_smiles_async(smiles)
                if not cid:
                    return AdapterResult(success=False, data=None, error="CID required for patents")

                patents = await self._get_patents_async(cid)
                if patents is None:
                    return AdapterResult(success=False, data=None, error="Failed to fetch patents")

                result_data = {
                    "cid": cid,
                    "num_patents": len(patents),
                    "patents": patents[:20]  # Limit to first 20
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
                    "source": "pubchem",
                    "adapter_version": self.version,
                    "mode": mode
                }
            )

        except Exception as e:
            logger.error(f"PubChem: Unexpected error: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unexpected error: {str(e)}"
            )
