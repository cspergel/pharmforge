"""
CompTox Chemistry Dashboard Adapter - Access EPA chemistry and toxicity data
API Documentation: https://api-ccte.epa.gov/docs/
Database: https://comptox.epa.gov/dashboard/

The CompTox Dashboard provides access to 900,000+ chemicals with:
- Predicted toxicity endpoints (OPERA models)
- Exposure pathways and use categories
- Physicochemical properties
- Bioactivity assays (ToxCast/Tox21)
- Hazard classifications (GHS)
- Environmental fate predictions
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class CompToxAdapter(AdapterProtocol):
    """
    Adapter for EPA CompTox Chemistry Dashboard API
    Provides comprehensive toxicity and environmental chemistry data

    Features:
    - Chemical search by DTXSID, CAS, name, SMILES, InChIKey
    - QSAR toxicity predictions (OPERA models)
    - Exposure pathway analysis
    - Bioactivity data (ToxCast/Tox21)
    - GHS hazard classifications
    - Environmental fate predictions
    """

    # EPA CompTox API base URL
    BASE_URL = "https://api-ccte.epa.gov/chemical"

    def __init__(self):
        super().__init__(
            name="comptox",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.1,  # 10 requests/second max
                "timeout": 45,
                "max_retries": 3
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input for CompTox query

        Args:
            input_data: Query string or dict with query_type and query

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return "query" in input_data
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
                    logger.error(f"CompTox: All {max_retries} retry attempts failed: {e}")
                    return None

                wait_time = 2 ** attempt  # Exponential backoff
                logger.warning(f"CompTox: Retry {attempt + 1}/{max_retries} after {wait_time}s delay")
                await asyncio.sleep(wait_time)

        return None

    async def _search_chemical(
        self,
        query: str,
        query_type: str = "name"
    ) -> Optional[Dict[str, Any]]:
        """
        Search for a chemical in CompTox Dashboard

        Args:
            query: Search query (identifier or name)
            query_type: Type of query - "dtxsid", "cas", "name", "smiles", "inchikey"

        Returns:
            Chemical data or None on error
        """
        # Map query types to API endpoints
        endpoint_map = {
            "dtxsid": f"detail/search/by-dtxsid/{quote(query)}",
            "cas": f"detail/search/by-cas/{quote(query)}",
            "name": f"detail/search/by-name/{quote(query)}",
            "smiles": f"detail/search/by-smiles/{quote(query)}",
            "inchikey": f"detail/search/by-inchikey/{quote(query)}"
        }

        endpoint = endpoint_map.get(query_type, f"detail/search/by-name/{quote(query)}")
        url = f"{self.BASE_URL}/{endpoint}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 45))

        headers = {
            'User-Agent': 'PharmForge/1.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        logger.info(f"CompTox: Chemical not found for {query_type}: {query}")
                        return None
                    elif response.status != 200:
                        text = await response.text()
                        logger.warning(f"CompTox API error {response.status}: {text[:200]}")
                        return None

                    data = await response.json()
                    logger.info(f"CompTox: Found chemical for {query_type}: {query}")
                    return data

        return await self._retry_request(make_request)

    async def _get_toxicity_data(self, dtxsid: str) -> Optional[Dict[str, Any]]:
        """
        Get predicted toxicity endpoints for a chemical

        Args:
            dtxsid: DTXSID identifier

        Returns:
            Toxicity predictions or None on error
        """
        # Note: CompTox API endpoints for toxicity may vary
        # Using property search which includes OPERA predictions
        url = f"{self.BASE_URL}/property/search/by-dtxsid/{quote(dtxsid)}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 45))

        headers = {
            'User-Agent': 'PharmForge/1.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        logger.debug(f"CompTox: No toxicity data for {dtxsid}")
                        return {}
                    elif response.status != 200:
                        logger.warning(f"CompTox toxicity request failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        return await self._retry_request(make_request)

    async def _get_bioactivity_data(self, dtxsid: str) -> Optional[Dict[str, Any]]:
        """
        Get bioactivity assay data for a chemical

        Args:
            dtxsid: DTXSID identifier

        Returns:
            Bioactivity data or None on error
        """
        url = f"{self.BASE_URL}/bioactivity/search/by-dtxsid/{quote(dtxsid)}"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 45))

        headers = {
            'User-Agent': 'PharmForge/1.0 (Drug Discovery Platform; +https://pharmforge.org)',
            'Accept': 'application/json'
        }

        async def make_request():
            async with aiohttp.ClientSession() as session:
                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if response.status == 404:
                        logger.debug(f"CompTox: No bioactivity data for {dtxsid}")
                        return {}
                    elif response.status != 200:
                        logger.warning(f"CompTox bioactivity request failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        return await self._retry_request(make_request)

    def _extract_properties(self, chemical_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract physicochemical properties from chemical data

        Args:
            chemical_data: Raw API response

        Returns:
            Formatted properties dict
        """
        properties = {}

        # Common property fields in CompTox API
        prop_map = {
            "MOLECULAR_WEIGHT": "molecular_weight",
            "LOGP": "logp",
            "WATER_SOLUBILITY": "water_solubility",
            "VAPOR_PRESSURE": "vapor_pressure",
            "MELTING_POINT": "melting_point",
            "BOILING_POINT": "boiling_point",
            "HENRYS_LAW_CONSTANT": "henrys_law_constant"
        }

        for api_key, our_key in prop_map.items():
            if api_key in chemical_data:
                properties[our_key] = chemical_data[api_key]

        return properties

    def _extract_toxicity(self, toxicity_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract and format toxicity predictions

        Args:
            toxicity_data: Raw toxicity API response

        Returns:
            Formatted toxicity dict
        """
        if not toxicity_data:
            return {
                "predicted_endpoints": [],
                "qsar_predictions": {}
            }

        # Extract OPERA model predictions
        predicted_endpoints = []
        qsar_predictions = {}

        # Common OPERA endpoints
        opera_endpoints = {
            "OPERA_LogP": "logp",
            "OPERA_MP": "melting_point",
            "OPERA_BP": "boiling_point",
            "OPERA_LogOH": "oh_rate_constant",
            "OPERA_BCF": "bioconcentration_factor",
            "OPERA_BioDeg": "biodegradation",
            "OPERA_HL": "half_life",
            "OPERA_KOA": "air_water_partition",
            "OPERA_Ames": "mutagenicity",
            "OPERA_RT": "developmental_toxicity"
        }

        for api_key, endpoint_name in opera_endpoints.items():
            if api_key in toxicity_data:
                value = toxicity_data[api_key]
                predicted_endpoints.append({
                    "endpoint": endpoint_name,
                    "value": value,
                    "source": "OPERA model"
                })

        # Binary QSAR predictions
        if "OPERA_Ames" in toxicity_data:
            qsar_predictions["mutagenicity"] = "positive" if toxicity_data["OPERA_Ames"] else "negative"
        if "OPERA_RT" in toxicity_data:
            qsar_predictions["developmental_toxicity"] = "positive" if toxicity_data["OPERA_RT"] else "negative"

        return {
            "predicted_endpoints": predicted_endpoints,
            "qsar_predictions": qsar_predictions
        }

    def _extract_bioactivity(self, bioactivity_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract bioactivity information

        Args:
            bioactivity_data: Raw bioactivity API response

        Returns:
            Formatted bioactivity dict
        """
        if not bioactivity_data:
            return {
                "num_assays": 0,
                "active_assays": 0,
                "gene_targets": []
            }

        # Extract assay counts and targets
        assays = bioactivity_data.get("assays", [])

        num_assays = len(assays)
        active_assays = sum(1 for a in assays if a.get("activity", "").lower() == "active")

        # Extract unique gene targets
        gene_targets = list(set(
            a.get("gene_target")
            for a in assays
            if a.get("gene_target")
        ))

        return {
            "num_assays": num_assays,
            "active_assays": active_assays,
            "gene_targets": gene_targets[:20]  # Limit to top 20
        }

    def _extract_hazard(self, chemical_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract hazard classification information

        Args:
            chemical_data: Raw API response

        Returns:
            Formatted hazard dict
        """
        hazard = {
            "ghs_classification": [],
            "ecological_hazard": "unknown",
            "environmental_fate": {}
        }

        # GHS classifications
        if "ghs_classifications" in chemical_data:
            hazard["ghs_classification"] = chemical_data["ghs_classifications"]

        # Environmental fate
        if "biodegradation" in chemical_data:
            biodeg = chemical_data["biodegradation"]
            if biodeg and biodeg > 0.5:
                hazard["environmental_fate"]["biodegradation"] = "readily biodegradable"
                hazard["environmental_fate"]["persistence"] = "low"
            else:
                hazard["environmental_fate"]["biodegradation"] = "not readily biodegradable"
                hazard["environmental_fate"]["persistence"] = "high"

        return hazard

    def _extract_exposure(self, chemical_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract exposure pathway information

        Args:
            chemical_data: Raw API response

        Returns:
            Formatted exposure dict
        """
        exposure = {
            "consumer_products": [],
            "exposure_pathways": [],
            "use_categories": []
        }

        # Extract use and exposure info
        if "consumer_uses" in chemical_data:
            exposure["consumer_products"] = chemical_data["consumer_uses"]

        if "exposure_pathways" in chemical_data:
            exposure["exposure_pathways"] = chemical_data["exposure_pathways"]

        if "use_categories" in chemical_data:
            exposure["use_categories"] = chemical_data["use_categories"]

        return exposure

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute CompTox chemistry and toxicity lookup

        Args:
            input_data: Query string or dict with parameters:
                       - Simple: "aspirin" or "DTXSID7020182"
                       - Dict: {
                           "query": "aspirin",
                           "query_type": "name",  # or "dtxsid", "cas", "smiles", "inchikey"
                           "include_toxicity": true,
                           "include_bioactivity": true,
                           "include_exposure": true,
                           "include_hazards": true
                         }
            **kwargs: Additional parameters:
                     - include_toxicity: bool (default: True)
                     - include_bioactivity: bool (default: True)
                     - include_properties: bool (default: True)
                     - include_exposure: bool (default: True)
                     - include_hazards: bool (default: True)

        Returns:
            AdapterResult containing comprehensive EPA CompTox data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a query string or dict with 'query' field"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.1)
        await asyncio.sleep(rate_delay)

        # Parse input and auto-detect query type
        if isinstance(input_data, str):
            query = input_data
            # Auto-detect query type
            if query.startswith("DTXSID"):
                query_type = "dtxsid"
            elif "InChI=" in query:
                query_type = "inchikey"
            elif any(c in query for c in ["=", "#", "@", "(", ")"]):  # Likely SMILES
                query_type = "smiles"
            elif query.replace("-", "").isdigit() or (query.count("-") == 2 and all(p.isdigit() for p in query.split("-"))):
                query_type = "cas"  # CAS number format
            else:
                query_type = "name"  # Default to name search
        else:
            query = input_data["query"]
            query_type = input_data.get("query_type", "name")

        # Search for chemical
        chemical_data = await self._search_chemical(query, query_type)

        if chemical_data is None:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to find chemical in CompTox: {query}",
                metadata={"query": query, "query_type": query_type}
            )

        # Extract DTXSID for additional queries
        dtxsid = chemical_data.get("dtxsid")

        if not dtxsid:
            return AdapterResult(
                success=False,
                data=None,
                error="No DTXSID found in CompTox response",
                metadata={"query": query}
            )

        # Build result with basic chemical info
        result = {
            "chemical": {
                "dtxsid": dtxsid,
                "preferred_name": chemical_data.get("preferredName") or chemical_data.get("casrn"),
                "casrn": chemical_data.get("casrn"),
                "smiles": chemical_data.get("smiles"),
                "inchi": chemical_data.get("inchi"),
                "inchikey": chemical_data.get("inchikey"),
                "molecular_formula": chemical_data.get("molecularFormula"),
                "molecular_weight": chemical_data.get("molecularWeight")
            },
            "warnings": []
        }

        # Extract properties
        if kwargs.get("include_properties", True):
            result["properties"] = self._extract_properties(chemical_data)

        # Fetch and extract toxicity data
        if kwargs.get("include_toxicity", True):
            toxicity_data = await self._get_toxicity_data(dtxsid)
            await asyncio.sleep(0.1)  # Rate limit

            if toxicity_data:
                result["toxicity"] = self._extract_toxicity(toxicity_data)
            else:
                result["toxicity"] = {"predicted_endpoints": [], "qsar_predictions": {}}
                result["warnings"].append("Toxicity data unavailable")

        # Fetch and extract bioactivity data
        if kwargs.get("include_bioactivity", True):
            bioactivity_data = await self._get_bioactivity_data(dtxsid)
            await asyncio.sleep(0.1)

            if bioactivity_data:
                result["bioactivity"] = self._extract_bioactivity(bioactivity_data)
            else:
                result["bioactivity"] = {"num_assays": 0, "active_assays": 0, "gene_targets": []}
                result["warnings"].append("Bioactivity data unavailable")

        # Extract exposure info
        if kwargs.get("include_exposure", True):
            result["exposure"] = self._extract_exposure(chemical_data)

        # Extract hazard info
        if kwargs.get("include_hazards", True):
            result["hazard"] = self._extract_hazard(chemical_data)

        return AdapterResult(
            success=True,
            data=result,
            cache_hit=False,
            metadata={
                "query": query,
                "query_type": query_type,
                "dtxsid": dtxsid,
                "source": "epa_comptox",
                "adapter_version": self.version
            }
        )
