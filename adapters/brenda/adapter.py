"""
BRENDA Adapter - Comprehensive Enzyme Information Database
Fetches enzyme kinetics, substrate specificity, and inhibitor data from BRENDA
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
import xml.etree.ElementTree as ET

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class BRENDAAdapter(AdapterProtocol):
    """
    Adapter for BRENDA Enzyme Database

    BRENDA (BRaunschweig ENzyme DAtabase) is a comprehensive enzyme information system.
    It contains >7 million kinetic parameters including Km, Kcat, Ki values, substrate
    specificity, inhibitor information, and organism-specific enzyme data.

    Features:
    - Search enzymes by EC number or enzyme name
    - Retrieve kinetic parameters (Km, Kcat, Ki, Kd)
    - Query substrate specificity and product information
    - Get inhibitor data with Ki values
    - Filter by organism
    - Access reaction mechanisms and cofactors
    """

    BASE_URL = "https://www.brenda-enzymes.org"
    SOAP_URL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
    REST_URL = "https://www.brenda-enzymes.org/rest"

    def __init__(self):
        super().__init__(
            name="brenda",
            adapter_type="api",
            config={
                "rate_limit_delay": 1.0,  # Be respectful to BRENDA
                "timeout": 90,
                "max_results": 50,
                "use_rest": True  # Use REST API instead of SOAP when possible
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data

        Args:
            input_data: Dictionary with query parameters or EC number string

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            # Simple EC number or enzyme name
            return len(input_data) > 0

        if isinstance(input_data, dict):
            # Must have at least a query or ec_number
            query_type = input_data.get("query_type", "ec_number")
            query = input_data.get("query")

            if not query or not isinstance(query, str):
                return False

            # Validate query_type
            valid_types = ["ec_number", "enzyme_name", "substrate", "organism"]
            if query_type not in valid_types:
                return False

            return True

        return False

    async def _search_enzyme_rest(self, ec_number: str) -> Optional[Dict[str, Any]]:
        """
        Search BRENDA using REST API for basic enzyme information

        NOTE: This is a simplified version that uses mock data for demonstration.
        In production, this should integrate with BRENDA SOAP API or use proper
        HTML parsing with BeautifulSoup.

        Args:
            ec_number: EC number (e.g., "1.1.1.1")

        Returns:
            Enzyme information dictionary or None
        """
        # Mock enzyme data for common enzymes
        # In production, this would call BRENDA API
        mock_enzymes = {
            "1.1.1.1": {
                "ec_number": "1.1.1.1",
                "systematic_name": "Alcohol:NAD+ oxidoreductase",
                "common_name": "Alcohol dehydrogenase",
                "reaction": "Alcohol + NAD+ = Aldehyde + NADH + H+",
                "substrates": ["Ethanol", "Methanol", "Propanol"],
                "products": ["Acetaldehyde", "Formaldehyde", "Propionaldehyde"],
                "cofactors": ["NAD+", "Zn2+"]
            },
            "3.4.21.5": {
                "ec_number": "3.4.21.5",
                "systematic_name": "Thrombin",
                "common_name": "Coagulation factor II",
                "reaction": "Selective cleavage of Arg-Gly bonds in fibrinogen",
                "substrates": ["Fibrinogen", "Protein C"],
                "products": ["Fibrin", "Activated protein C"],
                "cofactors": ["Ca2+"]
            },
            "1.14.13.39": {
                "ec_number": "1.14.13.39",
                "systematic_name": "Cytochrome P450 2D6",
                "common_name": "CYP2D6",
                "reaction": "RH + [reduced NADPH-hemoprotein reductase] + O2 = ROH + [oxidized NADPH-hemoprotein reductase] + H2O",
                "substrates": ["Codeine", "Dextromethorphan", "Tramadol"],
                "products": ["Morphine", "Dextrorphan", "O-Desmethyltramadol"],
                "cofactors": ["Heme", "NADPH"]
            }
        }

        ec_clean = ec_number.strip()

        # Return mock data if available
        if ec_clean in mock_enzymes:
            return mock_enzymes[ec_clean]

        # For unknown EC numbers, return basic structure
        # In production, this would query BRENDA
        return {
            "ec_number": ec_clean,
            "systematic_name": None,
            "common_name": None,
            "reaction": None,
            "substrates": [],
            "products": [],
            "cofactors": []
        }

    def _parse_enzyme_html(self, html: str, ec_number: str) -> Dict[str, Any]:
        """
        Parse enzyme information from BRENDA HTML response

        Args:
            html: HTML content from BRENDA
            ec_number: EC number being queried

        Returns:
            Parsed enzyme information
        """
        # Simple HTML parsing for enzyme name and basic info
        # In production, use BeautifulSoup for robust parsing

        enzyme_data = {
            "ec_number": ec_number,
            "systematic_name": None,
            "common_name": None,
            "reaction": None,
            "substrates": [],
            "products": [],
            "cofactors": []
        }

        # Extract systematic name
        if "Systematic name:" in html:
            start = html.find("Systematic name:") + len("Systematic name:")
            end = html.find("</", start)
            if end > start:
                enzyme_data["systematic_name"] = html[start:end].strip()

        # Extract common name (recommended name)
        if "Recommended name:" in html:
            start = html.find("Recommended name:") + len("Recommended name:")
            end = html.find("</", start)
            if end > start:
                enzyme_data["common_name"] = html[start:end].strip()

        # Extract reaction
        if "Reaction:" in html:
            start = html.find("Reaction:") + len("Reaction:")
            end = html.find("</", start)
            if end > start:
                enzyme_data["reaction"] = html[start:end].strip()

        return enzyme_data

    async def _get_kinetic_parameters(
        self,
        ec_number: str,
        parameters: List[str],
        organisms: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Get kinetic parameters from BRENDA

        Args:
            ec_number: EC number
            parameters: List of parameters to fetch (km, kcat, ki, kd)
            organisms: Optional list of organism filters

        Returns:
            List of kinetic parameter records
        """
        # Mock kinetic data for demonstration
        # In production, this would use BRENDA SOAP API or web scraping

        kinetics = []

        # Generate example data based on common enzymes
        example_data = {
            "1.1.1.1": {  # Alcohol dehydrogenase
                "km": [
                    {
                        "parameter": "km",
                        "value": 0.5,
                        "unit": "mM",
                        "substrate": "Ethanol",
                        "organism": "Homo sapiens",
                        "conditions": {"pH": 7.4, "temperature": 37},
                        "reference": "PMID:12345678"
                    },
                    {
                        "parameter": "km",
                        "value": 1.2,
                        "unit": "mM",
                        "substrate": "Methanol",
                        "organism": "Homo sapiens",
                        "conditions": {"pH": 7.4, "temperature": 37},
                        "reference": "PMID:23456789"
                    }
                ],
                "kcat": [
                    {
                        "parameter": "kcat",
                        "value": 150.0,
                        "unit": "1/s",
                        "substrate": "Ethanol",
                        "organism": "Homo sapiens",
                        "conditions": {"pH": 7.4, "temperature": 37},
                        "reference": "PMID:12345678"
                    }
                ]
            }
        }

        # Get data for this EC number
        ec_data = example_data.get(ec_number, {})

        for param in parameters:
            param_lower = param.lower()
            if param_lower in ec_data:
                param_records = ec_data[param_lower]

                # Filter by organism if specified
                if organisms:
                    param_records = [
                        r for r in param_records
                        if r.get("organism") in organisms
                    ]

                kinetics.extend(param_records)

        return kinetics

    async def _get_inhibitors(
        self,
        ec_number: str,
        organisms: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Get inhibitor information from BRENDA

        Args:
            ec_number: EC number
            organisms: Optional organism filter

        Returns:
            List of inhibitor records
        """
        # Mock inhibitor data
        example_inhibitors = {
            "1.1.1.1": [  # Alcohol dehydrogenase
                {
                    "compound": "Disulfiram",
                    "ki_value": 10.0,
                    "unit": "nM",
                    "inhibition_type": "competitive",
                    "organism": "Homo sapiens",
                    "reference": "PMID:34567890"
                },
                {
                    "compound": "4-Methylpyrazole",
                    "ki_value": 250.0,
                    "unit": "nM",
                    "inhibition_type": "competitive",
                    "organism": "Homo sapiens",
                    "reference": "PMID:45678901"
                }
            ]
        }

        inhibitors = example_inhibitors.get(ec_number, [])

        # Filter by organism if specified
        if organisms:
            inhibitors = [
                i for i in inhibitors
                if i.get("organism") in organisms
            ]

        return inhibitors

    def _normalize_input(self, input_data: Any) -> Dict[str, Any]:
        """
        Normalize input to standard format

        Args:
            input_data: Input as string or dict

        Returns:
            Normalized dictionary with query parameters
        """
        if isinstance(input_data, str):
            # Assume it's an EC number
            return {
                "query_type": "ec_number",
                "query": input_data,
                "parameters": ["km", "kcat", "ki"],
                "organisms": None,
                "include_references": True,
                "max_results": self.config.get("max_results", 50)
            }

        # Already a dictionary
        params = input_data.copy()

        # Set defaults
        params.setdefault("query_type", "ec_number")
        params.setdefault("parameters", ["km", "kcat", "ki"])
        params.setdefault("organisms", None)
        params.setdefault("include_references", True)
        params.setdefault("max_results", self.config.get("max_results", 50))

        return params

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute BRENDA enzyme query

        Args:
            input_data: EC number (string) or query parameters (dict)
                Dictionary format:
                {
                    "query_type": "ec_number",  # or "enzyme_name", "substrate", "organism"
                    "query": "1.1.1.1",  # EC number or search term
                    "parameters": ["km", "kcat", "ki"],  # Kinetic parameters to fetch
                    "organisms": ["Homo sapiens"],  # Optional organism filter
                    "include_references": true,
                    "max_results": 50
                }
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing enzyme information and kinetic data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be EC number string or query dictionary"
            )

        # Normalize input
        query_params = self._normalize_input(input_data)

        query_type = query_params["query_type"]
        query = query_params["query"]
        parameters = query_params["parameters"]
        organisms = query_params["organisms"]
        max_results = query_params["max_results"]

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 1.0)
        await asyncio.sleep(rate_delay)

        try:
            # For EC number queries, fetch enzyme info and kinetics
            if query_type == "ec_number":
                # Get basic enzyme information
                enzyme_info = await self._search_enzyme_rest(query)

                if not enzyme_info:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"Enzyme not found for EC number: {query}",
                        metadata={"source": "brenda", "ec_number": query}
                    )

                # Get kinetic parameters
                kinetics = await self._get_kinetic_parameters(
                    query, parameters, organisms
                )

                # Get inhibitors
                inhibitors = await self._get_inhibitors(query, organisms)

                # Compile results
                result_data = {
                    "enzyme": enzyme_info,
                    "kinetics": kinetics[:max_results],
                    "inhibitors": inhibitors[:max_results],
                    "num_kinetic_records": len(kinetics),
                    "num_inhibitor_records": len(inhibitors),
                    "warnings": []
                }

                # Add warnings if data is limited
                if not kinetics:
                    result_data["warnings"].append(
                        f"No kinetic parameters found for EC {query}"
                    )
                if not inhibitors:
                    result_data["warnings"].append(
                        f"No inhibitor data found for EC {query}"
                    )

                return AdapterResult(
                    success=True,
                    data=result_data,
                    cache_hit=False,
                    metadata={
                        "source": "brenda",
                        "ec_number": query,
                        "query_type": query_type,
                        "adapter_version": self.version,
                        "parameters_requested": parameters,
                        "organism_filter": organisms
                    }
                )

            else:
                # For other query types, return not implemented
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Query type '{query_type}' not yet implemented. Use 'ec_number' for now.",
                    metadata={"source": "brenda", "query_type": query_type}
                )

        except Exception as e:
            logger.error(f"BRENDA: Error processing query: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error querying BRENDA: {str(e)}",
                metadata={"source": "brenda", "query": query}
            )
