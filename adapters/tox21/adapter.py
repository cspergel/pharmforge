"""
Tox21 Adapter for PharmForge

Provides access to NIH Toxicology in the 21st Century (Tox21) screening data.
Tox21 is a collaboration between NIH, EPA, and FDA to develop predictive models
for toxicity screening using high-throughput assays.

This adapter queries Tox21 data via PubChem's BioAssay database, which hosts
the complete Tox21 screening data.

Key Features:
    - 12,000+ compounds tested across 50+ toxicity assays
    - Nuclear receptor binding assays
    - Stress response pathway assays
    - Cytotoxicity and genotoxicity endpoints
    - AC50 values (activity concentrations)
    - Public domain data

References:
    - Tox21: https://tripod.nih.gov/tox21/
    - PubChem BioAssay: https://pubchem.ncbi.nlm.nih.gov/
    - Tox21 Data: https://doi.org/10.1093/toxsci/kfy131
"""

import aiohttp
import asyncio
import logging
from typing import Dict, Any, List, Optional
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class Tox21Adapter(AdapterProtocol):
    """
    Tox21 adapter for toxicity screening data.

    Queries Tox21 data via PubChem BioAssay API to retrieve toxicity screening
    results for compounds. Tox21 focuses on:

    Nuclear Receptors (NR):
        - AR: Androgen Receptor (agonist/antagonist)
        - ER: Estrogen Receptor Alpha (agonist/antagonist)
        - ER-LBD: Estrogen Receptor Ligand Binding Domain
        - AR-LBD: Androgen Receptor Ligand Binding Domain
        - AhR: Aryl Hydrocarbon Receptor
        - Aromatase: CYP19A1 inhibition
        - PPAR-gamma: Peroxisome Proliferator-Activated Receptor Gamma

    Stress Response (SR):
        - ARE: Antioxidant Response Element
        - ATAD5: DNA damage response
        - HSE: Heat Shock Response Element
        - MMP: Mitochondrial Membrane Potential
        - p53: Tumor suppressor activation

    Other Endpoints:
        - Cell viability (cytotoxicity)
        - Genotoxicity
        - Mitochondrial toxicity

    Data includes activity calls (active/inactive), AC50 values, efficacy,
    and curve quality metrics.
    """

    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    # Major Tox21 assay categories with PubChem AID mappings
    TOX21_ASSAYS = {
        # Nuclear Receptors
        "NR-AR": {
            "name": "Androgen Receptor",
            "aids": [743040, 743053],  # Agonist, Antagonist
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor"
        },
        "NR-ER": {
            "name": "Estrogen Receptor Alpha",
            "aids": [743077, 743079],  # Agonist, Antagonist
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor"
        },
        "NR-ER-LBD": {
            "name": "Estrogen Receptor Ligand Binding Domain",
            "aids": [743078, 743080],
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor"
        },
        "NR-AR-LBD": {
            "name": "Androgen Receptor Ligand Binding Domain",
            "aids": [743054, 743055],
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor"
        },
        "NR-AhR": {
            "name": "Aryl Hydrocarbon Receptor",
            "aids": [743122],
            "endpoint": "metabolic_activation",
            "category": "nuclear_receptor"
        },
        "NR-Aromatase": {
            "name": "Aromatase (CYP19A1) Inhibition",
            "aids": [743139],
            "endpoint": "endocrine_disruption",
            "category": "nuclear_receptor"
        },
        "NR-PPAR-gamma": {
            "name": "PPAR Gamma",
            "aids": [743140],
            "endpoint": "metabolic_disruption",
            "category": "nuclear_receptor"
        },

        # Stress Response
        "SR-ARE": {
            "name": "Antioxidant Response Element",
            "aids": [743219],
            "endpoint": "oxidative_stress",
            "category": "stress_response"
        },
        "SR-ATAD5": {
            "name": "DNA Damage Response (ATAD5)",
            "aids": [720516],
            "endpoint": "genotoxicity",
            "category": "stress_response"
        },
        "SR-HSE": {
            "name": "Heat Shock Response Element",
            "aids": [743225],
            "endpoint": "stress_response",
            "category": "stress_response"
        },
        "SR-MMP": {
            "name": "Mitochondrial Membrane Potential",
            "aids": [743254],
            "endpoint": "mitochondrial_toxicity",
            "category": "stress_response"
        },
        "SR-p53": {
            "name": "p53 Tumor Suppressor Activation",
            "aids": [743212],
            "endpoint": "dna_damage",
            "category": "stress_response"
        }
    }

    def __init__(self):
        super().__init__(
            name="tox21",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.3,  # PubChem rate limit: ~5 req/sec
                "timeout": 45,
                "max_retries": 3
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data for Tox21 queries.

        Accepts:
            - SMILES strings
            - PubChem CIDs (as integers or strings)
            - Dictionary with query_type and query fields

        Args:
            input_data: Input to validate

        Returns:
            True if valid, False otherwise
        """
        # Dictionary input with query details
        if isinstance(input_data, dict):
            if "query" not in input_data:
                return False
            query_type = input_data.get("query_type", "smiles")
            query = input_data.get("query")

            if query_type == "smiles":
                return self._validate_smiles(query)
            elif query_type == "cid":
                return isinstance(query, (int, str)) and str(query).isdigit()
            elif query_type in ["cas", "name"]:
                return isinstance(query, str) and len(query) > 0
            return False

        # Simple string input (assume SMILES)
        if isinstance(input_data, str):
            return self._validate_smiles(input_data)

        # Integer input (assume CID)
        if isinstance(input_data, int):
            return True

        return False

    def _validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES string."""
        if not smiles or not isinstance(smiles, str):
            return False
        if len(smiles) == 0:
            return False
        # Basic SMILES character validation
        valid_chars = set("CNOPSFClBrI[]()=@#+-0123456789cnops/\\")
        return all(c in valid_chars or c.isspace() for c in smiles)

    async def _get_cid_from_smiles(self, session: aiohttp.ClientSession, smiles: str) -> Optional[int]:
        """
        Convert SMILES to PubChem CID.

        Args:
            session: aiohttp session
            smiles: SMILES string

        Returns:
            PubChem CID or None if not found
        """
        try:
            encoded_smiles = quote(smiles)
            url = f"{self.BASE_URL}/compound/smiles/{encoded_smiles}/cids/JSON"

            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))
            async with session.get(url, timeout=timeout) as response:
                if response.status == 200:
                    data = await response.json()
                    cids = data.get("IdentifierList", {}).get("CID", [])
                    if cids:
                        cid = cids[0]
                        logger.info(f"Tox21: Converted SMILES to CID {cid}")
                        return cid
                else:
                    logger.warning(f"Tox21: Could not convert SMILES to CID (status {response.status})")
                    return None
        except Exception as e:
            logger.error(f"Tox21: Error converting SMILES to CID: {e}")
            return None

    async def _get_cid_from_name(self, session: aiohttp.ClientSession, name: str) -> Optional[int]:
        """
        Convert compound name to PubChem CID.

        Args:
            session: aiohttp session
            name: Compound name or CAS number

        Returns:
            PubChem CID or None if not found
        """
        try:
            encoded_name = quote(name)
            url = f"{self.BASE_URL}/compound/name/{encoded_name}/cids/JSON"

            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))
            async with session.get(url, timeout=timeout) as response:
                if response.status == 200:
                    data = await response.json()
                    cids = data.get("IdentifierList", {}).get("CID", [])
                    if cids:
                        cid = cids[0]
                        logger.info(f"Tox21: Converted name '{name}' to CID {cid}")
                        return cid
                else:
                    logger.warning(f"Tox21: Could not convert name to CID (status {response.status})")
                    return None
        except Exception as e:
            logger.error(f"Tox21: Error converting name to CID: {e}")
            return None

    async def _query_bioassay(
        self,
        session: aiohttp.ClientSession,
        cid: int,
        aid: int
    ) -> Optional[Dict[str, Any]]:
        """
        Query a specific Tox21 bioassay for a compound.

        Args:
            session: aiohttp session
            cid: PubChem Compound ID
            aid: PubChem Assay ID

        Returns:
            Assay results dictionary or None
        """
        try:
            # Query bioassay data
            url = f"{self.BASE_URL}/assay/aid/{aid}/cid/{cid}/JSON"

            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))
            async with session.get(url, timeout=timeout) as response:
                if response.status == 200:
                    data = await response.json()

                    # Extract assay data
                    tables = data.get("PC_AssaySubmit", [{}])[0].get("data", [])
                    if not tables:
                        return None

                    # Parse activity data
                    result = {
                        "aid": aid,
                        "cid": cid,
                        "activity": "inactive",
                        "ac50": None,
                        "efficacy": None,
                        "curve_class": None
                    }

                    # Look for activity outcome
                    for table in tables:
                        tid_name = table.get("tid", {}).get("name", "")
                        value = table.get("value", {})

                        # Activity call
                        if "Phenotype" in tid_name or "Activity" in tid_name:
                            if isinstance(value, dict):
                                outcome = value.get("sval", "").lower()
                                if "active" in outcome or "agonist" in outcome or "antagonist" in outcome:
                                    result["activity"] = "active"
                            elif isinstance(value, str):
                                if "active" in value.lower():
                                    result["activity"] = "active"

                        # AC50 value
                        if "AC50" in tid_name or "ac50" in tid_name.lower():
                            if isinstance(value, dict):
                                result["ac50"] = value.get("fval")
                            elif isinstance(value, (int, float)):
                                result["ac50"] = float(value)

                        # Efficacy
                        if "Efficacy" in tid_name or "efficacy" in tid_name.lower():
                            if isinstance(value, dict):
                                result["efficacy"] = value.get("fval")
                            elif isinstance(value, (int, float)):
                                result["efficacy"] = float(value)

                        # Curve class
                        if "Curve" in tid_name and "Class" in tid_name:
                            if isinstance(value, dict):
                                result["curve_class"] = value.get("fval")
                            elif isinstance(value, (int, float)):
                                result["curve_class"] = float(value)

                    return result

                elif response.status == 404:
                    # Compound not tested in this assay
                    return None
                else:
                    logger.warning(f"Tox21: Error querying AID {aid} for CID {cid} (status {response.status})")
                    return None

        except Exception as e:
            logger.error(f"Tox21: Error querying bioassay AID {aid}: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Tox21 toxicity screening query.

        Args:
            input_data: SMILES string, CID, or dict with query details
            **kwargs: Additional parameters:
                - assays: List of assay IDs to query (default: all)
                - include_inactive: Include inactive results (default: True)
                - ac50_threshold: Filter by AC50 threshold in uM (default: None)

        Returns:
            AdapterResult containing Tox21 screening data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input data for Tox21 query"
            )

        # Parse input
        if isinstance(input_data, dict):
            query_type = input_data.get("query_type", "smiles")
            query = input_data.get("query")
            assay_filter = input_data.get("assays", None)
            include_inactive = input_data.get("include_inactive", True)
            ac50_threshold = input_data.get("ac50_threshold", None)
        else:
            query = input_data
            query_type = "cid" if isinstance(input_data, int) else "smiles"
            assay_filter = kwargs.get("assays", None)
            include_inactive = kwargs.get("include_inactive", True)
            ac50_threshold = kwargs.get("ac50_threshold", None)

        try:
            connector = aiohttp.TCPConnector(limit=10)
            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 45))

            async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
                # Get CID
                cid = None
                if query_type == "cid":
                    cid = int(query)
                elif query_type == "smiles":
                    await asyncio.sleep(self.config.get("rate_limit_delay", 0.3))
                    cid = await self._get_cid_from_smiles(session, query)
                elif query_type in ["name", "cas"]:
                    await asyncio.sleep(self.config.get("rate_limit_delay", 0.3))
                    cid = await self._get_cid_from_name(session, query)

                if not cid:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"Could not find PubChem CID for query: {query}"
                    )

                logger.info(f"Tox21: Querying toxicity data for CID {cid}")

                # Determine which assays to query
                if assay_filter and isinstance(assay_filter, list):
                    assays_to_query = {k: v for k, v in self.TOX21_ASSAYS.items() if k in assay_filter}
                else:
                    assays_to_query = self.TOX21_ASSAYS

                # Query all assays
                toxicity_results = []
                active_count = 0
                inactive_count = 0
                toxicity_flags = set()

                for assay_id, assay_info in assays_to_query.items():
                    for aid in assay_info["aids"]:
                        await asyncio.sleep(self.config.get("rate_limit_delay", 0.3))

                        result = await self._query_bioassay(session, cid, aid)

                        if result:
                            # Apply filters
                            if not include_inactive and result["activity"] == "inactive":
                                continue

                            if ac50_threshold and result.get("ac50"):
                                if result["ac50"] > ac50_threshold:
                                    continue

                            # Determine confidence based on curve class
                            confidence = "unknown"
                            if result.get("curve_class"):
                                cc = result["curve_class"]
                                if cc in [1.1, 1.2, 2.1, 2.2]:
                                    confidence = "high"
                                elif cc in [3, 4]:
                                    confidence = "medium"
                                else:
                                    confidence = "low"

                            assay_result = {
                                "assay_id": assay_id,
                                "assay_name": assay_info["name"],
                                "endpoint": assay_info["endpoint"],
                                "category": assay_info["category"],
                                "pubchem_aid": aid,
                                "activity": result["activity"],
                                "ac50": result.get("ac50"),
                                "efficacy": result.get("efficacy"),
                                "curve_class": result.get("curve_class"),
                                "confidence": confidence
                            }

                            toxicity_results.append(assay_result)

                            # Count activities
                            if result["activity"] == "active":
                                active_count += 1
                                toxicity_flags.add(assay_info["endpoint"])
                            else:
                                inactive_count += 1

                # Calculate overall risk
                total_assays = len(toxicity_results)
                if total_assays == 0:
                    overall_risk = "unknown"
                elif active_count == 0:
                    overall_risk = "low"
                elif active_count / total_assays < 0.2:
                    overall_risk = "low"
                elif active_count / total_assays < 0.5:
                    overall_risk = "medium"
                else:
                    overall_risk = "high"

                # Build result
                result_data = {
                    "compound": {
                        "cid": cid,
                        "query": query,
                        "query_type": query_type
                    },
                    "toxicity_results": toxicity_results,
                    "summary": {
                        "total_assays": total_assays,
                        "active_assays": active_count,
                        "inactive_assays": inactive_count,
                        "toxicity_flags": list(toxicity_flags),
                        "overall_risk": overall_risk
                    },
                    "warnings": []
                }

                # Add warnings
                if "endocrine_disruption" in toxicity_flags:
                    result_data["warnings"].append("Potential endocrine disruption activity detected")
                if "genotoxicity" in toxicity_flags:
                    result_data["warnings"].append("Potential genotoxicity detected")
                if "mitochondrial_toxicity" in toxicity_flags:
                    result_data["warnings"].append("Potential mitochondrial toxicity detected")

                logger.info(f"Tox21: Query complete - {active_count} active / {total_assays} total assays")

                return AdapterResult(
                    success=True,
                    data=result_data,
                    metadata={
                        "adapter_name": self.name,
                        "version": self.version,
                        "source": "Tox21 via PubChem BioAssay",
                        "cid": cid,
                        "assays_queried": list(assays_to_query.keys())
                    }
                )

        except Exception as e:
            logger.error(f"Tox21 query failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=f"Tox21 query failed: {str(e)}"
            )

    def get_assay_info(self, assay_id: str) -> Optional[Dict[str, Any]]:
        """
        Get information about a specific Tox21 assay.

        Args:
            assay_id: Tox21 assay identifier (e.g., "NR-AR", "SR-HSE")

        Returns:
            Assay information dictionary or None
        """
        return self.TOX21_ASSAYS.get(assay_id)

    def list_assays(self, category: Optional[str] = None) -> List[str]:
        """
        List available Tox21 assays.

        Args:
            category: Filter by category ("nuclear_receptor" or "stress_response")

        Returns:
            List of assay IDs
        """
        if category:
            return [
                assay_id for assay_id, info in self.TOX21_ASSAYS.items()
                if info["category"] == category
            ]
        return list(self.TOX21_ASSAYS.keys())

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata.

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "NIH Tox21 toxicity screening data via PubChem BioAssay",
            "data_source": "Tox21 (NIH/EPA/FDA collaboration)",
            "coverage": {
                "compounds": "12,000+",
                "assays": "50+",
                "categories": ["nuclear_receptor", "stress_response"]
            },
            "assays": {
                "total": len(self.TOX21_ASSAYS),
                "nuclear_receptor": len(self.list_assays("nuclear_receptor")),
                "stress_response": len(self.list_assays("stress_response")),
                "available": list(self.TOX21_ASSAYS.keys())
            },
            "endpoints": [
                "endocrine_disruption",
                "metabolic_activation",
                "metabolic_disruption",
                "oxidative_stress",
                "genotoxicity",
                "stress_response",
                "mitochondrial_toxicity",
                "dna_damage"
            ],
            "config": self.config,
            "references": [
                "https://tripod.nih.gov/tox21/",
                "https://pubchem.ncbi.nlm.nih.gov/",
                "Huang et al. (2018) Nucleic Acids Res. doi:10.1093/nar/gky1033"
            ],
            "citation": "Huang, R., et al. (2018). The NCATS Tox21 Program: Enabling technologies for high-throughput screening and data integration. Frontiers in Public Health, 6, 190."
        }
