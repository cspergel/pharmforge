"""
FDA FAERS Adapter - Analyze adverse events from FDA's openFDA API
API Documentation: https://open.fda.gov/apis/drug/event/
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
import math
from collections import Counter
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class FDAFAERSAdapter(AdapterProtocol):
    """
    Adapter for FDA FAERS (Adverse Event Reporting System) via openFDA API
    Searches for adverse events, calculates safety signals, and provides pharmacovigilance data
    """

    BASE_URL = "https://api.fda.gov/drug/event.json"

    def __init__(self):
        super().__init__(
            name="fda_faers",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.3,  # openFDA rate limit: 240 requests/minute
                "timeout": 60,
                "max_results": 1000,
                "limit_per_request": 100
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid drug name or active ingredient

        Args:
            input_data: Drug name string or dict with search parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return "drug" in input_data or "ingredient" in input_data
        return False

    async def _search_adverse_events(
        self,
        drug_name: str,
        search_field: str = "patient.drug.medicinalproduct"
    ) -> Optional[Dict[str, Any]]:
        """
        Search for adverse events associated with a drug

        Args:
            drug_name: Name of the drug or active ingredient
            search_field: Field to search (medicinalproduct, openfda.generic_name, etc.)

        Returns:
            API response with adverse events or None on error
        """
        # Build search query
        search_query = f'{search_field}:"{drug_name}"'

        params = {
            "search": search_query,
            "limit": self.config.get("limit_per_request", 100)
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.BASE_URL, params=params, timeout=timeout) as response:
                    if response.status == 404:
                        logger.info(f"FDA FAERS: No events found for {drug_name}")
                        return {"results": [], "meta": {"results": {"total": 0}}}
                    elif response.status != 200:
                        logger.warning(f"FDA FAERS search failed: {response.status}")
                        text = await response.text()
                        logger.debug(f"Response: {text[:500]}")
                        return None

                    data = await response.json()
                    return data

        except asyncio.TimeoutError:
            logger.error(f"FDA FAERS: Timeout searching for {drug_name}")
            return None
        except Exception as e:
            logger.error(f"FDA FAERS: Error searching for {drug_name}: {e}")
            return None

    async def _count_adverse_events(
        self,
        drug_name: str,
        count_field: str = "patient.reaction.reactionmeddrapt.exact"
    ) -> Optional[Dict[str, Any]]:
        """
        Get count statistics for adverse events

        Args:
            drug_name: Name of the drug
            count_field: Field to count by (reaction type, outcomes, etc.)

        Returns:
            Count statistics or None on error
        """
        search_query = f'patient.drug.medicinalproduct:"{drug_name}"'

        params = {
            "search": search_query,
            "count": count_field
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.BASE_URL, params=params, timeout=timeout) as response:
                    if response.status == 404:
                        return {"results": []}
                    elif response.status != 200:
                        logger.warning(f"FDA FAERS count failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        except Exception as e:
            logger.error(f"FDA FAERS: Error counting events for {drug_name}: {e}")
            return None

    def _extract_event_info(self, event: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract key information from an adverse event report

        Args:
            event: Raw event data from API

        Returns:
            Cleaned and structured event information
        """
        patient = event.get("patient", {})

        # Extract reactions
        reactions = patient.get("reaction", [])
        reaction_list = [
            {
                "term": r.get("reactionmeddrapt"),
                "outcome": r.get("reactionoutcome")
            }
            for r in reactions
        ]

        # Extract drug information
        drugs = patient.get("drug", [])
        drug_list = []
        for drug in drugs:
            drug_list.append({
                "name": drug.get("medicinalproduct"),
                "role": drug.get("drugcharacterization"),  # 1=suspect, 2=concomitant, 3=interacting
                "indication": drug.get("drugindication")
            })

        # Extract demographics
        info = {
            "report_id": event.get("safetyreportid"),
            "receive_date": event.get("receivedate"),
            "serious": event.get("serious", 0),
            "age": patient.get("patientonsetage"),
            "age_unit": patient.get("patientonsetageunit"),
            "sex": patient.get("patientsex"),  # 1=Male, 2=Female
            "weight_kg": patient.get("patientweight"),
            "reactions": reaction_list,
            "drugs": drug_list,
            "outcomes": [patient.get("seriousnessdeath"),
                        patient.get("seriousnesshospitalization"),
                        patient.get("seriousnesslifethreatening"),
                        patient.get("seriousnessdisabling")]
        }

        return info

    def _calculate_safety_metrics(
        self,
        drug_events: int,
        reaction_counts: Dict[str, int],
        total_events: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        Calculate pharmacovigilance safety metrics

        Args:
            drug_events: Total number of events for the drug
            reaction_counts: Dictionary of reaction types and their counts
            total_events: Total events in database (for PRR calculation)

        Returns:
            Dictionary of safety metrics
        """
        if drug_events == 0:
            return {
                "total_events": 0,
                "unique_reactions": 0,
                "top_reactions": [],
                "severity_indicators": {}
            }

        # Top reactions
        top_reactions = sorted(
            reaction_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]

        metrics = {
            "total_events": drug_events,
            "unique_reactions": len(reaction_counts),
            "top_reactions": [
                {"term": term, "count": count, "percentage": (count / drug_events) * 100}
                for term, count in top_reactions
            ]
        }

        # Calculate Proportional Reporting Ratio (PRR) if total events known
        # PRR = (a/b) / (c/d) where:
        # a = events with drug and reaction
        # b = events with drug
        # c = events with reaction (not including drug)
        # d = all events
        if total_events and total_events > drug_events:
            prr_values = {}
            for reaction, count in top_reactions[:5]:  # Top 5 for PRR
                a = count
                b = drug_events
                # Estimate c and d (simplified - real calculation needs more data)
                # This is a placeholder to show the concept
                estimated_reaction_total = count * 10  # Rough estimate
                c = max(estimated_reaction_total - a, 1)
                d = total_events - b

                if c > 0 and d > 0:
                    prr = (a / b) / (c / d)
                    # PRR > 2 with chi-square test often indicates signal
                    prr_values[reaction] = round(prr, 2)

            metrics["prr_signals"] = prr_values

        return metrics

    def _analyze_demographics(self, events: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Analyze demographic distribution of adverse events

        Args:
            events: List of processed event records

        Returns:
            Demographics summary
        """
        sex_counts = Counter()
        age_groups = Counter()
        serious_count = 0

        for event in events:
            # Sex distribution (1=Male, 2=Female)
            sex = event.get("sex")
            if sex == "1":
                sex_counts["Male"] += 1
            elif sex == "2":
                sex_counts["Female"] += 1
            else:
                sex_counts["Unknown"] += 1

            # Age groups (approximate)
            age = event.get("age")
            if age:
                try:
                    age_val = float(age)
                    if age_val < 18:
                        age_groups["<18"] += 1
                    elif age_val < 45:
                        age_groups["18-44"] += 1
                    elif age_val < 65:
                        age_groups["45-64"] += 1
                    else:
                        age_groups["65+"] += 1
                except (ValueError, TypeError):
                    age_groups["Unknown"] += 1

            # Serious events
            if event.get("serious") == 1:
                serious_count += 1

        return {
            "sex_distribution": dict(sex_counts),
            "age_distribution": dict(age_groups),
            "serious_events": serious_count,
            "serious_percentage": (serious_count / len(events) * 100) if events else 0
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute FDA FAERS adverse event search

        Args:
            input_data: Drug name string or dict with search parameters
                       Examples:
                       - "aspirin"
                       - {"drug": "aspirin", "ingredient": "acetylsalicylic acid"}
            **kwargs: Additional parameters:
                     - search_field: Field to search (default: medicinalproduct)
                     - include_demographics: bool - Include demographic analysis
                     - calculate_metrics: bool - Calculate safety metrics

        Returns:
            AdapterResult containing adverse event data and safety metrics
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a drug name or search parameters dict"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.3)
        await asyncio.sleep(rate_delay)

        # Extract drug name
        if isinstance(input_data, str):
            drug_name = input_data
        else:
            drug_name = input_data.get("drug") or input_data.get("ingredient")

        # Search for adverse events
        event_results = await self._search_adverse_events(drug_name)

        if event_results is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search FDA FAERS",
                metadata={"drug": drug_name}
            )

        # Extract events
        events_raw = event_results.get("results", [])
        total_count = event_results.get("meta", {}).get("results", {}).get("total", len(events_raw))

        if not events_raw:
            return AdapterResult(
                success=True,
                data={
                    "events": [],
                    "total_count": 0,
                    "demographics": {},
                    "safety_metrics": self._calculate_safety_metrics(0, {})
                },
                metadata={"drug": drug_name, "source": "fda_faers"}
            )

        # Process events
        events = [self._extract_event_info(event) for event in events_raw]

        # Get reaction count statistics
        reaction_count_results = await self._count_adverse_events(drug_name)
        reaction_counts = {}
        if reaction_count_results:
            for item in reaction_count_results.get("results", []):
                reaction_counts[item.get("term")] = item.get("count", 0)

        # Calculate safety metrics
        safety_metrics = self._calculate_safety_metrics(
            total_count,
            reaction_counts
        )

        # Analyze demographics if requested
        demographics = {}
        if kwargs.get("include_demographics", True):
            demographics = self._analyze_demographics(events)

        result_data = {
            "events": events[:100],  # Limit returned events for performance
            "total_count": total_count,
            "safety_metrics": safety_metrics,
            "demographics": demographics
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "drug": drug_name,
                "source": "fda_faers",
                "adapter_version": self.version,
                "api": "openFDA"
            }
        )
