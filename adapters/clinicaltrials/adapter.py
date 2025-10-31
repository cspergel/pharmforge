"""
ClinicalTrials.gov Adapter - Search and analyze clinical trials
API Documentation: https://clinicaltrials.gov/api/v2/studies
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ClinicalTrialsAdapter(AdapterProtocol):
    """
    Adapter for ClinicalTrials.gov API v2
    Searches for clinical trials, study details, outcomes, and results
    """

    BASE_URL = "https://clinicaltrials.gov/api/v2"

    def __init__(self):
        super().__init__(
            name="clinicaltrials",
            adapter_type="api",
            config={
                "rate_limit_delay": 1.0,  # Be respectful to the API
                "timeout": 60,
                "max_studies": 50,
                "page_size": 20
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid search query

        Args:
            input_data: Dictionary with search parameters or drug name string

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            # Must have at least one search parameter
            valid_keys = {"drug", "condition", "nct_id", "term", "phase", "status"}
            return any(key in input_data for key in valid_keys)
        return False

    async def _search_studies(self, query_params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Search ClinicalTrials.gov for studies

        Args:
            query_params: Dictionary of search parameters

        Returns:
            API response with study list or None on error
        """
        url = f"{self.BASE_URL}/studies"

        # Build query parameters
        params = {
            "format": "json",
            "pageSize": self.config.get("page_size", 20)
        }

        # Add search filters
        filter_parts = []

        if "drug" in query_params:
            filter_parts.append(f"AREA[InterventionName]{query_params['drug']}")
        if "condition" in query_params:
            filter_parts.append(f"AREA[Condition]{query_params['condition']}")
        if "phase" in query_params:
            filter_parts.append(f"AREA[Phase]{query_params['phase']}")
        if "status" in query_params:
            filter_parts.append(f"AREA[OverallStatus]{query_params['status']}")
        if "nct_id" in query_params:
            # Direct NCT ID search
            filter_parts.append(f"AREA[NCTId]{query_params['nct_id']}")
        if "term" in query_params:
            # General term search
            filter_parts.append(f"SEARCH[Study]{query_params['term']}")

        if filter_parts:
            params["query.term"] = " AND ".join(filter_parts)

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.warning(f"ClinicalTrials.gov search failed: {response.status}")
                        text = await response.text()
                        logger.debug(f"Response: {text[:500]}")
                        return None

                    data = await response.json()
                    return data

        except asyncio.TimeoutError:
            logger.error(f"ClinicalTrials.gov: Timeout during search")
            return None
        except Exception as e:
            logger.error(f"ClinicalTrials.gov: Error during search: {e}")
            return None

    async def _get_study_details(self, nct_id: str) -> Optional[Dict[str, Any]]:
        """
        Get detailed information for a specific study

        Args:
            nct_id: NCT identifier (e.g., "NCT12345678")

        Returns:
            Detailed study information or None on error
        """
        url = f"{self.BASE_URL}/studies/{nct_id}"
        params = {"format": "json"}

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.warning(f"ClinicalTrials.gov study detail failed: {response.status}")
                        return None

                    data = await response.json()
                    return data

        except Exception as e:
            logger.error(f"ClinicalTrials.gov: Error fetching study {nct_id}: {e}")
            return None

    def _extract_study_info(self, study: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract key information from a study record

        Args:
            study: Raw study data from API

        Returns:
            Cleaned and structured study information
        """
        protocol_section = study.get("protocolSection", {})
        id_module = protocol_section.get("identificationModule", {})
        status_module = protocol_section.get("statusModule", {})
        design_module = protocol_section.get("designModule", {})
        arms_module = protocol_section.get("armsInterventionsModule", {})
        outcomes_module = protocol_section.get("outcomesModule", {})

        # Extract basic info
        info = {
            "nct_id": id_module.get("nctId"),
            "title": id_module.get("briefTitle"),
            "status": status_module.get("overallStatus"),
            "phase": design_module.get("phases", []),
            "enrollment": status_module.get("enrollmentInfo", {}).get("count"),
            "start_date": status_module.get("startDateStruct", {}).get("date"),
            "completion_date": status_module.get("completionDateStruct", {}).get("date"),
        }

        # Extract interventions (drugs)
        interventions = arms_module.get("interventions", [])
        info["interventions"] = [
            {
                "type": i.get("type"),
                "name": i.get("name"),
                "description": i.get("description")
            }
            for i in interventions
        ]

        # Extract conditions
        conditions_module = protocol_section.get("conditionsModule", {})
        info["conditions"] = conditions_module.get("conditions", [])

        # Extract primary outcomes
        primary_outcomes = outcomes_module.get("primaryOutcomes", [])
        info["primary_outcomes"] = [
            {
                "measure": o.get("measure"),
                "description": o.get("description"),
                "timeFrame": o.get("timeFrame")
            }
            for o in primary_outcomes
        ]

        # Extract results if available
        results_section = study.get("resultsSection", {})
        if results_section:
            info["has_results"] = True
            outcome_measures = results_section.get("outcomeMeasuresModule", {})
            info["outcome_measures"] = outcome_measures.get("outcomeMeasures", [])
        else:
            info["has_results"] = False

        return info

    def _summarize_studies(self, studies: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Generate summary statistics across multiple studies

        Args:
            studies: List of processed study records

        Returns:
            Summary statistics
        """
        if not studies:
            return {
                "total_studies": 0,
                "by_phase": {},
                "by_status": {},
                "total_enrollment": 0
            }

        # Count by phase
        phase_counts = {}
        status_counts = {}
        total_enrollment = 0
        has_results_count = 0

        for study in studies:
            # Count phases
            phases = study.get("phase", [])
            if isinstance(phases, list):
                for phase in phases:
                    phase_counts[phase] = phase_counts.get(phase, 0) + 1

            # Count status
            status = study.get("status", "Unknown")
            status_counts[status] = status_counts.get(status, 0) + 1

            # Sum enrollment
            enrollment = study.get("enrollment")
            if enrollment:
                total_enrollment += enrollment

            # Count results
            if study.get("has_results"):
                has_results_count += 1

        return {
            "total_studies": len(studies),
            "by_phase": phase_counts,
            "by_status": status_counts,
            "total_enrollment": total_enrollment,
            "studies_with_results": has_results_count
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute ClinicalTrials.gov search

        Args:
            input_data: Drug name string or dict with search parameters
                       Examples:
                       - "aspirin"
                       - {"drug": "aspirin", "condition": "heart disease"}
                       - {"nct_id": "NCT12345678"}
            **kwargs: Additional parameters:
                     - get_details: bool - Fetch full details for each study
                     - max_studies: int - Maximum number of studies to return

        Returns:
            AdapterResult containing trial information
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a drug name or search parameters dict"
            )

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 1.0)
        await asyncio.sleep(rate_delay)

        # Convert string input to dict
        if isinstance(input_data, str):
            query_params = {"drug": input_data}
        else:
            query_params = input_data

        # Search for studies
        search_results = await self._search_studies(query_params)

        if search_results is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search ClinicalTrials.gov",
                metadata={"query": query_params}
            )

        # Extract studies from response
        studies_raw = search_results.get("studies", [])

        if not studies_raw:
            return AdapterResult(
                success=True,
                data={
                    "studies": [],
                    "summary": self._summarize_studies([]),
                    "total_count": 0
                },
                metadata={"query": query_params, "source": "clinicaltrials.gov"}
            )

        # Process each study
        studies = []
        max_studies = kwargs.get("max_studies", self.config.get("max_studies", 50))

        for study_raw in studies_raw[:max_studies]:
            study_info = self._extract_study_info(study_raw)
            studies.append(study_info)

        # Generate summary
        summary = self._summarize_studies(studies)

        result_data = {
            "studies": studies,
            "summary": summary,
            "total_count": search_results.get("totalCount", len(studies))
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "query": query_params,
                "source": "clinicaltrials.gov",
                "adapter_version": self.version,
                "api_version": "v2"
            }
        )
