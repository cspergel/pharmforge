"""
ChEMBL Adapter - Fetches bioactivity data from ChEMBL REST API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ChEMBLAdapter(AdapterProtocol):
    """
    Adapter for ChEMBL REST API
    Searches for known bioactivity data for compounds
    """

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self):
        super().__init__(
            name="chembl",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,  # 2 requests/second (be nice to ChEMBL)
                "timeout": 60,
                "max_activities": 100
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string or InChIKey

        Args:
            input_data: SMILES string or InChIKey to validate

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        return True

    async def _search_by_smiles_async(self, smiles: str) -> Optional[List[Dict[str, Any]]]:
        """
        Search ChEMBL for bioactivity data by SMILES

        Args:
            smiles: SMILES string of the molecule

        Returns:
            List of bioactivity records or None if not found
        """
        from urllib.parse import quote

        # First, search for the molecule by SMILES
        encoded_smiles = quote(smiles)
        molecule_url = f"{self.BASE_URL}/molecule.json"
        params = {"molecule_structures__canonical_smiles__flexmatch": smiles, "limit": 1}

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                # Find molecule in ChEMBL
                async with session.get(molecule_url, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.warning(f"ChEMBL molecule search failed: {response.status}")
                        return None

                    data = await response.json()
                    molecules = data.get("molecules", [])

                    if not molecules:
                        logger.info(f"ChEMBL: No molecule found for SMILES: {smiles}")
                        return []

                    chembl_id = molecules[0].get("molecule_chembl_id")
                    logger.info(f"ChEMBL: Found molecule {chembl_id} for {smiles}")

                # Get bioactivity data for this molecule
                activity_url = f"{self.BASE_URL}/activity.json"
                activity_params = {
                    "molecule_chembl_id": chembl_id,
                    "limit": self.config.get("max_activities", 100)
                }

                async with session.get(activity_url, params=activity_params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.warning(f"ChEMBL activity search failed: {response.status}")
                        return None

                    data = await response.json()
                    activities = data.get("activities", [])
                    logger.info(f"ChEMBL: Found {len(activities)} activities for {chembl_id}")

                    return activities

        except asyncio.TimeoutError:
            logger.error(f"ChEMBL: Timeout searching for {smiles}")
            return None
        except Exception as e:
            logger.error(f"ChEMBL: Error searching for {smiles}: {e}")
            return None

    def _summarize_activities(self, activities: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Summarize bioactivity data into useful metrics

        Args:
            activities: List of activity records from ChEMBL

        Returns:
            Dictionary of summarized bioactivity metrics
        """
        if not activities:
            return {
                "num_activities": 0,
                "num_targets": 0,
                "assay_types": [],
                "best_ic50_nm": None,
                "best_ki_nm": None
            }

        # Count unique targets
        targets = set()
        assay_types = set()
        ic50_values = []
        ki_values = []

        for activity in activities:
            if activity.get("target_chembl_id"):
                targets.add(activity["target_chembl_id"])

            if activity.get("assay_type"):
                assay_types.add(activity["assay_type"])

            # Collect binding affinity values
            standard_type = (activity.get("standard_type") or "").upper()
            standard_value = activity.get("standard_value")
            standard_units = (activity.get("standard_units") or "").upper()

            if standard_value and standard_units == "NM":  # nanomolar
                try:
                    value_nm = float(standard_value)
                    if standard_type == "IC50":
                        ic50_values.append(value_nm)
                    elif standard_type in ["KI", "KD"]:
                        ki_values.append(value_nm)
                except (ValueError, TypeError):
                    pass

        return {
            "num_activities": len(activities),
            "num_targets": len(targets),
            "assay_types": list(assay_types),
            "best_ic50_nm": min(ic50_values) if ic50_values else None,
            "best_ki_nm": min(ki_values) if ki_values else None,
            "targets": list(targets)[:10]  # Top 10 targets
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the ChEMBL bioactivity search

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing bioactivity data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string"
            )

        smiles = input_data

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Search ChEMBL
        activities = await self._search_by_smiles_async(smiles)

        if activities is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search ChEMBL",
                metadata={"source": "chembl", "smiles": smiles}
            )

        # Summarize the results
        summary = self._summarize_activities(activities)

        # Include both summary AND raw bioactivities for downstream use
        result_data = {
            **summary,
            "bioactivities": activities  # Include raw data for tests/validation
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "chembl",
                "smiles": smiles,
                "adapter_version": self.version,
                "raw_activities_count": len(activities) if activities else 0
            }
        )
