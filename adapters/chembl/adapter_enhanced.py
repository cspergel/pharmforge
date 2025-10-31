"""
Enhanced ChEMBL Adapter - Complete ChEMBL API integration
Adds target queries, drug info, mechanism of action, and advanced filtering
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ChEMBLEnhancedAdapter(AdapterProtocol):
    """
    Enhanced adapter for ChEMBL REST API

    Capabilities:
    - Bioactivity data search (existing)
    - Target information and queries (new)
    - Drug/clinical candidate lookup (new)
    - Mechanism of action data (new)
    - Similarity search (new)
    - Substructure search (new)
    - Advanced activity filtering (new)
    """

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    def __init__(self):
        super().__init__(
            name="chembl_enhanced",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,
                "timeout": 90,
                "max_activities": 200,
                "similarity_threshold": 70  # 70% similarity
            }
        )
        self.version = "2.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input

        Args:
            input_data: SMILES, ChEMBL ID, or dict

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data) > 0
        elif isinstance(input_data, dict):
            return any(k in input_data for k in ["smiles", "chembl_id", "target_id"])
        return False

    async def _search_molecule_by_smiles_async(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Search for molecule in ChEMBL by SMILES

        Args:
            smiles: SMILES string

        Returns:
            Molecule data or None
        """
        url = f"{self.BASE_URL}/molecule.json"
        params = {"molecule_structures__canonical_smiles__flexmatch": smiles, "limit": 1}
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        molecules = data.get("molecules", [])
                        if molecules:
                            logger.info(f"ChEMBL: Found molecule {molecules[0].get('molecule_chembl_id')}")
                            return molecules[0]
                    return None
        except Exception as e:
            logger.error(f"ChEMBL: Error searching molecule: {e}")
            return None

    async def _get_molecule_by_id_async(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """
        Get molecule details by ChEMBL ID

        Args:
            chembl_id: ChEMBL molecule ID

        Returns:
            Molecule data or None
        """
        url = f"{self.BASE_URL}/molecule/{chembl_id}.json"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    return None
        except Exception as e:
            logger.error(f"ChEMBL: Error fetching molecule: {e}")
            return None

    async def _get_activities_async(
        self,
        chembl_id: str,
        activity_type: Optional[str] = None,
        target_type: Optional[str] = None
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get bioactivity data for a molecule

        Args:
            chembl_id: ChEMBL molecule ID
            activity_type: Filter by activity type (IC50, Ki, Kd, etc.)
            target_type: Filter by target type (PROTEIN, ORGANISM, etc.)

        Returns:
            List of activities or None
        """
        url = f"{self.BASE_URL}/activity.json"
        params = {
            "molecule_chembl_id": chembl_id,
            "limit": self.config.get("max_activities", 200)
        }

        if activity_type:
            params["standard_type"] = activity_type
        if target_type:
            params["target_type"] = target_type

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        activities = data.get("activities", [])
                        logger.info(f"ChEMBL: Found {len(activities)} activities for {chembl_id}")
                        return activities
                    return []
        except Exception as e:
            logger.error(f"ChEMBL: Error fetching activities: {e}")
            return None

    async def _get_target_async(self, target_id: str) -> Optional[Dict[str, Any]]:
        """
        Get target information

        Args:
            target_id: ChEMBL target ID

        Returns:
            Target data or None
        """
        url = f"{self.BASE_URL}/target/{target_id}.json"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        logger.info(f"ChEMBL: Retrieved target {target_id}")
                        return data
                    return None
        except Exception as e:
            logger.error(f"ChEMBL: Error fetching target: {e}")
            return None

    async def _get_drug_info_async(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """
        Get drug/clinical candidate information

        Args:
            chembl_id: ChEMBL molecule ID

        Returns:
            Drug information or None
        """
        url = f"{self.BASE_URL}/drug.json"
        params = {"molecule_chembl_id": chembl_id}
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        drugs = data.get("drugs", [])
                        if drugs:
                            logger.info(f"ChEMBL: Found drug info for {chembl_id}")
                            return drugs[0]
                    return None
        except Exception as e:
            logger.error(f"ChEMBL: Error fetching drug info: {e}")
            return None

    async def _get_mechanism_async(self, chembl_id: str) -> Optional[List[Dict[str, Any]]]:
        """
        Get mechanism of action data

        Args:
            chembl_id: ChEMBL molecule ID

        Returns:
            List of mechanisms or None
        """
        url = f"{self.BASE_URL}/mechanism.json"
        params = {"molecule_chembl_id": chembl_id}
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        mechanisms = data.get("mechanisms", [])
                        logger.info(f"ChEMBL: Found {len(mechanisms)} mechanisms for {chembl_id}")
                        return mechanisms
                    return []
        except Exception as e:
            logger.error(f"ChEMBL: Error fetching mechanisms: {e}")
            return None

    async def _similarity_search_async(
        self,
        smiles: str,
        similarity: int = 70
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Perform similarity search

        Args:
            smiles: SMILES string
            similarity: Similarity threshold (0-100)

        Returns:
            List of similar molecules or None
        """
        url = f"{self.BASE_URL}/similarity/{smiles}/{similarity}.json"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        molecules = data.get("molecules", [])
                        logger.info(f"ChEMBL: Found {len(molecules)} similar molecules")
                        return molecules
                    return []
        except Exception as e:
            logger.error(f"ChEMBL: Error in similarity search: {e}")
            return None

    async def _substructure_search_async(self, smiles: str) -> Optional[List[Dict[str, Any]]]:
        """
        Perform substructure search

        Args:
            smiles: SMILES substructure query

        Returns:
            List of molecules containing substructure or None
        """
        url = f"{self.BASE_URL}/substructure/{smiles}.json"
        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 90))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        molecules = data.get("molecules", [])
                        logger.info(f"ChEMBL: Found {len(molecules)} molecules with substructure")
                        return molecules
                    return []
        except Exception as e:
            logger.error(f"ChEMBL: Error in substructure search: {e}")
            return None

    def _filter_activities(
        self,
        activities: List[Dict[str, Any]],
        max_value: Optional[float] = None,
        min_value: Optional[float] = None,
        relation: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Filter activities by value and relation

        Args:
            activities: List of activity records
            max_value: Maximum activity value (nM)
            min_value: Minimum activity value (nM)
            relation: Relation type ('=', '<', '>', '<=', '>=')

        Returns:
            Filtered activities
        """
        filtered = []

        for activity in activities:
            value = activity.get("standard_value")
            units = activity.get("standard_units", "").upper()
            rel = activity.get("standard_relation", "=")

            # Only process nM values
            if not value or units != "NM":
                continue

            try:
                val_float = float(value)

                # Apply filters
                if max_value is not None and val_float > max_value:
                    continue
                if min_value is not None and val_float < min_value:
                    continue
                if relation and rel != relation:
                    continue

                filtered.append(activity)
            except (ValueError, TypeError):
                continue

        return filtered

    def _summarize_activities(self, activities: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Summarize bioactivity data

        Args:
            activities: List of activities

        Returns:
            Summary statistics
        """
        if not activities:
            return {
                "num_activities": 0,
                "num_targets": 0,
                "assay_types": [],
                "best_ic50_nm": None,
                "best_ki_nm": None,
                "best_kd_nm": None
            }

        targets = set()
        assay_types = set()
        ic50_values = []
        ki_values = []
        kd_values = []

        for activity in activities:
            if activity.get("target_chembl_id"):
                targets.add(activity["target_chembl_id"])

            if activity.get("assay_type"):
                assay_types.add(activity["assay_type"])

            standard_type = (activity.get("standard_type") or "").upper()
            standard_value = activity.get("standard_value")
            standard_units = (activity.get("standard_units") or "").upper()

            if standard_value and standard_units == "NM":
                try:
                    value_nm = float(standard_value)
                    if standard_type == "IC50":
                        ic50_values.append(value_nm)
                    elif standard_type == "KI":
                        ki_values.append(value_nm)
                    elif standard_type == "KD":
                        kd_values.append(value_nm)
                except (ValueError, TypeError):
                    pass

        return {
            "num_activities": len(activities),
            "num_targets": len(targets),
            "assay_types": list(assay_types),
            "best_ic50_nm": min(ic50_values) if ic50_values else None,
            "best_ki_nm": min(ki_values) if ki_values else None,
            "best_kd_nm": min(kd_values) if kd_values else None,
            "targets": list(targets)[:20]
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute enhanced ChEMBL query

        Args:
            input_data: SMILES, ChEMBL ID, or dict
            **kwargs: Additional parameters:
                - mode: "bioactivity", "target", "drug", "mechanism", "similarity", "substructure"
                - activity_type: str (IC50, Ki, Kd, etc.)
                - target_type: str (PROTEIN, etc.)
                - max_value: float (filter activities)
                - min_value: float
                - similarity: int (0-100)

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
            # Check if it's a ChEMBL ID or SMILES
            if input_data.startswith("CHEMBL"):
                chembl_id = input_data
                smiles = None
                target_id = None
            elif input_data.startswith("TARGET:"):
                target_id = input_data.replace("TARGET:", "")
                chembl_id = None
                smiles = None
            else:
                smiles = input_data
                chembl_id = None
                target_id = None
        else:
            smiles = input_data.get("smiles")
            chembl_id = input_data.get("chembl_id")
            target_id = input_data.get("target_id")

        # Get mode
        mode = kwargs.get("mode", "bioactivity")

        # Rate limiting
        await asyncio.sleep(self.config.get("rate_limit_delay", 0.5))

        # Execute based on mode
        try:
            if mode == "bioactivity":
                # Get ChEMBL ID if we have SMILES
                if not chembl_id and smiles:
                    molecule = await self._search_molecule_by_smiles_async(smiles)
                    if not molecule:
                        return AdapterResult(success=False, data=None, error="Molecule not found in ChEMBL")
                    chembl_id = molecule.get("molecule_chembl_id")

                # Get activities
                activities = await self._get_activities_async(
                    chembl_id,
                    activity_type=kwargs.get("activity_type"),
                    target_type=kwargs.get("target_type")
                )

                if activities is None:
                    return AdapterResult(success=False, data=None, error="Failed to fetch activities")

                # Filter activities
                activities = self._filter_activities(
                    activities,
                    max_value=kwargs.get("max_value"),
                    min_value=kwargs.get("min_value"),
                    relation=kwargs.get("relation")
                )

                summary = self._summarize_activities(activities)

                result_data = {
                    "chembl_id": chembl_id,
                    **summary,
                    "bioactivities": activities
                }

            elif mode == "target" and target_id:
                target_data = await self._get_target_async(target_id)
                if not target_data:
                    return AdapterResult(success=False, data=None, error="Target not found")

                result_data = {
                    "target_id": target_id,
                    "target_type": target_data.get("target_type"),
                    "organism": target_data.get("organism"),
                    "pref_name": target_data.get("pref_name"),
                    "target_components": target_data.get("target_components", [])
                }

            elif mode == "drug":
                if not chembl_id and smiles:
                    molecule = await self._search_molecule_by_smiles_async(smiles)
                    if molecule:
                        chembl_id = molecule.get("molecule_chembl_id")

                if not chembl_id:
                    return AdapterResult(success=False, data=None, error="ChEMBL ID required")

                drug_info = await self._get_drug_info_async(chembl_id)

                result_data = {
                    "chembl_id": chembl_id,
                    "is_drug": drug_info is not None,
                    "drug_info": drug_info
                }

            elif mode == "mechanism":
                if not chembl_id and smiles:
                    molecule = await self._search_molecule_by_smiles_async(smiles)
                    if molecule:
                        chembl_id = molecule.get("molecule_chembl_id")

                if not chembl_id:
                    return AdapterResult(success=False, data=None, error="ChEMBL ID required")

                mechanisms = await self._get_mechanism_async(chembl_id)

                result_data = {
                    "chembl_id": chembl_id,
                    "num_mechanisms": len(mechanisms) if mechanisms else 0,
                    "mechanisms": mechanisms
                }

            elif mode == "similarity" and smiles:
                similarity = kwargs.get("similarity", self.config.get("similarity_threshold", 70))
                molecules = await self._similarity_search_async(smiles, similarity)

                if molecules is None:
                    return AdapterResult(success=False, data=None, error="Similarity search failed")

                result_data = {
                    "query_smiles": smiles,
                    "similarity_threshold": similarity,
                    "num_results": len(molecules),
                    "molecules": molecules
                }

            elif mode == "substructure" and smiles:
                molecules = await self._substructure_search_async(smiles)

                if molecules is None:
                    return AdapterResult(success=False, data=None, error="Substructure search failed")

                result_data = {
                    "query_smiles": smiles,
                    "num_results": len(molecules),
                    "molecules": molecules
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
                    "source": "chembl",
                    "adapter_version": self.version,
                    "mode": mode
                }
            )

        except Exception as e:
            logger.error(f"ChEMBL: Unexpected error: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unexpected error: {str(e)}"
            )
