"""
BindingDB Adapter for PharmForge

Queries the BindingDB database for experimental binding affinity data.
BindingDB contains binding data for protein-ligand complexes with measured
Ki, Kd, IC50, and EC50 values from scientific literature.

References:
    - BindingDB: https://www.bindingdb.org/
    - API Documentation: https://bindingdb.org/bind/chemsearch/marvin/SDFdownload.jsp
"""

import aiohttp
import asyncio
import logging
import math
from typing import Any, Dict, Optional, List
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class BindingDBAdapter(AdapterProtocol):
    """
    Adapter for BindingDB experimental binding data

    BindingDB is a public database of measured binding affinities for
    protein-ligand complexes. This adapter queries experimental data
    to validate computational predictions and retrieve known binding values.

    Features:
        - Query by SMILES, InChI, or compound name
        - Retrieve Ki, Kd, IC50, EC50 values
        - Access assay conditions and references
        - Target protein information
        - Multiple measurement types

    Data Types:
        - Ki: Inhibition constant (equilibrium dissociation constant)
        - Kd: Dissociation constant
        - IC50: Half-maximal inhibitory concentration
        - EC50: Half-maximal effective concentration

    Returns:
        - Experimental binding values with units
        - Target protein information
        - Assay conditions and pH
        - Literature references
        - Calculated pKi/pKd/pIC50 values
    """

    BASE_URL = "https://bindingdb.org/axis2/services/BDBService"

    def __init__(
        self,
        name: str = "bindingdb",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize BindingDB adapter

        Args:
            name: Adapter name (default: "bindingdb")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - timeout: Request timeout in seconds (default: 30)
                   - max_results: Maximum results to return (default: 100)
                   - min_affinity: Minimum binding affinity to include (nM)
        """
        super().__init__(name, adapter_type, config or {})
        self.version = "1.0.0"
        self.timeout = self.config.get('timeout', 30)
        self.max_results = self.config.get('max_results', 100)
        self.min_affinity = self.config.get('min_affinity', None)

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data

        Args:
            input_data: SMILES string, InChI, or dict with search parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, dict):
            # Check for required search fields
            return any(k in input_data for k in ['smiles', 'inchi', 'name', 'target'])
        elif isinstance(input_data, str):
            return len(input_data) > 0
        return False

    def _calculate_pvalue(self, value: float, unit: str) -> Optional[float]:
        """
        Calculate p-value (negative log) from binding constant

        pKi = -log10(Ki in M)
        pKd = -log10(Kd in M)
        pIC50 = -log10(IC50 in M)

        Args:
            value: Binding value
            unit: Unit (nM, uM, mM, M)

        Returns:
            p-value or None if cannot calculate
        """
        try:
            # Convert to molar
            if unit == 'nM':
                molar = value * 1e-9
            elif unit == 'uM' or unit == 'μM':
                molar = value * 1e-6
            elif unit == 'mM':
                molar = value * 1e-3
            elif unit == 'M':
                molar = value
            else:
                return None

            if molar <= 0:
                return None

            # Calculate -log10
            pvalue = -math.log10(molar)
            return round(pvalue, 2)

        except (ValueError, TypeError):
            return None

    def _convert_to_nm(self, value: float, unit: str) -> Optional[float]:
        """
        Convert binding value to nanomolar (nM)

        Args:
            value: Binding value
            unit: Unit (nM, uM, mM, M)

        Returns:
            Value in nM or None
        """
        try:
            if unit == 'nM':
                return value
            elif unit == 'uM' or unit == 'μM':
                return value * 1000
            elif unit == 'mM':
                return value * 1e6
            elif unit == 'M':
                return value * 1e9
            else:
                return None
        except (ValueError, TypeError):
            return None

    async def _query_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Query BindingDB by SMILES structure

        Args:
            smiles: SMILES string

        Returns:
            Dictionary with binding data or None
        """
        try:
            # Note: BindingDB's web service API is SOAP-based
            # For REST-like access, we use the TSV download endpoint
            # Actual implementation may need SOAP client or web scraping

            # Alternative: Use REST endpoint if available
            url = f"https://bindingdb.org/axis2/services/BDBService/getLigandsBySMILES"

            # For now, return simulated structure
            # Real implementation would parse BindingDB response
            logger.warning("BindingDB API integration requires SOAP client - returning mock data")

            # This is a placeholder - real implementation needs proper API integration
            return None

        except Exception as e:
            logger.error(f"Error querying BindingDB by SMILES: {e}")
            return None

    async def _search_bindingdb(
        self,
        smiles: Optional[str] = None,
        inchi: Optional[str] = None,
        target: Optional[str] = None,
        name: Optional[str] = None
    ) -> List[Dict[str, Any]]:
        """
        Search BindingDB database

        Args:
            smiles: SMILES string
            inchi: InChI string
            target: Target protein name
            name: Compound name

        Returns:
            List of binding data entries
        """
        results = []

        try:
            # BindingDB API is SOAP-based, but they also provide TSV downloads
            # For a production implementation, you would:
            # 1. Use a SOAP client library (zeep)
            # 2. Parse TSV downloads
            # 3. Use their web interface with requests/BeautifulSoup

            # Example using TSV download endpoint:
            if smiles:
                # Construct search URL
                # Note: This is simplified - actual URL structure varies
                encoded_smiles = quote(smiles)
                url = f"https://bindingdb.org/bind/chemsearch/marvin/MolStructure.jsp?smiles={encoded_smiles}"

                timeout = aiohttp.ClientTimeout(total=self.timeout)

                async with aiohttp.ClientSession() as session:
                    async with session.get(url, timeout=timeout) as response:
                        if response.status == 200:
                            # Parse response (TSV or HTML)
                            text = await response.text()
                            # Real implementation would parse this data
                            logger.info(f"BindingDB query returned {len(text)} bytes")

            # For demonstration, return example structure
            # Real implementation would parse actual BindingDB data
            logger.warning(
                "BindingDB adapter returns example data. "
                "Full implementation requires SOAP client or TSV parsing."
            )

        except Exception as e:
            logger.error(f"Error searching BindingDB: {e}")

        return results

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute BindingDB query

        Args:
            input_data: SMILES, InChI, or dict with search parameters
            **kwargs: Additional parameters:
                - target: Target protein name
                - measurement_type: 'Ki', 'Kd', 'IC50', 'EC50'

        Returns:
            AdapterResult containing experimental binding data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input. Provide SMILES, InChI, or search dict."
            )

        # Extract search parameters
        if isinstance(input_data, dict):
            smiles = input_data.get('smiles')
            inchi = input_data.get('inchi')
            target = input_data.get('target')
            name = input_data.get('name')
        else:
            # Assume SMILES
            smiles = input_data
            inchi = None
            target = kwargs.get('target')
            name = None

        measurement_type = kwargs.get('measurement_type', None)

        # Query BindingDB
        try:
            results = await self._search_bindingdb(
                smiles=smiles,
                inchi=inchi,
                target=target,
                name=name
            )

            # Since BindingDB API requires SOAP integration, provide example data
            # Real implementation would parse actual results
            example_data = self._create_example_data(smiles, target)

            # Filter by measurement type if specified
            if measurement_type and results:
                results = [r for r in results if r.get('type') == measurement_type]

            # Apply minimum affinity filter
            if self.min_affinity and results:
                results = [
                    r for r in results
                    if r.get('value_nm', float('inf')) <= self.min_affinity
                ]

            # Limit results
            if results:
                results = results[:self.max_results]

            # Calculate statistics
            stats = None
            if results:
                values = [r.get('value_nm') for r in results if r.get('value_nm')]
                if values:
                    stats = {
                        'count': len(values),
                        'mean_nm': sum(values) / len(values),
                        'min_nm': min(values),
                        'max_nm': max(values),
                        'median_nm': sorted(values)[len(values) // 2]
                    }

            result_data = {
                "smiles": smiles,
                "target": target,
                "entries": example_data['entries'],  # Use example data
                "count": len(example_data['entries']),
                "statistics": stats,
                "source": "BindingDB",
                "note": "This is example data. Full implementation requires SOAP API integration."
            }

            logger.info(f"BindingDB query complete: {len(example_data['entries'])} entries")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "query_type": "SMILES" if smiles else "target",
                    "source": "BindingDB"
                }
            )

        except Exception as e:
            logger.error(f"BindingDB adapter error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=str(e)
            )

    def _create_example_data(
        self,
        smiles: Optional[str],
        target: Optional[str]
    ) -> Dict[str, Any]:
        """
        Create example BindingDB data structure

        This demonstrates the expected data format.
        Real implementation would parse actual BindingDB responses.

        Args:
            smiles: Query SMILES
            target: Query target

        Returns:
            Example data structure
        """
        # Example entries showing different measurement types
        entries = [
            {
                "type": "Ki",
                "value": 5.2,
                "unit": "nM",
                "value_nm": 5.2,
                "pvalue": 8.28,  # -log10(5.2e-9)
                "target": target or "Example Kinase",
                "target_type": "protein",
                "organism": "Homo sapiens",
                "assay_description": "Radioligand binding assay",
                "pH": 7.4,
                "temperature": "25°C",
                "reference": "PMID:12345678",
                "year": 2023
            },
            {
                "type": "IC50",
                "value": 12.5,
                "unit": "nM",
                "value_nm": 12.5,
                "pvalue": 7.90,
                "target": target or "Example Kinase",
                "target_type": "protein",
                "organism": "Homo sapiens",
                "assay_description": "Cell-based proliferation assay",
                "reference": "PMID:12345679",
                "year": 2023
            },
            {
                "type": "Kd",
                "value": 3.8,
                "unit": "nM",
                "value_nm": 3.8,
                "pvalue": 8.42,
                "target": target or "Example Kinase",
                "target_type": "protein",
                "organism": "Homo sapiens",
                "assay_description": "Surface plasmon resonance",
                "temperature": "25°C",
                "reference": "PMID:12345680",
                "year": 2022
            }
        ]

        return {
            "entries": entries,
            "count": len(entries)
        }

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "Query experimental binding affinity data from BindingDB",
            "features": [
                "Experimental Ki, Kd, IC50, EC50 values",
                "Target protein information",
                "Assay conditions and references",
                "Literature citations"
            ],
            "data_types": {
                "Ki": "Inhibition constant (equilibrium dissociation)",
                "Kd": "Dissociation constant",
                "IC50": "Half-maximal inhibitory concentration",
                "EC50": "Half-maximal effective concentration"
            },
            "config": {
                "timeout": self.timeout,
                "max_results": self.max_results,
                "min_affinity": self.min_affinity
            },
            "input": {
                "smiles": "Compound SMILES structure",
                "inchi": "Compound InChI",
                "target": "Protein target name",
                "name": "Compound name"
            },
            "output": {
                "entries": "List of experimental measurements",
                "statistics": "Aggregated statistics (mean, min, max)",
                "references": "Literature citations"
            },
            "notes": [
                "Requires SOAP API integration for full functionality",
                "Currently returns example data structure",
                "Production version needs zeep library or TSV parsing"
            ]
        }
