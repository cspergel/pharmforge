"""
PubChem Adapter - Fetches molecular properties from PubChem REST API
"""
from typing import Any, Dict, Optional
import aiohttp
import asyncio
import logging

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class PubChemAdapter(AdapterProtocol):
    """
    Adapter for PubChem REST API
    Fetches molecular properties for compounds by SMILES
    """

    BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def __init__(self):
        super().__init__(
            name="pubchem",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.2,  # 5 requests/second
                "timeout": 30
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string

        Args:
            input_data: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        # Basic SMILES validation - contains valid characters
        valid_chars = set("CNOPSFClBrI[]()=@#+-0123456789cnops")
        return all(c in valid_chars for c in input_data.replace(" ", ""))

    async def _fetch_properties_async(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Fetch molecular properties from PubChem API

        Args:
            smiles: SMILES string of the molecule

        Returns:
            Dictionary of molecular properties or None if not found
        """
        # URL encode the SMILES string
        from urllib.parse import quote
        encoded_smiles = quote(smiles)

        # Request specific properties
        properties = [
            "MolecularWeight",
            "XLogP",
            "TPSA",
            "HBondDonorCount",
            "HBondAcceptorCount",
            "RotatableBondCount",
            "CanonicalSMILES",
            "InChI",
            "InChIKey"
        ]

        url = f"{self.BASE_URL}/compound/smiles/{encoded_smiles}/property/{','.join(properties)}/JSON"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

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
                        logger.error(f"PubChem API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"PubChem: Timeout fetching properties for {smiles}")
            return None
        except Exception as e:
            logger.error(f"PubChem: Error fetching properties for {smiles}: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the PubChem property lookup

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing molecular properties
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
        rate_delay = self.config.get("rate_limit_delay", 0.2)
        await asyncio.sleep(rate_delay)

        # Fetch properties
        props = await self._fetch_properties_async(smiles)

        if props is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to fetch properties from PubChem",
                metadata={"source": "pubchem", "smiles": smiles}
            )

        # Format the result - convert numeric fields to float
        def to_float(val):
            """Safely convert value to float, return None if not possible"""
            if val is None:
                return None
            try:
                return float(val)
            except (ValueError, TypeError):
                return None

        def to_int(val):
            """Safely convert value to int, return None if not possible"""
            if val is None:
                return None
            try:
                return int(val)
            except (ValueError, TypeError):
                return None

        result_data = {
            "molecular_weight": to_float(props.get("MolecularWeight")),
            "logp": to_float(props.get("XLogP")),
            "tpsa": to_float(props.get("TPSA")),
            "h_bond_donors": to_int(props.get("HBondDonorCount")),
            "h_bond_acceptors": to_int(props.get("HBondAcceptorCount")),
            "rotatable_bonds": to_int(props.get("RotatableBondCount")),
            "canonical_smiles": props.get("CanonicalSMILES"),
            "inchi": props.get("InChI"),
            "inchikey": props.get("InChIKey")
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "pubchem",
                "smiles": smiles,
                "adapter_version": self.version
            }
        )
