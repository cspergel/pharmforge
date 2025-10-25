"""
RDKit Local Adapter - Local molecular property calculations
Provides fallback when PubChem API is unavailable
"""
from typing import Any, Dict, Optional
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, Crippen
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - install with: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class RDKitAdapter(AdapterProtocol):
    """
    Adapter for local RDKit calculations
    Computes molecular properties locally (no API calls required)
    """

    def __init__(self):
        super().__init__(
            name="rdkit_local",
            adapter_type="local",
            config={
                "timeout": 10  # Local calculations should be fast
            }
        )
        self.version = "1.0.0"

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string

        Args:
            input_data: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE:
            return False

        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False

        # Try to parse with RDKit
        try:
            mol = Chem.MolFromSmiles(input_data)
            return mol is not None
        except Exception:
            return False

    def _calculate_properties(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Calculate molecular properties using RDKit

        Args:
            smiles: SMILES string

        Returns:
            Dictionary of calculated properties
        """
        if not RDKIT_AVAILABLE:
            raise Exception("RDKit is not installed")

        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                logger.warning(f"RDKit: Could not parse SMILES: {smiles}")
                return None

            # Calculate properties
            properties = {
                # Basic properties
                "molecular_weight": Descriptors.MolWt(mol),
                "logp": Crippen.MolLogP(mol),
                "tpsa": Descriptors.TPSA(mol),

                # Lipinski Rule of 5
                "h_bond_donors": Descriptors.NumHDonors(mol),
                "h_bond_acceptors": Descriptors.NumHAcceptors(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),

                # Counts
                "num_atoms": mol.GetNumAtoms(),
                "num_heavy_atoms": Descriptors.HeavyAtomCount(mol),
                "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
                "num_rings": Descriptors.RingCount(mol),

                # Canonical SMILES from RDKit
                "canonical_smiles": Chem.MolToSmiles(mol),

                # Lipinski violations
                "lipinski_violations": self._count_lipinski_violations(mol)
            }

            return properties

        except Exception as e:
            logger.error(f"RDKit: Error calculating properties for {smiles}: {e}")
            return None

    def _count_lipinski_violations(self, mol) -> int:
        """
        Count number of Lipinski Rule of 5 violations

        Args:
            mol: RDKit molecule object

        Returns:
            Number of violations (0-4)
        """
        violations = 0

        # Rule 1: Molecular weight <= 500
        if Descriptors.MolWt(mol) > 500:
            violations += 1

        # Rule 2: LogP <= 5
        if Crippen.MolLogP(mol) > 5:
            violations += 1

        # Rule 3: H-bond donors <= 5
        if Descriptors.NumHDonors(mol) > 5:
            violations += 1

        # Rule 4: H-bond acceptors <= 10
        if Descriptors.NumHAcceptors(mol) > 10:
            violations += 1

        return violations

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute local RDKit property calculations

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters

        Returns:
            AdapterResult containing calculated properties
        """
        # Check if RDKit is available
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string or could not parse with RDKit"
            )

        smiles = input_data

        # Calculate properties (synchronous, but fast)
        props = self._calculate_properties(smiles)

        if props is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to calculate properties with RDKit",
                metadata={"source": "rdkit_local", "smiles": smiles}
            )

        return AdapterResult(
            success=True,
            data=props,
            cache_hit=False,
            metadata={
                "source": "rdkit_local",
                "smiles": smiles,
                "adapter_version": self.version,
                "computation_type": "local"
            }
        )
