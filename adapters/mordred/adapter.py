"""
Mordred Adapter - Comprehensive molecular descriptor calculator
Provides 1800+ molecular descriptors using the Mordred library
"""
from typing import Any, Dict, Optional
import logging

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - required for Mordred. Install with: pip install rdkit")

try:
    from mordred import Calculator, descriptors
    MORDRED_AVAILABLE = True
except ImportError:
    MORDRED_AVAILABLE = False
    logging.warning("Mordred not available - install with: pip install mordred")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class MordredAdapter(AdapterProtocol):
    """
    Adapter for Mordred molecular descriptor calculations
    Computes 1800+ molecular descriptors locally (no API calls required)
    """

    def __init__(self):
        super().__init__(
            name="mordred",
            adapter_type="local",
            config={
                "timeout": 30,  # Descriptor calculation can take time for large molecules
                "ignore_3d": True  # Set to False if you want 3D descriptors (requires conformer generation)
            }
        )
        self.version = "1.0.0"
        self.calculator = None

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Mordred requires RDKit. Install with: pip install rdkit")
        elif not MORDRED_AVAILABLE:
            logger.error("Mordred is not installed! Install with: pip install mordred")
        else:
            # Initialize the Mordred calculator with all descriptors
            try:
                self.calculator = Calculator(descriptors, ignore_3D=self.config.get("ignore_3d", True))
                logger.info(f"Mordred calculator initialized with {len(self.calculator.descriptors)} descriptors")
            except Exception as e:
                logger.error(f"Failed to initialize Mordred calculator: {e}")
                self.calculator = None

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string

        Args:
            input_data: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not RDKIT_AVAILABLE or not MORDRED_AVAILABLE:
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

    def _calculate_descriptors(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Calculate molecular descriptors using Mordred

        Args:
            smiles: SMILES string

        Returns:
            Dictionary of calculated descriptors
        """
        if not RDKIT_AVAILABLE or not MORDRED_AVAILABLE:
            raise Exception("RDKit and Mordred are required but not installed")

        if self.calculator is None:
            raise Exception("Mordred calculator is not initialized")

        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                logger.warning(f"Mordred: Could not parse SMILES: {smiles}")
                return None

            # Calculate all descriptors
            logger.info(f"Mordred: Calculating {len(self.calculator.descriptors)} descriptors for {smiles}")
            result = self.calculator(mol)

            # Convert result to dictionary, handling missing/error values
            descriptors_dict = {}
            error_count = 0
            missing_count = 0

            for desc_name, desc_value in zip(self.calculator.descriptors, result):
                # Get the descriptor name as string
                name = str(desc_name)

                # Handle different types of descriptor values
                if desc_value is None:
                    missing_count += 1
                    descriptors_dict[name] = None
                elif isinstance(desc_value, Exception):
                    error_count += 1
                    descriptors_dict[name] = None
                else:
                    # Convert to native Python types for JSON serialization
                    try:
                        # Handle numpy types and other numeric types
                        if hasattr(desc_value, 'item'):
                            descriptors_dict[name] = desc_value.item()
                        else:
                            descriptors_dict[name] = float(desc_value) if isinstance(desc_value, (int, float)) else desc_value
                    except (ValueError, TypeError):
                        descriptors_dict[name] = None
                        error_count += 1

            logger.info(f"Mordred: Calculated {len(descriptors_dict) - error_count - missing_count} valid descriptors "
                       f"({error_count} errors, {missing_count} missing)")

            # Add metadata about the calculation
            metadata = {
                "total_descriptors": len(descriptors_dict),
                "valid_descriptors": len(descriptors_dict) - error_count - missing_count,
                "error_count": error_count,
                "missing_count": missing_count,
                "canonical_smiles": Chem.MolToSmiles(mol)
            }

            return {
                "descriptors": descriptors_dict,
                "metadata": metadata
            }

        except Exception as e:
            logger.error(f"Mordred: Error calculating descriptors for {smiles}: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Mordred descriptor calculations

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters
                - ignore_3d: Override config to ignore 3D descriptors (default: True)

        Returns:
            AdapterResult containing calculated descriptors
        """
        # Check if dependencies are available
        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Mordred requires RDKit. Install with: pip install rdkit"
            )

        if not MORDRED_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="Mordred is not installed. Install with: pip install mordred"
            )

        if self.calculator is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Mordred calculator failed to initialize"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string or could not parse with RDKit"
            )

        smiles = input_data

        # Calculate descriptors (synchronous, but may be slow for large molecules)
        try:
            result = self._calculate_descriptors(smiles)

            if result is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to calculate descriptors with Mordred",
                    metadata={"source": "mordred", "smiles": smiles}
                )

            return AdapterResult(
                success=True,
                data=result,
                cache_hit=False,
                metadata={
                    "source": "mordred",
                    "smiles": smiles,
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "ignore_3d": self.config.get("ignore_3d", True),
                    **result.get("metadata", {})
                }
            )

        except Exception as e:
            logger.error(f"Mordred: Unexpected error during execution: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unexpected error: {str(e)}",
                metadata={"source": "mordred", "smiles": smiles}
            )
