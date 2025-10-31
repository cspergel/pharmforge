"""
ADMET-AI Adapter for PharmForge

Provides 49 ADMET property predictions using ADMET-AI Chemprop models.
This is a clean, simple implementation focusing on ADMET-AI only.

Future enhancements (custom models, ensembles) will be added after validation.
"""

import hashlib
import logging
from typing import Dict, Any, List, Optional

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ADMETaiAdapter(AdapterProtocol):
    """
    ADMET-AI adapter providing ML-based ADMET predictions.

    Uses pre-trained Chemprop-RDKit models from the ADMET-AI library.
    Covers 49 ADMET properties across absorption, distribution, metabolism,
    excretion, and toxicity endpoints.

    Reference: https://github.com/swansonk14/admet_ai
    """

    def __init__(self, name: str = "admet_ai", adapter_type: str = "ml", config: Optional[Dict[str, Any]] = None):
        """
        Initialize ADMET-AI adapter.

        Args:
            name: Adapter name (default: "admet_ai")
            adapter_type: Adapter type (default: "ml")
            config: Optional configuration dictionary. Supported keys:
                   - properties: List of specific properties to predict (None = all 49 properties)
        """
        super().__init__(name, adapter_type, config)
        self.version = "1.0.0"

        # Properties to predict (None = all 49 properties)
        self.properties = self.config.get('properties', None)

        # Lazy-load ADMET-AI model (heavy import)
        self._model = None

    @property
    def model(self):
        """Lazy-load ADMET-AI model on first use."""
        if self._model is None:
            try:
                import argparse
                import torch
                import numpy as np
                import numpy.core.multiarray
                from admet_ai import ADMETModel

                # Fix for PyTorch 2.6+ weights_only=True default
                # ADMET AI models contain various Python objects that need to be allowlisted
                # We trust these models from the official ADMET AI library
                logger.info("Configuring PyTorch safe globals for ADMET-AI...")

                # Add all required safe globals for ADMET AI model loading
                safe_globals = [
                    argparse.Namespace,
                    numpy.core.multiarray._reconstruct,
                    numpy.ndarray,
                    numpy.dtype,
                    numpy.dtypes.Float64DType,
                ]

                torch.serialization.add_safe_globals(safe_globals)
                logger.info(f"Added {len(safe_globals)} safe globals for model loading")

                logger.info("Loading ADMET-AI model...")
                self._model = ADMETModel()
                logger.info("✓ ADMET-AI model loaded successfully")
            except ImportError as e:
                logger.error(f"Failed to import admet_ai: {e}")
                raise ImportError(
                    "ADMET-AI library not installed. "
                    "Install with: pip install admet-ai"
                ) from e
            except Exception as e:
                logger.error(f"Failed to load ADMET-AI model: {e}")
                raise

        return self._model

    def validate_input(self, smiles: str) -> bool:
        """
        Validate SMILES string input.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not smiles or not isinstance(smiles, str):
            return False

        # Basic SMILES validation with RDKit
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.warning(f"SMILES validation error: {e}")
            return False

    async def execute(self, smiles: str, **params) -> AdapterResult:
        """
        Execute ADMET predictions for a molecule.

        Args:
            smiles: SMILES string of the molecule
            **params: Additional parameters (e.g., specific properties to predict)

        Returns:
            AdapterResult containing predictions for all ADMET properties
        """
        try:
            # Validate input
            if not self.validate_input(smiles):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid SMILES string"
                )

            # Get properties to predict (from params or instance default)
            properties = params.get('properties', self.properties)

            # Generate cache key
            cache_key = self.generate_cache_key(smiles, properties=properties)

            # Make predictions using ADMET-AI
            logger.info(f"Predicting ADMET properties for: {smiles[:50]}...")
            predictions = self.model.predict(smiles=smiles)

            # Filter to requested properties if specified
            if properties:
                predictions = {
                    k: v for k, v in predictions.items()
                    if k in properties
                }

            # Separate raw predictions from percentiles
            raw_predictions = {
                k: v for k, v in predictions.items()
                if not k.endswith('_drugbank_approved_percentile')
            }

            percentiles = {
                k: v for k, v in predictions.items()
                if k.endswith('_drugbank_approved_percentile')
            }

            # Create result data structure
            result_data = {
                "smiles": smiles,
                "properties": raw_predictions,
                "percentiles": percentiles,
                "property_count": len(raw_predictions),
                "model": "ADMET-AI v1.4.0",
                "reference": "Trained on TDC ADMET datasets"
            }

            logger.info(f"✓ Predicted {len(raw_predictions)} ADMET properties")

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "cache_key": cache_key,
                    "version": self.version
                }
            )

        except Exception as e:
            logger.error(f"ADMET-AI prediction failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name
                }
            )

    def generate_cache_key(
        self,
        smiles: str,
        properties: Optional[List[str]] = None
    ) -> str:
        """
        Generate deterministic cache key for predictions.

        Args:
            smiles: SMILES string
            properties: Optional list of properties being predicted

        Returns:
            SHA256 hash as cache key
        """
        # Canonicalize SMILES for consistent caching
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            pass  # Use original SMILES if canonicalization fails

        # Include properties in cache key if specified
        props_str = ','.join(sorted(properties)) if properties else 'all'

        # Create cache key
        key_string = f"admet_ai_v1.0.0:{smiles}:{props_str}"
        return hashlib.sha256(key_string.encode()).hexdigest()

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
            "description": "ADMET property predictions using ADMET-AI",
            "properties": {
                "total_properties": 49,
                "categories": [
                    "Absorption",
                    "Distribution",
                    "Metabolism",
                    "Excretion",
                    "Toxicity",
                    "Physicochemical"
                ],
                "model": "Chemprop-RDKit",
                "reference": "https://github.com/swansonk14/admet_ai"
            },
            "config": {
                "properties_filter": self.properties or "all"
            }
        }

    # Helper methods for accessing specific property groups

    def get_absorption_properties(self) -> List[str]:
        """Get list of absorption-related properties."""
        return [
            "Caco2_Wang",
            "HIA_Hou",
            "Pgp_Broccatelli",
            "Bioavailability_Ma",
            "PAMPA_NCATS",
            "Solubility_AqSolDB",
            "HydrationFreeEnergy_FreeSolv"
        ]

    def get_distribution_properties(self) -> List[str]:
        """Get list of distribution-related properties."""
        return [
            "BBB_Martins",
            "PPBR_AZ",
            "VDss_Lombardo",
            "Lipophilicity_AstraZeneca"
        ]

    def get_metabolism_properties(self) -> List[str]:
        """Get list of metabolism-related properties."""
        return [
            "CYP1A2_Veith",
            "CYP2C19_Veith",
            "CYP2C9_Veith",
            "CYP2D6_Veith",
            "CYP3A4_Veith",
            "CYP2C9_Substrate_CarbonMangels",
            "CYP2D6_Substrate_CarbonMangels",
            "CYP3A4_Substrate_CarbonMangels"
        ]

    def get_excretion_properties(self) -> List[str]:
        """Get list of excretion-related properties."""
        return [
            "Clearance_Hepatocyte_AZ",
            "Clearance_Microsome_AZ",
            "Half_Life_Obach"
        ]

    def get_toxicity_properties(self) -> List[str]:
        """Get list of toxicity-related properties."""
        return [
            "AMES",
            "BBB_Martins",
            "hERG",
            "DILI",
            "LD50_Zhu",
            "Carcinogens_Lagunin",
            "ClinTox",
            "Skin_Reaction",
            "NR-AR",
            "NR-AR-LBD",
            "NR-AhR",
            "NR-Aromatase",
            "NR-ER",
            "NR-ER-LBD",
            "NR-PPAR-gamma",
            "SR-ARE",
            "SR-ATAD5",
            "SR-HSE",
            "SR-MMP",
            "SR-p53"
        ]

    def get_physicochemical_properties(self) -> List[str]:
        """Get list of physicochemical descriptor properties."""
        return [
            "molecular_weight",
            "logP",
            "tpsa",
            "hydrogen_bond_acceptors",
            "hydrogen_bond_donors",
            "stereo_centers",
            "QED",
            "Lipinski"
        ]
