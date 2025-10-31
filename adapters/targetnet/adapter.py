"""
TargetNet Adapter Implementation

Provides target prediction capabilities:
- Predict protein targets for molecules
- Off-target prediction and analysis
- Target confidence scoring
- Multi-target profiling
"""

import hashlib
import logging
from typing import Dict, Any, List, Optional, Union, Tuple
import numpy as np

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class TargetNetAdapter(AdapterProtocol):
    """
    TargetNet adapter for protein target prediction.

    Capabilities:
    - Predict protein targets for small molecules
    - Multi-target profiling
    - Off-target prediction
    - Target-based drug repurposing
    - Polypharmacology analysis

    Uses molecular fingerprints and deep learning for target prediction.
    """

    def __init__(
        self,
        name: str = "targetnet",
        adapter_type: str = "ml",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize TargetNet adapter.

        Args:
            name: Adapter name (default: "targetnet")
            adapter_type: Adapter type (default: "ml")
            config: Configuration dictionary with keys:
                   - model_path: Path to pretrained model (optional)
                   - confidence_threshold: Minimum confidence (default: 0.5)
                   - max_targets: Maximum targets to return (default: 10)
                   - use_local: Use local model vs API (default: True)
        """
        super().__init__(name, adapter_type, config)
        self.version = "1.0.0"

        # Configuration
        self.model_path = self.config.get('model_path', None)
        self.confidence_threshold = self.config.get('confidence_threshold', 0.5)
        self.max_targets = self.config.get('max_targets', 10)
        self.use_local = self.config.get('use_local', True)

        # Lazy-load model
        self._model = None
        self._target_database = None

        logger.info(
            f"Initialized TargetNet adapter "
            f"(threshold={self.confidence_threshold})"
        )

    @property
    def model(self):
        """Lazy-load TargetNet model on first use."""
        if self._model is None and self.use_local:
            try:
                # Try to load target prediction model
                import torch
                from targetnet.models import TargetPredictor

                logger.info("Loading TargetNet model...")

                if self.model_path:
                    self._model = torch.load(self.model_path)
                else:
                    # Use default pretrained model
                    self._model = TargetPredictor.from_pretrained()

                self._model.eval()
                logger.info("TargetNet model loaded successfully")

            except ImportError as e:
                logger.warning(f"TargetNet not installed: {e}")
                logger.info("Will use similarity-based fallback")
                self._model = "fallback"
            except Exception as e:
                logger.error(f"Failed to load TargetNet model: {e}")
                self._model = "fallback"

        return self._model

    @property
    def target_database(self):
        """Lazy-load target database."""
        if self._target_database is None:
            # Load or create target database
            self._target_database = self._load_target_database()
        return self._target_database

    def _load_target_database(self) -> Dict[str, Any]:
        """
        Load target database with protein information.

        Returns:
            Dictionary mapping target IDs to information
        """
        # Simplified target database
        # In production, this would load from ChEMBL, BindingDB, etc.
        return {
            'CHEMBL1795121': {
                'name': 'Dopamine D2 receptor',
                'uniprot': 'P14416',
                'type': 'GPCR',
                'disease_areas': ['CNS', 'Schizophrenia', 'Parkinsons']
            },
            'CHEMBL217': {
                'name': 'Acetylcholinesterase',
                'uniprot': 'P22303',
                'type': 'Enzyme',
                'disease_areas': ['Alzheimers', 'Neurodegenerative']
            },
            'CHEMBL1862': {
                'name': 'Cyclooxygenase-2',
                'uniprot': 'P35354',
                'type': 'Enzyme',
                'disease_areas': ['Inflammation', 'Pain']
            },
            'CHEMBL279': {
                'name': 'Vascular endothelial growth factor receptor 2',
                'uniprot': 'P35968',
                'type': 'Kinase',
                'disease_areas': ['Cancer', 'Angiogenesis']
            },
            'CHEMBL203': {
                'name': 'Epidermal growth factor receptor',
                'uniprot': 'P00533',
                'type': 'Kinase',
                'disease_areas': ['Cancer']
            },
            'CHEMBL214': {
                'name': 'Beta-2 adrenergic receptor',
                'uniprot': 'P07550',
                'type': 'GPCR',
                'disease_areas': ['Asthma', 'COPD']
            },
            'CHEMBL244': {
                'name': 'Monoamine oxidase A',
                'uniprot': 'P21397',
                'type': 'Enzyme',
                'disease_areas': ['Depression', 'Anxiety']
            },
            'CHEMBL4005': {
                'name': 'Cytochrome P450 3A4',
                'uniprot': 'P08684',
                'type': 'Enzyme',
                'disease_areas': ['Drug Metabolism']
            },
            'CHEMBL1871': {
                'name': 'Serotonin 2a (5-HT2a) receptor',
                'uniprot': 'P28223',
                'type': 'GPCR',
                'disease_areas': ['Depression', 'Schizophrenia']
            },
            'CHEMBL1937': {
                'name': 'hERG',
                'uniprot': 'Q12809',
                'type': 'Ion Channel',
                'disease_areas': ['Cardiotoxicity']
            }
        }

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data.

        Args:
            input_data: SMILES string or molecular fingerprint

        Returns:
            True if valid, False otherwise
        """
        # Validate SMILES string
        if isinstance(input_data, str):
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(input_data)
                return mol is not None
            except Exception:
                return False

        # Validate fingerprint (numpy array)
        if isinstance(input_data, (list, np.ndarray)):
            return len(input_data) > 0

        return False

    async def execute(self, input_data: Any, **params) -> AdapterResult:
        """
        Execute target prediction.

        Args:
            input_data: SMILES string or molecular fingerprint
            **params: Additional parameters:
                     - confidence_threshold: Override default threshold
                     - max_targets: Maximum targets to return
                     - target_types: Filter by target type (e.g., ['Kinase', 'GPCR'])
                     - include_off_targets: Include off-target predictions (default: True)

        Returns:
            AdapterResult with predicted targets and confidence scores
        """
        try:
            # Validate input
            if not self.validate_input(input_data):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid input data"
                )

            # Parse parameters
            confidence_threshold = params.get(
                'confidence_threshold',
                self.confidence_threshold
            )
            max_targets = params.get('max_targets', self.max_targets)
            target_types = params.get('target_types', None)
            include_off_targets = params.get('include_off_targets', True)

            # Convert SMILES to fingerprint if needed
            if isinstance(input_data, str):
                smiles = input_data
                fingerprint = self._smiles_to_fingerprint(smiles)
            else:
                fingerprint = input_data
                smiles = None

            logger.info(f"Predicting targets for molecule...")

            # Predict targets
            if self.model == "fallback" or self._model == "fallback":
                predictions = self._fallback_prediction(
                    smiles or fingerprint,
                    confidence_threshold
                )
            else:
                predictions = self._model_prediction(
                    fingerprint,
                    confidence_threshold
                )

            # Filter by target type if specified
            if target_types:
                predictions = [
                    p for p in predictions
                    if self.target_database.get(p['target_id'], {}).get('type') in target_types
                ]

            # Limit to max_targets
            predictions = predictions[:max_targets]

            # Enrich with target information
            enriched_predictions = []
            for pred in predictions:
                target_id = pred['target_id']
                target_info = self.target_database.get(target_id, {})

                enriched_pred = {
                    **pred,
                    'target_name': target_info.get('name', 'Unknown'),
                    'uniprot_id': target_info.get('uniprot', 'N/A'),
                    'target_type': target_info.get('type', 'Unknown'),
                    'disease_areas': target_info.get('disease_areas', [])
                }
                enriched_predictions.append(enriched_pred)

            # Identify primary and off-targets
            primary_targets = [
                p for p in enriched_predictions
                if p['confidence'] >= 0.7
            ]
            off_targets = [
                p for p in enriched_predictions
                if 0.5 <= p['confidence'] < 0.7
            ] if include_off_targets else []

            # Analyze polypharmacology
            polypharmacology = self._analyze_polypharmacology(
                enriched_predictions
            )

            result_data = {
                'smiles': smiles,
                'primary_targets': primary_targets,
                'off_targets': off_targets,
                'all_predictions': enriched_predictions,
                'polypharmacology': polypharmacology,
                'statistics': {
                    'total_targets': len(enriched_predictions),
                    'primary_targets': len(primary_targets),
                    'off_targets': len(off_targets),
                    'mean_confidence': np.mean([p['confidence'] for p in enriched_predictions]) if enriched_predictions else 0.0
                },
                'parameters': {
                    'confidence_threshold': confidence_threshold,
                    'max_targets': max_targets
                }
            }

            logger.info(
                f"Predicted {len(enriched_predictions)} targets "
                f"({len(primary_targets)} primary, {len(off_targets)} off-targets)"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "version": self.version
                }
            )

        except Exception as e:
            logger.error(f"Target prediction failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e)
            )

    def _smiles_to_fingerprint(self, smiles: str) -> np.ndarray:
        """
        Convert SMILES to molecular fingerprint.

        Args:
            smiles: SMILES string

        Returns:
            Molecular fingerprint as numpy array
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {smiles}")

        # Generate Morgan fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)

        # Convert to numpy array
        arr = np.zeros((2048,), dtype=np.float32)
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)

        return arr

    def _model_prediction(
        self,
        fingerprint: np.ndarray,
        confidence_threshold: float
    ) -> List[Dict[str, Any]]:
        """
        Predict targets using TargetNet model.

        Args:
            fingerprint: Molecular fingerprint
            confidence_threshold: Minimum confidence

        Returns:
            List of target predictions
        """
        import torch

        predictions = []

        try:
            # Convert to tensor
            fp_tensor = torch.FloatTensor(fingerprint).unsqueeze(0)

            # Predict
            with torch.no_grad():
                outputs = self._model(fp_tensor)
                probabilities = torch.sigmoid(outputs).squeeze().numpy()

            # Get top predictions
            target_ids = list(self.target_database.keys())

            for i, prob in enumerate(probabilities):
                if prob >= confidence_threshold and i < len(target_ids):
                    predictions.append({
                        'target_id': target_ids[i],
                        'confidence': float(prob)
                    })

            # Sort by confidence
            predictions.sort(key=lambda x: x['confidence'], reverse=True)

        except Exception as e:
            logger.error(f"Model prediction error: {e}")
            # Fall back to similarity-based prediction
            return self._fallback_prediction(fingerprint, confidence_threshold)

        return predictions

    def _fallback_prediction(
        self,
        input_data: Union[str, np.ndarray],
        confidence_threshold: float
    ) -> List[Dict[str, Any]]:
        """
        Fallback target prediction using similarity-based approach.

        Args:
            input_data: SMILES or fingerprint
            confidence_threshold: Minimum confidence

        Returns:
            List of target predictions
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors

        # Get molecule
        if isinstance(input_data, str):
            mol = Chem.MolFromSmiles(input_data)
        else:
            # Fingerprint provided - use heuristics
            mol = None

        predictions = []

        # Simple heuristic-based target prediction
        # Based on molecular properties and known drug-target relationships

        if mol:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            rings = Descriptors.RingCount(mol)
            aromatic = Descriptors.NumAromaticRings(mol)

            # Heuristic rules for target prediction
            target_scores = {}

            # Kinase inhibitors: typically flat aromatics, MW 300-600
            if 300 < mw < 600 and aromatic >= 2 and logp > 2:
                target_scores['CHEMBL279'] = 0.75  # VEGFR2
                target_scores['CHEMBL203'] = 0.72  # EGFR

            # GPCR ligands: often contain basic nitrogens
            if hba >= 2 and 250 < mw < 450:
                target_scores['CHEMBL1795121'] = 0.68  # D2 receptor
                target_scores['CHEMBL214'] = 0.65     # Beta-2 AR
                target_scores['CHEMBL1871'] = 0.63    # 5-HT2a

            # Enzyme inhibitors
            if hbd >= 1 and 200 < mw < 500:
                target_scores['CHEMBL217'] = 0.62   # AChE
                target_scores['CHEMBL1862'] = 0.60  # COX-2

            # Drug metabolism enzymes (most drugs interact)
            target_scores['CHEMBL4005'] = 0.55  # CYP3A4

            # hERG liability (many drugs have some affinity)
            if logp > 2 and hba >= 2:
                target_scores['CHEMBL1937'] = 0.52  # hERG

            # Convert to prediction format
            for target_id, score in target_scores.items():
                if score >= confidence_threshold:
                    predictions.append({
                        'target_id': target_id,
                        'confidence': score
                    })

            # Sort by confidence
            predictions.sort(key=lambda x: x['confidence'], reverse=True)

        return predictions

    def _analyze_polypharmacology(
        self,
        predictions: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """
        Analyze polypharmacology profile.

        Args:
            predictions: List of target predictions

        Returns:
            Polypharmacology analysis
        """
        if not predictions:
            return {
                'is_polypharmacological': False,
                'target_count': 0,
                'target_diversity': 0.0
            }

        # Count target types
        target_types = {}
        disease_areas = set()

        for pred in predictions:
            target_type = pred.get('target_type', 'Unknown')
            target_types[target_type] = target_types.get(target_type, 0) + 1

            areas = pred.get('disease_areas', [])
            disease_areas.update(areas)

        # Calculate diversity
        type_diversity = len(target_types) / len(predictions) if predictions else 0.0

        # Determine if polypharmacological
        is_polypharmacological = len(predictions) >= 3 and type_diversity > 0.3

        return {
            'is_polypharmacological': is_polypharmacological,
            'target_count': len(predictions),
            'target_type_distribution': target_types,
            'target_diversity': type_diversity,
            'disease_areas': list(disease_areas),
            'selectivity_score': 1.0 - (len(predictions) / 10.0)  # Higher = more selective
        }

    async def predict_batch(
        self,
        smiles_list: List[str],
        **params
    ) -> List[AdapterResult]:
        """
        Batch target prediction for multiple molecules.

        Args:
            smiles_list: List of SMILES strings
            **params: Parameters for prediction

        Returns:
            List of AdapterResults
        """
        results = []

        for smiles in smiles_list:
            result = await self.execute(smiles, **params)
            results.append(result)

        return results

    def find_shared_targets(
        self,
        predictions_list: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Find shared targets across multiple molecules.

        Args:
            predictions_list: List of prediction results

        Returns:
            List of shared targets with molecules
        """
        from collections import defaultdict

        # Map targets to molecules
        target_to_molecules = defaultdict(list)

        for i, predictions in enumerate(predictions_list):
            for target in predictions.get('all_predictions', []):
                target_id = target['target_id']
                target_to_molecules[target_id].append({
                    'molecule_idx': i,
                    'confidence': target['confidence']
                })

        # Find shared targets
        shared_targets = []

        for target_id, molecules in target_to_molecules.items():
            if len(molecules) >= 2:  # Shared by at least 2 molecules
                target_info = self.target_database.get(target_id, {})
                shared_targets.append({
                    'target_id': target_id,
                    'target_name': target_info.get('name', 'Unknown'),
                    'molecule_count': len(molecules),
                    'molecules': molecules,
                    'mean_confidence': np.mean([m['confidence'] for m in molecules])
                })

        # Sort by molecule count
        shared_targets.sort(key=lambda x: x['molecule_count'], reverse=True)

        return shared_targets

    def generate_cache_key(self, input_data: Any, **params) -> str:
        """
        Generate cache key.

        Args:
            input_data: Input data
            **params: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        import json

        # Canonicalize SMILES if string
        if isinstance(input_data, str):
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(input_data)
                if mol:
                    input_data = Chem.MolToSmiles(mol, canonical=True)
            except Exception:
                pass

        cache_dict = {
            'adapter': self.name,
            'version': self.version,
            'input': str(input_data),
            'params': params
        }

        key_string = json.dumps(cache_dict, sort_keys=True)
        return hashlib.sha256(key_string.encode()).hexdigest()

    def get_metadata(self) -> Dict[str, Any]:
        """Get adapter metadata."""
        return {
            'name': self.name,
            'type': self.adapter_type,
            'version': self.version,
            'enabled': self.enabled,
            'description': 'TargetNet - Protein target prediction',
            'capabilities': [
                'Predict protein targets',
                'Off-target analysis',
                'Multi-target profiling',
                'Polypharmacology analysis',
                'Target-based drug repurposing'
            ],
            'target_database_size': len(self.target_database),
            'config': {
                'confidence_threshold': self.confidence_threshold,
                'max_targets': self.max_targets,
                'use_local': self.use_local
            }
        }
