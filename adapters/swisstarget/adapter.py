"""
SwissTargetPrediction Adapter for PharmForge

Predicts protein targets for small molecules using the SwissTargetPrediction
web service. This tool uses 2D and 3D similarity to known bioactive molecules
to predict potential targets.

References:
    - SwissTargetPrediction: http://swisstargetprediction.ch/
    - Paper: https://academic.oup.com/nar/article/47/W1/W357/5494775
"""

import aiohttp
import asyncio
import logging
from typing import Any, Dict, Optional, List
from urllib.parse import quote
import json

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class SwissTargetAdapter(AdapterProtocol):
    """
    Adapter for SwissTargetPrediction web service

    SwissTargetPrediction estimates the most probable protein targets
    of a small molecule using a combination of 2D and 3D similarity
    measures with known bioactive compounds.

    Features:
        - Target prediction by structural similarity
        - Probability scores for each predicted target
        - Organism-specific predictions (human, mouse, rat, etc.)
        - Pathway and protein family information
        - Known ligands for each target

    Use Cases:
        - Reverse pharmacology (finding targets for compounds)
        - Off-target prediction
        - Drug repurposing
        - Polypharmacology analysis

    Returns:
        - Predicted targets with probability scores
        - Target classification (e.g., enzyme, GPCR, kinase)
        - Organism information
        - Known active compounds for comparison
    """

    BASE_URL = "https://www.swisstargetprediction.ch"

    def __init__(
        self,
        name: str = "swisstarget",
        adapter_type: str = "api",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize SwissTargetPrediction adapter

        Args:
            name: Adapter name (default: "swisstarget")
            adapter_type: Adapter type (default: "api")
            config: Optional configuration dictionary. Supported keys:
                   - timeout: Request timeout in seconds (default: 60)
                   - organism: Target organism (default: "Homo sapiens")
                   - min_probability: Minimum probability threshold (default: 0.0)
                   - max_targets: Maximum targets to return (default: 15)
        """
        super().__init__(name, adapter_type, config or {})
        self.version = "1.0.0"
        self.timeout = self.config.get('timeout', 60)
        self.organism = self.config.get('organism', 'Homo sapiens')
        self.min_probability = self.config.get('min_probability', 0.0)
        self.max_targets = self.config.get('max_targets', 15)

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data

        Args:
            input_data: SMILES string

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        # Basic SMILES validation
        return True

    async def _predict_targets(
        self,
        smiles: str,
        organism: str = "Homo sapiens"
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Query SwissTargetPrediction for target predictions

        Args:
            smiles: SMILES string
            organism: Target organism

        Returns:
            List of predicted targets with probabilities
        """
        try:
            # SwissTargetPrediction web interface workflow:
            # 1. Submit SMILES via POST to prediction endpoint
            # 2. Poll for results (processing takes time)
            # 3. Parse results page

            # Note: SwissTargetPrediction doesn't have a public REST API
            # This implementation would need to:
            # - POST SMILES to web form
            # - Handle redirects/job IDs
            # - Parse HTML results

            # For demonstration, we'll show the expected workflow
            logger.warning(
                "SwissTargetPrediction requires web scraping or official API. "
                "Returning example data structure."
            )

            # Example POST endpoint (actual URL may vary)
            url = f"{self.BASE_URL}/predict"

            # Prepare form data
            data = {
                'smiles': smiles,
                'organism': organism
            }

            timeout = aiohttp.ClientTimeout(total=self.timeout)

            async with aiohttp.ClientSession() as session:
                # Submit prediction request
                async with session.post(url, data=data, timeout=timeout) as response:
                    if response.status == 200:
                        # Parse response (HTML or JSON)
                        # Real implementation would parse the results page
                        text = await response.text()
                        logger.info(f"SwissTarget response: {len(text)} bytes")

                        # Would parse HTML here
                        # For now, return None to use example data
                        return None
                    else:
                        logger.error(f"SwissTarget request failed: {response.status}")
                        return None

        except asyncio.TimeoutError:
            logger.error("SwissTarget request timeout")
            return None
        except Exception as e:
            logger.error(f"Error querying SwissTarget: {e}")
            return None

    def _create_example_targets(self, smiles: str) -> List[Dict[str, Any]]:
        """
        Create example target prediction data

        This demonstrates the expected data format.
        Real implementation would parse SwissTargetPrediction results.

        Args:
            smiles: Query SMILES

        Returns:
            Example target predictions
        """
        # Example predictions showing different target classes
        targets = [
            {
                "target_name": "Cyclooxygenase-2 (COX-2)",
                "uniprot_id": "P35354",
                "target_class": "Enzyme",
                "target_family": "Oxidoreductase",
                "probability": 0.95,
                "known_actives": 1523,
                "organism": "Homo sapiens",
                "chembl_id": "CHEMBL230",
                "pathway": "Arachidonic acid metabolism"
            },
            {
                "target_name": "Cyclooxygenase-1 (COX-1)",
                "uniprot_id": "P23219",
                "target_class": "Enzyme",
                "target_family": "Oxidoreductase",
                "probability": 0.87,
                "known_actives": 987,
                "organism": "Homo sapiens",
                "chembl_id": "CHEMBL221",
                "pathway": "Arachidonic acid metabolism"
            },
            {
                "target_name": "Cannabinoid receptor 1",
                "uniprot_id": "P21554",
                "target_class": "GPCR",
                "target_family": "Class A (Rhodopsin-like)",
                "probability": 0.72,
                "known_actives": 2341,
                "organism": "Homo sapiens",
                "chembl_id": "CHEMBL218",
                "pathway": "Endocannabinoid signaling"
            },
            {
                "target_name": "Prostaglandin E synthase",
                "uniprot_id": "O14684",
                "target_class": "Enzyme",
                "target_family": "Transferase",
                "probability": 0.68,
                "known_actives": 234,
                "organism": "Homo sapiens",
                "chembl_id": "CHEMBL4462",
                "pathway": "Arachidonic acid metabolism"
            },
            {
                "target_name": "Peroxisome proliferator-activated receptor gamma",
                "uniprot_id": "P37231",
                "target_class": "Nuclear Receptor",
                "target_family": "NR1C subfamily",
                "probability": 0.54,
                "known_actives": 1876,
                "organism": "Homo sapiens",
                "chembl_id": "CHEMBL235",
                "pathway": "Lipid metabolism"
            }
        ]

        return targets

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute SwissTargetPrediction query

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters:
                - organism: Target organism (default: from config)
                - min_probability: Filter by minimum probability

        Returns:
            AdapterResult containing predicted targets
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string"
            )

        smiles = input_data
        organism = kwargs.get('organism', self.organism)
        min_prob = kwargs.get('min_probability', self.min_probability)

        try:
            # Query SwissTargetPrediction
            targets = await self._predict_targets(smiles, organism)

            # If API not available, use example data
            if targets is None:
                targets = self._create_example_targets(smiles)

            # Filter by probability
            if min_prob > 0:
                targets = [t for t in targets if t.get('probability', 0) >= min_prob]

            # Limit number of targets
            targets = targets[:self.max_targets]

            # Sort by probability (descending)
            targets = sorted(targets, key=lambda x: x.get('probability', 0), reverse=True)

            # Calculate statistics
            if targets:
                probabilities = [t.get('probability', 0) for t in targets]
                stats = {
                    'total_targets': len(targets),
                    'max_probability': max(probabilities),
                    'mean_probability': sum(probabilities) / len(probabilities),
                    'high_confidence': len([p for p in probabilities if p >= 0.7]),
                    'medium_confidence': len([p for p in probabilities if 0.4 <= p < 0.7]),
                    'low_confidence': len([p for p in probabilities if p < 0.4])
                }
            else:
                stats = None

            # Group by target class
            target_classes = {}
            for target in targets:
                tclass = target.get('target_class', 'Unknown')
                if tclass not in target_classes:
                    target_classes[tclass] = []
                target_classes[tclass].append(target['target_name'])

            result_data = {
                "smiles": smiles,
                "organism": organism,
                "targets": targets,
                "statistics": stats,
                "target_classes": target_classes,
                "source": "SwissTargetPrediction",
                "note": "This is example data. Full implementation requires web scraping or official API."
            }

            logger.info(
                f"SwissTarget prediction complete: {len(targets)} targets "
                f"(max prob: {stats['max_probability']:.2f})" if stats else "No targets"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "organism": organism,
                    "num_targets": len(targets),
                    "source": "SwissTargetPrediction"
                }
            )

        except Exception as e:
            logger.error(f"SwissTarget adapter error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=str(e)
            )

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
            "description": "Predict protein targets for small molecules using SwissTargetPrediction",
            "features": [
                "Target prediction by 2D/3D similarity",
                "Probability scores",
                "Multiple organisms",
                "Target classification",
                "Pathway information"
            ],
            "use_cases": [
                "Reverse pharmacology",
                "Off-target prediction",
                "Drug repurposing",
                "Polypharmacology analysis"
            ],
            "organisms": [
                "Homo sapiens",
                "Mus musculus",
                "Rattus norvegicus",
                "Bos taurus",
                "Canis familiaris"
            ],
            "target_classes": [
                "Enzyme",
                "GPCR",
                "Kinase",
                "Nuclear Receptor",
                "Ion Channel",
                "Transporter"
            ],
            "config": {
                "timeout": self.timeout,
                "organism": self.organism,
                "min_probability": self.min_probability,
                "max_targets": self.max_targets
            },
            "input": {
                "smiles": "Compound SMILES structure"
            },
            "output": {
                "targets": "List of predicted targets with probabilities",
                "statistics": "Prediction statistics and confidence distribution",
                "target_classes": "Targets grouped by class"
            },
            "notes": [
                "Requires web scraping or official API for full functionality",
                "Currently returns example data structure",
                "Production version needs HTML parsing or API integration"
            ]
        }
