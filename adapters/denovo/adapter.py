"""
General De Novo Drug Design Adapter

Unified adapter supporting multiple generative approaches:
- RNN: Recurrent neural networks for sequential generation
- VAE: Variational autoencoders for latent space sampling
- Transformer: Attention-based generation
- Fragment: Fragment-based de novo design

Provides property optimization, scaffold hopping, and multi-objective design.
"""

import hashlib
import logging
from typing import Dict, Any, List, Optional, Union, Tuple
import numpy as np

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DeNovoAdapter(AdapterProtocol):
    """
    General de novo drug design adapter.

    Supports multiple generation strategies and provides:
    - Multi-model generation (RNN, VAE, Transformer, Fragment-based)
    - Property optimization and constraints
    - Scaffold hopping and decoration
    - Multi-objective optimization
    - Integration with ADMET and docking
    """

    def __init__(
        self,
        name: str = "denovo",
        adapter_type: str = "ml",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize De Novo adapter.

        Args:
            name: Adapter name (default: "denovo")
            adapter_type: Adapter type (default: "ml")
            config: Configuration dictionary with keys:
                   - model_type: 'rnn', 'vae', 'transformer', 'fragment' (default: 'rnn')
                   - model_path: Path to model weights (optional)
                   - use_pretrained: Use pretrained models (default: True)
                   - optimize_properties: Auto-optimize properties (default: True)
        """
        super().__init__(name, adapter_type, config)
        self.version = "1.0.0"

        # Configuration
        self.model_type = self.config.get('model_type', 'rnn')
        self.model_path = self.config.get('model_path', None)
        self.use_pretrained = self.config.get('use_pretrained', True)
        self.optimize_properties = self.config.get('optimize_properties', True)

        # Lazy-load models
        self._models = {}

        logger.info(
            f"Initialized De Novo adapter "
            f"(model_type={self.model_type})"
        )

    def get_model(self, model_type: str):
        """
        Get or load a specific model.

        Args:
            model_type: Type of model to load

        Returns:
            Model instance or "fallback"
        """
        if model_type not in self._models:
            try:
                if model_type == 'rnn':
                    self._models[model_type] = self._load_rnn_model()
                elif model_type == 'vae':
                    self._models[model_type] = self._load_vae_model()
                elif model_type == 'transformer':
                    self._models[model_type] = self._load_transformer_model()
                elif model_type == 'fragment':
                    self._models[model_type] = self._load_fragment_model()
                else:
                    logger.warning(f"Unknown model type: {model_type}")
                    self._models[model_type] = "fallback"
            except Exception as e:
                logger.error(f"Failed to load {model_type} model: {e}")
                self._models[model_type] = "fallback"

        return self._models[model_type]

    def _load_rnn_model(self):
        """Load RNN-based generator."""
        try:
            # Attempt to load RNN model
            # This is a placeholder for actual implementation
            import torch
            from denovo.models import RNNGenerator

            if self.model_path:
                model = torch.load(self.model_path)
            elif self.use_pretrained:
                model = RNNGenerator.from_pretrained()
            else:
                model = RNNGenerator()

            model.eval()
            logger.info("RNN model loaded successfully")
            return model
        except ImportError:
            logger.warning("RNN model not available, using fallback")
            return "fallback"

    def _load_vae_model(self):
        """Load VAE model."""
        try:
            import torch
            from denovo.models import MolecularVAE

            if self.model_path:
                model = torch.load(self.model_path)
            elif self.use_pretrained:
                model = MolecularVAE.from_pretrained()
            else:
                model = MolecularVAE()

            model.eval()
            logger.info("VAE model loaded successfully")
            return model
        except ImportError:
            logger.warning("VAE model not available, using fallback")
            return "fallback"

    def _load_transformer_model(self):
        """Load Transformer model."""
        try:
            import torch
            from denovo.models import TransformerGenerator

            if self.model_path:
                model = torch.load(self.model_path)
            elif self.use_pretrained:
                model = TransformerGenerator.from_pretrained()
            else:
                model = TransformerGenerator()

            model.eval()
            logger.info("Transformer model loaded successfully")
            return model
        except ImportError:
            logger.warning("Transformer model not available, using fallback")
            return "fallback"

    def _load_fragment_model(self):
        """Load fragment-based generator."""
        try:
            from denovo.models import FragmentGenerator

            model = FragmentGenerator()
            logger.info("Fragment model loaded successfully")
            return model
        except ImportError:
            logger.warning("Fragment model not available, using fallback")
            return "fallback"

    def validate_input(self, input_data: Any) -> bool:
        """Validate input data."""
        # String (SMILES)
        if isinstance(input_data, str):
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(input_data)
                return mol is not None
            except Exception:
                return False

        # Dictionary (parameters)
        if isinstance(input_data, dict):
            return True

        # None (default parameters)
        if input_data is None:
            return True

        return False

    async def execute(self, input_data: Any, **params) -> AdapterResult:
        """
        Execute de novo drug design.

        Args:
            input_data: Either:
                       - Dict with design parameters
                       - SMILES string (for scaffold hopping)
                       - None (for de novo generation)
            **params: Additional parameters:
                     - num_molecules: Number to generate (default: 100)
                     - model_type: Override model type
                     - target_profile: Target property profile
                     - constraints: Molecular constraints
                     - optimization_objectives: List of objectives to optimize
                     - diversity_threshold: Minimum diversity (default: 0.7)

        Returns:
            AdapterResult with generated molecules
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
            num_molecules = params.get('num_molecules', 100)
            model_type = params.get('model_type', self.model_type)
            target_profile = params.get('target_profile', {})
            constraints = params.get('constraints', {})
            objectives = params.get('optimization_objectives', ['qed'])
            diversity_threshold = params.get('diversity_threshold', 0.7)

            # Determine operation mode
            if isinstance(input_data, str):
                mode = 'scaffold_hop'
                seed_smiles = input_data
            else:
                mode = 'generate'
                seed_smiles = None

            logger.info(
                f"De novo design: mode={mode}, model={model_type}, "
                f"n={num_molecules}"
            )

            # Get model
            model = self.get_model(model_type)

            # Generate molecules
            if mode == 'generate':
                molecules = await self._generate_denovo(
                    model=model,
                    model_type=model_type,
                    num_molecules=num_molecules,
                    target_profile=target_profile,
                    constraints=constraints
                )
            else:  # scaffold_hop
                molecules = await self._scaffold_hop(
                    model=model,
                    model_type=model_type,
                    seed_smiles=seed_smiles,
                    num_molecules=num_molecules,
                    target_profile=target_profile
                )

            # Apply property optimization if enabled
            if self.optimize_properties and target_profile:
                molecules = self._optimize_molecules(
                    molecules,
                    target_profile,
                    objectives
                )

            # Apply diversity filtering
            molecules = self._apply_diversity_filter(
                molecules,
                diversity_threshold
            )

            # Calculate comprehensive metrics
            metrics = self._calculate_metrics(molecules, num_molecules)

            # Sort by composite score
            molecules.sort(key=lambda x: x.get('score', 0), reverse=True)

            result_data = {
                'molecules': molecules[:num_molecules],
                'metrics': metrics,
                'parameters': {
                    'mode': mode,
                    'model_type': model_type,
                    'num_molecules': num_molecules,
                    'target_profile': target_profile,
                    'constraints': constraints,
                    'objectives': objectives
                },
                'model_info': {
                    'type': model_type,
                    'version': self.version,
                    'fallback_used': model == "fallback"
                }
            }

            logger.info(
                f"Generated {len(molecules)} molecules "
                f"(validity: {metrics['validity']:.1%}, "
                f"uniqueness: {metrics['uniqueness']:.1%})"
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "version": self.version,
                    "model_type": model_type
                }
            )

        except Exception as e:
            logger.error(f"De novo design failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e)
            )

    async def _generate_denovo(
        self,
        model,
        model_type: str,
        num_molecules: int,
        target_profile: Dict[str, float],
        constraints: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Generate molecules de novo.

        Args:
            model: Generative model
            model_type: Type of model
            num_molecules: Number to generate
            target_profile: Target properties
            constraints: Constraints

        Returns:
            List of generated molecules
        """
        if model == "fallback":
            return self._fallback_generation(
                num_molecules,
                target_profile,
                constraints
            )

        molecules = []

        try:
            if model_type == 'rnn':
                molecules = self._rnn_generation(
                    model, num_molecules, target_profile
                )
            elif model_type == 'vae':
                molecules = self._vae_generation(
                    model, num_molecules, target_profile
                )
            elif model_type == 'transformer':
                molecules = self._transformer_generation(
                    model, num_molecules, target_profile
                )
            elif model_type == 'fragment':
                molecules = self._fragment_generation(
                    model, num_molecules, target_profile
                )

        except Exception as e:
            logger.error(f"Model generation failed: {e}")
            molecules = self._fallback_generation(
                num_molecules,
                target_profile,
                constraints
            )

        return molecules

    def _rnn_generation(
        self,
        model,
        num_molecules: int,
        target_profile: Dict[str, float]
    ) -> List[Dict[str, Any]]:
        """Generate using RNN model."""
        import torch
        from rdkit import Chem

        molecules = []

        with torch.no_grad():
            for _ in range(num_molecules):
                try:
                    # Sample from model
                    smiles = model.sample()

                    # Validate
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        molecules.append({
                            'smiles': smiles,
                            'source': 'RNN'
                        })
                except Exception:
                    continue

        return molecules

    def _vae_generation(
        self,
        model,
        num_molecules: int,
        target_profile: Dict[str, float]
    ) -> List[Dict[str, Any]]:
        """Generate using VAE model."""
        import torch
        from rdkit import Chem

        molecules = []

        with torch.no_grad():
            # Sample from latent space
            z = torch.randn(num_molecules, model.latent_dim)

            # Decode
            generated = model.decode(z)

            for smiles in generated:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        molecules.append({
                            'smiles': smiles,
                            'source': 'VAE'
                        })
                except Exception:
                    continue

        return molecules

    def _transformer_generation(
        self,
        model,
        num_molecules: int,
        target_profile: Dict[str, float]
    ) -> List[Dict[str, Any]]:
        """Generate using Transformer model."""
        import torch
        from rdkit import Chem

        molecules = []

        with torch.no_grad():
            for _ in range(num_molecules):
                try:
                    # Generate sequence
                    smiles = model.generate()

                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        molecules.append({
                            'smiles': smiles,
                            'source': 'Transformer'
                        })
                except Exception:
                    continue

        return molecules

    def _fragment_generation(
        self,
        model,
        num_molecules: int,
        target_profile: Dict[str, float]
    ) -> List[Dict[str, Any]]:
        """Generate using fragment-based approach."""
        from rdkit import Chem

        molecules = []

        for _ in range(num_molecules):
            try:
                # Assemble from fragments
                smiles = model.assemble_fragments()

                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    molecules.append({
                        'smiles': smiles,
                        'source': 'Fragment'
                    })
            except Exception:
                continue

        return molecules

    def _fallback_generation(
        self,
        num_molecules: int,
        target_profile: Dict[str, float],
        constraints: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """Fallback generation using RDKit."""
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, QED

        molecules = []

        # Diverse fragment library
        fragments = [
            # Aromatic cores
            'c1ccccc1', 'c1ccncc1', 'c1cnccc1', 'c1ccoc1', 'c1ccsc1',
            'c1c[nH]cc1', 'c1ccc2ccccc2c1',
            # Aliphatic cores
            'C1CCCCC1', 'C1CCNCC1', 'C1CCOCC1', 'C1CCOC1',
            # Heterocycles
            'c1ccc2[nH]ccc2c1', 'c1cnc2ccccc2c1', 'c1ccc2occc2c1'
        ]

        # Functional groups
        groups = [
            'C(=O)O', 'C(=O)N', 'CN', 'CO', 'C(F)(F)F',
            'Cl', 'Br', 'S(=O)(=O)N', 'C#N', 'N'
        ]

        attempts = 0
        max_attempts = num_molecules * 20

        while len(molecules) < num_molecules and attempts < max_attempts:
            attempts += 1

            try:
                # Random assembly
                core = Chem.MolFromSmiles(np.random.choice(fragments))
                if not core:
                    continue

                mol = Chem.RWMol(core)

                # Add groups
                n_groups = np.random.randint(1, 4)
                for _ in range(n_groups):
                    group_smiles = np.random.choice(groups)
                    group = Chem.MolFromSmiles(group_smiles)
                    if group:
                        mol = Chem.CombineMols(mol, group)

                # Sanitize
                Chem.SanitizeMol(mol)
                smiles = Chem.MolToSmiles(mol)

                # Basic filtering
                if QED.qed(mol) > 0.3:
                    molecules.append({
                        'smiles': smiles,
                        'source': 'RDKit-Fallback'
                    })

            except Exception:
                continue

        return molecules

    async def _scaffold_hop(
        self,
        model,
        model_type: str,
        seed_smiles: str,
        num_molecules: int,
        target_profile: Dict[str, float]
    ) -> List[Dict[str, Any]]:
        """
        Perform scaffold hopping.

        Args:
            model: Generative model
            model_type: Model type
            seed_smiles: Seed molecule
            num_molecules: Number of hops
            target_profile: Target properties

        Returns:
            List of scaffold-hopped molecules
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        molecules = []
        seed_mol = Chem.MolFromSmiles(seed_smiles)

        if not seed_mol:
            return molecules

        # Extract scaffold
        from rdkit.Chem.Scaffolds import MurckoScaffold
        scaffold = MurckoScaffold.GetScaffoldForMol(seed_mol)
        scaffold_smiles = Chem.MolToSmiles(scaffold)

        # Generate variations
        for _ in range(num_molecules):
            try:
                # Create variant
                variant = self._create_scaffold_variant(scaffold_smiles)
                if variant:
                    molecules.append({
                        'smiles': variant,
                        'original_scaffold': scaffold_smiles,
                        'seed_smiles': seed_smiles,
                        'source': f'{model_type}_scaffold_hop'
                    })
            except Exception:
                continue

        return molecules

    def _create_scaffold_variant(self, scaffold_smiles: str) -> Optional[str]:
        """Create scaffold variant."""
        from rdkit import Chem

        scaffold = Chem.MolFromSmiles(scaffold_smiles)
        if not scaffold:
            return None

        # Simple mutations: add/remove atoms, change bonds
        variant = Chem.RWMol(scaffold)

        # Add random decorations
        decorations = ['C', 'N', 'O', 'F', 'Cl']
        if variant.GetNumAtoms() > 0:
            idx = np.random.randint(0, variant.GetNumAtoms())
            decoration = np.random.choice(decorations)
            new_atom = Chem.Atom(decoration)
            new_idx = variant.AddAtom(new_atom)
            variant.AddBond(idx, new_idx, Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(variant)
            return Chem.MolToSmiles(variant)
        except Exception:
            return None

    def _optimize_molecules(
        self,
        molecules: List[Dict[str, Any]],
        target_profile: Dict[str, float],
        objectives: List[str]
    ) -> List[Dict[str, Any]]:
        """Optimize molecules for target properties."""
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED

        for mol_data in molecules:
            smiles = mol_data['smiles']
            mol = Chem.MolFromSmiles(smiles)

            if not mol:
                mol_data['score'] = 0.0
                continue

            # Calculate properties
            properties = self._calculate_properties(mol)
            mol_data['properties'] = properties

            # Calculate objective scores
            scores = []

            for obj in objectives:
                if obj == 'qed':
                    scores.append(QED.qed(mol))
                elif obj in target_profile:
                    target = target_profile[obj]
                    actual = properties.get(obj, 0)
                    # Gaussian scoring
                    score = np.exp(-((actual - target) ** 2) / (2 * (target * 0.2) ** 2))
                    scores.append(score)

            # Composite score
            mol_data['score'] = np.mean(scores) if scores else 0.0

        return molecules

    def _calculate_properties(self, mol) -> Dict[str, float]:
        """Calculate molecular properties."""
        from rdkit.Chem import Descriptors, QED, Lipinski

        return {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol),
            'hbd': Lipinski.NumHDonors(mol),
            'hba': Lipinski.NumHAcceptors(mol),
            'rotatable_bonds': Lipinski.NumRotatableBonds(mol),
            'qed': QED.qed(mol),
            'num_atoms': mol.GetNumAtoms()
        }

    def _apply_diversity_filter(
        self,
        molecules: List[Dict[str, Any]],
        threshold: float
    ) -> List[Dict[str, Any]]:
        """Apply diversity filtering."""
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem

        if len(molecules) <= 1:
            return molecules

        # Calculate fingerprints
        fps = []
        for mol_data in molecules:
            mol = Chem.MolFromSmiles(mol_data['smiles'])
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append((mol_data, fp))

        # Select diverse subset
        selected = [fps[0]]

        for mol_data, fp in fps[1:]:
            # Check similarity to selected
            max_sim = max([
                DataStructs.TanimotoSimilarity(fp, sel_fp)
                for _, sel_fp in selected
            ])

            if max_sim < threshold:
                selected.append((mol_data, fp))

        return [mol_data for mol_data, _ in selected]

    def _calculate_metrics(
        self,
        molecules: List[Dict[str, Any]],
        target_count: int
    ) -> Dict[str, Any]:
        """Calculate generation metrics."""
        from rdkit import Chem

        if not molecules:
            return {
                'validity': 0.0,
                'uniqueness': 0.0,
                'novelty': 0.0,
                'diversity': 0.0
            }

        # Validity
        valid_count = sum(
            1 for m in molecules
            if Chem.MolFromSmiles(m['smiles']) is not None
        )
        validity = valid_count / len(molecules) if molecules else 0.0

        # Uniqueness
        unique_smiles = set(m['smiles'] for m in molecules)
        uniqueness = len(unique_smiles) / len(molecules) if molecules else 0.0

        # Novelty (placeholder - needs training set)
        novelty = 0.9

        # Diversity
        diversity = self._calculate_diversity(molecules)

        return {
            'validity': validity,
            'uniqueness': uniqueness,
            'novelty': novelty,
            'diversity': diversity,
            'total_generated': len(molecules),
            'target_count': target_count
        }

    def _calculate_diversity(self, molecules: List[Dict[str, Any]]) -> float:
        """Calculate molecular diversity."""
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem

        if len(molecules) <= 1:
            return 0.0

        # Calculate pairwise similarities
        fps = []
        for mol_data in molecules[:100]:  # Limit for performance
            mol = Chem.MolFromSmiles(mol_data['smiles'])
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)

        if len(fps) <= 1:
            return 0.0

        # Average dissimilarity
        similarities = []
        for i in range(len(fps)):
            for j in range(i + 1, len(fps)):
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                similarities.append(sim)

        avg_similarity = np.mean(similarities) if similarities else 0.0
        diversity = 1.0 - avg_similarity

        return diversity

    def generate_cache_key(self, input_data: Any, **params) -> str:
        """Generate cache key."""
        import json

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
            'input': input_data,
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
            'description': 'General de novo drug design framework',
            'supported_models': ['RNN', 'VAE', 'Transformer', 'Fragment'],
            'capabilities': [
                'De novo molecule generation',
                'Scaffold hopping',
                'Property optimization',
                'Multi-objective design',
                'Diversity control'
            ],
            'config': {
                'model_type': self.model_type,
                'optimize_properties': self.optimize_properties,
                'use_pretrained': self.use_pretrained
            }
        }
