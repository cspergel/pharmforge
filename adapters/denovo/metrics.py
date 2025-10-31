"""
Molecular Generation Metrics

Comprehensive metrics for evaluating de novo generated molecules:
- Validity: Percentage of valid molecules
- Uniqueness: Percentage of unique structures
- Novelty: Percentage of molecules not in training set
- Diversity: Structural diversity of generated set
- Drug-likeness: QED and other drug-like properties
- Synthetic Accessibility: SA score
"""

import logging
from typing import List, Dict, Any, Optional, Set, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class GenerationMetrics:
    """
    Comprehensive metrics for molecular generation evaluation.

    Calculates standard benchmarking metrics used in de novo design:
    - Validity, Uniqueness, Novelty (VUN metrics)
    - Diversity (internal and external)
    - Drug-likeness scores
    - Property distributions
    """

    def __init__(
        self,
        reference_set: Optional[List[str]] = None,
        calculate_synthetic_accessibility: bool = True
    ):
        """
        Initialize metrics calculator.

        Args:
            reference_set: Reference SMILES set for novelty calculation
            calculate_synthetic_accessibility: Calculate SA scores (default: True)
        """
        self.reference_set = set(reference_set) if reference_set else set()
        self.calculate_sa = calculate_synthetic_accessibility

    def calculate_all_metrics(
        self,
        generated_smiles: List[str],
        verbose: bool = True
    ) -> Dict[str, Any]:
        """
        Calculate all metrics for generated molecules.

        Args:
            generated_smiles: List of generated SMILES strings
            verbose: Print detailed results (default: True)

        Returns:
            Dictionary containing all metrics
        """
        if not generated_smiles:
            logger.warning("Empty generation set")
            return self._empty_metrics()

        # Basic metrics
        validity = self.calculate_validity(generated_smiles)
        valid_smiles = [s for s in generated_smiles if self._is_valid(s)]

        uniqueness = self.calculate_uniqueness(valid_smiles)
        unique_smiles = list(set(valid_smiles))

        novelty = self.calculate_novelty(unique_smiles)

        # Diversity
        diversity_internal = self.calculate_diversity(unique_smiles)

        # Drug-likeness
        druglikeness = self.calculate_druglikeness(unique_smiles)

        # Property statistics
        property_stats = self.calculate_property_statistics(unique_smiles)

        # Synthetic accessibility
        sa_stats = None
        if self.calculate_sa:
            sa_stats = self.calculate_sa_statistics(unique_smiles)

        # Compile results
        metrics = {
            'validity': validity,
            'uniqueness': uniqueness,
            'novelty': novelty,
            'diversity_internal': diversity_internal,
            'druglikeness': druglikeness,
            'property_statistics': property_stats,
            'sa_statistics': sa_stats,
            'total_generated': len(generated_smiles),
            'valid_count': len(valid_smiles),
            'unique_count': len(unique_smiles)
        }

        if verbose:
            self._print_metrics(metrics)

        return metrics

    def calculate_validity(self, smiles_list: List[str]) -> float:
        """
        Calculate validity: percentage of valid molecules.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Validity ratio (0-1)
        """
        if not smiles_list:
            return 0.0

        valid_count = sum(1 for s in smiles_list if self._is_valid(s))
        return valid_count / len(smiles_list)

    def calculate_uniqueness(self, smiles_list: List[str]) -> float:
        """
        Calculate uniqueness: percentage of unique structures.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Uniqueness ratio (0-1)
        """
        if not smiles_list:
            return 0.0

        unique_count = len(set(smiles_list))
        return unique_count / len(smiles_list)

    def calculate_novelty(self, smiles_list: List[str]) -> float:
        """
        Calculate novelty: percentage not in reference set.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Novelty ratio (0-1)
        """
        if not smiles_list:
            return 0.0

        if not self.reference_set:
            logger.warning("No reference set provided, assuming all novel")
            return 1.0

        # Canonicalize for comparison
        canonical_generated = set()
        for smiles in smiles_list:
            canonical = self._canonicalize(smiles)
            if canonical:
                canonical_generated.add(canonical)

        # Count novel molecules
        novel_count = sum(
            1 for s in canonical_generated
            if s not in self.reference_set
        )

        return novel_count / len(canonical_generated) if canonical_generated else 0.0

    def calculate_diversity(
        self,
        smiles_list: List[str],
        method: str = 'tanimoto'
    ) -> float:
        """
        Calculate internal diversity of molecule set.

        Args:
            smiles_list: List of SMILES strings
            method: Diversity metric ('tanimoto', 'scaffold', 'mcs')

        Returns:
            Diversity score (0-1, higher = more diverse)
        """
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem

        if len(smiles_list) <= 1:
            return 0.0

        # Limit for performance
        sample = smiles_list[:200] if len(smiles_list) > 200 else smiles_list

        # Calculate fingerprints
        fps = []
        for smiles in sample:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)

        if len(fps) <= 1:
            return 0.0

        # Calculate pairwise similarities
        similarities = []
        for i in range(len(fps)):
            for j in range(i + 1, len(fps)):
                sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                similarities.append(sim)

        # Diversity = 1 - average similarity
        avg_similarity = np.mean(similarities) if similarities else 0.0
        diversity = 1.0 - avg_similarity

        return diversity

    def calculate_druglikeness(
        self,
        smiles_list: List[str]
    ) -> Dict[str, float]:
        """
        Calculate drug-likeness metrics.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Dictionary with drug-likeness scores
        """
        from rdkit import Chem
        from rdkit.Chem import QED, Lipinski

        qed_scores = []
        lipinski_pass = 0

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            # QED score
            qed_scores.append(QED.qed(mol))

            # Lipinski rule of 5
            mw = mol.GetProp('MolWt') if mol.HasProp('MolWt') else Chem.Descriptors.MolWt(mol)
            logp = Chem.Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)

            violations = 0
            if mw > 500:
                violations += 1
            if logp > 5:
                violations += 1
            if hbd > 5:
                violations += 1
            if hba > 10:
                violations += 1

            if violations <= 1:  # Allow 1 violation
                lipinski_pass += 1

        return {
            'mean_qed': np.mean(qed_scores) if qed_scores else 0.0,
            'median_qed': np.median(qed_scores) if qed_scores else 0.0,
            'lipinski_pass_rate': lipinski_pass / len(smiles_list) if smiles_list else 0.0
        }

    def calculate_property_statistics(
        self,
        smiles_list: List[str]
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate property distribution statistics.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Dictionary of property statistics
        """
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski

        properties = {
            'molecular_weight': [],
            'logp': [],
            'tpsa': [],
            'hbd': [],
            'hba': [],
            'rotatable_bonds': [],
            'num_rings': [],
            'num_aromatic_rings': []
        }

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            properties['molecular_weight'].append(Descriptors.MolWt(mol))
            properties['logp'].append(Descriptors.MolLogP(mol))
            properties['tpsa'].append(Descriptors.TPSA(mol))
            properties['hbd'].append(Lipinski.NumHDonors(mol))
            properties['hba'].append(Lipinski.NumHAcceptors(mol))
            properties['rotatable_bonds'].append(Lipinski.NumRotatableBonds(mol))
            properties['num_rings'].append(Lipinski.RingCount(mol))
            properties['num_aromatic_rings'].append(Descriptors.NumAromaticRings(mol))

        # Calculate statistics
        stats = {}
        for prop, values in properties.items():
            if values:
                stats[prop] = {
                    'mean': float(np.mean(values)),
                    'std': float(np.std(values)),
                    'min': float(np.min(values)),
                    'max': float(np.max(values)),
                    'median': float(np.median(values))
                }
            else:
                stats[prop] = {
                    'mean': 0.0,
                    'std': 0.0,
                    'min': 0.0,
                    'max': 0.0,
                    'median': 0.0
                }

        return stats

    def calculate_sa_statistics(
        self,
        smiles_list: List[str]
    ) -> Dict[str, float]:
        """
        Calculate synthetic accessibility statistics.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            SA score statistics
        """
        try:
            from adapters.custom.synthesis.sascore import calculateScore
            from rdkit import Chem

            sa_scores = []

            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    try:
                        sa_score = calculateScore(mol)
                        sa_scores.append(sa_score)
                    except Exception:
                        continue

            if sa_scores:
                return {
                    'mean_sa': float(np.mean(sa_scores)),
                    'median_sa': float(np.median(sa_scores)),
                    'std_sa': float(np.std(sa_scores)),
                    'easy_to_synthesize': sum(1 for s in sa_scores if s <= 3) / len(sa_scores)
                }
            else:
                return self._empty_sa_stats()

        except ImportError:
            logger.warning("SA score calculation not available")
            return self._empty_sa_stats()

    def compare_to_reference(
        self,
        generated_smiles: List[str],
        reference_smiles: List[str]
    ) -> Dict[str, Any]:
        """
        Compare generated molecules to reference set.

        Args:
            generated_smiles: Generated SMILES
            reference_smiles: Reference SMILES

        Returns:
            Comparison metrics
        """
        # Calculate metrics for both sets
        gen_metrics = self.calculate_property_statistics(generated_smiles)
        ref_metrics = self.calculate_property_statistics(reference_smiles)

        # Calculate property distribution similarity
        distribution_similarity = {}

        for prop in gen_metrics.keys():
            gen_mean = gen_metrics[prop]['mean']
            gen_std = gen_metrics[prop]['std']
            ref_mean = ref_metrics[prop]['mean']
            ref_std = ref_metrics[prop]['std']

            # Wasserstein distance approximation
            mean_diff = abs(gen_mean - ref_mean)
            std_diff = abs(gen_std - ref_std)

            # Normalized similarity (0-1, higher = more similar)
            similarity = 1.0 / (1.0 + mean_diff + std_diff)
            distribution_similarity[prop] = similarity

        return {
            'distribution_similarity': distribution_similarity,
            'generated_metrics': gen_metrics,
            'reference_metrics': ref_metrics
        }

    def _is_valid(self, smiles: str) -> bool:
        """Check if SMILES is valid."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None and mol.GetNumAtoms() > 0
        except Exception:
            return False

    def _canonicalize(self, smiles: str) -> Optional[str]:
        """Canonicalize SMILES."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            pass
        return None

    def _empty_metrics(self) -> Dict[str, Any]:
        """Return empty metrics."""
        return {
            'validity': 0.0,
            'uniqueness': 0.0,
            'novelty': 0.0,
            'diversity_internal': 0.0,
            'druglikeness': {
                'mean_qed': 0.0,
                'median_qed': 0.0,
                'lipinski_pass_rate': 0.0
            },
            'property_statistics': {},
            'sa_statistics': self._empty_sa_stats(),
            'total_generated': 0,
            'valid_count': 0,
            'unique_count': 0
        }

    def _empty_sa_stats(self) -> Dict[str, float]:
        """Return empty SA statistics."""
        return {
            'mean_sa': 0.0,
            'median_sa': 0.0,
            'std_sa': 0.0,
            'easy_to_synthesize': 0.0
        }

    def _print_metrics(self, metrics: Dict[str, Any]):
        """Print formatted metrics."""
        print("\n" + "=" * 60)
        print("MOLECULAR GENERATION METRICS")
        print("=" * 60)

        print(f"\nBasic Metrics:")
        print(f"  Validity:    {metrics['validity']:.1%}")
        print(f"  Uniqueness:  {metrics['uniqueness']:.1%}")
        print(f"  Novelty:     {metrics['novelty']:.1%}")
        print(f"  Diversity:   {metrics['diversity_internal']:.1%}")

        print(f"\nDrug-likeness:")
        dl = metrics['druglikeness']
        print(f"  Mean QED:    {dl['mean_qed']:.3f}")
        print(f"  Lipinski:    {dl['lipinski_pass_rate']:.1%}")

        if metrics.get('sa_statistics'):
            sa = metrics['sa_statistics']
            print(f"\nSynthetic Accessibility:")
            print(f"  Mean SA:     {sa['mean_sa']:.2f}")
            print(f"  Easy to synthesize: {sa['easy_to_synthesize']:.1%}")

        print(f"\nCounts:")
        print(f"  Total:       {metrics['total_generated']}")
        print(f"  Valid:       {metrics['valid_count']}")
        print(f"  Unique:      {metrics['unique_count']}")

        print("=" * 60 + "\n")


def calculate_molecular_similarity(
    smiles1: str,
    smiles2: str,
    method: str = 'tanimoto'
) -> float:
    """
    Calculate similarity between two molecules.

    Args:
        smiles1: First SMILES string
        smiles2: Second SMILES string
        method: Similarity method ('tanimoto', 'dice')

    Returns:
        Similarity score (0-1)
    """
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem

    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if not mol1 or not mol2:
        return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)

    if method == 'tanimoto':
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    elif method == 'dice':
        return DataStructs.DiceSimilarity(fp1, fp2)
    else:
        return DataStructs.TanimotoSimilarity(fp1, fp2)
