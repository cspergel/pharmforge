# Synthetic Accessibility Score (SAScore) Evaluator
# Adapted from RDKit contrib: Ertl & Schuffenhauer, J. Cheminform. 2009
import gzip
import math
import os
import pickle
from typing import Any, Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


class SynthesisAccessibilityEvaluator:
    """
    Computes Synthetic Accessibility Score (SAscore) from RDKit.
    Score is between 1 (easy) and 10 (hard). Combines fragment score, penalties, and symmetry.

    Reference: Ertl P, Schuffenhauer A. J. Cheminform. 2009, 1:8
    """

    def __init__(self, fscore_path: str = None):
        # Default path: adapters/custom/resources/fpscores.pkl.gz
        if fscore_path is None:
            import os
            current_dir = os.path.dirname(os.path.abspath(__file__))
            fscore_path = os.path.join(current_dir, "..", "resources", "fpscores.pkl.gz")
        self.fscore_path = os.path.abspath(fscore_path)
        self._fscores = self._load_fragment_scores()
        self._fpgen = GetMorganGenerator(radius=2)

    def _load_fragment_scores(self) -> Dict[int, float]:
        if not os.path.exists(self.fscore_path):
            raise FileNotFoundError(f"Fragment score file not found: {self.fscore_path}")
        try:
            with gzip.open(self.fscore_path, "rb") as f:
                data = pickle.load(f)
        except OSError:
            with open(self.fscore_path, "rb") as f:
                data = pickle.load(f)

        scores = {}
        for entry in data:
            for fragment in entry[1:]:
                scores[fragment] = float(entry[0])
        print(f"[SAScore] Loaded {len(scores)} fragment scores.")
        return scores

    def evaluate(self, smiles: str, **kwargs) -> Dict[str, Any]:
        """
        Evaluate synthetic accessibility for a single molecule.

        Args:
            smiles: SMILES string of molecule

        Returns:
            Dictionary with score, tags, reason, summary, raw data
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol or mol.GetNumAtoms() == 0:
            return {
                "score": 10.0,  # Maximum difficulty
                "tags": ["InvalidMolecule"],
                "reason": "Invalid SMILES string",
                "summary": {},
                "raw": {}
            }

        try:
            sascore, fallback_count, total_fragments = self._calculate_sascore(mol)
            tag = (
                "EasySynthesis"
                if sascore < 3
                else "HardToSynthesize"
                if sascore > 7
                else "ModerateSynthesis"
            )
            difficulty = "low" if sascore < 3 else "high" if sascore > 7 else "moderate"

            return {
                "score": round(sascore, 3),
                "tags": [tag],
                "reason": f"SAscore = {sascore:.2f} (1 = easy, 10 = hard)",
                "summary": {
                    "sascore": round(sascore, 3),
                    "synthesis_difficulty": difficulty,
                },
                "raw": {
                    "sascore": sascore,
                    "synthesis_difficulty": difficulty,
                    "fallback_fragments": fallback_count,
                    "total_fragments": total_fragments,
                },
            }

        except Exception as e:
            return {
                "score": 10.0,
                "tags": ["EvaluationError"],
                "reason": f"SAScore computation failed: {str(e)}",
                "summary": {},
                "raw": {}
            }

    def evaluate_batch(self, smiles_list: List[str], **kwargs) -> List[Dict[str, Any]]:
        """Evaluate multiple molecules"""
        return [self.evaluate(smi, **kwargs) for smi in smiles_list]

    def _calculate_sascore(self, mol) -> (float, int, int):
        sfp = self._fpgen.GetSparseCountFingerprint(mol)
        nze = sfp.GetNonzeroElements()
        n_fragments = sum(nze.values())

        score_sum = 0.0
        fallback_count = 0
        for fid, count in nze.items():
            if fid not in self._fscores:
                fallback_count += 1
            score_sum += self._fscores.get(fid, -4) * count

        score1 = score_sum / max(n_fragments, 1)

        n_atoms = mol.GetNumAtoms()
        n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_bridge = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        n_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        n_macro = sum(1 for ring in mol.GetRingInfo().AtomRings() if len(ring) > 8)

        score2 = (
            -(n_atoms**1.005 - n_atoms)
            - math.log10(n_chiral + 1)
            - math.log10(n_spiro + 1)
            - math.log10(n_bridge + 1)
            - (math.log10(2) if n_macro > 0 else 0)
        )

        num_bits = len(nze)
        score3 = math.log(n_atoms / num_bits) * 0.5 if n_atoms > num_bits > 0 else 0

        raw_score = score1 + score2 + score3

        min_score, max_score = -4.0, 2.5
        sascore = 11.0 - ((raw_score - min_score + 1) / (max_score - min_score)) * 9.0
        if sascore > 8:
            sascore = 8 + math.log(sascore + 1 - 9)

        return float(min(max(sascore, 1.0), 10.0)), fallback_count, len(nze)

    def get_property_names(self) -> List[str]:
        return ["sascore", "synthesis_difficulty"]

    def get_evaluator_info(self) -> Dict[str, Any]:
        return {
            "name": "Synthesis Accessibility Score (SAscore)",
            "description": "Estimates the synthetic difficulty of a molecule based on fragment occurrence and structural complexity.",
            "source": "Ertl & Schuffenhauer, J. Cheminform. 2009",
            "domain": "Synthesis Feasibility",
            "type": "rule-based",
            "range": "1 (easy) to 10 (hard)",
            "reference": "https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8",
        }

    def validate_parameters(self, parameters: Dict[str, Any]) -> bool:
        return True
