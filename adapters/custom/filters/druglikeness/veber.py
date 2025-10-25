# Veber Rules (Oral Bioavailability)
# Based on: Veber et al., J. Med. Chem. 2002
from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import Descriptors


class VeberEvaluator:
    """
    Assesses compliance with Veber rules for oral bioavailability.
    Rules: Rotatable bonds ≤ 10, TPSA ≤ 140, Total H-bonds ≤ 12

    Reference: Veber DF, et al. J. Med. Chem. 2002, 45, 12, 2615-2623
    """

    def evaluate(self, smiles: str, **kwargs) -> Dict[str, Any]:
        """
        Evaluate a molecule for Veber rules compliance.

        Args:
            smiles: SMILES string of molecule

        Returns:
            Dictionary with score, pass/fail, tags, reason, summary, raw data
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol or mol.GetNumAtoms() == 0:
            return {
                "score": 0.0,
                "pass": False,
                "tags": ["InvalidMolecule"],
                "reason": "Invalid SMILES string",
                "summary": {},
                "raw": {}
            }

        try:
            rotatable = Descriptors.NumRotatableBonds(mol)
            psa = Descriptors.TPSA(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            hb_total = hbd + hba

            rules = {
                "rotatable_bonds": rotatable <= 10,
                "polar_surface_area": psa <= 140,
                "hydrogen_bonds": hb_total <= 12,
            }

            violations = sum(not passed for passed in rules.values())
            pass_veber = violations == 0
            score = self._calculate_normalized_score(rotatable, psa)

            bioavailability = (
                "high" if violations == 0 else "moderate" if violations == 1 else "low"
            )
            tag = (
                "PassVeber" if violations == 0 else
                "PartialVeber" if violations == 1 else
                "FailVeber"
            )

            return {
                "score": score,
                "pass": pass_veber,
                "tags": [tag],
                "reason": f"{violations} Veber rule violation(s)",
                "summary": {
                    "veber_score": score,
                    "oral_bioavailability": bioavailability,
                    "veber_violations": violations,
                },
                "raw": {
                    "rotatable_bonds": rotatable,
                    "polar_surface_area": round(psa, 2),
                    "h_bond_donors": hbd,
                    "h_bond_acceptors": hba,
                    "total_h_bonds": hb_total,
                    "rules_check": rules,
                },
            }

        except Exception as e:
            return {
                "score": 0.0,
                "pass": False,
                "tags": ["EvaluationError"],
                "reason": f"Veber evaluation failed: {str(e)}",
                "summary": {},
                "raw": {}
            }

    def evaluate_batch(self, smiles_list: List[str], **kwargs) -> List[Dict[str, Any]]:
        """Evaluate multiple molecules"""
        return [self.evaluate(smi, **kwargs) for smi in smiles_list]

    def _calculate_normalized_score(self, rotatable: int, psa: float) -> float:
        rot_score = max(0.0, min(1.0, (10 - rotatable) / 10))
        psa_score = max(0.0, min(1.0, (140 - psa) / 140))
        return round((rot_score + psa_score) / 2, 3)

    def get_property_names(self) -> List[str]:
        return [
            "veber_score",
            "oral_bioavailability",
            "veber_violations",
            "rotatable_bonds",
            "polar_surface_area",
            "total_h_bonds",
        ]

    def get_evaluator_info(self) -> Dict[str, Any]:
        return {
            "name": "Veber Rules",
            "description": "Assesses compliance with Veber rules (TPSA, rotatable bonds, H-bonds) for oral drugs.",
            "source": "Veber et al., J. Med. Chem. 2002",
            "domain": "Oral Bioavailability",
            "type": "rule-based",
            "reference": "https://pubs.acs.org/doi/10.1021/jm020017n",
        }

    def validate_parameters(self, parameters: Dict[str, Any]) -> bool:
        return True
