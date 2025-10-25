from typing import Any, Dict, List
from rdkit.Chem import Descriptors

from pharmforge_core.core.molecule import UnifiedMolecule
from pharmforge_core.evaluators.evaluator_registry import EvaluatorComponent
from pharmforge_core.utils.decorators.unified_molecule_support import support_unified_molecule
from pharmforge_core.evaluators.base import register_evaluator


@register_evaluator
class VeberEvaluator(EvaluatorComponent):
    name = "veber"

    @classmethod
    def get_component_name(cls) -> str:
        return cls.name

    @classmethod
    def get_component_type(cls) -> str:
        return "rule"

    @support_unified_molecule(expected="mol", allow_batch=False, trace_output=True)
    def evaluate(self, mol, **kwargs) -> Dict[str, Any]:
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
            "classifications": [bioavailability],
            "evaluator_metadata": {
                "name": self.name,
                "type": "rule",
                "domain": "Oral Bioavailability",
                "description": "Assesses compliance with Veber rules (TPSA, rotatable bonds, H-bonds) for oral drugs.",
                "thresholds": {
                    "rotatable_bonds_max": 10,
                    "TPSA_max": 140,
                    "total_H_bonds_max": 12,
                },
            },
        }

    def _calculate_normalized_score(self, rotatable: int, psa: float) -> float:
        rot_score = max(0.0, min(1.0, (10 - rotatable) / 10))
        psa_score = max(0.0, min(1.0, (140 - psa) / 140))
        return round((rot_score + psa_score) / 2, 3)

    def validate_parameters(self, parameters: Dict[str, Any]) -> bool:
        return True

    def get_property_names(self) -> List[str]:
        return [
            "veber_score",
            "oral_bioavailability",
            "veber_violations",
            "rotatable_bonds",
            "polar_surface_area",
            "total_h_bonds",
        ]
