from typing import Any, Dict, List

from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from pharmforge_core.core.molecule import UnifiedMolecule
from pharmforge_core.evaluators.evaluator_registry import EvaluatorComponent
from pharmforge_core.utils.decorators.unified_molecule_support import support_unified_molecule
from pharmforge_core.evaluators.base import register_evaluator


@register_evaluator
class PAINSEvaluator(EvaluatorComponent):
    name = "pains"

    @classmethod
    def get_component_name(cls) -> str:
        return cls.name

    @classmethod
    def get_component_type(cls) -> str:
        return "rule"

    def __init__(self):
        super().__init__()
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        self.catalog = FilterCatalog.FilterCatalog(params)

    @support_unified_molecule(expected="mol", allow_batch=False, trace_output=True)
    def evaluate(self, mol, **kwargs) -> Dict[str, Any]:
        matches = self.catalog.GetMatches(mol)
        alert_count = len(matches)
        score = 1.0 if alert_count == 0 else max(0.0, 1.0 - (alert_count * 0.2))

        passed = alert_count == 0
        tags = ["pains_pass"] if passed else ["pains_fail"] + [f"pains_alert_{i+1}" for i in range(alert_count)]
        label = "pass" if passed else "fail"

        alert_names = [m.GetDescription() for m in matches]

        return {
            "score": round(score, 3),
            "pass": passed,
            "tags": tags,
            "reason": "No PAINS alerts detected" if passed else f"{alert_count} PAINS alert(s) found",
            "summary": {
                "pains_alerts": alert_count,
                "class_label": label,
            },
            "raw": {
                "alerts_detected": alert_count,
                "alert_descriptions": alert_names,
                "pass_pains": passed,
            },
            "classifications": [label],
            "evaluator_metadata": {
                "name": self.name,
                "type": "rule",
                "domain": "Toxicity / Screening Artifact",
                "description": "Detects Pan-Assay Interference Compounds (PAINS) using RDKit catalog.",
                "thresholds": {
                    "max_allowed_alerts": 0,
                    "penalty_per_alert": 0.2,
                },
            },
        }

    def validate_parameters(self, parameters: Dict[str, Any]) -> bool:
        return True

    def get_property_names(self) -> List[str]:
        return ["pains_alerts", "pass_pains"]
