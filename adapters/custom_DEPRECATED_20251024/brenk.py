# pharmforge_core/evaluators/rules/brenk.py

from typing import Any, Dict, List
from rdkit import Chem
from rdkit.Chem import FilterCatalog

from pharmforge_core.core.molecule import UnifiedMolecule
from pharmforge_core.evaluators.evaluator_registry import EvaluatorComponent
from pharmforge_core.utils.decorators.unified_molecule_support import support_unified_molecule
from pharmforge_core.evaluators.base import register_evaluator


@register_evaluator
class BrenkAlertEvaluator(EvaluatorComponent):
    name = "brenk"

    @classmethod
    def get_component_name(cls) -> str:
        return cls.name

    @classmethod
    def get_component_type(cls) -> str:
        return "rule"

    def __init__(self):
        super().__init__()
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        self.catalog = FilterCatalog.FilterCatalog(params)
        self._initialize_alert_categories()

    def _initialize_alert_categories(self):
        self.categories = {
            "reactivity": ["michael", "aldehyde", "quinone", "peroxide"],
            "metabolic": ["nitro", "sulfoxide", "phosphor"],
            "toxicity": ["halo", "epoxide", "aziridine", "metal"],
            "aggregation": ["crown", "polycycle", "benzene"],
        }

    @support_unified_molecule(expected="mol", allow_batch=False, trace_output=True)
    def evaluate(self, mol: Chem.Mol, **kwargs) -> Dict[str, Any]:
        matches = self.catalog.GetMatches(mol)
        alerts = []
        for match in matches:
            desc = match.GetDescription()
            alerts.append({
                "name": desc,
                "category": self._categorize_alert(desc),
                "severity": self._assess_severity(desc),
            })

        categorized = self._categorize_alerts(alerts)
        risk = self._assess_risk(alerts)
        score = round(1.0 - (risk["score"] / 3.0), 3)

        return {
            "score": score,
            "pass": risk["level"] == "low",
            "tags": [f"brenk_{risk['level']}_risk"],
            "reason": risk["description"],
            "summary": {
                "brenk_alert_count": len(alerts),
                "brenk_risk_level": risk["level"],
                "class_label": risk["level"],
                "score": score,
            },
            "raw": {
                "alerts": alerts,
                "categories": categorized,
                "risk_assessment": risk,
            },
            "classifications": [risk["level"]],
            "evaluator_metadata": {
                "name": self.name,
                "type": "rule",
                "domain": "Toxicophore Alerts",
                "description": "Detects structural alerts using Brenk filter (ChemMedChem 2008).",
            },
        }

    def validate_parameters(self, parameters: Dict[str, Any]) -> bool:
        return True

    def get_property_names(self) -> List[str]:
        return ["brenk_alert_count", "brenk_risk_level", "score"]

    def _categorize_alert(self, description: str) -> str:
        description = description.lower()
        for category, keywords in self.categories.items():
            if any(keyword in description for keyword in keywords):
                return category
        return "other"

    def _assess_severity(self, description: str) -> int:
        description = description.lower()
        if any(k in description for k in ["peroxide", "epoxide", "aziridine", "thiol", "michael"]):
            return 3
        elif any(k in description for k in ["aldehyde", "nitro", "sulfoxide", "benzene", "halogen"]):
            return 2
        return 1

    def _categorize_alerts(self, alerts: List[Dict]) -> Dict:
        categories = {}
        for alert in alerts:
            category = alert["category"]
            if category not in categories:
                categories[category] = {"count": 0, "max_severity": 0, "alerts": []}
            categories[category]["count"] += 1
            categories[category]["max_severity"] = max(categories[category]["max_severity"], alert["severity"])
            categories[category]["alerts"].append(alert["name"])
        return categories

    def _assess_risk(self, alerts: List[Dict]) -> Dict:
        if not alerts:
            return {"level": "low", "score": 0, "description": "No structural alerts found"}

        max_severity = max(alert["severity"] for alert in alerts)
        num_high_risk = sum(1 for alert in alerts if alert["severity"] == 3)

        if num_high_risk > 1 or max_severity == 3:
            return {
                "level": "high",
                "score": 3,
                "description": "Multiple high-risk or severe structural alerts",
            }
        elif num_high_risk == 1 or len(alerts) > 2:
            return {
                "level": "medium",
                "score": 2,
                "description": "One high-risk or multiple lower-risk alerts",
            }
        else:
            return {"level": "low", "score": 1, "description": "Only minor structural alerts"}
