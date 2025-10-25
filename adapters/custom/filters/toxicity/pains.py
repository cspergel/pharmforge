# PAINS (Pan-Assay Interference Compounds) Filter
# Adapted from RDKit FilterCatalog: Baell & Holloway, J. Med. Chem. 2010
from typing import Any, Dict, List

from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams


class PAINSEvaluator:
    """
    Detects Pan-Assay Interference Compounds (PAINS) using RDKit filter catalog.
    PAINS are compounds that frequently produce false positives in biochemical assays.

    Reference: Baell JB, Holloway GA. J. Med. Chem. 2010, 53, 7, 2719-2740
    """

    def __init__(self):
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        self.catalog = FilterCatalog.FilterCatalog(params)

    def evaluate(self, smiles: str, **kwargs) -> Dict[str, Any]:
        """
        Evaluate a molecule for PAINS alerts.

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
            }

        except Exception as e:
            return {
                "score": 0.0,
                "pass": False,
                "tags": ["EvaluationError"],
                "reason": f"PAINS evaluation failed: {str(e)}",
                "summary": {},
                "raw": {}
            }

    def evaluate_batch(self, smiles_list: List[str], **kwargs) -> List[Dict[str, Any]]:
        """Evaluate multiple molecules"""
        return [self.evaluate(smi, **kwargs) for smi in smiles_list]

    def get_property_names(self) -> List[str]:
        return ["pains_alerts", "pass_pains"]

    def get_evaluator_info(self) -> Dict[str, Any]:
        return {
            "name": "PAINS Filter",
            "description": "Detects Pan-Assay Interference Compounds (PAINS) using RDKit catalog.",
            "source": "Baell & Holloway, J. Med. Chem. 2010",
            "domain": "Toxicity / Screening Artifact",
            "type": "rule-based",
            "reference": "https://pubs.acs.org/doi/10.1021/jm901137j",
        }

    def validate_parameters(self, parameters: Dict[str, Any]) -> bool:
        return True
