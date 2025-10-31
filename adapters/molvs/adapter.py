"""
MolVS Adapter - Molecule validation and standardization
Provides molecule standardization, validation, and structure cleanup
"""
from typing import Any, Dict, Optional, List
import logging

try:
    from molvs import Standardizer, validate_smiles
    from molvs.fragment import LargestFragmentChooser
    from molvs.charge import Uncharger
    from molvs.tautomer import TautomerCanonicalizer
    MOLVS_AVAILABLE = True
except ImportError:
    MOLVS_AVAILABLE = False
    logging.warning("MolVS not available - install with: pip install molvs")

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - install with: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class MolVSAdapter(AdapterProtocol):
    """
    Adapter for MolVS molecule validation and standardization
    Standardizes molecules, validates structures, and identifies chemical issues
    """

    def __init__(self):
        super().__init__(
            name="molvs",
            adapter_type="local",
            config={
                "timeout": 10,  # Local calculations should be fast
                "include_validation": True,
                "include_tautomers": True
            }
        )
        self.version = "1.0.0"

        if not MOLVS_AVAILABLE:
            logger.error("MolVS is not installed! Install with: pip install molvs")

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

        # Initialize MolVS components if available
        if MOLVS_AVAILABLE:
            self.standardizer = Standardizer()
            self.fragment_chooser = LargestFragmentChooser()
            self.uncharger = Uncharger()
            self.tautomer_canonicalizer = TautomerCanonicalizer()

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid SMILES string

        Args:
            input_data: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not MOLVS_AVAILABLE or not RDKIT_AVAILABLE:
            return False

        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False

        # Try to parse with RDKit
        try:
            mol = Chem.MolFromSmiles(input_data)
            return mol is not None
        except Exception:
            return False

    def _standardize_molecule(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Standardize a molecule using MolVS

        Args:
            smiles: Input SMILES string

        Returns:
            Dictionary containing standardization results
        """
        if not MOLVS_AVAILABLE or not RDKIT_AVAILABLE:
            raise Exception("MolVS or RDKit is not installed")

        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                logger.warning(f"MolVS: Could not parse SMILES: {smiles}")
                return None

            results = {
                "original_smiles": smiles,
                "original_mol_valid": True,
            }

            # Step 1: Get largest fragment (removes salts, solvents)
            try:
                mol_fragment = self.fragment_chooser.choose(mol)
                fragment_smiles = Chem.MolToSmiles(mol_fragment)
                results["fragment_smiles"] = fragment_smiles
                results["fragments_removed"] = fragment_smiles != Chem.MolToSmiles(mol)
            except Exception as e:
                logger.warning(f"MolVS: Fragment selection failed: {e}")
                mol_fragment = mol
                results["fragment_smiles"] = smiles
                results["fragments_removed"] = False

            # Step 2: Neutralize charges
            try:
                mol_uncharged = self.uncharger.uncharge(mol_fragment)
                uncharged_smiles = Chem.MolToSmiles(mol_uncharged)
                results["uncharged_smiles"] = uncharged_smiles
                results["charges_neutralized"] = uncharged_smiles != results["fragment_smiles"]
            except Exception as e:
                logger.warning(f"MolVS: Uncharging failed: {e}")
                mol_uncharged = mol_fragment
                results["uncharged_smiles"] = results["fragment_smiles"]
                results["charges_neutralized"] = False

            # Step 3: Full standardization (functional groups, bonds, etc.)
            try:
                mol_standardized = self.standardizer.standardize(mol_uncharged)
                standardized_smiles = Chem.MolToSmiles(mol_standardized)
                results["standardized_smiles"] = standardized_smiles
            except Exception as e:
                logger.warning(f"MolVS: Standardization failed: {e}")
                results["standardized_smiles"] = results["uncharged_smiles"]
                mol_standardized = mol_uncharged

            # Step 4: Canonical tautomer
            if self.config.get("include_tautomers", True):
                try:
                    mol_tautomer = self.tautomer_canonicalizer.canonicalize(mol_standardized)
                    tautomer_smiles = Chem.MolToSmiles(mol_tautomer)
                    results["canonical_tautomer_smiles"] = tautomer_smiles
                    results["tautomer_changed"] = tautomer_smiles != results["standardized_smiles"]
                except Exception as e:
                    logger.warning(f"MolVS: Tautomer canonicalization failed: {e}")
                    results["canonical_tautomer_smiles"] = results["standardized_smiles"]
                    results["tautomer_changed"] = False
                    mol_tautomer = mol_standardized
            else:
                results["canonical_tautomer_smiles"] = results["standardized_smiles"]
                results["tautomer_changed"] = False
                mol_tautomer = mol_standardized

            # Final canonical SMILES
            results["final_smiles"] = Chem.MolToSmiles(mol_tautomer)
            results["structure_changed"] = results["final_smiles"] != smiles

            return results

        except Exception as e:
            logger.error(f"MolVS: Error standardizing {smiles}: {e}")
            return None

    def _validate_molecule(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Validate molecule structure and identify issues

        Args:
            smiles: SMILES string to validate

        Returns:
            Dictionary containing validation results
        """
        if not MOLVS_AVAILABLE or not RDKIT_AVAILABLE:
            raise Exception("MolVS or RDKit is not installed")

        try:
            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                return {
                    "valid": False,
                    "issues": ["Could not parse SMILES string"],
                    "num_issues": 1
                }

            # Use MolVS validate_smiles function
            validation_issues = []

            try:
                # validate_smiles returns a list of validation messages
                issues = validate_smiles(smiles)

                if issues:
                    for issue in issues:
                        validation_issues.append(str(issue))

            except Exception as e:
                logger.warning(f"MolVS: Validation check failed: {e}")

            # Additional RDKit-based checks
            num_atoms = mol.GetNumAtoms()
            num_heavy_atoms = mol.GetNumHeavyAtoms()

            # Check for common issues
            if num_atoms == 0:
                validation_issues.append("No atoms in molecule")

            if num_heavy_atoms == 0:
                validation_issues.append("No heavy atoms in molecule")

            # Check for disconnected fragments
            fragments = Chem.GetMolFrags(mol, asMols=True)
            if len(fragments) > 1:
                validation_issues.append(f"Multiple fragments detected ({len(fragments)} fragments)")

            # Check for radical electrons
            for atom in mol.GetAtoms():
                if atom.GetNumRadicalElectrons() > 0:
                    validation_issues.append(f"Radical electron detected on atom {atom.GetSymbol()}")
                    break

            return {
                "valid": len(validation_issues) == 0,
                "issues": validation_issues,
                "num_issues": len(validation_issues),
                "num_atoms": num_atoms,
                "num_heavy_atoms": num_heavy_atoms,
                "num_fragments": len(fragments)
            }

        except Exception as e:
            logger.error(f"MolVS: Error validating {smiles}: {e}")
            return {
                "valid": False,
                "issues": [f"Validation error: {str(e)}"],
                "num_issues": 1
            }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute MolVS standardization and validation

        Args:
            input_data: SMILES string
            **kwargs: Additional parameters
                - include_validation: bool - Include validation results (default: True)

        Returns:
            AdapterResult containing standardized SMILES and validation info
        """
        # Check if MolVS and RDKit are available
        if not MOLVS_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="MolVS is not installed. Install with: pip install molvs"
            )

        if not RDKIT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="RDKit is not installed. Install with: pip install rdkit"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid SMILES string or could not parse with RDKit"
            )

        smiles = input_data

        # Standardize molecule
        standardization = self._standardize_molecule(smiles)

        if standardization is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to standardize molecule with MolVS",
                metadata={"source": "molvs", "smiles": smiles}
            )

        result_data = {
            "standardization": standardization
        }

        # Optionally include validation
        include_validation = kwargs.get("include_validation", self.config.get("include_validation", True))

        if include_validation:
            # Validate the original SMILES
            validation = self._validate_molecule(smiles)
            result_data["validation"] = validation

            # Also validate the standardized SMILES
            final_smiles = standardization.get("final_smiles", smiles)
            if final_smiles != smiles:
                standardized_validation = self._validate_molecule(final_smiles)
                result_data["standardized_validation"] = standardized_validation

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "molvs",
                "smiles": smiles,
                "adapter_version": self.version,
                "computation_type": "local",
                "final_smiles": standardization.get("final_smiles", smiles),
                "structure_changed": standardization.get("structure_changed", False)
            }
        )
