"""
Datamol Adapter - Modern molecular manipulation toolkit
Provides SMILES standardization, molecular standardization, conformer generation,
clustering, and diversity analysis using the datamol library
"""
from typing import Any, Dict, List, Optional, Union
import logging

try:
    import datamol as dm
    DATAMOL_AVAILABLE = True
except ImportError:
    DATAMOL_AVAILABLE = False
    logging.warning("Datamol not available - install with: pip install datamol")

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit not available - install with: pip install rdkit")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DatamolAdapter(AdapterProtocol):
    """
    Adapter for datamol molecular manipulation toolkit
    Provides SMILES handling, standardization, conformer generation, clustering, and diversity analysis

    Supported operations:
    - standardize: Standardize molecules (default)
    - to_smiles: Convert to various SMILES formats
    - conformers: Generate 3D conformers
    - cluster: Cluster molecules by similarity
    - diversity: Calculate molecular diversity
    """

    def __init__(self):
        super().__init__(
            name="datamol",
            adapter_type="local",
            config={
                "timeout": 120,  # Some operations can be slow for large datasets
                "default_operation": "standardize",
                "supported_operations": ["standardize", "to_smiles", "conformers", "cluster", "diversity"]
            }
        )
        self.version = "1.0.0"

        if not DATAMOL_AVAILABLE:
            logger.error("Datamol is not installed! Install with: pip install datamol")

        if not RDKIT_AVAILABLE:
            logger.error("RDKit is not installed! Install with: pip install rdkit")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is valid for datamol operations

        Accepts:
        - Single SMILES string
        - Dictionary with 'smiles' key
        - Dictionary with 'smiles_list' key for batch operations

        Args:
            input_data: Input data to validate

        Returns:
            True if valid, False otherwise
        """
        if not DATAMOL_AVAILABLE or not RDKIT_AVAILABLE:
            return False

        # Single SMILES string
        if isinstance(input_data, str):
            if len(input_data) == 0:
                return False
            try:
                mol = dm.to_mol(input_data)
                return mol is not None
            except Exception:
                return False

        # Dictionary with 'smiles' key
        if isinstance(input_data, dict):
            if "smiles" in input_data:
                smiles = input_data["smiles"]
                if not isinstance(smiles, str) or len(smiles) == 0:
                    return False
                try:
                    mol = dm.to_mol(smiles)
                    return mol is not None
                except Exception:
                    return False

            # Dictionary with 'smiles_list' key for batch operations
            if "smiles_list" in input_data:
                smiles_list = input_data["smiles_list"]
                if not isinstance(smiles_list, list) or len(smiles_list) == 0:
                    return False
                if not all(isinstance(s, str) for s in smiles_list):
                    return False
                return True

        return False

    def _standardize_molecule(self, smiles: str, **kwargs) -> Optional[Dict[str, Any]]:
        """
        Standardize a molecule using datamol

        Args:
            smiles: Input SMILES string
            **kwargs: Additional standardization options

        Returns:
            Dictionary containing standardization results
        """
        if not DATAMOL_AVAILABLE:
            raise Exception("Datamol is not installed")

        try:
            mol = dm.to_mol(smiles)

            if mol is None:
                logger.warning(f"Datamol: Could not parse SMILES: {smiles}")
                return None

            # Standardize using datamol
            standardized_mol = dm.standardize_mol(mol)

            if standardized_mol is None:
                logger.warning(f"Datamol: Standardization failed for: {smiles}")
                return None

            results = {
                "original_smiles": smiles,
                "standardized_smiles": dm.to_smiles(standardized_mol),
                "structure_changed": dm.to_smiles(standardized_mol) != smiles,
                "num_atoms": standardized_mol.GetNumAtoms(),
                "num_heavy_atoms": standardized_mol.GetNumHeavyAtoms(),
                "molecular_weight": dm.descriptors.mw(standardized_mol),
                "canonical_smiles": dm.to_smiles(standardized_mol, canonical=True)
            }

            return results

        except Exception as e:
            logger.error(f"Datamol: Error standardizing {smiles}: {e}")
            return None

    def _to_smiles_operation(self, smiles: str, **kwargs) -> Optional[Dict[str, Any]]:
        """
        Convert molecule to various SMILES formats

        Args:
            smiles: Input SMILES string
            **kwargs: SMILES conversion options
                - canonical: bool (default: True)
                - isomeric: bool (default: True)
                - kekulize: bool (default: False)
                - ordered: bool (default: False)

        Returns:
            Dictionary containing various SMILES representations
        """
        if not DATAMOL_AVAILABLE:
            raise Exception("Datamol is not installed")

        try:
            mol = dm.to_mol(smiles)

            if mol is None:
                logger.warning(f"Datamol: Could not parse SMILES: {smiles}")
                return None

            # Get various SMILES representations
            results = {
                "original_smiles": smiles,
                "canonical_smiles": dm.to_smiles(mol, canonical=True),
                "isomeric_smiles": dm.to_smiles(mol, isomeric=True),
                "non_isomeric_smiles": dm.to_smiles(mol, isomeric=False),
                "kekulized_smiles": dm.to_smiles(mol, kekulize=True) if not kwargs.get("no_kekulize", False) else None,
                "inchi": dm.to_inchi(mol),
                "inchikey": dm.to_inchikey(mol),
                "selfies": dm.to_selfies(mol) if hasattr(dm, 'to_selfies') else None
            }

            return results

        except Exception as e:
            logger.error(f"Datamol: Error converting SMILES {smiles}: {e}")
            return None

    def _generate_conformers(self, smiles: str, **kwargs) -> Optional[Dict[str, Any]]:
        """
        Generate 3D conformers for a molecule

        Args:
            smiles: Input SMILES string
            **kwargs: Conformer generation options
                - n_confs: int - Number of conformers to generate (default: 10)
                - minimize_energy: bool - Minimize conformer energies (default: True)
                - align_conformers: bool - Align conformers to each other (default: True)

        Returns:
            Dictionary containing conformer information
        """
        if not DATAMOL_AVAILABLE:
            raise Exception("Datamol is not installed")

        try:
            mol = dm.to_mol(smiles)

            if mol is None:
                logger.warning(f"Datamol: Could not parse SMILES: {smiles}")
                return None

            # Generate conformers
            n_confs = kwargs.get("n_confs", 10)
            minimize_energy = kwargs.get("minimize_energy", True)
            align_conformers = kwargs.get("align_conformers", True)

            mol_with_confs = dm.conformers.generate(
                mol,
                n_confs=n_confs,
                minimize_energy=minimize_energy,
                align_conformers=align_conformers
            )

            if mol_with_confs is None or mol_with_confs.GetNumConformers() == 0:
                logger.warning(f"Datamol: Could not generate conformers for: {smiles}")
                return None

            # Get conformer energies if available
            conformer_energies = []
            if minimize_energy and hasattr(mol_with_confs, 'GetPropsAsDict'):
                props = mol_with_confs.GetPropsAsDict()
                for i in range(mol_with_confs.GetNumConformers()):
                    energy_key = f"Energy_{i}"
                    if energy_key in props:
                        conformer_energies.append(float(props[energy_key]))

            results = {
                "original_smiles": smiles,
                "num_conformers": mol_with_confs.GetNumConformers(),
                "n_confs_requested": n_confs,
                "minimized": minimize_energy,
                "aligned": align_conformers,
                "conformer_energies": conformer_energies if conformer_energies else None,
                "lowest_energy": min(conformer_energies) if conformer_energies else None
            }

            return results

        except Exception as e:
            logger.error(f"Datamol: Error generating conformers for {smiles}: {e}")
            return None

    def _cluster_molecules(self, smiles_list: List[str], **kwargs) -> Optional[Dict[str, Any]]:
        """
        Cluster molecules by structural similarity

        Args:
            smiles_list: List of SMILES strings
            **kwargs: Clustering options
                - n_clusters: int - Number of clusters (default: auto-detect)
                - method: str - Clustering method (default: "butina")
                - cutoff: float - Distance cutoff for Butina clustering (default: 0.35)

        Returns:
            Dictionary containing cluster assignments and statistics
        """
        if not DATAMOL_AVAILABLE:
            raise Exception("Datamol is not installed")

        try:
            # Convert SMILES to molecules
            mols = [dm.to_mol(s) for s in smiles_list]
            mols = [m for m in mols if m is not None]

            if len(mols) == 0:
                logger.warning("Datamol: No valid molecules to cluster")
                return None

            if len(mols) < 2:
                logger.warning("Datamol: Need at least 2 molecules for clustering")
                return {
                    "num_molecules": len(mols),
                    "num_clusters": 1,
                    "cluster_labels": [0] * len(mols),
                    "cluster_sizes": [len(mols)],
                    "method": "single"
                }

            # Generate fingerprints
            fps = [dm.to_fp(mol) for mol in mols]

            # Perform clustering using Butina algorithm (datamol's default)
            method = kwargs.get("method", "butina")
            cutoff = kwargs.get("cutoff", 0.35)

            # Datamol's cluster_mols function
            if hasattr(dm, 'cluster_mols'):
                clusters = dm.cluster_mols(mols, cutoff=cutoff)

                # Convert clusters to labels
                cluster_labels = [0] * len(mols)
                for cluster_id, cluster_indices in enumerate(clusters):
                    for idx in cluster_indices:
                        cluster_labels[idx] = cluster_id

                # Calculate cluster sizes
                cluster_sizes = [len(cluster) for cluster in clusters]

                results = {
                    "num_molecules": len(mols),
                    "num_clusters": len(clusters),
                    "cluster_labels": cluster_labels,
                    "cluster_sizes": cluster_sizes,
                    "method": method,
                    "cutoff": cutoff,
                    "smiles_list": smiles_list[:len(mols)]  # Only valid molecules
                }

                return results
            else:
                # Fallback to simple fingerprint-based clustering
                from sklearn.cluster import KMeans
                import numpy as np

                # Convert fingerprints to numpy array
                fp_array = np.array([list(fp) for fp in fps])

                n_clusters = kwargs.get("n_clusters", min(5, len(mols) // 2))
                n_clusters = max(2, min(n_clusters, len(mols)))

                kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
                cluster_labels = kmeans.fit_predict(fp_array).tolist()

                # Calculate cluster sizes
                cluster_sizes = [cluster_labels.count(i) for i in range(n_clusters)]

                results = {
                    "num_molecules": len(mols),
                    "num_clusters": n_clusters,
                    "cluster_labels": cluster_labels,
                    "cluster_sizes": cluster_sizes,
                    "method": "kmeans",
                    "smiles_list": smiles_list[:len(mols)]
                }

                return results

        except Exception as e:
            logger.error(f"Datamol: Error clustering molecules: {e}")
            return None

    def _calculate_diversity(self, smiles_list: List[str], **kwargs) -> Optional[Dict[str, Any]]:
        """
        Calculate molecular diversity metrics

        Args:
            smiles_list: List of SMILES strings
            **kwargs: Diversity calculation options

        Returns:
            Dictionary containing diversity metrics
        """
        if not DATAMOL_AVAILABLE:
            raise Exception("Datamol is not installed")

        try:
            # Convert SMILES to molecules
            mols = [dm.to_mol(s) for s in smiles_list]
            mols = [m for m in mols if m is not None]

            if len(mols) == 0:
                logger.warning("Datamol: No valid molecules to analyze")
                return None

            if len(mols) < 2:
                return {
                    "num_molecules": len(mols),
                    "mean_similarity": 1.0,
                    "mean_diversity": 0.0,
                    "min_similarity": 1.0,
                    "max_similarity": 1.0,
                    "pairwise_similarities": [1.0]
                }

            # Generate fingerprints
            fps = [dm.to_fp(mol) for mol in mols]

            # Calculate pairwise similarities
            from rdkit import DataStructs
            import numpy as np

            n_mols = len(fps)
            similarities = []

            for i in range(n_mols):
                for j in range(i + 1, n_mols):
                    sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    similarities.append(sim)

            similarities = np.array(similarities)
            diversities = 1.0 - similarities

            results = {
                "num_molecules": len(mols),
                "mean_similarity": float(np.mean(similarities)),
                "mean_diversity": float(np.mean(diversities)),
                "std_diversity": float(np.std(diversities)),
                "min_similarity": float(np.min(similarities)),
                "max_similarity": float(np.max(similarities)),
                "median_similarity": float(np.median(similarities)),
                "diversity_score": float(np.mean(diversities) / (np.std(diversities) + 1e-10)),
                "num_pairwise_comparisons": len(similarities),
                "smiles_list": smiles_list[:len(mols)]
            }

            return results

        except Exception as e:
            logger.error(f"Datamol: Error calculating diversity: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute datamol operations

        Args:
            input_data: SMILES string or dictionary with 'smiles' or 'smiles_list'
            **kwargs: Additional parameters
                - operation: str - Operation to perform (standardize, to_smiles, conformers, cluster, diversity)
                - Additional operation-specific parameters

        Returns:
            AdapterResult containing operation results
        """
        # Check if datamol and RDKit are available
        if not DATAMOL_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="Datamol is not installed. Install with: pip install datamol"
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
                error="Invalid input data. Expected SMILES string or dict with 'smiles' or 'smiles_list'"
            )

        # Get operation
        operation = kwargs.get("operation", self.config["default_operation"]).lower()

        if operation not in self.config["supported_operations"]:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Unsupported operation: {operation}. Supported: {self.config['supported_operations']}"
            )

        # Extract SMILES
        if isinstance(input_data, str):
            smiles = input_data
            smiles_list = None
        elif isinstance(input_data, dict):
            smiles = input_data.get("smiles")
            smiles_list = input_data.get("smiles_list")
        else:
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input format"
            )

        # Execute operation
        try:
            if operation == "standardize":
                if smiles is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="standardize operation requires a single SMILES string"
                    )
                result_data = self._standardize_molecule(smiles, **kwargs)

            elif operation == "to_smiles":
                if smiles is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="to_smiles operation requires a single SMILES string"
                    )
                result_data = self._to_smiles_operation(smiles, **kwargs)

            elif operation == "conformers":
                if smiles is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="conformers operation requires a single SMILES string"
                    )
                result_data = self._generate_conformers(smiles, **kwargs)

            elif operation == "cluster":
                if smiles_list is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="cluster operation requires a list of SMILES (use 'smiles_list' key)"
                    )
                result_data = self._cluster_molecules(smiles_list, **kwargs)

            elif operation == "diversity":
                if smiles_list is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="diversity operation requires a list of SMILES (use 'smiles_list' key)"
                    )
                result_data = self._calculate_diversity(smiles_list, **kwargs)

            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Operation {operation} not implemented"
                )

            if result_data is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Failed to execute {operation} operation",
                    metadata={
                        "source": "datamol",
                        "operation": operation
                    }
                )

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "datamol",
                    "adapter_version": self.version,
                    "computation_type": "local",
                    "operation": operation
                }
            )

        except Exception as e:
            logger.error(f"Datamol: Error executing {operation}: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error executing {operation}: {str(e)}",
                metadata={
                    "source": "datamol",
                    "operation": operation
                }
            )
