"""
ChemPlot Adapter - Chemical space visualization using dimensionality reduction
Provides PCA, t-SNE, and UMAP projections for chemical space analysis
"""
from typing import Any, Dict, List, Optional
import logging
import numpy as np

try:
    from chemplot import Plotter
    CHEMPLOT_AVAILABLE = True
except ImportError:
    CHEMPLOT_AVAILABLE = False
    logging.warning("ChemPlot not available - install with: pip install chemplot")

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class ChemPlotAdapter(AdapterProtocol):
    """
    Adapter for chemical space visualization using ChemPlot
    Provides dimensionality reduction (PCA, t-SNE, UMAP) for SMILES lists
    Returns coordinates for frontend plotting, diversity metrics, and cluster assignments
    """

    def __init__(self):
        super().__init__(
            name="chemplot",
            adapter_type="local",
            config={
                "timeout": 120,  # Dimensionality reduction can be slow for large datasets
                "default_method": "pca",
                "supported_methods": ["pca", "tsne", "umap"]
            }
        )
        self.version = "1.0.0"

        if not CHEMPLOT_AVAILABLE:
            logger.error("ChemPlot is not installed! Install with: pip install chemplot")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a dictionary with a list of SMILES strings

        Args:
            input_data: Expected to be a dict with 'smiles_list' key

        Returns:
            True if valid, False otherwise
        """
        if not CHEMPLOT_AVAILABLE:
            return False

        # Check if input is a dictionary
        if not isinstance(input_data, dict):
            logger.warning("ChemPlot: Input must be a dictionary")
            return False

        # Check if smiles_list key exists
        if "smiles_list" not in input_data:
            logger.warning("ChemPlot: Input must contain 'smiles_list' key")
            return False

        smiles_list = input_data["smiles_list"]

        # Check if smiles_list is a list
        if not isinstance(smiles_list, list):
            logger.warning("ChemPlot: 'smiles_list' must be a list")
            return False

        # Check if list is not empty
        if len(smiles_list) == 0:
            logger.warning("ChemPlot: 'smiles_list' cannot be empty")
            return False

        # Check if all elements are strings
        if not all(isinstance(s, str) for s in smiles_list):
            logger.warning("ChemPlot: All elements in 'smiles_list' must be strings")
            return False

        # Minimum requirement: at least 2 molecules for meaningful visualization
        if len(smiles_list) < 2:
            logger.warning("ChemPlot: Need at least 2 molecules for visualization")
            return False

        return True

    def _calculate_diversity_metrics(self, coordinates: np.ndarray) -> Dict[str, float]:
        """
        Calculate diversity metrics from dimensionality-reduced coordinates

        Args:
            coordinates: Numpy array of shape (n_samples, n_dimensions)

        Returns:
            Dictionary containing diversity metrics
        """
        try:
            # Calculate pairwise distances
            from scipy.spatial.distance import pdist, squareform

            distances = pdist(coordinates, metric='euclidean')
            distance_matrix = squareform(distances)

            metrics = {
                "mean_distance": float(np.mean(distances)),
                "std_distance": float(np.std(distances)),
                "min_distance": float(np.min(distances)),
                "max_distance": float(np.max(distances)),
                "median_distance": float(np.median(distances)),
                "diversity_score": float(np.mean(distances) / (np.std(distances) + 1e-10))  # Normalized diversity
            }

            return metrics
        except Exception as e:
            logger.error(f"ChemPlot: Error calculating diversity metrics: {e}")
            return {
                "mean_distance": 0.0,
                "std_distance": 0.0,
                "min_distance": 0.0,
                "max_distance": 0.0,
                "median_distance": 0.0,
                "diversity_score": 0.0
            }

    def _perform_clustering(self, coordinates: np.ndarray, n_clusters: int = None) -> List[int]:
        """
        Perform clustering on dimensionality-reduced coordinates

        Args:
            coordinates: Numpy array of shape (n_samples, n_dimensions)
            n_clusters: Number of clusters (auto-detected if None)

        Returns:
            List of cluster IDs for each molecule
        """
        try:
            from sklearn.cluster import KMeans
            from sklearn.metrics import silhouette_score

            n_samples = coordinates.shape[0]

            # Auto-detect optimal number of clusters if not specified
            if n_clusters is None:
                # Try different numbers of clusters and use silhouette score
                max_k = min(10, n_samples // 2)  # Don't try more than 10 or n/2 clusters

                if max_k < 2:
                    # Not enough samples for clustering
                    return [0] * n_samples

                best_k = 2
                best_score = -1

                for k in range(2, max_k + 1):
                    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
                    labels = kmeans.fit_predict(coordinates)
                    score = silhouette_score(coordinates, labels)

                    if score > best_score:
                        best_score = score
                        best_k = k

                n_clusters = best_k
                logger.info(f"ChemPlot: Auto-detected {n_clusters} clusters (silhouette score: {best_score:.3f})")

            # Perform final clustering
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            cluster_labels = kmeans.fit_predict(coordinates)

            return cluster_labels.tolist()

        except Exception as e:
            logger.error(f"ChemPlot: Error performing clustering: {e}")
            # Return all molecules in one cluster as fallback
            return [0] * coordinates.shape[0]

    def _project_chemical_space(
        self,
        smiles_list: List[str],
        method: str = "pca",
        n_components: int = 2
    ) -> Optional[np.ndarray]:
        """
        Project chemical space using specified dimensionality reduction method

        Args:
            smiles_list: List of SMILES strings
            method: Dimensionality reduction method ('pca', 'tsne', or 'umap')
            n_components: Number of dimensions to reduce to (default: 2)

        Returns:
            Numpy array of coordinates or None on failure
        """
        if not CHEMPLOT_AVAILABLE:
            raise Exception("ChemPlot is not installed")

        try:
            # Create ChemPlot plotter
            cp = Plotter.from_smiles(smiles_list, target=None, sim_type="structural")

            # Perform dimensionality reduction
            if method.lower() == "pca":
                cp.pca()
            elif method.lower() == "tsne":
                cp.tsne()
            elif method.lower() == "umap":
                cp.umap()
            else:
                logger.warning(f"ChemPlot: Unknown method '{method}', falling back to PCA")
                cp.pca()

            # Extract coordinates from the plotter
            # ChemPlot stores the reduced coordinates in df_2_dim
            coordinates = cp.df_2_dim[["PC-1", "PC-2"]].values if method.lower() == "pca" else \
                         cp.df_2_dim[["t-SNE-1", "t-SNE-2"]].values if method.lower() == "tsne" else \
                         cp.df_2_dim[["UMAP-1", "UMAP-2"]].values

            return coordinates

        except Exception as e:
            logger.error(f"ChemPlot: Error during projection: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute chemical space visualization

        Args:
            input_data: Dictionary with 'smiles_list' key containing list of SMILES
            **kwargs: Additional parameters:
                - method: Dimensionality reduction method ('pca', 'tsne', 'umap')
                - n_components: Number of dimensions (default: 2)
                - n_clusters: Number of clusters (auto-detected if None)

        Returns:
            AdapterResult containing coordinates, cluster assignments, and diversity metrics
        """
        # Check if ChemPlot is available
        if not CHEMPLOT_AVAILABLE:
            return AdapterResult(
                success=False,
                data=None,
                error="ChemPlot is not installed. Install with: pip install chemplot"
            )

        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input format. Expected: {'smiles_list': ['CCO', 'CC(=O)O', ...]}"
            )

        smiles_list = input_data["smiles_list"]
        method = kwargs.get("method", self.config["default_method"]).lower()
        n_components = kwargs.get("n_components", 2)
        n_clusters = kwargs.get("n_clusters", None)

        # Validate method
        if method not in self.config["supported_methods"]:
            logger.warning(f"ChemPlot: Method '{method}' not supported, using PCA")
            method = "pca"

        # Perform dimensionality reduction
        coordinates = self._project_chemical_space(smiles_list, method, n_components)

        if coordinates is None:
            return AdapterResult(
                success=False,
                data=None,
                error=f"Failed to perform {method.upper()} projection",
                metadata={
                    "source": "chemplot",
                    "method": method,
                    "n_molecules": len(smiles_list)
                }
            )

        # Calculate diversity metrics
        diversity_metrics = self._calculate_diversity_metrics(coordinates)

        # Perform clustering
        cluster_labels = self._perform_clustering(coordinates, n_clusters)

        # Prepare result data
        result_data = {
            "coordinates": coordinates.tolist(),  # Convert numpy array to list for JSON serialization
            "cluster_labels": cluster_labels,
            "diversity_metrics": diversity_metrics,
            "method": method.upper(),
            "n_molecules": len(smiles_list),
            "n_dimensions": n_components,
            "n_clusters": len(set(cluster_labels)),
            "smiles_list": smiles_list  # Include original SMILES for reference
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "chemplot",
                "adapter_version": self.version,
                "computation_type": "local",
                "method": method,
                "n_molecules": len(smiles_list),
                "n_clusters": len(set(cluster_labels))
            }
        )
