"""
MDAnalysis Adapter for PharmForge

Provides molecular dynamics trajectory analysis using the MDAnalysis toolkit.

MDAnalysis is a powerful library for analyzing molecular dynamics simulations,
supporting multiple trajectory formats (DCD, XTC, TRR, etc.) and topology files.

Features:
- RMSD (Root Mean Square Deviation) calculation
- RMSF (Root Mean Square Fluctuation) calculation
- Radius of gyration tracking
- Hydrogen bond analysis
- Distance analysis
- Dihedral angle analysis
- Support for multiple trajectory formats

Reference: https://www.mdanalysis.org/
Papers:
- Michaud-Agrawal et al., J. Comput. Chem. 2011, 32, 2319-2327
- Gowers et al., Proc. of the 15th Python in Science Conf. 2016, 98-105
"""

import hashlib
import logging
import json
from typing import Dict, Any, Optional, List, Tuple
import asyncio
from pathlib import Path

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)

# Try to import MDAnalysis
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align, distances
    try:
        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
        HBA_AVAILABLE = True
    except ImportError:
        # Older versions use different import
        try:
            from MDAnalysis.analysis.hbonds.hbond_analysis import HydrogenBondAnalysis
            HBA_AVAILABLE = True
        except ImportError:
            HBA_AVAILABLE = False
            logger.warning("HydrogenBondAnalysis not available in this MDAnalysis version")

    MDANALYSIS_AVAILABLE = True
    logger.info(f"MDAnalysis version: {mda.__version__}")
except ImportError:
    MDANALYSIS_AVAILABLE = False
    HBA_AVAILABLE = False
    logger.warning("MDAnalysis not available - install with: pip install MDAnalysis")


class MDAnalysisAdapter(AdapterProtocol):
    """
    Adapter for MDAnalysis trajectory analysis.

    Analyzes molecular dynamics simulation trajectories to compute structural
    properties and dynamics over time. Useful for:
    - Validating MD simulation results
    - Analyzing protein stability and flexibility
    - Tracking conformational changes
    - Understanding molecular behavior

    Features:
    - Multiple trajectory format support (DCD, XTC, TRR, NetCDF, etc.)
    - RMSD/RMSF calculations
    - Radius of gyration tracking
    - Hydrogen bond analysis
    - Distance and angle measurements
    - Comprehensive statistical analysis
    """

    def __init__(
        self,
        name: str = "mdanalysis",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize MDAnalysis adapter.

        Args:
            name: Adapter name (default: "mdanalysis")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - selection: Atom selection string (default: "protein and name CA")
                   - compute_rmsd: Calculate RMSD (default: True)
                   - compute_rmsf: Calculate RMSF (default: True)
                   - compute_rg: Calculate radius of gyration (default: True)
                   - compute_hbonds: Analyze hydrogen bonds (default: False)
                   - align_trajectory: Align trajectory before analysis (default: True)
                   - reference_frame: Frame to use as reference (default: 0)
                   - step: Step size for trajectory iteration (default: 1)
        """
        default_config = {
            "selection": "protein and name CA",  # Default to C-alpha atoms
            "compute_rmsd": True,
            "compute_rmsf": True,
            "compute_rg": True,
            "compute_hbonds": False,  # More expensive
            "align_trajectory": True,
            "reference_frame": 0,  # First frame as reference
            "step": 1,  # Analyze every frame
            "hbond_distance": 3.0,  # Angstroms
            "hbond_angle": 150.0,  # Degrees
        }

        # Merge with provided config
        merged_config = {**default_config, **(config or {})}

        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        if not MDANALYSIS_AVAILABLE:
            logger.error("MDAnalysis is not installed! Install with: pip install MDAnalysis")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input contains required trajectory and topology paths.

        Args:
            input_data: Dictionary with 'topology' and 'trajectory' file paths,
                       or just 'topology' for single-frame analysis

        Returns:
            True if valid, False otherwise
        """
        if not MDANALYSIS_AVAILABLE:
            logger.error("MDAnalysis is not installed")
            return False

        if not isinstance(input_data, dict):
            logger.error("Input must be a dictionary with 'topology' and optional 'trajectory' keys")
            return False

        if 'topology' not in input_data:
            logger.error("Input must contain 'topology' key with path to topology file")
            return False

        # Check if files exist
        topology_path = Path(input_data['topology'])
        if not topology_path.exists():
            logger.error(f"Topology file does not exist: {topology_path}")
            return False

        # Trajectory is optional (can analyze single structure)
        if 'trajectory' in input_data:
            trajectory_path = Path(input_data['trajectory'])
            if not trajectory_path.exists():
                logger.error(f"Trajectory file does not exist: {trajectory_path}")
                return False

        return True

    def _load_universe(
        self,
        topology: str,
        trajectory: Optional[str] = None
    ) -> Any:
        """
        Load MDAnalysis Universe from topology and trajectory files.

        Args:
            topology: Path to topology file (PDB, PSF, GRO, etc.)
            trajectory: Optional path to trajectory file (DCD, XTC, TRR, etc.)

        Returns:
            MDAnalysis Universe object
        """
        if not MDANALYSIS_AVAILABLE:
            raise ImportError("MDAnalysis is not installed")

        try:
            if trajectory:
                universe = mda.Universe(topology, trajectory)
                logger.info(f"Loaded Universe: {universe.atoms.n_atoms} atoms, "
                          f"{universe.trajectory.n_frames} frames")
            else:
                universe = mda.Universe(topology)
                logger.info(f"Loaded Universe: {universe.atoms.n_atoms} atoms, single frame")

            return universe

        except Exception as e:
            logger.error(f"Failed to load Universe: {e}")
            raise

    def _calculate_rmsd(
        self,
        universe: Any,
        selection: str,
        reference_frame: int = 0
    ) -> Dict[str, Any]:
        """
        Calculate RMSD over trajectory.

        Args:
            universe: MDAnalysis Universe
            selection: Atom selection string
            reference_frame: Frame index to use as reference

        Returns:
            Dictionary with RMSD results
        """
        try:
            # Select atoms
            select = universe.select_atoms(selection)

            if select.n_atoms == 0:
                logger.warning(f"No atoms selected with: {selection}")
                return {"error": "No atoms selected"}

            # Align if configured
            if self.config.get("align_trajectory", True):
                align.AlignTraj(
                    universe,
                    universe,
                    select=selection,
                    in_memory=False
                ).run()
                logger.info("Aligned trajectory to reference")

            # Calculate RMSD
            rmsd_analysis = rms.RMSD(
                select,
                select,
                ref_frame=reference_frame
            )
            rmsd_analysis.run(step=self.config.get("step", 1))

            # Extract results
            rmsd_values = rmsd_analysis.results.rmsd[:, 2]  # Column 2 is RMSD
            time_values = rmsd_analysis.results.rmsd[:, 1]  # Column 1 is time

            import numpy as np

            result = {
                "n_frames": len(rmsd_values),
                "rmsd_mean": float(np.mean(rmsd_values)),
                "rmsd_std": float(np.std(rmsd_values)),
                "rmsd_min": float(np.min(rmsd_values)),
                "rmsd_max": float(np.max(rmsd_values)),
                "rmsd_values": [float(v) for v in rmsd_values],
                "time_values": [float(t) for t in time_values],
                "selection": selection,
                "n_atoms": select.n_atoms
            }

            logger.info(f"RMSD: {result['rmsd_mean']:.3f} ± {result['rmsd_std']:.3f} Å "
                       f"(min: {result['rmsd_min']:.3f}, max: {result['rmsd_max']:.3f})")

            return result

        except Exception as e:
            logger.error(f"RMSD calculation failed: {e}")
            return {"error": str(e)}

    def _calculate_rmsf(
        self,
        universe: Any,
        selection: str
    ) -> Dict[str, Any]:
        """
        Calculate RMSF (per-residue or per-atom fluctuations).

        Args:
            universe: MDAnalysis Universe
            selection: Atom selection string

        Returns:
            Dictionary with RMSF results
        """
        try:
            select = universe.select_atoms(selection)

            if select.n_atoms == 0:
                logger.warning(f"No atoms selected with: {selection}")
                return {"error": "No atoms selected"}

            # Align trajectory first
            if self.config.get("align_trajectory", True):
                align.AlignTraj(
                    universe,
                    universe,
                    select=selection,
                    in_memory=False
                ).run()

            # Calculate RMSF
            rmsf_analysis = rms.RMSF(select)
            rmsf_analysis.run(step=self.config.get("step", 1))

            import numpy as np

            rmsf_values = rmsf_analysis.results.rmsf

            result = {
                "n_atoms": len(rmsf_values),
                "rmsf_mean": float(np.mean(rmsf_values)),
                "rmsf_std": float(np.std(rmsf_values)),
                "rmsf_min": float(np.min(rmsf_values)),
                "rmsf_max": float(np.max(rmsf_values)),
                "rmsf_values": [float(v) for v in rmsf_values],
                "selection": selection
            }

            logger.info(f"RMSF: {result['rmsf_mean']:.3f} ± {result['rmsf_std']:.3f} Å "
                       f"(min: {result['rmsf_min']:.3f}, max: {result['rmsf_max']:.3f})")

            return result

        except Exception as e:
            logger.error(f"RMSF calculation failed: {e}")
            return {"error": str(e)}

    def _calculate_radius_of_gyration(
        self,
        universe: Any,
        selection: str
    ) -> Dict[str, Any]:
        """
        Calculate radius of gyration over trajectory.

        Args:
            universe: MDAnalysis Universe
            selection: Atom selection string

        Returns:
            Dictionary with Rg results
        """
        try:
            select = universe.select_atoms(selection)

            if select.n_atoms == 0:
                logger.warning(f"No atoms selected with: {selection}")
                return {"error": "No atoms selected"}

            rg_values = []
            time_values = []

            # Iterate through trajectory
            for ts in universe.trajectory[::self.config.get("step", 1)]:
                rg = select.radius_of_gyration()
                rg_values.append(rg)
                time_values.append(ts.time)

            import numpy as np
            rg_array = np.array(rg_values)

            result = {
                "n_frames": len(rg_values),
                "rg_mean": float(np.mean(rg_array)),
                "rg_std": float(np.std(rg_array)),
                "rg_min": float(np.min(rg_array)),
                "rg_max": float(np.max(rg_array)),
                "rg_values": [float(v) for v in rg_values],
                "time_values": [float(t) for t in time_values],
                "selection": selection,
                "n_atoms": select.n_atoms
            }

            logger.info(f"Radius of gyration: {result['rg_mean']:.3f} ± {result['rg_std']:.3f} Å "
                       f"(min: {result['rg_min']:.3f}, max: {result['rg_max']:.3f})")

            return result

        except Exception as e:
            logger.error(f"Radius of gyration calculation failed: {e}")
            return {"error": str(e)}

    def _analyze_hydrogen_bonds(
        self,
        universe: Any
    ) -> Dict[str, Any]:
        """
        Analyze hydrogen bonds in trajectory.

        Args:
            universe: MDAnalysis Universe

        Returns:
            Dictionary with hydrogen bond analysis results
        """
        if not HBA_AVAILABLE:
            return {"error": "HydrogenBondAnalysis not available in this MDAnalysis version"}

        try:
            # Create hydrogen bond analysis
            hbonds = HydrogenBondAnalysis(
                universe=universe,
                donors_sel="protein",
                hydrogens_sel="protein",
                acceptors_sel="protein",
                d_a_cutoff=self.config.get("hbond_distance", 3.0),
                d_h_a_angle_cutoff=self.config.get("hbond_angle", 150.0)
            )

            # Run analysis
            hbonds.run(step=self.config.get("step", 1))

            # Extract results
            n_hbonds = len(hbonds.results.hbonds)

            result = {
                "n_hbonds_total": n_hbonds,
                "n_frames": universe.trajectory.n_frames,
                "distance_cutoff": self.config.get("hbond_distance", 3.0),
                "angle_cutoff": self.config.get("hbond_angle", 150.0)
            }

            if n_hbonds > 0:
                import numpy as np
                # Count unique hydrogen bonds
                unique_hbonds = len(np.unique(hbonds.results.hbonds[:, [0, 1, 2]], axis=0))
                result["n_unique_hbonds"] = unique_hbonds
                result["avg_hbonds_per_frame"] = n_hbonds / universe.trajectory.n_frames

            logger.info(f"Hydrogen bonds: {n_hbonds} total, {result.get('n_unique_hbonds', 0)} unique")

            return result

        except Exception as e:
            logger.error(f"Hydrogen bond analysis failed: {e}")
            return {"error": str(e)}

    def _get_system_info(self, universe: Any) -> Dict[str, Any]:
        """
        Extract basic system information from Universe.

        Args:
            universe: MDAnalysis Universe

        Returns:
            Dictionary with system information
        """
        try:
            info = {
                "n_atoms": universe.atoms.n_atoms,
                "n_residues": universe.atoms.n_residues,
                "n_segments": universe.atoms.n_segments,
            }

            # Trajectory information
            if hasattr(universe, 'trajectory'):
                info["n_frames"] = universe.trajectory.n_frames
                info["total_time"] = float(universe.trajectory.totaltime)
                info["dt"] = float(universe.trajectory.dt)

            # Residue composition
            if hasattr(universe.atoms, 'resnames'):
                resname_counts = {}
                for resname in universe.atoms.resnames:
                    resname_counts[resname] = resname_counts.get(resname, 0) + 1
                info["residue_composition"] = resname_counts

            return info

        except Exception as e:
            logger.warning(f"Could not extract all system info: {e}")
            return {"n_atoms": universe.atoms.n_atoms}

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute MDAnalysis trajectory analysis.

        Args:
            input_data: Dictionary containing:
                       - 'topology': Path to topology file (required)
                       - 'trajectory': Path to trajectory file (optional)
            **kwargs: Additional parameters:
                     - selection: Override default atom selection
                     - compute_rmsd: Override RMSD calculation
                     - compute_rmsf: Override RMSF calculation
                     - compute_rg: Override Rg calculation
                     - compute_hbonds: Override H-bond analysis

        Returns:
            AdapterResult containing trajectory analysis results
        """
        try:
            # Check if MDAnalysis is available
            if not MDANALYSIS_AVAILABLE:
                return AdapterResult(
                    success=False,
                    data={},
                    error="MDAnalysis is not installed. Install with: pip install MDAnalysis"
                )

            # Validate input
            if not self.validate_input(input_data):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid input. Expected dict with 'topology' and optional 'trajectory' keys"
                )

            # Extract file paths
            topology = input_data['topology']
            trajectory = input_data.get('trajectory')

            # Get configuration (allow runtime override)
            selection = kwargs.get('selection', self.config.get('selection', 'protein and name CA'))
            compute_rmsd = kwargs.get('compute_rmsd', self.config.get('compute_rmsd', True))
            compute_rmsf = kwargs.get('compute_rmsf', self.config.get('compute_rmsf', True))
            compute_rg = kwargs.get('compute_rg', self.config.get('compute_rg', True))
            compute_hbonds = kwargs.get('compute_hbonds', self.config.get('compute_hbonds', False))

            logger.info(f"Analyzing trajectory: {topology}")
            if trajectory:
                logger.info(f"  Trajectory: {trajectory}")
            logger.info(f"  Selection: {selection}")

            # Run analysis in thread pool (MDAnalysis operations are synchronous)
            loop = asyncio.get_event_loop()
            result_data = await loop.run_in_executor(
                None,
                self._run_analysis,
                topology,
                trajectory,
                selection,
                compute_rmsd,
                compute_rmsf,
                compute_rg,
                compute_hbonds
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "version": self.version,
                    "topology": topology,
                    "trajectory": trajectory,
                    "selection": selection
                }
            )

        except ImportError as e:
            logger.error(f"Required library not available: {e}")
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name,
                    "installation_help": "pip install MDAnalysis"
                }
            )

        except Exception as e:
            logger.error(f"MDAnalysis trajectory analysis failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    def _run_analysis(
        self,
        topology: str,
        trajectory: Optional[str],
        selection: str,
        compute_rmsd: bool,
        compute_rmsf: bool,
        compute_rg: bool,
        compute_hbonds: bool
    ) -> Dict[str, Any]:
        """
        Run the full analysis workflow (synchronous method for thread pool).

        Args:
            topology: Path to topology file
            trajectory: Optional path to trajectory file
            selection: Atom selection string
            compute_rmsd: Whether to compute RMSD
            compute_rmsf: Whether to compute RMSF
            compute_rg: Whether to compute radius of gyration
            compute_hbonds: Whether to analyze hydrogen bonds

        Returns:
            Dictionary with analysis results
        """
        # Load Universe
        logger.info("Loading Universe...")
        universe = self._load_universe(topology, trajectory)

        # Get system information
        system_info = self._get_system_info(universe)

        # Initialize results
        results = {
            "system": system_info,
            "selection": selection
        }

        # Only run trajectory analyses if we have a trajectory
        has_trajectory = trajectory is not None and universe.trajectory.n_frames > 1

        # Calculate RMSD
        if compute_rmsd and has_trajectory:
            logger.info("Calculating RMSD...")
            results["rmsd"] = self._calculate_rmsd(
                universe,
                selection,
                self.config.get("reference_frame", 0)
            )

        # Calculate RMSF
        if compute_rmsf and has_trajectory:
            logger.info("Calculating RMSF...")
            results["rmsf"] = self._calculate_rmsf(universe, selection)

        # Calculate radius of gyration
        if compute_rg:
            logger.info("Calculating radius of gyration...")
            results["radius_of_gyration"] = self._calculate_radius_of_gyration(
                universe,
                selection
            )

        # Analyze hydrogen bonds
        if compute_hbonds and has_trajectory:
            logger.info("Analyzing hydrogen bonds...")
            results["hydrogen_bonds"] = self._analyze_hydrogen_bonds(universe)

        # Add summary statistics
        if not has_trajectory:
            results["note"] = "Single frame analysis - trajectory-dependent metrics not computed"

        logger.info("Analysis complete!")

        return results

    def generate_cache_key(self, input_data: Any, **kwargs) -> str:
        """
        Generate deterministic cache key for MDAnalysis analysis.

        Args:
            input_data: Dictionary with topology and trajectory paths
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        # Include file paths and important parameters in cache key
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "topology": input_data.get('topology'),
            "trajectory": input_data.get('trajectory'),
            "selection": kwargs.get('selection', self.config.get('selection')),
            "compute_rmsd": kwargs.get('compute_rmsd', self.config.get('compute_rmsd')),
            "compute_rmsf": kwargs.get('compute_rmsf', self.config.get('compute_rmsf')),
            "compute_rg": kwargs.get('compute_rg', self.config.get('compute_rg')),
            "compute_hbonds": kwargs.get('compute_hbonds', self.config.get('compute_hbonds')),
            "step": self.config.get('step'),
        }

        cache_string = json.dumps(cache_dict, sort_keys=True)
        cache_key = hashlib.sha256(cache_string.encode()).hexdigest()
        return cache_key

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata.

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "description": "Molecular dynamics trajectory analysis using MDAnalysis",
            "capabilities": {
                "rmsd_calculation": True,
                "rmsf_calculation": True,
                "radius_of_gyration": True,
                "hydrogen_bond_analysis": HBA_AVAILABLE,
                "multiple_formats": True,
                "alignment": True
            },
            "supported_formats": {
                "topology": ["PDB", "PSF", "GRO", "TPR", "PRMTOP", "TOP"],
                "trajectory": ["DCD", "XTC", "TRR", "NetCDF", "H5MD", "LAMMPS", "CRD"]
            },
            "config": {
                "selection": self.config.get("selection"),
                "compute_rmsd": self.config.get("compute_rmsd"),
                "compute_rmsf": self.config.get("compute_rmsf"),
                "compute_rg": self.config.get("compute_rg"),
                "compute_hbonds": self.config.get("compute_hbonds"),
                "align_trajectory": self.config.get("align_trajectory"),
                "step": self.config.get("step")
            },
            "references": [
                {
                    "paper": "Michaud-Agrawal et al., J. Comput. Chem. 2011, 32, 2319-2327",
                    "doi": "10.1002/jcc.21787"
                },
                {
                    "paper": "Gowers et al., Proc. of the 15th Python in Science Conf. 2016, 98-105",
                    "doi": "10.25080/Majora-629e541a-00e"
                }
            ],
            "website": "https://www.mdanalysis.org/",
            "documentation": "https://docs.mdanalysis.org/",
            "computational_cost": {
                "rmsd": "Fast - linear with trajectory length",
                "rmsf": "Fast - linear with trajectory length",
                "radius_of_gyration": "Fast - linear with trajectory length",
                "hydrogen_bonds": "Moderate - depends on system size and trajectory length"
            }
        }
