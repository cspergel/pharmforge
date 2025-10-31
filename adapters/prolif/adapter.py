"""
ProLIF (Protein-Ligand Interaction Fingerprints) Adapter for PharmForge

Analyzes protein-ligand interactions from docking poses or MD trajectories
to identify interaction types and frequencies.

ProLIF computes interaction fingerprints for protein-ligand complexes, detecting:
- Hydrogen bonds (H-bonds)
- Hydrophobic contacts
- Pi-stacking interactions
- Pi-cation interactions
- Salt bridges
- Halogen bonds
- Metal coordination
- Water bridges

Features:
- Multiple frame support (MD trajectories)
- Residue-level interaction analysis
- Interaction frequency calculations
- Support for multiple input formats
- Compatible with MDAnalysis Universe objects

Reference: https://prolif.readthedocs.io/
"""

import logging
import asyncio
from typing import Any, Dict, Optional, List, Tuple
from pathlib import Path
import json
import hashlib

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)

# Try to import ProLIF and dependencies
try:
    import prolif as plf
    PROLIF_AVAILABLE = True
    logger.info(f"ProLIF version: {plf.__version__}")
except ImportError:
    PROLIF_AVAILABLE = False
    plf = None
    logger.warning("ProLIF not available - install with: pip install prolif")

# Try to import MDAnalysis (required by ProLIF)
try:
    import MDAnalysis as mda
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False
    mda = None
    logger.warning("MDAnalysis not available - install with: pip install MDAnalysis")

# Try to import pandas (used for results processing)
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    pd = None
    logger.warning("Pandas not available - install with: pip install pandas")


class ProLIFAdapter(AdapterProtocol):
    """
    Adapter for ProLIF protein-ligand interaction fingerprint analysis.

    Analyzes protein-ligand complexes to identify and quantify interaction types.
    Supports both single-frame analysis (docking poses) and multi-frame analysis
    (MD trajectories).

    Features:
    - Comprehensive interaction detection (H-bonds, hydrophobic, pi-stacking, etc.)
    - Residue-level interaction mapping
    - Interaction frequency calculation
    - Support for multiple frames/trajectories
    - Integration with MDAnalysis Universe objects
    """

    # Available interaction types in ProLIF
    INTERACTION_TYPES = [
        "HBDonor",           # Hydrogen bond (protein as donor)
        "HBAcceptor",        # Hydrogen bond (protein as acceptor)
        "PiStacking",        # Pi-stacking interactions
        "PiCation",          # Pi-cation interactions
        "CationPi",          # Cation-pi interactions
        "Hydrophobic",       # Hydrophobic contacts
        "Anionic",           # Salt bridges (anionic)
        "Cationic",          # Salt bridges (cationic)
        "XBDonor",           # Halogen bond (protein as donor)
        "XBAcceptor",        # Halogen bond (protein as acceptor)
        "MetalDonor",        # Metal coordination (protein as donor)
        "MetalAcceptor",     # Metal coordination (protein as acceptor)
    ]

    def __init__(
        self,
        name: str = "prolif",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize ProLIF adapter.

        Args:
            name: Adapter name (default: "prolif")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - protein_selection: MDAnalysis selection for protein (default: "protein")
                   - ligand_selection: MDAnalysis selection for ligand (default: "resname LIG")
                   - interactions: List of interaction types to compute (default: all)
                   - distance_cutoff: Distance cutoff for interactions in Angstroms (default: varies by type)
                   - compute_frequency: Calculate interaction frequencies (default: True)
        """
        default_config = {
            "protein_selection": "protein",
            "ligand_selection": "resname LIG",
            "interactions": None,  # None means use all available
            "compute_frequency": True,
            "distance_cutoff": None,  # Use ProLIF defaults
        }

        # Merge with provided config
        merged_config = {**default_config, **(config or {})}

        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        if not PROLIF_AVAILABLE:
            logger.error("ProLIF is not installed! Install with: pip install prolif")
        if not MDANALYSIS_AVAILABLE:
            logger.error("MDAnalysis is not installed! Install with: pip install MDAnalysis")
        if not PANDAS_AVAILABLE:
            logger.warning("Pandas is not installed - some features may be limited")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input contains required protein and ligand file paths or Universe.

        Args:
            input_data: Dictionary with file paths or MDAnalysis Universe:
                       Option 1: {'protein_file': path, 'ligand_file': path}
                       Option 2: {'complex_file': path, 'ligand_selection': selection}
                       Option 3: {'universe': mda.Universe, 'ligand_selection': selection}

        Returns:
            True if valid, False otherwise
        """
        if not PROLIF_AVAILABLE or not MDANALYSIS_AVAILABLE:
            logger.error("ProLIF or MDAnalysis is not installed")
            return False

        if not isinstance(input_data, dict):
            logger.error("Input must be a dictionary")
            return False

        # Check for different input formats
        has_separate_files = 'protein_file' in input_data and 'ligand_file' in input_data
        has_complex_file = 'complex_file' in input_data
        has_universe = 'universe' in input_data

        if not (has_separate_files or has_complex_file or has_universe):
            logger.error("Input must contain either 'protein_file' and 'ligand_file', "
                        "'complex_file', or 'universe'")
            return False

        # Validate ligand selection is provided when needed
        if (has_complex_file or has_universe) and 'ligand_selection' not in input_data:
            logger.error("'ligand_selection' is required when using 'complex_file' or 'universe'")
            return False

        # Check if files exist
        if has_separate_files:
            protein_path = Path(input_data['protein_file'])
            ligand_path = Path(input_data['ligand_file'])
            if not protein_path.exists():
                logger.error(f"Protein file does not exist: {protein_path}")
                return False
            if not ligand_path.exists():
                logger.error(f"Ligand file does not exist: {ligand_path}")
                return False

        if has_complex_file:
            complex_path = Path(input_data['complex_file'])
            if not complex_path.exists():
                logger.error(f"Complex file does not exist: {complex_path}")
                return False

        return True

    def _load_structures(
        self,
        input_data: Dict[str, Any]
    ) -> Tuple[Optional[Any], Optional[Any]]:
        """
        Load protein and ligand structures from input data.

        Args:
            input_data: Dictionary with file paths or Universe

        Returns:
            Tuple of (protein_molecule, ligand_molecule) as ProLIF Molecule objects
        """
        try:
            # Case 1: Separate protein and ligand files
            if 'protein_file' in input_data and 'ligand_file' in input_data:
                logger.info(f"Loading protein from: {input_data['protein_file']}")
                logger.info(f"Loading ligand from: {input_data['ligand_file']}")

                protein = plf.Molecule.from_file(input_data['protein_file'])
                ligand = plf.Molecule.from_file(input_data['ligand_file'])

                return protein, ligand

            # Case 2: Complex file with ligand selection
            elif 'complex_file' in input_data:
                logger.info(f"Loading complex from: {input_data['complex_file']}")

                # Load as MDAnalysis Universe
                universe = mda.Universe(input_data['complex_file'])

                protein_sel = self.config.get('protein_selection', 'protein')
                ligand_sel = input_data.get('ligand_selection', self.config.get('ligand_selection'))

                logger.info(f"Protein selection: {protein_sel}")
                logger.info(f"Ligand selection: {ligand_sel}")

                # Create ProLIF molecules from selections
                protein = plf.Molecule.from_mda(universe.select_atoms(protein_sel))
                ligand = plf.Molecule.from_mda(universe.select_atoms(ligand_sel))

                return protein, ligand

            # Case 3: MDAnalysis Universe directly
            elif 'universe' in input_data:
                universe = input_data['universe']

                protein_sel = self.config.get('protein_selection', 'protein')
                ligand_sel = input_data.get('ligand_selection', self.config.get('ligand_selection'))

                logger.info(f"Using provided Universe")
                logger.info(f"Protein selection: {protein_sel}")
                logger.info(f"Ligand selection: {ligand_sel}")

                protein = plf.Molecule.from_mda(universe.select_atoms(protein_sel))
                ligand = plf.Molecule.from_mda(universe.select_atoms(ligand_sel))

                return protein, ligand

            else:
                logger.error("Invalid input format")
                return None, None

        except Exception as e:
            logger.error(f"Failed to load structures: {e}")
            return None, None

    def _compute_fingerprint(
        self,
        protein: Any,
        ligand: Any,
        n_frames: int = 1
    ) -> Optional[Any]:
        """
        Compute interaction fingerprint using ProLIF.

        Args:
            protein: ProLIF Molecule object for protein
            ligand: ProLIF Molecule object for ligand
            n_frames: Number of frames to analyze (for trajectories)

        Returns:
            ProLIF Fingerprint object
        """
        try:
            # Get interaction types to compute
            interactions = self.config.get('interactions')
            if interactions is None:
                # Use all available interactions
                fp = plf.Fingerprint()
            else:
                # Use specified interactions
                fp = plf.Fingerprint(interactions)

            logger.info(f"Computing fingerprint for {n_frames} frame(s)")
            logger.info(f"Interactions: {fp.interactions if hasattr(fp, 'interactions') else 'all'}")

            # Run fingerprint calculation
            if n_frames == 1:
                # Single frame
                fp.run_from_iterable([ligand], protein)
            else:
                # Multiple frames (trajectory)
                # Note: ligand should be an iterable of Molecule objects, one per frame
                fp.run_from_iterable(ligand, protein)

            logger.info(f"Fingerprint computation complete")
            return fp

        except Exception as e:
            logger.error(f"Fingerprint computation failed: {e}")
            return None

    def _analyze_fingerprint(
        self,
        fp: Any,
        compute_frequency: bool = True
    ) -> Dict[str, Any]:
        """
        Analyze fingerprint and extract interaction data.

        Args:
            fp: ProLIF Fingerprint object
            compute_frequency: Whether to compute interaction frequencies

        Returns:
            Dictionary with interaction analysis results
        """
        try:
            import numpy as np

            # Convert fingerprint to DataFrame if pandas is available
            if PANDAS_AVAILABLE:
                df = fp.to_dataframe()
                logger.info(f"Fingerprint DataFrame shape: {df.shape}")
            else:
                df = None

            # Get interaction counts per frame
            results = {
                "n_frames": len(fp.ifp),
                "interactions": {},
                "residue_pairs": [],
            }

            # Analyze each frame
            all_interactions = []
            residue_pairs_set = set()

            for frame_idx, frame_data in enumerate(fp.ifp):
                frame_interactions = {}

                for residue_key, interaction_dict in frame_data.items():
                    residue_id = str(residue_key)

                    for interaction_type, value in interaction_dict.items():
                        if value:  # Interaction is present
                            # Track interaction type
                            if interaction_type not in frame_interactions:
                                frame_interactions[interaction_type] = []
                            frame_interactions[interaction_type].append(residue_id)

                            # Track unique residue pairs
                            residue_pairs_set.add((residue_id, interaction_type))

                all_interactions.append(frame_interactions)

            # Compute overall statistics
            interaction_counts = {}
            for frame_interactions in all_interactions:
                for interaction_type, residues in frame_interactions.items():
                    if interaction_type not in interaction_counts:
                        interaction_counts[interaction_type] = 0
                    interaction_counts[interaction_type] += len(residues)

            results["interaction_counts"] = interaction_counts

            # Compute frequencies if requested
            if compute_frequency and results["n_frames"] > 0:
                frequencies = {}
                for interaction_type, count in interaction_counts.items():
                    frequencies[interaction_type] = count / results["n_frames"]
                results["interaction_frequencies"] = frequencies

            # Convert residue pairs to list
            results["residue_pairs"] = [
                {"residue": res_id, "interaction_type": int_type}
                for res_id, int_type in sorted(residue_pairs_set)
            ]
            results["n_unique_interactions"] = len(residue_pairs_set)

            # Add summary statistics
            total_interactions = sum(interaction_counts.values())
            results["total_interactions"] = total_interactions
            results["avg_interactions_per_frame"] = (
                total_interactions / results["n_frames"] if results["n_frames"] > 0 else 0
            )

            # Most common interaction types
            if interaction_counts:
                sorted_interactions = sorted(
                    interaction_counts.items(),
                    key=lambda x: x[1],
                    reverse=True
                )
                results["most_common_interactions"] = [
                    {"type": itype, "count": count}
                    for itype, count in sorted_interactions[:5]
                ]

            logger.info(f"Found {total_interactions} total interactions across {results['n_frames']} frame(s)")
            logger.info(f"Unique interaction types: {list(interaction_counts.keys())}")

            return results

        except Exception as e:
            logger.error(f"Fingerprint analysis failed: {e}")
            return {"error": str(e)}

    def _run_analysis(
        self,
        input_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Run the full ProLIF analysis workflow (synchronous method for thread pool).

        Args:
            input_data: Dictionary with structure file paths or Universe

        Returns:
            Dictionary with analysis results
        """
        # Load structures
        logger.info("Loading structures...")
        protein, ligand = self._load_structures(input_data)

        if protein is None or ligand is None:
            return {"error": "Failed to load protein and ligand structures"}

        # Determine number of frames
        n_frames = 1
        if hasattr(ligand, '__len__'):
            n_frames = len(ligand)
        elif 'universe' in input_data:
            universe = input_data['universe']
            if hasattr(universe, 'trajectory'):
                n_frames = universe.trajectory.n_frames

        # Compute fingerprint
        logger.info("Computing interaction fingerprint...")
        fp = self._compute_fingerprint(protein, ligand, n_frames)

        if fp is None:
            return {"error": "Failed to compute interaction fingerprint"}

        # Analyze fingerprint
        logger.info("Analyzing interactions...")
        results = self._analyze_fingerprint(
            fp,
            compute_frequency=self.config.get('compute_frequency', True)
        )

        # Add metadata
        results["input_type"] = self._get_input_type(input_data)
        results["protein_selection"] = self.config.get('protein_selection', 'protein')
        results["ligand_selection"] = input_data.get(
            'ligand_selection',
            self.config.get('ligand_selection', 'resname LIG')
        )

        logger.info("ProLIF analysis complete!")

        return results

    def _get_input_type(self, input_data: Dict[str, Any]) -> str:
        """Determine the input type from input_data."""
        if 'protein_file' in input_data and 'ligand_file' in input_data:
            return "separate_files"
        elif 'complex_file' in input_data:
            return "complex_file"
        elif 'universe' in input_data:
            return "mdanalysis_universe"
        else:
            return "unknown"

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute ProLIF interaction fingerprint analysis.

        Args:
            input_data: Dictionary containing structure information:
                       - Option 1: {'protein_file': path, 'ligand_file': path}
                       - Option 2: {'complex_file': path, 'ligand_selection': selection}
                       - Option 3: {'universe': mda.Universe, 'ligand_selection': selection}
            **kwargs: Additional parameters:
                     - protein_selection: Override protein selection
                     - ligand_selection: Override ligand selection
                     - interactions: Override interaction types to compute

        Returns:
            AdapterResult containing interaction fingerprint results
        """
        try:
            # Check if dependencies are available
            if not PROLIF_AVAILABLE:
                return AdapterResult(
                    success=False,
                    data={},
                    error="ProLIF is not installed. Install with: pip install prolif"
                )

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
                    error="Invalid input. See documentation for supported input formats."
                )

            # Update config with runtime overrides
            if 'protein_selection' in kwargs:
                self.config['protein_selection'] = kwargs['protein_selection']
            if 'ligand_selection' in kwargs:
                input_data['ligand_selection'] = kwargs['ligand_selection']
            if 'interactions' in kwargs:
                self.config['interactions'] = kwargs['interactions']

            logger.info("Starting ProLIF interaction analysis")

            # Run analysis in thread pool (ProLIF/MDAnalysis operations are synchronous)
            loop = asyncio.get_event_loop()
            result_data = await loop.run_in_executor(
                None,
                self._run_analysis,
                input_data
            )

            # Check for errors
            if "error" in result_data:
                return AdapterResult(
                    success=False,
                    data=result_data,
                    error=result_data["error"]
                )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "version": self.version,
                    "prolif_version": plf.__version__ if PROLIF_AVAILABLE else None,
                    "input_type": result_data.get("input_type", "unknown")
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
                    "installation_help": "pip install prolif MDAnalysis"
                }
            )

        except Exception as e:
            logger.error(f"ProLIF analysis failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={"adapter_name": self.name}
            )

    def generate_cache_key(self, input_data: Any, **kwargs) -> str:
        """
        Generate deterministic cache key for ProLIF analysis.

        Args:
            input_data: Dictionary with structure file paths or Universe
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        # Create cache dict with relevant parameters
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "protein_selection": self.config.get('protein_selection'),
            "ligand_selection": input_data.get('ligand_selection', self.config.get('ligand_selection')),
            "interactions": self.config.get('interactions'),
            "compute_frequency": self.config.get('compute_frequency'),
        }

        # Add file paths if present
        if 'protein_file' in input_data:
            cache_dict['protein_file'] = str(input_data['protein_file'])
        if 'ligand_file' in input_data:
            cache_dict['ligand_file'] = str(input_data['ligand_file'])
        if 'complex_file' in input_data:
            cache_dict['complex_file'] = str(input_data['complex_file'])

        # Universe objects can't be easily hashed, so use a placeholder
        if 'universe' in input_data:
            cache_dict['has_universe'] = True

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
            "description": "Protein-ligand interaction fingerprint analysis using ProLIF",
            "capabilities": {
                "hydrogen_bonds": True,
                "hydrophobic_contacts": True,
                "pi_stacking": True,
                "pi_cation": True,
                "salt_bridges": True,
                "halogen_bonds": True,
                "metal_coordination": True,
                "trajectory_support": True,
                "residue_level_analysis": True,
                "frequency_calculation": True,
            },
            "supported_interactions": self.INTERACTION_TYPES,
            "supported_input_formats": {
                "separate_files": "protein_file and ligand_file paths",
                "complex_file": "single PDB/MOL2 file with ligand_selection",
                "mdanalysis_universe": "MDAnalysis Universe object with ligand_selection",
            },
            "config": {
                "protein_selection": self.config.get("protein_selection"),
                "ligand_selection": self.config.get("ligand_selection"),
                "interactions": self.config.get("interactions"),
                "compute_frequency": self.config.get("compute_frequency"),
            },
            "requirements": [
                "ProLIF (pip install prolif)",
                "MDAnalysis (pip install MDAnalysis)",
                "Pandas (pip install pandas, optional)",
            ],
            "references": [
                {
                    "paper": "ProLIF: Protein-Ligand Interaction Fingerprints",
                    "url": "https://prolif.readthedocs.io/",
                }
            ],
            "computational_cost": {
                "single_frame": "Fast - seconds for typical protein-ligand complex",
                "trajectory": "Moderate - linear with number of frames",
            },
        }
