"""
gmx_MMPBSA Adapter for PharmForge

Calculate binding free energies from molecular dynamics trajectories using
Molecular Mechanics Poisson-Boltzmann Surface Area (MM/PBSA) and
Molecular Mechanics Generalized Born Surface Area (MM/GBSA) methods.

gmx_MMPBSA is a Python package for calculating binding free energies from
GROMACS molecular dynamics trajectories. It's a wrapper around AMBER's
MMPBSA.py tool optimized for GROMACS workflows.

Features:
- MM/PBSA and MM/GBSA free energy calculations
- Per-residue energy decomposition
- Support for multiple trajectory formats (.xtc, .trr)
- Compatible with AMBER, GROMACS, CHARMM, NAMD topologies
- Detailed energy component breakdown

Reference: https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA
Paper: Valdes-Tresanco et al., J. Chem. Theory Comput. 2021, 17(10), 6281-6291
"""

import os
import tempfile
import asyncio
import logging
import json
import subprocess
from pathlib import Path
from typing import Any, Dict, Optional, List, Tuple
import shutil

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class GmxMMPBSAAdapter(AdapterProtocol):
    """
    Adapter for gmx_MMPBSA binding free energy calculations.

    Calculates binding free energies from molecular dynamics trajectories
    using MM/PBSA or MM/GBSA methods. Provides detailed energy decomposition
    including:
    - Van der Waals energy
    - Electrostatic energy
    - Polar solvation energy
    - Non-polar solvation (SASA)
    - Per-residue decomposition (optional)

    Requirements:
        - gmx_MMPBSA (pip install gmx_MMPBSA)
        - GROMACS (for trajectory processing)
        - AmberTools (installed with gmx_MMPBSA)

    Workflow:
        1. Validate input files (topology, trajectory, index)
        2. Generate gmx_MMPBSA input file
        3. Run MM/PBSA or MM/GBSA calculation
        4. Parse output for energy components
        5. Return structured results with binding free energy
    """

    def __init__(
        self,
        name: str = "gmx_mmpbsa",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize gmx_MMPBSA adapter.

        Args:
            name: Adapter name (default: "gmx_mmpbsa")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - method: "pb" or "gb" (default: "pb")
                   - startframe: First frame to analyze (default: 1)
                   - endframe: Last frame to analyze (default: -1, all frames)
                   - interval: Frame interval (default: 1)
                   - decomp: Enable per-residue decomposition (default: False)
                   - igb: GB model (1-8, default: 5 for GB)
                   - saltcon: Salt concentration in M (default: 0.15)
                   - temperature: Temperature in K (default: 298.15)
                   - gmx_path: Path to GROMACS binary (default: "gmx")
                   - mmpbsa_path: Path to gmx_MMPBSA binary (default: "gmx_MMPBSA")
        """
        default_config = {
            "method": "pb",  # Options: "pb", "gb"
            "startframe": 1,
            "endframe": -1,  # -1 means all frames
            "interval": 1,
            "decomp": False,  # Per-residue decomposition
            "igb": 5,  # GB model (for GB method)
            "saltcon": 0.15,  # Salt concentration (M)
            "temperature": 298.15,  # Temperature (K)
            "gmx_path": "gmx",
            "mmpbsa_path": "gmx_MMPBSA",
        }

        # Merge with provided config
        merged_config = {**default_config, **(config or {})}

        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Validate method
        if self.config["method"] not in ["pb", "gb"]:
            raise ValueError(f"Invalid method: {self.config['method']}. Must be 'pb' or 'gb'")

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data structure.

        Args:
            input_data: Dictionary with required keys:
                       - topology_file: Path to topology file
                       - trajectory_file: Path to trajectory file
                       - index_file: Path to GROMACS index file
                       - receptor_group: Receptor group name
                       - ligand_group: Ligand group name

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, dict):
            logger.error("Input must be a dictionary")
            return False

        required_keys = ["topology_file", "trajectory_file", "index_file",
                        "receptor_group", "ligand_group"]

        for key in required_keys:
            if key not in input_data:
                logger.error(f"Missing required key: {key}")
                return False

        # Validate file paths
        for file_key in ["topology_file", "trajectory_file", "index_file"]:
            file_path = input_data[file_key]
            if not isinstance(file_path, str):
                logger.error(f"{file_key} must be a string path")
                return False

        # Validate group names
        for group_key in ["receptor_group", "ligand_group"]:
            if not isinstance(input_data[group_key], str):
                logger.error(f"{group_key} must be a string")
                return False

        return True

    def _check_dependencies(self) -> Tuple[bool, str]:
        """
        Check if required dependencies are installed.

        Returns:
            Tuple of (success, error_message)
        """
        try:
            # Check gmx_MMPBSA
            result = subprocess.run(
                [self.config["mmpbsa_path"], "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode != 0:
                return False, f"gmx_MMPBSA not found at: {self.config['mmpbsa_path']}"

            logger.info(f"gmx_MMPBSA version: {result.stdout.strip()}")

        except FileNotFoundError:
            return False, (
                "gmx_MMPBSA not installed. Install with:\n"
                "  conda install -c conda-forge ambertools\n"
                "  pip install gmx_MMPBSA"
            )
        except subprocess.TimeoutExpired:
            return False, "gmx_MMPBSA version check timed out"

        try:
            # Check GROMACS (optional but recommended)
            result = subprocess.run(
                [self.config["gmx_path"], "--version"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                logger.info("GROMACS found")
            else:
                logger.warning("GROMACS not found - some features may be limited")

        except FileNotFoundError:
            logger.warning("GROMACS not found - processing may be limited")
        except subprocess.TimeoutExpired:
            logger.warning("GROMACS version check timed out")

        return True, ""

    def _verify_files(self, input_data: Dict[str, Any]) -> Tuple[bool, str]:
        """
        Verify that input files exist and are readable.

        Args:
            input_data: Input dictionary with file paths

        Returns:
            Tuple of (success, error_message)
        """
        for file_key in ["topology_file", "trajectory_file", "index_file"]:
            file_path = input_data[file_key]

            if not os.path.exists(file_path):
                return False, f"File not found: {file_path}"

            if not os.path.isfile(file_path):
                return False, f"Not a file: {file_path}"

            if os.path.getsize(file_path) == 0:
                return False, f"File is empty: {file_path}"

        logger.info("All input files verified")
        return True, ""

    def _generate_input_file(
        self,
        input_data: Dict[str, Any],
        output_path: str
    ) -> bool:
        """
        Generate gmx_MMPBSA input file.

        Args:
            input_data: Input parameters
            output_path: Path to write input file

        Returns:
            True if successful, False otherwise
        """
        try:
            method = self.config["method"]

            # Base input configuration
            input_content = f"""# gmx_MMPBSA input file generated by PharmForge
# Method: MM/{method.upper()}SA

&general
startframe={self.config['startframe']},
endframe={self.config['endframe']},
interval={self.config['interval']},
"""

            # Add method-specific parameters
            if method == "pb":
                input_content += f"""
/

&pb
istrng={self.config['saltcon']},
fillratio=4.0,
"""
            else:  # gb
                input_content += f"""
/

&gb
igb={self.config['igb']},
saltcon={self.config['saltcon']},
"""

            # Add decomposition if requested
            if self.config.get("decomp", False):
                input_content += """
/

&decomp
idecomp=1,
dec_verbose=3,
"""

            input_content += "/\n"

            # Write input file
            with open(output_path, 'w') as f:
                f.write(input_content)

            logger.info(f"Generated gmx_MMPBSA input file: {output_path}")
            return True

        except Exception as e:
            logger.error(f"Failed to generate input file: {e}")
            return False

    async def _run_mmpbsa(
        self,
        input_file: str,
        topology: str,
        trajectory: str,
        index: str,
        receptor_group: str,
        ligand_group: str,
        work_dir: str
    ) -> Tuple[bool, str, str]:
        """
        Run gmx_MMPBSA calculation.

        Args:
            input_file: Path to gmx_MMPBSA input file
            topology: Path to topology file
            trajectory: Path to trajectory file
            index: Path to index file
            receptor_group: Receptor group name
            ligand_group: Ligand group name
            work_dir: Working directory for calculation

        Returns:
            Tuple of (success, stdout, stderr)
        """
        try:
            # Build command
            cmd = [
                self.config["mmpbsa_path"],
                "-i", input_file,
                "-cs", topology,  # Complex structure
                "-ct", trajectory,  # Complex trajectory
                "-n", index,  # Index file
                "-cg", receptor_group, ligand_group,  # Complex groups
                "-o", "FINAL_RESULTS_MMPBSA.dat",
                "-eo", "FINAL_RESULTS_MMPBSA.csv",
            ]

            logger.info(f"Running gmx_MMPBSA: {' '.join(cmd)}")

            # Run gmx_MMPBSA as subprocess
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=work_dir
            )

            stdout, stderr = await process.communicate()

            stdout_text = stdout.decode() if stdout else ""
            stderr_text = stderr.decode() if stderr else ""

            if process.returncode != 0:
                logger.error(f"gmx_MMPBSA failed with return code {process.returncode}")
                logger.error(f"STDERR: {stderr_text}")
                return False, stdout_text, stderr_text

            logger.info("gmx_MMPBSA calculation completed successfully")
            return True, stdout_text, stderr_text

        except Exception as e:
            logger.error(f"Error running gmx_MMPBSA: {e}")
            return False, "", str(e)

    def _parse_results(self, work_dir: str) -> Optional[Dict[str, Any]]:
        """
        Parse gmx_MMPBSA output files.

        Args:
            work_dir: Working directory containing output files

        Returns:
            Dictionary with parsed results, or None on failure
        """
        try:
            # Look for output files
            dat_file = os.path.join(work_dir, "FINAL_RESULTS_MMPBSA.dat")
            csv_file = os.path.join(work_dir, "FINAL_RESULTS_MMPBSA.csv")

            if not os.path.exists(dat_file):
                logger.error("Output file not found: FINAL_RESULTS_MMPBSA.dat")
                return None

            # Parse DAT file for detailed results
            results = {
                "binding_free_energy": None,
                "delta_h": None,
                "delta_s": None,
                "components": {},
                "per_residue_decomposition": [],
                "warnings": []
            }

            with open(dat_file, 'r') as f:
                lines = f.readlines()

            # Parse energy components
            # Look for DELTA TOTAL section
            parsing_delta = False
            for i, line in enumerate(lines):
                if "DELTA TOTAL" in line or "DELTA Energy" in line:
                    parsing_delta = True
                    continue

                if parsing_delta:
                    # Parse energy line: Energy Component = Value +/- Std
                    if "=" in line and "+/-" in line:
                        parts = line.split("=")
                        if len(parts) == 2:
                            component = parts[0].strip()
                            value_str = parts[1].split("+/-")[0].strip()

                            try:
                                value = float(value_str)

                                # Map component names
                                if "VDWAALS" in component or "van der Waals" in component:
                                    results["components"]["van_der_waals"] = value
                                elif "EEL" in component or "Electrostatic" in component:
                                    results["components"]["electrostatic"] = value
                                elif "EPB" in component or "Polar" in component or "PB" in component:
                                    results["components"]["polar_solvation"] = value
                                elif "ENPOLAR" in component or "SASA" in component or "Non-Polar" in component:
                                    results["components"]["sasa"] = value
                                elif "DELTA G" in component or "TOTAL" in component:
                                    results["binding_free_energy"] = value
                                elif "DELTA H" in component:
                                    results["delta_h"] = value
                                elif "-T" in component and "DELTA S" in component:
                                    # This is -T*Delta S
                                    results["delta_s"] = value

                            except ValueError:
                                continue

                    # Stop parsing after empty line
                    if line.strip() == "":
                        break

            # Parse CSV file if available (may have additional data)
            if os.path.exists(csv_file):
                try:
                    import csv
                    with open(csv_file, 'r') as f:
                        reader = csv.DictReader(f)
                        csv_data = list(reader)

                    # CSV may contain frame-by-frame data
                    if csv_data:
                        logger.info(f"Found {len(csv_data)} frames in CSV output")

                except Exception as e:
                    logger.warning(f"Could not parse CSV file: {e}")

            # Validate results
            if results["binding_free_energy"] is None:
                logger.warning("Could not extract binding free energy from output")
                results["warnings"].append("Binding free energy not found in output")

            logger.info(f"Parsed results: ΔG = {results['binding_free_energy']} kcal/mol")

            return results

        except Exception as e:
            logger.error(f"Error parsing results: {e}")
            return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute gmx_MMPBSA binding free energy calculation.

        Args:
            input_data: Dictionary with required keys:
                       - topology_file: Path to topology file (.prmtop, .gro, etc.)
                       - trajectory_file: Path to trajectory file (.xtc, .trr, etc.)
                       - index_file: Path to GROMACS index file (.ndx)
                       - receptor_group: Receptor group name from index
                       - ligand_group: Ligand group name from index
            **kwargs: Additional parameters (override config):
                     - method: "pb" or "gb"
                     - startframe, endframe, interval
                     - decomp: Enable decomposition

        Returns:
            AdapterResult containing binding free energy and energy components
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input data. See logs for details."
            )

        # Check dependencies
        deps_ok, deps_error = self._check_dependencies()
        if not deps_ok:
            return AdapterResult(
                success=False,
                data=None,
                error=deps_error
            )

        # Verify files exist
        files_ok, files_error = self._verify_files(input_data)
        if not files_ok:
            return AdapterResult(
                success=False,
                data=None,
                error=files_error
            )

        # Extract input parameters
        topology_file = input_data["topology_file"]
        trajectory_file = input_data["trajectory_file"]
        index_file = input_data["index_file"]
        receptor_group = input_data["receptor_group"]
        ligand_group = input_data["ligand_group"]

        # Create temporary working directory
        temp_dir = tempfile.mkdtemp(prefix="gmx_mmpbsa_")

        try:
            logger.info(f"Working directory: {temp_dir}")
            logger.info(f"Topology: {topology_file}")
            logger.info(f"Trajectory: {trajectory_file}")
            logger.info(f"Index: {index_file}")
            logger.info(f"Receptor: {receptor_group}, Ligand: {ligand_group}")

            # Generate input file
            input_file = os.path.join(temp_dir, "mmpbsa.in")
            if not self._generate_input_file(input_data, input_file):
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to generate gmx_MMPBSA input file"
                )

            # Run gmx_MMPBSA
            success, stdout, stderr = await self._run_mmpbsa(
                input_file=input_file,
                topology=topology_file,
                trajectory=trajectory_file,
                index=index_file,
                receptor_group=receptor_group,
                ligand_group=ligand_group,
                work_dir=temp_dir
            )

            if not success:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"gmx_MMPBSA calculation failed: {stderr}"
                )

            # Parse results
            results = self._parse_results(temp_dir)
            if results is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to parse gmx_MMPBSA output"
                )

            # Add metadata
            results["method"] = self.config["method"].upper()
            results["frames_analyzed"] = {
                "start": self.config["startframe"],
                "end": self.config["endframe"],
                "interval": self.config["interval"]
            }
            results["receptor_group"] = receptor_group
            results["ligand_group"] = ligand_group

            logger.info(
                f"✓ gmx_MMPBSA complete: "
                f"ΔG = {results['binding_free_energy']:.2f} kcal/mol"
            )

            return AdapterResult(
                success=True,
                data=results,
                metadata={
                    "adapter_name": self.name,
                    "adapter_version": self.version,
                    "method": self.config["method"],
                    "temperature": self.config["temperature"],
                    "saltcon": self.config["saltcon"]
                }
            )

        except Exception as e:
            logger.error(f"gmx_MMPBSA adapter error: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data=None,
                error=str(e)
            )

        finally:
            # Cleanup temporary directory
            try:
                shutil.rmtree(temp_dir)
                logger.info(f"Cleaned up temporary directory: {temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to cleanup temp directory: {e}")

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
            "description": "Calculate binding free energies using MM/PBSA and MM/GBSA methods",
            "capabilities": {
                "mm_pbsa": True,
                "mm_gbsa": True,
                "energy_decomposition": True,
                "per_residue_decomposition": True,
                "trajectory_analysis": True
            },
            "requirements": [
                "gmx_MMPBSA (pip install gmx_MMPBSA)",
                "AmberTools (conda install -c conda-forge ambertools)",
                "GROMACS (optional, for preprocessing)",
                "Topology file (.prmtop, .gro, etc.)",
                "Trajectory file (.xtc, .trr, etc.)",
                "Index file (.ndx)"
            ],
            "config": {
                "method": self.config["method"],
                "startframe": self.config["startframe"],
                "endframe": self.config["endframe"],
                "interval": self.config["interval"],
                "decomp": self.config["decomp"],
                "saltcon": self.config["saltcon"],
                "temperature": self.config["temperature"]
            },
            "reference": {
                "paper": "Valdes-Tresanco et al., J. Chem. Theory Comput. 2021, 17(10), 6281-6291",
                "github": "https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA",
                "documentation": "https://valdes-tresanco-ms.github.io/gmx_MMPBSA/"
            },
            "input_format": {
                "topology_file": "Path to topology file",
                "trajectory_file": "Path to MD trajectory file",
                "index_file": "Path to GROMACS index file",
                "receptor_group": "Receptor group name from index",
                "ligand_group": "Ligand group name from index"
            },
            "output_format": {
                "binding_free_energy": "kcal/mol (lower is better)",
                "delta_h": "Enthalpy change (kcal/mol)",
                "delta_s": "Entropy change (-T*ΔS, kcal/mol)",
                "components": {
                    "van_der_waals": "Van der Waals energy",
                    "electrostatic": "Electrostatic energy",
                    "polar_solvation": "Polar solvation energy",
                    "sasa": "Non-polar solvation (SASA)"
                }
            },
            "computational_cost": {
                "typical_runtime": "1-30 minutes depending on trajectory length",
                "memory_usage": "1-4 GB",
                "parallelization": "Limited (mainly single-threaded)"
            }
        }
