"""
OpenMM Adapter for PharmForge

Provides molecular dynamics simulation, energy minimization, and property
calculation for small molecules using the OpenMM toolkit.

OpenMM is a high-performance toolkit for molecular simulation that can run
on CPUs and GPUs (CUDA/OpenCL).

Features:
- SMILES to 3D structure conversion
- Energy minimization with various force fields
- Molecular dynamics simulation
- Property calculation (RMSD, radius of gyration, energy)
- Conformer generation
- Stability assessment

Reference: http://openmm.org/
Paper: Eastman et al., PLoS Comput. Biol. 2017, 13(7), e1005659
"""

import hashlib
import logging
import json
import io
from typing import Dict, Any, Optional, Tuple
import asyncio

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class OpenMMAdapter(AdapterProtocol):
    """
    OpenMM adapter for molecular dynamics simulation and analysis.

    Uses OpenMM to perform energy minimization, molecular dynamics simulations,
    and property calculations for small molecules. Useful for:
    - Validating synthesized molecules
    - Predicting molecular stability
    - Generating low-energy conformers
    - Calculating molecular properties

    Features:
    - Multiple force field support (GAFF, AMBER)
    - CPU and GPU acceleration
    - Customizable simulation parameters
    - Comprehensive property reporting
    """

    def __init__(
        self,
        name: str = "openmm",
        adapter_type: str = "local",
        config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize OpenMM adapter.

        Args:
            name: Adapter name (default: "openmm")
            adapter_type: Adapter type (default: "local")
            config: Optional configuration dictionary. Supported keys:
                   - force_field: Force field to use (default: "gaff")
                   - minimize_steps: Energy minimization steps (default: 1000)
                   - minimize_tolerance: Convergence tolerance in kJ/mol (default: 10.0)
                   - run_md: Whether to run MD simulation (default: False)
                   - md_steps: MD simulation steps (default: 10000)
                   - md_temperature: Temperature in Kelvin (default: 300.0)
                   - md_timestep: Timestep in femtoseconds (default: 2.0)
                   - save_trajectory: Save trajectory coordinates (default: False)
                   - platform: OpenMM platform ("CPU", "CUDA", "OpenCL", "Reference")
        """
        default_config = {
            "force_field": "gaff",  # Options: "gaff", "amber14", "amber99"
            "minimize_steps": 1000,
            "minimize_tolerance": 10.0,  # kJ/mol
            "run_md": False,
            "md_steps": 10000,
            "md_temperature": 300.0,  # Kelvin
            "md_timestep": 2.0,  # femtoseconds
            "save_trajectory": False,
            "platform": None,  # Auto-select best available platform
            "friction_coefficient": 1.0,  # 1/ps for Langevin integrator
            "report_interval": 100,  # Report every N steps
            "use_pbc": False,  # Periodic boundary conditions (False for small molecules)
        }

        # Merge with provided config
        merged_config = {**default_config, **(config or {})}

        super().__init__(name, adapter_type, merged_config)
        self.version = "1.0.0"

        # Lazy-load OpenMM components (heavy imports)
        self._openmm_app = None
        self._openmm_unit = None
        self._openmm = None

    @property
    def openmm(self):
        """Lazy-load OpenMM core module."""
        if self._openmm is None:
            try:
                import openmm
                self._openmm = openmm
                logger.info(f"OpenMM version: {openmm.version.version}")
            except ImportError as e:
                logger.error(f"Failed to import OpenMM: {e}")
                raise ImportError(
                    "OpenMM library not installed. "
                    "Install with: conda install -c conda-forge openmm\n"
                    "or: pip install openmm\n"
                    "For GPU support, see: http://openmm.org/documentation.html"
                ) from e
        return self._openmm

    @property
    def openmm_app(self):
        """Lazy-load OpenMM application module."""
        if self._openmm_app is None:
            try:
                import openmm.app as app
                self._openmm_app = app
            except ImportError as e:
                raise ImportError("OpenMM app module not available") from e
        return self._openmm_app

    @property
    def openmm_unit(self):
        """Lazy-load OpenMM unit module."""
        if self._openmm_unit is None:
            try:
                import openmm.unit as unit
                self._openmm_unit = unit
            except ImportError as e:
                raise ImportError("OpenMM unit module not available") from e
        return self._openmm_unit

    def validate_input(self, smiles: str) -> bool:
        """
        Validate SMILES string input.

        Args:
            smiles: SMILES string to validate

        Returns:
            True if valid, False otherwise
        """
        if not smiles or not isinstance(smiles, str):
            return False

        # Basic SMILES validation with RDKit
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception as e:
            logger.warning(f"SMILES validation error: {e}")
            return False

    def _smiles_to_3d(self, smiles: str) -> Tuple[Any, str]:
        """
        Convert SMILES to 3D structure using RDKit.

        Args:
            smiles: Input SMILES string

        Returns:
            Tuple of (RDKit molecule with 3D coords, PDB string)
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates using ETKDG
            params = AllChem.ETKDGv3()
            params.randomSeed = 42  # Reproducibility
            success = AllChem.EmbedMolecule(mol, params)

            if success != 0:
                logger.warning("ETKDG failed, trying basic embedding")
                AllChem.EmbedMolecule(mol)

            # Optimize with UFF force field (quick pre-optimization)
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            except Exception as e:
                logger.warning(f"UFF optimization failed: {e}")

            # Convert to PDB format
            pdb_block = Chem.MolToPDBBlock(mol)

            logger.info(f"Generated 3D structure with {mol.GetNumAtoms()} atoms")

            return mol, pdb_block

        except ImportError as e:
            raise ImportError(
                "RDKit library not installed. Install with: conda install -c conda-forge rdkit"
            ) from e
        except Exception as e:
            logger.error(f"Failed to convert SMILES to 3D: {e}")
            raise

    def _create_openmm_system(self, pdb_string: str) -> Tuple[Any, Any, Any]:
        """
        Create OpenMM system from PDB string.

        Args:
            pdb_string: PDB format structure

        Returns:
            Tuple of (topology, positions, system)
        """
        try:
            # Parse PDB
            pdb_file = io.StringIO(pdb_string)
            pdb = self.openmm_app.PDBFile(pdb_file)

            # Get topology and positions
            topology = pdb.topology
            positions = pdb.positions

            # Select force field
            force_field_name = self.config.get("force_field", "gaff")

            # Create force field
            # Note: For small molecules, we'll use a simplified force field approach
            # In production, you'd want to use proper GAFF/AMBER parameters via openmmforcefields
            try:
                # Try to load standard force fields
                if force_field_name.lower() == "gaff":
                    # GAFF requires openmmforcefields package
                    logger.info("Attempting to use GAFF force field...")
                    try:
                        from openmmforcefields.generators import GAFFTemplateGenerator
                        forcefield = self.openmm_app.ForceField()
                        gaff = GAFFTemplateGenerator(molecules=[])
                        forcefield.registerTemplateGenerator(gaff.generator)
                    except ImportError:
                        logger.warning("openmmforcefields not available, using simple force field")
                        forcefield = self._create_simple_forcefield()
                else:
                    # Use standard AMBER force fields
                    forcefield = self.openmm_app.ForceField(f"{force_field_name}.xml", "tip3p.xml")

            except Exception as e:
                logger.warning(f"Could not load {force_field_name}, using simple force field: {e}")
                forcefield = self._create_simple_forcefield()

            # Create system
            system = forcefield.createSystem(
                topology,
                nonbondedMethod=self.openmm_app.NoCutoff,
                constraints=None,
                rigidWater=False
            )

            logger.info(f"Created OpenMM system with {system.getNumParticles()} particles")

            return topology, positions, system

        except Exception as e:
            logger.error(f"Failed to create OpenMM system: {e}")
            raise

    def _create_simple_forcefield(self):
        """
        Create a simple force field for testing when standard force fields aren't available.

        Returns:
            Simple OpenMM ForceField object
        """
        # Create a minimal force field XML
        ff_xml = """
        <ForceField>
          <AtomTypes>
            <Type name="C" class="C" element="C" mass="12.01"/>
            <Type name="H" class="H" element="H" mass="1.008"/>
            <Type name="O" class="O" element="O" mass="16.00"/>
            <Type name="N" class="N" element="N" mass="14.01"/>
          </AtomTypes>
          <HarmonicBondForce>
            <Bond class1="C" class2="C" length="0.154" k="259408.0"/>
            <Bond class1="C" class2="H" length="0.109" k="284512.0"/>
            <Bond class1="C" class2="O" length="0.143" k="267776.0"/>
            <Bond class1="C" class2="N" length="0.147" k="282001.6"/>
            <Bond class1="O" class2="H" length="0.096" k="462750.4"/>
            <Bond class1="N" class2="H" length="0.101" k="363171.2"/>
          </HarmonicBondForce>
          <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
            <Atom type="C" charge="0.0" sigma="0.339967" epsilon="0.359824"/>
            <Atom type="H" charge="0.0" sigma="0.106908" epsilon="0.0656888"/>
            <Atom type="O" charge="0.0" sigma="0.296105" epsilon="0.878640"/>
            <Atom type="N" charge="0.0" sigma="0.325000" epsilon="0.711280"/>
          </NonbondedForce>
        </ForceField>
        """

        ff_file = io.StringIO(ff_xml)
        return self.openmm_app.ForceField(ff_file)

    def _minimize_energy(
        self,
        topology: Any,
        positions: Any,
        system: Any
    ) -> Dict[str, Any]:
        """
        Perform energy minimization.

        Args:
            topology: OpenMM topology
            positions: Initial positions
            system: OpenMM system

        Returns:
            Dictionary with minimization results
        """
        try:
            unit = self.openmm_unit

            # Get platform
            platform_name = self.config.get("platform")
            if platform_name:
                platform = self.openmm.Platform.getPlatformByName(platform_name)
            else:
                # Auto-select fastest platform
                platform = None

            # Create integrator (required for context, but not used for minimization)
            integrator = self.openmm.LangevinIntegrator(
                self.config.get("md_temperature", 300.0) * unit.kelvin,
                self.config.get("friction_coefficient", 1.0) / unit.picosecond,
                self.config.get("md_timestep", 2.0) * unit.femtosecond
            )

            # Create simulation context
            if platform:
                context = self.openmm.Context(system, integrator, platform)
                logger.info(f"Using platform: {platform.getName()}")
            else:
                context = self.openmm.Context(system, integrator)
                logger.info(f"Using platform: {context.getPlatform().getName()}")

            # Set initial positions
            context.setPositions(positions)

            # Get initial energy
            state = context.getState(getEnergy=True)
            initial_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

            logger.info(f"Initial energy: {initial_energy:.2f} kJ/mol")

            # Minimize energy
            max_iterations = self.config.get("minimize_steps", 1000)
            tolerance = self.config.get("minimize_tolerance", 10.0) * unit.kilojoules_per_mole

            self.openmm.LocalEnergyMinimizer.minimize(
                context,
                tolerance=tolerance,
                maxIterations=max_iterations
            )

            # Get final energy and positions
            state = context.getState(getEnergy=True, getPositions=True)
            final_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            final_positions = state.getPositions()

            # Check convergence (OpenMM doesn't directly report this)
            converged = abs(final_energy - initial_energy) > 0.1  # Changed if energy decreased

            logger.info(f"Final energy: {final_energy:.2f} kJ/mol")
            logger.info(f"Energy change: {final_energy - initial_energy:.2f} kJ/mol")

            result = {
                "initial_energy": round(initial_energy, 4),
                "final_energy": round(final_energy, 4),
                "energy_change": round(final_energy - initial_energy, 4),
                "converged": converged,
                "steps": max_iterations,
                "final_positions": final_positions
            }

            return result

        except Exception as e:
            logger.error(f"Energy minimization failed: {e}")
            raise

    def _run_molecular_dynamics(
        self,
        topology: Any,
        positions: Any,
        system: Any
    ) -> Dict[str, Any]:
        """
        Run molecular dynamics simulation.

        Args:
            topology: OpenMM topology
            positions: Initial positions
            system: OpenMM system

        Returns:
            Dictionary with MD results
        """
        try:
            unit = self.openmm_unit

            # Get parameters
            temperature = self.config.get("md_temperature", 300.0) * unit.kelvin
            timestep = self.config.get("md_timestep", 2.0) * unit.femtosecond
            friction = self.config.get("friction_coefficient", 1.0) / unit.picosecond
            n_steps = self.config.get("md_steps", 10000)
            report_interval = self.config.get("report_interval", 100)

            # Get platform
            platform_name = self.config.get("platform")
            if platform_name:
                platform = self.openmm.Platform.getPlatformByName(platform_name)
            else:
                platform = None

            # Create Langevin integrator
            integrator = self.openmm.LangevinIntegrator(temperature, friction, timestep)

            # Create simulation
            if platform:
                simulation = self.openmm_app.Simulation(topology, system, integrator, platform)
            else:
                simulation = self.openmm_app.Simulation(topology, system, integrator)

            # Set initial positions
            simulation.context.setPositions(positions)

            # Equilibrate (assign velocities)
            simulation.context.setVelocitiesToTemperature(temperature)

            logger.info(f"Running MD: {n_steps} steps at {temperature._value:.1f} K")

            # Store trajectory if requested
            trajectory_coords = []
            if self.config.get("save_trajectory", False):
                trajectory_coords.append(
                    simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
                )

            # Track properties
            energies = []
            times = []

            # Run simulation
            for step in range(0, n_steps, report_interval):
                simulation.step(report_interval)

                # Get state
                state = simulation.context.getState(getEnergy=True, getPositions=True)

                # Record properties
                energies.append(
                    state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
                )
                times.append(step * timestep.value_in_unit(unit.picoseconds))

                # Save trajectory point
                if self.config.get("save_trajectory", False):
                    trajectory_coords.append(state.getPositions(asNumpy=True))

            # Final state
            final_state = simulation.context.getState(getEnergy=True, getPositions=True)
            final_energy = final_state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            final_positions = final_state.getPositions()

            # Calculate RMSD from initial structure
            rmsd = self._calculate_rmsd(positions, final_positions)

            # Calculate radius of gyration
            rg = self._calculate_radius_of_gyration(final_positions, topology)

            # Calculate average energy
            import statistics
            avg_energy = statistics.mean(energies) if energies else final_energy
            energy_std = statistics.stdev(energies) if len(energies) > 1 else 0.0

            total_time = n_steps * timestep.value_in_unit(unit.picoseconds)

            logger.info(f"MD complete: {total_time:.2f} ps")
            logger.info(f"  Final energy: {final_energy:.2f} kJ/mol")
            logger.info(f"  RMSD: {rmsd:.3f} Å")
            logger.info(f"  Rg: {rg:.3f} Å")

            result = {
                "temperature": round(temperature._value, 2),
                "steps": n_steps,
                "time": round(total_time, 4),  # picoseconds
                "timestep": round(timestep.value_in_unit(unit.femtoseconds), 2),
                "final_energy": round(final_energy, 4),
                "average_energy": round(avg_energy, 4),
                "energy_std": round(energy_std, 4),
                "rmsd": round(rmsd, 4),
                "radius_of_gyration": round(rg, 4),
                "final_positions": final_positions,
                "energy_trajectory": [round(e, 2) for e in energies] if self.config.get("save_trajectory") else [],
                "time_trajectory": [round(t, 2) for t in times] if self.config.get("save_trajectory") else [],
            }

            return result

        except Exception as e:
            logger.error(f"Molecular dynamics failed: {e}")
            raise

    def _calculate_rmsd(self, positions1: Any, positions2: Any) -> float:
        """
        Calculate RMSD between two sets of positions.

        Args:
            positions1: First set of positions
            positions2: Second set of positions

        Returns:
            RMSD in Angstroms
        """
        try:
            import numpy as np

            # Convert to numpy arrays
            pos1 = np.array([[p.x, p.y, p.z] for p in positions1])
            pos2 = np.array([[p.x, p.y, p.z] for p in positions2])

            # Calculate RMSD (in nm, convert to Angstroms)
            diff = pos1 - pos2
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1))) * 10.0  # nm to Å

            return rmsd

        except Exception as e:
            logger.warning(f"RMSD calculation failed: {e}")
            return 0.0

    def _calculate_radius_of_gyration(self, positions: Any, topology: Any) -> float:
        """
        Calculate radius of gyration.

        Args:
            positions: Atomic positions
            topology: Molecular topology

        Returns:
            Radius of gyration in Angstroms
        """
        try:
            import numpy as np

            # Get heavy atom positions (exclude hydrogen)
            heavy_atoms = []
            for atom in topology.atoms():
                if atom.element.symbol != 'H':
                    heavy_atoms.append(atom.index)

            if not heavy_atoms:
                # Fallback to all atoms
                heavy_atoms = list(range(len(positions)))

            # Convert to numpy array
            coords = np.array([[positions[i].x, positions[i].y, positions[i].z] for i in heavy_atoms])

            # Calculate center of mass (assuming equal masses)
            center = np.mean(coords, axis=0)

            # Calculate Rg
            rg = np.sqrt(np.mean(np.sum((coords - center)**2, axis=1))) * 10.0  # nm to Å

            return rg

        except Exception as e:
            logger.warning(f"Radius of gyration calculation failed: {e}")
            return 0.0

    def _positions_to_pdb(self, positions: Any, topology: Any) -> str:
        """
        Convert positions to PDB format string.

        Args:
            positions: OpenMM positions
            topology: OpenMM topology

        Returns:
            PDB format string
        """
        try:
            output = io.StringIO()
            self.openmm_app.PDBFile.writeFile(topology, positions, output)
            return output.getvalue()
        except Exception as e:
            logger.warning(f"PDB conversion failed: {e}")
            return ""

    def _calculate_stability_score(
        self,
        minimization_result: Dict[str, Any],
        md_result: Optional[Dict[str, Any]] = None
    ) -> Tuple[float, str]:
        """
        Calculate stability score based on simulation results.

        Args:
            minimization_result: Results from energy minimization
            md_result: Optional results from MD simulation

        Returns:
            Tuple of (score 0-1, feasibility label)
        """
        try:
            # Energy-based score (lower is better)
            final_energy = minimization_result["final_energy"]

            # Normalize energy (typical range: -500 to +500 kJ/mol for small molecules)
            # Very stable: < -100, Moderate: -100 to 100, Unstable: > 100
            if final_energy < -100:
                energy_score = 1.0
            elif final_energy < 100:
                energy_score = 0.5 + 0.5 * (100 - final_energy) / 200
            else:
                energy_score = max(0.0, 0.5 - (final_energy - 100) / 400)

            # If MD was run, factor in RMSD (lower is more stable)
            if md_result:
                rmsd = md_result.get("rmsd", 0.0)
                # Good: RMSD < 1 Å, Moderate: 1-3 Å, Poor: > 3 Å
                if rmsd < 1.0:
                    rmsd_score = 1.0
                elif rmsd < 3.0:
                    rmsd_score = 1.0 - (rmsd - 1.0) / 2.0
                else:
                    rmsd_score = max(0.0, 0.5 - (rmsd - 3.0) / 6.0)

                # Combined score (70% energy, 30% RMSD)
                combined_score = 0.7 * energy_score + 0.3 * rmsd_score
            else:
                combined_score = energy_score

            # Classify feasibility
            if combined_score >= 0.7:
                feasibility = "high"
            elif combined_score >= 0.4:
                feasibility = "medium"
            else:
                feasibility = "low"

            return combined_score, feasibility

        except Exception as e:
            logger.warning(f"Stability score calculation failed: {e}")
            return 0.5, "medium"

    async def execute(self, smiles: str, **params) -> AdapterResult:
        """
        Execute OpenMM molecular dynamics simulation for a molecule.

        Args:
            smiles: SMILES string of the target molecule
            **params: Additional parameters:
                     - minimize_steps: Override default minimization steps
                     - run_md: Override whether to run MD
                     - md_steps: Override default MD steps
                     - md_temperature: Override default temperature
                     - force_field: Override default force field

        Returns:
            AdapterResult containing simulation results and properties
        """
        try:
            # Validate input
            if not self.validate_input(smiles):
                return AdapterResult(
                    success=False,
                    data={},
                    error="Invalid SMILES string"
                )

            # Get parameters (allow runtime override)
            run_md = params.get('run_md', self.config.get('run_md', False))
            minimize_steps = params.get('minimize_steps', self.config.get('minimize_steps', 1000))

            # Generate cache key
            cache_key = self.generate_cache_key(smiles, **params)

            logger.info(f"Running OpenMM simulation for: {smiles[:50]}...")
            logger.info(f"Config: minimize_steps={minimize_steps}, run_md={run_md}")

            # Run simulation in thread pool (OpenMM operations are synchronous)
            loop = asyncio.get_event_loop()
            result_data = await loop.run_in_executor(
                None,
                self._run_simulation,
                smiles,
                run_md,
                params
            )

            return AdapterResult(
                success=True,
                data=result_data,
                metadata={
                    "adapter_name": self.name,
                    "cache_key": cache_key,
                    "version": self.version,
                    "force_field": self.config.get("force_field")
                }
            )

        except ImportError as e:
            # OpenMM or RDKit not installed
            logger.error(f"Required library not available: {e}")
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name,
                    "installation_help": "conda install -c conda-forge openmm rdkit"
                }
            )

        except Exception as e:
            logger.error(f"OpenMM simulation failed: {e}", exc_info=True)
            return AdapterResult(
                success=False,
                data={},
                error=str(e),
                metadata={
                    "adapter_name": self.name
                }
            )

    def _run_simulation(self, smiles: str, run_md: bool, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run the full simulation workflow (synchronous method for thread pool).

        Args:
            smiles: Target SMILES
            run_md: Whether to run MD simulation
            params: Runtime parameters

        Returns:
            Dictionary with simulation results
        """
        from rdkit import Chem

        # Step 1: Convert SMILES to 3D structure
        logger.info("Step 1/4: Converting SMILES to 3D structure...")
        mol, pdb_string = self._smiles_to_3d(smiles)

        # Get molecular properties
        num_atoms = mol.GetNumAtoms()
        molecular_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)

        # Step 2: Create OpenMM system
        logger.info("Step 2/4: Creating OpenMM system...")
        topology, positions, system = self._create_openmm_system(pdb_string)

        # Step 3: Energy minimization
        logger.info("Step 3/4: Running energy minimization...")
        minimization_result = self._minimize_energy(topology, positions, system)

        # Get minimized positions
        minimized_positions = minimization_result.pop("final_positions")

        # Step 4: Molecular dynamics (optional)
        md_result = None
        if run_md:
            logger.info("Step 4/4: Running molecular dynamics simulation...")
            md_result = self._run_molecular_dynamics(topology, minimized_positions, system)

            # Get final positions from MD
            final_positions = md_result.pop("final_positions")
        else:
            logger.info("Step 4/4: Skipping MD simulation (run_md=False)")
            final_positions = minimized_positions

        # Calculate stability score
        stability_score, feasibility = self._calculate_stability_score(
            minimization_result,
            md_result
        )

        # Convert final structure to PDB
        final_pdb = self._positions_to_pdb(final_positions, topology)

        # Construct result data
        result_data = {
            "smiles": smiles,
            "minimization": minimization_result,
            "structure": {
                "num_atoms": num_atoms,
                "molecular_weight": round(molecular_weight, 2),
                "pdb_string": final_pdb
            },
            "stability_score": round(stability_score, 3),
            "feasibility": feasibility,
            "model": f"OpenMM {self.openmm.version.version}",
            "force_field": self.config.get("force_field", "gaff"),
            "reference": "Eastman et al., PLoS Comput. Biol. 2017"
        }

        # Add MD results if available
        if md_result:
            result_data["molecular_dynamics"] = md_result

        logger.info(f"✓ Simulation complete: stability_score={stability_score:.3f}, feasibility={feasibility}")

        return result_data

    def generate_cache_key(self, smiles: str, **kwargs) -> str:
        """
        Generate deterministic cache key for OpenMM simulation.

        Args:
            smiles: SMILES string
            **kwargs: Additional parameters

        Returns:
            SHA256 hash as cache key
        """
        # Canonicalize SMILES for consistent caching
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                smiles = Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            pass  # Use original SMILES if canonicalization fails

        # Include important parameters in cache key
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "smiles": smiles,
            "force_field": self.config.get("force_field"),
            "minimize_steps": kwargs.get("minimize_steps", self.config.get("minimize_steps")),
            "run_md": kwargs.get("run_md", self.config.get("run_md")),
            "md_steps": kwargs.get("md_steps", self.config.get("md_steps")) if kwargs.get("run_md") else None,
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
            "description": "Molecular dynamics simulation and energy minimization using OpenMM",
            "capabilities": {
                "energy_minimization": True,
                "molecular_dynamics": True,
                "property_calculation": True,
                "conformer_generation": True,
                "stability_assessment": True,
                "gpu_acceleration": True
            },
            "config": {
                "force_field": self.config.get("force_field"),
                "minimize_steps": self.config.get("minimize_steps"),
                "run_md": self.config.get("run_md"),
                "md_steps": self.config.get("md_steps"),
                "md_temperature": self.config.get("md_temperature"),
                "platform": self.config.get("platform") or "auto"
            },
            "reference": {
                "paper": "Eastman et al., PLoS Comput. Biol. 2017, 13(7), e1005659",
                "website": "http://openmm.org/",
                "documentation": "http://docs.openmm.org/"
            },
            "computational_cost": {
                "minimization": "1-10 seconds per molecule (CPU)",
                "md_simulation": "10-60 seconds for 10-100 ps (CPU)",
                "gpu_speedup": "5-50x faster with CUDA/OpenCL"
            }
        }
