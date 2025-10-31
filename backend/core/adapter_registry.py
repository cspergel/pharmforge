"""
Adapter Registry Initialization
Registers all available adapters at startup
"""
import logging
from backend.core.adapters.protocol import registry

logger = logging.getLogger(__name__)


def register_all_adapters():
    """
    Register all available adapters (75 adapters)
    Called at application startup

    Latest additions:
    - xTB (quantum chemistry - semiempirical QM)
    - IntAct (protein-protein interactions)
    - ORD (Open Reaction Database - forward synthesis)
    - RNAcentral (RNA databases)
    - COCONUT (natural products)
    - HMDB (metabolomics)
    - Tox21 (toxicity)
    - CompTox (EPA chemistry)
    - SAbDab (antibody structures)
    - ImmuneBuilder (antibody prediction)
    """
    logger.info("Registering all PharmForge adapters...")

    # Import all adapter classes
    # Chemical Databases
    from adapters.pubchem.adapter import PubChemAdapter
    from adapters.chembl.adapter import ChEMBLAdapter
    from adapters.chemspider.adapter import ChemSpiderAdapter
    from adapters.drugcentral.adapter import DrugCentralAdapter
    from adapters.zinc_fragments.adapter import ZINCFragmentsAdapter
    from adapters.surechembl.adapter import SureChEMBLAdapter
    from adapters.coconut.adapter import COCONUTAdapter

    # Target & Disease
    from adapters.opentargets.adapter import OpenTargetsAdapter
    from adapters.disgenet.adapter import DisGeNETAdapter
    from adapters.uniprot.adapter import UniProtAdapter
    from adapters.string_db.adapter import StringDBAdapter
    from adapters.biogrid.adapter import BioGRIDAdapter
    from adapters.intact.adapter import IntActAdapter

    # Protein Structure
    from adapters.alphafold.adapter import AlphaFoldAdapter
    from adapters.rcsb_pdb.adapter import RCSBPDBAdapter
    from adapters.pdbe.adapter import PDBEAdapter
    from adapters.swissmodel.adapter import SwissModelAdapter
    from adapters.pdb_redo.adapter import PDBRedoAdapter
    from adapters.sabdab.adapter import SAbDabAdapter
    from adapters.immunebuilder.adapter import ImmuneBuilderAdapter

    # Docking & Binding
    from adapters.vina.adapter import VinaAdapter
    from adapters.gnina.adapter import GNINAAdapter
    from adapters.diffdock.adapter import DiffDockAdapter
    from adapters.bindingdb.adapter import BindingDBAdapter
    from adapters.openmm.adapter import OpenMMAdapter
    from adapters.gmx_mmpbsa.adapter import GmxMMPBSAAdapter
    from adapters.molscrub.adapter import MolScrubAdapter

    # ADMET & Properties
    from adapters.admet_ai.adapter import ADMETaiAdapter
    from adapters.rdkit_local.adapter import RDKitAdapter
    from adapters.pkcsm.adapter import PkCSMAdapter
    from adapters.swisstarget.adapter import SwissTargetAdapter
    from adapters.targetnet.adapter import TargetNetAdapter
    from adapters.olorenchemengine.adapter import OlorenChemEngineAdapter
    from adapters.comptox.adapter import CompToxAdapter
    from adapters.tox21.adapter import Tox21Adapter
    from adapters.xtb.adapter import xTBAdapter

    # Retrosynthesis & De Novo
    from adapters.aizynthfinder.adapter import AiZynthFinderAdapter
    from adapters.llm_retrosynthesis.adapter import LLMRetrosynthesisAdapter
    from adapters.reinvent.adapter import REINVENTAdapter
    from adapters.molgan.adapter import MolGANAdapter
    from adapters.denovo.adapter import DeNovoAdapter
    from adapters.ord.adapter import ORDAdapter

    # Clinical & Safety
    from adapters.clinicaltrials.adapter import ClinicalTrialsAdapter
    from adapters.fda_faers.adapter import FDAFAERSAdapter

    # Literature & Patents
    from adapters.pubmed.adapter import PubMedAdapter
    from adapters.europepmc.adapter import EuropePMCAdapter
    from adapters.google_patents.adapter import GooglePatentsAdapter
    from adapters.lens.adapter import LensAdapter

    # Genomics & Expression
    from adapters.geo.adapter import GEOAdapter
    from adapters.gtex.adapter import GTExAdapter
    from adapters.reactome.adapter import ReactomeAdapter
    from adapters.kegg.adapter import KEGGAdapter
    from adapters.brenda.adapter import BRENDAAdapter

    # Metabolomics
    from adapters.hmdb.adapter import HMDBAdapter

    # RNA Databases
    from adapters.rnacentral.adapter import RNAcentralAdapter

    # Molecular Features & ML Round 1 (from awesome-drug-discovery)
    from adapters.mordred.adapter import MordredAdapter
    from adapters.molvs.adapter import MolVSAdapter
    from adapters.deepchem.adapter import DeepChemAdapter
    from adapters.scikit_mol.adapter import ScikitMolAdapter
    from adapters.chemprop.adapter import ChempropAdapter
    from adapters.mdanalysis.adapter import MDAnalysisAdapter
    from adapters.openbabel.adapter import OpenBabelAdapter
    from adapters.meeko.adapter import MeekoAdapter

    # Molecular Features & ML Round 2 (GNN, Featurization, Visualization, Optimization)
    from adapters.torchdrug.adapter import TorchDrugAdapter
    from adapters.dgllifesci.adapter import DGLLifeSciAdapter
    from adapters.molfeat.adapter import MolFeatAdapter
    from adapters.chemplot.adapter import ChemPlotAdapter
    from adapters.optuna.adapter import OptunaAdapter

    # Molecular Features & ML Round 3 (Interaction Analysis, Modern Toolkit, 3D Viz)
    from adapters.prolif.adapter import ProLIFAdapter
    from adapters.datamol.adapter import DatamolAdapter
    from adapters.py3dmol.adapter import Py3DmolAdapter
    from adapters.pymol.adapter import PyMOLAdapter

    # Molecular Features & ML Round 4 (AutoML & Feature Engineering)
    from adapters.chemml.adapter import ChemMLAdapter
    from adapters.tpot.adapter import TPOTAdapter
    from adapters.autosklearn.adapter import AutoSklearnAdapter

    # Create adapter instances - using error handling for each
    adapters_to_register = []

    adapter_classes = [
        # Chemical Databases
        (PubChemAdapter, "PubChem"),
        (ChEMBLAdapter, "ChEMBL"),
        (ChemSpiderAdapter, "ChemSpider"),
        (DrugCentralAdapter, "DrugCentral"),
        (ZINCFragmentsAdapter, "ZINC Fragments"),
        (SureChEMBLAdapter, "SureChEMBL"),
        (COCONUTAdapter, "COCONUT"),
        # Target & Disease
        (OpenTargetsAdapter, "Open Targets"),
        (DisGeNETAdapter, "DisGeNET"),
        (UniProtAdapter, "UniProt"),
        (StringDBAdapter, "STRING DB"),
        (BioGRIDAdapter, "BioGRID"),
        (IntActAdapter, "IntAct"),
        # Protein Structure
        (AlphaFoldAdapter, "AlphaFold"),
        (RCSBPDBAdapter, "RCSB PDB"),
        (PDBEAdapter, "PDBe"),
        (SwissModelAdapter, "SWISS-MODEL"),
        (PDBRedoAdapter, "PDB-REDO"),
        (SAbDabAdapter, "SAbDab"),
        (ImmuneBuilderAdapter, "ImmuneBuilder"),
        # Docking & Binding
        (VinaAdapter, "AutoDock Vina"),
        (GNINAAdapter, "GNINA"),
        (DiffDockAdapter, "DiffDock"),
        (BindingDBAdapter, "BindingDB"),
        (OpenMMAdapter, "OpenMM"),
        (GmxMMPBSAAdapter, "gmx_MMPBSA"),
        (MolScrubAdapter, "MolScrub"),
        # ADMET & Properties
        (ADMETaiAdapter, "ADMET-ai"),
        (RDKitAdapter, "RDKit"),
        (PkCSMAdapter, "pkCSM"),
        (SwissTargetAdapter, "SwissTargetPrediction"),
        (TargetNetAdapter, "TargetNet"),
        (OlorenChemEngineAdapter, "Oloren ChemEngine"),
        (CompToxAdapter, "EPA CompTox"),
        (Tox21Adapter, "Tox21"),
        (xTBAdapter, "xTB (Quantum Chemistry)"),
        # Retrosynthesis & De Novo
        (AiZynthFinderAdapter, "AiZynthFinder"),
        (LLMRetrosynthesisAdapter, "LLM Retrosynthesis"),
        (REINVENTAdapter, "REINVENT"),
        (MolGANAdapter, "MolGAN"),
        (DeNovoAdapter, "De Novo Design"),
        (ORDAdapter, "Open Reaction Database"),
        # Clinical & Safety
        (ClinicalTrialsAdapter, "ClinicalTrials.gov"),
        (FDAFAERSAdapter, "FDA FAERS"),
        # Literature & Patents
        (PubMedAdapter, "PubMed"),
        (EuropePMCAdapter, "Europe PMC"),
        (GooglePatentsAdapter, "Google Patents"),
        (LensAdapter, "Lens.org"),
        # Genomics & Expression
        (GEOAdapter, "GEO"),
        (GTExAdapter, "GTEx"),
        (ReactomeAdapter, "Reactome"),
        (KEGGAdapter, "KEGG"),
        (BRENDAAdapter, "BRENDA"),
        # Metabolomics
        (HMDBAdapter, "HMDB"),
        # RNA Databases
        (RNAcentralAdapter, "RNAcentral"),
        # Molecular Features & ML Round 1
        (MordredAdapter, "Mordred"),
        (MolVSAdapter, "MolVS"),
        (DeepChemAdapter, "DeepChem"),
        (ScikitMolAdapter, "scikit-mol"),
        (ChempropAdapter, "Chemprop"),
        (MDAnalysisAdapter, "MDAnalysis"),
        (OpenBabelAdapter, "OpenBabel"),
        (MeekoAdapter, "Meeko"),
        # Molecular Features & ML Round 2
        (TorchDrugAdapter, "TorchDrug"),
        (DGLLifeSciAdapter, "DGL-LifeSci"),
        (MolFeatAdapter, "MolFeat"),
        (ChemPlotAdapter, "ChemPlot"),
        (OptunaAdapter, "Optuna"),
        # Molecular Features & ML Round 3
        (ProLIFAdapter, "ProLIF"),
        (DatamolAdapter, "Datamol"),
        (Py3DmolAdapter, "py3Dmol"),
        (PyMOLAdapter, "PyMOL"),
        # Molecular Features & ML Round 4
        (ChemMLAdapter, "ChemML"),
        (TPOTAdapter, "TPOT"),
        (AutoSklearnAdapter, "Auto-sklearn"),
    ]

    # Instantiate and register each adapter
    for adapter_class, display_name in adapter_classes:
        try:
            adapter_instance = adapter_class()
            adapters_to_register.append(adapter_instance)
            logger.info(f"  Created {display_name} adapter instance")
        except Exception as e:
            logger.error(f"  Failed to create {display_name} adapter: {e}")

    # Register all successfully created adapters
    registered_count = 0
    failed_count = 0

    for adapter in adapters_to_register:
        try:
            registry.register(adapter)
            logger.info(f"  ✓ Registered: {adapter.name} (type: {adapter.adapter_type})")
            registered_count += 1
        except Exception as e:
            logger.error(f"  ✗ Failed to register {adapter.name}: {e}")
            failed_count += 1

    logger.info(f"\n{'='*70}")
    logger.info(f"Adapter Registry Initialization Complete")
    logger.info(f"  Total registered: {registered_count} adapters")
    logger.info(f"  Failed: {failed_count} adapters")
    logger.info(f"{'='*70}\n")

    if registered_count > 0:
        logger.info(f"Available adapters: {', '.join(sorted(registry.list_adapters()))}")


def get_adapter_status():
    """
    Get status of all registered adapters

    Returns:
        Dictionary with adapter names and their status
    """
    status = {}

    for adapter_name in registry.list_adapters():
        adapter = registry.get(adapter_name)
        if adapter:
            status[adapter_name] = {
                "enabled": adapter.enabled,
                "type": adapter.adapter_type,
                "version": adapter.version
            }

    return status


def assert_required_adapters():
    """
    Assert that required adapters are registered.

    Raises:
        AssertionError: If a required adapter is not registered
    """
    required = ["vina_docking", "tdc_admet", "aizynthfinder"]

    missing_adapters = []
    for name in required:
        adapter = registry.get(name)
        if adapter is None:
            missing_adapters.append(name)

    if missing_adapters:
        raise AssertionError(
            f"Required adapters not registered: {', '.join(missing_adapters)}. "
            f"Available adapters: {', '.join(registry.list_adapters())}"
        )

    logger.info(f"✓ All required adapters registered: {', '.join(required)}")
