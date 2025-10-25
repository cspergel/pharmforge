"""
Adapter Registry Initialization
Registers all available adapters at startup
"""
import logging
from backend.core.adapters.protocol import registry

logger = logging.getLogger(__name__)


def register_all_adapters():
    """
    Register all available adapters
    Called at application startup
    """
    logger.info("Registering adapters...")

    # Import adapters
    from adapters.pubchem.adapter import PubChemAdapter
    from adapters.chembl.adapter import ChEMBLAdapter
    from adapters.rdkit_local.adapter import RDKitAdapter
    from adapters.admet_ai.adapter import ADMETaiAdapter
    from adapters.aizynthfinder.adapter import AiZynthFinderAdapter
    from adapters.vina.adapter import VinaAdapter

    # Register adapters
    adapters_to_register = [
        PubChemAdapter(),
        ChEMBLAdapter(),
        RDKitAdapter(),
        ADMETaiAdapter(),
        AiZynthFinderAdapter(),
        VinaAdapter(),
    ]

    for adapter in adapters_to_register:
        try:
            registry.register(adapter)
            logger.info(f"✓ Registered adapter: {adapter.name} (type: {adapter.adapter_type})")
        except Exception as e:
            logger.error(f"✗ Failed to register adapter {adapter.name}: {e}")

    logger.info(f"Registry initialized with {len(registry.list_adapters())} adapters")
    logger.info(f"Available adapters: {', '.join(registry.list_adapters())}")


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
