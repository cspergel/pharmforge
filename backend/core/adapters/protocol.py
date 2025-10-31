"""
Adapter Protocol - Base class for all PharmForge adapters
Defines the interface that all adapters must implement
"""
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional, List
from dataclasses import dataclass
import hashlib
import json
import logging
from backend.core.cache import get_cache

logger = logging.getLogger(__name__)


@dataclass
class AdapterResult:
    """
    Standard result format returned by all adapters
    """
    success: bool
    data: Any
    error: Optional[str] = None
    cache_hit: bool = False
    metadata: Optional[Dict[str, Any]] = None


class AdapterProtocol(ABC):
    """
    Base class for all adapters
    All adapters must inherit from this class and implement the required methods
    """

    def __init__(self, name: str, adapter_type: str, config: Optional[Dict[str, Any]] = None):
        """
        Initialize the adapter

        Args:
            name: Unique identifier for this adapter (e.g., "pubchem", "vina_docking")
            adapter_type: Type of adapter ("api", "local", "ml")
            config: Optional configuration dictionary
        """
        self.name = name
        self.adapter_type = adapter_type
        self.config = config or {}
        self.enabled = True
        self.version = "0.1.0"

        logger.info(f"Initialized adapter: {self.name} (type: {self.adapter_type})")

    @abstractmethod
    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the adapter logic

        Args:
            input_data: Primary input data (e.g., SMILES string, compound ID)
            **kwargs: Additional adapter-specific parameters

        Returns:
            AdapterResult containing the execution results
        """
        pass

    @abstractmethod
    def validate_input(self, input_data: Any) -> bool:
        """
        Validate input data before execution

        Args:
            input_data: Input to validate

        Returns:
            True if valid, False otherwise
        """
        pass

    def generate_cache_key(self, input_data: Any, **kwargs) -> str:
        """
        Generate a deterministic cache key for this adapter execution

        Args:
            input_data: Primary input data
            **kwargs: Additional parameters that affect the result

        Returns:
            SHA256 hash as cache key
        """
        # Create a dictionary with all inputs
        cache_dict = {
            "adapter": self.name,
            "version": self.version,
            "input": input_data,
            "params": kwargs
        }

        # Convert to JSON string (sorted keys for determinism)
        cache_string = json.dumps(cache_dict, sort_keys=True)

        # Generate SHA256 hash
        cache_key = hashlib.sha256(cache_string.encode()).hexdigest()

        return cache_key

    def get_metadata(self) -> Dict[str, Any]:
        """
        Get adapter metadata

        Returns:
            Dictionary containing adapter information
        """
        return {
            "name": self.name,
            "type": self.adapter_type,
            "version": self.version,
            "enabled": self.enabled,
            "config": self.config
        }

    async def __call__(self, input_data: Any, use_cache: bool = True, **kwargs) -> AdapterResult:
        """
        Make the adapter callable with caching support
        Allows usage like: result = await adapter(smiles)

        Args:
            input_data: Primary input data
            use_cache: Whether to use cache (default: True)
            **kwargs: Additional adapter-specific parameters

        Returns:
            AdapterResult with cache_hit flag set appropriately
        """
        cache = get_cache()

        # Generate cache key
        cache_key = self.generate_cache_key(input_data, **kwargs)

        # Try to get from cache if enabled
        if use_cache and cache.enabled:
            cached_data = cache.get(cache_key)
            if cached_data is not None:
                logger.debug(f"Cache hit for {self.name}: {cache_key[:16]}...")
                return AdapterResult(
                    success=True,
                    data=cached_data,
                    cache_hit=True,
                    metadata={"cache_key": cache_key}
                )

        # Execute adapter logic
        result = await self.execute(input_data, **kwargs)

        # Cache successful results
        if use_cache and result.success and cache.enabled:
            cache.set(cache_key, result.data)
            logger.debug(f"Cached result for {self.name}: {cache_key[:16]}...")

        # Add cache metadata
        if result.metadata is None:
            result.metadata = {}
        result.metadata["cache_key"] = cache_key

        return result


class AdapterRegistry:
    """
    Registry for managing all available adapters
    """

    def __init__(self):
        self._adapters: Dict[str, AdapterProtocol] = {}
        logger.info("Adapter registry initialized")

    def register(self, adapter: AdapterProtocol):
        """
        Register an adapter

        Args:
            adapter: AdapterProtocol instance to register
        """
        if adapter.name in self._adapters:
            logger.warning(f"Adapter {adapter.name} already registered, overwriting")

        self._adapters[adapter.name] = adapter
        logger.info(f"Registered adapter: {adapter.name}")

    def get(self, name: str) -> Optional[AdapterProtocol]:
        """
        Get an adapter by name

        Args:
            name: Adapter name

        Returns:
            AdapterProtocol instance or None if not found
        """
        return self._adapters.get(name)

    def list_adapters(self) -> List[str]:
        """
        List all registered adapter names

        Returns:
            List of adapter names
        """
        return list(self._adapters.keys())

    def get_all_metadata(self) -> List[Dict[str, Any]]:
        """
        Get metadata for all registered adapters

        Returns:
            List of metadata dictionaries
        """
        return [adapter.get_metadata() for adapter in self._adapters.values()]


# Global adapter registry instance
registry = AdapterRegistry()
