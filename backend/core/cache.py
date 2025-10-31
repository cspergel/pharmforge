"""
Redis Cache Manager for PharmForge

Provides a caching layer for adapter results with:
- Configurable TTL (Time To Live)
- JSON serialization/deserialization
- Graceful degradation if Redis unavailable
- Cache hit/miss logging
"""
import os
import json
import logging
from typing import Any, Optional
import redis
from redis.exceptions import RedisError, ConnectionError

logger = logging.getLogger(__name__)


class CacheManager:
    """
    Manages Redis caching for PharmForge adapters

    Features:
    - JSON serialization of Python objects
    - Configurable TTL per key
    - Graceful degradation when Redis unavailable
    - Cache statistics logging
    """

    def __init__(
        self,
        redis_url: Optional[str] = None,
        default_ttl: int = 86400,  # 24 hours in seconds
        enabled: bool = True
    ):
        """
        Initialize cache manager

        Args:
            redis_url: Redis connection URL (defaults to REDIS_URL env var)
            default_ttl: Default time-to-live in seconds (default: 24 hours)
            enabled: Whether caching is enabled (defaults to ENABLE_CACHING env var)
        """
        self.redis_url = redis_url or os.getenv("REDIS_URL", "redis://localhost:6379/0")
        self.default_ttl = default_ttl

        # Check if caching is explicitly disabled via environment variable
        enable_caching_env = os.getenv("ENABLE_CACHING", "true").lower()
        self.enabled = enabled and enable_caching_env in ("true", "1", "yes")

        self.redis_client: Optional[redis.Redis] = None
        self._connect()

        # Statistics
        self.hits = 0
        self.misses = 0
        self.errors = 0

    def _connect(self):
        """
        Establish connection to Redis
        Falls back to disabled mode if connection fails
        """
        if not self.enabled:
            logger.info("Cache is disabled via configuration")
            return

        try:
            self.redis_client = redis.from_url(
                self.redis_url,
                decode_responses=True,  # Auto-decode bytes to strings
                socket_connect_timeout=5,
                socket_timeout=5
            )
            # Test connection
            self.redis_client.ping()
            logger.info(f"Connected to Redis at {self.redis_url}")
        except (RedisError, ConnectionError) as e:
            logger.warning(f"Failed to connect to Redis: {e}. Caching disabled.")
            self.redis_client = None
            self.enabled = False

    def get(self, key: str) -> Optional[Any]:
        """
        Get value from cache

        Args:
            key: Cache key

        Returns:
            Cached value (deserialized from JSON) or None if not found
        """
        if not self.enabled or self.redis_client is None:
            return None

        try:
            value = self.redis_client.get(key)
            if value is None:
                self.misses += 1
                logger.debug(f"Cache MISS: {key[:16]}...")
                return None

            # Deserialize JSON
            result = json.loads(value)
            self.hits += 1
            logger.debug(f"Cache HIT: {key[:16]}...")
            return result

        except (RedisError, json.JSONDecodeError) as e:
            self.errors += 1
            logger.error(f"Cache GET error for key {key[:16]}...: {e}")
            return None

    def set(
        self,
        key: str,
        value: Any,
        ttl: Optional[int] = None
    ) -> bool:
        """
        Set value in cache

        Args:
            key: Cache key
            value: Value to cache (must be JSON-serializable)
            ttl: Time-to-live in seconds (uses default_ttl if None)

        Returns:
            True if successfully cached, False otherwise
        """
        if not self.enabled or self.redis_client is None:
            return False

        try:
            # Serialize to JSON
            serialized = json.dumps(value)

            # Use default TTL if not specified
            expire_seconds = ttl if ttl is not None else self.default_ttl

            # Set with expiration
            self.redis_client.setex(key, expire_seconds, serialized)
            logger.debug(f"Cache SET: {key[:16]}... (TTL: {expire_seconds}s)")
            return True

        except (RedisError, TypeError, ValueError) as e:
            self.errors += 1
            logger.error(f"Cache SET error for key {key[:16]}...: {e}")
            return False

    def delete(self, key: str) -> bool:
        """
        Delete key from cache

        Args:
            key: Cache key to delete

        Returns:
            True if deleted, False otherwise
        """
        if not self.enabled or self.redis_client is None:
            return False

        try:
            result = self.redis_client.delete(key)
            logger.debug(f"Cache DELETE: {key[:16]}... (existed: {result > 0})")
            return result > 0

        except RedisError as e:
            self.errors += 1
            logger.error(f"Cache DELETE error for key {key[:16]}...: {e}")
            return False

    def clear(self, pattern: str = "*") -> int:
        """
        Clear cache keys matching pattern

        Args:
            pattern: Redis key pattern (default: "*" for all keys)

        Returns:
            Number of keys deleted
        """
        if not self.enabled or self.redis_client is None:
            return 0

        try:
            # Find all matching keys
            keys = self.redis_client.keys(pattern)
            if not keys:
                logger.info(f"No keys found matching pattern: {pattern}")
                return 0

            # Delete all matching keys
            deleted = self.redis_client.delete(*keys)
            logger.info(f"Cache CLEAR: deleted {deleted} keys matching '{pattern}'")
            return deleted

        except RedisError as e:
            self.errors += 1
            logger.error(f"Cache CLEAR error for pattern {pattern}: {e}")
            return 0

    def get_stats(self) -> dict:
        """
        Get cache statistics

        Returns:
            Dictionary with cache statistics
        """
        total_requests = self.hits + self.misses
        hit_rate = (self.hits / total_requests * 100) if total_requests > 0 else 0.0

        return {
            "enabled": self.enabled,
            "hits": self.hits,
            "misses": self.misses,
            "errors": self.errors,
            "total_requests": total_requests,
            "hit_rate": f"{hit_rate:.2f}%"
        }

    def reset_stats(self):
        """Reset cache statistics"""
        self.hits = 0
        self.misses = 0
        self.errors = 0
        logger.info("Cache statistics reset")

    def healthcheck(self) -> bool:
        """
        Check if Redis connection is healthy

        Returns:
            True if Redis is accessible, False otherwise
        """
        if not self.enabled or self.redis_client is None:
            return False

        try:
            self.redis_client.ping()
            return True
        except RedisError:
            return False


# Global cache instance
# This is initialized when the module is imported
_cache_instance: Optional[CacheManager] = None


def get_cache() -> CacheManager:
    """
    Get the global cache instance (singleton pattern)

    Returns:
        CacheManager instance
    """
    global _cache_instance
    if _cache_instance is None:
        _cache_instance = CacheManager()
    return _cache_instance


def reset_cache():
    """
    Reset the global cache instance
    Useful for testing
    """
    global _cache_instance
    _cache_instance = None
