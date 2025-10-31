"""
Health check endpoint for PharmForge API
Returns 200 when healthy, 503 when degraded
"""
from fastapi import APIRouter, Response, Depends
from sqlalchemy import text
from sqlalchemy.orm import Session
from redis import Redis
import json
import os
import logging
from datetime import datetime, timezone

logger = logging.getLogger(__name__)
router = APIRouter()

def check_database(db: Session) -> bool:
    """Check if database is reachable"""
    try:
        db.execute(text("SELECT 1"))
        return True
    except Exception as e:
        logger.error(f"Database health check failed: {e}")
        return False

def check_redis() -> bool:
    """Check if Redis is reachable"""
    try:
        redis_url = os.getenv("REDIS_URL", "redis://localhost:6379/0")
        r = Redis.from_url(redis_url, socket_connect_timeout=2)
        r.ping()
        return True
    except Exception as e:
        logger.error(f"Redis health check failed: {e}")
        return False

@router.get("/health")
def health(db: Session = Depends(lambda: None)):
    """
    Health check endpoint
    Returns:
        - 200 if all systems healthy
        - 503 if any system degraded

    Response includes:
        - status: "ok" or "degraded"
        - timestamp: ISO 8601 timestamp
        - db: boolean indicating database health
        - redis: boolean indicating Redis health
    """
    # Check database health
    db_ok = False
    if db is not None:
        db_ok = check_database(db)
    else:
        # Fallback: create temporary session
        try:
            from ..db.database import SessionLocal
            temp_db = SessionLocal()
            try:
                db_ok = check_database(temp_db)
            finally:
                temp_db.close()
        except Exception as e:
            logger.error(f"Failed to create database session: {e}")
            db_ok = False

    # Check Redis health
    redis_ok = check_redis()

    checks = {
        "db": db_ok,
        "redis": redis_ok
    }

    healthy = all(checks.values())
    status_code = 200 if healthy else 503

    response_body = {
        "status": "ok" if healthy else "degraded",
        "timestamp": datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
        **checks
    }

    return Response(
        content=json.dumps(response_body),
        media_type="application/json",
        status_code=status_code
    )


@router.get("/")
def root():
    """Root health check - simple ping"""
    return {"status": "ok", "service": "PharmForge API"}


@router.get("/cache/stats")
def cache_stats():
    """
    Get cache statistics

    Returns:
        Cache hit/miss statistics and health status
    """
    from backend.core.cache import get_cache

    cache = get_cache()
    stats = cache.get_stats()
    healthy = cache.healthcheck()

    return {
        "healthy": healthy,
        **stats
    }
