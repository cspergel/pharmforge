"""
Health check endpoint for PharmForge API
Returns 200 when healthy, 503 when degraded
"""
from fastapi import APIRouter, Response
from sqlalchemy import text
from redis import Redis
import json
import os
import logging

logger = logging.getLogger(__name__)
router = APIRouter()

def check_database() -> bool:
    """Check if database is reachable"""
    try:
        from ..db.database import engine
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
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
def health():
    """
    Health check endpoint
    Returns:
        - 200 if all systems healthy
        - 503 if any system degraded
    """
    db_ok = check_database()
    redis_ok = check_redis()

    checks = {
        "db": db_ok,
        "redis": redis_ok
    }

    healthy = all(checks.values())
    status_code = 200 if healthy else 503

    response_body = {
        "status": "ok" if healthy else "degraded",
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
