"""
Celery Application Configuration
Handles async task execution for pipeline runs
"""
import os
from celery import Celery
from kombu import serialization
import logging

logger = logging.getLogger(__name__)

# Get environment variables
CELERY_BROKER_URL = os.getenv("CELERY_BROKER_URL", os.getenv("REDIS_URL", "redis://localhost:6379/0"))
CELERY_RESULT_BACKEND = os.getenv("CELERY_RESULT_BACKEND", os.getenv("REDIS_URL", "redis://localhost:6379/0"))

# Create Celery application
celery_app = Celery(
    "pharmforge",
    broker=CELERY_BROKER_URL,
    backend=CELERY_RESULT_BACKEND,
    include=["backend.core.tasks"]  # Auto-discover tasks
)

# Configure Celery
celery_app.conf.update(
    # Serialization
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",

    # Timezone
    timezone="UTC",
    enable_utc=True,

    # Time limits
    task_time_limit=3600,  # 1 hour hard limit
    task_soft_time_limit=3000,  # 50 minutes soft limit

    # Results
    result_expires=86400,  # 24 hours
    result_persistent=True,

    # Task execution
    task_acks_late=True,
    task_reject_on_worker_lost=True,

    # Worker
    worker_prefetch_multiplier=1,  # Only fetch one task at a time
    worker_max_tasks_per_child=100,  # Restart worker after 100 tasks

    # Logging
    worker_log_format="[%(asctime)s: %(levelname)s/%(processName)s] %(message)s",
    worker_task_log_format="[%(asctime)s: %(levelname)s/%(processName)s][%(task_name)s(%(task_id)s)] %(message)s",
)

logger.info(f"Celery app configured with broker: {CELERY_BROKER_URL}")
logger.info(f"Celery result backend: {CELERY_RESULT_BACKEND}")

# Auto-discover tasks from backend.core module
celery_app.autodiscover_tasks(["backend.core"], force=True)
