"""
PharmForge FastAPI Application
Main entry point for the backend API
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create FastAPI application
app = FastAPI(
    title="PharmForge",
    description="AI-powered drug discovery workflow orchestrator",
    version="0.1.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # TODO: Restrict in production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Import routers
from .api import health, runs, stock_routes
from .api.adapters import router as adapters_router

# Include routers
app.include_router(health.router, tags=["health"])
app.include_router(adapters_router, prefix="/api/adapters", tags=["adapters"])
app.include_router(runs.router, prefix="/api/v1", tags=["runs"])
app.include_router(stock_routes.router, tags=["stock"])

@app.on_event("startup")
async def startup_event():
    """Initialize services on startup"""
    logger.info("PharmForge API starting up...")

    # Register all adapters
    from .core.adapter_registry import register_all_adapters
    register_all_adapters()

    # TODO: Initialize database connection pool
    # TODO: Initialize Redis connection pool

    logger.info("PharmForge API ready!")

@app.on_event("shutdown")
async def shutdown_event():
    """Cleanup on shutdown"""
    logger.info("PharmForge API shutting down...")
    # TODO: Close database connections
    # TODO: Close Redis connections
    logger.info("PharmForge API stopped")

@app.get("/")
async def root():
    """Root endpoint"""
    return {
        "name": "PharmForge API",
        "version": "0.1.0",
        "status": "running",
        "docs": "/docs"
    }
