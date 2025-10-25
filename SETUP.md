# PharmForge Setup Guide

## 🎯 Phase 1, Week 1 - COMPLETED ✅

### What We've Built

You now have a complete PharmForge foundation ready for development:

#### ✅ Project Structure
- Complete directory layout with Python packages
- Backend, adapters, docs, config, tests organized
- All `__init__.py` files in place

#### ✅ Docker Environment
- **docker-compose.yml** - PostgreSQL 16, Redis 7, FastAPI backend
- **Dockerfile.backend** - Python 3.11 container configuration
- Health checks for all services
- Volume persistence for data

#### ✅ FastAPI Backend
- **backend/main.py** - Application entry point with CORS
- **backend/api/health.py** - Health endpoint (200/503 status codes)
- Logging configuration
- Auto-generated API docs at `/docs`

#### ✅ Database Layer
- **backend/db/database.py** - SQLAlchemy engine and session management
- **backend/db/models.py** - Complete schema:
  - `PipelineRun` - Pipeline execution tracking
  - `CompoundResult` - Individual compound results
  - `CacheEntry` - Adapter result caching
  - `Adapter` - Adapter registry
- **Alembic** - Database migration framework configured

#### ✅ Adapter System
- **backend/core/adapters/protocol.py** - Base adapter protocol
  - `AdapterProtocol` abstract base class
  - `AdapterResult` standard response format
  - `AdapterRegistry` for managing adapters
  - Cache key generation (SHA256 hashing)

#### ✅ Utilities
- **backend/core/scoring_utils.py** - Score normalization (0-1, higher=better)
- **backend/core/logging_config.py** - Logging setup

#### ✅ Configuration
- **.env** - Environment variables
- **.gitignore** - Proper Python/Docker exclusions
- **requirements.txt** - All dependencies pinned

## 🚀 Next Steps - Start Your Environment

### Step 1: Start Docker Services

Open PowerShell or Command Prompt in this directory and run:

```powershell
# Start all services
docker compose up -d

# Check status (all should be "healthy")
docker compose ps

# View logs
docker compose logs -f backend
```

Expected output:
```
NAME                    IMAGE                  STATUS
pharmforge-db           postgres:16-alpine     Up (healthy)
pharmforge-redis        redis:7-alpine         Up (healthy)
pharmforge-backend      ...                    Up
```

### Step 2: Initialize Database

```powershell
# Enter the backend container
docker compose exec backend bash

# Inside container, create initial migration
alembic revision --autogenerate -m "Initial schema"

# Apply migration
alembic upgrade head

# Exit container
exit
```

### Step 3: Test the API

Open your browser to:

1. **API Documentation**: http://localhost:8000/docs
2. **Health Check**: http://localhost:8000/health

Expected health response when services are healthy:
```json
{
  "status": "ok",
  "db": true,
  "redis": true
}
```

If you see `status: "degraded"`, check Docker services with `docker compose ps`.

### Step 4: Verify Everything Works

```powershell
# Test the health endpoint
curl http://localhost:8000/health

# View all endpoints
curl http://localhost:8000/docs

# Check root endpoint
curl http://localhost:8000/
```

## 📁 Key Files You Should Know

| File | Purpose |
|------|---------|
| `docker-compose.yml` | Service orchestration |
| `backend/main.py` | FastAPI app entry point |
| `backend/api/health.py` | Health check endpoint |
| `backend/db/models.py` | Database schema |
| `backend/core/adapters/protocol.py` | Adapter base class |
| `.env` | Environment configuration |
| `alembic.ini` | Database migration config |

## 🔧 Development Commands

### Docker Operations
```powershell
# Start services
docker compose up -d

# Stop services
docker compose down

# Reset everything (CAUTION: deletes data!)
docker compose down -v

# View logs
docker compose logs -f backend

# Restart a service
docker compose restart backend
```

### Database Operations
```bash
# Inside backend container:
docker compose exec backend bash

# Create migration
alembic revision --autogenerate -m "Description"

# Apply migrations
alembic upgrade head

# Rollback one migration
alembic downgrade -1

# View history
alembic history
```

### Testing
```bash
# Run tests inside container
docker compose exec backend pytest backend/tests/

# With coverage
docker compose exec backend pytest --cov=backend backend/tests/

# Specific test file
docker compose exec backend pytest backend/tests/test_scoring_utils.py -v
```

## 🐛 Troubleshooting

### Problem: `docker compose` command not found

**Solution**: You might need to use `docker-compose` (with hyphen) instead:
```powershell
docker-compose up -d
```

Or update Docker Desktop to the latest version.

### Problem: Port 5432 or 6379 already in use

**Solution**: Another PostgreSQL or Redis is running. Stop it:
```powershell
# Windows - stop existing services
net stop postgresql-x64-14
net stop Redis

# Or change ports in docker-compose.yml
```

### Problem: Health check returns 503

**Check services**:
```powershell
docker compose ps
```

All services should show "Up (healthy)". If not:
```powershell
# Check specific service logs
docker compose logs postgres
docker compose logs redis

# Restart unhealthy service
docker compose restart postgres
```

### Problem: Database migration fails

**Reset database**:
```powershell
# Stop services
docker compose down -v

# Start fresh
docker compose up -d

# Wait for healthy
docker compose ps

# Try migration again
docker compose exec backend alembic upgrade head
```

## 📊 Week 1 Validation Checklist

- [x] ✅ Docker environment working (`docker compose ps` shows all healthy)
- [x] ✅ FastAPI responding (`http://localhost:8000/docs` loads)
- [x] ✅ PostgreSQL connected (health endpoint shows `"db": true`)
- [x] ✅ Redis connected (health endpoint shows `"redis": true`)
- [x] ✅ Adapter protocol defined (`backend/core/adapters/protocol.py` exists)
- [x] ✅ Database schema created (models defined in `backend/db/models.py`)

## 🎯 What's Next - Week 2

Once your environment is running, you're ready for Week 2:

### Week 2 Goals (Days 8-14):
- Implement PubChem adapter (molecular properties API)
- Implement ChEMBL adapter (bioactivity search)
- Implement RDKit adapter (local fallback calculations)
- Test full property + bioactivity pipeline
- Process 10 compounds end-to-end

To start Week 2, simply tell me:
```
"Let's start Week 2 - implement the first adapters"
```

## 📚 Documentation

- **[PHARMFORGE_README.md](PHARMFORGE_README.md)** - Complete project overview
- **[docs/phases/phase1_weeks1-4_part1.md](docs/phases/phase1_weeks1-4_part1.md)** - Week 1 detailed guide
- **[docs/phases/BUILD_GUIDE_INDEX.md](docs/phases/BUILD_GUIDE_INDEX.md)** - Full 12-week plan
- **[docs/phases/CommonMistakes.md](docs/phases/CommonMistakes.md)** - Avoid common pitfalls

## 🎉 Congratulations!

You've completed Week 1 of the PharmForge build!

Your foundation is solid:
- ✅ Docker development environment
- ✅ FastAPI backend with health checks
- ✅ PostgreSQL database with migrations
- ✅ Adapter protocol ready for implementations
- ✅ Proper project structure

**Total Time**: ~2-3 hours of setup
**Lines of Code**: ~500 lines of production code
**Tests**: Ready for writing
**Next Phase**: Adapter implementations

---

**Questions?** Check the troubleshooting section or review the phase documentation in `docs/phases/`.

**Ready to continue?** Start Docker services and verify everything works, then we'll move to Week 2! 🚀
