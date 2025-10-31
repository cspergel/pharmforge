# PharmForge Deployment Guide

**Version:** 0.3.0
**Last Updated:** October 26, 2025
**Status:** Phase 3 Implementation

---

## Table of Contents

1. [Local Development Setup](#local-development-setup)
2. [Docker Compose Deployment](#docker-compose-deployment)
3. [Environment Variables](#environment-variables)
4. [GPU Configuration](#gpu-configuration)
5. [AWS Cloud Deployment](#aws-cloud-deployment)
6. [Production Checklist](#production-checklist)
7. [Monitoring & Health Checks](#monitoring--health-checks)
8. [Troubleshooting](#troubleshooting)

---

## Local Development Setup

### Prerequisites

- **Docker** 24.0+ and **Docker Compose** 2.20+
- **Python** 3.9+ (for local development)
- **Git** for version control
- **NVIDIA GPU** (optional, for GPU-accelerated adapters)
- **16GB+ RAM** recommended

### Quick Start

```bash
# Clone the repository
git clone https://github.com/your-org/pharmforge.git
cd pharmforge/claude-code-agents-wizard-v2

# Copy environment template
cp .env.example .env

# Edit .env with your settings
nano .env

# Start all services
docker-compose up -d

# Check logs
docker-compose logs -f

# Verify services are healthy
curl http://localhost:8000/health
```

### Service URLs

- **Backend API:** http://localhost:8000
- **API Documentation:** http://localhost:8000/docs
- **Frontend UI:** http://localhost:8501
- **Redis:** localhost:6379
- **PostgreSQL:** localhost:5432

---

## Docker Compose Deployment

### Architecture

PharmForge uses a multi-container architecture:

```
┌─────────────────────────────────────────────────────────┐
│                    Docker Compose                        │
├─────────────────────────────────────────────────────────┤
│  ┌──────────┐  ┌──────────┐  ┌──────────┐             │
│  │ Backend  │  │ Frontend │  │  Celery  │             │
│  │ (FastAPI)│  │(Streamlit)│  │ Worker   │             │
│  │  :8000   │  │  :8501   │  │          │             │
│  └────┬─────┘  └──────────┘  └────┬─────┘             │
│       │                            │                    │
│  ┌────▼────────────────────────────▼─────┐             │
│  │         Redis (Cache + Queue)          │             │
│  │              :6379                     │             │
│  └────────────────────────────────────────┘             │
│                                                          │
│  ┌────────────────────────────────────────┐             │
│  │     PostgreSQL (Database)              │             │
│  │              :5432                     │             │
│  └────────────────────────────────────────┘             │
└─────────────────────────────────────────────────────────┘
```

### docker-compose.yml Configuration

```yaml
version: '3.8'

services:
  # PostgreSQL Database
  db:
    image: postgres:15-alpine
    container_name: pharmforge-db
    environment:
      POSTGRES_DB: ${POSTGRES_DB:-pharmforge}
      POSTGRES_USER: ${POSTGRES_USER:-pharmforge}
      POSTGRES_PASSWORD: ${POSTGRES_PASSWORD:-changeme}
    volumes:
      - postgres_data:/var/lib/postgresql/data
    ports:
      - "5432:5432"
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U pharmforge"]
      interval: 10s
      timeout: 5s
      retries: 5

  # Redis Cache + Message Queue
  redis:
    image: redis:7-alpine
    container_name: pharmforge-redis
    command: redis-server --appendonly yes
    volumes:
      - redis_data:/data
    ports:
      - "6379:6379"
    healthcheck:
      test: ["CMD", "redis-cli", "ping"]
      interval: 10s
      timeout: 5s
      retries: 5

  # FastAPI Backend
  backend:
    build:
      context: .
      dockerfile: backend/Dockerfile
    container_name: pharmforge-backend
    environment:
      - DATABASE_URL=postgresql://${POSTGRES_USER}:${POSTGRES_PASSWORD}@db:5432/${POSTGRES_DB}
      - REDIS_URL=redis://redis:6379/0
      - OPENAI_API_KEY=${OPENAI_API_KEY}
      - BIOGRID_ACCESS_KEY=${BIOGRID_ACCESS_KEY}
      - ENV=development
    volumes:
      - ./backend:/app/backend
      - ./adapters:/app/adapters
      - ./cache:/app/cache
      - ./models:/app/models
    ports:
      - "8000:8000"
    depends_on:
      db:
        condition: service_healthy
      redis:
        condition: service_healthy
    command: uvicorn backend.main:app --host 0.0.0.0 --port 8000 --reload

  # Celery Worker
  worker:
    build:
      context: .
      dockerfile: backend/Dockerfile
    container_name: pharmforge-worker
    environment:
      - DATABASE_URL=postgresql://${POSTGRES_USER}:${POSTGRES_PASSWORD}@db:5432/${POSTGRES_DB}
      - REDIS_URL=redis://redis:6379/0
      - OPENAI_API_KEY=${OPENAI_API_KEY}
      - BIOGRID_ACCESS_KEY=${BIOGRID_ACCESS_KEY}
    volumes:
      - ./backend:/app/backend
      - ./adapters:/app/adapters
      - ./cache:/app/cache
      - ./models:/app/models
    depends_on:
      - db
      - redis
    command: celery -A backend.celery_app worker --loglevel=info

  # Streamlit Frontend
  frontend:
    build:
      context: ./frontend
      dockerfile: Dockerfile
    container_name: pharmforge-frontend
    environment:
      - PHARMFORGE_API_URL=http://backend:8000
    volumes:
      - ./frontend:/app
    ports:
      - "8501:8501"
    depends_on:
      - backend
    command: streamlit run app.py --server.port=8501 --server.address=0.0.0.0

volumes:
  postgres_data:
  redis_data:
```

### Starting Services

```bash
# Start all services in detached mode
docker-compose up -d

# Start specific service
docker-compose up -d backend

# View logs
docker-compose logs -f backend

# Stop all services
docker-compose down

# Stop and remove volumes (WARNING: deletes data)
docker-compose down -v
```

### Scaling Workers

```bash
# Scale Celery workers to 4 instances
docker-compose up -d --scale worker=4

# Check worker status
docker-compose exec worker celery -A backend.celery_app inspect active
```

---

## Environment Variables

### .env File Template

Create a `.env` file in the project root:

```bash
# Database Configuration
POSTGRES_DB=pharmforge
POSTGRES_USER=pharmforge
POSTGRES_PASSWORD=your_secure_password_here

# Redis Configuration
REDIS_URL=redis://redis:6379/0

# API Keys (Optional - most adapters work without)
OPENAI_API_KEY=sk-your-openai-key-here
BIOGRID_ACCESS_KEY=your-biogrid-key-here
GOOGLE_CSE_ID=your-google-cse-id-here

# Application Settings
ENV=development
LOG_LEVEL=INFO
CACHE_TTL=86400

# GPU Support (if available)
CUDA_VISIBLE_DEVICES=0
```

### Required vs Optional Keys

| Variable | Required? | Purpose | How to Get |
|----------|-----------|---------|------------|
| `POSTGRES_PASSWORD` | ✅ Yes | Database security | Set your own |
| `OPENAI_API_KEY` | ⚠️ Optional | LLM retrosynthesis | [OpenAI Platform](https://platform.openai.com/) |
| `BIOGRID_ACCESS_KEY` | ⚠️ Optional | Protein interactions | [BioGRID Registration](https://webservice.thebiogrid.org/) (free) |
| `GOOGLE_CSE_ID` | ⚠️ Optional | Google Patents search | [Google CSE](https://programmablesearchengine.google.com/) |

**Note:** 36 out of 39 adapters work without any API keys!

---

## GPU Configuration

### Enabling GPU Support

PharmForge supports GPU acceleration for:
- **AutoDock Vina** (docking)
- **GNINA** (scoring)
- **DiffDock** (ML docking)
- **OpenMM** (MD simulations)

### Docker GPU Setup

1. **Install NVIDIA Container Toolkit:**

```bash
# Ubuntu/Debian
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/nvidia-docker/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/nvidia-docker/$distribution/nvidia-docker.list | \
  sudo tee /etc/apt/sources.list.d/nvidia-docker.list

sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

2. **Update docker-compose.yml:**

```yaml
services:
  backend:
    # ... existing config ...
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
```

3. **Verify GPU Access:**

```bash
docker-compose exec backend nvidia-smi
```

### GPU Requirements

| Adapter | Min VRAM | Recommended |
|---------|----------|-------------|
| Vina | N/A (CPU-only fallback) | Any NVIDIA GPU |
| GNINA | 4GB | 8GB |
| DiffDock | 8GB | 12GB+ |
| OpenMM | 2GB | 8GB |

---

## AWS Cloud Deployment

### Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│                     AWS Cloud                            │
├─────────────────────────────────────────────────────────┤
│                                                           │
│  ┌────────────────────────────────────────────┐         │
│  │   Application Load Balancer (ALB)          │         │
│  │   pharmforge-alb.us-east-1.elb.amazonaws...│         │
│  └──────────────┬──────────────┬───────────────┘         │
│                 │              │                          │
│  ┌──────────────▼────┐  ┌──────▼──────────────┐         │
│  │  ECS Service      │  │  ECS Service         │         │
│  │  (Backend)        │  │  (Worker)            │         │
│  │  Fargate Tasks    │  │  Fargate Tasks       │         │
│  └──────────┬────────┘  └──────┬───────────────┘         │
│             │                  │                          │
│  ┌──────────▼──────────────────▼───────────────┐         │
│  │      RDS PostgreSQL (Managed)                │         │
│  │      pharmforge-db.xxxxx.us-east-1.rds...   │         │
│  └──────────────────────────────────────────────┘         │
│                                                           │
│  ┌──────────────────────────────────────────────┐         │
│  │   ElastiCache Redis (Managed)                │         │
│  │   pharmforge-redis.xxxxx.cache.amazonaws...  │         │
│  └──────────────────────────────────────────────┘         │
│                                                           │
│  ┌──────────────────────────────────────────────┐         │
│  │   S3 Bucket (Results Storage)                │         │
│  │   pharmforge-results                         │         │
│  └──────────────────────────────────────────────┘         │
│                                                           │
│  ┌──────────────────────────────────────────────┐         │
│  │   CloudWatch (Monitoring & Logs)             │         │
│  └──────────────────────────────────────────────┘         │
└─────────────────────────────────────────────────────────┘
```

### Prerequisites

- AWS account with admin access
- AWS CLI configured (`aws configure`)
- Terraform 1.5+ installed
- Docker images built and tagged

### Step 1: Terraform Setup

```bash
cd infra/terraform

# Initialize Terraform
terraform init

# Review planned changes
terraform plan

# Apply infrastructure
terraform apply

# Note the outputs (ALB DNS, ECR URLs)
terraform output
```

### Step 2: Build and Push Docker Images

```bash
# Get ECR login
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin <ECR_URL>

# Build backend image
docker build -t pharmforge/backend:latest -f backend/Dockerfile .

# Tag for ECR
docker tag pharmforge/backend:latest <ECR_URL>/pharmforge/backend:latest

# Push to ECR
docker push <ECR_URL>/pharmforge/backend:latest

# Repeat for worker and frontend images
```

### Step 3: Update ECS Service

```bash
# Force new deployment with updated images
aws ecs update-service \
  --cluster pharmforge-cluster \
  --service pharmforge-backend \
  --force-new-deployment \
  --region us-east-1

# Monitor deployment
aws ecs describe-services \
  --cluster pharmforge-cluster \
  --services pharmforge-backend \
  --region us-east-1
```

### Step 4: Configure Secrets

```bash
# Store secrets in AWS Secrets Manager
aws secretsmanager create-secret \
  --name pharmforge/openai_key \
  --secret-string "sk-your-key-here" \
  --region us-east-1

# Update ECS task definition to reference secrets
# (See Terraform configuration for details)
```

### Infrastructure Costs (Estimated)

| Service | Configuration | Monthly Cost |
|---------|--------------|--------------|
| ECS Fargate (Backend) | 0.25 vCPU, 0.5GB RAM, 2 tasks | ~$15 |
| ECS Fargate (Worker) | 0.5 vCPU, 1GB RAM, 2 tasks | ~$30 |
| RDS PostgreSQL | db.t3.micro, 20GB | ~$15 |
| ElastiCache Redis | cache.t3.micro | ~$12 |
| ALB | 1 load balancer | ~$20 |
| S3 | 100GB storage, minimal transfer | ~$2 |
| **Total** | | **~$94/month** |

**Free Tier:** First year eligible for significant discounts (~50% reduction).

---

## Production Checklist

### Security

- [ ] Change default passwords in `.env`
- [ ] Enable SSL/TLS (ALB HTTPS listener)
- [ ] Configure AWS WAF for API protection
- [ ] Set up IAM roles with least privilege
- [ ] Enable CloudTrail audit logging
- [ ] Rotate API keys regularly
- [ ] Configure VPC security groups (restrict ports)
- [ ] Enable RDS encryption at rest
- [ ] Use AWS Secrets Manager (not .env files)

### Performance

- [ ] Configure RDS connection pooling (max 20 connections)
- [ ] Set Redis maxmemory policy (allkeys-lru)
- [ ] Enable CloudFront CDN for frontend
- [ ] Configure auto-scaling policies (CPU > 70%)
- [ ] Set up read replicas for RDS (if needed)
- [ ] Enable ECS task auto-scaling
- [ ] Configure ALB health checks (30s interval)

### Monitoring

- [ ] Set up CloudWatch dashboards
- [ ] Configure alarms (CPU, memory, errors)
- [ ] Enable X-Ray tracing
- [ ] Set up SNS notifications for alerts
- [ ] Configure log retention (7-30 days)
- [ ] Integrate Sentry for error tracking
- [ ] Set up uptime monitoring (Pingdom/UptimeRobot)

### Backup & DR

- [ ] Enable RDS automated backups (7-day retention)
- [ ] Configure S3 versioning
- [ ] Set up cross-region replication (optional)
- [ ] Document restore procedures
- [ ] Test backup restoration monthly
- [ ] Export critical data to cold storage

### Compliance

- [ ] Review data residency requirements
- [ ] Implement GDPR compliance (if applicable)
- [ ] Configure audit logging
- [ ] Document data retention policies
- [ ] Set up disaster recovery plan

---

## Monitoring & Health Checks

### Health Endpoint

```bash
# Check service health
curl http://localhost:8000/health

# Example response (healthy)
{
  "status": "ok",
  "checks": {
    "db": "ok",
    "redis": "ok",
    "adapters": "ok"
  }
}

# Example response (degraded)
{
  "status": "degraded",
  "checks": {
    "db": "ok",
    "redis": "failed",
    "adapters": "ok"
  }
}
```

**HTTP Status Codes:**
- `200` - All systems healthy
- `503` - One or more systems degraded

### CloudWatch Metrics

Key metrics to monitor:

| Metric | Threshold | Action |
|--------|-----------|--------|
| CPU Utilization | > 80% | Scale up tasks |
| Memory Utilization | > 85% | Scale up memory |
| DB Connections | > 15/20 | Investigate leaks |
| 5XX Errors | > 1% | Check logs |
| Response Time | > 2s | Optimize queries |
| Cache Hit Rate | < 50% | Review cache TTL |

### Example CloudWatch Dashboard

```json
{
  "widgets": [
    {
      "type": "metric",
      "properties": {
        "metrics": [
          ["AWS/ECS", "CPUUtilization", {"stat": "Average"}],
          [".", "MemoryUtilization", {"stat": "Average"}]
        ],
        "period": 300,
        "stat": "Average",
        "region": "us-east-1",
        "title": "ECS Resource Utilization"
      }
    }
  ]
}
```

### Logging

```bash
# View backend logs
docker-compose logs -f backend

# View worker logs
docker-compose logs -f worker

# View all logs
docker-compose logs -f

# AWS CloudWatch logs
aws logs tail /ecs/pharmforge-backend --follow
```

---

## Troubleshooting

### Service Won't Start

**Symptom:** Container exits immediately

```bash
# Check logs
docker-compose logs backend

# Common issues:
# 1. Database not ready
#    Solution: Wait for db health check
# 2. Missing environment variables
#    Solution: Check .env file
# 3. Port conflict
#    Solution: Change ports in docker-compose.yml
```

### Database Connection Errors

**Symptom:** `FATAL: password authentication failed`

```bash
# Reset database
docker-compose down -v
docker-compose up -d db

# Wait for health check
docker-compose ps

# Run migrations
docker-compose exec backend alembic upgrade head
```

### Redis Connection Errors

**Symptom:** `Error 111 connecting to redis:6379. Connection refused`

```bash
# Check Redis status
docker-compose ps redis

# Restart Redis
docker-compose restart redis

# Test connection
docker-compose exec redis redis-cli ping
```

### GPU Not Detected

**Symptom:** `CUDA not available` errors

```bash
# Check NVIDIA driver
nvidia-smi

# Check Docker GPU support
docker run --rm --gpus all nvidia/cuda:11.8.0-base-ubuntu22.04 nvidia-smi

# Verify container GPU access
docker-compose exec backend nvidia-smi
```

### High Memory Usage

**Symptom:** Container OOM killed

```bash
# Check resource usage
docker stats

# Increase memory limits in docker-compose.yml
services:
  backend:
    deploy:
      resources:
        limits:
          memory: 4G

# Restart services
docker-compose up -d
```

### Slow Pipeline Execution

**Symptom:** Pipelines taking > 5 minutes

```bash
# Check cache hit rate
curl http://localhost:8000/api/v1/stats/cache

# Clear cache if corrupted
docker-compose exec redis redis-cli FLUSHALL

# Check worker status
docker-compose exec worker celery -A backend.celery_app inspect active

# Scale workers
docker-compose up -d --scale worker=4
```

### AWS Deployment Issues

**Symptom:** ECS tasks failing to start

```bash
# Check task logs
aws ecs describe-tasks \
  --cluster pharmforge-cluster \
  --tasks <task-id> \
  --region us-east-1

# Check ECR image availability
aws ecr describe-images \
  --repository-name pharmforge/backend \
  --region us-east-1

# Verify IAM permissions
aws sts get-caller-identity
```

---

## Support & Resources

### Documentation
- [User Guide](USER_GUIDE.md)
- [API Documentation](http://localhost:8000/docs)
- [Phase 3 Plan](PHASE3_IMPLEMENTATION_PLAN.md)

### Community
- GitHub Issues: https://github.com/your-org/pharmforge/issues
- Email: support@pharmforge.org

### Emergency Contacts
- Database issues: Check PostgreSQL logs first
- API issues: Check backend logs and health endpoint
- Deployment issues: Review Terraform state and AWS console

---

**Version:** 0.3.0
**Last Updated:** October 26, 2025
**Maintained by:** PharmForge Team
