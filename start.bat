@echo off
REM PharmForge Startup Script for Windows Command Prompt

echo.
echo ========================================
echo  PharmForge Development Environment
echo ========================================
echo.

echo Checking Docker...
docker --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker not found. Please install Docker Desktop.
    pause
    exit /b 1
)
echo [OK] Docker is installed
echo.

echo Starting Docker services...
docker compose up -d
if errorlevel 1 (
    echo [ERROR] Failed to start services
    pause
    exit /b 1
)
echo.

echo Waiting for services to start (15 seconds)...
timeout /t 15 /nobreak >nul
echo.

echo Checking service status...
docker compose ps
echo.

echo Testing health endpoint...
curl -s http://localhost:8000/health
echo.
echo.

echo ========================================
echo  Next Steps:
echo ========================================
echo.
echo 1. Initialize database:
echo    docker compose exec backend alembic upgrade head
echo.
echo 2. Visit API docs:
echo    http://localhost:8000/docs
echo.
echo 3. View logs:
echo    docker compose logs -f backend
echo.
echo ========================================
echo.
pause
