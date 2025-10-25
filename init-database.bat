@echo off
REM PharmForge Database Initialization Script

echo.
echo ========================================
echo  PharmForge Database Initialization
echo ========================================
echo.

echo Step 1: Creating initial migration...
echo.
docker compose exec -w /app/backend backend alembic revision --autogenerate -m "Initial schema"
if errorlevel 1 (
    echo [ERROR] Failed to create migration
    pause
    exit /b 1
)
echo.
echo [OK] Migration created successfully
echo.

echo Step 2: Applying migration to database...
echo.
docker compose exec -w /app/backend backend alembic upgrade head
if errorlevel 1 (
    echo [ERROR] Failed to apply migration
    pause
    exit /b 1
)
echo.
echo [OK] Database initialized successfully
echo.

echo ========================================
echo  Database Setup Complete!
echo ========================================
echo.
echo Your database now has these tables:
echo   - pipeline_runs
echo   - compound_results
echo   - cache_entries
echo   - adapters
echo   - alembic_version
echo.
echo ========================================
echo  What's Next?
echo ========================================
echo.
echo 1. Visit API docs:
echo    http://localhost:8000/docs
echo.
echo 2. Test health endpoint:
echo    http://localhost:8000/health
echo.
echo 3. Ready for Week 2?
echo    Tell Claude: "Let's start Week 2"
echo.
echo ========================================
echo.
pause
