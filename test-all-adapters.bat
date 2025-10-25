@echo off
REM Test All Adapters (PubChem + ChEMBL)

echo.
echo ========================================
echo  Testing All PharmForge Adapters
echo ========================================
echo.

echo Running adapter tests...
echo (This may take a minute due to API calls)
echo.
docker compose exec backend python /app/backend/test-adapters.py

echo.
echo ========================================
echo  Tests Complete!
echo ========================================
echo.
pause
