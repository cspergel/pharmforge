@echo off
REM Test PubChem Adapter

echo.
echo ========================================
echo  Testing PubChem Adapter
echo ========================================
echo.

echo Running PubChem adapter test...
docker compose exec backend python /app/backend/test-pubchem.py

echo.
echo ========================================
echo  Test Complete!
echo ========================================
echo.
pause
