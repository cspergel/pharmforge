@echo off
echo Restarting backend container...
docker compose restart backend
echo.
echo Waiting for backend to start (5 seconds)...
timeout /t 5 /nobreak >nul
echo.
echo Testing health endpoint...
curl http://localhost:8000/health
echo.
echo.
pause
