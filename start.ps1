# PharmForge Startup Script
# Run this in PowerShell to start your development environment

Write-Host "🚀 Starting PharmForge Development Environment..." -ForegroundColor Cyan
Write-Host ""

# Check if Docker is running
Write-Host "Checking Docker..." -ForegroundColor Yellow
docker --version
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Docker not found. Please install Docker Desktop." -ForegroundColor Red
    exit 1
}
Write-Host "✅ Docker is installed" -ForegroundColor Green
Write-Host ""

# Start services
Write-Host "Starting Docker services..." -ForegroundColor Yellow
docker compose up -d
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Failed to start services" -ForegroundColor Red
    exit 1
}
Write-Host ""

# Wait for services to be healthy
Write-Host "Waiting for services to be healthy (15 seconds)..." -ForegroundColor Yellow
Start-Sleep -Seconds 15
Write-Host ""

# Check service status
Write-Host "Checking service status..." -ForegroundColor Yellow
docker compose ps
Write-Host ""

# Check health endpoint
Write-Host "Testing health endpoint..." -ForegroundColor Yellow
try {
    $response = Invoke-RestMethod -Uri "http://localhost:8000/health" -ErrorAction SilentlyContinue
    if ($response.status -eq "ok") {
        Write-Host "✅ Health check passed!" -ForegroundColor Green
        Write-Host "   Database: $($response.db)" -ForegroundColor Green
        Write-Host "   Redis: $($response.redis)" -ForegroundColor Green
    } else {
        Write-Host "⚠️  Health check returned degraded status" -ForegroundColor Yellow
        Write-Host $response
    }
} catch {
    Write-Host "⚠️  Backend not responding yet. It may still be starting up." -ForegroundColor Yellow
    Write-Host "   Try: curl http://localhost:8000/health" -ForegroundColor Yellow
}
Write-Host ""

Write-Host "📝 Next Steps:" -ForegroundColor Cyan
Write-Host "   1. Initialize database:" -ForegroundColor White
Write-Host "      docker compose exec backend alembic upgrade head" -ForegroundColor Gray
Write-Host ""
Write-Host "   2. Visit API docs:" -ForegroundColor White
Write-Host "      http://localhost:8000/docs" -ForegroundColor Gray
Write-Host ""
Write-Host "   3. View logs:" -ForegroundColor White
Write-Host "      docker compose logs -f backend" -ForegroundColor Gray
Write-Host ""
Write-Host "✅ PharmForge is ready!" -ForegroundColor Green
