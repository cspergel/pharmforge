# PharmForge Frontend - Quick Start

## 30-Second Start

```bash
# From project root
cd C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2

# Start everything
docker-compose up -d

# Open browser
http://localhost:8501
```

That's it! Frontend is running and connected to backend.

---

## Manual Setup (No Docker)

```bash
# 1. Install Python 3.10+
python --version

# 2. Navigate to frontend
cd frontend

# 3. Install dependencies
pip install -r requirements.txt

# 4. Set API URL (Windows)
set PHARMFORGE_API_URL=http://localhost:8000

# 5. Run Streamlit
streamlit run streamlit_app.py
```

---

## Verify It's Working

### Check Health
1. Open http://localhost:8501
2. Sidebar should show "✅ System Online"
3. Should see "DB: healthy | Redis: healthy"

### Create Test Run
1. Click "🚀 New Run"
2. Enter: "Design EGFR inhibitors"
3. Click "🚀 Launch Pipeline"
4. Should see success message with run_id

### View Results
1. Click "📊 My Runs"
2. Find your run
3. Click "📊 View Results" (when completed)
4. Should see ranked compounds with visualizations

---

## Troubleshooting

### Problem: Can't connect to backend

```bash
# Check backend is running
docker-compose ps

# Restart backend
docker-compose restart backend

# Check logs
docker-compose logs backend
```

### Problem: Import errors

```bash
# Reinstall dependencies
pip install -r frontend/requirements.txt --force-reinstall

# Verify RDKit
python -c "from rdkit import Chem; print('RDKit OK')"
```

### Problem: Port 8501 in use

```bash
# Option 1: Kill process
netstat -ano | findstr :8501
taskkill /PID <PID> /F

# Option 2: Change port
streamlit run streamlit_app.py --server.port 8502
```

---

## Key Features

- **Natural Language Input:** Just describe what you want
- **Real-time Progress:** Watch pipeline execute step-by-step
- **Interactive Results:** Pareto plots, sortable tables, molecule viewer
- **Export:** Download results as CSV
- **Adapters:** 23 production-ready integrations

---

## Next Steps

👉 Read full guide: `FRONTEND_SETUP_GUIDE.md`
👉 API documentation: http://localhost:8000/docs
👉 Design spec: `docs/phases/PharmForge_Frontend_Design_Spec.md`

---

**Need Help?**
- Check `FRONTEND_SETUP_GUIDE.md` for detailed troubleshooting
- Review docker-compose logs: `docker-compose logs -f`
- Verify all services running: `docker-compose ps`
