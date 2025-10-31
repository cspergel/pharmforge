# PharmForge Frontend - Streamlit MVP

This is the Streamlit web interface for PharmForge, providing a user-friendly UI for drug discovery workflows.

## Quick Start

### With Docker Compose (Recommended)
```bash
# From project root
docker-compose up -d frontend

# Access UI
http://localhost:8501
```

### Standalone Development
```bash
# Install dependencies
pip install streamlit==1.28.0 requests pandas plotly numpy

# Set API URL
export PHARMFORGE_API_URL=http://localhost:8000

# Run Streamlit
streamlit run streamlit_app.py
```

## Features

### Pages

1. **Home (🏠)**
   - Quick action cards for different input modes
   - Recent runs overview
   - System health monitoring

2. **New Run (🚀)**
   - Natural Language Query input
   - Batch Upload (CSV)
   - Evolution Mode (coming soon)
   - Advanced pipeline options

3. **My Runs (📊)**
   - Run listing with filters
   - Real-time progress tracking
   - Status indicators
   - Results viewer

4. **Adapters (🔌)**
   - System health dashboard
   - Adapter status monitoring
   - Configuration management

5. **Analytics (📈)**
   - Usage statistics
   - Popular targets
   - Performance metrics

### Design System

- **Primary Color:** #00bcd4 (Teal)
- **Font:** Inter (Google Fonts)
- **Components:** Custom styled buttons, inputs, cards
- **Responsive:** Streamlit's built-in grid system

## Configuration

### Environment Variables

- `PHARMFORGE_API_URL`: Backend API URL (default: http://backend:8000)

### Streamlit Config

Located in `.streamlit/config.toml`:
- Theme colors
- Server settings
- Port configuration (8501)

## Development

### File Structure
```
frontend/
├── streamlit_app.py      # Main application
├── .streamlit/
│   └── config.toml       # Streamlit configuration
└── README.md             # This file
```

### Adding New Pages

1. Create new function: `def show_new_page():`
2. Add route in `main()` function
3. Add navigation item in sidebar radio

### Customizing Design

Edit the CSS in the `st.markdown()` block at the top of `streamlit_app.py`:
- Colors: Update CSS variables
- Typography: Modify font sizes and weights
- Components: Add new component styles

## API Integration

### Backend Endpoints Used

- `GET /health` - System health check
- `POST /runs/` - Create new run
- `GET /runs/` - List all runs
- `GET /runs/{id}/results` - Fetch run results

### Making API Calls

```python
import requests

response = requests.get(f"{API_URL}/endpoint", timeout=5)
if response.status_code == 200:
    data = response.json()
```

## Testing

### Manual Testing
1. Start services: `docker-compose up -d`
2. Open http://localhost:8501
3. Test each page navigation
4. Try creating a run
5. Check adapter status

### Integration Testing
```bash
# Test frontend-backend connectivity
docker exec pharmforge-frontend python -c "
import requests
resp = requests.get('http://backend:8000/health')
print('Status:', resp.status_code)
"
```

## Troubleshooting

### Frontend Not Starting
```bash
# Check logs
docker-compose logs frontend

# Rebuild container
docker-compose up -d --build frontend
```

### Cannot Connect to Backend
```bash
# Verify backend is running
docker-compose ps backend

# Check network connectivity
docker exec pharmforge-frontend ping backend
```

### Port Already in Use
```bash
# Check what's using port 8501
netstat -ano | findstr :8501

# Change port in docker-compose.yml:
ports:
  - "8502:8501"  # External:Internal
```

## Migration to Next.js

This Streamlit MVP will eventually be replaced with a full Next.js application. The Next.js version will include:

- Command palette (⌘K)
- 3D molecule viewer
- Real-time WebSocket updates
- Dark mode
- Advanced visualizations
- Keyboard shortcuts
- Full component library

See `docs/PharmForge_Frontend_Design_Spec.md` for complete Next.js roadmap.

## Dependencies

Core:
- streamlit==1.28.0
- requests==2.31.0
- pandas==2.1.3
- plotly==5.18.0
- numpy==1.26.2

See Dockerfile.frontend for complete list.

## License

MIT License - See project root LICENSE file

## Support

For issues or questions:
1. Check FRONTEND_VALIDATION.md for common issues
2. Review docker-compose logs
3. Consult Phase 2 documentation in docs/phases/
