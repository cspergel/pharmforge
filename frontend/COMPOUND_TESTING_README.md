# Compound Testing Page

## Overview

The Compound Testing page is a critical feature of PharmForge that allows users to quickly test individual compounds across multiple adapters. It provides a simple, streamlined interface for validating SMILES strings and executing adapter workflows.

## Location

- **File:** `frontend/pages/compound_testing.py`
- **Navigation:** Accessible via sidebar "🧪 Compound Testing" or home page "Test Compound" button

## Features

### 1. SMILES Input Section

- **Text input** for SMILES string entry
- **Validate button** to check SMILES validity using RDKit
- **Example button** for quick testing with Aspirin (CC(=O)Oc1ccccc1C(=O)O)
- **2D molecule preview** when SMILES is validated
- **Example library** with 6 common drug molecules:
  - Caffeine
  - Ibuprofen
  - Penicillin
  - Metformin
  - Atorvastatin
  - Sildenafil

### 2. Adapter Selection

- **Categorized adapters** grouped by type:
  - Docking (AutoDock Vina, GNINA)
  - ADMET (TDC ADMET, pkCSM)
  - Retrosynthesis (AiZynthFinder, LLM Retrosynthesis)
  - Novelty (ChEMBL Similarity)
  - Database (ZINC, DrugCentral)

- **Select All / Deselect All** buttons per category
- **Status indicators** for each adapter (✅ healthy, ⚠️ degraded, ❌ offline)
- **Info buttons** showing adapter descriptions
- **Checkboxes** for individual adapter selection

### 3. Parameter Customization (Optional)

- **Collapsible section** for advanced users
- **Parameter inputs** for each selected adapter
- **Default values** pre-filled
- **Type-specific inputs**:
  - Integer sliders with min/max
  - Float inputs with precision
  - Select dropdowns for choices
  - Text inputs for strings

### 4. Execute Button

- **Smart validation**: Disabled until SMILES is valid and adapters selected
- **Progress tracking**: Shows which adapter is currently running
- **Real-time updates**: Progress bar and status text
- **Error handling**: Gracefully handles API failures with mock data fallback

### 5. Results Display

#### Summary Section
- Total adapters executed
- Successful vs failed count
- Total execution time

#### Results Table
- Adapter name
- Status with emoji
- Execution time
- Preview of results

#### Detailed Results
- Expandable cards for each adapter
- Full JSON output for successful runs
- Error messages for failed runs
- Execution time per adapter

#### Export Options
- **CSV export**: Tabular format for spreadsheet analysis
- **JSON export**: Full data structure with metadata
- Timestamped filenames

## API Integration

### Backend Endpoints Used

1. **List Adapters**: `GET /api/v1/adapters`
   - Fetches available adapters and their metadata
   - Fallback to mock data if unavailable

2. **Test Adapter**: `POST /api/v1/adapters/test`
   - Executes single adapter on SMILES string
   - Returns adapter-specific results

### Mock Data Support

For development without backend connection, the page includes:
- Mock adapter list with 9 common adapters
- Mock result generation based on adapter type
- Realistic execution time simulation

## User Experience

### Workflow
1. Enter SMILES string (or use example)
2. Validate SMILES
3. View 2D molecule structure
4. Select adapters to run
5. (Optional) Customize parameters
6. Click "Run Selected Adapters"
7. View results in multiple formats
8. Export data for further analysis

### Error Handling
- Invalid SMILES detection
- Connection error fallback
- Adapter timeout handling
- Missing parameter validation

### Visual Feedback
- Progress bar during execution
- Status messages
- Color-coded results (✅ success, ❌ error)
- Loading spinners

## Code Structure

### Main Function
```python
render_compound_testing_page()
```
Main entry point that orchestrates the entire page.

### Helper Functions

- `validate_smiles(smiles)`: RDKit-based SMILES validation
- `group_adapters_by_category(adapters)`: Categorize adapters
- `get_status_icon(status)`: Map status to emoji
- `format_execution_time(seconds)`: Human-readable time format
- `export_results_csv(results)`: CSV export
- `export_results_json(results)`: JSON export

### Execution Functions

- `execute_adapters(...)`: Main execution loop with progress tracking
- `generate_mock_result(...)`: Mock data for development
- `display_results(...)`: Results rendering

## Session State

The page uses Streamlit session state for:
- `smiles`: Current SMILES string
- `selected_adapters`: Set of selected adapter names
- `results`: Execution results
- `smiles_valid`: Validation status

## Dependencies

- `streamlit`: UI framework
- `rdkit`: SMILES validation and molecule rendering
- `pandas`: Data manipulation
- `components.api_client`: Backend communication
- `components.molecule_viewer`: 2D structure display

## Usage Example

```python
# In streamlit_app.py
from pages.compound_testing import render_compound_testing_page

# Add to navigation
page = st.radio("Go to", [..., "🧪 Compound Testing", ...])

# Route to page
if page == "🧪 Compound Testing":
    render_compound_testing_page()
```

## Future Enhancements

1. **Batch compound testing**: Test multiple SMILES at once
2. **Result comparison**: Compare results across compounds
3. **Visualization**: Charts and plots for adapter results
4. **History**: Save and recall previous tests
5. **Favorites**: Save commonly used adapter combinations
6. **Real-time updates**: WebSocket for long-running adapters
7. **3D viewer**: Interactive 3D molecule structures
8. **Parameter presets**: Save/load parameter configurations

## Testing

To test the page without backend:
1. Run Streamlit app: `streamlit run streamlit_app.py`
2. Navigate to "🧪 Compound Testing"
3. Enter SMILES or use example
4. Select adapters
5. Click "Run Selected Adapters"
6. Mock data will be displayed

## Performance

- Page loads in <1s
- SMILES validation: <100ms
- Mock adapter execution: ~500ms per adapter
- Real adapter execution: varies (2-30s typical)
- Results display: <500ms

## Accessibility

- Clear labels for all inputs
- Help text and tooltips
- Keyboard navigation support
- Color-blind friendly icons (emoji + text)
- Mobile-responsive layout

## Security

- SMILES validation prevents code injection
- API client handles authentication
- No sensitive data stored in session state
- Export files contain only results data
