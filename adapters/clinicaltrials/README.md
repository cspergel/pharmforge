# ClinicalTrials.gov Adapter

Comprehensive adapter for searching and analyzing clinical trials using the ClinicalTrials.gov API v2.

## Features

- **Trial Search**: Search by drug name, condition, NCT ID, phase, and status
- **Study Details**: Retrieve comprehensive trial information including:
  - Trial status and phases
  - Enrollment numbers
  - Start and completion dates
  - Primary outcomes
  - Interventions and conditions
- **Summary Statistics**: Aggregate data across multiple trials
- **Results Analysis**: Access published trial results when available

## API Documentation

- Base URL: https://clinicaltrials.gov/api/v2
- Documentation: https://clinicaltrials.gov/api/
- Rate Limit: 1 request/second (configured, be respectful)

## Installation

No additional dependencies beyond the base requirements.

## Usage

### Basic Search by Drug Name

```python
from adapters.clinicaltrials.adapter import ClinicalTrialsAdapter

adapter = ClinicalTrialsAdapter()

# Search for trials by drug name
result = await adapter.execute("pembrolizumab")

if result.success:
    print(f"Found {result.data['total_count']} trials")
    summary = result.data['summary']
    print(f"Trials with results: {summary['studies_with_results']}")
    print(f"By phase: {summary['by_phase']}")
```

### Advanced Search with Multiple Parameters

```python
# Search by drug and condition
result = await adapter.execute({
    "drug": "metformin",
    "condition": "diabetes",
    "phase": "Phase 3",
    "status": "Completed"
})

# Direct NCT ID lookup
result = await adapter.execute({
    "nct_id": "NCT12345678"
})
```

### Analyzing Trial Data

```python
result = await adapter.execute("aspirin")

if result.success:
    studies = result.data['studies']

    for study in studies[:5]:
        print(f"NCT ID: {study['nct_id']}")
        print(f"Title: {study['title']}")
        print(f"Status: {study['status']}")
        print(f"Phase: {study['phase']}")
        print(f"Enrollment: {study['enrollment']}")
        print(f"Primary Outcomes:")
        for outcome in study['primary_outcomes']:
            print(f"  - {outcome['measure']}")
```

## Input Parameters

### Simple String Input
- **drug_name** (str): Drug or intervention name to search

### Dictionary Input
- **drug** (str): Drug/intervention name
- **condition** (str): Medical condition or disease
- **nct_id** (str): NCT identifier for direct lookup
- **term** (str): General search term
- **phase** (str): Trial phase (e.g., "Phase 1", "Phase 2", "Phase 3")
- **status** (str): Overall status (e.g., "Completed", "Active", "Recruiting")

### Keyword Arguments
- **max_studies** (int): Maximum number of studies to return (default: 50)

## Output Format

```python
{
    "studies": [
        {
            "nct_id": "NCT12345678",
            "title": "Study Title",
            "status": "Completed",
            "phase": ["Phase 3"],
            "enrollment": 500,
            "start_date": "2020-01-01",
            "completion_date": "2023-01-01",
            "interventions": [
                {
                    "type": "Drug",
                    "name": "Drug Name",
                    "description": "Description"
                }
            ],
            "conditions": ["Condition 1", "Condition 2"],
            "primary_outcomes": [
                {
                    "measure": "Outcome measure",
                    "description": "Description",
                    "timeFrame": "24 weeks"
                }
            ],
            "has_results": true
        }
    ],
    "summary": {
        "total_studies": 100,
        "by_phase": {
            "Phase 1": 10,
            "Phase 2": 30,
            "Phase 3": 50,
            "Phase 4": 10
        },
        "by_status": {
            "Completed": 60,
            "Active": 20,
            "Recruiting": 20
        },
        "total_enrollment": 50000,
        "studies_with_results": 45
    },
    "total_count": 150
}
```

## Example Workflows

### 1. Drug Development Pipeline Analysis

```python
# Analyze all trials for a drug across phases
result = await adapter.execute("durvalumab")

summary = result.data['summary']
print("Development Pipeline:")
for phase in ['Phase 1', 'Phase 2', 'Phase 3', 'Phase 4']:
    count = summary['by_phase'].get(phase, 0)
    print(f"  {phase}: {count} trials")
```

### 2. Competitive Analysis

```python
# Compare clinical development for competing drugs
drugs = ["pembrolizumab", "nivolumab", "atezolizumab"]

for drug in drugs:
    result = await adapter.execute(drug)
    if result.success:
        summary = result.data['summary']
        print(f"{drug}:")
        print(f"  Total trials: {summary['total_studies']}")
        print(f"  Phase 3: {summary['by_phase'].get('Phase 3', 0)}")
        print(f"  With results: {summary['studies_with_results']}")
```

### 3. Indication Discovery

```python
# Find all conditions being tested for a drug
result = await adapter.execute("metformin")

conditions = set()
for study in result.data['studies']:
    conditions.update(study['conditions'])

print(f"Conditions being studied: {len(conditions)}")
for condition in sorted(conditions):
    print(f"  - {condition}")
```

## Error Handling

```python
result = await adapter.execute("nonexistent_drug")

if not result.success:
    print(f"Error: {result.error}")
else:
    if result.data['total_count'] == 0:
        print("No trials found")
```

## Rate Limiting

The adapter implements a 1-second delay between requests by default. Adjust in config:

```python
adapter = ClinicalTrialsAdapter()
adapter.config['rate_limit_delay'] = 2.0  # 2 seconds between requests
```

## Caching

Results are automatically cached using the AdapterProtocol caching mechanism:

```python
# First call - hits API
result1 = await adapter.execute("aspirin")  # cache_hit=False

# Second call - returns cached result
result2 = await adapter.execute("aspirin")  # cache_hit=True
```

## Notes

- The API returns a maximum of 1000 studies per search
- Use pagination for large result sets (not yet implemented)
- Some fields may be None if not available for a study
- Results section is only available for trials that have published results

## Troubleshooting

### No results found
- Try simpler search terms
- Check spelling of drug/condition names
- Use generic names instead of brand names

### Timeout errors
- Increase timeout in config: `adapter.config['timeout'] = 120`
- Check internet connection
- API may be temporarily unavailable

## Future Enhancements

- [ ] Pagination support for large result sets
- [ ] Advanced filtering options
- [ ] Statistical analysis of trial outcomes
- [ ] Export to common formats (CSV, JSON)
- [ ] Integration with other clinical databases
