# FDA FAERS Adapter

Comprehensive adapter for analyzing adverse drug events using the FDA Adverse Event Reporting System (FAERS) via openFDA API.

## Features

- **Adverse Event Search**: Query adverse events by drug name or active ingredient
- **Safety Metrics**: Calculate pharmacovigilance metrics including:
  - Reporting frequencies
  - Proportional Reporting Ratio (PRR)
  - Top adverse reactions
- **Demographics Analysis**: Patient demographics and serious event analysis
- **Signal Detection**: Identify potential safety signals

## API Documentation

- Base URL: https://api.fda.gov/drug/event.json
- Documentation: https://open.fda.gov/apis/drug/event/
- Rate Limit: 240 requests/minute (4 per second)

## Installation

No additional dependencies beyond the base requirements.

## Usage

### Basic Adverse Event Search

```python
from adapters.fda_faers.adapter import FDAFAERSAdapter

adapter = FDAFAERSAdapter()

# Search for adverse events
result = await adapter.execute("aspirin")

if result.success:
    metrics = result.data['safety_metrics']
    print(f"Total events: {metrics['total_events']}")
    print(f"Unique reactions: {metrics['unique_reactions']}")

    # Top adverse reactions
    for reaction in metrics['top_reactions'][:5]:
        print(f"{reaction['term']}: {reaction['count']} ({reaction['percentage']:.1f}%)")
```

### Demographics Analysis

```python
result = await adapter.execute("warfarin", include_demographics=True)

if result.success:
    demographics = result.data['demographics']

    print(f"Serious events: {demographics['serious_events']}")
    print(f"Serious percentage: {demographics['serious_percentage']:.1f}%")

    print("\nSex distribution:")
    for sex, count in demographics['sex_distribution'].items():
        print(f"  {sex}: {count}")

    print("\nAge distribution:")
    for age_group, count in demographics['age_distribution'].items():
        print(f"  {age_group}: {count}")
```

### Safety Signal Detection

```python
result = await adapter.execute("drug_name")

if result.success:
    metrics = result.data['safety_metrics']

    # Check for signals using PRR
    if 'prr_signals' in metrics:
        print("Potential Safety Signals (PRR > 2.0):")
        for reaction, prr in metrics['prr_signals'].items():
            if prr > 2.0:
                print(f"  {reaction}: PRR = {prr}")
```

## Input Parameters

### Simple String Input
- **drug_name** (str): Drug name or active ingredient

### Dictionary Input
- **drug** (str): Drug name
- **ingredient** (str): Active ingredient name

### Keyword Arguments
- **include_demographics** (bool): Include demographic analysis (default: True)
- **calculate_metrics** (bool): Calculate safety metrics (default: True)
- **search_field** (str): Field to search (default: "patient.drug.medicinalproduct")

## Output Format

```python
{
    "events": [
        {
            "report_id": "12345678",
            "receive_date": "20230101",
            "serious": 1,
            "age": "65",
            "age_unit": "801",  # Years
            "sex": "2",  # 1=Male, 2=Female
            "weight_kg": 70.0,
            "reactions": [
                {
                    "term": "Nausea",
                    "outcome": "1"  # 1=Recovered
                }
            ],
            "drugs": [
                {
                    "name": "ASPIRIN",
                    "role": "1",  # 1=Suspect, 2=Concomitant, 3=Interacting
                    "indication": "PAIN"
                }
            ],
            "outcomes": [1, 0, 0, 0]  # Death, Hospitalization, Life-threatening, Disabling
        }
    ],
    "total_count": 50000,
    "safety_metrics": {
        "total_events": 50000,
        "unique_reactions": 150,
        "top_reactions": [
            {
                "term": "Nausea",
                "count": 5000,
                "percentage": 10.0
            }
        ],
        "prr_signals": {
            "Serious Reaction": 2.5
        }
    },
    "demographics": {
        "sex_distribution": {
            "Male": 20000,
            "Female": 28000,
            "Unknown": 2000
        },
        "age_distribution": {
            "<18": 1000,
            "18-44": 10000,
            "45-64": 20000,
            "65+": 15000,
            "Unknown": 4000
        },
        "serious_events": 25000,
        "serious_percentage": 50.0
    }
}
```

## Pharmacovigilance Metrics

### Proportional Reporting Ratio (PRR)

The PRR is calculated as:

```
PRR = (a/b) / (c/d)

where:
a = number of reports with drug AND reaction
b = number of reports with drug
c = number of reports with reaction (excluding drug)
d = total number of reports

Signal threshold: PRR ≥ 2.0
```

### Reporting Odds Ratio (ROR)

```python
# Calculate ROR for a specific drug-reaction pair
def calculate_ror(a, b, c, d):
    """
    a = drug+reaction
    b = drug only
    c = reaction only
    d = neither
    """
    ror = (a * d) / (b * c)
    return ror
```

## Example Workflows

### 1. Comparative Safety Analysis

```python
adapter = FDAFAERSAdapter()

drugs = ["ibuprofen", "naproxen", "celecoxib"]

for drug in drugs:
    result = await adapter.execute(drug)
    if result.success:
        metrics = result.data['safety_metrics']
        demographics = result.data['demographics']

        print(f"\n{drug}:")
        print(f"  Total events: {metrics['total_events']:,}")
        print(f"  Serious events: {demographics['serious_percentage']:.1f}%")
        print(f"  Top reaction: {metrics['top_reactions'][0]['term']}")
```

### 2. Signal Detection Workflow

```python
# Screen for safety signals
result = await adapter.execute("new_drug")

if result.success:
    metrics = result.data['safety_metrics']

    # Apply signal detection criteria
    # - N ≥ 3 cases
    # - PRR ≥ 2.0
    # - Chi-square ≥ 4

    signals = []
    for reaction in metrics['top_reactions']:
        if reaction['count'] >= 3 and reaction['percentage'] > 5.0:
            signals.append(reaction)

    if signals:
        print(f"Potential signals detected: {len(signals)}")
        for signal in signals:
            print(f"  - {signal['term']}: {signal['count']} cases")
```

### 3. Age-Specific Risk Assessment

```python
result = await adapter.execute("medication", include_demographics=True)

if result.success:
    demographics = result.data['demographics']
    age_dist = demographics['age_distribution']

    total = sum(age_dist.values())

    print("Age-specific risk:")
    for age_group, count in age_dist.items():
        percentage = (count / total * 100) if total > 0 else 0
        print(f"  {age_group}: {count} ({percentage:.1f}%)")

    # Check for pediatric or geriatric concerns
    pediatric_pct = (age_dist.get('<18', 0) / total * 100) if total > 0 else 0
    geriatric_pct = (age_dist.get('65+', 0) / total * 100) if total > 0 else 0

    if pediatric_pct > 20:
        print(f"\n⚠ HIGH pediatric exposure: {pediatric_pct:.1f}%")
    if geriatric_pct > 40:
        print(f"\n⚠ HIGH geriatric exposure: {geriatric_pct:.1f}%")
```

### 4. Drug-Drug Interaction Analysis

```python
result = await adapter.execute("warfarin")

if result.success:
    events = result.data['events']

    # Analyze concomitant medications
    concomitant_drugs = {}

    for event in events:
        drugs = event.get('drugs', [])
        for drug in drugs:
            if drug.get('role') == '2':  # Concomitant
                name = drug.get('name')
                if name:
                    concomitant_drugs[name] = concomitant_drugs.get(name, 0) + 1

    print("Most common concomitant medications:")
    for drug, count in sorted(concomitant_drugs.items(), key=lambda x: x[1], reverse=True)[:10]:
        print(f"  {drug}: {count}")
```

## Understanding FAERS Data

### Serious Outcomes

Events are considered serious if they result in:
1. Death
2. Hospitalization (initial or prolonged)
3. Life-threatening situation
4. Disability or permanent damage
5. Congenital anomaly
6. Required intervention to prevent permanent damage

### Limitations

- **Underreporting**: Not all adverse events are reported
- **Reporting bias**: More likely to report serious or unusual events
- **Causality**: Reports don't prove causation
- **Data quality**: Self-reported data may have inaccuracies
- **Duplicate reports**: Same event may be reported multiple times

## Error Handling

```python
result = await adapter.execute("drug_name")

if not result.success:
    print(f"Error: {result.error}")
elif result.data['total_count'] == 0:
    print("No adverse events found")
else:
    # Process data
    pass
```

## Rate Limiting

The adapter implements automatic rate limiting (240 requests/minute):

```python
adapter = FDAFAERSAdapter()
adapter.config['rate_limit_delay'] = 0.3  # 300ms between requests
```

## Best Practices

1. **Signal Detection**: Use multiple metrics (PRR, ROR, IC) for robust signal detection
2. **Clinical Context**: Always interpret findings in clinical context
3. **Time Trends**: Monitor temporal patterns in reporting
4. **Denominator Data**: Consider exposure data when available
5. **Expert Review**: Consult pharmacovigilance experts for interpretation

## Notes

- Results are limited to 100 events per request (configurable)
- Use count queries for aggregate statistics
- openFDA updates weekly
- Data from 2004 onwards
- US-centric but includes international reports

## Troubleshooting

### Rate limit errors
- Reduce request frequency
- Implement exponential backoff

### No results
- Try alternative drug names (brand vs generic)
- Check spelling
- Drug may have no reported events

### Inconsistent data
- FAERS contains real-world data with inherent variability
- Validate findings with additional sources

## Future Enhancements

- [ ] Batch processing for multiple drugs
- [ ] Temporal trend analysis
- [ ] Advanced signal detection algorithms (IC, BCPNN)
- [ ] Export to regulatory formats (E2B, MedWatch)
- [ ] Integration with drug exposure databases
