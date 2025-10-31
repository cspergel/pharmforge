# DrugCentral Adapter

Comprehensive adapter for accessing drug information, targets, indications, and pharmacology from the DrugCentral database.

## Features

- **Drug Search**: Search by drug name, structure, or InChIKey
- **Target Information**: Retrieve drug targets, genes, and mechanisms of action
- **Approved Indications**: Access FDA-approved uses and conditions
- **Pharmacology**: Get ADME (Absorption, Distribution, Metabolism, Excretion) data
- **Chemical Properties**: Molecular properties and drug-like characteristics

## API Documentation

- Base URL: https://drugcentral.org/api/v1
- Documentation: https://drugcentral.org/
- Database: Curated by NCATS/NIH

## Installation

No additional dependencies beyond the base requirements.

## Usage

### Basic Drug Search

```python
from adapters.drugcentral.adapter import DrugCentralAdapter

adapter = DrugCentralAdapter()

# Search for a drug
result = await adapter.execute("aspirin")

if result.success:
    drugs = result.data['drugs']
    for drug in drugs:
        print(f"Name: {drug['name']}")
        print(f"Generic: {drug['generic_name']}")
        print(f"Targets: {drug['num_targets']}")
        print(f"Indications: {drug['num_indications']}")
```

### Get Comprehensive Drug Information

```python
# Get all available information
result = await adapter.execute(
    "metformin",
    include_targets=True,
    include_indications=True,
    include_pharmacology=True,
    include_properties=True
)

if result.success:
    drug = result.data['drugs'][0]

    # Drug targets
    print("\nTargets:")
    for target in drug['targets']:
        print(f"  - {target['gene']} ({target['protein']})")
        print(f"    Action: {target['action']}")
        print(f"    MOA: {target['moa']}")

    # Approved indications
    print("\nIndications:")
    for indication in drug['indications']:
        print(f"  - {indication['condition']}")

    # Pharmacology
    if 'pharmacology' in drug:
        pharm = drug['pharmacology']
        print("\nPharmacokinetics:")
        print(f"  Absorption: {pharm['absorption']}")
        print(f"  Half-life: {pharm['half_life']}")
        print(f"  Metabolism: {pharm['metabolism']}")
```

### Search by InChIKey

```python
# Search using InChIKey
result = await adapter.execute({
    "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
})
```

## Input Parameters

### Simple String Input
- **drug_name** (str): Drug name to search

### Dictionary Input
- **drug** (str): Drug name
- **structure** (str): Chemical structure (SMILES or InChI)
- **inchikey** (str): InChIKey identifier

### Keyword Arguments
- **include_targets** (bool): Fetch target information (default: True)
- **include_indications** (bool): Fetch approved indications (default: True)
- **include_pharmacology** (bool): Fetch pharmacology data (default: True)
- **include_properties** (bool): Fetch chemical properties (default: True)

## Output Format

```python
{
    "drugs": [
        {
            "drug_id": "123",
            "name": "Aspirin",
            "generic_name": "acetylsalicylic acid",
            "inn": "Acetylsalicylic Acid",
            "cas_number": "50-78-2",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O",
            "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
            "num_targets": 5,
            "primary_moa": "Cyclooxygenase inhibitor",
            "targets": [
                {
                    "gene": "PTGS1",
                    "protein": "Prostaglandin-endoperoxide synthase 1",
                    "action": "inhibitor",
                    "moa": "Cyclooxygenase inhibitor",
                    "act_value": 10.0,
                    "act_type": "IC50",
                    "organism": "Homo sapiens"
                }
            ],
            "num_indications": 10,
            "indications": [
                {
                    "condition": "Pain",
                    "snomed_id": "22253000",
                    "umls_cui": "C0030193",
                    "efo_id": "EFO_0001422"
                }
            ],
            "pharmacology": {
                "absorption": "Rapid and complete from GI tract",
                "distribution": "Widely distributed; crosses placenta",
                "metabolism": "Hepatic via conjugation and hydroxylation",
                "excretion": "Renal",
                "half_life": "15-20 minutes",
                "protein_binding": "80-90%"
            },
            "properties": {
                "molecular_weight": 180.16,
                "logp": 1.19,
                "psa": 63.6,
                "hbd": 1,
                "hba": 4,
                "rotatable_bonds": 3,
                "formula": "C9H8O4"
            }
        }
    ],
    "total_found": 1
}
```

## Understanding Drug Targets

### Target Actions
- **inhibitor**: Blocks target activity
- **agonist**: Activates target
- **antagonist**: Blocks target activation
- **modulator**: Alters target activity
- **substrate**: Processed by target enzyme

### Mechanism of Action (MOA)
Primary biological mechanism:
- Enzyme inhibitor
- Receptor agonist/antagonist
- Ion channel modulator
- Transporter inhibitor
- etc.

## Example Workflows

### 1. Multi-Target Drug Analysis

```python
adapter = DrugCentralAdapter()

result = await adapter.execute("imatinib", include_targets=True)

if result.success:
    drug = result.data['drugs'][0]
    targets = drug['targets']

    print(f"{drug['name']} targets {len(targets)} proteins:")

    # Group by action
    actions = {}
    for target in targets:
        action = target['action']
        if action not in actions:
            actions[action] = []
        actions[action].append(target['gene'])

    for action, genes in actions.items():
        print(f"  {action}: {', '.join(genes)}")
```

### 2. Indication Analysis

```python
# Analyze approved indications for a drug
result = await adapter.execute("metformin", include_indications=True)

if result.success:
    drug = result.data['drugs'][0]
    indications = drug['indications']

    print(f"{drug['name']} approved for {len(indications)} indications:")

    # Group by therapeutic area (simplified)
    for indication in indications:
        condition = indication['condition']
        umls_cui = indication['umls_cui']
        print(f"  - {condition} (UMLS: {umls_cui})")
```

### 3. Drug-Like Properties Assessment

```python
# Assess Lipinski's Rule of Five
result = await adapter.execute("compound_name", include_properties=True)

if result.success:
    drug = result.data['drugs'][0]

    if 'properties' in drug:
        props = drug['properties']

        violations = 0
        print("Lipinski's Rule of Five:")

        # MW ≤ 500
        if props['molecular_weight'] > 500:
            violations += 1
            print(f"  ⚠ MW: {props['molecular_weight']} (>500)")
        else:
            print(f"  ✓ MW: {props['molecular_weight']}")

        # LogP ≤ 5
        if props['logp'] > 5:
            violations += 1
            print(f"  ⚠ LogP: {props['logp']} (>5)")
        else:
            print(f"  ✓ LogP: {props['logp']}")

        # HBD ≤ 5
        if props['hbd'] > 5:
            violations += 1
            print(f"  ⚠ H-bond donors: {props['hbd']} (>5)")
        else:
            print(f"  ✓ H-bond donors: {props['hbd']}")

        # HBA ≤ 10
        if props['hba'] > 10:
            violations += 1
            print(f"  ⚠ H-bond acceptors: {props['hba']} (>10)")
        else:
            print(f"  ✓ H-bond acceptors: {props['hba']}")

        print(f"\nLipinski violations: {violations}/4")
        if violations <= 1:
            print("✓ Drug-like properties")
```

### 4. Pharmacokinetic Profile

```python
result = await adapter.execute("warfarin", include_pharmacology=True)

if result.success:
    drug = result.data['drugs'][0]

    if 'pharmacology' in drug:
        pk = drug['pharmacology']

        print("Pharmacokinetic Profile:")
        print(f"  Absorption: {pk['absorption']}")
        print(f"  Distribution: {pk['distribution']}")
        print(f"  Metabolism: {pk['metabolism']}")
        print(f"  Excretion: {pk['excretion']}")
        print(f"  Half-life: {pk['half_life']}")
        print(f"  Protein binding: {pk['protein_binding']}")
```

### 5. Target-Based Drug Discovery

```python
# Find all drugs targeting a specific gene
target_gene = "EGFR"

# Search for the gene (simplified - would need specific endpoint)
result = await adapter.execute(target_gene)

if result.success:
    print(f"Drugs targeting {target_gene}:")

    for drug in result.data['drugs']:
        if 'targets' in drug:
            for target in drug['targets']:
                if target['gene'] == target_gene:
                    print(f"  - {drug['name']}")
                    print(f"    Action: {target['action']}")
                    if target['act_value']:
                        print(f"    {target['act_type']}: {target['act_value']}")
```

## Data Quality

DrugCentral is curated by NCATS and includes:
- FDA-approved drugs
- Experimental drugs in clinical trials
- Manually curated data from literature and databases
- Regular updates from FDA Orange Book and other sources

### Data Sources
- FDA Orange Book
- ChEMBL
- PubChem
- UniProt
- OMIM
- SNOMED CT
- UMLS

## Chemical Identifiers

### Supported Formats
- **Drug names**: Brand and generic names
- **InChIKey**: Standard chemical identifier
- **CAS Number**: Chemical Abstracts Service registry number
- **SMILES**: Simplified molecular-input line-entry system

## Therapeutic Applications

### Drug Repurposing
Use target and indication data to identify repurposing opportunities:

```python
# Find drugs with similar targets
result1 = await adapter.execute("drug1", include_targets=True)
result2 = await adapter.execute("drug2", include_targets=True)

targets1 = set(t['gene'] for t in result1.data['drugs'][0]['targets'])
targets2 = set(t['gene'] for t in result2.data['drugs'][0]['targets'])

common_targets = targets1 & targets2
print(f"Common targets: {common_targets}")
```

## Error Handling

```python
result = await adapter.execute("unknown_drug")

if not result.success:
    print(f"Error: {result.error}")
elif result.data['total_found'] == 0:
    print("Drug not found in DrugCentral")
else:
    # Process data
    pass
```

## Rate Limiting

```python
adapter = DrugCentralAdapter()
adapter.config['rate_limit_delay'] = 0.5  # 500ms between requests
```

## Notes

- Returns top 5 results for ambiguous searches
- Some fields may be None if data not available
- Pharmacology data may be limited for older drugs
- Target information includes experimental data

## Limitations

- API coverage may be incomplete
- Some endpoints may require authentication
- Pharmacokinetic data quality varies
- Not all drugs in FDA database are included

## Integration with Other Adapters

### With Reactome
```python
# Get drug targets, then analyze pathways
drug_result = await drugcentral_adapter.execute("drug", include_targets=True)
target_genes = [t['gene'] for t in drug_result.data['drugs'][0]['targets']]

pathway_result = await reactome_adapter.execute(target_genes)
```

### With FAERS
```python
# Compare approved indications with adverse events
drug_result = await drugcentral_adapter.execute("drug", include_indications=True)
safety_result = await faers_adapter.execute("drug")
```

## Future Enhancements

- [ ] Structure-based search
- [ ] Drug-drug interaction data
- [ ] Dosing information
- [ ] Patent and regulatory status
- [ ] Bioavailability predictions
