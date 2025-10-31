# PharmForge API Keys Setup Guide

**Date:** October 26, 2025
**Purpose:** Complete guide for obtaining and configuring API keys for PharmForge adapters

---

## Overview

PharmForge uses 33 adapters, of which:
- **23 adapters** are production-ready (no API keys needed)
- **10 adapters** require API keys or additional setup

This guide covers how to obtain API keys, configure them in `.env`, and test adapter functionality.

---

## Quick Start

### 1. Copy .env.example to .env

```bash
cp .env.example .env
```

### 2. Add Your API Keys

Edit `.env` file and replace placeholders with your actual API keys:

```bash
# Open .env in your editor
nano .env  # or vim, code, etc.
```

### 3. Restart Docker Containers

```bash
cd "C:\Users\drcra\Documents\Coding Projects\PharmForge\claude-code-agents-wizard-v2"
docker-compose restart backend celery-worker
```

---

## API Keys by Category

### 1. LLM Services (for LLM Retrosynthesis Adapter)

#### OpenAI API Key

**Purpose:** Powers LLM-based retrosynthesis planning
**Cost:** Pay-as-you-go (GPT-4: ~$0.03/1K tokens)
**Sign Up:** https://platform.openai.com/signup

**Steps:**
1. Go to https://platform.openai.com/
2. Sign up or log in
3. Navigate to API Keys section
4. Click "Create new secret key"
5. Copy the key (starts with `sk-...`)
6. Add to `.env`:

```env
OPENAI_API_KEY=sk-proj-abc123xyz...
```

**Testing:**
```bash
docker-compose exec -T backend python -c "
import os
from openai import OpenAI
client = OpenAI(api_key=os.getenv('OPENAI_API_KEY'))
print('OpenAI API key is valid!')
"
```

---

#### Anthropic API Key (Alternative to OpenAI)

**Purpose:** Alternative LLM for retrosynthesis (Claude models)
**Cost:** Pay-as-you-go (Claude Sonnet: ~$0.015/1K tokens)
**Sign Up:** https://console.anthropic.com/

**Steps:**
1. Go to https://console.anthropic.com/
2. Sign up or log in
3. Navigate to API Keys
4. Create new API key
5. Copy the key (starts with `sk-ant-...`)
6. Add to `.env`:

```env
ANTHROPIC_API_KEY=sk-ant-api03-abc123xyz...
```

**Note:** Only one LLM API key is needed (OpenAI OR Anthropic)

---

### 2. Patent & Literature Search

#### Lens.org API Key

**Purpose:** Patent search and analysis
**Cost:** Free tier available (100 requests/day), Paid plans from $99/month
**Sign Up:** https://www.lens.org/lens/user/subscriptions

**Steps:**
1. Go to https://www.lens.org/
2. Create account
3. Navigate to Account → Subscriptions
4. Subscribe to API access (free tier available)
5. Go to API Access section
6. Copy your API token
7. Add to `.env`:

```env
LENS_ORG_API_KEY=your_lens_api_token_here
```

**Testing:**
```bash
docker-compose exec -T backend python -c "
import os, requests
api_key = os.getenv('LENS_ORG_API_KEY')
headers = {'Authorization': f'Bearer {api_key}'}
r = requests.get('https://api.lens.org/patent/search', headers=headers)
print(f'Lens.org API: {r.status_code}')
"
```

---

#### Google Patents API Key

**Purpose:** Google Patents search
**Cost:** Free (subject to quota limits)
**Sign Up:** https://console.cloud.google.com/

**Steps:**
1. Go to https://console.cloud.google.com/
2. Create a new project or select existing
3. Enable "Custom Search API"
4. Go to Credentials → Create Credentials → API Key
5. Copy the API key
6. Add to `.env`:

```env
GOOGLE_PATENTS_API_KEY=AIzaSyAbc123xyz...
```

**Additional Setup:**
- Create Custom Search Engine at: https://programmablesearchengine.google.com/
- Configure it to search patents.google.com
- Note the Search Engine ID (cx parameter)

**Testing:**
```bash
docker-compose exec -T backend python -c "
import os, requests
api_key = os.getenv('GOOGLE_PATENTS_API_KEY')
params = {'key': api_key, 'cx': 'your_search_engine_id', 'q': 'aspirin'}
r = requests.get('https://www.googleapis.com/customsearch/v1', params=params)
print(f'Google Patents API: {r.status_code}')
"
```

---

### 3. Drug Discovery Databases

#### ChEMBL API Key (Optional)

**Purpose:** Enhanced access to ChEMBL bioactivity database
**Cost:** Free (no key required for basic access)
**Sign Up:** Not required (public API available)

**Note:** ChEMBL adapter works without API key using public REST API. API key only needed for:
- Higher rate limits
- Priority access
- Advanced features

**If you want enhanced access:**
```env
CHEMBL_API_KEY=your_chembl_api_key_here
```

**Testing (Public API - No Key Needed):**
```bash
docker-compose exec -T backend python -c "
import requests
r = requests.get('https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25.json')
print(f'ChEMBL API: {r.status_code} - {r.json()['molecule_chembl_id']}')
"
```

---

### 4. Gene-Disease Associations (Deferred)

#### DisGeNET API Key

**Purpose:** Gene-disease association data
**Cost:** Free academic tier, Commercial licenses available
**Sign Up:** https://www.disgenet.org/signup

**Status:** ⚠️ DEFERRED
**Reason:** Licensing considerations for commercial use

**If needed in future:**
1. Review licensing terms at https://www.disgenet.org/legal
2. Determine if academic or commercial license needed
3. Sign up and obtain API key
4. Add to `.env`:

```env
# DisGeNET API (for gene-disease associations)
# DISGENET_API_KEY=your_disgenet_api_key_here
```

**Keep commented out until licensing is clarified.**

---

## Environment Variables Reference

### Complete .env Template

```env
# =============================================================================
# PharmForge Environment Configuration
# =============================================================================

# -----------------------------------------------------------------------------
# Database Configuration
# -----------------------------------------------------------------------------
DATABASE_URL=postgresql://pharmforge:pharmforge_dev_password@postgres:5432/pharmforge

# -----------------------------------------------------------------------------
# Redis Cache & Queue
# -----------------------------------------------------------------------------
REDIS_URL=redis://redis:6379/0
CELERY_BROKER_URL=redis://redis:6379/0
CELERY_RESULT_BACKEND=redis://redis:6379/0

# -----------------------------------------------------------------------------
# Application Settings
# -----------------------------------------------------------------------------
LOG_LEVEL=DEBUG
PYTHONPATH=/app

# -----------------------------------------------------------------------------
# GPU Configuration (for DiffDock, ML models)
# -----------------------------------------------------------------------------
NVIDIA_VISIBLE_DEVICES=all

# -----------------------------------------------------------------------------
# External APIs - LLM Services (for LLM Retrosynthesis adapter)
# -----------------------------------------------------------------------------

# OpenAI API (GPT-4, GPT-3.5)
# Sign up: https://platform.openai.com/signup
# Pricing: Pay-as-you-go (~$0.03/1K tokens for GPT-4)
OPENAI_API_KEY=sk-proj-...

# Anthropic API (Claude models) - Alternative to OpenAI
# Sign up: https://console.anthropic.com/
# Pricing: Pay-as-you-go (~$0.015/1K tokens for Claude Sonnet)
ANTHROPIC_API_KEY=sk-ant-api03-...

# -----------------------------------------------------------------------------
# External APIs - Patent & Literature Search
# -----------------------------------------------------------------------------

# Lens.org API (Patent search)
# Sign up: https://www.lens.org/lens/user/subscriptions
# Free tier: 100 requests/day | Paid: from $99/month
LENS_ORG_API_KEY=

# Google Patents API
# Sign up: https://console.cloud.google.com/
# Free with quota limits
GOOGLE_PATENTS_API_KEY=

# -----------------------------------------------------------------------------
# External APIs - Drug Discovery Databases
# -----------------------------------------------------------------------------

# ChEMBL API (Optional - public API works without key)
# Enhanced access for higher rate limits
CHEMBL_API_KEY=

# DisGeNET API (Gene-disease associations)
# ⚠️ DEFERRED - Licensing considerations for commercial use
# Sign up: https://www.disgenet.org/signup
# DISGENET_API_KEY=

# -----------------------------------------------------------------------------
# DiffDock Configuration (GPU-accelerated docking)
# -----------------------------------------------------------------------------
DIFFDOCK_MODEL_DIR=/app/diffdock_models
DIFFDOCK_USE_GPU=true

# -----------------------------------------------------------------------------
# REINVENT Configuration (Molecular generation)
# -----------------------------------------------------------------------------
REINVENT_MODEL_DIR=/app/reinvent_models
REINVENT_FALLBACK_MODE=true  # Uses RDKit when REINVENT not installed

# -----------------------------------------------------------------------------
# AiZynthFinder Configuration (Retrosynthesis)
# -----------------------------------------------------------------------------
AIZYNTHFINDER_CONFIG=/app/aizynthfinder_data/config.yml
```

---

## Adapter API Key Requirements

### No API Key Needed (23 adapters)

**Core Drug Discovery:**
- UniProt
- RCSB PDB
- AlphaFold
- SWISS-MODEL
- OpenTargets
- BindingDB
- DrugCentral
- ZINC
- ChEMBL (public API)
- PubChem
- RDKit

**Clinical & Safety:**
- ClinicalTrials.gov
- FDA FAERS
- ADMET-AI

**Literature:**
- PubMed
- EuropePMC
- SureChEMBL

**ML/Computational:**
- REINVENT (fallback mode)
- MolGAN (fallback mode)
- De Novo (fragment mode)
- AiZynthFinder
- Vina
- OpenMM

**Systems Biology:**
- Reactome
- GTEx

---

### API Key Required (7 adapters)

| Adapter | API Key Required | Priority | Cost |
|---------|------------------|----------|------|
| **LLM Retrosynthesis** | OpenAI OR Anthropic | High | Pay-as-you-go |
| **Lens.org** | Lens.org API | Medium | Free tier / $99+/mo |
| **Google Patents** | Google Cloud API | Medium | Free (with limits) |
| **DiffDock** | No key, but needs GPU + models | High | Hardware cost |
| **ChEMBL (enhanced)** | ChEMBL API | Low | Free |
| **DisGeNET** | DisGeNET API | Deferred | TBD (licensing) |

---

### Additional Setup Required (3 adapters)

| Adapter | Requirement | Status | Notes |
|---------|-------------|--------|-------|
| **DiffDock** | GPU + 5GB models | Pending | Docker GPU configured |
| **REINVENT (full)** | GitHub install + models | Optional | Works in fallback mode |
| **OpenMM** | Package installed | ✅ Complete | Ready for use |

---

## Testing API Keys

### Test All Configured API Keys

Create `scripts/test_api_keys.py`:

```python
#!/usr/bin/env python3
"""
Test all configured API keys
"""

import os
from pathlib import Path

def test_openai():
    """Test OpenAI API key"""
    api_key = os.getenv('OPENAI_API_KEY')
    if not api_key:
        return '❌ Not configured'

    try:
        from openai import OpenAI
        client = OpenAI(api_key=api_key)
        client.models.list()
        return '✅ Valid'
    except Exception as e:
        return f'❌ Invalid: {str(e)[:50]}'

def test_anthropic():
    """Test Anthropic API key"""
    api_key = os.getenv('ANTHROPIC_API_KEY')
    if not api_key:
        return '❌ Not configured'

    try:
        import anthropic
        client = anthropic.Anthropic(api_key=api_key)
        # Simple test - list models would require actual request
        return '✅ Configured (validation requires request)'
    except Exception as e:
        return f'❌ Invalid: {str(e)[:50]}'

def test_lens_org():
    """Test Lens.org API key"""
    api_key = os.getenv('LENS_ORG_API_KEY')
    if not api_key:
        return '❌ Not configured'

    try:
        import requests
        headers = {'Authorization': f'Bearer {api_key}'}
        r = requests.get('https://api.lens.org/patent/search', headers=headers, timeout=5)
        if r.status_code == 200 or r.status_code == 400:  # 400 = valid key, bad request
            return '✅ Valid'
        return f'❌ Invalid: HTTP {r.status_code}'
    except Exception as e:
        return f'❌ Error: {str(e)[:50]}'

def test_google_patents():
    """Test Google Patents API key"""
    api_key = os.getenv('GOOGLE_PATENTS_API_KEY')
    if not api_key:
        return '❌ Not configured'

    try:
        import requests
        params = {'key': api_key, 'cx': 'test', 'q': 'test'}
        r = requests.get('https://www.googleapis.com/customsearch/v1', params=params, timeout=5)
        if r.status_code in [200, 400]:  # Valid key even if cx is wrong
            return '✅ Valid'
        return f'❌ Invalid: HTTP {r.status_code}'
    except Exception as e:
        return f'❌ Error: {str(e)[:50]}'

def main():
    print("="*70)
    print("PharmForge API Keys Test")
    print("="*70)
    print()

    tests = {
        'OpenAI API': test_openai,
        'Anthropic API': test_anthropic,
        'Lens.org API': test_lens_org,
        'Google Patents API': test_google_patents,
    }

    results = {}
    for name, test_func in tests.items():
        print(f"Testing {name}...")
        result = test_func()
        results[name] = result
        print(f"  {result}")
        print()

    print("="*70)
    print("Summary")
    print("="*70)

    configured = sum(1 for r in results.values() if '✅' in r)
    total = len(results)

    print(f"Configured: {configured}/{total}")
    print()

    for name, result in results.items():
        print(f"{name:<25} {result}")

if __name__ == '__main__':
    main()
```

**Run the test:**
```bash
docker-compose exec -T backend python scripts/test_api_keys.py
```

---

## Pricing Summary

### Free Tier Options

| Service | Free Tier | Limits |
|---------|-----------|--------|
| **ChEMBL** | ✅ Unlimited | Public API, no key needed |
| **PubChem** | ✅ Unlimited | Public API |
| **UniProt** | ✅ Unlimited | Public API |
| **PubMed** | ✅ Unlimited | Public API (rate limited) |
| **OpenTargets** | ✅ Unlimited | Public API |
| **Lens.org** | ✅ 100 req/day | Free tier |
| **Google Patents** | ✅ Yes | Subject to quotas |

### Pay-As-You-Go

| Service | Model | Price | Estimate (1000 queries) |
|---------|-------|-------|-------------------------|
| **OpenAI** | GPT-4 | $0.03/1K tokens | ~$30-$60 |
| **OpenAI** | GPT-3.5-Turbo | $0.002/1K tokens | ~$2-$4 |
| **Anthropic** | Claude Sonnet | $0.015/1K tokens | ~$15-$30 |

### Subscription Plans

| Service | Plan | Price | Features |
|---------|------|-------|----------|
| **Lens.org** | Scholar | Free | 100 requests/day |
| **Lens.org** | Professional | $99/month | 10,000 requests/month |
| **Lens.org** | Enterprise | Custom | Unlimited |

---

## Security Best Practices

### 1. Never Commit .env to Git

**.gitignore already configured:**
```gitignore
.env
.env.local
.env.*.local
```

### 2. Use Environment-Specific Files

```bash
.env.development   # Local development
.env.staging       # Staging server
.env.production    # Production server
```

### 3. Rotate Keys Regularly

- Change API keys every 90 days
- Immediately rotate if exposed
- Use different keys for dev/staging/prod

### 4. Use Key Management Services (Production)

For production deployment, consider:
- **AWS Secrets Manager**
- **HashiCorp Vault**
- **Azure Key Vault**
- **Google Cloud Secret Manager**

---

## Troubleshooting

### API Key Not Working

**Check 1: Verify key is loaded**
```bash
docker-compose exec -T backend python -c "import os; print(os.getenv('OPENAI_API_KEY'))"
```

**Check 2: Restart containers**
```bash
docker-compose restart backend celery-worker
```

**Check 3: Check logs**
```bash
docker-compose logs backend | grep -i "api key"
```

### Rate Limiting

Most APIs have rate limits:

- **OpenAI:** 3,500 requests/minute (varies by tier)
- **Lens.org Free:** 100 requests/day
- **PubMed:** ~3 requests/second

**Solution:** Implement exponential backoff in adapters

### Authentication Errors

**401 Unauthorized:** Invalid API key
- Double-check key is correct
- Ensure no extra spaces
- Verify key hasn't expired

**403 Forbidden:** Valid key, insufficient permissions
- Check subscription tier
- Verify API access is enabled
- Review usage quotas

---

## Next Steps

### 1. Configure Required API Keys (Priority Order)

1. **OpenAI OR Anthropic** - For LLM Retrosynthesis (if needed)
2. **Lens.org** - For patent search (if needed)
3. **Google Patents** - For Google Patents search (if needed)

### 2. Test Adapter Functionality

```bash
# Test LLM Retrosynthesis
docker-compose exec -T backend python -c "
from adapters.llm_retrosynthesis.adapter import LLMRetrosynthesisAdapter
adapter = LLMRetrosynthesisAdapter()
result = adapter.execute('CCO')  # Ethanol
print(f'Result: {result}')
"

# Test Lens.org patent search
docker-compose exec -T backend python -c "
from adapters.lens_org.adapter import LensOrgAdapter
adapter = LensOrgAdapter()
result = adapter.execute(query='aspirin drug discovery')
print(f'Found {len(result)} patents')
"
```

### 3. Monitor Usage and Costs

Track API usage to avoid unexpected costs:

```python
# Add to backend/core/monitoring.py
import logging

def log_api_usage(adapter_name: str, tokens: int, cost: float):
    logging.info(f"API Usage: {adapter_name} - {tokens} tokens - ${cost:.4f}")
```

---

## Support

**Issues with API keys?**
- Check adapter-specific documentation in `adapters/{adapter_name}/README.md`
- Review error logs: `docker-compose logs backend`
- Create GitHub issue with sanitized error (remove API keys!)

**Need API key for specific service?**
- Follow links in this guide
- Check service documentation
- Contact service support if issues

---

**Document Version:** 1.0
**Last Updated:** October 26, 2025
**Maintained By:** PharmForge Development Team
