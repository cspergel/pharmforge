"""
Example usage of ClinicalTrials.gov Adapter
Demonstrates how to query clinical trial data
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

from adapters.clinicaltrials.adapter import ClinicalTrialsAdapter


async def example_basic_search():
    """Example: Basic drug name search"""
    print("\n" + "="*80)
    print("Example 1: Basic Drug Search - Aspirin")
    print("="*80)

    adapter = ClinicalTrialsAdapter()

    # Search for aspirin trials
    result = await adapter.execute("aspirin")

    if result.success:
        data = result.data
        print(f"\nTotal trials found: {data['total_count']}")
        print(f"Studies returned: {len(data['studies'])}")

        # Print summary
        summary = data['summary']
        print(f"\nSummary:")
        print(f"  Total enrollment: {summary['total_enrollment']}")
        print(f"  Studies with results: {summary['studies_with_results']}")

        print(f"\nBy Phase:")
        for phase, count in summary['by_phase'].items():
            print(f"  {phase}: {count}")

        print(f"\nBy Status:")
        for status, count in summary['by_status'].items():
            print(f"  {status}: {count}")

        # Show first few trials
        print(f"\nFirst 3 trials:")
        for i, study in enumerate(data['studies'][:3], 1):
            print(f"\n  {i}. {study['title']}")
            print(f"     NCT ID: {study['nct_id']}")
            print(f"     Status: {study['status']}")
            print(f"     Phase: {', '.join(study['phase']) if study['phase'] else 'N/A'}")
            print(f"     Enrollment: {study['enrollment']}")
            if study['interventions']:
                print(f"     Interventions: {', '.join([i['name'] for i in study['interventions'][:3]])}")
    else:
        print(f"Error: {result.error}")


async def example_advanced_search():
    """Example: Advanced search with multiple parameters"""
    print("\n" + "="*80)
    print("Example 2: Advanced Search - Cancer + Phase 3")
    print("="*80)

    adapter = ClinicalTrialsAdapter()

    # Search for cancer trials in phase 3
    search_params = {
        "condition": "cancer",
        "phase": "Phase 3",
        "status": "Completed"
    }

    result = await adapter.execute(search_params, max_studies=10)

    if result.success:
        data = result.data
        print(f"\nTotal trials found: {data['total_count']}")
        print(f"\nShowing {len(data['studies'])} completed Phase 3 cancer trials:")

        for i, study in enumerate(data['studies'], 1):
            print(f"\n{i}. {study['title']}")
            print(f"   NCT ID: {study['nct_id']}")
            print(f"   Enrollment: {study['enrollment']}")
            print(f"   Completion Date: {study.get('completion_date', 'N/A')}")

            if study['interventions']:
                drugs = [i['name'] for i in study['interventions'] if i['type'] == 'DRUG']
                if drugs:
                    print(f"   Drug Interventions: {', '.join(drugs[:3])}")

            if study['primary_outcomes']:
                print(f"   Primary Outcome: {study['primary_outcomes'][0]['measure']}")
    else:
        print(f"Error: {result.error}")


async def example_specific_trial():
    """Example: Search by NCT ID"""
    print("\n" + "="*80)
    print("Example 3: Specific Trial by NCT ID")
    print("="*80)

    adapter = ClinicalTrialsAdapter()

    # Search for a specific trial (use a real NCT ID)
    search_params = {
        "nct_id": "NCT04280705"  # Example COVID-19 trial
    }

    result = await adapter.execute(search_params)

    if result.success:
        data = result.data
        if data['studies']:
            study = data['studies'][0]
            print(f"\nTrial Details:")
            print(f"  Title: {study['title']}")
            print(f"  NCT ID: {study['nct_id']}")
            print(f"  Status: {study['status']}")
            print(f"  Phase: {', '.join(study['phase']) if study['phase'] else 'N/A'}")
            print(f"  Enrollment: {study['enrollment']}")
            print(f"  Start Date: {study.get('start_date', 'N/A')}")
            print(f"  Completion Date: {study.get('completion_date', 'N/A')}")

            print(f"\n  Conditions:")
            for condition in study['conditions']:
                print(f"    - {condition}")

            print(f"\n  Interventions:")
            for intervention in study['interventions']:
                print(f"    - {intervention['name']} ({intervention['type']})")

            print(f"\n  Primary Outcomes:")
            for outcome in study['primary_outcomes']:
                print(f"    - {outcome['measure']}")
                print(f"      Timeframe: {outcome.get('timeFrame', 'N/A')}")
        else:
            print("No trial found with that NCT ID")
    else:
        print(f"Error: {result.error}")


async def example_drug_condition_combination():
    """Example: Search for specific drug in specific condition"""
    print("\n" + "="*80)
    print("Example 4: Drug + Condition Search - Pembrolizumab + Melanoma")
    print("="*80)

    adapter = ClinicalTrialsAdapter()

    search_params = {
        "drug": "pembrolizumab",
        "condition": "melanoma"
    }

    result = await adapter.execute(search_params, max_studies=5)

    if result.success:
        data = result.data
        print(f"\nTotal trials found: {data['total_count']}")

        summary = data['summary']
        print(f"\nPhase Distribution:")
        for phase, count in summary['by_phase'].items():
            print(f"  {phase}: {count}")

        print(f"\nRecent Trials:")
        for i, study in enumerate(data['studies'][:5], 1):
            print(f"\n{i}. {study['title'][:80]}...")
            print(f"   Status: {study['status']}")
            print(f"   Phase: {', '.join(study['phase']) if study['phase'] else 'N/A'}")
            print(f"   Enrollment: {study['enrollment']}")
            print(f"   Has Results: {'Yes' if study.get('has_results') else 'No'}")
    else:
        print(f"Error: {result.error}")


async def example_with_caching():
    """Example: Demonstrate caching"""
    print("\n" + "="*80)
    print("Example 5: Caching Demonstration")
    print("="*80)

    adapter = ClinicalTrialsAdapter()

    print("\nFirst call (will fetch from API)...")
    import time
    start = time.time()
    result1 = await adapter("ibuprofen", use_cache=True)
    time1 = time.time() - start

    print(f"Time: {time1:.2f}s")
    print(f"Cache hit: {result1.cache_hit}")
    print(f"Studies found: {len(result1.data['studies']) if result1.success else 0}")

    print("\nSecond call (should use cache)...")
    start = time.time()
    result2 = await adapter("ibuprofen", use_cache=True)
    time2 = time.time() - start

    print(f"Time: {time2:.2f}s")
    print(f"Cache hit: {result2.cache_hit}")
    print(f"Studies found: {len(result2.data['studies']) if result2.success else 0}")

    print(f"\nSpeedup: {time1/time2:.1f}x faster with cache")


async def main():
    """Run all examples"""
    print("\n" + "="*80)
    print("ClinicalTrials.gov Adapter - Example Usage")
    print("="*80)

    try:
        await example_basic_search()
        await asyncio.sleep(1)  # Be nice to the API

        await example_advanced_search()
        await asyncio.sleep(1)

        await example_specific_trial()
        await asyncio.sleep(1)

        await example_drug_condition_combination()
        await asyncio.sleep(1)

        await example_with_caching()

    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()

    print("\n" + "="*80)
    print("Examples completed!")
    print("="*80 + "\n")


if __name__ == "__main__":
    asyncio.run(main())
