"""
Example usage of FDA FAERS Adapter
Demonstrates how to query adverse event data and safety signals
"""
import asyncio
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

from adapters.fda_faers.adapter import FDAFAERSAdapter


async def example_basic_search():
    """Example: Basic drug name search"""
    print("\n" + "="*80)
    print("Example 1: Basic Adverse Event Search - Aspirin")
    print("="*80)

    adapter = FDAFAERSAdapter()

    # Search for aspirin adverse events
    result = await adapter.execute("aspirin")

    if result.success:
        data = result.data
        print(f"\nTotal adverse events found: {data['total_count']:,}")
        print(f"Events returned: {len(data['events'])}")

        # Print safety metrics
        metrics = data['safety_metrics']
        print(f"\nSafety Metrics:")
        print(f"  Unique reaction types: {metrics['unique_reactions']}")

        print(f"\n  Top 10 Adverse Reactions:")
        for i, reaction in enumerate(metrics['top_reactions'], 1):
            print(f"    {i}. {reaction['term']}")
            print(f"       Count: {reaction['count']:,} ({reaction['percentage']:.1f}%)")

        # Print demographics
        if 'demographics' in data and data['demographics']:
            demo = data['demographics']
            print(f"\nDemographics:")
            print(f"  Sex Distribution: {demo['sex_distribution']}")
            print(f"  Age Distribution: {demo['age_distribution']}")
            print(f"  Serious Events: {demo['serious_events']:,} ({demo['serious_percentage']:.1f}%)")
    else:
        print(f"Error: {result.error}")


async def example_specific_drug():
    """Example: Specific drug with detailed analysis"""
    print("\n" + "="*80)
    print("Example 2: Detailed Analysis - Metformin")
    print("="*80)

    adapter = FDAFAERSAdapter()

    result = await adapter.execute("metformin", include_demographics=True)

    if result.success:
        data = result.data
        print(f"\nTotal events: {data['total_count']:,}")

        # Safety metrics
        metrics = data['safety_metrics']
        print(f"\nTop 5 Most Common Reactions:")
        for i, reaction in enumerate(metrics['top_reactions'][:5], 1):
            print(f"\n  {i}. {reaction['term']}")
            print(f"     Occurrences: {reaction['count']:,}")
            print(f"     Percentage: {reaction['percentage']:.2f}%")

        # PRR signals if available
        if 'prr_signals' in metrics and metrics['prr_signals']:
            print(f"\n  Proportional Reporting Ratio (PRR) Signals:")
            print(f"  (PRR > 2.0 may indicate safety signal)")
            for reaction, prr in metrics['prr_signals'].items():
                signal = "⚠️ SIGNAL" if prr > 2.0 else "Normal"
                print(f"    {reaction}: {prr} ({signal})")

        # Show sample events
        print(f"\n  Sample Adverse Event Reports:")
        for i, event in enumerate(data['events'][:3], 1):
            print(f"\n  Event {i}:")
            print(f"    Report ID: {event['report_id']}")
            print(f"    Received: {event['receive_date']}")
            print(f"    Serious: {'Yes' if event['serious'] == 1 else 'No'}")
            print(f"    Patient: Age {event.get('age', 'Unknown')}, Sex {event.get('sex', 'Unknown')}")

            if event['reactions']:
                print(f"    Reactions:")
                for reaction in event['reactions'][:3]:
                    print(f"      - {reaction['term']}")

            if event['drugs']:
                print(f"    Drugs involved:")
                for drug in event['drugs'][:3]:
                    role_map = {1: "Suspect", 2: "Concomitant", 3: "Interacting"}
                    role = role_map.get(drug.get('role'), 'Unknown')
                    print(f"      - {drug.get('name')} ({role})")
    else:
        print(f"Error: {result.error}")


async def example_safety_comparison():
    """Example: Compare safety profiles of two drugs"""
    print("\n" + "="*80)
    print("Example 3: Safety Profile Comparison - Ibuprofen vs Naproxen")
    print("="*80)

    adapter = FDAFAERSAdapter()

    drugs = ["ibuprofen", "naproxen"]
    results = {}

    for drug in drugs:
        result = await adapter.execute(drug)
        if result.success:
            results[drug] = result.data
        await asyncio.sleep(0.5)  # Rate limiting

    # Compare the results
    print("\nComparison:")
    print(f"\n{'Metric':<30} {'Ibuprofen':<15} {'Naproxen':<15}")
    print("-" * 60)

    for drug in drugs:
        if drug in results:
            data = results[drug]
            if drug == "ibuprofen":
                total_events = f"{data['total_count']:,}"
                unique_reactions = data['safety_metrics']['unique_reactions']
                serious_pct = data['demographics'].get('serious_percentage', 0)
            else:
                print(f"{'Total Events':<30} {results['ibuprofen']['total_count']:,<15} {data['total_count']:,<15}")
                print(f"{'Unique Reactions':<30} {results['ibuprofen']['safety_metrics']['unique_reactions']:<15} {data['safety_metrics']['unique_reactions']:<15}")
                print(f"{'Serious Events %':<30} {results['ibuprofen']['demographics'].get('serious_percentage', 0):<15.1f} {data['demographics'].get('serious_percentage', 0):<15.1f}")

    # Compare top reactions
    print("\nTop 5 Reactions Comparison:")
    for i in range(5):
        ibu_reaction = results['ibuprofen']['safety_metrics']['top_reactions'][i]['term'] if len(results['ibuprofen']['safety_metrics']['top_reactions']) > i else 'N/A'
        nap_reaction = results['naproxen']['safety_metrics']['top_reactions'][i]['term'] if len(results['naproxen']['safety_metrics']['top_reactions']) > i else 'N/A'
        print(f"  {i+1}. {ibu_reaction:<30} | {nap_reaction}")


async def example_search_by_ingredient():
    """Example: Search by active ingredient"""
    print("\n" + "="*80)
    print("Example 4: Search by Active Ingredient - Acetaminophen")
    print("="*80)

    adapter = FDAFAERSAdapter()

    search_params = {
        "ingredient": "acetaminophen"
    }

    result = await adapter.execute(search_params)

    if result.success:
        data = result.data
        print(f"\nTotal events for acetaminophen: {data['total_count']:,}")

        metrics = data['safety_metrics']
        print(f"\nTop Adverse Reactions:")
        for i, reaction in enumerate(metrics['top_reactions'][:10], 1):
            print(f"  {i:2d}. {reaction['term']:<40} {reaction['count']:>6,} events ({reaction['percentage']:>5.1f}%)")

        # Demographics
        demo = data['demographics']
        print(f"\nDemographic Analysis:")
        print(f"  Sex Distribution:")
        for sex, count in demo['sex_distribution'].items():
            print(f"    {sex}: {count:,}")

        print(f"\n  Age Distribution:")
        for age_group, count in sorted(demo['age_distribution'].items()):
            print(f"    {age_group}: {count:,}")
    else:
        print(f"Error: {result.error}")


async def example_pharmacovigilance_report():
    """Example: Generate a pharmacovigilance-style report"""
    print("\n" + "="*80)
    print("Example 5: Pharmacovigilance Report - Warfarin")
    print("="*80)

    adapter = FDAFAERSAdapter()

    result = await adapter.execute("warfarin", include_demographics=True)

    if result.success:
        data = result.data
        metrics = data['safety_metrics']
        demo = data['demographics']

        print("\n" + "="*60)
        print("PHARMACOVIGILANCE SAFETY REPORT")
        print("Drug: Warfarin")
        print("Data Source: FDA FAERS (openFDA)")
        print("="*60)

        print(f"\n1. OVERALL ADVERSE EVENT SUMMARY")
        print(f"   Total Reports: {data['total_count']:,}")
        print(f"   Unique Reaction Types: {metrics['unique_reactions']}")
        print(f"   Serious Events: {demo['serious_events']:,} ({demo['serious_percentage']:.1f}%)")

        print(f"\n2. DEMOGRAPHIC PROFILE")
        print(f"   Sex Distribution:")
        for sex, count in demo['sex_distribution'].items():
            pct = (count / sum(demo['sex_distribution'].values())) * 100
            print(f"     {sex}: {count:,} ({pct:.1f}%)")

        print(f"\n   Age Distribution:")
        for age_group, count in sorted(demo['age_distribution'].items()):
            total_age = sum(demo['age_distribution'].values())
            pct = (count / total_age) * 100 if total_age > 0 else 0
            print(f"     {age_group}: {count:,} ({pct:.1f}%)")

        print(f"\n3. TOP ADVERSE REACTIONS (MedDRA Terms)")
        print(f"   {'Rank':<5} {'Reaction':<40} {'Count':<10} {'%':<8}")
        print(f"   " + "-"*63)
        for i, reaction in enumerate(metrics['top_reactions'][:15], 1):
            print(f"   {i:<5} {reaction['term']:<40} {reaction['count']:<10,} {reaction['percentage']:<8.2f}")

        print(f"\n4. SAMPLE CASE REPORTS")
        for i, event in enumerate(data['events'][:5], 1):
            print(f"\n   Case {i}: {event['report_id']}")
            print(f"     Received: {event['receive_date']}")
            print(f"     Serious: {'Yes' if event['serious'] == 1 else 'No'}")

            if event['reactions']:
                reactions_str = ", ".join([r['term'] for r in event['reactions'][:3]])
                print(f"     Reactions: {reactions_str}")

            if event['drugs']:
                suspect_drugs = [d['name'] for d in event['drugs'] if d.get('role') == 1]
                if suspect_drugs:
                    print(f"     Suspect Drugs: {', '.join(suspect_drugs[:3])}")

        print("\n" + "="*60)
        print("END OF REPORT")
        print("="*60)
    else:
        print(f"Error: {result.error}")


async def example_with_caching():
    """Example: Demonstrate caching"""
    print("\n" + "="*80)
    print("Example 6: Caching Demonstration")
    print("="*80)

    adapter = FDAFAERSAdapter()

    print("\nFirst call (will fetch from API)...")
    import time
    start = time.time()
    result1 = await adapter("lisinopril", use_cache=True)
    time1 = time.time() - start

    print(f"Time: {time1:.2f}s")
    print(f"Cache hit: {result1.cache_hit}")
    print(f"Events found: {result1.data['total_count']:,}" if result1.success else "Error")

    print("\nSecond call (should use cache)...")
    start = time.time()
    result2 = await adapter("lisinopril", use_cache=True)
    time2 = time.time() - start

    print(f"Time: {time2:.2f}s")
    print(f"Cache hit: {result2.cache_hit}")
    print(f"Events found: {result2.data['total_count']:,}" if result2.success else "Error")

    if time2 > 0:
        print(f"\nSpeedup: {time1/time2:.1f}x faster with cache")


async def main():
    """Run all examples"""
    print("\n" + "="*80)
    print("FDA FAERS Adapter - Example Usage")
    print("Adverse Event Reporting System Analysis")
    print("="*80)

    try:
        await example_basic_search()
        await asyncio.sleep(0.5)  # Rate limiting

        await example_specific_drug()
        await asyncio.sleep(0.5)

        await example_safety_comparison()
        await asyncio.sleep(0.5)

        await example_search_by_ingredient()
        await asyncio.sleep(0.5)

        await example_pharmacovigilance_report()
        await asyncio.sleep(0.5)

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
