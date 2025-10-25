"""
Test script to verify Celery integration works
Run this after starting docker-compose up
"""
import requests
import time
import json

API_BASE = "http://localhost:8000/api/v1"

def test_pipeline_run():
    """Test creating and monitoring a pipeline run"""
    print("=" * 60)
    print("Testing PharmForge Celery Pipeline Integration")
    print("=" * 60)

    # 1. Create a new run
    print("\n1. Creating new pipeline run...")
    create_data = {
        "input_smiles": ["CCO", "CC(=O)O"],  # Ethanol and Acetic acid
        "user_id": "test_user"
    }

    response = requests.post(f"{API_BASE}/runs", json=create_data)

    if response.status_code != 201:
        print(f"❌ Failed to create run: {response.status_code}")
        print(response.text)
        return

    run_data = response.json()
    run_id = run_data["run_id"]
    print(f"✓ Created run: {run_id}")
    print(f"  Status: {run_data['status']}")

    # 2. Poll for completion
    print(f"\n2. Monitoring run {run_id}...")
    max_polls = 60  # 5 minutes max
    poll_interval = 5  # seconds

    for i in range(max_polls):
        time.sleep(poll_interval)

        response = requests.get(f"{API_BASE}/runs/{run_id}")
        if response.status_code != 200:
            print(f"❌ Failed to get run status: {response.status_code}")
            return

        run_data = response.json()
        status = run_data["status"]

        print(f"  Poll {i+1}: Status = {status}")

        if status == "completed":
            print("\n✓ Pipeline completed successfully!")
            print(f"\nResults preview:")
            if run_data.get("results"):
                results = run_data["results"]
                print(f"  - Total compounds processed: {results.get('total_processed', 0)}")
                print(f"  - Adapter sequence: {results.get('adapter_sequence', [])}")

                if results.get("compounds"):
                    print(f"\n  First compound results:")
                    first = results["compounds"][0]
                    print(f"    SMILES: {first['smiles']}")
                    if first.get("properties"):
                        props = first["properties"]
                        print(f"    MW: {props.get('molecular_weight', 'N/A')}")
                        print(f"    LogP: {props.get('logp', 'N/A')}")
                        print(f"    TPSA: {props.get('tpsa', 'N/A')}")
            break

        elif status == "failed":
            print(f"\n❌ Pipeline failed!")
            print(f"Error: {run_data.get('error_message', 'Unknown error')}")
            break

        elif i == max_polls - 1:
            print(f"\n⚠ Timeout waiting for pipeline completion")

    # 3. List all runs
    print(f"\n3. Listing all runs...")
    response = requests.get(f"{API_BASE}/runs?page=1&page_size=5")

    if response.status_code != 200:
        print(f"❌ Failed to list runs: {response.status_code}")
        return

    list_data = response.json()
    print(f"✓ Found {list_data['total']} total runs")
    print(f"  Showing {len(list_data['runs'])} on page {list_data['page']}")

    for run in list_data["runs"]:
        print(f"  - {run['run_id']}: {run['status']} ({len(run['input_smiles'])} compounds)")

    print("\n" + "=" * 60)
    print("Test completed!")
    print("=" * 60)


if __name__ == "__main__":
    try:
        test_pipeline_run()
    except requests.exceptions.ConnectionError:
        print("❌ Cannot connect to API. Make sure docker-compose is running:")
        print("   docker-compose up -d")
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
