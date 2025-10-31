"""
Example usage of the py3Dmol Adapter
Demonstrates various use cases for 3D molecular visualization
"""
import asyncio
from adapters.py3dmol import Py3DmolAdapter


async def example_1_simple_smiles():
    """Example 1: Visualize a simple molecule from SMILES"""
    print("\n=== Example 1: Simple SMILES Visualization ===")

    adapter = Py3DmolAdapter()

    # Aspirin
    result = await adapter.execute(
        {"structure": "CC(=O)Oc1ccccc1C(=O)O"},
        style="stick",
        width=500,
        height=400
    )

    if result.success:
        print(f"Structure type: {result.data['structure_type']}")
        print(f"Format: {result.data['structure_format']}")
        print(f"HTML length: {len(result.data['html'])} characters")
        print("Visualization HTML generated successfully!")

        # Save to file for testing
        with open("aspirin_visualization.html", "w") as f:
            f.write(result.data['html'])
        print("Saved to: aspirin_visualization.html")
    else:
        print(f"Error: {result.error}")


async def example_2_drug_molecules():
    """Example 2: Visualize multiple drug molecules"""
    print("\n=== Example 2: Drug Molecules ===")

    adapter = Py3DmolAdapter()

    drugs = {
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Paracetamol": "CC(=O)Nc1ccc(O)cc1"
    }

    for name, smiles in drugs.items():
        result = await adapter.execute(
            {"structure": smiles},
            style="sphere",
            width=400,
            height=300,
            background_color="#f5f5f5"
        )

        if result.success:
            print(f"{name}: {result.data['structure_type']} ({len(result.data['html'])} chars)")

            # Save individual visualizations
            filename = f"{name.lower()}_viz.html"
            with open(filename, "w") as f:
                f.write(result.data['html'])
            print(f"  Saved to: {filename}")
        else:
            print(f"{name}: Error - {result.error}")


async def example_3_pdb_file():
    """Example 3: Visualize from PDB file"""
    print("\n=== Example 3: PDB File Visualization ===")

    adapter = Py3DmolAdapter()

    # Example with a PDB string (simulating file content)
    pdb_string = """
ATOM      1  N   ALA A   1      -8.778   1.234   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      -7.483   0.573   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1      -7.365  -0.660  -0.897  1.00  0.00           C
ATOM      4  O   ALA A   1      -8.349  -1.357  -1.124  1.00  0.00           O
ATOM      5  CB  ALA A   1      -7.166   0.193   1.443  1.00  0.00           C
"""

    result = await adapter.execute(
        {"structure": pdb_string},
        style="cartoon",
        width=600,
        height=600
    )

    if result.success:
        print(f"Structure type: {result.data['structure_type']}")
        print(f"Format: {result.data['structure_format']}")
        print(f"Viewer ID: {result.data['viewer_id']}")
        print("PDB visualization generated successfully!")
    else:
        print(f"Error: {result.error}")


async def example_4_different_styles():
    """Example 4: Compare different visual styles"""
    print("\n=== Example 4: Different Visual Styles ===")

    adapter = Py3DmolAdapter()
    smiles = "c1ccccc1"  # Benzene

    styles = ["stick", "sphere", "line", "cross"]

    for style in styles:
        result = await adapter.execute(
            {"structure": smiles},
            style=style,
            width=300,
            height=300
        )

        if result.success:
            print(f"{style.capitalize()} style: Success")

            # Save with style-specific filename
            with open(f"benzene_{style}.html", "w") as f:
                f.write(result.data['html'])
            print(f"  Saved to: benzene_{style}.html")
        else:
            print(f"{style.capitalize()} style: {result.error}")


async def example_5_complex_molecule():
    """Example 5: Visualize a complex molecule"""
    print("\n=== Example 5: Complex Molecule ===")

    adapter = Py3DmolAdapter()

    # Penicillin G
    penicillin_smiles = "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"

    result = await adapter.execute(
        {"structure": penicillin_smiles},
        style="stick",
        width=700,
        height=500,
        background_color="white"
    )

    if result.success:
        print(f"Penicillin G visualization:")
        print(f"  Structure type: {result.data['structure_type']}")
        print(f"  Format: {result.data['structure_format']}")
        print(f"  Dimensions: {result.data['width']}x{result.data['height']}")

        # Save with full HTML structure
        html_full = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Penicillin G - 3D Visualization</title>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f0f0f0;
        }}
        .container {{
            max-width: 800px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            text-align: center;
        }}
        .info {{
            margin-top: 20px;
            padding: 10px;
            background-color: #e8f4f8;
            border-radius: 5px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Penicillin G - 3D Molecular Structure</h1>
        <div style="text-align: center;">
            {result.data['html']}
        </div>
        <div class="info">
            <h3>Molecule Information</h3>
            <p><strong>SMILES:</strong> {penicillin_smiles}</p>
            <p><strong>Structure Type:</strong> {result.data['structure_type']}</p>
            <p><strong>Visualization Style:</strong> {result.data['style']}</p>
        </div>
    </div>
</body>
</html>
"""
        with open("penicillin_full.html", "w") as f:
            f.write(html_full)
        print("  Saved to: penicillin_full.html")
    else:
        print(f"Error: {result.error}")


async def example_6_error_handling():
    """Example 6: Error handling"""
    print("\n=== Example 6: Error Handling ===")

    adapter = Py3DmolAdapter()

    # Test invalid inputs
    test_cases = [
        ({"structure": ""}, "Empty structure"),
        ({"structure": "INVALID_SMILES_###"}, "Invalid SMILES"),
        ({"invalid_key": "data"}, "Missing structure key"),
        ("not_a_dict", "Not a dictionary")
    ]

    for input_data, description in test_cases:
        result = await adapter.execute(input_data)
        status = "Success" if result.success else f"Error: {result.error}"
        print(f"{description}: {status}")


async def example_7_metadata():
    """Example 7: Accessing metadata"""
    print("\n=== Example 7: Metadata Information ===")

    adapter = Py3DmolAdapter()

    # Get adapter metadata
    metadata = adapter.get_metadata()
    print("Adapter Metadata:")
    print(f"  Name: {metadata['name']}")
    print(f"  Type: {metadata['type']}")
    print(f"  Version: {metadata['version']}")
    print(f"  Enabled: {metadata['enabled']}")
    print(f"  Supported styles: {metadata['config']['supported_styles']}")

    # Execute and get result metadata
    result = await adapter.execute(
        {"structure": "CCO"},  # Ethanol
        style="stick"
    )

    if result.success:
        print("\nResult Metadata:")
        for key, value in result.metadata.items():
            print(f"  {key}: {value}")


async def main():
    """Run all examples"""
    print("=" * 60)
    print("py3Dmol Adapter - Usage Examples")
    print("=" * 60)

    await example_1_simple_smiles()
    await example_2_drug_molecules()
    await example_3_pdb_file()
    await example_4_different_styles()
    await example_5_complex_molecule()
    await example_6_error_handling()
    await example_7_metadata()

    print("\n" + "=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    asyncio.run(main())
