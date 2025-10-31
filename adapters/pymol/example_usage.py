"""
PyMOL Adapter - Example Usage
Demonstrates various visualization workflows with the PyMOL adapter
"""
import asyncio
import os
from pathlib import Path
from adapters.pymol import PyMOLAdapter


async def example_basic_smiles():
    """Example 1: Basic SMILES visualization"""
    print("\n=== Example 1: Basic SMILES Visualization ===")

    adapter = PyMOLAdapter()

    # Aspirin molecule
    smiles = "CC(=O)Oc1ccccc1C(=O)O"

    result = await adapter.execute(
        {"structure": smiles},
        style="sticks",
        color_scheme="by_element",
        ray_trace=True,
        width=800,
        height=600
    )

    if result.success:
        print(f"✓ Image saved to: {result.data['image_path']}")
        print(f"  Render time: {result.data['render_time']}s")
        print(f"  Image size: {result.data['file_size']} bytes")
        print(f"  Atoms: {result.data['num_atoms']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_multiple_styles():
    """Example 2: Compare different visualization styles"""
    print("\n=== Example 2: Multiple Visualization Styles ===")

    adapter = PyMOLAdapter()

    # Caffeine molecule
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    styles = ["sticks", "spheres", "lines", "surface"]

    for style in styles:
        result = await adapter.execute(
            {"structure": smiles},
            style=style,
            color_scheme="by_element",
            ray_trace=False,  # Fast preview
            width=600,
            height=600,
            output_path=f"caffeine_{style}.png"
        )

        if result.success:
            print(f"✓ {style.capitalize()}: {result.data['image_path']}")
        else:
            print(f"✗ {style.capitalize()}: {result.error}")


async def example_color_schemes():
    """Example 3: Different color schemes"""
    print("\n=== Example 3: Color Schemes ===")

    adapter = PyMOLAdapter()

    smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen

    color_schemes = ["by_element", "rainbow", "spectrum", "cyan"]

    for scheme in color_schemes:
        result = await adapter.execute(
            {"structure": smiles},
            style="sticks",
            color_scheme=scheme,
            ray_trace=True,
            width=800,
            height=600,
            output_path=f"ibuprofen_{scheme}.png"
        )

        if result.success:
            print(f"✓ {scheme}: {result.data['image_path']}")
        else:
            print(f"✗ {scheme}: {result.error}")


async def example_ray_tracing_comparison():
    """Example 4: Ray tracing quality comparison"""
    print("\n=== Example 4: Ray Tracing Comparison ===")

    adapter = PyMOLAdapter()

    smiles = "CC(=O)Nc1ccc(cc1)O"  # Paracetamol

    # Without ray tracing (fast)
    result_fast = await adapter.execute(
        {"structure": smiles},
        ray_trace=False,
        output_path="paracetamol_fast.png"
    )

    # With ray tracing (high quality)
    result_quality = await adapter.execute(
        {"structure": smiles},
        ray_trace=True,
        output_path="paracetamol_quality.png"
    )

    if result_fast.success and result_quality.success:
        print(f"✓ Fast render: {result_fast.data['render_time']}s")
        print(f"✓ Quality render: {result_quality.data['render_time']}s")
        print(f"  Speedup: {result_quality.data['render_time'] / result_fast.data['render_time']:.1f}x faster without ray tracing")
    else:
        print("✗ Error in rendering")


async def example_batch_processing():
    """Example 5: Batch processing multiple molecules"""
    print("\n=== Example 5: Batch Processing ===")

    adapter = PyMOLAdapter()

    # Common pain relief drugs
    drugs = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "Paracetamol": "CC(=O)Nc1ccc(cc1)O",
        "Naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O"
    }

    structures = [{"structure": smiles} for smiles in drugs.values()]

    results = adapter.render_batch(
        structures,
        style="sticks",
        color_scheme="by_element",
        ray_trace=True,
        width=600,
        height=600
    )

    for name, result in zip(drugs.keys(), results):
        if result.success:
            print(f"✓ {name}: {result.data['render_time']}s")
        else:
            print(f"✗ {name}: {result.error}")


async def example_high_resolution():
    """Example 6: High-resolution publication-quality image"""
    print("\n=== Example 6: High-Resolution Render ===")

    adapter = PyMOLAdapter()

    # Complex molecule - Doxorubicin (anticancer drug)
    smiles = "CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)(C(=O)CO)O)N)O"

    result = await adapter.execute(
        {"structure": smiles},
        style="sticks",
        color_scheme="by_element",
        background="white",
        ray_trace=True,
        width=1920,   # Full HD width
        height=1080,  # Full HD height
        output_path="doxorubicin_hd.png"
    )

    if result.success:
        print(f"✓ High-res image: {result.data['image_path']}")
        print(f"  Resolution: {result.data['width']}x{result.data['height']}")
        print(f"  File size: {result.data['file_size'] / 1024:.1f} KB")
        print(f"  Render time: {result.data['render_time']}s")
    else:
        print(f"✗ Error: {result.error}")


async def example_pdb_file():
    """Example 7: Loading from PDB file"""
    print("\n=== Example 7: PDB File Visualization ===")

    adapter = PyMOLAdapter()

    # This example assumes you have a PDB file
    # For demo purposes, we'll create a simple PDB string
    pdb_content = """ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00  0.00           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00  0.00           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00  0.00           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00  0.00           O
ATOM      5  N   ALA A   2      -7.656   2.923   3.155  1.00  0.00           N
ATOM      6  CA  ALA A   2      -6.522   2.038   2.831  1.00  0.00           C
ATOM      7  C   ALA A   2      -5.241   2.537   3.427  1.00  0.00           C
ATOM      8  O   ALA A   2      -4.978   3.742   3.426  1.00  0.00           O
END"""

    result = await adapter.execute(
        {"structure": pdb_content},
        format="pdb",
        style="cartoon",
        color_scheme="rainbow",
        ray_trace=True,
        width=800,
        height=600
    )

    if result.success:
        print(f"✓ PDB visualization: {result.data['image_path']}")
        print(f"  Style: {result.data['style']}")
    else:
        print(f"✗ Error: {result.error}")


async def example_custom_backgrounds():
    """Example 8: Different background colors"""
    print("\n=== Example 8: Custom Backgrounds ===")

    adapter = PyMOLAdapter()

    smiles = "c1ccccc1"  # Benzene

    backgrounds = ["white", "black", "gray"]

    for bg in backgrounds:
        result = await adapter.execute(
            {"structure": smiles},
            style="sticks",
            color_scheme="by_element",
            background=bg,
            ray_trace=True,
            width=600,
            height=600,
            output_path=f"benzene_bg_{bg}.png"
        )

        if result.success:
            print(f"✓ Background {bg}: {result.data['image_path']}")
        else:
            print(f"✗ Background {bg}: {result.error}")


async def example_error_handling():
    """Example 9: Error handling"""
    print("\n=== Example 9: Error Handling ===")

    adapter = PyMOLAdapter()

    # Invalid SMILES
    result = await adapter.execute(
        {"structure": "INVALID_SMILES_STRING"},
        style="sticks"
    )

    if not result.success:
        print(f"✓ Error correctly handled: {result.error}")

    # Missing structure key
    result = await adapter.execute(
        {"invalid_key": "value"}
    )

    if not result.success:
        print(f"✓ Validation error handled: {result.error}")


async def example_cache_usage():
    """Example 10: Caching behavior"""
    print("\n=== Example 10: Caching ===")

    adapter = PyMOLAdapter()

    smiles = "CCO"  # Ethanol

    # First render (no cache)
    result1 = await adapter(
        {"structure": smiles},
        use_cache=True
    )

    if result1.success:
        print(f"✓ First render: {result1.data['render_time']}s (cache_hit: {result1.cache_hit})")

    # Second render (should hit cache)
    result2 = await adapter(
        {"structure": smiles},
        use_cache=True
    )

    if result2.success:
        print(f"✓ Second render: cached (cache_hit: {result2.cache_hit})")
        if result2.cache_hit:
            print("  → Cache working correctly!")


async def main():
    """Run all examples"""
    print("=" * 60)
    print("PyMOL Adapter - Example Usage")
    print("=" * 60)

    try:
        # Run examples
        await example_basic_smiles()
        await example_multiple_styles()
        await example_color_schemes()
        await example_ray_tracing_comparison()
        await example_batch_processing()
        await example_high_resolution()
        await example_pdb_file()
        await example_custom_backgrounds()
        await example_error_handling()
        await example_cache_usage()

        print("\n" + "=" * 60)
        print("All examples completed!")
        print("=" * 60)

    except Exception as e:
        print(f"\n✗ Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # Run the examples
    asyncio.run(main())
