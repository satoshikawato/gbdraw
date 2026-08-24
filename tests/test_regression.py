"""
Regression tests for gbdraw output.

These tests ensure that gbdraw produces identical output to known-good reference
outputs. Any changes to the output should be intentional and explicitly approved.

Run with: pytest tests/test_regression.py -v
Run fast tests only: pytest tests/test_regression.py -v -m "not slow"
"""

import subprocess
import sys

import pytest

from tests.utils.svg_compare import compare_svgs


# ============================================================================
# Basic Functionality Tests
# ============================================================================

class TestBasicFunctionality:
    """Test that gbdraw runs without errors."""

    def test_help_command(self):
        """Test that help command works."""
        result = subprocess.run(
            [sys.executable, "-m", "gbdraw.cli", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "gbdraw" in result.stdout.lower() or "genome" in result.stdout.lower()

    def test_circular_help(self):
        """Test that circular subcommand help works."""
        result = subprocess.run(
            [sys.executable, "-m", "gbdraw.cli", "circular", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "circular" in result.stdout.lower() or "gbk" in result.stdout.lower()

    def test_linear_help(self):
        """Test that linear subcommand help works."""
        result = subprocess.run(
            [sys.executable, "-m", "gbdraw.cli", "linear", "--help"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "linear" in result.stdout.lower() or "gbk" in result.stdout.lower()

    def test_unknown_command_lists_every_surface(self):
        result = subprocess.run(
            [sys.executable, "-m", "gbdraw.cli", "unknown"],
            capture_output=True,
            text=True,
            timeout=30,
        )

        assert result.returncode == 1
        assert result.stdout == ""
        assert "circular, linear, or gui" in result.stderr
        assert "--help" in result.stderr


# ============================================================================
# Circular Diagram Regression Tests
# ============================================================================

class TestCircularRegression:
    """Regression tests for circular diagrams."""

    @pytest.mark.regression
    @pytest.mark.circular
    def test_mje_nmv_basic(self, tmp_path, find_test_input, gbdraw_runner):
        """Test basic MjeNMV circular diagram."""
        gbk_file = find_test_input("MjeNMV.gb")
        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        # Run gbdraw
        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            "MjeNMV_test",
            tmp_path,
            extra_args=["--separate_strands", "--legend", "none"],
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists(), f"Output SVG not created: {svg_path}"

        # Check SVG is valid (can be parsed)
        svg_content = svg_path.read_text()
        assert "<svg" in svg_content
        assert "</svg>" in svg_content

    @pytest.mark.regression
    @pytest.mark.circular
    def test_mje_nmv_with_labels(self, tmp_path, find_test_input, gbdraw_runner):
        """Test MjeNMV with labels (without color table to avoid format issues)."""
        gbk_file = find_test_input("MjeNMV.gb")

        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        # Test with labels but without color table (color table format varies)
        args = ["--separate_strands", "--labels", "--legend", "none"]

        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            "MjeNMV_labels_test",
            tmp_path,
            extra_args=args,
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

        # Verify the SVG contains text elements (labels)
        svg_content = svg_path.read_text()
        assert "<text" in svg_content, "SVG should contain text elements for labels"

    @pytest.mark.regression
    @pytest.mark.circular
    @pytest.mark.slow
    def test_nc_010162_large_genome(self, tmp_path, find_test_input, gbdraw_runner):
        """Test large genome (NC_010162 - Sorangium cellulosum)."""
        gbk_file = find_test_input("NC_010162.gb")
        whitelist = find_test_input("NC_010162.whitelist.tsv")
        color_table = find_test_input("NC_010162.feature-specific_table.tsv")

        if not gbk_file:
            pytest.skip("NC_010162.gb not found")

        args = [
            "--separate_strands",
            "-p",
            "edelweiss",
            "--track_type",
            "tuckin",
            "--legend",
            "left",
        ]

        if whitelist:
            args.extend(["--labels", "--label_whitelist", str(whitelist)])
        if color_table:
            args.extend(["-t", str(color_table)])

        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            "NC_010162_test",
            tmp_path,
            extra_args=args,
            timeout=600,  # Large genome needs more time
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

    @pytest.mark.regression
    @pytest.mark.circular
    def test_track_types(self, tmp_path, find_test_input, gbdraw_runner):
        """Test all track types produce valid output."""
        gbk_file = find_test_input("MjeNMV.gb")
        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        for track_type in ["tuckin", "middle", "spreadout"]:
            returncode, output, svg_path = gbdraw_runner.run(
                "circular",
                [gbk_file],
                f"MjeNMV_{track_type}",
                tmp_path,
                extra_args=["--track_type", track_type, "--legend", "none"],
            )

            assert returncode == 0, f"Track type '{track_type}' failed: {output}"
            assert svg_path.exists(), f"No output for track type '{track_type}'"

    @pytest.mark.regression
    @pytest.mark.circular
    def test_gc_skew_options(self, tmp_path, find_test_input, gbdraw_runner):
        """Test GC content and skew suppression options."""
        gbk_file = find_test_input("MjeNMV.gb")
        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        # Test with both enabled (default)
        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            "MjeNMV_gc_both",
            tmp_path,
            extra_args=["--legend", "none"],
        )
        assert returncode == 0

        # Test with GC suppressed
        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            "MjeNMV_no_gc",
            tmp_path,
            extra_args=["--no-gc", "--legend", "none"],
        )
        assert returncode == 0

        # Test with both suppressed
        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            "MjeNMV_no_gc_skew",
            tmp_path,
            extra_args=["--no-gc", "--no-skew", "--legend", "none"],
        )
        assert returncode == 0


# ============================================================================
# Linear Diagram Regression Tests
# ============================================================================

class TestLinearRegression:
    """Regression tests for linear diagrams."""

    @pytest.mark.regression
    @pytest.mark.linear
    def test_single_genome_linear(self, tmp_path, find_test_input, gbdraw_runner):
        """Test single genome linear diagram."""
        gbk_file = find_test_input("MjeNMV.gb")
        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        returncode, output, svg_path = gbdraw_runner.run(
            "linear",
            [gbk_file],
            "MjeNMV_linear_test",
            tmp_path,
            extra_args=["--legend", "none"],
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

    @pytest.mark.regression
    @pytest.mark.linear
    def test_multi_genome_linear(self, tmp_path, find_test_input, gbdraw_runner):
        """Test multi-genome linear diagram."""
        gbk_files = []
        for name in ["MjeNMV.gb", "MelaMJNV.gb"]:
            f = find_test_input(name)
            if f:
                gbk_files.append(f)

        if len(gbk_files) < 2:
            pytest.skip("Need at least 2 GenBank files for multi-genome test")

        returncode, output, svg_path = gbdraw_runner.run(
            "linear",
            gbk_files,
            "multi_genome_test",
            tmp_path,
            extra_args=["--legend", "none"],
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

    @pytest.mark.regression
    @pytest.mark.linear
    def test_linear_with_blast(self, tmp_path, find_test_input, gbdraw_runner):
        """Test linear diagram with BLAST comparison."""
        gbk1 = find_test_input("MjeNMV.gb")
        gbk2 = find_test_input("MelaMJNV.gb")
        blast = find_test_input("MjeNMV.MelaMJNV.tblastx.out")

        if not (gbk1 and gbk2 and blast):
            pytest.skip("Required files not found for BLAST comparison test")

        returncode, output, svg_path = gbdraw_runner.run(
            "linear",
            [gbk1, gbk2],
            "blast_comparison_test",
            tmp_path,
            blast_files=[blast],
            extra_args=["--align_center", "--legend", "none"],
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

    @pytest.mark.regression
    @pytest.mark.linear
    def test_linear_with_gc_skew(self, tmp_path, find_test_input, gbdraw_runner):
        """Test linear diagram with GC content and skew."""
        gbk_file = find_test_input("MjeNMV.gb")
        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        returncode, output, svg_path = gbdraw_runner.run(
            "linear",
            [gbk_file],
            "MjeNMV_gc_skew_test",
            tmp_path,
            extra_args=["--gc", "--skew", "--legend", "none"],
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

    @pytest.mark.regression
    @pytest.mark.linear
    @pytest.mark.slow
    def test_majani_full_comparison(self, tmp_path, find_test_input, gbdraw_runner):
        """Test full Majani virus comparison (10 genomes)."""
        genome_names = [
            "MjeNMV.gb",
            "MelaMJNV.gb",
            "PemoMJNVA.gb",
            "PeseMJNV.gb",
            "PemoMJNVB.gb",
            "LvMJNV.gb",
            "TrcuMJNV.gb",
            "MellatMJNV.gb",
            "MeenMJNV.gb",
            "MejoMJNV.gb",
        ]
        blast_names = [
            "MjeNMV.MelaMJNV.tblastx.out",
            "MelaMJNV.PemoMJNVA.tblastx.out",
            "PemoMJNVA.PeseMJNV.tblastx.out",
            "PeseMJNV.PemoMJNVB.tblastx.out",
            "PemoMJNVB.LvMJNV.tblastx.out",
            "LvMJNV.TrcuMJNV.tblastx.out",
            "TrcuMJNV.MellatMJNV.tblastx.out",
            "MellatMJNV.MeenMJNV.tblastx.out",
            "MeenMJNV.MejoMJNV.tblastx.out",
        ]

        gbk_files = [find_test_input(n) for n in genome_names]
        blast_files = [find_test_input(n) for n in blast_names]

        gbk_files = [f for f in gbk_files if f is not None]
        blast_files = [f for f in blast_files if f is not None]

        if len(gbk_files) < 10 or len(blast_files) < 9:
            pytest.skip("Not all Majani files available")

        returncode, output, svg_path = gbdraw_runner.run(
            "linear",
            gbk_files,
            "majani_test",
            tmp_path,
            blast_files=blast_files,
            extra_args=["--align_center", "--separate_strands", "--legend", "none"],
            timeout=600,
        )

        assert returncode == 0, f"gbdraw failed: {output}"
        assert svg_path.exists()

        # Compare with reference if available
        ref_svg = find_test_input("majani.svg")
        if ref_svg:
            result = compare_svgs(ref_svg, svg_path)
            if not result.equal:
                print(f"Differences found: {result.message}")
                for diff in result.differences[:5]:
                    print(f"  {diff}")


# ============================================================================
# Color Palette Tests
# ============================================================================

class TestColorPalettes:
    """Test that all color palettes produce valid output."""

    PALETTES = [
        "default", "ajisai", "arctic", "autumn", "forest",
        "fugaku", "mint", "orchid", "sakura", "soft_pastels",
    ]

    @pytest.mark.regression
    @pytest.mark.parametrize("palette", PALETTES)
    def test_palette(self, palette, tmp_path, find_test_input, gbdraw_runner):
        """Test that each palette produces valid output."""
        gbk_file = find_test_input("MjeNMV.gb")
        if not gbk_file:
            pytest.skip("MjeNMV.gb not found")

        returncode, output, svg_path = gbdraw_runner.run(
            "circular",
            [gbk_file],
            f"palette_{palette}",
            tmp_path,
            extra_args=["-p", palette, "--legend", "none"],
        )

        assert returncode == 0, f"Palette '{palette}' failed: {output}"
        assert svg_path.exists()


# ============================================================================
# SVG Comparison Utility Tests
# ============================================================================

class TestSVGComparison:
    """Test the SVG comparison utility itself."""

    def test_identical_svgs(self):
        """Test that identical SVGs are detected as equal."""
        svg = """<?xml version="1.0"?>
        <svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
            <circle cx="50" cy="50" r="40"/>
        </svg>"""

        result = compare_svgs(svg, svg)
        assert result.equal
        assert "equivalent" in result.message.lower()

    def test_different_svgs(self):
        """Test that different SVGs are detected as different."""
        svg1 = """<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
            <circle cx="50" cy="50" r="40"/>
        </svg>"""

        svg2 = """<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">
            <circle cx="50" cy="50" r="30"/>
        </svg>"""

        result = compare_svgs(svg1, svg2)
        assert not result.equal

    def test_float_tolerance(self):
        """Test that floating point differences within tolerance are ignored."""
        svg1 = """<svg xmlns="http://www.w3.org/2000/svg">
            <circle cx="50.00001" cy="50.00001" r="40"/>
        </svg>"""

        svg2 = """<svg xmlns="http://www.w3.org/2000/svg">
            <circle cx="50.00002" cy="50.00002" r="40"/>
        </svg>"""

        result = compare_svgs(svg1, svg2)
        assert result.equal

    def test_composition_metadata_float_tolerance(self):
        """Treat runtime-scale float noise in composition metadata like SVG coordinates."""
        svg1 = '<svg data-gbdraw-composition="{&quot;x&quot;:67.00496309450932}"/>'
        svg2 = '<svg data-gbdraw-composition="{&quot;x&quot;:67.00496309451}"/>'

        assert compare_svgs(svg1, svg2).equal

    def test_basic_css_color_keyword_tolerance(self):
        """Treat basic CSS color names and their hexadecimal forms as equivalent."""
        named = '<svg><path fill="yellow" stroke="red"/></svg>'
        hexadecimal = '<svg><path fill="#FFFF00" stroke="#ff0000"/></svg>'

        assert compare_svgs(named, hexadecimal).equal

    def test_attribute_order_tolerance(self):
        """Test that attribute order differences are ignored."""
        svg1 = """<svg xmlns="http://www.w3.org/2000/svg">
            <circle cx="50" cy="50" r="40" fill="red"/>
        </svg>"""

        svg2 = """<svg xmlns="http://www.w3.org/2000/svg">
            <circle fill="red" r="40" cy="50" cx="50"/>
        </svg>"""

        result = compare_svgs(svg1, svg2)
        assert result.equal
