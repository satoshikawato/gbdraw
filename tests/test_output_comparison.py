"""
Output comparison tests for validating refactoring.

These tests compare current output against known reference outputs to ensure
refactoring doesn't change the visual output of diagrams.

Usage:
    # Generate reference outputs first (run once before refactoring)
    pytest tests/test_output_comparison.py::TestGenerateReferences -v

    # Run comparison tests (run after each refactoring step)
    pytest tests/test_output_comparison.py::TestOutputComparison -v
"""

import hashlib
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pytest

from tests.utils.svg_compare import compare_svgs


# Reference output directory
TESTS_DIR = Path(__file__).parent
REFERENCE_DIR = TESTS_DIR / "reference_outputs"
EXAMPLES_DIR = TESTS_DIR.parent / "examples"


def get_file_hash(path: Path) -> str:
    """Get SHA256 hash of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_input(filename: str) -> Optional[Path]:
    """Find an input file in examples or external directories."""
    search_dirs = [
        EXAMPLES_DIR,
        Path("/home/kawato/study/2025-09-17_gbdraw_test"),
        Path("/home/kawato/study/2025-05-15_gbdraw"),
    ]
    for d in search_dirs:
        if d.exists():
            f = d / filename
            if f.exists():
                return f
    return None


def run_gbdraw_circular(gbk_file: Path, output_dir: Path, output_name: str, extra_args: list[str]) -> Path:
    """Run gbdraw circular and return output path."""
    cmd = [
        sys.executable, "-m", "gbdraw.cli", "circular",
        "--gbk", str(gbk_file),
        "-o", output_name,
        "-f", "svg",
    ] + extra_args

    result = subprocess.run(cmd, cwd=str(output_dir), capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise RuntimeError(f"gbdraw failed: {result.stderr}")

    return output_dir / f"{output_name}.svg"


def run_gbdraw_linear(gbk_files: list[Path], output_dir: Path, output_name: str,
                      extra_args: list[str], blast_files: list[Path] = None) -> Path:
    """Run gbdraw linear and return output path."""
    cmd = [
        sys.executable, "-m", "gbdraw.cli", "linear",
        "--gbk", *[str(f) for f in gbk_files],
        "-o", output_name,
        "-f", "svg",
    ]
    if blast_files:
        cmd.extend(["-b", *[str(f) for f in blast_files]])
    cmd.extend(extra_args)

    result = subprocess.run(cmd, cwd=str(output_dir), capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        raise RuntimeError(f"gbdraw failed: {result.stderr}")

    return output_dir / f"{output_name}.svg"


# Test case definitions - these define the canonical test cases for regression testing
TEST_CASES = {
    # Circular test cases
    "circular_basic": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--separate_strands", "--legend", "none"],
    },
    "circular_with_labels": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--separate_strands", "--labels", "--legend", "none"],
    },
    "circular_tuckin": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--track_type", "tuckin", "--legend", "none"],
    },
    "circular_middle": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--track_type", "middle", "--legend", "none"],
    },
    "circular_spreadout": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--track_type", "spreadout", "--legend", "none"],
    },
    "circular_no_gc": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--suppress_gc", "--legend", "none"],
    },
    "circular_no_skew": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--suppress_skew", "--legend", "none"],
    },
    "circular_no_gc_skew": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--suppress_gc", "--suppress_skew", "--legend", "none"],
    },
    # Linear test cases
    "linear_basic": {
        "type": "linear",
        "gbk": ["MjeNMV.gb"],
        "args": ["--legend", "none"],
    },
    "linear_with_gc_skew": {
        "type": "linear",
        "gbk": ["MjeNMV.gb"],
        "args": ["--show_gc", "--show_skew", "--legend", "none"],
    },
    "linear_separate_strands": {
        "type": "linear",
        "gbk": ["MjeNMV.gb"],
        "args": ["--separate_strands", "--legend", "none"],
    },
    "linear_multi_genome": {
        "type": "linear",
        "gbk": ["MjeNMV.gb", "MelaMJNV.gb"],
        "args": ["--legend", "none"],
    },
    "linear_with_blast": {
        "type": "linear",
        "gbk": ["MjeNMV.gb", "MelaMJNV.gb"],
        "blast": ["MjeNMV.MelaMJNV.tblastx.out"],
        "args": ["--align_center", "--legend", "none"],
    },
}


class TestGenerateReferences:
    """Generate reference outputs for comparison testing.

    Run this ONCE before starting refactoring to create baseline outputs.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        """Ensure reference directory exists."""
        REFERENCE_DIR.mkdir(exist_ok=True)

    @pytest.mark.parametrize("name,config", TEST_CASES.items())
    def test_generate_reference(self, name, config, tmp_path):
        """Generate a reference output for a test case."""
        if config["type"] == "circular":
            gbk = find_input(config["gbk"])
            if not gbk:
                pytest.skip(f"Input file {config['gbk']} not found")

            output_svg = run_gbdraw_circular(gbk, tmp_path, name, config["args"])
        else:
            gbk_files = [find_input(g) for g in config["gbk"]]
            if not all(gbk_files):
                pytest.skip(f"Some input files not found: {config['gbk']}")

            blast_files = None
            if "blast" in config:
                blast_files = [find_input(b) for b in config["blast"]]
                if not all(blast_files):
                    pytest.skip(f"BLAST files not found: {config['blast']}")

            output_svg = run_gbdraw_linear(gbk_files, tmp_path, name, config["args"], blast_files)

        # Copy to reference directory
        ref_path = REFERENCE_DIR / f"{name}.svg"
        shutil.copy(output_svg, ref_path)

        print(f"Generated reference: {ref_path}")
        print(f"  Hash: {get_file_hash(ref_path)[:16]}...")


class TestOutputComparison:
    """Compare current output against reference outputs.

    Run this after each refactoring step to ensure no output changes.
    """

    @pytest.mark.parametrize("name,config", TEST_CASES.items())
    def test_compare_output(self, name, config, tmp_path):
        """Compare current output against reference."""
        ref_path = REFERENCE_DIR / f"{name}.svg"
        if not ref_path.exists():
            pytest.skip(f"Reference file not found: {ref_path}. Run TestGenerateReferences first.")

        if config["type"] == "circular":
            gbk = find_input(config["gbk"])
            if not gbk:
                pytest.skip(f"Input file {config['gbk']} not found")

            output_svg = run_gbdraw_circular(gbk, tmp_path, name, config["args"])
        else:
            gbk_files = [find_input(g) for g in config["gbk"]]
            if not all(gbk_files):
                pytest.skip(f"Some input files not found: {config['gbk']}")

            blast_files = None
            if "blast" in config:
                blast_files = [find_input(b) for b in config["blast"]]
                if not all(blast_files):
                    pytest.skip(f"BLAST files not found: {config['blast']}")

            output_svg = run_gbdraw_linear(gbk_files, tmp_path, name, config["args"], blast_files)

        # Compare with reference
        result = compare_svgs(ref_path, output_svg)

        if not result.equal:
            # Save the actual output for inspection
            actual_path = REFERENCE_DIR / f"{name}.actual.svg"
            shutil.copy(output_svg, actual_path)

            diff_msg = "\n".join(result.differences[:10])
            pytest.fail(
                f"Output differs from reference!\n"
                f"Reference: {ref_path}\n"
                f"Actual: {actual_path}\n"
                f"Differences:\n{diff_msg}"
            )


class TestQuickValidation:
    """Quick validation that gbdraw produces valid output.

    These tests don't compare against references, just verify output is valid SVG.
    """

    def test_circular_produces_svg(self, tmp_path):
        """Verify circular command produces valid SVG."""
        gbk = find_input("MjeNMV.gb")
        if not gbk:
            pytest.skip("MjeNMV.gb not found")

        output = run_gbdraw_circular(gbk, tmp_path, "test", ["--legend", "none"])

        content = output.read_text()
        assert "<svg" in content
        assert "</svg>" in content
        assert "xmlns" in content

    def test_linear_produces_svg(self, tmp_path):
        """Verify linear command produces valid SVG."""
        gbk = find_input("MjeNMV.gb")
        if not gbk:
            pytest.skip("MjeNMV.gb not found")

        output = run_gbdraw_linear([gbk], tmp_path, "test", ["--legend", "none"])

        content = output.read_text()
        assert "<svg" in content
        assert "</svg>" in content
        assert "xmlns" in content
