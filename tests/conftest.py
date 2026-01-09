"""
Pytest configuration and fixtures for gbdraw regression tests.

This module provides shared fixtures for test input files, reference outputs,
and helper functions for running gbdraw commands and comparing outputs.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Generator, Optional

import pytest

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "examples"
TESTS_DIR = PROJECT_ROOT / "tests"
REFERENCE_OUTPUTS_DIR = TESTS_DIR / "reference_outputs"
TEST_INPUTS_DIR = TESTS_DIR / "test_inputs"

# External test directories (if available)
EXTERNAL_TEST_DIR_1 = Path("/home/kawato/study/2025-09-17_gbdraw_test")
EXTERNAL_TEST_DIR_2 = Path("/home/kawato/study/2025-05-15_gbdraw")


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "regression: marks tests as regression tests"
    )
    config.addinivalue_line(
        "markers", "circular: marks tests for circular diagrams"
    )
    config.addinivalue_line(
        "markers", "linear: marks tests for linear diagrams"
    )


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Return the project root directory."""
    return PROJECT_ROOT


@pytest.fixture(scope="session")
def examples_dir() -> Path:
    """Return the examples directory."""
    return EXAMPLES_DIR


@pytest.fixture(scope="session")
def reference_outputs_dir() -> Path:
    """Return the reference outputs directory."""
    return REFERENCE_OUTPUTS_DIR


@pytest.fixture(scope="session")
def test_inputs_dir() -> Path:
    """Return the test inputs directory."""
    return TEST_INPUTS_DIR


@pytest.fixture(scope="session")
def external_test_dir_1() -> Optional[Path]:
    """Return the first external test directory if available."""
    if EXTERNAL_TEST_DIR_1.exists():
        return EXTERNAL_TEST_DIR_1
    return None


@pytest.fixture(scope="session")
def external_test_dir_2() -> Optional[Path]:
    """Return the second external test directory if available."""
    if EXTERNAL_TEST_DIR_2.exists():
        return EXTERNAL_TEST_DIR_2
    return None


@pytest.fixture(scope="function")
def temp_output_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture(scope="session")
def gbdraw_command() -> list[str]:
    """Return the command to invoke gbdraw."""
    return [sys.executable, "-m", "gbdraw.cli"]


class GbdrawRunner:
    """Helper class to run gbdraw commands and capture output."""

    def __init__(self, base_command: list[str]):
        self.base_command = base_command

    def run_circular(
        self,
        gbk_files: list[Path],
        output_prefix: str,
        output_dir: Path,
        extra_args: Optional[list[str]] = None,
        timeout: int = 120,
    ) -> tuple[int, str, Path]:
        """
        Run gbdraw circular command.

        Args:
            gbk_files: List of GenBank files
            output_prefix: Output file prefix
            output_dir: Directory for output files
            extra_args: Additional command-line arguments
            timeout: Command timeout in seconds

        Returns:
            Tuple of (return_code, stdout+stderr, output_svg_path)
        """
        cmd = self.base_command + ["circular"]
        cmd.extend(["--gbk"] + [str(f) for f in gbk_files])
        cmd.extend(["-o", str(output_dir / output_prefix)])
        cmd.extend(["-f", "svg"])

        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(output_dir),
        )

        output_svg = output_dir / f"{output_prefix}.svg"
        return result.returncode, result.stdout + result.stderr, output_svg

    def run_linear(
        self,
        gbk_files: list[Path],
        output_prefix: str,
        output_dir: Path,
        blast_files: Optional[list[Path]] = None,
        extra_args: Optional[list[str]] = None,
        timeout: int = 120,
    ) -> tuple[int, str, Path]:
        """
        Run gbdraw linear command.

        Args:
            gbk_files: List of GenBank files
            output_prefix: Output file prefix
            output_dir: Directory for output files
            blast_files: Optional BLAST comparison files
            extra_args: Additional command-line arguments
            timeout: Command timeout in seconds

        Returns:
            Tuple of (return_code, stdout+stderr, output_svg_path)
        """
        cmd = self.base_command + ["linear"]
        cmd.extend(["--gbk"] + [str(f) for f in gbk_files])
        cmd.extend(["-o", str(output_dir / output_prefix)])
        cmd.extend(["-f", "svg"])

        if blast_files:
            cmd.extend(["-b"] + [str(f) for f in blast_files])

        if extra_args:
            cmd.extend(extra_args)

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(output_dir),
        )

        output_svg = output_dir / f"{output_prefix}.svg"
        return result.returncode, result.stdout + result.stderr, output_svg


@pytest.fixture(scope="session")
def gbdraw_runner(gbdraw_command) -> GbdrawRunner:
    """Return a GbdrawRunner instance."""
    return GbdrawRunner(gbdraw_command)


# Test case definitions for regression testing
# Each entry defines: (name, type, inputs, args, reference_output)

CIRCULAR_TEST_CASES = [
    # Basic circular diagrams
    {
        "name": "mje_nmv_basic",
        "gbk": "MjeNMV.gb",
        "args": ["--separate_strands"],
        "reference": "MjeNMV_basic.svg",
        "description": "Basic circular diagram with separate strands",
    },
    {
        "name": "mje_nmv_labels",
        "gbk": "MjeNMV.gb",
        "args": ["--separate_strands", "--show_labels", "-t", "feature_specific_color_table.txt"],
        "reference": "MjeNMV_feature_specifc_colors_with_labels.svg",
        "description": "Circular diagram with labels and feature-specific colors",
    },
    {
        "name": "nc_010162_edelweiss",
        "gbk": "NC_010162.gb",
        "args": [
            "--separate_strands", "-p", "edelweiss", "--track_type", "tuckin",
            "--show_labels", "--label_whitelist", "NC_010162.whitelist.tsv",
            "-t", "NC_010162.feature-specific_table.tsv", "--legend", "left"
        ],
        "reference": "NC_010162_edelweiss.svg",
        "description": "Large genome with whitelist labels and custom colors",
    },
]

LINEAR_TEST_CASES = [
    # Basic linear diagrams
    {
        "name": "majani_linear",
        "gbk": ["MjeNMV.gb", "MelaMJNV.gb", "PemoMJNVA.gb", "PeseMJNV.gb",
                "PemoMJNVB.gb", "LvMJNV.gb", "TrcuMJNV.gb", "MellatMJNV.gb", "MeenMJNV.gb", "MejoMJNV.gb"],
        "blast": [
            "MjeNMV.MelaMJNV.tblastx.out", "MelaMJNV.PemoMJNVA.tblastx.out",
            "PemoMJNVA.PeseMJNV.tblastx.out", "PeseMJNV.PemoMJNVB.tblastx.out",
            "PemoMJNVB.LvMJNV.tblastx.out", "LvMJNV.TrcuMJNV.tblastx.out",
            "TrcuMJNV.MellatMJNV.tblastx.out", "MellatMJNV.MeenMJNV.tblastx.out",
            "MeenMJNV.MejoMJNV.tblastx.out"
        ],
        "args": ["--align_center", "--separate_strands"],
        "reference": "majani.svg",
        "description": "Multi-genome linear comparison",
    },
]


@pytest.fixture(scope="session")
def circular_test_cases() -> list[dict]:
    """Return circular test case definitions."""
    return CIRCULAR_TEST_CASES


@pytest.fixture(scope="session")
def linear_test_cases() -> list[dict]:
    """Return linear test case definitions."""
    return LINEAR_TEST_CASES


def get_test_input_path(filename: str) -> Optional[Path]:
    """
    Find a test input file across available directories.

    Args:
        filename: Name of the file to find

    Returns:
        Path to the file if found, None otherwise
    """
    # Check in order: test_inputs, examples, external dirs
    search_dirs = [
        TEST_INPUTS_DIR,
        EXAMPLES_DIR,
        EXTERNAL_TEST_DIR_1,
        EXTERNAL_TEST_DIR_2,
    ]

    for dir_path in search_dirs:
        if dir_path and dir_path.exists():
            file_path = dir_path / filename
            if file_path.exists():
                return file_path

    return None


def get_reference_output_path(filename: str) -> Optional[Path]:
    """
    Find a reference output file across available directories.

    Args:
        filename: Name of the file to find

    Returns:
        Path to the file if found, None otherwise
    """
    # Check in order: reference_outputs, examples, external dirs
    search_dirs = [
        REFERENCE_OUTPUTS_DIR,
        EXAMPLES_DIR,
        EXTERNAL_TEST_DIR_1,
        EXTERNAL_TEST_DIR_2,
    ]

    for dir_path in search_dirs:
        if dir_path and dir_path.exists():
            file_path = dir_path / filename
            if file_path.exists():
                return file_path

    return None


@pytest.fixture(scope="session")
def find_test_input():
    """Return a function to find test input files."""
    return get_test_input_path


@pytest.fixture(scope="session")
def find_reference_output():
    """Return a function to find reference output files."""
    return get_reference_output_path
