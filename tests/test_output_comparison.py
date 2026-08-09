"""
Output comparison tests for validating refactoring.

These tests compare current output against known reference outputs to ensure
refactoring doesn't change the visual output of diagrams.

Usage:
    # Run comparison tests (safe for normal test runs)
    pytest tests/test_output_comparison.py::TestOutputComparison -v

    # Intentionally update tracked references
    pytest tests/test_output_comparison.py::TestGenerateReferences \
        --update-reference-outputs -v
"""

import hashlib
import re
import shutil
import tempfile
from pathlib import Path

import pytest

from tests.utils.svg_compare import compare_svgs


# Reference output directory
TESTS_DIR = Path(__file__).parent
REFERENCE_DIR = TESTS_DIR / "reference_outputs"


def _strip_additive_semantic_attributes(svg_text: str) -> str:
    """Ignore non-visual metadata in geometry reference comparisons."""

    def strip_definition_group_metadata(match: re.Match[str]) -> str:
        tag = match.group(0)
        role_match = re.search(r'\sdata-gbdraw-role="([^"]*)"', tag)
        if role_match is None or role_match.group(1) not in {
            "record-definition",
            "record-definition-row",
            "plot-title",
        }:
            return tag
        for attribute in ("data-gbdraw-role", "data-gbdraw-definition-part"):
            tag = re.sub(rf'\s{attribute}="[^"]*"', "", tag)
        return tag

    svg_text = re.sub(r"<g\b[^>]*>", strip_definition_group_metadata, svg_text)
    for attribute in (
        "data-gbdraw-match-id",
        "data-gbdraw-record-id",
        "data-gbdraw-record-index",
        "data-gbdraw-slot-id",
        "data-gbdraw-slot-renderer",
    ):
        svg_text = re.sub(rf'\s{attribute}="[^"]*"', "", svg_text)
    return svg_text


def get_file_hash(path: Path) -> str:
    """Get SHA256 hash of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_reference_atomically(source: Path, destination: Path) -> None:
    """Replace a reference from a temporary file in the same directory."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            with source.open("rb") as source_handle:
                shutil.copyfileobj(source_handle, temporary)
        temporary_path.replace(destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


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
    "circular_radial_labels": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": [
            "--labels",
            "both",
            "--label_placement",
            "radial",
            "--label_rendering",
            "external_only",
            "--legend",
            "none",
        ],
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
        "args": ["--no-gc", "--legend", "none"],
    },
    "circular_no_skew": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--no-skew", "--legend", "none"],
    },
    "circular_no_gc_skew": {
        "type": "circular",
        "gbk": "MjeNMV.gb",
        "args": ["--no-gc", "--no-skew", "--legend", "none"],
    },
    "circular_hepatoplasma_repeat_underlay": {
        "type": "circular",
        "gbk": "AP027078.gb",
        "args": ["--no-gc", "--no-skew", "--legend", "none"],
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
        "args": ["--gc", "--skew", "--legend", "none"],
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
    "linear_hepatoplasma_repeat_underlay": {
        "type": "linear",
        "gbk": ["AP027078.gb"],
        "args": ["--legend", "none"],
    },
}


def _render_case(name, config, output_dir, find_test_input, gbdraw_runner) -> Path:
    input_names = config["gbk"]
    if isinstance(input_names, str):
        input_names = [input_names]
    gbk_files = [find_test_input(filename) for filename in input_names]
    if not all(gbk_files):
        pytest.skip(f"Some input files not found: {input_names}")

    blast_files = [find_test_input(filename) for filename in config.get("blast", ())]
    if not all(blast_files):
        pytest.skip(f"Some BLAST files not found: {config['blast']}")

    returncode, output, output_svg = gbdraw_runner.run(
        config["type"],
        gbk_files,
        name,
        output_dir,
        blast_files=blast_files,
        extra_args=config["args"],
    )
    if returncode != 0:
        raise RuntimeError(f"gbdraw failed: {output}")
    return output_svg


@pytest.mark.reference_generation
class TestGenerateReferences:
    """Generate reference outputs for comparison testing.

    Run only after approving an intentional SVG change.
    """

    @pytest.mark.parametrize("name,config", TEST_CASES.items())
    def test_generate_reference(
        self, name, config, tmp_path, find_test_input, gbdraw_runner
    ):
        """Generate a reference output for a test case."""
        output_svg = _render_case(
            name, config, tmp_path, find_test_input, gbdraw_runner
        )

        # Replace only after generation succeeds, without exposing a partial fixture.
        ref_path = REFERENCE_DIR / f"{name}.svg"
        replace_reference_atomically(output_svg, ref_path)

        print(f"Generated reference: {ref_path}")
        print(f"  Hash: {get_file_hash(ref_path)[:16]}...")


class TestOutputComparison:
    """Compare current output against reference outputs.

    Run this after each refactoring step to ensure no output changes.
    """

    @pytest.mark.parametrize("name,config", TEST_CASES.items())
    def test_compare_output(
        self, name, config, tmp_path, find_test_input, gbdraw_runner
    ):
        """Compare current output against reference."""
        ref_path = REFERENCE_DIR / f"{name}.svg"
        if not ref_path.exists():
            pytest.skip(
                f"Reference file not found: {ref_path}. Run TestGenerateReferences first."
            )

        output_svg = _render_case(
            name, config, tmp_path, find_test_input, gbdraw_runner
        )

        # Compare with reference
        result = compare_svgs(
            _strip_additive_semantic_attributes(ref_path.read_text(encoding="utf-8")),
            _strip_additive_semantic_attributes(output_svg.read_text(encoding="utf-8")),
        )

        if not result.equal:
            # Keep diagnostics outside the tracked reference directory.
            actual_path = tmp_path / f"{name}.actual.svg"
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

    def test_circular_produces_svg(self, tmp_path, find_test_input, gbdraw_runner):
        """Verify circular command produces valid SVG."""
        output = _render_case(
            "test",
            {"type": "circular", "gbk": "MjeNMV.gb", "args": ["--legend", "none"]},
            tmp_path,
            find_test_input,
            gbdraw_runner,
        )

        content = output.read_text()
        assert "<svg" in content
        assert "</svg>" in content
        assert "xmlns" in content

    def test_linear_produces_svg(self, tmp_path, find_test_input, gbdraw_runner):
        """Verify linear command produces valid SVG."""
        output = _render_case(
            "test",
            {"type": "linear", "gbk": ["MjeNMV.gb"], "args": ["--legend", "none"]},
            tmp_path,
            find_test_input,
            gbdraw_runner,
        )

        content = output.read_text()
        assert "<svg" in content
        assert "</svg>" in content
        assert "xmlns" in content
