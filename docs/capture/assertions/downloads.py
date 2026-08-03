"""Mechanical checks for documentation-owned downloads."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from assertions.svg_semantics import (
    assert_finished_circular_svg,
    assert_finished_linear_svg,
    assert_first_linear_svg,
    assert_gui_circular_layout_svg,
    assert_gui_losatn_svg,
    inspect_svg_file,
)


EXPECTED_FIRST_CIRCULAR_SVG = "human_mitochondrion.svg"
EXPECTED_FIRST_LINEAR_SVG = "lambda_linear.svg"
EXPECTED_GUI_INPUTS_SVG = "lambda_gff3.svg"
EXPECTED_GUI_LOSATN_TSV = "lambda-de3.losatn.tsv"
EXPECTED_GUI_LOSATN_SVG = "lambda-de3-losatn.svg"
EXPECTED_GUI_CIRCULAR_LAYOUT_SVG = "multi_record_circular.svg"
EXPECTED_GUI_LOSATN_TSV_SHA256 = (
    "703b0ac749669152a8e6d5fa6fb246cf5973cb1a3e0ca9db2f346ee33628317c"
)
MINIMUM_STATIC_SVG_BYTES = 10_000


def assert_first_circular_svg_download(path: Path) -> dict[str, Any]:
    """Validate the promised static SVG filename, XML, semantics, and safety."""

    if path.name != EXPECTED_FIRST_CIRCULAR_SVG:
        raise AssertionError(
            f"Expected {EXPECTED_FIRST_CIRCULAR_SVG}, downloaded {path.name}"
        )
    size = path.stat().st_size
    if size < MINIMUM_STATIC_SVG_BYTES:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {size} bytes")

    report = inspect_svg_file(path)
    if report.get("rootTag") != "svg":
        raise AssertionError(f"Downloaded XML root is not SVG: {report.get('rootTag')}")
    assert_finished_circular_svg(report)
    return {"filename": path.name, "bytes": size, "semantics": report}


def assert_first_linear_svg_download(path: Path) -> dict[str, Any]:
    """Validate the Lambda filename, XML, layout semantics, and safety."""

    if path.name != EXPECTED_FIRST_LINEAR_SVG:
        raise AssertionError(f"Expected {EXPECTED_FIRST_LINEAR_SVG}, downloaded {path.name}")
    size = path.stat().st_size
    if size < MINIMUM_STATIC_SVG_BYTES:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {size} bytes")

    report = inspect_svg_file(path)
    if report.get("rootTag") != "svg":
        raise AssertionError(f"Downloaded XML root is not SVG: {report.get('rootTag')}")
    assert_finished_linear_svg(report)
    return {"filename": path.name, "bytes": size, "semantics": report}


def assert_gui_inputs_svg_download(path: Path) -> dict[str, Any]:
    """Validate the whole-record GFF3 + FASTA static SVG download."""

    if path.name != EXPECTED_GUI_INPUTS_SVG:
        raise AssertionError(f"Expected {EXPECTED_GUI_INPUTS_SVG}, downloaded {path.name}")
    size = path.stat().st_size
    if size < MINIMUM_STATIC_SVG_BYTES:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {size} bytes")

    report = inspect_svg_file(path)
    if report.get("rootTag") != "svg":
        raise AssertionError(f"Downloaded XML root is not SVG: {report.get('rootTag')}")
    assert_first_linear_svg(report)
    return {"filename": path.name, "bytes": size, "semantics": report}


def assert_gui_losatn_tsv_download(
    path: Path,
    reference_path: Path,
) -> dict[str, Any]:
    """Require browser LOSATN output to reproduce the qualified TSV bytes."""

    if path.name != EXPECTED_GUI_LOSATN_TSV:
        raise AssertionError(f"Expected {EXPECTED_GUI_LOSATN_TSV}, downloaded {path.name}")
    contents = path.read_bytes()
    if contents != reference_path.read_bytes():
        raise AssertionError("Browser LOSATN TSV differs from the qualified reference bytes")
    digest = hashlib.sha256(contents).hexdigest()
    if digest != EXPECTED_GUI_LOSATN_TSV_SHA256:
        raise AssertionError(f"Unexpected LOSATN TSV checksum: {digest}")

    rows = [line.split("\t") for line in contents.decode("utf-8").splitlines() if line]
    if len(rows) != 6 or any(len(row) != 12 for row in rows):
        raise AssertionError("Qualified LOSATN TSV must contain six 12-column rows")
    if {(row[0], row[1]) for row in rows} != {
        ("NC_001416.1", "NC_042057.1")
    }:
        raise AssertionError("LOSATN TSV endpoints do not match whole Lambda and DE3")
    if any(max(int(row[6]), int(row[7])) > 48_502 for row in rows):
        raise AssertionError("LOSATN TSV contains a query endpoint outside whole Lambda")
    if any(max(int(row[8]), int(row[9])) > 42_925 for row in rows):
        raise AssertionError("LOSATN TSV contains a subject endpoint outside whole DE3")
    return {
        "filename": path.name,
        "bytes": len(contents),
        "sha256": digest,
        "rows": len(rows),
    }


def assert_gui_losatn_svg_download(path: Path) -> dict[str, Any]:
    """Validate the static two-record LOSATN comparison SVG."""

    if path.name != EXPECTED_GUI_LOSATN_SVG:
        raise AssertionError(f"Expected {EXPECTED_GUI_LOSATN_SVG}, downloaded {path.name}")
    size = path.stat().st_size
    if size < MINIMUM_STATIC_SVG_BYTES:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {size} bytes")

    report = inspect_svg_file(path)
    if report.get("rootTag") != "svg":
        raise AssertionError(f"Downloaded XML root is not SVG: {report.get('rootTag')}")
    assert_gui_losatn_svg(report)
    return {"filename": path.name, "bytes": size, "semantics": report}


def assert_gui_circular_layout_svg_download(path: Path) -> dict[str, Any]:
    """Validate the static four-record mitochondrial grid SVG."""

    if path.name != EXPECTED_GUI_CIRCULAR_LAYOUT_SVG:
        raise AssertionError(
            f"Expected {EXPECTED_GUI_CIRCULAR_LAYOUT_SVG}, downloaded {path.name}"
        )
    size = path.stat().st_size
    if size < MINIMUM_STATIC_SVG_BYTES:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {size} bytes")

    report = inspect_svg_file(path)
    if report.get("rootTag") != "svg":
        raise AssertionError(f"Downloaded XML root is not SVG: {report.get('rootTag')}")
    assert_gui_circular_layout_svg(report)
    return {"filename": path.name, "bytes": size, "semantics": report}
