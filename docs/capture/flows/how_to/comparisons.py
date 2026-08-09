"""Capture focused protein-comparison GUI How-to journeys."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from playwright.sync_api import BrowserType

from config import (
    GUI_LOSATP_COLLINEAR_SCREENSHOT_NAMES,
    GUI_LOSATP_GROUPS_HOW_TO_SCREENSHOT_NAMES,
)
from flows.bgc_losatp import BgcLosatpResult, capture_bgc_losatp
from flows.hepatoplasmataceae_losatp import (
    HepatoplasmataceaeCollinearResult,
    capture_hepatoplasmataceae_collinear,
)
from flows.how_to.nucleotide_comparisons import (
    capture_gui_circular_rings,
    capture_gui_tlosatx,
    capture_gui_uploaded_comparison,
)
from flows.web_capture import assert_output_paths


__all__ = [
    "capture_gui_circular_rings",
    "capture_gui_losatp_collinear",
    "capture_gui_losatp_groups_how_to",
    "capture_gui_tlosatx",
    "capture_gui_uploaded_comparison",
]


def capture_gui_losatp_groups_how_to(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> BgcLosatpResult:
    """Run H-GUI-07 as the canonical five-record Similarity groups task."""

    assert_output_paths(
        output_paths,
        GUI_LOSATP_GROUPS_HOW_TO_SCREENSHOT_NAMES,
        "H-GUI-07",
    )
    return capture_bgc_losatp(
        browser_type,
        base_url,
        output_paths,
        download_dir,
        scenario_id="H-GUI-07",
        mode="orthogroup",
        output_prefix="bgc_similarity_groups",
        raw_tsv_name="bgc_similarity_groups.tsv",
        screenshot_names={
            "settings": "group-settings.png",
            "result": "group-result.png",
            "popup": "group-popup.png",
        },
        include_plain_result=False,
        download_member_fasta=False,
        separate_strands=False,
    )


def capture_gui_losatp_collinear(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> HepatoplasmataceaeCollinearResult:
    """Run H-GUI-08 with all-record Hepatoplasmataceae evidence."""

    assert_output_paths(
        output_paths,
        GUI_LOSATP_COLLINEAR_SCREENSHOT_NAMES,
        "H-GUI-08",
    )
    return capture_hepatoplasmataceae_collinear(
        browser_type,
        base_url,
        output_paths,
        download_dir,
        output_prefix="hepatoplasmataceae_collinear",
        screenshot_names={
            "settings": "collinear-settings.png",
            "result": "collinear-result.png",
            "detail": "collinear-detail.png",
            "popup": "block-popup.png",
        },
    )
