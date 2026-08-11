"""Capture the five-record LOSATP Similarity groups tutorial."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType

from assertions.svg_semantics import (
    assert_gui_bgc_gallery_definition_presentation,
    assert_gui_bgc_similarity_groups_svg,
)
from config import GUI_LOSATP_GROUPS_SCREENSHOT_NAMES
from flows.bgc_losatp import BgcLosatpResult, capture_bgc_losatp
from flows.web_capture import assert_output_paths


def capture_gui_losatp_groups(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> BgcLosatpResult:
    """Run T-GUI-04 from five complete BGC records through the visible UI."""

    assert_output_paths(
        output_paths,
        GUI_LOSATP_GROUPS_SCREENSHOT_NAMES,
        "T-GUI-04",
    )
    def assert_final(report: dict[str, Any]) -> None:
        assert_gui_bgc_similarity_groups_svg(report)
        assert_gui_bgc_gallery_definition_presentation(report)

    return capture_bgc_losatp(
        browser_type,
        base_url,
        output_paths,
        download_dir,
        scenario_id="T-GUI-04",
        mode="orthogroup",
        output_prefix="bgc_losatp_groups",
        raw_tsv_name="bgc_losatp_groups.tsv",
        screenshot_names={
            "input": "01-input-ready.png",
            "plain": "02-first-diagram.png",
            "settings": "03-losatp-settings.png",
            "align": "04-align-og1.png",
            "result": "05-comparison-result.png",
            "popup": "06-match-popup.png",
        },
        include_plain_result=True,
        download_member_fasta=False,
        separate_strands=False,
        align_orthogroup_id="og_1",
        match_gallery_definitions=True,
        assert_final=assert_final,
    )
