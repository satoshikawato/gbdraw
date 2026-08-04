"""Capture the LOSATP Collinear project tutorial."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from playwright.sync_api import BrowserType

from config import GUI_LOSATP_TUTORIAL_COLLINEAR_SCREENSHOT_NAMES
from flows.hepatoplasmataceae_losatp import (
    HepatoplasmataceaeCollinearResult,
    capture_hepatoplasmataceae_collinear,
)
from flows.web_capture import assert_output_paths


def capture_gui_losatp_tutorial_collinear(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> HepatoplasmataceaeCollinearResult:
    """Find and inspect conserved-order blocks across five genomes."""

    assert_output_paths(
        output_paths,
        GUI_LOSATP_TUTORIAL_COLLINEAR_SCREENSHOT_NAMES,
        "T-GUI-08",
    )
    return capture_hepatoplasmataceae_collinear(
        browser_type,
        base_url,
        output_paths,
        download_dir,
        output_prefix="losatp_collinear",
        screenshot_names={
            "input": "01-input-ready.png",
            "plain": "02-first-diagram.png",
            "settings": "03-collinear-settings.png",
            "result": "04-collinear-result.png",
            "popup": "05-block-popup.png",
        },
        starting_project="gallery-collinear",
    )


__all__ = ["capture_gui_losatp_tutorial_collinear"]
