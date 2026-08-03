"""Visible-UI setup shared by finished human mitochondrial GUI workflows."""

from __future__ import annotations

from typing import Any

from Bio import SeqIO
from playwright.sync_api import Page, expect

from assertions.svg_semantics import (
    assert_finished_circular_svg,
    inspect_first_circular_svg,
)
from config import (
    FIRST_CIRCULAR_FIXTURE_PATH,
    FIRST_CIRCULAR_FIXTURE_SHA256,
    FIRST_CIRCULAR_FIXTURE_SIZE,
)
from flows.web_capture import (
    assert_fixture_identity,
    generate_and_inspect,
)


HUMAN_RECORD_ID = "NC_012920.1"
HUMAN_RECORD_LENGTH = 16_569
HUMAN_SPECIES_MARKUP = "<i>Homo sapiens</i>"


def assert_human_mitochondrion_fixture() -> dict[str, Any]:
    """Require the unchanged complete circular RefSeq fixture."""

    assert_fixture_identity(
        FIRST_CIRCULAR_FIXTURE_PATH,
        expected_size=FIRST_CIRCULAR_FIXTURE_SIZE,
        expected_sha256=FIRST_CIRCULAR_FIXTURE_SHA256,
    )
    records = list(SeqIO.parse(FIRST_CIRCULAR_FIXTURE_PATH, "genbank"))
    if len(records) != 1:
        raise AssertionError("HmmtDNA.gbk must contain exactly one complete record")
    record = records[0]
    if (
        record.id != HUMAN_RECORD_ID
        or len(record) != HUMAN_RECORD_LENGTH
        or str(record.annotations.get("topology", "")).lower() != "circular"
        or "complete" not in record.description.lower()
    ):
        raise AssertionError(
            "Finished-diagram workflows require complete circular NC_012920.1"
        )
    return {
        "record_id": record.id,
        "length": len(record),
        "topology": record.annotations.get("topology"),
        "description": record.description,
    }


def load_raw_human_circular(page: Page, *, output_prefix: str) -> None:
    """Load HmmtDNA through the public Circular GenBank controls."""

    circular = page.get_by_role("button", name="Circular", exact=True)
    circular.click()
    expect(circular).to_have_attribute("aria-pressed", "true")
    genbank = page.get_by_role("radio", name="GenBank", exact=True)
    genbank.check()
    expect(genbank).to_be_checked()
    page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
        FIRST_CIRCULAR_FIXTURE_PATH
    )
    expect(
        page.get_by_role(
            "group", name="GenBank/DDBJ File selection", exact=True
        )
    ).to_contain_text(FIRST_CIRCULAR_FIXTURE_PATH.name)
    apply_finished_human_settings(page, output_prefix=output_prefix)


def apply_finished_human_settings(page: Page, *, output_prefix: str) -> None:
    """Apply the accepted T-GUI-01 finished state through visible controls."""

    prefix = page.get_by_label("Output Prefix", exact=True)
    prefix.fill(output_prefix)
    expect(prefix).to_have_value(output_prefix)

    species = page.get_by_label("Species", exact=True)
    species.fill(HUMAN_SPECIES_MARKUP)
    expect(species).to_have_value(HUMAN_SPECIES_MARKUP)

    track_preset = page.get_by_label("Track Preset", exact=True)
    track_preset.select_option("middle")
    expect(track_preset).to_have_value("middle")
    separate_strands = page.get_by_label("Separate Strands", exact=True)
    separate_strands.check()
    expect(separate_strands).to_be_checked()
    hide_gc_content = page.get_by_label("Hide GC Content", exact=True)
    hide_gc_content.uncheck()
    expect(hide_gc_content).not_to_be_checked()
    hide_gc_skew = page.get_by_label("Hide GC Skew", exact=True)
    hide_gc_skew.uncheck()
    expect(hide_gc_skew).not_to_be_checked()

    title_and_legend = page.get_by_label("Title & Legend", exact=True)
    title_and_legend.click()
    legend_position = page.get_by_label("Legend Position", exact=True)
    legend_position.select_option("right")
    expect(legend_position).to_have_value("right")
    title_and_legend.click()

    labels = page.get_by_label("Labels", exact=True)
    labels.click()
    label_mode = page.get_by_label("Label Mode", exact=True)
    label_mode.select_option("out")
    expect(label_mode).to_have_value("out")
    labels.click()


def generate_finished_human_diagram(page: Page) -> dict[str, Any]:
    """Generate and validate the accepted finished Circular diagram."""

    return generate_and_inspect(
        page,
        inspect_first_circular_svg,
        assert_finished_circular_svg,
    )


def fit_finished_human_preview(page: Page) -> None:
    """Reach the documented 60 percent preview zoom with public controls."""

    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    for _ in range(10):
        if "60%" in reset_zoom.inner_text():
            return
        zoom_out.click()
    raise AssertionError("Could not reach the documented 60% Circular preview zoom")


def pan_finished_human_preview_left(page: Page, *, distance_ratio: float = 0.25) -> None:
    """Pan the public preview so the right-side legend remains in frame."""

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    box = result_region.bounding_box()
    if box is None:
        raise AssertionError("Could not resolve the Result Preview bounds")
    y = box["y"] + (box["height"] * 0.72)
    page.mouse.move(box["x"] + (box["width"] * 0.86), y)
    page.mouse.down()
    page.mouse.move(
        box["x"] + (box["width"] * (0.86 - distance_ratio)),
        y,
        steps=10,
    )
    page.mouse.up()
    page.wait_for_timeout(250)


__all__ = [
    "HUMAN_RECORD_ID",
    "HUMAN_RECORD_LENGTH",
    "HUMAN_SPECIES_MARKUP",
    "apply_finished_human_settings",
    "assert_human_mitochondrion_fixture",
    "fit_finished_human_preview",
    "generate_finished_human_diagram",
    "load_raw_human_circular",
    "pan_finished_human_preview_left",
]
