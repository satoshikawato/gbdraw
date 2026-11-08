"""Numbered labels of Circular scalar axes keep clear of each other and of the data under them.

A Percent GC content axis reserves the room its numbered labels take past its band, as a depth axis does, so the
labels of two neighbouring axes never overlap (OV-2764). A numbered label sits on its track's data and is drawn on
a halo, so it stays readable where the data reach it (OV-2817).
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from itertools import combinations
from pathlib import Path

import pytest

from gbdraw.circular import circular_main
from gbdraw.labels.label_style import HALO_COLOR, HALO_WIDTH_EM

SVG = "{http://www.w3.org/2000/svg}"
XLINK_HREF = "{http://www.w3.org/1999/xlink}href"
NUMBER = re.compile(r"-?[0-9]+(?:\.[0-9]+)?[x%]")
DATA = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "tutorial-data" / "depth-1kb"


def _tutorial_figure(tmp_path: Path, *extra: str) -> ET.Element:
    """The T-CLI-05 / T-GUI-12 tutorial figure (a 52 px depth track at 0-80x right above a 42 px Percent GC content
    track at 10-55 %) with numbered axis labels."""

    circular_main([
        "--gbk", str(DATA / "AP027133.gb"),
        "--depth_track", str(DATA / "AP027133.DRR394922.depth-1kb.tsv"),
        "--depth_window", "1", "--depth_step", "1000", "--depth_min", "0", "--depth_max", "80",
        "--no_depth_log_scale", "--show_depth_axis", "--show_depth_ticks",
        "--depth_large_tick_interval", "20", "--depth_small_tick_interval", "10",
        "--gc", "--skew", "--window", "1000", "--step", "1000",
        "--gc_content_mode", "percent", "--gc_content_min_percent", "10", "--gc_content_max_percent", "55",
        "--show_gc_content_axis", "--show_gc_content_ticks",
        "--gc_content_large_tick_interval", "10", "--gc_content_small_tick_interval", "5",
        "--circular_track_slot", "ticks:ticks@side=outside,tick_label_layout=label_out_tick_in",
        "--circular_track_slot", "features:features@side=overlay,lane_direction=split",
        "--circular_track_slot", "depth_1:depth@track_index=0,side=inside,w=52px",
        "--circular_track_slot", "gc_content:dinucleotide_content@side=inside,w=42px,nt=GC",
        "--circular_track_slot", "gc_skew:dinucleotide_skew@side=inside,w=34px,nt=GC",
        "--labels", "none", "--legend", "none",
        "--depth_axis_labels", "numbered", "--gc_content_axis_labels", "numbered",
        *extra,
        "-o", str(tmp_path / "figure"), "-f", "svg",
    ])
    return ET.parse(tmp_path / "figure.svg").getroot()


def _numbered_labels(root: ET.Element) -> list[ET.Element]:
    labels = [text for text in root.iter(f"{SVG}text") if NUMBER.fullmatch(text.text or "")]
    assert labels
    return labels


@pytest.mark.circular
@pytest.mark.parametrize("gc_font", [(), ("--gc_content_tick_font_size", "11")], ids=["auto", "explicit"])
def test_numbered_labels_of_a_depth_axis_and_the_percent_gc_axis_below_it_do_not_overlap(
    tmp_path: Path, gc_font: tuple[str, ...]
) -> None:
    labels = [
        (text.text or "", float(text.get("y", "nan")), float(text.get("font-size", "nan")))
        for text in _numbered_labels(_tutorial_figure(tmp_path, *gc_font))
    ]
    assert {"0x", "80x", "10%", "55%"} <= {text for text, _, _ in labels}
    # Every label starts right of its axis at the top of the circle, centred on its tick: two labels overlap when
    # their centres are closer than half of each label line.
    for (text_a, y_a, size_a), (text_b, y_b, size_b) in combinations(labels, 2):
        assert abs(y_a - y_b) >= 0.5 * (size_a + size_b) - 1e-6, (text_a, text_b, labels)


@pytest.mark.circular
def test_each_numbered_circular_axis_label_is_drawn_on_a_halo_just_before_it(tmp_path: Path) -> None:
    root = _tutorial_figure(tmp_path)
    for group in root.iter(f"{SVG}g"):
        children = list(group)
        for index, child in enumerate(children):
            if child.tag != f"{SVG}text" or not NUMBER.fullmatch(child.text or ""):
                continue
            halo = children[index - 1]
            assert halo.tag == f"{SVG}use", child.text
            assert (halo.get("href") or halo.get(XLINK_HREF)) == f"#{child.get('id')}"
            assert halo.get("fill") == "none"
            assert halo.get("stroke") == HALO_COLOR
            assert float(halo.get("stroke-width", "nan")) == pytest.approx(
                HALO_WIDTH_EM * float(child.get("font-size", "nan"))
            )
    # Ids are unique within the figure, so each halo outlines its own label.
    ids = [element.get("id") for element in root.iter() if element.get("id")]
    assert len(ids) == len(set(ids))
