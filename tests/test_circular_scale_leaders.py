"""OV-3301: feature labels, their leader lines and the scale's numbers.

With a scale standing between the feature track and the outer labels (a ticks slot outside the features, its labels
outside the ticks), the outer labels and their leader elbows stood on the numbers, and the leaders ran through them
("8 kbp", "9 kbp" of NC_002333.2 in docs H-CLI-03). The figures here are drawn by the command line and measured
from the SVG alone: each number's letters as an annular sector (its label path's radius and the letters' measured
ink), each label's letters as a box, each leader as its two lines, each backing as the sector its path draws.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pytest

from gbdraw import circular as circular_cli_module
from gbdraw.core.text import text_ink, text_ink_offsets

TUTORIAL_DATA = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "tutorial-data"
ZEBRAFISH = TUTORIAL_DATA / "metazoan-mitochondria-four" / "NC_002333.2.gb"
HUMAN = TUTORIAL_DATA / "human-mitochondrion" / "HmmtDNA.gbk"
PRIORITY = TUTORIAL_DATA / "shared" / "cds_gene_qualifier_priority.tsv"
SCALE_FIRST = ["--circular_track_order", "ticks,features,gc_content,gc_skew", "--circular_track_axis_index", "1"]
EXPLICIT_OUT = [
    "--circular_track_slot", "ticks:ticks@tick_label_layout=label_out_tick_in",
    "--circular_track_slot", "features:features",
    "--circular_track_slot", "gc_content:dinucleotide_content",
    "--circular_track_axis_index", "1",
]


def _draw(tmp_path: Path, gbk: Path, *args: str) -> str:
    out = tmp_path / "figure"
    circular_cli_module.circular_main(
        ["--gbk", str(gbk), "--qualifier_priority", str(PRIORITY), *args, "-o", str(out), "-f", "svg"]
    )
    return (tmp_path / "figure.svg").read_text(encoding="utf-8")


def _group(svg: str, marker: str) -> tuple[str, str]:
    start = svg.index(marker)
    tag_start = svg.rindex("<g ", 0, start)
    tag_end = svg.index(">", start) + 1
    depth, at = 1, tag_end
    while depth:
        opening, closing = svg.find("<g", at), svg.find("</g>", at)
        if opening != -1 and opening < closing:
            depth, at = depth + 1, opening + 2
        else:
            depth, at = depth - 1, closing + 4
    return svg[tag_start:tag_end], svg[tag_end:at]


def _wrapped(angle: float) -> float:
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


class _Sector:
    """An annular sector around the circle centre: ``inner`` to ``outer`` radius, ``middle`` +/- ``half`` radians."""

    def __init__(self, inner: float, outer: float, low: float, high: float) -> None:
        self.inner, self.outer = inner, outer
        self.middle = math.atan2(math.sin(low) + math.sin(high), math.cos(low) + math.cos(high))
        self.half = 0.5 * abs(_wrapped(high - low))

    def holds(self, point: tuple[float, float], pad: float = 0.0) -> bool:
        radius = math.hypot(*point)
        if not self.inner - pad <= radius <= self.outer + pad:
            return False
        return abs(_wrapped(math.atan2(point[1], point[0]) - self.middle)) * radius <= self.half * radius + pad

    def covers(self, other: _Sector) -> bool:
        return (
            self.inner <= other.inner + 1e-6
            and self.outer >= other.outer - 1e-6
            and abs(_wrapped(other.middle - self.middle)) + other.half <= self.half + 1e-6
        )


def _numbers(ticks_body: str) -> list[tuple[str, _Sector]]:
    """Each scale number's letters: its path's radius and ends, the text centred on it, the ink measured."""

    numbers = []
    for match in re.finditer(
        r'<path d="M ([-\d.]+),([-\d.]+) A ([\d.]+),[\d.]+ 0 0 \d ([-\d.]+),([-\d.]+)"[^>]*'
        r'id="circular_tick_label_path_[^"]+"[^>]*/><text[^>]*><textPath dominant-baseline="([^"]+)" '
        r'font-family="([^"]+)" font-size="([\d.]+)"[^>]*>([^<]+)<',
        ticks_body,
    ):
        x1, y1, radius, x2, y2 = map(float, match.group(1, 2, 3, 4, 5))
        baseline, family, size, text = match.group(6), match.group(7).replace("&apos;", "'"), float(match.group(8)), match.group(9)
        ink = text_ink(text, family)
        offsets = text_ink_offsets(text, family, size, dominant_baseline=baseline, text_anchor="middle")
        assert ink is not None and offsets is not None
        _left, _right, top, bottom = offsets
        # The letters stand outward of the path on both halves: above an upper path, below a reversed lower one.
        inner, outer = (radius - bottom, radius - top) if baseline == "text-after-edge" else (radius + top, radius + bottom)
        middle = math.atan2(math.sin(math.atan2(y1, x1)) + math.sin(math.atan2(y2, x2)), math.cos(math.atan2(y1, x1)) + math.cos(math.atan2(y2, x2)))
        half = 0.5 * (ink.right - ink.left) * size / radius
        centre = middle + (0.5 * (ink.left + ink.right) - 0.5 * ink.advance) * size / radius * (1 if baseline == "text-after-edge" else -1)
        numbers.append((text, _Sector(inner, outer, centre - half, centre + half)))
    return numbers


def _backings(ticks_body: str) -> list[_Sector]:
    sectors = []
    for match in re.finditer(
        r'<path d="M ([-\d.]+),([-\d.]+) A ([\d.]+),[\d.]+ 0 \d 1 ([-\d.]+),([-\d.]+) L [-\d.]+,[-\d.]+ '
        r'A ([\d.]+),[\d.]+ 0 \d 0 [-\d.]+,[-\d.]+ Z" fill="#ffffff" fill-opacity="[\d.]+"',
        ticks_body,
    ):
        x1, y1, outer, x2, y2, inner = map(float, match.groups())
        sectors.append(_Sector(inner, outer, math.atan2(y1, x1), math.atan2(y2, x2)))
    return sectors


def _leaders(svg: str) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    _tag, body = _group(svg, 'id="label_leaders"')
    lines = []
    for match in re.finditer(r"<line ([^>]*)/>", body):
        attrs = dict(re.findall(r'([\w-]+)="([^"]*)"', match.group(1)))
        lines.append(((float(attrs["x1"]), float(attrs["y1"])), (float(attrs["x2"]), float(attrs["y2"]))))
    return lines


def _label_boxes(svg: str) -> list[tuple[str, tuple[float, float, float, float]]]:
    _tag, body = _group(svg, 'id="label_text"')
    boxes = []
    for match in re.finditer(r"<text ([^>]*)>([^<]*)</text>", body):
        attrs = dict(re.findall(r'([\w-]+)="([^"]*)"', match.group(1)))
        offsets = text_ink_offsets(
            match.group(2), attrs["font-family"].replace("&apos;", "'"), float(attrs["font-size"]),
            dominant_baseline=attrs.get("dominant-baseline", "auto"), text_anchor=attrs.get("text-anchor", "start"),
        )
        assert offsets is not None
        x, y = float(attrs["x"]), float(attrs["y"])
        boxes.append((match.group(2), (x + offsets[0], y + offsets[2], x + offsets[1], y + offsets[3])))
    return boxes


def _meets(sector: _Sector, start: tuple[float, float], end: tuple[float, float], pad: float) -> bool:
    steps = max(2, int(math.hypot(end[0] - start[0], end[1] - start[1]) * 4))
    return any(
        sector.holds((start[0] + (end[0] - start[0]) * i / steps, start[1] + (end[1] - start[1]) * i / steps), pad)
        for i in range(steps + 1)
    )


def _scale_and_labels(svg: str) -> dict[str, object]:
    """What the leaders and labels do to the numbers: the numbers a leader crosses or touches without a backing, the labels
    whose letters meet a number, the leader elbows in the numbers' band, and the clearances kept."""

    ticks_tag, ticks_body = _group(svg, 'data-gbdraw-slot-renderer="ticks"')
    leaders_tag, _ = _group(svg, 'id="label_leaders"')
    assert re.findall(r"translate\([^)]*\)", ticks_tag) == re.findall(r"translate\([^)]*\)", leaders_tag)
    numbers = _numbers(ticks_body)
    assert numbers
    backings = _backings(ticks_body)
    leaders = _leaders(svg)
    # A leader (1 px wide) whose stroke reaches the letters or comes within a backing's padding (0.18 em) of them.
    pad = 0.5 + 0.18 * float(re.findall(r'font-size="([\d.]+)"', ticks_body)[0])
    crossed = {text for text, sector in numbers for start, end in leaders if _meets(sector, start, end, pad)}
    unbacked = sorted(
        text for text, sector in numbers if text in crossed and not any(backing.covers(sector) for backing in backings)
    )
    on_numbers = sorted(
        {
            (label, text)
            for label, (x0, y0, x1, y1) in _label_boxes(svg)
            for text, sector in numbers
            if any(sector.holds((x0 + (x1 - x0) * i / 12, y0 + (y1 - y0) * j / 4)) for i in range(13) for j in range(5))
        }
    )
    band = (min(sector.inner for _t, sector in numbers), max(sector.outer for _t, sector in numbers))
    elbows = [start for start, _end in leaders[::2]]
    return {
        "unbacked": unbacked,
        "crossed": sorted(crossed),
        "on_numbers": on_numbers,
        "elbows_in_band": sum(band[0] <= math.hypot(*point) <= band[1] for point in elbows),
        "band": band,
    }


@pytest.mark.circular
@pytest.mark.parametrize(
    ("gbk", "args"),
    [
        (ZEBRAFISH, [*SCALE_FIRST, "--labels", "out", "--label_font_size", "10"]),
        (HUMAN, [*SCALE_FIRST, "--labels", "out"]),
        (ZEBRAFISH, [*EXPLICIT_OUT, "--labels", "out", "--scale_interval", "2000"]),
        (HUMAN, [*EXPLICIT_OUT, "--labels", "out", "--label_font_size", "10", "--scale_interval", "500"]),
    ],
    ids=["zebrafish-scale-first", "human-scale-first", "zebrafish-explicit-out-2kbp", "human-explicit-out-500bp"],
)
def test_outer_labels_and_their_leaders_keep_clear_of_a_scale_outside_the_features(tmp_path, gbk, args) -> None:
    """OV-3301: outer labels and their leader elbows stand past the numbers of a scale outside the features, and a
    number a leader runs through is drawn in front of it on a backing, as a region label is (CIR-RGN-FRONT)."""

    measured = _scale_and_labels(_draw(tmp_path, gbk, *args))

    assert measured["on_numbers"] == []
    assert measured["elbows_in_band"] == 0
    assert measured["unbacked"] == []


@pytest.mark.circular
def test_inner_leaders_through_the_default_scale_pass_behind_backed_numbers(tmp_path) -> None:
    """The same rule inside the features: with ``--labels both`` the inner labels' leaders cross the default scale
    ("5 kbp" and "8 kbp" of NC_002333.2), and the numbers they cross are backed."""

    measured = _scale_and_labels(_draw(tmp_path, ZEBRAFISH, "--labels", "both"))

    assert measured["crossed"]
    assert measured["unbacked"] == []
    assert measured["on_numbers"] == []
