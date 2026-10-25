"""Labels, ticks, the scale and definitions lay out the box their text is drawn in (OV-165 part 2).

Layout reserves ``drawn_text_box``: the drawn advance, or the ink's right edge
when a last glyph paints past it, and the face's line box. A consumer that
needs room beyond the drawn box adds its own named margin (Circular feature
labels: ``CIRCULAR_LABEL_GAP_EM`` of the width). The line box is the same for
every string of a face, so two labels of one face and size are equally tall
whatever letters they have (OV-1743).
"""

from __future__ import annotations

import math
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group

from gbdraw.api.requests import RecordPresentation
from gbdraw.canvas.linear import LinearCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig, LinearRenderProfile
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.core.text import drawn_text_box, get_font_vertical_metrics, text_ink, text_ink_offsets
from gbdraw.diagrams.linear.assemble import _axis_ruler_extents
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.labels.circular import prepare_label_list
from gbdraw.labels.linear import linear_record_labels
from gbdraw.render.drawers.circular.labels import LabelDrawer as CircularLabelDrawer
from gbdraw.render.groups.linear.definition import DefinitionGroup, LinearDefinitionStyle
from gbdraw.render.groups.linear.length_bar import RULER_LABEL_OFFSET, RULER_TICK_LENGTH, LengthBarGroup
from gbdraw.svg.circular_ticks import resolve_circular_tick_label_geometry

# The must-pass cases (OV-165 PR-2 survey): an inward radial ``multipart enzyme`` lane, HmmtDNA's embedded
# sRNA, and a centre name line that stays one line (rule L never wraps it). They pass with the point-mode reserve and keep passing.
from tests.test_circular_definition_fit import (  # noqa: F401
    test_a_name_line_is_never_wrapped as test_guard_name_line_stays_one_line,
)
from tests.test_circular_label_placement import (  # noqa: F401
    test_hmmtdna_label_override_keeps_current_embedded_flags as test_guard_hmmtdna_embedded_srna,
)
from tests.test_feature_placement import (  # noqa: F401
    test_final_svg_isolated_secondary_lane_keeps_empty_main as test_guard_multipart_enzyme_lanes,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SESSIONS = REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions"
SANS = "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif"
FAMILIES = (SANS, "Times New Roman, serif", "Courier New, monospace")
# Captions with and without ascenders and descenders, with 0 to 5 interior spaces.
CAPTIONS = (
    "CDS",
    "acen",
    "Agy",
    "dnaK",
    "GC skew (-)",
    "type IV pilin protein",
    "a b c d e",
    "i l i l i l",
    "hypothetical protein",
)
SIZES = (8.0, 14.0, 24.0)


def _config(overrides: dict[str, object]) -> GbdrawConfig:
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {"labels.filtering.blacklist_keywords": [], **overrides},
    )
    return GbdrawConfig.from_dict(config_dict)


def _feature(feature_id: str, start: int, end: int, label: str) -> FeatureObject:
    location = [FeatureLocationPart("block", "001", "positive", start, end, True)]
    return FeatureObject(
        feature_id=feature_id,
        location=location,
        is_directional=True,
        color="#54bcf8",
        note="",
        label_text=label,
        coordinates=location,
        type="CDS",
        qualifiers={"product": [label]},
    )


def _features(captions: tuple[str, ...]) -> dict[str, FeatureObject]:
    # Ten-bp features far apart: every caption is an external label.
    return {f"f{index}": _feature(f"f{index}", 1000 * index + 1, 1000 * index + 10, caption)
            for index, caption in enumerate(captions)}


def _line_box(family: str, size: float) -> float:
    metrics = get_font_vertical_metrics(family)
    assert metrics is not None
    return (metrics.ascent + metrics.descent) * size


def _ink_right(text: str, family: str, size: float) -> float:
    ink = text_ink(text, family)
    assert ink is not None
    return ink.right * size


def _circular_candidates(family: str, size: float):
    cache: dict[str, object] = {}
    prepare_label_list(
        _features(CAPTIONS),
        1000 * len(CAPTIONS),
        390,
        0.19,
        CircularRenderProfile(_config({"objects.text.font_family": family, "labels.circular.scope": "outer"})),
        label_font_size=size,
        _candidate_cache=cache,  # type: ignore[arg-type]
    )
    candidates = cache["candidates"]
    assert isinstance(candidates, tuple) and len(candidates) == len(CAPTIONS)
    return candidates


# 1, 2: Circular feature labels.
@pytest.mark.parametrize("family", FAMILIES)
def test_circular_label_box_contains_its_drawn_ink(family: str) -> None:
    misses = [
        (candidate.text, candidate.width_px, _ink_right(candidate.text, family, 14.0))
        for candidate in _circular_candidates(family, 14.0)
        if candidate.width_px < _ink_right(candidate.text, family, 14.0) - 1e-9
    ]
    assert not misses


@pytest.mark.parametrize("family", FAMILIES)
def test_circular_label_height_is_the_face_line_box(family: str) -> None:
    heights = {candidate.text: candidate.height_px for candidate in _circular_candidates(family, 14.0)}
    assert heights == pytest.approx({caption: _line_box(family, 14.0) for caption in CAPTIONS})


@pytest.mark.parametrize("family", FAMILIES)
def test_circular_label_width_is_the_drawn_box_plus_the_label_gap(family: str) -> None:
    from gbdraw.labels.label_style import CIRCULAR_LABEL_GAP_EM

    widths = {candidate.text: candidate.width_px for candidate in _circular_candidates(family, 14.0)}
    assert widths == pytest.approx(
        {caption: drawn_text_box(caption, family, 14.0)[0] + CIRCULAR_LABEL_GAP_EM * 14.0 for caption in CAPTIONS}
    )


# 3: an embedded label's text path (OV-1743).
def _embedded_path_radius(text: str) -> float:
    cfg = _config({})
    family = cfg.objects.text.font_family
    label = {
        "feature_id": f"embedded-{text}",
        "stable_id": f"embedded-{text}",
        "input_order": 0,
        "label_text": text,
        "track_id": 0,
        "strand": "positive",
        "middle": 250.0,
        "start": 200.0,
        "end": 300.0,
        "feature_center_radius_px": 300.0,
        "width_px": drawn_text_box(text, family, 14.0)[0],
        "height_px": _line_box(family, 14.0),
        "is_embedded": True,
        "font_size": 14.0,
        "font_family": family,
    }
    group = CircularLabelDrawer(profile=CircularRenderProfile(cfg)).draw(
        label, Group(id="test"), 1000, float(cfg.canvas.circular.radius), float(cfg.canvas.circular.track_ratio)
    )
    match = re.search(r'd="M [^"]*A ?([-\d.e+]+),', group.tostring())
    assert match is not None
    return float(match.group(1))


def test_embedded_label_path_radius_does_not_depend_on_its_letters() -> None:
    assert _embedded_path_radius("acen") == pytest.approx(_embedded_path_radius("Agy"))


# 4: Circular tick labels.
def test_circular_tick_label_reserves_its_drawn_box() -> None:
    geometry = resolve_circular_tick_label_geometry(
        center_radius_px=100.0,
        total_len=50_000,
        size="large",
        tick=10_000,
        label_text="10 kbp",
        font_size=14.0,
        font_family=SANS,
        track_type="tuckin",
        strandedness=True,
        dpi=96,
        label_side="outside",
        tick_side="inside",
        tick_length_px=10.0,
        tick_width=2.0,
        length_threshold=int(_config({}).labels.length_threshold.circular),
    )
    assert geometry.bbox_width_px == pytest.approx(drawn_text_box("10 kbp", SANS, 14.0)[0])
    assert geometry.bbox_height_px == pytest.approx(_line_box(SANS, 14.0))
    assert geometry.radial_outer_px - geometry.radial_inner_px == pytest.approx(_line_box(SANS, 14.0))


# 5: the Linear ruler and scale bar.
def test_linear_axis_ruler_band_holds_the_line_box_of_its_labels() -> None:
    cfg = _config({})
    canvas_config = SimpleNamespace(length_param="short")
    context = SimpleNamespace(
        axis_ruler_enabled=True, feature_track_layout="above", profile=SimpleNamespace(config=cfg)
    )
    font_size = cfg.objects.scale.ruler_label_font_size.for_length_param("short")
    above, below = _axis_ruler_extents(canvas_config, context)  # type: ignore[arg-type]
    expected = max(
        0.5 * float(cfg.objects.scale.stroke_width),
        float(RULER_TICK_LENGTH),
        float(RULER_LABEL_OFFSET) + _line_box(cfg.objects.text.font_family, font_size),
    )
    assert (above, below) == (0.0, pytest.approx(expected))


def _length_bar(style: str) -> LengthBarGroup:
    cfg = _config({"objects.scale.style": style, "objects.scale.interval": 2000})
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1, longest_genome=10_000, profile=LinearRenderProfile(cfg), legend="none"
    )
    return LengthBarGroup(
        float(canvas_config.alignment_width),
        10_000,
        canvas_config.length_param,
        scale=cfg.objects.scale,
        font_family=cfg.objects.text.font_family,
        dpi=cfg.canvas.dpi,
    )


def test_scale_bar_label_reserves_its_drawn_width_and_the_bar_gap() -> None:
    bar = _length_bar("bar")
    width, height = drawn_text_box(bar.label_text, bar.font_family, bar.font_size, bar.font_weight)
    half_stroke = 0.5 * float(bar.stroke_width)
    assert bar.local_bounds.min_x == pytest.approx(bar.start_x - 10.0 - width - half_stroke)
    assert bar.local_bounds.max_y == pytest.approx(max(half_stroke, 0.5 * height))


def test_ruler_labels_reserve_their_line_box_below_the_ruler() -> None:
    ruler = _length_bar("ruler")
    height = _line_box(ruler.font_family, ruler.ruler_label_font_size)
    assert ruler.local_bounds.max_y == pytest.approx(
        max(RULER_TICK_LENGTH + 0.5 * float(ruler.stroke_width), RULER_LABEL_OFFSET + height)
    )


# 6: the Linear definition column.
def _record() -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * 2500), id="rec1", name="rec1", description="Escherichia coli strain K-12")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["organism"] = "Escherichia coli"
    record.features.append(SeqFeature(FeatureLocation(0, 300, strand=1), type="CDS", qualifiers={"product": ["dnaK"]}))
    return record


def test_linear_definition_lines_are_a_line_box_tall_and_interval_apart() -> None:
    cfg = _config({})
    record = _record()
    group = DefinitionGroup(
        record,
        "short",
        style=LinearDefinitionStyle.of(cfg),
        record_presentation=RecordPresentation(label="<i>Escherichia coli</i> K-12", subtitle="plasmid pK"),
    )
    lines = group.definition_lines
    assert len(lines) >= 2
    heights = [line.height for line in lines]
    assert heights == pytest.approx(
        [
            _line_box(group._style_for_line(line.kind).font_family, group._style_for_line(line.kind).font_size)
            for line in lines
        ]
    )
    assert group.definition_bounding_box_height == pytest.approx(sum(heights) + group.interval * (len(lines) - 1))


# 7: Linear feature labels.
def test_linear_label_height_is_the_face_line_box() -> None:
    cfg = _config({"labels.linear.scope": "all"})
    record_labels = linear_record_labels(
        _features(CAPTIONS), 1000 * len(CAPTIONS), 2000, 20, False, "middle", None, LinearRenderProfile(cfg),
        label_font_size=14.0,
    )
    heights = {label["label_text"]: label["height_px"] for label in record_labels.labels}
    assert heights == pytest.approx({caption: _line_box(cfg.objects.text.font_family, 14.0) for caption in CAPTIONS})
    assert record_labels.max_label_height == pytest.approx(_line_box(cfg.objects.text.font_family, 14.0))


# 8: corpus property over every face and size.
def test_every_label_and_tick_layout_width_holds_the_drawn_ink() -> None:
    misses = []
    for family in FAMILIES:
        for size in SIZES:
            for candidate in _circular_candidates(family, size):
                if candidate.width_px < _ink_right(candidate.text, family, size) - 1e-9:
                    misses.append(("circular label", candidate.text, family, size))
            cfg = _config({"objects.text.font_family": family, "labels.linear.scope": "all"})
            linear = linear_record_labels(
                _features(CAPTIONS), 1000 * len(CAPTIONS), 2000, 20, False, "middle", None,
                LinearRenderProfile(cfg), label_font_size=size,
            )
            for label in linear.labels:
                if label["width_px"] < _ink_right(label["label_text"], family, size) - 1e-9:
                    misses.append(("linear label", label["label_text"], family, size))
            for text in ("0", "10 kbp", "1.5 Mbp", "250 bp"):
                geometry = resolve_circular_tick_label_geometry(
                    center_radius_px=100.0, total_len=50_000, size="large", tick=10_000, label_text=text,
                    font_size=size, font_family=family, track_type="tuckin", strandedness=True, dpi=96,
                    label_side="outside", tick_side="inside", tick_length_px=10.0, tick_width=2.0,
                    length_threshold=int(cfg.labels.length_threshold.circular),
                )
                if geometry.bbox_width_px < _ink_right(text, family, size) - 1e-9:
                    misses.append(("tick label", text, family, size))
    assert not misses, f"{len(misses)} layout widths are narrower than their ink, e.g. {misses[:5]}"


# 11: the drawn labels of Gallery-quality figures neither overlap nor go missing.
_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_SVG = "{http://www.w3.org/2000/svg}"
Affine = tuple[float, float, float, float, float, float]


def _compose(a: Affine, b: Affine) -> Affine:
    return (
        a[0] * b[0] + a[2] * b[1],
        a[1] * b[0] + a[3] * b[1],
        a[0] * b[2] + a[2] * b[3],
        a[1] * b[2] + a[3] * b[3],
        a[0] * b[4] + a[2] * b[5] + a[4],
        a[1] * b[4] + a[3] * b[5] + a[5],
    )


def _transform(text: str | None) -> Affine:
    matrix: Affine = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    for name, args in re.findall(r"(\w+)\s*\(([^)]*)\)", text or ""):
        values = [float(value) for value in re.findall(_NUM, args)]
        if name == "translate":
            step: Affine = (1.0, 0.0, 0.0, 1.0, values[0], values[1] if len(values) > 1 else 0.0)
        elif name == "rotate":
            cos, sin = math.cos(math.radians(values[0])), math.sin(math.radians(values[0]))
            cx, cy = (values[1], values[2]) if len(values) == 3 else (0.0, 0.0)
            step = (cos, sin, -sin, cos, cx - cos * cx + sin * cy, cy - sin * cx - cos * cy)
        else:
            raise AssertionError(f"unexpected transform {name}")
        matrix = _compose(matrix, step)
    return matrix


def _drawn_label_quads(svg: Path) -> list[list[tuple[float, float]]]:
    """The ink of every drawn external feature label (not on a text path), as a quadrilateral on the canvas."""
    quads: list[list[tuple[float, float]]] = []

    def walk(element: ET.Element, matrix: Affine) -> None:
        if element.tag == f"{_SVG}defs":
            return
        matrix = _compose(matrix, _transform(element.get("transform")))
        if (
            element.tag == f"{_SVG}text"
            and element.get("data-label-feature-id")
            and not any(child.tag == f"{_SVG}textPath" for child in element.iter())
        ):
            offsets = text_ink_offsets(
                "".join(element.itertext()),
                element.get("font-family", ""),
                float(re.findall(_NUM, element.get("font-size", "12"))[0]),
                dominant_baseline=element.get("dominant-baseline", "auto"),
                text_anchor=element.get("text-anchor", "start"),
            )
            assert offsets is not None
            left, right, top, bottom = offsets
            x, y = float(element.get("x", 0.0)), float(element.get("y", 0.0))
            quads.append([
                (matrix[0] * px + matrix[2] * py + matrix[4], matrix[1] * px + matrix[3] * py + matrix[5])
                for px, py in ((x + left, y + top), (x + right, y + top), (x + right, y + bottom), (x + left, y + bottom))
            ])
        for child in element:
            walk(child, matrix)

    walk(ET.parse(svg).getroot(), (1.0, 0.0, 0.0, 1.0, 0.0, 0.0))
    return quads


def _overlap_depth(p: list[tuple[float, float]], q: list[tuple[float, float]]) -> float:
    """How deep two convex quadrilaterals overlap along their separating axes (<= 0: apart)."""
    depth = math.inf
    for quad in (p, q):
        for (x1, y1), (x2, y2) in zip(quad, quad[1:] + quad[:1]):
            norm = math.hypot(x2 - x1, y2 - y1) or 1.0
            axis = (-(y2 - y1) / norm, (x2 - x1) / norm)
            pp = [axis[0] * x + axis[1] * y for x, y in p]
            qq = [axis[0] * x + axis[1] * y for x, y in q]
            depth = min(depth, min(max(pp), max(qq)) - max(min(pp), min(qq)))
    return depth


GALLERY_CASES = {
    "tobacco": (["--session", str(SESSIONS / "tobacco-chloroplast.gbdraw-session.json")], 118),
    "hmmtdna-serif": (
        [
            "--gbk", str(REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"), "--separate_strands",
            "--track_type", "middle", "--labels", "out", "--species", "<i>Homo sapiens</i>",
            "--font_family", "Liberation Serif",
        ],
        32,
    ),
    "hmmtdna-basic": (["--session", str(SESSIONS / "HmmtDNA_basic_circular.gbdraw-session.json")], 32),
}


@pytest.mark.circular
@pytest.mark.parametrize("case", sorted(GALLERY_CASES))
def test_gallery_labels_keep_their_count_and_do_not_overlap(case: str, tmp_path: Path) -> None:
    args, expected_count = GALLERY_CASES[case]
    output = tmp_path / case
    subprocess.run(
        [sys.executable, "-m", "gbdraw.cli", "circular", *args, "-f", "svg", "--overwrite", "-o", str(output)],
        cwd=REPO_ROOT, check=True, capture_output=True, text=True,
    )
    quads = _drawn_label_quads(output.with_suffix(".svg"))
    overlaps = [
        (i, j) for i in range(len(quads)) for j in range(i + 1, len(quads)) if _overlap_depth(quads[i], quads[j]) > 1e-3
    ]
    assert (len(quads), overlaps) == (expected_count, [])
