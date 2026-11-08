"""Scalar axes (GC content in percent, depth) carry one range label by default (D-16); numbered tick labels never
overlap (OV-2763); the Percent band is light (OV-2759)."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api.assembly_inputs import CircularAssemblyInputs, LinearAssemblyInputs
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_circular_diagram_from_record, assemble_linear_diagram_from_records
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml

SVG = "{http://www.w3.org/2000/svg}"

# The T-CLI-05 tutorial's GC content axis: 10-55 % with a large tick every 10 %.
DOCS_PERCENT_AXIS = {
    "objects.gc_content.mode": "percent",
    "objects.gc_content.min_percent": 10,
    "objects.gc_content.max_percent": 55,
    "objects.gc_content.large_tick_interval": 10,
}
NUMBERED = {"objects.gc_content.axis_labels": "numbered", "objects.depth.axis_labels": "numbered"}


def _make_record(length: int = 4000) -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * (length // 4)), id="rec1", name="rec1")
    record.annotations["molecule_type"] = "DNA"
    record.features = [SeqFeature(FeatureLocation(10, 50, strand=1), type="CDS", qualifiers={"product": ["p"]})]
    return record


def _depth_table(length: int = 4000) -> pd.DataFrame:
    positions = list(range(1, length + 1, 100))
    return pd.DataFrame({"reference_name": ["rec1"] * len(positions), "position": positions,
                         "depth": [30.0] * len(positions)})


def _axis_labels(svg: str, unit: str) -> list[tuple[str, float, float]]:
    """The tick labels of the one scalar axis whose labels end with ``unit``: text, y of the label centre,
    font size. A slot-placed track renames its axis group, so the labels are found by their axis lines."""

    axis_label = re.compile(rf"-?[0-9]+(?:\.[0-9]+)?{re.escape(unit)}")
    labels = [
        (text.text or "", float(text.get("y", "nan")), float(text.get("font-size", "nan")))
        for group in ET.fromstring(svg).iter(f"{SVG}g")
        if group.find(f"{SVG}line") is not None
        for text in group.findall(f"{SVG}text")
        if axis_label.fullmatch(text.text or "")
    ]
    assert labels
    return labels


def _assert_no_overlap(labels: list[tuple[str, float, float]]) -> None:
    centres = sorted(y for _, y, _ in labels)
    (font_size,) = {size for _, _, size in labels}
    for lower, upper in zip(centres, centres[1:]):
        assert upper - lower >= font_size - 1e-6, labels


@pytest.mark.parametrize(
    ("positions", "kept"),
    [
        # The upper bound is kept first, then the lower bound, then interior labels from the bottom up.
        ((0.0, 3.0, 10.0, 12.0, 20.0, 21.5), (0, 2, 5)),
        # The lower bound yields to the upper bound when they are closer than one label line.
        ((0.0, 2.0), (1,)),
        # Uneven spacing (a log-scaled axis): only the labels a line apart from every kept one stay.
        ((0.0, 9.0, 13.0, 15.5, 17.0, 18.0), (0, 1, 5)),
        ((0.0, 5.0, 10.0), (0, 1, 2)),
        ((4.0,), (0,)),
        ((), ()),
    ],
)
def test_the_label_rule_keeps_the_upper_bound_then_labels_a_line_apart_from_the_bottom_up(
    positions: tuple[float, ...], kept: tuple[int, ...]
) -> None:
    from gbdraw.layout.scalar_axis import scalar_axis_kept_labels

    assert scalar_axis_kept_labels(positions, 5.0) == kept


@pytest.mark.circular
def test_circular_percent_axis_of_the_docs_figure_draws_its_bound_and_no_overlapping_labels() -> None:
    svg = assemble_circular_diagram_from_record(
        _make_record(),
        CircularAssemblyInputs(
            legend="none",
            circular_track_slots=["gc_content:dinucleotide_content@side=inside,w=42px,nt=GC"],
            cfg=apply_config_overrides(None, {"canvas.show_skew": False, **DOCS_PERCENT_AXIS, **NUMBERED}),
            window=100,
            step=100,
        ),
    ).tostring()

    labels = _axis_labels(svg, "%")
    texts = [text for text, _, _ in labels]
    assert "55%" in texts
    assert "10%" in texts
    _assert_no_overlap(labels)


@pytest.mark.linear
def test_linear_percent_axis_of_the_docs_figure_draws_its_bound_and_no_overlapping_labels() -> None:
    svg = assemble_linear_diagram_from_records(
        [_make_record()],
        LinearAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {"canvas.show_gc": True, "canvas.show_skew": False,
                                              **DOCS_PERCENT_AXIS, **NUMBERED}),
            window=100,
            step=100,
        ),
    ).tostring()

    labels = _axis_labels(svg, "%")
    texts = [text for text, _, _ in labels]
    assert "55%" in texts
    _assert_no_overlap(labels)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_depth_maximum_off_the_interval_grid_keeps_its_label_and_drops_the_next(mode: str) -> None:
    overrides = {
        "canvas.show_gc": False,
        "canvas.show_skew": False,
        "objects.depth.normalize": False,
        "objects.depth.min_depth": 0,
        "objects.depth.max_depth": 62,
        "objects.depth.large_tick_interval": 20,
        **NUMBERED,
    }
    if mode == "circular":
        svg = assemble_circular_diagram_from_record(
            _make_record(),
            CircularAssemblyInputs(legend="none", depth_table=_depth_table(), window=100, step=100,
                                   cfg=apply_config_overrides(None, overrides)),
        ).tostring()
    else:
        svg = assemble_linear_diagram_from_records(
            [_make_record()],
            LinearAssemblyInputs(legend="none", depth_table=_depth_table(), window=100, step=100,
                                 cfg=apply_config_overrides(None, overrides)),
        ).tostring()

    labels = _axis_labels(svg, "x")
    texts = [text for text, _, _ in labels]
    assert "62x" in texts
    assert "60x" not in texts
    assert "0x" in texts
    _assert_no_overlap(labels)


def test_the_percent_band_is_light_gray_by_default() -> None:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    assert cfg.objects.gc_content.percent_background_color == "#e5e5e5"

    svg = assemble_circular_diagram_from_record(
        _make_record(),
        CircularAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {"canvas.show_skew": False, "objects.gc_content.mode": "percent"}),
            window=100,
            step=100,
        ),
    ).tostring()
    background = next(
        element for element in ET.fromstring(svg).iter(f"{SVG}circle")
        if element.get("id") == "gc_content_percent_background"
    )
    assert background.get("stroke") == "#e5e5e5"


def _range_labels(root: ET.Element) -> list[ET.Element]:
    """Every text or text-along-a-path element whose content is a range label ``[...]``."""

    return [
        element
        for element in root.iter()
        if element.tag in (f"{SVG}text", f"{SVG}textPath") and (element.text or "").startswith("[")
    ]


def _depth_overrides(**extra: object) -> dict[str, object]:
    return {
        "canvas.show_gc": False,
        "canvas.show_skew": False,
        "objects.depth.min_depth": 0,
        "objects.depth.max_depth": 62,
        "objects.depth.large_tick_interval": 20,
        **extra,
    }


def test_the_range_label_prints_bounds_like_the_tick_labels_with_an_en_dash() -> None:
    from gbdraw.layout.scalar_axis import scalar_axis_range_label

    assert scalar_axis_range_label(10, 55, "%") == "[10\u201355%]"
    assert scalar_axis_range_label(0, 74.5, "x") == "[0\u201374.5x]"
    assert scalar_axis_range_label(1, 1000, "x", log_scale=True) == "[1\u20131000x, log]"


@pytest.mark.circular
def test_circular_axes_draw_every_tick_and_one_range_label_along_the_circle_by_default() -> None:
    svg = assemble_circular_diagram_from_record(
        _make_record(),
        CircularAssemblyInputs(
            legend="none",
            depth_table=_depth_table(),
            circular_track_slots=[
                "depth_1:depth@track_index=0,side=inside,w=52px",
                "gc_content:dinucleotide_content@side=inside,w=42px,nt=GC",
            ],
            cfg=apply_config_overrides(None, {"canvas.show_skew": False, **DOCS_PERCENT_AXIS, **_depth_overrides(
                **{"canvas.show_gc": True})}),
            window=100,
            step=100,
        ),
    ).tostring()

    root = ET.fromstring(svg)
    ranges = _range_labels(root)
    assert sorted(element.text or "" for element in ranges) == ["[0\u201362x]", "[10\u201355%]"]
    # Text along a path, on its alphabetic baseline: CairoSVG ignores dominant-baseline there.
    assert all(element.tag == f"{SVG}textPath" and element.get("dominant-baseline") is None for element in ranges)
    assert not any(
        re.fullmatch(r"[0-9.]+[x%]", text.text or "") for text in root.iter(f"{SVG}text")
    )
    # Every large tick is still drawn, the bounds included: 10, 20, .., 50, 55 % and 0, 20, 40, 60, 62 x, plus the
    # axis lines.
    axis_lines = [
        len(group.findall(f"{SVG}line")) for group in root.iter(f"{SVG}g")
        if group.find(f"{SVG}line") is not None and group.find(f"{SVG}text") is not None
    ]
    assert sorted(axis_lines) == [1 + 5, 1 + 6]


@pytest.mark.linear
@pytest.mark.parametrize(("log_scale", "expected"), [(False, "[0\u201362x]"), (True, "[0\u201362x, log]")])
def test_linear_depth_axis_draws_one_range_label_left_of_the_axis_centred_on_the_track(
    log_scale: bool, expected: str
) -> None:
    svg = assemble_linear_diagram_from_records(
        [_make_record()],
        LinearAssemblyInputs(legend="none", depth_table=_depth_table(), window=100, step=100,
                             cfg=apply_config_overrides(None, _depth_overrides(
                                 **{"objects.depth.normalize": log_scale}))),
    ).tostring()

    root = ET.fromstring(svg)
    (label,) = _range_labels(root)
    assert label.text == expected
    assert label.get("text-anchor") == "end"
    axis_group = next(
        group for group in root.iter(f"{SVG}g")
        if any(text is label for text in group.findall(f"{SVG}text"))
    )
    vertical = next(
        line for line in axis_group.findall(f"{SVG}line") if line.get("x1") == line.get("x2")
    )
    top, bottom = float(vertical.get("y1", "nan")), float(vertical.get("y2", "nan"))
    assert float(label.get("y", "nan")) == pytest.approx(0.5 * (top + bottom))
    assert float(label.get("x", "nan")) < float(vertical.get("x1", "nan"))


@pytest.mark.circular
def test_a_circular_depth_axis_with_the_range_label_reserves_no_room_outside_its_band() -> None:
    from types import SimpleNamespace

    from gbdraw.layout.circular_depth_axis import resolve_depth_axis_footprint

    def footprint(axis_labels: str) -> float:
        config = SimpleNamespace(show_axis=True, show_ticks=True, tick_font_size=10.0, axis_labels=axis_labels)
        return resolve_depth_axis_footprint(config, 40.0).radial_outer_extra_px

    assert footprint("numbered") == pytest.approx(5.0)
    assert footprint("range") < 1.0


V0145_SCALAR_AXES_SESSION = (
    Path(__file__).parent / "fixtures" / "sessions" / "REGEX-scalar-axes.v50-cli0145.gbdraw-session.json.gz"
)


@pytest.mark.circular
def test_a_session_written_before_the_label_style_keeps_its_numbered_ticks(tmp_path: Path) -> None:
    """A full config the v0.14.5 CLI wrote has no ``axis_labels``; its axes drew numbers, and replay keeps them."""

    from gbdraw.circular import circular_main

    circular_main(["--session", str(V0145_SCALAR_AXES_SESSION), "-o", str(tmp_path / "replay"), "-f", "svg"])
    svg = (tmp_path / "replay.svg").read_text(encoding="utf-8")

    assert [text for text, _, _ in _axis_labels(svg, "x")] == ["0x", "20x", "40x", "60x"]
    assert [text for text, _, _ in _axis_labels(svg, "%")] == ["20%", "30%", "40%", "50%", "60%", "70%"]
    assert _range_labels(ET.fromstring(svg)) == []
