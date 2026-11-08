"""Scalar-axis tick labels (GC content in percent, depth) never overlap; the Percent band is light (OV-2759, OV-2763)."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET

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
            cfg=apply_config_overrides(None, {"canvas.show_skew": False, **DOCS_PERCENT_AXIS}),
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
                                              **DOCS_PERCENT_AXIS}),
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
