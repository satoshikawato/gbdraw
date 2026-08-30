from __future__ import annotations

from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.layout.spatial import Aabb
from gbdraw.diagrams.circular.positioning import place_definition_group_on_canvas
from gbdraw.render.groups.circular.definition import DefinitionGroup
from svgwrite.container import Group


def _definition(
    *,
    plot_title: str | None = None,
    profile: str = "full",
    record_label: str = "",
    record_subtitle: str = "",
) -> DefinitionGroup:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    render_profile = CircularRenderProfile(cfg)
    record = SeqRecord(Seq("ATGC" * 250), id="NC_TEST.1")
    if record_label:
        record.annotations["gbdraw_record_label"] = record_label
    if record_subtitle:
        record.annotations["gbdraw_record_subtitle"] = record_subtitle
    canvas = CircularCanvasConfigurator("test", render_profile, "none", record)
    return DefinitionGroup(
        record,
        canvas,
        cfg=cfg,
        species="<i>Testus boundsii</i>",
        strain="strain A",
        plot_title=plot_title,
        definition_profile=profile,
        definition_group_id=("plot_title" if profile == "shared_common" else None),
    )


def test_circular_definition_exposes_authoritative_local_bounds() -> None:
    definition = _definition()

    assert isinstance(definition.local_bounds, Aabb)
    assert definition.local_bounds.width > 0.0
    assert definition.local_bounds.height > 0.0
    for element in definition.get_group().elements:
        x = float(element.attribs["x"])
        y = float(element.attribs["y"])
        assert definition.local_bounds.min_x < x < definition.local_bounds.max_x
        assert definition.local_bounds.min_y < y < definition.local_bounds.max_y


def test_circular_plot_title_bounds_are_centered_in_local_coordinates() -> None:
    title = _definition(
        plot_title="Comparative <i>mitochondrial</i> rings",
        profile="shared_common",
    )

    assert title.get_group().attribs["data-gbdraw-role"] == "plot-title"
    assert title.local_bounds.min_x == pytest.approx(-title.local_bounds.max_x)
    assert title.local_bounds.width > 100.0
    assert title.local_bounds.height > 0.0


def test_circular_definition_local_bounds_center_on_record_axis() -> None:
    group = Group(debug=False)
    local_bounds = Aabb(-40.0, -10.0, 60.0, 50.0)
    setattr(group, "_gbdraw_local_bounds", local_bounds)

    place_definition_group_on_canvas(
        group,
        SimpleNamespace(offset_x=400.0, offset_y=300.0),
    )

    translate_x, translate_y = getattr(group, "_gbdraw_plot_translation")
    assert translate_x + 0.5 * (local_bounds.min_x + local_bounds.max_x) == pytest.approx(400.0)
    assert translate_y + 0.5 * (local_bounds.min_y + local_bounds.max_y) == pytest.approx(300.0)


def test_circular_record_label_and_subtitle_override_inferred_title_lines() -> None:
    definition = _definition(
        record_label="長い <i>環状ゲノム</i> ラベル",
        record_subtitle="Selected record subtitle",
    )

    svg = definition.get_group().tostring()

    assert "長い " in svg
    assert "環状ゲノム" in svg
    assert 'font-style="italic"' in svg
    assert "Selected record subtitle" in svg
    assert "Testus boundsii" not in svg
    assert "strain A" not in svg
    assert definition.local_bounds.width > 100.0


def test_empty_circular_record_title_values_preserve_inference() -> None:
    definition = _definition()
    svg = definition.get_group().tostring()

    assert "Testus boundsii" in svg
    assert "strain A" in svg
