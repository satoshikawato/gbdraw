from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from gbdraw.config.models import (
    CircularRenderProfile,
    GbdrawConfig,
    LinearRenderProfile,
)
from gbdraw.config.toml import load_config_toml
from gbdraw.api.config import apply_config_overrides
from gbdraw.core.text import drawn_text_box
from gbdraw.configurators import LegendDrawingConfigurator, LegendMeasurement
from gbdraw.configurators.legend import measure_legend
from gbdraw.layout.spatial import Aabb, union_aabbs
from gbdraw.render.groups.circular.legend import LegendGroup as CircularLegendGroup
from gbdraw.render.groups.linear.legend import LegendGroup as LinearLegendGroup


def _config() -> GbdrawConfig:
    return GbdrawConfig.from_dict(
        load_config_toml("gbdraw.data", "config.toml")
    )


def _record() -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * 250), id="legend-test")
    record.annotations["molecule_type"] = "DNA"
    return record


def _configurator(profile, canvas_config) -> LegendDrawingConfigurator:
    return LegendDrawingConfigurator(
        profile=profile,
        canvas_config=canvas_config,
    )


def _legend_table() -> dict[str, dict[str, object]]:
    return {
        "CDS": {
            "type": "solid",
            "fill": "#54bcf8",
            "stroke": "none",
            "width": 0,
        },
        "Pairwise match identity": {
            "type": "gradient",
            "min_color": "#ffffff",
            "max_color": "#000000",
            "stroke": "none",
            "width": 0,
            "min_value": 70,
            "max_value": 100,
        },
    }


def _circular_legend_table() -> dict[str, dict[str, object]]:
    """A feature row and one similarity-ring row, as the Circular assembler writes them."""

    return {
        "CDS": _legend_table()["CDS"],
        "Reference A": {
            "type": "ring",
            "source_index": 0,
            "track_index": 1,
            "fill": None,
            "stroke": "none",
            "width": 0,
            "min_color": "#f0f1f5",
            "max_color": "#8b9cc1",
            "min_value": 70,
            "max_value": 100,
        },
    }


def test_circular_legend_measurement_carries_the_render_layout() -> None:
    profile = CircularRenderProfile(_config())
    canvas_config = CircularCanvasConfigurator(
        output_prefix="out",
        profile=profile,
        legend="right",
        gb_record=_record(),
    )
    configurator = _configurator(profile, canvas_config)

    measurement = measure_legend(
        configurator,
        _circular_legend_table(),
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )

    assert isinstance(measurement, LegendMeasurement)
    assert measurement.circular_layout is not None
    assert measurement.legend_width == measurement.circular_layout.width
    assert measurement.legend_height == measurement.circular_layout.height
    assert measurement.local_bounds.min_x == pytest.approx(0.0)
    assert measurement.local_bounds.min_y == pytest.approx(
        -0.5 * measurement.color_rect_size
    )
    assert measurement.local_bounds.max_x == pytest.approx(
        measurement.legend_width
    )
    assert measurement.local_bounds.max_y == pytest.approx(
        measurement.legend_height - (0.5 * measurement.color_rect_size)
    )
    assert measurement.has_gradient is True
    assert not hasattr(configurator, "legend_width")

    group = CircularLegendGroup(measurement, _circular_legend_table())
    assert group.layout is measurement.circular_layout


@pytest.mark.parametrize("legend_position", ["bottom", "right"])
@pytest.mark.parametrize(
    "legend_table",
    [
        {
            "CDS": {
                "type": "solid",
                "fill": "#54bcf8",
                "stroke": "none",
                "width": 0,
            },
        },
        {
            "Pairwise match identity": {
                "type": "gradient",
                "min_color": "#ffffff",
                "max_color": "#000000",
                "stroke": "none",
                "width": 0,
                "min_value": 70,
                "max_value": 100,
            },
        },
        {
            "Collinear": {
                "type": "gradient",
                "min_color": "#ffffff",
                "max_color": "#00ffff",
                "stroke": "none",
                "width": 0,
                "min_value": 70,
                "max_value": 100,
            },
            "Inverted": {
                "type": "gradient",
                "min_color": "#ffffff",
                "max_color": "#ff0000",
                "stroke": "none",
                "width": 0,
                "min_value": 70,
                "max_value": 100,
            },
        },
        _legend_table(),
    ],
    ids=["solid-only", "single-gradient-only", "multi-gradient-only", "mixed"],
)
def test_linear_legend_measurement_is_the_renderer_geometry_authority(
    legend_position: str,
    legend_table: dict[str, dict[str, object]],
) -> None:
    profile = LinearRenderProfile(_config())
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1_000,
        profile=profile,
        legend=legend_position,
    )
    configurator = _configurator(profile, canvas_config)

    measurement = measure_legend(
        configurator,
        legend_table,
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )
    group = LinearLegendGroup(measurement, legend_table)

    assert measurement.circular_layout is None
    assert measurement.linear_layout is not None
    assert group.layout is measurement.linear_layout
    assert group.get_horizontal_dimensions() == pytest.approx(
        (
            measurement.linear_layout.horizontal.width,
            measurement.linear_layout.horizontal.height,
        )
    )
    assert group.get_vertical_dimensions() == pytest.approx(
        (
            measurement.linear_layout.vertical.width,
            measurement.linear_layout.vertical.height,
        )
    )
    assert group.legend_width == pytest.approx(
        measurement.linear_layout.active.width
    )
    assert group.legend_height == pytest.approx(
        measurement.linear_layout.active.height
    )
    assert measurement.legend_width == pytest.approx(group.legend_width)
    assert measurement.legend_height == pytest.approx(group.legend_height)
    assert measurement.local_bounds.min_x == pytest.approx(0.0)
    assert measurement.local_bounds.min_y == pytest.approx(0.0)
    assert measurement.local_bounds.max_x == pytest.approx(group.legend_width)
    assert measurement.local_bounds.max_y == pytest.approx(group.legend_height)
    with pytest.raises(FrozenInstanceError):
        measurement.linear_layout.active.width = 1


@pytest.mark.parametrize("feature_row", [False, True], ids=["gradient-only", "with-feature-row"])
@pytest.mark.parametrize("legend_position", ["left", "right", "top", "bottom"])
def test_single_gradient_title_wider_than_its_bar_stays_inside_the_legend(
    legend_position: str, feature_row: bool
) -> None:
    """The Linear block is as wide as its title, with the bar centred under it (OV-1640).

    The Circular similarity-ring block has no gradient title; its long names are
    tested in ``tests/test_circular_conservation.py`` (OV-1261, OV-1611, OV-1612).
    """
    font_size = 32.0
    cfg = apply_config_overrides(
        load_config_toml("gbdraw.data", "config.toml"),
        {
            "objects.legends.font_size.short": font_size,
            "objects.legends.font_size.long": font_size,
        },
    )
    table: dict[str, dict[str, object]] = {
        **({"CDS": {"type": "solid", "fill": "#54bcf8", "stroke": "none", "width": 0}} if feature_row else {}),
        "Pairwise match identity": {
            "type": "gradient",
            "min_color": "#ffffff",
            "max_color": "#000000",
            "stroke": "none",
            "width": 0,
            "min_value": 70,
            "max_value": 100,
        },
    }
    profile = LinearRenderProfile(cfg)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1, longest_genome=1_000, profile=profile, legend=legend_position
    )
    configurator = _configurator(profile, canvas_config)
    measurement = measure_legend(
        configurator, table, placement=canvas_config.legend_position, wrap_width=canvas_config.total_width
    )
    title_width = drawn_text_box("Pairwise match identity", configurator.font_family, font_size)[0]
    assert measurement.linear_layout is not None
    active = measurement.linear_layout.active
    gradient = active.gradient
    assert gradient is not None
    entry_x = active.gradient_x
    (single,) = gradient.entries
    title_x, bar_x = single.title_x, single.bar_x
    assert (gradient.min_label_x, gradient.max_label_x) == (bar_x, bar_x + gradient.bar_width)
    assert title_width > gradient.bar_width
    assert gradient.width == title_width
    assert bar_x == (title_width - gradient.bar_width) / 2.0
    assert title_x == title_width / 2.0
    # The title's measured extent lies inside the Legend's local bounds.
    assert measurement.local_bounds.min_x <= entry_x + title_x - (title_width / 2.0)
    assert entry_x + title_x + (title_width / 2.0) <= measurement.local_bounds.max_x


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_empty_legend_measurement_is_an_immutable_zero_value(mode: str) -> None:
    cfg = _config()
    if mode == "circular":
        profile = CircularRenderProfile(cfg)
        canvas_config = CircularCanvasConfigurator(
            output_prefix="out",
            profile=profile,
            legend="right",
            gb_record=_record(),
        )
    else:
        profile = LinearRenderProfile(cfg)
        canvas_config = LinearCanvasConfigurator(
            num_of_entries=1,
            longest_genome=1_000,
            profile=profile,
            legend="right",
        )
    configurator = _configurator(profile, canvas_config)

    measurement = measure_legend(
        configurator,
        {},
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )

    assert measurement.legend_width == 0
    assert measurement.legend_height == 0
    assert measurement.total_feature_legend_width == 0
    assert measurement.pairwise_legend_width == 0
    assert measurement.num_of_lines == 0
    assert measurement.num_of_columns == 0
    assert measurement.num_of_items_per_line == 0
    assert measurement.has_gradient is False
    assert measurement.circular_layout is None
    assert measurement.linear_layout is None
    assert measurement.local_bounds.min_x == 0.0
    assert measurement.local_bounds.min_y == 0.0
    assert measurement.local_bounds.max_x == 0.0
    assert measurement.local_bounds.max_y == 0.0
    with pytest.raises(FrozenInstanceError):
        measurement.legend_width = 1


def test_legend_mode_dispatch_uses_the_render_profile() -> None:
    profile = LinearRenderProfile(_config())
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1_000,
        profile=profile,
        legend="right",
    )
    configurator = _configurator(profile, canvas_config)
    measurement = measure_legend(
        configurator,
        _legend_table(),
        placement="right",
        wrap_width=canvas_config.total_width,
    )

    assert measurement.circular_layout is None
    assert measurement.legend_height > 0


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_legend_local_bounds_include_stroke_at_its_painted_position(mode: str) -> None:
    cfg = _config()
    if mode == "circular":
        profile = CircularRenderProfile(cfg)
        canvas_config = CircularCanvasConfigurator(
            output_prefix="out",
            profile=profile,
            legend="right",
            gb_record=_record(),
        )
    else:
        profile = LinearRenderProfile(cfg)
        canvas_config = LinearCanvasConfigurator(
            num_of_entries=1,
            longest_genome=1_000,
            profile=profile,
            legend="right",
        )
    legend_table = _legend_table()
    legend_table["CDS"]["stroke"] = "#222222"
    legend_table["CDS"]["width"] = 6
    measurement = measure_legend(
        _configurator(profile, canvas_config),
        legend_table,
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )

    rect_size = float(measurement.color_rect_size)
    if mode == "circular":
        assert measurement.circular_layout is not None
        layout = measurement.circular_layout
        entry = layout.solid_entries[0]
        base = Aabb(
            0.0,
            -0.5 * rect_size,
            float(layout.width),
            float(layout.height) - (0.5 * rect_size),
        )
        rect_x = float(entry.rect_x)
        rect_y = float(entry.rect_y)
    else:
        assert measurement.linear_layout is not None
        layout = measurement.linear_layout.active
        entry = layout.feature.entries[0]
        base = Aabb(0.0, 0.0, float(layout.width), float(layout.height))
        rect_x = float(layout.feature_x) + float(entry.rect_x)
        rect_y = float(layout.feature_y) + float(entry.rect_y)
    stroked_rect = Aabb(
        rect_x,
        rect_y - (0.5 * rect_size),
        rect_x + rect_size,
        rect_y + (0.5 * rect_size),
    ).expanded(3.0)
    assert measurement.local_bounds == union_aabbs((base, stroked_rect))


@pytest.mark.parametrize("wrap_width", [-1.0, float("inf"), float("nan")])
def test_legend_measurement_rejects_invalid_wrap_width(wrap_width: float) -> None:
    profile = CircularRenderProfile(_config())
    canvas_config = CircularCanvasConfigurator(
        output_prefix="out",
        profile=profile,
        legend="right",
        gb_record=_record(),
    )

    with pytest.raises(ValueError, match="wrap width"):
        measure_legend(
            _configurator(profile, canvas_config),
            _legend_table(),
            placement="right",
            wrap_width=wrap_width,
        )
