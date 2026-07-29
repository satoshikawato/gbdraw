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
from gbdraw.configurators import LegendDrawingConfigurator, LegendMeasurement
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
        color_table=None,
        default_colors=None,
        selected_features_set=["CDS"],
        profile=profile,
        gc_config=None,
        skew_config=None,
        feature_config=None,
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

    measurement = configurator.measure_legend(
        _legend_table(),
        canvas_config,
    )

    assert isinstance(measurement, LegendMeasurement)
    assert measurement.circular_layout is not None
    assert measurement.legend_width == measurement.circular_layout.width
    assert measurement.legend_height == measurement.circular_layout.height
    assert measurement.has_gradient is True
    assert not hasattr(configurator, "legend_width")

    group = CircularLegendGroup(
        canvas_config,
        measurement,
        _legend_table(),
    )
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
            },
            "Inverted": {
                "type": "gradient",
                "min_color": "#ffffff",
                "max_color": "#ff0000",
                "stroke": "none",
                "width": 0,
                "min_value": 70,
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

    measurement = configurator.measure_legend(
        legend_table,
        canvas_config,
    )
    group = LinearLegendGroup(
        canvas_config,
        measurement,
        legend_table,
        cfg=profile.config,
    )

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
    with pytest.raises(FrozenInstanceError):
        measurement.linear_layout.active.width = 1


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

    measurement = configurator.measure_legend({}, canvas_config)

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
    with pytest.raises(FrozenInstanceError):
        measurement.legend_width = 1


def test_legend_mode_dispatch_does_not_use_the_canvas_class_name() -> None:
    profile = LinearRenderProfile(_config())
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1_000,
        profile=profile,
        legend="right",
    )
    configurator = _configurator(profile, canvas_config)
    misleading_canvas_type = type(
        "CircularCanvasConfigurator",
        (),
        {
            "legend_position": "right",
            "total_width": canvas_config.total_width,
        },
    )

    measurement = configurator.measure_legend(
        _legend_table(),
        misleading_canvas_type(),
    )

    assert measurement.circular_layout is None
    assert measurement.legend_height > 0
