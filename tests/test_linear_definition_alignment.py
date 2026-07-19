from __future__ import annotations

import copy
from pathlib import Path
import re
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.api import assemble_linear_diagram_from_records
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.core import text as text_module
from gbdraw.diagrams.linear import precalc as linear_precalc
from gbdraw.diagrams.linear.builders import add_record_definition_group
from gbdraw.render.export import save_figure
from gbdraw.render.groups.linear import DefinitionGroup


def _record(label: str = "Record A", record_id: str = "record_a") -> SeqRecord:
    record = SeqRecord(Seq("A" * 100), id=record_id)
    record.annotations["gbdraw_record_label"] = label
    return record


def _canvas_config(
    *,
    keep_definition_left_aligned: bool,
    definition_gap: float = 20,
    horizontal_offset: float = 120,
) -> SimpleNamespace:
    return SimpleNamespace(
        canvas_padding=50,
        definition_gap=definition_gap,
        horizontal_offset=horizontal_offset,
        keep_definition_left_aligned=keep_definition_left_aligned,
        length_param="short",
    )


def _definition_translate_x(canvas: Drawing) -> float:
    definition_group = next(
        element for element in canvas.elements if element.attribs.get("id") == "record_a_definition"
    )
    transform = str(definition_group.attribs["transform"])
    return float(transform.split("translate(", 1)[1].split(",", 1)[0])


def _definition_text_anchors(canvas: Drawing) -> set[str]:
    definition_group = next(
        element for element in canvas.elements if element.attribs.get("id") == "record_a_definition"
    )
    return {
        str(element.attribs.get("text-anchor"))
        for element in definition_group.elements
        if getattr(element, "elementname", "") == "text"
    }


def _definition_width(record: SeqRecord, config_dict: dict, canvas_config: SimpleNamespace) -> float:
    return DefinitionGroup(record, config_dict, canvas_config).definition_bounding_box_width


def _definition_only_config(*, font_weight: str = "normal", definition_gap: float = 20) -> dict:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    config_dict["canvas"]["show_gc"] = False
    config_dict["canvas"]["show_skew"] = False
    config_dict["canvas"]["show_depth"] = False
    config_dict["canvas"]["show_labels"] = False
    config_dict["canvas"]["linear"]["keep_definition_left_aligned"] = True
    config_dict["canvas"]["linear"]["definition_gap"] = definition_gap
    definition_cfg = config_dict["objects"]["definition"]["linear"]
    definition_cfg["font_weight"] = font_weight
    definition_cfg.pop("line_styles", None)
    definition_cfg["show_replicon"] = False
    definition_cfg["show_accession"] = False
    definition_cfg["show_length"] = False
    return config_dict


def _bbox_in_svg_space_script() -> str:
    return """
    async ({ definitionId, recordId }) => {
      await document.fonts.ready;
      const svg = document.querySelector('svg');
      const definition = document.getElementById(definitionId);
      const record = document.getElementById(recordId);
      if (!svg || !definition || !record) {
        throw new Error('Expected SVG, definition group, and record group to exist');
      }
      function bboxInSvgSpace(element) {
        const bbox = element.getBBox();
        const matrix = element.getCTM();
        const point = svg.createSVGPoint();
        const corners = [
          [bbox.x, bbox.y],
          [bbox.x + bbox.width, bbox.y],
          [bbox.x, bbox.y + bbox.height],
          [bbox.x + bbox.width, bbox.y + bbox.height],
        ].map(([x, y]) => {
          point.x = x;
          point.y = y;
          const transformed = point.matrixTransform(matrix);
          return { x: transformed.x, y: transformed.y };
        });
        return {
          left: Math.min(...corners.map((point) => point.x)),
          right: Math.max(...corners.map((point) => point.x)),
          top: Math.min(...corners.map((point) => point.y)),
          bottom: Math.max(...corners.map((point) => point.y)),
        };
      }
      const definitionBox = bboxInSvgSpace(definition);
      const recordBox = bboxInSvgSpace(record);
      return {
        gap: recordBox.left - definitionBox.right,
        definitionRight: definitionBox.right,
        recordLeft: recordBox.left,
      };
    }
    """


def _linear_definition_canvas(
    label: str,
    *,
    font_weight: str = "normal",
    definition_gap: float = 20,
    output_prefix: str = "linear_definition_gap",
) -> Drawing:
    record = _record(label, record_id="record_a")
    return assemble_linear_diagram_from_records(
        [record],
        config_dict=_definition_only_config(
            font_weight=font_weight,
            definition_gap=definition_gap,
        ),
        selected_features_set=[],
        output_prefix=output_prefix,
        legend="none",
    )


@pytest.mark.linear
def test_linear_definition_group_follows_record_offset_by_default() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        _record(),
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=_canvas_config(keep_definition_left_aligned=False),
        config_dict=config_dict,
        max_def_width=0,
    )
    add_record_definition_group(
        canvas_b,
        _record(),
        record_offset_y=10,
        record_offset_x=30,
        canvas_config=_canvas_config(keep_definition_left_aligned=False),
        config_dict=config_dict,
        max_def_width=0,
    )

    assert _definition_translate_x(canvas_b) == pytest.approx(_definition_translate_x(canvas_a) + 30)


@pytest.mark.linear
def test_linear_definition_group_can_stay_in_left_column() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = _canvas_config(keep_definition_left_aligned=True)
    short_record = _record("Short")
    long_record = _record("A much longer definition label", "record_a")
    max_def_width = max(
        _definition_width(short_record, config_dict, canvas_config),
        _definition_width(long_record, config_dict, canvas_config),
    )
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        short_record,
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=max_def_width,
    )
    add_record_definition_group(
        canvas_b,
        long_record,
        record_offset_y=10,
        record_offset_x=45,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=max_def_width,
    )

    assert _definition_translate_x(canvas_a) == pytest.approx(_definition_translate_x(canvas_b))
    assert _definition_text_anchors(canvas_a) == {"start"}
    assert _definition_text_anchors(canvas_b) == {"start"}


@pytest.mark.linear
def test_locked_linear_definition_column_uses_rendered_svg_text_width(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = _canvas_config(keep_definition_left_aligned=True)

    def overestimated_bbox(*_args, **_kwargs):
        return 100.0, 12.0

    def rendered_svg_bbox(*_args, **_kwargs):
        return 60.0, 12.0

    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_bbox_dimensions",
        overestimated_bbox,
    )
    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_svg_bbox_dimensions",
        rendered_svg_bbox,
    )

    definition_width = _definition_width(_record("Wide Label"), config_dict, canvas_config)
    canvas = Drawing()
    add_record_definition_group(
        canvas,
        _record("Wide Label"),
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=definition_width,
    )

    definition_x = _definition_translate_x(canvas)
    gap = canvas_config.horizontal_offset - (definition_x + definition_width)

    assert definition_width == pytest.approx(60.0)
    assert gap == pytest.approx(canvas_config.definition_gap)


@pytest.mark.linear
def test_locked_linear_definition_column_uses_configured_definition_gap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = _canvas_config(
        keep_definition_left_aligned=True,
        definition_gap=32,
    )

    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_bbox_dimensions",
        lambda *_args, **_kwargs: (64.0, 12.0),
    )
    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_svg_bbox_dimensions",
        lambda *_args, **_kwargs: (64.0, 12.0),
    )

    definition_width = _definition_width(_record("Wide Label"), config_dict, canvas_config)
    canvas = Drawing()
    add_record_definition_group(
        canvas,
        _record("Wide Label"),
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=definition_width,
    )

    definition_x = _definition_translate_x(canvas)
    gap = canvas_config.horizontal_offset - (definition_x + definition_width)

    assert gap == pytest.approx(32.0)


@pytest.mark.linear
def test_locked_linear_definition_column_follows_global_horizontal_shift() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    record = _record()
    base_config = _canvas_config(keep_definition_left_aligned=True, horizontal_offset=120)
    shifted_config = _canvas_config(keep_definition_left_aligned=True, horizontal_offset=150)
    definition_width = _definition_width(record, config_dict, base_config)
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        record,
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=base_config,
        config_dict=config_dict,
        max_def_width=definition_width,
    )
    add_record_definition_group(
        canvas_b,
        record,
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=shifted_config,
        config_dict=config_dict,
        max_def_width=definition_width,
    )

    assert _definition_translate_x(canvas_b) == pytest.approx(_definition_translate_x(canvas_a) + 30)
    shifted_gap = shifted_config.horizontal_offset - (_definition_translate_x(canvas_b) + definition_width)
    assert shifted_gap == pytest.approx(shifted_config.definition_gap)


@pytest.mark.linear
def test_locked_linear_definition_column_gap_uses_records_column_left() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    canvas_config = _canvas_config(keep_definition_left_aligned=True, horizontal_offset=120)
    record = _record("A much longer definition label")
    definition_width = _definition_width(record, config_dict, canvas_config)
    leftmost_record_offset_x = -10.0
    canvas = Drawing()

    add_record_definition_group(
        canvas,
        record,
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=canvas_config,
        config_dict=config_dict,
        max_def_width=definition_width - leftmost_record_offset_x,
    )

    records_column_left = canvas_config.horizontal_offset + leftmost_record_offset_x
    definition_right = _definition_translate_x(canvas) + definition_width
    assert records_column_left - definition_right == pytest.approx(canvas_config.definition_gap)


@pytest.mark.linear
def test_precalculated_max_definition_width_is_ceiled(monkeypatch: pytest.MonkeyPatch) -> None:
    records = [_record("A", "record_a"), _record("B", "record_b")]
    widths = {"record_a": 12.01, "record_b": 7.5}

    class FakeDefinitionGroup:
        def __init__(self, record, *_args, **_kwargs):
            self.definition_bounding_box_width = widths[record.id]
            self.definition_bounding_box_height = 10.0

    monkeypatch.setattr(linear_precalc, "DefinitionGroup", FakeDefinitionGroup)

    max_width, heights, half_heights = linear_precalc._precalculate_definition_metrics(
        records,
        load_config_toml("gbdraw.data", "config.toml"),
        _canvas_config(keep_definition_left_aligned=True),
    )

    assert max_width == 13
    assert heights == [10.0, 10.0]
    assert half_heights == [5.0, 5.0]


@pytest.mark.linear
def test_linear_definition_gap_defaults_to_twenty_for_legacy_config() -> None:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    del config_dict["canvas"]["linear"]["definition_gap"]

    cfg = GbdrawConfig.from_dict(config_dict)

    assert cfg.canvas.linear.definition_gap == pytest.approx(20.0)


@pytest.mark.linear
def test_linear_definition_gap_reads_explicit_config_value() -> None:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    config_dict["canvas"]["linear"]["definition_gap"] = 34

    cfg = GbdrawConfig.from_dict(config_dict)

    assert cfg.canvas.linear.definition_gap == pytest.approx(34.0)


@pytest.mark.linear
@pytest.mark.parametrize(("track_layout", "direction"), [("above", -1), ("below", 1)])
def test_linear_definition_band_matches_resolved_feature_center(
    track_layout: str,
    direction: int,
) -> None:
    record = _record("Definition follows resolved feature lane")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(FeatureLocation(20, 80, strand=1), type="CDS")
    ]
    config_dict = _definition_only_config()
    config_dict["canvas"]["linear"]["track_layout"] = track_layout
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    definition_height = DefinitionGroup(
        record,
        config_dict,
        canvas_config,
        cfg=cfg,
    ).definition_bounding_box_height

    drawing = assemble_linear_diagram_from_records(
        [record],
        config_dict=config_dict,
        selected_features_set=["CDS"],
        legend="none",
    )
    elements = {
        element.attribs.get("id"): element
        for element in drawing.elements
        if element.attribs.get("id")
    }

    def translate_y(element) -> float:
        match = re.search(
            r"translate\([^,]+,([-+0-9.eE]+)\)",
            str(element.attribs["transform"]),
        )
        assert match is not None
        return float(match.group(1))

    axis_y = translate_y(elements["record_a"])
    definition_center_y = translate_y(elements["record_a_definition"])
    half_height = 0.5 * definition_height
    canvas_band = drawing._gbdraw_track_slot_geometry["records"][0]["canvasBand"]

    assert (definition_center_y - axis_y) * direction > 0.0
    assert canvas_band["absoluteTopPx"] <= definition_center_y - half_height
    assert canvas_band["absoluteBottomPx"] >= definition_center_y + half_height


@pytest.mark.linear
def test_liberation_sans_bundled_fonts_preferred_over_system_fontconfig(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    text_module._font_path_cache.clear()
    monkeypatch.setattr(
        text_module.subprocess,
        "run",
        lambda *_args, **_kwargs: pytest.fail("fc-match should not be used for bundled Liberation Sans"),
    )

    cases = [
        ("normal", "normal", "LiberationSans-Regular.ttf"),
        ("bold", "normal", "LiberationSans-Bold.ttf"),
        ("normal", "italic", "LiberationSans-Italic.ttf"),
        ("bold", "italic", "LiberationSans-BoldItalic.ttf"),
    ]
    for font_weight, font_style, expected_name in cases:
        font_path = text_module._resolve_font_path(
            "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif",
            allow_system=True,
            font_weight=font_weight,
            font_style=font_style,
        )

        assert font_path is not None
        assert Path(font_path).name == expected_name


@pytest.mark.linear
def test_definition_plain_lines_use_definition_font_weight(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dict = _definition_only_config(font_weight="bold")
    config_dict["objects"]["definition"]["linear"]["show_accession"] = True
    calls: list[tuple[str, str, str]] = []

    def fake_bbox(text, _font_family, _font_size, _dpi, *, font_weight="normal", font_style="normal"):
        calls.append((str(text), str(font_weight), str(font_style)))
        return 10.0, 12.0

    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_bbox_dimensions",
        fake_bbox,
    )
    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_svg_bbox_dimensions",
        fake_bbox,
    )

    DefinitionGroup(_record("", "record_a"), config_dict, _canvas_config(keep_definition_left_aligned=True))

    assert ("record_a", "bold", "normal") in calls


@pytest.mark.linear
def test_definition_mixed_content_width_sums_style_aware_parts(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_dict = _definition_only_config(font_weight="bold")
    widths = {
        ("Alpha ", "bold", "normal"): 10.0,
        ("Beta", "bold", "italic"): 22.0,
        (" Gamma", "bold", "normal"): 30.0,
    }
    calls: list[tuple[str, str, str]] = []

    def fake_bbox(text, _font_family, _font_size, _dpi, *, font_weight="normal", font_style="normal"):
        key = (str(text), str(font_weight), str(font_style))
        calls.append(key)
        return widths.get(key, 1.0), 12.0

    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_bbox_dimensions",
        fake_bbox,
    )
    monkeypatch.setattr(
        "gbdraw.render.groups.linear.definition.calculate_svg_bbox_dimensions",
        fake_bbox,
    )

    definition_group = DefinitionGroup(
        _record("Alpha <i>Beta</i> Gamma"),
        config_dict,
        _canvas_config(keep_definition_left_aligned=True),
    )

    assert definition_group.definition_bounding_box_width == pytest.approx(62.0)
    assert ("Alpha ", "bold", "normal") in calls
    assert ("Beta", "bold", "italic") in calls
    assert (" Gamma", "bold", "normal") in calls


@pytest.mark.linear
@pytest.mark.parametrize(
    ("label", "font_weight"),
    [
        ("Plain definition", "normal"),
        ("Plain <i>italic definition</i> tail", "normal"),
        ("Bold definition", "bold"),
    ],
)
def test_browser_rendered_definition_gap_is_at_least_configured_gap(
    label: str,
    font_weight: str,
) -> None:
    playwright_sync_api = pytest.importorskip(
        "playwright.sync_api",
        reason="playwright is not available in this environment",
    )

    canvas = _linear_definition_canvas(label, font_weight=font_weight)
    svg_source = canvas.tostring()

    with playwright_sync_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1800, "height": 600})
        page.set_content(svg_source)
        result = page.evaluate(
            _bbox_in_svg_space_script(),
            {"definitionId": "record_a_definition", "recordId": "record_a"},
        )
        browser.close()

    assert result["gap"] >= 20.0


@pytest.mark.linear
def test_linear_definition_gap_svg_converts_with_cairosvg(tmp_path: Path) -> None:
    canvas = _linear_definition_canvas(
        "Plain <i>italic definition</i> tail",
        output_prefix=str(tmp_path / "linear_definition_gap"),
    )

    save_figure(canvas, ["svg", "png", "pdf", "eps", "ps"])

    for suffix in (".svg", ".png", ".pdf", ".eps", ".ps"):
        output_path = tmp_path / f"linear_definition_gap{suffix}"
        assert output_path.exists()
        assert output_path.stat().st_size > 0
