from __future__ import annotations

import copy
from pathlib import Path
import re
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.api import LinearMultiRecordOptions, save_figure_to
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.toml import load_config_toml
from gbdraw.core import text as text_module
from gbdraw.diagrams.linear import precalc as linear_precalc
from gbdraw.diagrams.linear import assemble as linear_assemble
from gbdraw.diagrams.linear.builders import add_record_definition_group
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
    config_dict: dict | None = None,
) -> SimpleNamespace:
    cfg = GbdrawConfig.from_dict(
        config_dict
        if config_dict is not None
        else load_config_toml("gbdraw.data", "config.toml")
    )
    return SimpleNamespace(
        canvas_padding=50,
        definition_gap=definition_gap,
        horizontal_offset=horizontal_offset,
        keep_definition_left_aligned=keep_definition_left_aligned,
        length_param="short",
        profile=LinearRenderProfile(cfg),
    )


def _definition_translate_x(canvas: Drawing) -> float:
    definition_group = next(
        element for element in canvas.elements if element.attribs.get("id") == "record_a_definition"
    )
    translations = re.findall(
        r"translate\(\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)\s*\)",
        str(definition_group.attribs["transform"]),
    )
    assert translations
    return sum(float(x_value) for x_value, _y_value in translations)


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
    return DefinitionGroup(
        record,
        canvas_config,
        cfg=GbdrawConfig.from_dict(config_dict),
    ).definition_bounding_box_width


def _definition_only_config(*, font_weight: str = "normal", definition_gap: float = 20) -> dict:
    config_dict = copy.deepcopy(load_config_toml("gbdraw.data", "config.toml"))
    config_dict["canvas"]["show_gc"] = False
    config_dict["canvas"]["show_skew"] = False
    config_dict["canvas"]["show_depth"] = False
    config_dict["labels"]["linear"]["scope"] = "none"
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
      const record = Array.from(
        svg?.querySelectorAll('g[data-gbdraw-record-id]') || []
      ).find(
        (element) => element.getAttribute('data-gbdraw-record-id') === recordId
      ) || document.getElementById(recordId);
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
    cfg = GbdrawConfig.from_dict(
        _definition_only_config(
            font_weight=font_weight,
            definition_gap=definition_gap,
        )
    )
    return assemble_linear_diagram_from_records(
        [record],
        cfg=cfg,
        selected_features_set=[],
        output_prefix=output_prefix,
        legend="none",
    )


@pytest.mark.linear
def test_linear_definition_group_follows_record_offset_by_default() -> None:
    canvas_a = Drawing()
    canvas_b = Drawing()

    add_record_definition_group(
        canvas_a,
        _record(),
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=_canvas_config(keep_definition_left_aligned=False),
        max_def_width=0,
    )
    add_record_definition_group(
        canvas_b,
        _record(),
        record_offset_y=10,
        record_offset_x=30,
        canvas_config=_canvas_config(keep_definition_left_aligned=False),
        max_def_width=0,
    )

    assert _definition_translate_x(canvas_b) == pytest.approx(_definition_translate_x(canvas_a) + 30)


def _definition_row_canvas(
    row_sizes: tuple[int, int], locked: bool, align_center: bool,
    *, subtitles: bool = True, show_replicon: bool = False, text_anchor: str = "middle",
) -> Drawing:
    records = []
    positions = []
    for row, count in enumerate(row_sizes):
        for _ in range(count):
            record = _record(
                "Aeromonas hydrophila" if row == 0 else "Aeromonas sp.",
                f"record_{len(records)}",
            )
            if subtitles:
                record.annotations["gbdraw_record_subtitle"] = "A1" if row == 0 else "B"
            # Unequal row lengths make centered sequence alignment move one row.
            record.seq = Seq("ATGC" * (250 if row == 0 else 150))
            record.features = [SeqFeature(
                FeatureLocation(0, len(record)), type="source",
                qualifiers={"plasmid": [f"p{len(records)}"]},
            )]
            records.append(record)
            positions.append(f"#{len(records)}@{row + 1}")
    config = _definition_only_config()
    config["canvas"]["linear"].update(
        keep_definition_left_aligned=locked, align_center=align_center,
    )
    config["objects"]["definition"]["linear"].update(
        show_replicon=show_replicon, text_anchor=text_anchor,
    )
    return assemble_linear_diagram_from_records(
        records, cfg=GbdrawConfig.from_dict(config), selected_features_set=[],
        layout=LinearMultiRecordOptions(multi_record_positions=tuple(positions)),
        legend="none",
    )


@pytest.mark.parametrize("row_sizes", [(1, 1), (2, 2), (1, 2)])
@pytest.mark.parametrize("locked", [False, True])
@pytest.mark.parametrize("align_center", [False, True])
def test_definition_column_aligns_long_and_short_rows(
    row_sizes: tuple[int, int], locked: bool, align_center: bool,
) -> None:
    drawing = _definition_row_canvas(row_sizes, locked, align_center)
    root = ET.fromstring(drawing.tostring())
    ns = {"s": "http://www.w3.org/2000/svg"}
    headings = [g for g in root.findall("s:g", ns)
                if g.find("s:text[@data-definition-line-kind='name']", ns) is not None]
    assert len(headings) == 2
    assert [["".join(t.itertext()) for t in g.findall("s:text", ns)] for g in headings] == [
        ["Aeromonas hydrophila", "A1"], ["Aeromonas sp.", "B"],
    ]
    assert {t.get("text-anchor") for g in headings for t in g.findall("s:text", ns)} == {
        "start" if locked else "middle"
    }

    def x(group: ET.Element) -> float:
        return sum(float(v) for v in re.findall(r"translate\(\s*([-+0-9.eE]+)", group.get("transform", "")))

    axes = [next(g for g in root.findall("s:g", ns)
                 if g.get("data-gbdraw-record-id") == rid and g.get("data-gbdraw-role") != "record-definition"
                 and g.get("data-gbdraw-role") != "record-definition-row")
            for rid in ("record_0", f"record_{row_sizes[0]}")]
    expected_shift = 0.0 if locked else x(axes[1]) - x(axes[0])
    assert x(headings[1]) - x(headings[0]) == pytest.approx(expected_shift)


@pytest.mark.parametrize("anchor", ["start", "end"])
@pytest.mark.parametrize("row_sizes", [(1, 1), (2, 2), (1, 2)])
@pytest.mark.parametrize("locked", [False, True])
def test_explicit_anchor_keeps_its_existing_single_record_scope(anchor, row_sizes, locked) -> None:
    drawing = _definition_row_canvas(row_sizes, locked, False, text_anchor=anchor)
    root = ET.fromstring(drawing.tostring())
    anchors = {text.get("text-anchor") for text in root.iter("{http://www.w3.org/2000/svg}text")
               if text.get("data-definition-line-kind")}
    assert anchors == {"start" if locked else anchor if row_sizes == (1, 1) else "middle"}


@pytest.mark.browser
@pytest.mark.parametrize("row_sizes", [(1, 1), (2, 2), (1, 2)])
@pytest.mark.parametrize("locked", [False, True])
@pytest.mark.parametrize("subtitles", [False, True])
def test_browser_definition_column_matches_collision_bounds(
    monkeypatch, row_sizes, locked, subtitles,
) -> None:
    playwright_api = pytest.importorskip("playwright.sync_api")
    measured = []
    original = linear_assemble._record_collision_bands

    def capture(**kwargs):
        bands = original(**kwargs)
        measured.append(bands)
        return bands

    monkeypatch.setattr(linear_assemble, "_record_collision_bands", capture)
    drawing = _definition_row_canvas(
        row_sizes, locked, True, subtitles=subtitles, show_replicon=True,
    )
    # The assembler's last pass uses final record positions. Compare those
    # domains to actual SVG text, without copying the placement formula.
    final_bands = measured[-sum(row_sizes):]
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page()
        page.set_content(drawing.tostring())
        result = page.evaluate("""async () => {
          await document.fonts.ready;
          const svg = document.querySelector('svg');
          const matrix = (el) => svg.getCTM().inverse().multiply(el.getCTM());
          const origin = (el) => new DOMPoint(0, 0).matrixTransform(matrix(el)).x;
          const box = (el) => {
            const b = el.getBBox(), m = matrix(el);
            const p = new DOMPoint(b.x, b.y).matrixTransform(m);
            const q = new DOMPoint(b.x + b.width, b.y + b.height).matrixTransform(m);
            return {left: p.x, right: q.x, top: p.y, bottom: q.y};
          };
          const defs = [...svg.querySelectorAll('g')].filter(g =>
            [...g.children].some(t => t.hasAttribute('data-definition-line-kind')));
          const axes = [...svg.querySelectorAll('g[data-gbdraw-record-id]')].filter(g =>
            !g.getAttribute('data-gbdraw-role')?.startsWith('record-definition'));
          return {
            width: svg.viewBox.baseVal.width || svg.width.baseVal.value,
            height: svg.viewBox.baseVal.height || svg.height.baseVal.value,
            axes: axes.map(g => ({id: g.dataset.gbdrawRecordId, x: origin(g)})),
            defs: defs.map(g => ({id: g.dataset.gbdrawRecordId, ...box(g),
              lines: [...g.querySelectorAll('text')].map(t => ({kind: t.dataset.definitionLineKind, ...box(t)}))
            }))
          };
        }""")
        browser.close()
    tolerance = 1.0  # Font glyph bearings may differ by at most one SVG pixel.
    assert len(result["axes"]) == sum(row_sizes)
    headings = []
    for index, bands in enumerate(final_bands):
        axis = next(a for a in result["axes"] if a["id"] == f"record_{index}")
        origin = axis["x"] - bands[0].x_start
        definitions = sorted((d for d in result["defs"] if d["id"] == f"record_{index}"), key=lambda d: d["left"])
        expected = sorted((b for b in bands if b.kind == "definition"), key=lambda b: b.x_start)
        assert len(definitions) == len(expected)
        for definition, band in zip(definitions, expected):
            assert definition["left"] == pytest.approx(origin + band.x_start, abs=tolerance)
            assert definition["right"] == pytest.approx(origin + band.x_end, abs=tolerance)
            assert definition["left"] >= 0
            assert definition["right"] <= result["width"]
            assert 0 <= definition["top"] < definition["bottom"] <= result["height"]
            names = [t for t in definition["lines"] if t["kind"] == "name"]
            if names:
                headings.append((names[0], axis["x"]))
                assert axis["x"] - definition["right"] >= 20 - tolerance
                for line in definition["lines"]:
                    coordinate = lambda b: b["left"] if locked else (b["left"] + b["right"]) / 2
                    assert coordinate(line) == pytest.approx(coordinate(names[0]), abs=tolerance)
    assert len(headings) == 2
    (first, axis_a), (second, axis_b) = headings
    coordinate = lambda b: b["left"] if locked else (b["left"] + b["right"]) / 2
    assert coordinate(second) - coordinate(first) == pytest.approx(0 if locked else axis_b - axis_a, abs=tolerance)


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
        max_def_width=max_def_width,
    )
    add_record_definition_group(
        canvas_b,
        long_record,
        record_offset_y=10,
        record_offset_x=45,
        canvas_config=canvas_config,
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
        max_def_width=definition_width,
    )
    add_record_definition_group(
        canvas_b,
        record,
        record_offset_y=10,
        record_offset_x=0,
        canvas_config=shifted_config,
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

    canvas_config = _canvas_config(keep_definition_left_aligned=True)
    max_width, measured_widths, heights = linear_precalc._precalculate_definition_metrics(
        records,
        canvas_config,
        cfg=canvas_config.profile.config,
    )

    assert max_width == 13
    assert heights == [10.0, 10.0]
    assert measured_widths == [12.01, 7.5]


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
        profile=LinearRenderProfile(cfg),
        legend="none",
    )
    definition_height = DefinitionGroup(
        record,
        canvas_config,
        cfg=cfg,
    ).definition_bounding_box_height

    drawing = assemble_linear_diagram_from_records(
        [record],
        cfg=cfg,
        selected_features_set=["CDS"],
        legend="none",
    )
    elements = {
        element.attribs.get("id"): element
        for element in drawing.elements
        if element.attribs.get("id")
    }

    def translate_y(element) -> float:
        translations = re.findall(
            r"translate\(\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)\s*\)",
            str(element.attribs["transform"]),
        )
        assert translations
        return sum(float(y_value) for _x_value, y_value in translations)

    record_group = next(
        element
        for element in drawing.elements
        if element.attribs.get("data-gbdraw-record-id") == "record_a"
    )
    axis_y = translate_y(record_group)
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

    DefinitionGroup(
        _record("", "record_a"),
        _canvas_config(
            keep_definition_left_aligned=True,
            config_dict=config_dict,
        ),
        cfg=GbdrawConfig.from_dict(config_dict),
    )

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
        _canvas_config(
            keep_definition_left_aligned=True,
            config_dict=config_dict,
        ),
        cfg=GbdrawConfig.from_dict(config_dict),
    )

    assert definition_group.definition_bounding_box_width == pytest.approx(62.0)
    assert ("Alpha ", "bold", "normal") in calls
    assert ("Beta", "bold", "italic") in calls
    assert (" Gamma", "bold", "normal") in calls


@pytest.mark.linear
@pytest.mark.browser
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

    save_figure_to(canvas, ["svg", "png", "pdf", "eps", "ps"])

    for suffix in (".svg", ".png", ".pdf", ".eps", ".ps"):
        output_path = tmp_path / f"linear_definition_gap{suffix}"
        assert output_path.exists()
        assert output_path.stat().st_size > 0
