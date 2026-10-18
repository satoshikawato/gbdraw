"""Per-Result Legend row facts reported to the Web (OV-63).

`drawn` is the Legend table the Result drew; `suppressed` is what the records can
name that the draft did not draw. A key in neither is stale.
"""

from __future__ import annotations

import inspect
import json
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    CircularBatchRequest,
    CircularDiagramOptions,
    CircularDiagramRequest,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularTrackOptions,
    ColorOptions,
    InMemoryRecordSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.api.request_render import render_request
from gbdraw.legend.row_facts import collect_legend_row_facts
from gbdraw.legend.sources import DepthLegendSource, NumericLegendSource, legend_slot_sources
from gbdraw.diagrams.circular.assemble import _sync_legend_table_for_circular_slots
from gbdraw.diagrams.linear.assemble import _sync_legend_table_for_linear_slots
from gbdraw.legend.table import LegendTable, _generated_legend_rows
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA, _decode_track_slots
from gbdraw.tracks.circular import parse_circular_track_slots
from gbdraw.tracks.linear import parse_linear_track_slots

_COLOR_COLUMNS = ["feature_type", "qualifier_key", "value", "color", "caption"]
_VISIBILITY_COLUMNS = ["record_id", "feature_type", "qualifier", "value", "action"]


def _record(record_id: str, *, repeats: bool = True) -> SeqRecord:
    record = SeqRecord(
        Seq("ATGC" * 500), id=record_id, annotations={"molecule_type": "DNA", "topology": "circular"}
    )
    record.features = [
        SeqFeature(SimpleLocation(100, 400, 1), type="CDS", qualifiers={"locus_tag": [f"{record_id}_1"]}),
        SeqFeature(SimpleLocation(600, 900, -1), type="CDS", qualifiers={"locus_tag": [f"{record_id}_2"]}),
        SeqFeature(SimpleLocation(0, 2000, 1), type="source", qualifiers={}),
    ]
    if repeats:
        record.features.append(
            SeqFeature(SimpleLocation(1200, 1300, 1), type="repeat_region", qualifiers={"rpt_type": ["tandem"]})
        )
    return record


def _color_table(*rows: tuple[str, str, str, str, str]) -> ColorOptions:
    return ColorOptions(color_table=pd.DataFrame(list(rows), columns=_COLOR_COLUMNS))


def _hide(record_id: str, feature_type: str) -> pd.DataFrame:
    return pd.DataFrame([[record_id, feature_type, "location", ".", "off"]], columns=_VISIBILITY_COLUMNS)


def _render(request, tmp_path: Path):
    result = render_request(request, include_feature_catalog=True)
    items = result.items if hasattr(result, "items") else (result,)
    out = []
    for index, item in enumerate(items):
        svg = Path(item.output_paths[0]).read_text(encoding="utf-8")
        facts = collect_legend_row_facts(item.drawing, result_index=index, result_name=Path(item.output_paths[0]).name)
        out.append((facts, list(dict.fromkeys(re.findall(r'data-legend-key="([^"]*)"', svg)))))
    return out


def _circular(tmp_path: Path, records=None, **options):
    output = RenderOutputRequest(output_prefix="facts", formats=("svg",), output_directory=tmp_path, overwrite=True)
    records = records or [_record("rec")]
    return CircularDiagramRequest(
        records=tuple(RecordInput(source=InMemoryRecordSource(r), record_key=r.id) for r in records),
        options=CircularDiagramOptions(**options),
        output=output,
    )


def test_drawn_rows_equal_the_svg_legend_keys_and_unrelated_keys_are_in_neither(tmp_path):
    [(facts, svg_keys)] = _render(_circular(tmp_path), tmp_path)
    assert facts["drawn"] == svg_keys
    assert {"CDS", "repeat_region"} <= set(facts["drawn"])
    assert not set(facts["drawn"]) & set(facts["suppressed"])
    # A type of the records with no drawn row is producible: its row may be absent.
    assert "source" in facts["suppressed"]
    # A caption no feature of the records can name is stale: in neither list.
    assert "Ghost" not in facts["drawn"] + facts["suppressed"]
    assert "no_such_type" not in facts["drawn"] + facts["suppressed"]


def test_hiding_every_feature_of_a_type_suppresses_its_row(tmp_path):
    [(facts, svg_keys)] = _render(
        _circular(tmp_path, feature_visibility_table=_hide("rec", "repeat_region")), tmp_path
    )
    assert facts["drawn"] == svg_keys
    assert "repeat_region" not in facts["drawn"]
    assert "repeat_region" in facts["suppressed"]
    assert "CDS" in facts["drawn"]


def test_a_track_whose_row_another_track_draws_reports_its_key_suppressed(tmp_path):
    """Two tracks of one label and one swatch share a row; the second keeps its key (0.13.0 Sessions hold it)."""
    slots = [
        "features:features",
        "gc_content:dinucleotide_content",
        "gc_content_2:dinucleotide_content",
        "gc_skew:dinucleotide_skew",
        "gc_skew_2:dinucleotide_skew@positive_color=#00ff00",
    ]
    [(facts, svg_keys)] = _render(
        _circular(tmp_path, tracks=CircularTrackOptions(circular_track_slots=slots)), tmp_path
    )
    assert facts["drawn"] == svg_keys
    assert {"GC content", "GC skew (+)", "GC skew (-)", "GC skew (+) (2)"} <= set(facts["drawn"])
    assert {"GC content (2)", "GC skew (-) (2)"} <= set(facts["suppressed"])
    assert not set(facts["drawn"]) & set(facts["suppressed"])
    svg = (tmp_path / "facts.svg").read_text(encoding="utf-8")
    names = dict(re.findall(r'data-legend-key="([^"]*)" data-legend-name="([^"]*)"', svg))
    # Only a row whose name is not its key carries the name.
    assert names == {"GC skew (+) (2)": "GC skew (+)"}
    for text in ["GC content", "GC skew (+) (1)", "GC skew (+) (2)", "GC skew (-)"]:
        assert f">{text}</text>" in svg
    assert ">GC content (1)</text>" not in svg and ">GC skew (-) (1)</text>" not in svg


def test_a_type_outside_the_selection_is_suppressed_not_stale(tmp_path):
    [(facts, _)] = _render(_circular(tmp_path, selected_features_set=("CDS",)), tmp_path)
    assert "repeat_region" not in facts["drawn"]
    assert "repeat_region" in facts["suppressed"]


def test_a_rule_that_recaptions_every_feature_suppresses_the_type_row(tmp_path):
    [(facts, svg_keys)] = _render(
        _circular(tmp_path, colors=_color_table(("CDS", "locus_tag", ".", "#c83366", "Zeta"))), tmp_path
    )
    assert facts["drawn"] == svg_keys
    assert "Zeta" in facts["drawn"]
    assert {"CDS", "other proteins"} <= set(facts["suppressed"])
    assert "other proteins" not in facts["drawn"]


def test_a_rule_row_whose_features_are_all_hidden_is_suppressed(tmp_path):
    [(facts, svg_keys)] = _render(
        _circular(
            tmp_path,
            colors=_color_table(("CDS", "locus_tag", ".", "#c83366", "Zeta")),
            feature_visibility_table=_hide("rec", "CDS"),
        ),
        tmp_path,
    )
    assert facts["drawn"] == svg_keys
    assert "Zeta" not in facts["drawn"]
    assert {"Zeta", "CDS"} <= set(facts["suppressed"])


def test_a_rule_that_matches_no_feature_names_no_row(tmp_path):
    [(facts, _)] = _render(
        _circular(tmp_path, colors=_color_table(("CDS", "locus_tag", "^never$", "#c83366", "Zeta"))), tmp_path
    )
    assert "Zeta" not in facts["drawn"] + facts["suppressed"]


def test_each_batch_result_reports_its_own_rows(tmp_path):
    records = [_record("one"), _record("two", repeats=False)]
    request = CircularBatchRequest(
        records=tuple(RecordInput(source=InMemoryRecordSource(r), record_key=r.id) for r in records),
        options=CircularDiagramOptions(),
        outputs=tuple(
            RenderOutputRequest(
                output_prefix=f"facts-{index}", formats=("svg",), output_directory=tmp_path, overwrite=True
            )
            for index in range(2)
        ),
    )
    (first, first_keys), (second, second_keys) = _render(request, tmp_path)
    assert first["drawn"] == first_keys and second["drawn"] == second_keys
    assert "repeat_region" in first["drawn"]
    assert "repeat_region" not in second["drawn"] + second["suppressed"]
    assert [first["resultIndex"], second["resultIndex"]] == [0, 1]


def test_a_multi_record_canvas_reports_one_set_of_rows(tmp_path):
    request = CircularDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(r), record_key=r.id)
            for r in (_record("one"), _record("two", repeats=False))
        ),
        options=CircularDiagramOptions(feature_visibility_table=_hide("one", "repeat_region")),
        layout=CircularMultiRecordOptions(),
        grouping="grid",
        output=RenderOutputRequest(output_prefix="canvas", formats=("svg",), output_directory=tmp_path, overwrite=True),
    )
    [(facts, svg_keys)] = _render(request, tmp_path)
    assert facts["drawn"] == svg_keys
    assert "CDS" in facts["drawn"]
    assert "repeat_region" in facts["suppressed"]


def test_linear_reports_its_rows(tmp_path):
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(r), record_key=r.id)
            for r in (_record("one"), _record("two"))
        ),
        options=LinearDiagramOptions(feature_visibility_table=pd.DataFrame(
            [["one", "repeat_region", "location", ".", "off"], ["two", "repeat_region", "location", ".", "off"]],
            columns=_VISIBILITY_COLUMNS,
        )),
        output=RenderOutputRequest(output_prefix="linear", formats=("svg",), output_directory=tmp_path, overwrite=True),
    )
    [(facts, svg_keys)] = _render(request, tmp_path)
    assert facts["drawn"] == svg_keys
    assert "repeat_region" in facts["suppressed"]
    assert "repeat_region" not in facts["drawn"]


def test_a_drawing_without_facts_reports_empty_lists():
    assert collect_legend_row_facts(object(), result_index=3, result_name="x.svg") == {
        "resultIndex": 3, "resultName": "x.svg", "drawn": [], "suppressed": [],
    }


def test_the_web_request_carries_the_facts_in_run_metadata(tmp_path):
    from gbdraw.session import build_session_document
    from gbdraw.web_support.request_render import render_embedded_canonical_web_request

    request = _circular(tmp_path, feature_visibility_table=_hide("rec", "repeat_region"))
    document = build_session_document(request).to_dict()
    response = render_embedded_canonical_web_request(
        document["renderRequest"], resources=document["resources"], workspace=tmp_path / "workspace"
    )
    [row] = response["metadata"]["legendRows"]
    assert row["resultName"] == response["results"][0]["name"]
    assert "repeat_region" in row["suppressed"] and "repeat_region" not in row["drawn"]
    keys = re.findall(r'data-legend-key="([^"]*)"', response["results"][0]["content"])
    assert list(dict.fromkeys(keys)) == row["drawn"]


def test_a_hidden_legend_draws_no_row_and_suppresses_what_the_records_name(tmp_path):
    request = CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(_record("rec")), record_key="rec"),),
        options=CircularDiagramOptions(output=CircularOutputOptions(legend="none")),
        output=RenderOutputRequest(output_prefix="none", formats=("svg",), output_directory=tmp_path, overwrite=True),
    )
    [(facts, svg_keys)] = _render(request, tmp_path)
    assert svg_keys == [] and facts["drawn"] == []
    assert {"CDS", "repeat_region"} <= set(facts["suppressed"])


def test_a_hidden_legend_suppresses_the_gc_rows_too(tmp_path):
    request = CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(_record("rec")), record_key="rec"),),
        options=CircularDiagramOptions(output=CircularOutputOptions(legend="none")),
        output=RenderOutputRequest(output_prefix="none", formats=("svg",), output_directory=tmp_path, overwrite=True),
    )
    [(facts, _)] = _render(request, tmp_path)
    assert {"GC content", "GC skew (+)", "GC skew (-)"} <= set(facts["suppressed"])


def test_a_multi_record_canvas_with_a_hidden_legend_still_reports_its_facts(tmp_path):
    request = CircularDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(r), record_key=r.id)
            for r in (_record("one"), _record("two", repeats=False))
        ),
        options=CircularDiagramOptions(output=CircularOutputOptions(legend="none")),
        layout=CircularMultiRecordOptions(),
        grouping="grid",
        output=RenderOutputRequest(output_prefix="canvas", formats=("svg",), output_directory=tmp_path, overwrite=True),
    )
    [(facts, svg_keys)] = _render(request, tmp_path)
    assert svg_keys == [] and facts["drawn"] == []
    assert {"CDS", "repeat_region", "GC content"} <= set(facts["suppressed"])


GENERATED_ROWS = json.loads(
    (Path(__file__).parent / "fixtures" / "legend_generated_rows.json").read_text(encoding="utf-8")
)["cases"]


@pytest.mark.parametrize("case", GENERATED_ROWS, ids=[case["name"] for case in GENERATED_ROWS])
def test_generated_legend_rows_are_the_rows_the_web_reads_from_their_keys(case):
    """The Web app reads these keys as Python's own rows (OV-307) and their palette keys (OV-309)."""
    def track(fills):
        return NumericLegendSource(case["dinucleotide"], fills[0], fills[1], "none", 0)

    default_colors = pd.DataFrame(
        [{"feature_type": key, "color": color} for key, color in (
            ("CDS", "#54bcf8"), ("tRNA", "#e8b441"), ("repeat_region", "#d3d3d3"), ("default", "#d3d3d3")
        )]
    )
    fills = _generated_legend_rows(
        {key: [tuple(row) for row in rows] for key, rows in case["feature_specific_colors"].items()},
        case["features_present"],
        default_colors,
        used_color_rules={tuple(row) for row in case["used_color_rules"]},
        default_used_features=set(case["default_used_features"]),
        gc=track(case["gc"]),
        skew=track(case["skew"]),
        depth=DepthLegendSource("", case["depth"], "none", 0),
        show_gc=case["show_gc"],
        show_skew=case["show_skew"],
        show_depth=case["show_depth"],
    )
    assert list(fills) == [key for key, _ in case["rows"]]


def test_generated_legend_rows_take_no_input_the_web_pin_lacks():
    """A new kind of generated row arrives as a new input; add it to the fixture and the Web reader."""
    assert list(inspect.signature(_generated_legend_rows).parameters) == [
        "feature_specific_colors", "features_present", "default_colors", "used_color_rules",
        "default_used_features", "gc", "skew", "depth", "show_gc", "show_skew", "show_depth",
        "block_stroke",
    ]


TRACK_SLOT_ROWS = json.loads(
    (Path(__file__).parent / "fixtures" / "legend_track_slot_rows.json").read_text(encoding="utf-8")
)["cases"]


@pytest.mark.parametrize("case", TRACK_SLOT_ROWS, ids=[case["name"] for case in TRACK_SLOT_ROWS])
def test_track_slot_legend_rows_are_the_rows_the_web_reads_from_the_slots(case):
    """The Web app reads each of these rows' slot and fill from the request's slot entries (OV-484)."""
    palette = case["palette"]

    def track(high, low):
        return NumericLegendSource(
            dinucleotide=case["dinucleotide"], high_fill_color=high, low_fill_color=low, stroke_color="none", stroke_width=0
        )

    gc_config = track(palette["gc_content"], palette["gc_content"])
    skew_config = track(palette["skew_high"], palette["skew_low"])
    # The request's entries (spec strings or slot objects), read as a render reads them.
    entries = list(_decode_track_slots(case["slots"], mode=case["mode"], schema=CANONICAL_REQUEST_SCHEMA, path="slots"))
    table = LegendTable({key: {"fill": fill} for key, fill in case.get("legend", [])})
    if case["mode"] == "circular":
        table = _sync_legend_table_for_circular_slots(
            table, slots=legend_slot_sources(parse_circular_track_slots(entries)), gc=gc_config, skew=skew_config,
            depth=None,
        )
    else:
        defaults = _generated_legend_rows(
            {}, [], pd.DataFrame(columns=["feature_type", "color"]), used_color_rules=None, default_used_features=set(),
            gc=gc_config, skew=skew_config, depth=None,
            show_gc=True, show_skew=True, show_depth=False,
        )
        for key, entry in defaults.items():
            table.add_default_track(key, {"fill": entry["fill"]})
        table = _sync_legend_table_for_linear_slots(
            table, slots=legend_slot_sources(parse_linear_track_slots(entries)), gc=gc_config, skew=skew_config,
            depth_tracks=None,
        )
    assert [[key, row["fill"]] for key, row in table.items()] == [[key, fill] for key, _, _, fill in case["rows"]]
    assert list(table.shared_keys) == case.get("shared", [])
