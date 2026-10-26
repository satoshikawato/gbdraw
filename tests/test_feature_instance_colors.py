"""PD-OI-069 revision 2: a `hash` color rule that names a feature instance
(`<biological feature ID>@<record key>`) colors that feature only; a rule that
names the stable hash keeps coloring every feature with that hash."""
from __future__ import annotations

import copy
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

import json

import pytest

from gbdraw import cli
from gbdraw.api import RenderOutputRequest
from gbdraw.api.options import ColorOptions, LinearDiagramOptions
from gbdraw.api.request_render import build_request_diagram, plan_request
from gbdraw.api.requests import InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.features.colors import get_color, preprocess_color_tables
from gbdraw.features.ids import compute_source_feature_hash, feature_instance, feature_instances
from gbdraw.features.overrides import FeatureOverride
from gbdraw.exceptions import ValidationError
from gbdraw.session import build_session_document
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from gbdraw.web_support.rule_matching import evaluate_rules_json, exact_value_pattern

REPO = Path(__file__).resolve().parents[1]
LAMBDA = REPO / "tests" / "test_inputs" / "NC_001416.gb"
FILL = "#ff00aa"


def _color_map(value: str) -> dict:
    rules = DataFrame([dict(feature_type="CDS", qualifier_key="hash", value=value, color=FILL, caption="mine")])
    return preprocess_color_tables(rules, DataFrame(columns=["feature_type", "color"]))[0]


def _keyed(record: SeqRecord, key: str) -> SeqRecord:
    copy = record[:]
    copy.id, copy.annotations = record.id, {**record.annotations, "gbdraw_record_key": key}
    return copy


def _fills(records: list[SeqRecord], value: str) -> list[list[bool]]:
    color_map = _color_map(value)
    fills = []
    for record in records:
        instance_of = feature_instances(record)
        fills.append([
            get_color(feature, color_map, {}, record_id=record.id, instance=instance_of(feature)) == FILL
            for feature in record.features if feature.type == "CDS"
        ])
    return fills


def test_an_instance_rule_colors_its_record_instance_and_a_hash_rule_both():
    record = SeqRecord(Seq("A" * 300), id="Chr1",
                       features=[SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS")])
    first, second = _keyed(record, "record-1"), _keyed(record, "record-2")
    feature_hash = compute_source_feature_hash(record.features[0], record)
    assert _fills([first, second], exact_value_pattern(feature_instance(feature_hash, "record-2"))) == [[False], [True]]
    assert _fills([first, second], feature_hash) == [[True], [True]]


def test_an_instance_rule_names_one_of_identical_features_of_a_record():
    twin = SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS")
    record = _keyed(SeqRecord(Seq("A" * 300), id="r", features=[twin, SeqFeature(twin.location, type="CDS")]), "k")
    feature_hash = compute_source_feature_hash(twin, record)
    instances = [feature_instances(record)(feature) for feature in record.features]
    assert instances == [f"{feature_hash}~0@k", f"{feature_hash}~1@k"]
    assert _fills([record], exact_value_pattern(instances[1])) == [[False, True]]
    # Anchored: `f@k` names neither identical feature, and `~1` is not `~10`.
    assert _fills([record], exact_value_pattern(f"{feature_hash}@k")) == [[False, False]]


def test_a_typed_instance_value_names_the_whole_instance_and_a_hash_keeps_its_search():
    # Eleven identical features in two records: `~1` is not `~10`, and `@k-1` is not `@k-10`.
    twin = SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS")
    record = SeqRecord(Seq("A" * 300), id="r", features=[SeqFeature(twin.location, type="CDS") for _ in range(11)])
    records = [_keyed(record, "k-1"), _keyed(record, "k-10")]
    feature_hash = compute_source_feature_hash(twin, record)
    one = [[index == 1 for index in range(11)], [False] * 11]
    assert _fills(records, f"{feature_hash}~1@k-1") == one
    assert _fills(records, exact_value_pattern(f"{feature_hash}~1@k-1")) == one
    assert _fills(records, f"{feature_hash}~1") == [[False] * 11, [False] * 11]
    # A part of the stable hash still names every feature that shares it.
    assert _fills(records, feature_hash[:5]) == [[True] * 11, [True] * 11]


def test_the_instances_of_a_record_follow_its_record_key_and_feature_list():
    feature = SeqFeature(SimpleLocation(10, 70, strand=1), type="CDS")
    record = _keyed(SeqRecord(Seq("A" * 300), id="r", features=[feature]), "k1")
    feature_hash = compute_source_feature_hash(feature, record)
    assert feature_instances(record)(record.features[0]) == f"{feature_hash}@k1"
    record.annotations["gbdraw_record_key"] = "k2"
    assert feature_instances(record)(record.features[0]) == f"{feature_hash}@k2"
    moved = SeqFeature(SimpleLocation(100, 160, strand=1), type="CDS")
    record.features = [moved]
    assert feature_instances(record)(moved) == f"{compute_source_feature_hash(moved, record)}@k2"
    # A feature outside the record's features has no instance.
    assert feature_instances(record)(feature) is None
    # A deep copy has new feature objects; it computes its own instances.
    copied = copy.deepcopy(record)
    assert feature_instances(copied)(copied.features[0]) == f"{compute_source_feature_hash(moved, record)}@k2"


def _run_cli(*args: str) -> None:
    argv = sys.argv
    sys.argv = ["gbdraw", *args]
    try:
        cli.main()
    finally:
        sys.argv = argv


def _filled_features(tmp_path: Path, value: str) -> tuple[list[str], int]:
    tmp_path = tmp_path / str(len(list(tmp_path.iterdir())))
    tmp_path.mkdir()
    table = tmp_path / "colors.tsv"
    table.write_text(f"CDS\thash\t{value}\t{FILL}\tmine\n", encoding="utf-8")
    _run_cli("linear", "--gbk", str(LAMBDA), str(LAMBDA), "-t", str(table),
             "-o", str(tmp_path / "out"), "-f", "svg")
    svg = (tmp_path / "out.svg").read_text(encoding="utf-8")
    filled = re.findall(r'<path [^>]*data-gbdraw-feature-id="([^"]+)"[^>]*fill="' + FILL, svg)
    # Each legend layout (horizontal, vertical) lists a row once.
    return filled, svg.count('data-legend-key="mine"') // len(re.findall(r'id="legend_', svg))


def test_the_cli_colors_the_feature_of_the_nth_record_input_by_its_record_key(tmp_path: Path):
    record = next(SeqIO.parse(LAMBDA, "genbank"))
    feature_hash = compute_source_feature_hash(next(f for f in record.features if f.type == "CDS"), record)
    # The CLI record key of the n-th record input is `record-<n>`.
    filled, legend_rows = _filled_features(tmp_path, exact_value_pattern(feature_instance(feature_hash, "record-2")))
    assert filled == [f"{feature_hash}_record_2"] and legend_rows == 1
    filled, legend_rows = _filled_features(tmp_path, feature_hash)
    assert sorted(filled) == [f"{feature_hash}_record_1", f"{feature_hash}_record_2"] and legend_rows == 1


def _source_twins_request(**options) -> LinearDiagramRequest:
    record = SeqRecord(Seq("ATG" * 40), id="dup", annotations={"molecule_type": "DNA", "topology": "linear"})
    record.features = [
        SeqFeature(SimpleLocation(0, 120, 1), type="source", qualifiers={"note": [note]}) for note in ("a", "b")
    ]
    return LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(record), record_key="one"),),
        options=LinearDiagramOptions(**options),
        output=RenderOutputRequest(output_prefix="twins", formats=("svg",)),
    )


def _web_render(tmp_path: Path, request: LinearDiagramRequest) -> tuple[list[dict], str]:
    document = build_session_document(request).to_dict()
    reply = render_embedded_canonical_web_request(
        document["renderRequest"], resources=document["resources"],
        workspace=tmp_path / str(len(list(tmp_path.iterdir()))),
    )
    assert reply["featureCatalog"] is not None
    return reply["featureCatalog"]["items"][0]["biologicalFeatures"], reply["results"][0]["content"]


def test_a_shown_source_feature_with_an_undrawn_identical_twin_colors_alike_live_and_after_generate(tmp_path: Path):
    # Two identical source features; the first is shown by its own Feature
    # visibility, the second stays undrawn. The catalog names them as the
    # source feature catalog does, so the instance the live preview colors is
    # the one Generate colors.
    first, second = (
        entry.biological_feature_id
        for entry in plan_request(_source_twins_request()).provenance[0].source_feature_catalog
    )
    shown = {"feature_overrides": (FeatureOverride("one", first, feature_visibility="on"),)}
    rows, _ = _web_render(tmp_path, _source_twins_request(**shown))
    assert [(row["biologicalFeatureId"], row["recordKey"]) for row in rows] == [(first, "one"), (second, "one")]

    value = exact_value_pattern(feature_instance(first, "one"))
    rule = {"feat": "source", "qual": "hash", "val": value, "color": FILL, "cap": "mine"}
    payloads = [
        {"type": "source", "record": row["record_id"], "qualifiers": row["qualifiers"], "label": "",
         "selector": {"hash": row["stableFeatureId"], "record_key": row["recordKey"],
                      "feature_id": row["biologicalFeatureId"], "location": None, "record_location": None}}
        for row in rows
    ]
    assert json.loads(evaluate_rules_json(json.dumps(payloads), json.dumps([rule])))["winners"] == [0, -1]

    table = DataFrame([dict(feature_type="source", qualifier_key="hash", value=value, color=FILL, caption="mine")])
    _, svg = _web_render(tmp_path, _source_twins_request(**shown, colors=ColorOptions(color_table=table)))
    filled = re.findall(r'<path [^>]*data-gbdraw-source-feature-index="(\d+)"[^>]*fill="' + FILL, svg)
    assert filled == ["0"]


def test_a_fill_and_a_label_edit_of_one_feature_both_name_its_record_key():
    # OV-1949 (open): a "This feature only" fill names the record key as a
    # per-feature label edit does, so a presentation change that gives the
    # record another key drops both. The Web sends no label row of a record key
    # the request lacks (Python refuses one); the fill row stays and colors nothing.
    record = SeqRecord(Seq("ATG" * 40), id="r", annotations={"molecule_type": "DNA", "topology": "linear"})
    record.features = [SeqFeature(SimpleLocation(20, 80, 1), type="CDS")]

    def request(record_key: str, **options) -> LinearDiagramRequest:
        return LinearDiagramRequest(
            records=(RecordInput(InMemoryRecordSource(copy.deepcopy(record)), record_key=record_key),),
            options=LinearDiagramOptions(selected_features_set=("CDS",), **options),
        )

    feature_id = plan_request(request("k1")).provenance[0].source_feature_catalog[0].biological_feature_id
    value = exact_value_pattern(feature_instance(feature_id, "k1"))
    fill = ColorOptions(color_table=DataFrame([dict(feature_type="CDS", qualifier_key="hash", value=value, color=FILL,
                                                    caption="mine")]))
    label = (FeatureOverride("k1", feature_id, label_visibility="on", label_text="Mine"),)
    svg = build_request_diagram(request("k1", colors=fill, feature_overrides=label)).drawing.tostring()
    assert FILL in svg and "Mine" in svg
    with pytest.raises(ValidationError, match="record key 'k1'"):
        build_request_diagram(request("k2", colors=fill, feature_overrides=label))
    assert FILL not in build_request_diagram(request("k2", colors=fill)).drawing.tostring()
