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

from gbdraw import cli
from gbdraw.features.colors import get_color, preprocess_color_tables
from gbdraw.features.ids import compute_source_feature_hash, feature_instance, feature_instances
from gbdraw.web_support.rule_matching import exact_value_pattern

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
