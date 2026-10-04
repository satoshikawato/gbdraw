"""Feature override table: per-feature edits as TSV rows (design Q4, PR-Q4-3)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import (
    build_request_diagram,
    plan_request,
    read_request_feature_override_table,
    resolve_request,
)
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.overrides import FeatureOverride
from gbdraw.session_request_codec import CanonicalRequestEncodingError, encode_canonical_request

INPUTS = Path(__file__).parent / "test_inputs"
COLUMNS = ("record", "feature_selector", "feature_visibility", "label_visibility", "label_text")
GENBANK = """\
LOCUS       dup                      120 bp    DNA     linear   UNK 01-JAN-1980
DEFINITION  Identical CDS fixture.
FEATURES             Location/Qualifiers
     CDS             21..80
                     /locus_tag="first"
     CDS             21..80
                     /locus_tag="second"
     CDS             91..117
                     /locus_tag="third"
ORIGIN
        1 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg
       61 atgatgatga tgatgatgat gatgatgatg atgatgatga tgatgatgat gatgatgatg
//
"""


def _record() -> SeqRecord:
    record = SeqRecord(Seq("ATG" * 40), id="dup", annotations={"molecule_type": "DNA", "topology": "linear"})
    record.features = [
        SeqFeature(SimpleLocation(start, end, 1), type="CDS", qualifiers={"locus_tag": [tag]})
        for start, end, tag in ((20, 80, "first"), (20, 80, "second"), (90, 117, "third"))
    ]
    return record


def _request(mode: str = "linear", **options):
    request_type, options_type = (
        (LinearDiagramRequest, LinearDiagramOptions) if mode == "linear"
        else (CircularDiagramRequest, CircularDiagramOptions)
    )
    options.setdefault("selected_features_set", ("CDS",))
    return request_type(
        records=(RecordInput(InMemoryRecordSource(_record()), record_key="one"),),
        options=options_type(**options),
    )


def _ids(mode: str = "linear") -> tuple[str, str, str]:
    catalog = plan_request(_request(mode)).provenance[0].source_feature_catalog
    return tuple(entry.biological_feature_id for entry in catalog)


def _table(*rows: dict) -> DataFrame:
    return DataFrame([{column: row.get(column, "") for column in COLUMNS} for row in rows])


def _labels(svg: str) -> list[str]:
    return [
        re.sub(r"<[^>]+>", "", body)
        for body in re.findall(r'<text[^>]*data-label-feature-id="[^"]+"[^>]*>(.*?)</text>', svg, re.S)
    ]


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_table_rows_materialize_to_the_exact_rows_they_name(mode):
    first, second, third = _ids(mode)
    assert second.endswith("~1")
    table = _table(
        {"record": "#1", "feature_selector": f"hash={second}", "label_visibility": "on",
         "label_text": ' Second "copy" '},
        {"feature_selector": "locus_tag=first", "feature_visibility": "OFF"},
        {"record": "dup", "feature_selector": "locus_tag=third", "feature_visibility": "exclude_matching",
         "label_visibility": "off"},
    )
    exact = (
        FeatureOverride("one", first, feature_visibility="off"),
        FeatureOverride("one", second, label_visibility="on", label_text=' Second "copy" '),
        FeatureOverride("one", third, feature_visibility="exclude_matching", label_visibility="off"),
    )
    plan = plan_request(_request(mode, feature_override_table=table))
    assert plan.request.options.feature_overrides == tuple(sorted(exact, key=lambda row: row.biological_feature_id))
    assert plan.request.options.feature_override_table is None
    assert plan.request.options.feature_override_table_file is None
    assert plan.inputs.feature_identity_notices == ()
    from_table = build_request_diagram(_request(mode, feature_override_table=table)).drawing.tostring()
    from_rows = build_request_diagram(_request(mode, feature_overrides=exact)).drawing.tostring()
    assert from_table == from_rows
    assert _labels(from_table) == [' Second "copy" ']


def test_table_file_is_utf8_with_bom_and_keeps_label_text_verbatim(tmp_path):
    _first, second, _third = _ids()
    path = tmp_path / "edits.tsv"
    path.write_text(
        "﻿" + "\t".join(COLUMNS) + "\n"
        f"#1\thash={second}\t\ton\t\"  tab-free \"\"quoted\"\" text \"\n",
        encoding="utf-8",
    )
    plan = plan_request(_request(feature_override_table_file=str(path)))
    assert plan.request.options.feature_overrides == (
        FeatureOverride("one", second, label_visibility="on", label_text='  tab-free "quoted" text '),
    )
    assert plan.request.options.feature_override_table_file is None


def test_edit_columns_are_optional_and_blank_cells_set_nothing():
    first, _second, _third = _ids()
    table = DataFrame([{"feature_selector": f"hash={first}", "label_text": "Only text"}])
    plan = plan_request(_request(feature_override_table=table))
    assert plan.request.options.feature_overrides == (FeatureOverride("one", first, label_text="Only text"),)


@pytest.mark.parametrize(
    ("table", "pattern", "row"),
    [
        (DataFrame([{"feature_selector": "locus_tag=first", "label_visibility": "on", "color": "red"}]),
         "unknown columns", None),
        (DataFrame([{"record": "#1", "label_visibility": "on"}]), "requires feature_selector column", None),
        (_table({"feature_selector": "locus_tag=first", "feature_visibility": "hidden"}),
         "row 2: feature_visibility must be on, off, or exclude_matching", 2),
        (_table({"feature_selector": "locus_tag=first", "label_visibility": "on"},
                {"feature_selector": "locus_tag=third", "label_visibility": "shown"}),
         "row 3: label_visibility must be on or off", 3),
        (_table({"feature_selector": "locus_tag=first"}), "row 2: .*at least one edit", 2),
        (_table({"feature_selector": "locus_tag=first", "label_text": "   "}), "row 2: label_text", 2),
        (_table({"feature_selector": "locus_tag=missing", "label_visibility": "on"}), "row 2: .*matched 0", 2),
        (_table({"feature_selector": "type=CDS", "label_visibility": "on"}), "row 2: .*matched 3", 2),
        (_table({"feature_selector": "locus_tag=first", "label_visibility": "on"},
                {"feature_selector": "locus_tag=first", "label_text": "again"}),
         "row 3: Duplicate resolved feature identity", 3),
    ],
)
def test_invalid_tables_name_the_table_and_row(table, pattern, row):
    with pytest.raises(ValidationError, match=pattern) as caught:
        plan_request(_request(feature_override_table=table))
    assert "Feature override table" in str(caught.value)
    assert caught.value.diagnostic["code"] == "TABLE_INVALID"
    assert caught.value.diagnostic.get("row") == row


def test_rows_naming_no_record_or_feature_are_listed_instead_of_failing():
    # The Web's Load Feature Edits TSV (design Q4 6.4, Owner Q3 = A) reports a
    # row whose record or feature the current records lack; every other defect
    # still rejects the table, as the CLI does.
    first, second, _third = _ids()
    table = _table(
        {"record": "#1", "feature_selector": f"hash={second}", "label_text": "Kept"},
        {"record": "#2", "feature_selector": f"hash={first}", "feature_visibility": "off"},
        {"record": "other", "feature_selector": f"hash={first}", "feature_visibility": "off"},
        {"record": "#1", "feature_selector": "locus_tag=missing", "label_visibility": "on"},
    )
    unmatched: list[int] = []
    rows = read_request_feature_override_table(_request(), table, unmatched=unmatched)
    assert rows == (FeatureOverride("one", second, label_text="Kept"),)
    assert unmatched == [3, 4, 5]
    for defect, pattern in (
        ({"feature_selector": "type=CDS", "label_visibility": "on"}, "row 2: .*matched 3"),
        ({"record": "#0", "feature_selector": f"hash={first}", "label_visibility": "on"}, "row 2: Record index"),
        ({"feature_selector": f"hash={first}", "feature_visibility": "hidden"}, "row 2: feature_visibility"),
    ):
        with pytest.raises(ValidationError, match=pattern) as caught:
            read_request_feature_override_table(_request(), _table(defect), unmatched=[])
        assert caught.value.diagnostic == {"code": "TABLE_INVALID", "row": 2}


def test_unreadable_table_file_is_a_typed_error(tmp_path):
    with pytest.raises(ValidationError, match="Cannot read Feature override table") as caught:
        plan_request(_request(feature_override_table_file=str(tmp_path / "missing.tsv")))
    assert caught.value.diagnostic["code"] == "TABLE_INVALID"


def test_exact_rows_table_and_file_are_mutually_exclusive(tmp_path):
    first, _second, _third = _ids()
    row = FeatureOverride("one", first, label_visibility="on")
    table = _table({"feature_selector": "locus_tag=first", "label_visibility": "on"})
    for kwargs in (
        {"feature_overrides": (row,), "feature_override_table": table},
        {"feature_overrides": (row,), "feature_override_table_file": "edits.tsv"},
        {"feature_override_table": table, "feature_override_table_file": "edits.tsv"},
    ):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            LinearDiagramOptions(**kwargs)
    with pytest.raises(ValidationError, match="feature_override_table must be a DataFrame"):
        LinearDiagramOptions(feature_override_table=[{"feature_selector": "x"}])
    with pytest.raises(ValidationError, match="feature_override_table_file must identify a file"):
        LinearDiagramOptions(feature_override_table_file=" ")


def test_canonical_encoding_requires_a_materialized_table():
    table = _table({"feature_selector": "locus_tag=first", "feature_visibility": "off"})
    with pytest.raises(CanonicalRequestEncodingError, match="override tables must be materialized"):
        encode_canonical_request(_request(feature_override_table=table))
    first, _second, _third = _ids()
    payload = encode_canonical_request(resolve_request(_request(feature_override_table=table))).payload
    assert payload["diagramOptions"]["featureOverrides"] == [
        FeatureOverride("one", first, feature_visibility="off").to_mapping()
    ]


@pytest.mark.parametrize(
    ("request_type", "options_type"),
    [(LinearDiagramRequest, LinearDiagramOptions), (CircularDiagramRequest, CircularDiagramOptions)],
)
def test_gff_gene_shown_by_a_table_row_loads_like_an_exact_row(request_type, options_type):
    # NC_013668.gff3 links each CDS to its gene through Parent; the type filter drops genes.
    source = RecordInput(GffFastaInputSource(INPUTS / "NC_013668.gff3", INPUTS / "NC_013668.fasta"),
                         record_key="g")

    def svg(**options):
        return build_request_diagram(request_type(records=(source,), options=options_type(
            selected_features_set=("CDS", "tRNA", "rRNA", "repeat_region"), **options,
        ))).drawing.tostring()

    gene = next(
        entry for entry in plan_request(request_type(records=(source,), options=options_type()))
        .provenance[0].source_feature_catalog if entry.feature_type == "gene"
    )
    exact = svg(feature_overrides=(FeatureOverride("g", gene.biological_feature_id, feature_visibility="on"),))
    from_table = svg(feature_override_table=DataFrame([
        {"feature_selector": f"hash={gene.biological_feature_id}", "feature_visibility": "on"}
    ]))
    assert from_table == exact != svg()


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_cli_feature_override_table_applies_rows(tmp_path, mode):
    from gbdraw.circular import circular_main
    from gbdraw.linear import linear_main

    genbank = tmp_path / "dup.gb"
    genbank.write_text(GENBANK, encoding="utf-8")
    first, second, _third = _ids(mode)
    table = tmp_path / "edits.tsv"
    table.write_text(
        "\t".join(COLUMNS) + "\n"
        f"#1\thash={first}\toff\t\t\n"
        f"#1\thash={second}\t\ton\tSecond copy\n",
        encoding="utf-8",
    )
    main = linear_main if mode == "linear" else circular_main
    outputs = {}
    for name, extra in (("plain", []), ("edited", ["--feature_override_table", str(table)])):
        prefix = tmp_path / name
        main(["--gbk", str(genbank), "-o", str(prefix), "-f", "svg", "-k", "CDS", *extra])
        outputs[name] = prefix.with_suffix(".svg").read_text(encoding="utf-8")
    drawn = {name: len(re.findall(r'data-gbdraw-feature-id="', svg)) for name, svg in outputs.items()}
    assert drawn["edited"] == drawn["plain"] - 1
    assert "Second copy" in _labels(outputs["edited"])
    assert "Second copy" not in outputs["plain"]


def test_hash_rows_resolve_without_scanning_the_catalog(monkeypatch):
    # A 4,318-row hash= table on MG1655 took 145 s with one catalog scan per row.
    from gbdraw.features.source import SourceFeatureIdentity

    _first, second, _third = _ids()
    monkeypatch.setattr(SourceFeatureIdentity, "matches", lambda *_args, **_kwargs: pytest.fail("scan"))
    table = _table({"record": "#1", "feature_selector": f"hash={second}", "label_visibility": "on"})
    plan = plan_request(_request(feature_override_table=table))
    assert plan.request.options.feature_overrides == (FeatureOverride("one", second, label_visibility="on"),)
