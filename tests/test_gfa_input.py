"""A GFA 1 assembly graph is a set of records: each segment (S line) with its sequence and no features (D-86, D-126).

One loader (``gbdraw.io.genome.load_gfa``, reading through ``gbdraw.io.gfa``)
owns the records, their topology and the refusals for the CLI, the Python API,
the typed request and the Web's Worker helper.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pytest

import gbdraw
from gbdraw.api import GfaInputSource
from gbdraw.api.input_summary import input_warnings
from gbdraw.api.record_planning import record_input_manifest_from_table, resolve_record_inputs
from gbdraw.api.requests import RecordCardinality, RecordInput
from gbdraw.circular import circular_main
from gbdraw.cli_utils.common import add_input_args, validate_input_args
from gbdraw.exceptions import ValidationError
from gbdraw.io.genome import load_gfa
from gbdraw.layout.record_transform import detected_topology
from gbdraw.linear import linear_main
from gbdraw.web_support.error_adapter import serialize_web_error

CIRCLE = "ACGTTGCA" * 20  # 160 bp; its last 10 bases repeat its first 10 in the graph below
GRAPH_INPUT = Path(__file__).resolve().parent / "test_inputs" / "gfa" / "HmmtDNA_lambda.gfa"


def _graph(path: Path, *, loop: str = "L\tcircle\t+\tcircle\t+\t10M", extra: str = "") -> Path:
    path.write_text(
        "H\tVN:Z:1.0\n"
        f"S\tcircle\t{CIRCLE}{CIRCLE[:10]}\n"
        "S\tleft\tGGGGCCCCAAAATTTT\n"
        "S\tright\tttttaaaaccccgggg\n"
        f"{loop}\n"
        "L\tleft\t+\tright\t+\t4M\n"
        f"{extra}"
    )
    return path


def _refusal(path: Path) -> ValidationError:
    with pytest.raises(ValidationError) as caught:
        load_gfa([str(path)])
    return caught.value


def test_each_segment_is_a_record_in_file_order_with_its_topology(tmp_path: Path) -> None:
    records = load_gfa([str(_graph(tmp_path / "graph.gfa"))])

    assert [(record.id, len(record), len(record.features)) for record in records] == [
        ("circle", 160, 0), ("left", 16, 0), ("right", 16, 0),
    ]
    # The self-loop's 10M overlap is dropped once: the circle has the genome's length.
    assert str(records[0].seq) == CIRCLE
    assert [detected_topology(record, "gfa") for record in records] == ["circular", "linear", "linear"]
    assert records[0].annotations["gbdraw_source_file"] == str(tmp_path / "graph.gfa")


@pytest.mark.parametrize(
    ("loop", "length", "topology"),
    [
        ("L\tcircle\t-\tcircle\t-\t10M", 160, "circular"),
        ("L\tcircle\t+\tcircle\t+\t0M", 170, "circular"),
        ("L\tcircle\t+\tcircle\t+\t*", 170, "circular"),
        ("L\tcircle\t+\tcircle\t+", 170, "circular"),
        # A link from the end to the end (an inversion) is no circle.
        ("L\tcircle\t+\tcircle\t-\t10M", 170, "linear"),
        # The same self-loop written both ways is one circle.
        ("L\tcircle\t+\tcircle\t+\t10M\nL\tcircle\t-\tcircle\t-\t10M", 160, "circular"),
    ],
)
def test_a_self_loop_drops_the_bases_its_overlap_repeats(tmp_path: Path, loop: str, length: int, topology: str) -> None:
    (circle, *_) = load_gfa([str(_graph(tmp_path / "graph.gfa", loop=loop))])

    assert (len(circle), circle.annotations["topology"]) == (length, topology)


@pytest.mark.parametrize(
    "loop",
    [
        "L\tcircle\t+\tcircle\t+\t5M1I4M",
        "L\tcircle\t+\tcircle\t+\t10M2D",
        "L\tcircle\t+\tcircle\t+\t170M",
        "L\tcircle\t+\tcircle\t+\t10M\nL\tcircle\t-\tcircle\t-\t12M",
    ],
)
def test_a_self_loop_overlap_gbdraw_cannot_drop_is_refused(tmp_path: Path, loop: str) -> None:
    error = _refusal(_graph(tmp_path / "graph.gfa", loop=loop))

    assert error.diagnostic == {"code": "INPUT_INVALID", "reason": "GFA_SELF_LOOP_OVERLAP", "recordName": "circle"}
    assert "circle" in str(error)


def test_a_segment_without_sequence_is_refused_with_its_name(tmp_path: Path) -> None:
    path = _graph(tmp_path / "graph.gfa", extra="S\tunknown\t*\tLN:i:500\n")

    error = _refusal(path)

    assert error.diagnostic == {"code": "INPUT_INVALID", "reason": "GFA_NO_SEQUENCE", "recordName": "unknown"}
    payload = serialize_web_error(error, operation="listSequenceRecords", stage="helper")
    assert payload["context"] == {"reason": "GFA_NO_SEQUENCE", "recordName": "unknown"}


@pytest.mark.parametrize(
    "text",
    [
        "H\tVN:Z:2.0\nS\ts1\t4\tACGT\n",
        "S\ts1\t4\tACGT\nE\te1\ts1+\ts1+\t0\t4$\t0\t4$\t4M\n",
    ],
)
def test_a_gfa_2_file_is_refused(tmp_path: Path, text: str) -> None:
    path = tmp_path / "graph.gfa"
    path.write_text(text)

    assert _refusal(path).diagnostic == {"code": "INPUT_UNREADABLE", "reason": "GFA_VERSION"}


@pytest.mark.parametrize(
    ("text", "diagnostic"),
    [
        ("H\tVN:Z:1.0\n", {"code": "NO_RECORDS", "reason": "GFA_NO_SEGMENTS"}),
        ("S\ts1\tACGT\nS\ts1\tGGCC\n", {"code": "INPUT_INVALID", "reason": "UNIQUE_IDS", "recordName": "s1"}),
        ("S\ts1\n", {"code": "INPUT_UNREADABLE", "reason": "GFA_LINE", "row": 1}),
        ("S\ts1\tACGT\nL\ts1\t+\n", {"code": "INPUT_UNREADABLE", "reason": "GFA_LINE", "row": 2}),
        (
            "S\ts1\tMKLVE\n",
            {"code": "INPUT_INVALID", "reason": "GFA_NOT_NUCLEOTIDE", "recordName": "s1", "recordPosition": 1,
             "codepoint": ord("L")},
        ),
    ],
)
def test_a_file_gbdraw_cannot_read_as_segments_is_refused(tmp_path: Path, text: str, diagnostic: dict) -> None:
    path = tmp_path / "graph.gfa"
    path.write_text(text)

    assert _refusal(path).diagnostic == diagnostic


def test_paths_and_walks_are_not_drawn_and_one_warning_names_them(tmp_path: Path) -> None:
    path = _graph(
        tmp_path / "graph.gfa",
        extra="P\tp1\tleft+,right+\t4M\nW\tsample\t1\tchr1\t0\t28\t>left>right\n",
    )
    inputs = (RecordInput(GfaInputSource(path), cardinality=RecordCardinality.ALL),)

    resolved = resolve_record_inputs(inputs, gff_candidate_features=None, gff_keep_all_features=True)

    assert [record.id for record in resolved.records] == ["circle", "left", "right"]
    assert [item.graph_lines_not_drawn for item in resolved.provenance] == [2, 2, 2]
    from gbdraw.api.input_summary import summarize_inputs

    summary = summarize_inputs(
        resolved.records, resolved.provenance, selected_features=[], feature_visibility_rules=None,
        specific_color_rules=None,
    )
    warnings = [item for item in input_warnings(summary) if item.code == "gfa_paths_not_drawn"]
    assert [item.message for item in warnings] == [
        "graph.gfa: its 2 path and walk lines (P, W) are not drawn; each segment (S) is drawn as a record."
    ]


def test_load_gfa_reads_gzip_and_applies_selector_and_reverse(tmp_path: Path) -> None:
    plain = _graph(tmp_path / "graph.gfa")
    packed = tmp_path / "graph.gfa.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()))

    (record,) = load_gfa([str(packed)], record_selectors=["left"], reverse_flags=[True])

    assert (record.id, str(record.seq)) == ("left", "AAAATTTTGGGGCCCC")


def test_read_gfa_is_the_beginner_reader(tmp_path: Path) -> None:
    assert "read_gfa" in gbdraw.__all__
    assert [record.id for record in gbdraw.read_gfa(_graph(tmp_path / "graph.gfa"))] == ["circle", "left", "right"]


def _parse_input_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    add_input_args(parser)
    parser.add_argument("--records_table", default=None)
    args = parser.parse_args(argv)
    validate_input_args(parser, args)
    return args


def test_cli_accepts_gfa_alone_and_refuses_it_with_other_inputs(capsys: pytest.CaptureFixture[str]) -> None:
    assert _parse_input_args(["--gfa", "a.gfa"]).gfa == ["a.gfa"]
    for other in (["--gbk", "a.gb"], ["--fasta", "a.fasta"]):
        with pytest.raises(SystemExit):
            _parse_input_args(["--gfa", "a.gfa", *other])
        assert "--gfa cannot be used with --gbk, --gff or --fasta" in capsys.readouterr().err


def test_cli_draws_a_gfa_in_both_modes(tmp_path: Path) -> None:
    circular_main(["--gfa", str(GRAPH_INPUT), "-o", str(tmp_path / "circular"), "-f", "svg"])
    linear_main(["--gfa", str(GRAPH_INPUT), "-o", str(tmp_path / "linear"), "-f", "svg"])

    assert "HmmtDNA" in (tmp_path / "circular_1.svg").read_text()
    linear = (tmp_path / "linear.svg").read_text()
    assert "lambda_left" in linear and "lambda_right" in linear


def test_records_table_row_with_a_gfa_path_is_a_gfa_source(tmp_path: Path) -> None:
    _graph(tmp_path / "graph.gfa")
    table = tmp_path / "records.tsv"
    table.write_text("gfa\trecord_id\ttopology\ngraph.gfa\tleft\tcircular\n")

    manifest = record_input_manifest_from_table(str(table))

    assert manifest.records[0].source == GfaInputSource(tmp_path / "graph.gfa")
    resolved = resolve_record_inputs(manifest.records, gff_candidate_features=None, gff_keep_all_features=True)
    assert [record.id for record in resolved.records] == ["left"]
    assert (resolved.provenance[0].source_kind, resolved.provenance[0].detected_topology) == ("gfa", "linear")
    assert resolved.provenance[0].display.is_circular is True


def test_records_table_refuses_mixing_gfa_rows_with_fasta_rows(tmp_path: Path) -> None:
    table = tmp_path / "records.tsv"
    table.write_text("gfa\tfasta\na.gfa\t\n\tb.fasta\n")

    with pytest.raises(ValidationError, match="cannot mix FASTA rows with GFA rows"):
        record_input_manifest_from_table(str(table))


def test_worker_helper_lists_a_gfa_through_the_gfa_loader(tmp_path: Path) -> None:
    source = (Path(__file__).resolve().parents[1] / "gbdraw/web/js/app/python-helpers.js").read_text()
    helpers: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0], helpers)
    list_records = helpers["list_sequence_records"]
    assert callable(list_records)

    listed = json.loads(list_records(str(_graph(tmp_path / "graph.gfa")), "gfa"))
    assert [(record["record_id"], record["record_length"], record["topology"]) for record in listed["records"]] == [
        ("circle", 160, "circular"), ("left", 16, "linear"), ("right", 16, "linear"),
    ]
    refused = json.loads(list_records(str(_graph(tmp_path / "bad.gfa", extra="S\tunknown\t*\n")), "gfa"))
    assert refused["error"]["context"] == {"reason": "GFA_NO_SEQUENCE", "recordName": "unknown"}


GFA_SESSION = Path(__file__).resolve().parent / "fixtures/sessions/cli-circular-gfa.v51.gbdraw"


def test_a_session_51_with_a_gfa_source_loads_and_replays(tmp_path: Path) -> None:
    from gbdraw.api import load_session_document
    from gbdraw.session import session_drawing_artifacts

    session = load_session_document(GFA_SESSION)
    drawing = session_drawing_artifacts(session, "circular").fields
    assert drawing["renderRequest"]["records"][0]["source"] == {"kind": "gfa", "gfaResourceId": "record-1-gfa"}

    circular_main(["--session", str(GFA_SESSION), "-o", str(tmp_path / "replay"), "-f", "svg"])
    assert "HmmtDNA" in (tmp_path / "replay.svg").read_text()
