"""A LOSATP comparison draws records that have no CDS protein (OV-3493).

A record without a CDS protein, such as a short contig of a draft assembly,
adds nothing to a protein search: it draws without protein links, and a
database it belongs to holds the other records' proteins, so no E-value
changes. Only a comparison whose displayed pairs all lack proteins on a side
has nothing to compare, and its refusal names a record without proteins.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest

import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    build_pairwise_protein_blastp_comparisons,
    extract_cds_proteins,
    extract_protein_identity_manifest,
    promote_legacy_protein_raw_cache_entries,
    proteins_to_fasta,
)
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.web_support.error_adapter import serialize_web_error
from tests.utils.session_drawing_view import drawing_view
from tests.utils.session_files import read_session_payload

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUTS = REPO_ROOT / "tests" / "test_inputs"

# A LOSAT that refuses what real LOSAT refuses: a search with no sequence on a
# side (an empty database is a LOSAT error, exit 3).
_FAKE = """#!{python}
import json, sys
args = sys.argv[1:]
if args in (["--version"], ["-version"]):
    print("losat 0.1.0")
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    sys.exit(0)
def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()
query = read(args[args.index("-query") + 1])
subject = read(args[args.index("-subject") + 1])
def ids(text):
    return [line[1:].split()[0] for line in text.splitlines() if line.startswith(">")]
with open({log!r}, "a", encoding="utf-8") as handle:
    handle.write(json.dumps({{"query": ids(query), "subject": ids(subject)}}) + "\\n")
if not ids(query) or not ids(subject):
    sys.stderr.write("BLAST engine error: Empty CBlastQueryVector\\n")
    sys.exit(3)
for query_id in ids(query):
    for subject_id in ids(subject):
        print("\\t".join([query_id, subject_id, "80.000", "10", "2", "0", "1", "10",
                         "1", "10", "1e-20", "90"]))
"""

MODES = ("pairwise", "similarity_groups", "collinear")


@pytest.fixture(autouse=True)
def _fresh_probe_caches(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(runtime_module, "_RUNTIME_VERSIONS", {})
    monkeypatch.setattr(runtime_module, "_TASK_VALUES", {}, raising=False)


def _contig(record_id: str) -> SeqRecord:
    """A contig with a sequence and no CDS."""

    record = SeqRecord(Seq("ACGTTGCA" * 120), id=record_id, name=record_id, description=f"{record_id} contig")
    record.annotations["molecule_type"] = "DNA"
    record.features = [SeqFeature(FeatureLocation(0, len(record.seq)), type="source")]
    return record


def _inputs(tmp_path: Path) -> Path:
    """A draft assembly whose middle contig has no CDS, and a single-record file."""

    draft = [
        *SeqIO.parse(INPUTS / "BGC0000708.gbk", "genbank"),
        _contig("NOCDS1"),
        *SeqIO.parse(INPUTS / "BGC0000709.gbk", "genbank"),
    ]
    SeqIO.write(draft, tmp_path / "draft.gbk", "genbank")
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_id\n"
        + "".join(f"{tmp_path / 'draft.gbk'}\t{record.id}\n" for record in draft)
        + f"{INPUTS / 'BGC0000711.gbk'}\tBGC0000711\n",
        encoding="utf-8",
    )
    return table


def _run(
    tmp_path: Path, name: str, argv: list[str], mode: str, *, save_session: bool = True
) -> tuple[Path, list[dict]]:
    log = tmp_path / name / "searches.jsonl"
    losat = tmp_path / name / "bin" / "losat"
    losat.parent.mkdir(parents=True, exist_ok=True)
    losat.write_text(_FAKE.format(python=sys.executable, log=str(log)), encoding="utf-8")
    losat.chmod(0o755)
    prefix = tmp_path / name / "out"
    linear_main([
        *argv, "--losat", "losatp", "--losatp_mode", mode, "--losat_bin", str(losat),
        "--losat_threads", "1", "-o", str(prefix), "-f", "svg", *(["--save_session"] if save_session else []),
    ])
    searches = [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()] if log.exists() else []
    return prefix, searches


@pytest.mark.linear
@pytest.mark.parametrize("mode", MODES)
def test_a_record_without_cds_draws_without_protein_links(tmp_path: Path, mode: str) -> None:
    table = _inputs(tmp_path)
    prefix, searches = _run(tmp_path, mode, ["--records_table", str(table)], mode)

    assert prefix.with_suffix(".svg").is_file()
    # No search reads a side without proteins; the record adds nothing to a database it is in.
    assert searches and all(search["query"] and search["subject"] for search in searches)
    session = drawing_view(read_session_payload(prefix.with_suffix(".gbdraw")))
    assert [record["selector"] for record in session["renderRequest"]["records"]][:3] == [
        {"kind": "recordId", "value": value} for value in ("BGC0000708", "NOCDS1", "BGC0000709")
    ]
    entries = session["losatCache"]["entries"]
    without = [entry for entry in entries if "record-2" in (entry["queryRecordInstanceKey"], entry["subjectRecordInstanceKey"])]
    others = [entry for entry in entries if entry not in without]
    # The pairs of the record without proteins found nothing; the others found their rows.
    assert without and all(entry["text"] == "" for entry in without)
    assert others and any(entry["text"] for entry in others)


def _proteinless_pair_records() -> list[SeqRecord]:
    return [next(SeqIO.parse(INPUTS / "BGC0000708.gbk", "genbank")), _contig("NOCDS1")]


@pytest.mark.linear
def test_a_comparison_whose_pairs_all_lack_proteins_names_a_record_without_them() -> None:
    records = _proteinless_pair_records()
    with pytest.raises(ValidationError) as caught:
        build_pairwise_protein_blastp_comparisons(records, protein_extraction=extract_cds_proteins(records))
    assert "NOCDS1" in str(caught.value)
    web = serialize_web_error(caught.value, operation="generate", stage="render")
    assert (web["code"], web["context"]["reason"], web["context"]["recordName"]) == (
        "COMPARISON_INPUT", "NO_PROTEINS", "NOCDS1",
    )


@pytest.mark.linear
@pytest.mark.parametrize("mode", MODES)
def test_the_command_line_refuses_only_a_comparison_with_nothing_to_compare(tmp_path: Path, mode: str) -> None:
    SeqIO.write([_contig("NOCDS1"), _contig("NOCDS2")], tmp_path / "contigs.gbk", "genbank")
    with pytest.raises(ValidationError, match="NOCDS1"):
        _run(tmp_path, mode, ["--gbk", str(tmp_path / "contigs.gbk")], mode)


@pytest.mark.linear
def test_the_web_helper_extracts_a_record_without_cds_as_no_proteins(tmp_path: Path) -> None:
    SeqIO.write([_contig("NOCDS1")], tmp_path / "contig.gbk", "genbank")
    extract = _helpers()["extract_cds_protein_fasta"]
    assert callable(extract)
    reply = json.loads(str(extract(str(tmp_path / "contig.gbk"), "genbank", None, None, None, "0", 0, "record-1")))
    assert "error" not in reply, reply
    assert (reply["fasta"], reply["protein_count"], reply["protein_map"], reply["record_id"]) == ("", 0, {}, "NOCDS1")
    assert reply["protein_set_hash"].startswith("sha256:")


def _helpers() -> dict[str, object]:
    """The Worker's Python helpers, run as the Worker runs them (python-helpers-source.test.mjs pins the text)."""

    helpers_js = (REPO_ROOT / "gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(helpers_js.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    return namespace


@pytest.mark.linear
@pytest.mark.parametrize("mode", MODES)
def test_a_figure_without_a_displayed_pair_draws(tmp_path: Path, mode: str) -> None:
    # Both records in one row: no pair is displayed, so there is nothing to refuse. (With --save_session, Pairwise
    # stops on dev too: OV-3511.)
    prefix, _searches = _run(tmp_path, mode, [
        "--gbk", str(INPUTS / "BGC0000708.gbk"), str(INPUTS / "BGC0000709.gbk"),
        "--multi_record_position", "#1@1", "--multi_record_position", "#2@1",
    ], mode, save_session=False)
    assert prefix.with_suffix(".svg").is_file()


@pytest.mark.linear
@pytest.mark.parametrize("mode", ["pairwise", "orthogroup", "collinear"])
def test_the_web_converter_refuses_a_comparison_with_nothing_to_compare(tmp_path: Path, mode: str) -> None:
    pairs = tmp_path / "pairs.json"
    raw = tmp_path / "raw.tsv"
    raw.write_text("", encoding="utf-8")
    pairs.write_text(json.dumps({
        "records": [
            {"recordIndex": index, "recordId": record_id, "proteinMap": {}, "viewTransform": {"length": 960}}
            for index, record_id in enumerate(("NOCDS1", "NOCDS2"))
        ],
        "pairs": [{
            "queryIndex": 0, "subjectIndex": 1, "pairIndex": 0, "displayPair": True, "cacheKey": "raw-0-1",
            "rawTsvOffset": 0, "rawTsvBytes": 0,
        }],
    }), encoding="utf-8")
    convert = _helpers()["convert_losatp_blastp_pairs_to_genomic_payload"]
    assert callable(convert)
    reply = json.loads(str(convert(str(pairs), str(raw), mode=mode)))
    error = reply.get("error")
    assert isinstance(error, dict), reply
    assert (error["code"], error["context"]["reason"], error["context"]["recordName"]) == (
        "COMPARISON_INPUT", "NO_PROTEINS", "NOCDS1",
    )


def _legacy_promotion_inputs() -> tuple[list[SeqRecord], object]:
    records = _proteinless_pair_records()
    return records, extract_protein_identity_manifest(records)


@pytest.mark.linear
def test_legacy_promotion_skips_a_record_without_proteins() -> None:
    records, extraction = _legacy_promotion_inputs()
    query = extraction.proteins_by_record[0]
    scan = promote_legacy_protein_raw_cache_entries(
        [],
        query_proteins=query,
        subject_proteins=[],
        query_fasta=proteins_to_fasta(query),
        subject_fasta="",
        identity_manifest=extraction.identity_manifest,
        expected_args=[],
    )
    assert (scan.promotion, scan.rejections) == (None, ())


@pytest.mark.linear
def test_the_web_legacy_promotion_skips_records_without_proteins() -> None:
    _records, extraction = _legacy_promotion_inputs()
    promote = _helpers()["promote_legacy_losatp_cache_candidates"]
    assert callable(promote)
    reply = json.loads(str(promote(
        "[]", "", "", "{}", "{}", json.dumps(extraction.identity_manifest), json.dumps({"args": []}),
    )))
    assert reply == {"status": "no-match", "rejections": []}


@pytest.mark.linear
def test_a_raw_cache_without_an_identity_manifest_is_a_validation_error() -> None:
    records = [next(SeqIO.parse(INPUTS / name, "genbank")) for name in ("BGC0000708.gbk", "BGC0000709.gbk")]
    with pytest.raises(ValidationError, match="identity manifest"):
        build_pairwise_protein_blastp_comparisons(records, losatp_cache=LosatpCacheManager())
