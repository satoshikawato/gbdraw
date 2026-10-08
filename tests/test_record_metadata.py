from __future__ import annotations

import json
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import (
    RecordSourceMetadata,
    format_inferred_definition,
    format_replicon_label,
    infer_record_source_metadata,
)
from gbdraw.exceptions import GbdrawError
from gbdraw.io.genome import load_gbks, load_gff_fasta

# The Web upload fast path reimplements this inference in
# gbdraw/web/js/services/record-discovery.js. Both suites assert the same table, so
# the two implementations cannot drift apart unnoticed.
_CASES = json.loads(
    (Path(__file__).parent / "fixtures" / "record_metadata_inference_cases.json").read_text(
        encoding="utf-8"
    )
)


@pytest.mark.parametrize("case", _CASES["definition"], ids=lambda case: case["name"])
def test_format_inferred_definition_matches_shared_cases(case: dict[str, str]) -> None:
    meta = RecordSourceMetadata(
        organism=case["organism"],
        strain=case["strain"],
        replicon=None,
        organelle=None,
    )
    assert format_inferred_definition(meta) == case["expected"]


@pytest.mark.parametrize("qualifiers,expected", [
    ({"chromosome": ["1"], "plasmid": ["p1"], "organelle": ["chloroplast"]}, "Chromosome 1"),
    ({"plasmid": ["p1"], "organelle": ["mitochondrion"]}, "p1"),
    ({"organelle": ["mitochondrion"]}, "Mitochondrion"),
    ({"organelle": ["plastid:chloroplast"]}, "Plastid:chloroplast"),
    ({}, ""),
])
def test_replicon_label_uses_source_qualifiers(qualifiers: dict, expected: str) -> None:
    record = SeqRecord(Seq("ATGC"), id="test", description="plasmid pDescription, complete genome")
    record.features.append(SeqFeature(SimpleLocation(0, 4), type="source", qualifiers=qualifiers))
    metadata = infer_record_source_metadata(record)
    assert format_replicon_label(metadata) == expected
    # Circular consumes these distinct fields; formatting a Linear label cannot mutate them.
    assert metadata.organelle == qualifiers.get("organelle", [None])[0]
    if not qualifiers.get("chromosome") and not qualifiers.get("plasmid"):
        assert metadata.replicon is None


def test_infer_record_source_metadata():
    record = SeqRecord(Seq("ATGC"), id="test", annotations={"organism": "Bacillus subtilis"})
    feature = SeqFeature(
        SimpleLocation(0, 4),
        type="source",
        qualifiers={"isolate": ["168"], "chromosome": ["1"]}
    )
    record.features.append(feature)
    meta = infer_record_source_metadata(record)
    assert meta.organism == "Bacillus subtilis"
    assert meta.strain == "168"
    assert meta.replicon == "Chromosome 1"
    assert meta.organelle is None



# G-D (Web GUI audit 2026-09-30): record discovery shares small synthetic inputs
# with tests/web/record-metadata-inference.test.mjs. The gbdraw.io.genome loader
# is the oracle; the Worker helper and the JS fast path must agree with it.
_REPO_ROOT = Path(__file__).resolve().parents[1]
_PYTHON_HELPERS_PATH = _REPO_ROOT / "gbdraw" / "web" / "js" / "app" / "python-helpers.js"


def _discovery_case_id(case: dict) -> str:
    return case["name"]


def _write_discovery_file(directory: Path, spec: dict) -> str:
    data = ("\n".join(spec["lines"]) + "\n").encode("utf-8")
    if spec.get("bom"):
        data = b"\xef\xbb\xbf" + data
    path = directory / spec["name"]
    path.write_bytes(data)
    return str(path)


def _write_discovery_inputs(directory: Path, case: dict) -> list[str]:
    if case["format"] == "genbank":
        return [_write_discovery_file(directory, case["files"]["source"])]
    return [
        _write_discovery_file(directory, case["files"]["gff"]),
        _write_discovery_file(directory, case["files"]["fasta"]),
    ]


def _loader_discovery(directory: Path, case: dict) -> dict:
    paths = _write_discovery_inputs(directory, case)
    try:
        records = load_gbks(paths) if case["format"] == "genbank" else load_gff_fasta([paths[0]], [paths[1]])
    except GbdrawError as error:
        return {"error": type(error).__name__}
    projected = []
    for record in records:
        entry = {"recordId": record.id, "recordLength": len(record.seq)}
        if case["format"] == "genbank":
            metadata = infer_record_source_metadata(record)
            entry["organism"] = metadata.organism or ""
            entry["inferredDefinition"] = format_inferred_definition(metadata)
        projected.append(entry)
    return {"records": projected}


@pytest.fixture(scope="module")
def _python_helpers_namespace() -> dict[str, object]:
    source = _PYTHON_HELPERS_PATH.read_text(encoding="utf-8")
    helper_source = source.split("export const PYTHON_HELPERS = `", 1)[1].rsplit("\n`;", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)
    return namespace


@pytest.mark.parametrize("case", _CASES["discovery"], ids=_discovery_case_id)
def test_loader_is_the_record_discovery_oracle(case: dict, tmp_path: Path) -> None:
    assert _loader_discovery(tmp_path, case) == case["expected"]


@pytest.mark.parametrize("case", _CASES["discovery"], ids=_discovery_case_id)
def test_worker_record_helper_matches_the_loader(
    case: dict, tmp_path: Path, _python_helpers_namespace: dict[str, object]
) -> None:
    paths = _write_discovery_inputs(tmp_path, case)
    if case["format"] == "genbank":
        payload = json.loads(_python_helpers_namespace["list_sequence_records"](paths[0], "genbank"))
    else:
        payload = json.loads(_python_helpers_namespace["list_gff_fasta_records"](paths[0], paths[1]))
    if case["expected"].get("error"):
        assert "error" in payload
        return
    assert "error" not in payload, payload
    observed = []
    for entry in payload["records"]:
        projected = {"recordId": entry["record_id"], "recordLength": entry["record_length"]}
        if case["format"] == "genbank":
            projected["organism"] = entry["organism"]
            projected["inferredDefinition"] = entry["inferred_definition"]
        observed.append(projected)
    assert observed == case["expected"]["records"]
