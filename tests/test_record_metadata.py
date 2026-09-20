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

# The Web upload fast path reimplements this inference in
# gbdraw/web/js/app/record-discovery.js. Both suites assert the same table, so
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

