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
    format_inferred_subtitle,
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


@pytest.mark.parametrize("case", _CASES["subtitle"], ids=lambda case: case["name"])
def test_format_inferred_subtitle_matches_shared_cases(case: dict[str, str]) -> None:
    meta = RecordSourceMetadata(
        organism=case["organism"],
        strain="",
        replicon=case["replicon"] or None,
        organelle=case["organelle"] or None,
    )
    assert format_inferred_subtitle(meta, case["description"]) == case["expected"]


@pytest.mark.parametrize(
    "case",
    [case for case in _CASES["subtitle"] if case["chromosome"] or case["plasmid"]],
    ids=lambda case: case["name"],
)
def test_shared_cases_agree_with_inferred_replicon(case: dict[str, str]) -> None:
    """The fixture's replicon must be what the qualifiers actually produce."""
    qualifiers = {
        key: [case[key]]
        for key in ("chromosome", "plasmid")
        if case[key]
    }
    record = SeqRecord(Seq("ATGC"), id="test")
    record.features.append(SeqFeature(SimpleLocation(0, 4), type="source", qualifiers=qualifiers))
    assert infer_record_source_metadata(record).replicon == case["replicon"]


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

