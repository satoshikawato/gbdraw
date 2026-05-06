from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity_units import build_collinearity_unit_index
from gbdraw.analysis.protein_colinearity import extract_cds_proteins
from gbdraw.exceptions import ValidationError


def _record(record_id: str, features: list[SeqFeature]) -> SeqRecord:
    record = SeqRecord(Seq("ATGAAATAG" * 20), id=record_id)
    record.features = features
    return record


def _cds(start: int, end: int, qualifiers: dict[str, list[str]]) -> SeqFeature:
    merged = {"translation": ["MK"]}
    merged.update(qualifiers)
    return SeqFeature(FeatureLocation(start, end, strand=1), type="CDS", qualifiers=merged)


@pytest.mark.linear
def test_auto_unit_mode_collapses_strong_locus_ids() -> None:
    record = _record(
        "record_a",
        [
            _cds(0, 9, {"locus_tag": ["locus_a"], "protein_id": ["p1"]}),
            _cds(12, 24, {"locus_tag": ["locus_a"], "protein_id": ["p2"]}),
            _cds(30, 39, {"locus_tag": ["locus_b"], "protein_id": ["p3"]}),
        ],
    )
    extraction = extract_cds_proteins([record])

    unit_index = build_collinearity_unit_index(extraction, records=[record], mode="auto")

    units = unit_index.units_by_record[0]
    assert [unit.unit_kind for unit in units] == ["locus", "locus"]
    assert units[0].locus_id == "locus_a"
    assert units[0].cds_members == ("p1", "p2")
    assert units[0].start == 0
    assert units[0].end == 24


@pytest.mark.linear
def test_gene_labels_do_not_drive_unit_collapse() -> None:
    record = _record(
        "record_a",
        [
            _cds(0, 9, {"gene": ["abc"], "protein_id": ["p1"]}),
            _cds(12, 24, {"gene": ["abc"], "protein_id": ["p2"]}),
        ],
    )
    extraction = extract_cds_proteins([record])

    unit_index = build_collinearity_unit_index(extraction, records=[record], mode="auto")

    units = unit_index.units_by_record[0]
    assert [unit.unit_kind for unit in units] == ["cds", "cds"]
    assert units[0].unit_id != units[1].unit_id
    assert "abc" in unit_index.ambiguous_aliases_by_record[0]


@pytest.mark.linear
def test_locus_unit_mode_rejects_product_only_cds() -> None:
    record = _record(
        "record_a",
        [_cds(0, 9, {"product": ["hypothetical protein"]})],
    )
    extraction = extract_cds_proteins([record])

    with pytest.raises(ValidationError, match="requires stable locus identifiers"):
        build_collinearity_unit_index(extraction, records=[record], mode="locus")


@pytest.mark.linear
def test_cds_unit_mode_never_collapses_same_locus() -> None:
    record = _record(
        "record_a",
        [
            _cds(0, 9, {"locus_tag": ["locus_a"], "protein_id": ["p1"]}),
            _cds(12, 24, {"locus_tag": ["locus_a"], "protein_id": ["p2"]}),
        ],
    )
    extraction = extract_cds_proteins([record])

    unit_index = build_collinearity_unit_index(extraction, records=[record], mode="cds")

    units = unit_index.units_by_record[0]
    assert [unit.unit_kind for unit in units] == ["cds", "cds"]
    assert [unit.representative_protein_id for unit in units] == ["p1", "p2"]
