"""One GFF3 parse serves every type filter, and filtering it in memory equals the type-filtered load."""

from __future__ import annotations

from pathlib import Path

import pytest
from BCBio import GFF

from gbdraw.api.record_planning import resolve_record_inputs
from gbdraw.api.requests import GffFastaInputSource, RecordCardinality, RecordInput
from gbdraw.core.record_metadata import _feature_source_index_map, _iter_source_features, _source_feature_index
from gbdraw.io.genome import (
    _normalize_gff3_multipart_features,
    filter_features_by_type,
    load_gff_fasta,
    parse_gff_fasta,
)

INPUTS = Path(__file__).parent / "test_inputs"
NESTED_GFF3 = """\
##gff-version 3
rec1\ttest\tgene\t1\t60\t.\t+\t.\tID=g1;Name=alpha
rec1\ttest\tmRNA\t1\t60\t.\t+\t.\tID=m1;Parent=g1
rec1\ttest\texon\t1\t20\t.\t+\t.\tID=e1;Parent=m1
rec1\ttest\texon\t41\t60\t.\t+\t.\tID=e2;Parent=m1
rec1\ttest\tCDS\t1\t20\t.\t+\t0\tID=c1;Parent=m1
rec1\ttest\tCDS\t41\t60\t.\t+\t1\tID=c1;Parent=m1
rec1\ttest\tgene\t70\t120\t.\t-\t.\tID=g2
rec1\ttest\ttRNA\t70\t120\t.\t-\t.\tID=t1;Parent=g2
rec1\ttest\tmisc_feature\t5\t9\t.\t+\t.\tNote=no id
rec2\ttest\tgene\t1\t30\t.\t+\t.\tID=g3
rec2\ttest\tCDS\t1\t30\t.\t+\t0\tID=c3;Parent=g3
"""
TYPE_SETS = (
    None, set(), {"CDS"}, {"gene", "CDS"}, {"mRNA", "exon"}, {"tRNA", "misc_feature"},
    {"polyA_site", "repeat_region"}, {"absent"},
)


@pytest.fixture(params=["NC_013668", "nested"])
def source(request, tmp_path) -> tuple[Path, Path]:
    if request.param == "NC_013668":
        return INPUTS / "NC_013668.gff3", INPUTS / "NC_013668.fasta"
    gff, fasta = tmp_path / "nested.gff3", tmp_path / "nested.fasta"
    gff.write_text(NESTED_GFF3)
    fasta.write_text(f">rec1\n{'ATG' * 40}\n>rec2\n{'ATG' * 10}\n")
    return gff, fasta


def _features(record) -> list[dict]:
    """Every attribute of each feature, in order, for exact comparison."""
    return [vars(feature) for feature in record.features]


def test_type_filters_of_one_parse_select_the_source_features_in_order(source):
    gff, fasta = source
    parsed = parse_gff_fasta(str(gff), str(fasta), source_feature_catalogs=[])
    nested_types = [[feature.type for feature in record.features] for record in parsed]
    reference = {record.id: _normalize_gff3_multipart_features(record) for record in GFF.parse(str(gff))}
    for types in (*TYPE_SETS, *reversed(TYPE_SETS)):  # a filter must leave the parse as it found it
        for record in parsed:
            kept = filter_features_by_type(record, types, source_indexes=_feature_source_index_map(record.features))
            expected = [
                (feature.type, feature.id, feature.location, feature.qualifiers, ordinal)
                for ordinal, feature in enumerate(_iter_source_features(reference[record.id].features))
                if types is None or feature.type in types
            ]
            assert [
                (feature.type, feature.id, feature.location, feature.qualifiers, _source_feature_index(feature))
                for feature in kept.features
            ] == expected, types
    assert [[feature.type for feature in record.features] for record in parsed] == nested_types


@pytest.mark.parametrize("types", TYPE_SETS, ids=lambda types: "all" if types is None else "+".join(sorted(types)) or "none")
def test_resolving_other_candidate_types_reuses_the_parse_and_equals_the_type_filtered_load(
    monkeypatch, source, types,
):
    gff, fasta = source
    inputs = (RecordInput(GffFastaInputSource(gff, fasta), record_key="r", cardinality=RecordCardinality.ALL),)
    parses, parse = [], GFF.parse
    monkeypatch.setattr(GFF, "parse", lambda *args, **kwargs: parses.append(args) or parse(*args, **kwargs))
    parsed_sources = {}
    for candidates in (None, {"CDS"}, types):
        resolved = resolve_record_inputs(
            inputs, gff_candidate_features=None if candidates is None else sorted(candidates),
            gff_keep_all_features=False, parsed_sources=parsed_sources,
        )
    typed = load_gff_fasta([str(gff)], [str(fasta)], selected_features_set=types, source_feature_catalogs=[])
    assert len(parses) == 2  # one for the three resolutions, one for the typed load
    assert [_features(record) for record in resolved.records] == [_features(record) for record in typed]
