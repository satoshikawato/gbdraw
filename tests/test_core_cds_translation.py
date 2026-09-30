from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature

from gbdraw.core.sequence import translate_cds


MG1655 = Path(__file__).resolve().parent / "test_inputs" / "MG1655.gbk"


def _cds(location: object, **qualifiers: str) -> SeqFeature:
    return SeqFeature(
        location,
        type="CDS",
        qualifiers={key: [value] for key, value in qualifiers.items()},
    )


@pytest.mark.parametrize(
    ("location", "qualifiers", "nucleotides", "expected"),
    [
        (FeatureLocation(0, 9, strand=1), {"transl_table": "11"}, "gtgaaatag", "MK"),
        (FeatureLocation(0, 9, strand=1), {}, "GTGAAATAG", "VK"),
        (FeatureLocation(0, 9, strand=1), {"transl_table": "11", "pseudo": ""}, "GTGAAATAG", "VK"),
        (FeatureLocation(0, 9, strand=1), {"transl_table": "11", "pseudogene": "unitary"}, "TTGAAATAG", "LK"),
        (FeatureLocation(BeforePosition(0), 9, strand=1), {"transl_table": "11"}, "GTGAAATAG", "VK"),
        (FeatureLocation(0, AfterPosition(9), strand=-1), {"transl_table": "11"}, "GTGAAATAG", "VK"),
        (FeatureLocation(0, AfterPosition(9), strand=1), {"transl_table": "11"}, "GTGAAATAG", "MK"),
        (
            CompoundLocation([FeatureLocation(10, 13, strand=-1), FeatureLocation(0, 6, strand=-1)]),
            {"transl_table": "11"},
            "TTGAAATAG",
            "MK",
        ),
        (FeatureLocation(0, 10, strand=1), {"transl_table": "11", "codon_start": "2"}, "AGTGAAATAG", "VK"),
        (FeatureLocation(0, 11, strand=1), {"transl_table": "11", "phase": "2"}, "AAGTGAAATAG", "VK"),
        (FeatureLocation(0, 11, strand=1), {"transl_table": "11", "codon_start": "3", "phase": "0"}, "AAATGAAATAG", "MK"),
        (FeatureLocation(0, 11, strand=1), {"transl_table": "11"}, "GTGAAATAGCC", "MK"),
        (FeatureLocation(0, 12, strand=1), {"transl_table": "11"}, "GTGNNNAARTAG", "MXK"),
        (FeatureLocation(0, 9, strand=1), {"transl_table": "11"}, "NTGAAATAG", "XK"),
    ],
)
def test_translate_cds_applies_the_insdc_start_codon_rule(
    location: object,
    qualifiers: dict[str, str],
    nucleotides: str,
    expected: str,
) -> None:
    assert translate_cds(_cds(location, **qualifiers), nucleotides) == expected


@pytest.mark.parametrize(
    ("qualifiers", "nucleotides", "message"),
    [
        ({"codon_start": "x"}, "ATGAAATAG", "codon_start is invalid: x"),
        ({"codon_start": "4"}, "ATGAAATAG", "codon_start is outside 1..3: 4"),
        ({"phase": "3"}, "ATGAAATAG", "phase is outside 0..2: 3"),
        ({"transl_table": "7"}, "ATGAAATAG", "transl_table is invalid: 7"),
        ({}, "AT", "coding sequence has no complete codon"),
    ],
)
def test_translate_cds_rejects_invalid_frames_and_tables(
    qualifiers: dict[str, str],
    nucleotides: str,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        translate_cds(_cds(FeatureLocation(0, len(nucleotides), strand=1), **qualifiers), nucleotides)


def test_translate_cds_can_require_whole_codons() -> None:
    feature = _cds(FeatureLocation(0, 11, strand=1), transl_table="11")

    with pytest.raises(ValueError, match="not divisible by 3"):
        translate_cds(feature, "GTGAAATAGCC", require_whole_codons=True)


@pytest.mark.slow
def test_translate_cds_matches_every_mg1655_translation() -> None:
    record = SeqIO.read(MG1655, "genbank")
    compared = 0
    mismatches: list[str] = []
    for feature in record.features:
        qualifiers = feature.qualifiers
        if (
            feature.type != "CDS"
            or "translation" not in qualifiers
            or "transl_except" in qualifiers
            or "pseudo" in qualifiers
        ):
            continue
        compared += 1
        if translate_cds(feature, feature.extract(record.seq)) != qualifiers["translation"][0]:
            mismatches.append(qualifiers["locus_tag"][0])

    assert compared > 4000
    assert mismatches == []
