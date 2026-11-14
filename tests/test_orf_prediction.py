"""The ORF owner: orfind.py's rules with the four D-09 changes (origin, 3'-partial, start set, nested sweep)."""

from __future__ import annotations

import hashlib
import io
import random
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.orf_prediction import (
    OrfPredictionParameters,
    PredictedOrf,
    predict_orfs,
    predicted_orf_features,
)
from gbdraw.exceptions import ValidationError
from gbdraw.version import __version__

INPUTS = Path(__file__).resolve().parent / "test_inputs"


def _spans(orfs: tuple[PredictedOrf, ...]) -> list[tuple[int, int, int, bool, str]]:
    return [
        (orf.start, orf.end, orf.strand, orf.partial_3prime, orf.protein)
        for orf in orfs
    ]


def _digest(orfs: tuple[PredictedOrf, ...]) -> str:
    text = "".join(
        f"{orf.start}\t{orf.end}\t{orf.strand}\t{orf.protein}\n" for orf in orfs
    )
    return hashlib.sha256(text.encode()).hexdigest()


# orfind.py (satoshikawato/bio_small_scripts) with -m 50, nested ORFs removed: the sha256 of its ORFs sorted by
# start, end, strand, one "start<TAB>end<TAB>strand<TAB>protein" line each (prh/logs/orfind-digests.log).
@pytest.mark.parametrize(
    ("name", "table", "count", "digest"),
    [
        (
            "NC_001416.gb",
            11,
            110,
            "171b25c57ed6f9acbf0197f339786708c1b377117a6b50660883f0a53cae7ffc",
        ),
        (
            "HmmtDNA.gbk",
            2,
            16,
            "8ad292626c01bef26d6259522dfffb365bb28d621bd837d4c1c46cbf9dcf6a46",
        ),
        (
            "NZ_CP046414.1.gbk",
            11,
            119,
            "9e41549c182047700d4c2aa64687eed91cbfaad145d33bb587318f6a124ca472",
        ),
    ],
)
def test_a_linear_reading_gives_orfinds_orfs(
    name: str, table: int, count: int, digest: str
) -> None:
    record = SeqIO.read(INPUTS / name, "genbank")
    orfs = predict_orfs(
        record.seq,
        topology="linear",
        parameters=OrfPredictionParameters(genetic_code=table),
    )
    assert len(orfs) == count
    assert _digest(orfs) == digest


def test_a_circular_plasmid_differs_from_orfind_only_over_the_origin() -> None:
    record = SeqIO.read(INPUTS / "NZ_CP046414.1.gbk", "genbank")
    parameters = OrfPredictionParameters(genetic_code=11)
    linear = {
        orf[:3]
        for orf in _spans(
            predict_orfs(record.seq, topology="linear", parameters=parameters)
        )
    }
    circular = {
        orf[:3]
        for orf in _spans(
            predict_orfs(record.seq, topology="circular", parameters=parameters)
        )
    }
    # orfind cut the ORF at the record end (no stop) and kept a minus-strand piece from base 1; read across the
    # origin, the ORF ends 490 bases after it and the minus-strand ORF over the origin lies inside it.
    assert circular - linear == {(67127, 69338, 1)}
    assert linear - circular == {(67127, 68848, 1), (0, 207, -1)}


ORIGIN_SEQUENCE = "CCCCCCCCCCCCCCCTAA" + "C" * 29 + "A" + "ATG" + "CCC" * 10
SHORT = OrfPredictionParameters(min_protein_length=5)


def test_an_orf_crosses_the_origin_of_a_circular_record() -> None:
    assert _spans(
        predict_orfs(ORIGIN_SEQUENCE, topology="circular", parameters=SHORT)
    ) == [(48, 99, 1, False, "M" + "P" * 15)]
    reverse = str(Seq(ORIGIN_SEQUENCE).reverse_complement())
    assert _spans(predict_orfs(reverse, topology="circular", parameters=SHORT)) == [
        (63, 114, -1, False, "M" + "P" * 15)
    ]


def test_an_orf_without_a_stop_codon_is_3prime_partial_on_a_linear_record() -> None:
    assert _spans(
        predict_orfs(ORIGIN_SEQUENCE, topology="linear", parameters=SHORT)
    ) == [(48, 81, 1, True, "M" + "P" * 10)]
    reverse = str(Seq(ORIGIN_SEQUENCE).reverse_complement())
    assert _spans(predict_orfs(reverse, topology="linear", parameters=SHORT)) == [
        (0, 33, -1, True, "M" + "P" * 10)
    ]


def test_a_frame_without_a_start_codon_gives_no_orf() -> None:
    # OV-3567: orfind read such a frame from its last codon, giving ORFs at bases -2 and 0..189 of a 186-base
    # record that contained, and so removed, the real ORF.
    sequence = "ATG" + "CCC" * 60 + "TAA"
    assert _spans(
        predict_orfs(sequence, topology="linear", parameters=OrfPredictionParameters())
    ) == [(0, 186, 1, False, "M" + "P" * 60)]


def test_orfs_with_the_same_span_are_both_kept() -> None:
    # OV-3568: each contains the other, so orfind removed both and returned no ORF.
    sequence = "ATGCCCGGGCCCATGGGGGGGGATGGATGGCCCCATGGGCCATGGCAT"
    orfs = predict_orfs(
        sequence,
        topology="linear",
        parameters=OrfPredictionParameters(min_protein_length=10),
    )
    assert [orf[:4] for orf in _spans(orfs)] == [(0, 48, -1, True), (0, 48, 1, True)]


def test_the_table_start_codons_add_alternative_starts() -> None:
    sequence = "CC" + "GTG" + "CCC" * 6 + "TAA" + "CC"
    methionine = OrfPredictionParameters(genetic_code=11, min_protein_length=5)
    table = OrfPredictionParameters(
        genetic_code=11, start_codons="table", min_protein_length=5
    )
    assert predict_orfs(sequence, topology="linear", parameters=methionine) == ()
    assert _spans(predict_orfs(sequence, topology="linear", parameters=table)) == [
        (2, 26, 1, False, "MPPPPPP")
    ]
    # In table 2 ATA is read as M, so it starts an ORF without the table's other starts (orfind).
    mitochondrial = OrfPredictionParameters(genetic_code=2, min_protein_length=5)
    assert _spans(
        predict_orfs(
            "ATA" + "CCC" * 6 + "TAA", topology="linear", parameters=mitochondrial
        )
    ) == [(0, 24, 1, False, "MPPPPPP")]


@pytest.mark.parametrize(
    ("table", "protein"),
    [(1, "MPPPPPP"), (4, "MPPPPPPWPPP"), (2, "MPPPPPPWPPP")],
)
def test_stop_codons_follow_the_genetic_code(table: int, protein: str) -> None:
    sequence = "ATG" + "CCC" * 6 + "TGA" + "CCC" * 3 + "TAG"
    parameters = OrfPredictionParameters(genetic_code=table, min_protein_length=5)
    assert [
        orf.protein
        for orf in predict_orfs(sequence, topology="linear", parameters=parameters)
    ] == [protein]


def test_table_2_stops_at_aga() -> None:
    sequence = "ATG" + "CCC" * 6 + "AGA" + "CCC" * 3 + "TAG"
    standard = OrfPredictionParameters(genetic_code=1, min_protein_length=5)
    mitochondrial = OrfPredictionParameters(genetic_code=2, min_protein_length=5)
    assert [
        orf.protein
        for orf in predict_orfs(sequence, topology="linear", parameters=standard)
    ] == ["MPPPPPPRPPP"]
    assert [
        orf.protein
        for orf in predict_orfs(sequence, topology="linear", parameters=mitochondrial)
    ] == ["MPPPPPP"]


def test_case_u_and_ambiguity_codes_read_as_biopython_translates_them() -> None:
    parameters = OrfPredictionParameters(min_protein_length=3)
    assert [
        orf.protein
        for orf in predict_orfs(
            "atgnnngcnccctar", topology="linear", parameters=parameters
        )
    ] == ["MXAP"]
    assert [
        orf.protein
        for orf in predict_orfs(
            "AUGCCCCCCCCCUAA", topology="linear", parameters=parameters
        )
    ] == ["MPPP"]
    assert [
        orf.protein
        for orf in predict_orfs(
            "ATGCC-CCCCCCTAA", topology="linear", parameters=parameters
        )
    ] == ["MXPP"]


def test_a_circular_frame_without_a_stop_codon_gives_no_orf() -> None:
    assert predict_orfs("ATG" + "CCC" * 30, topology="circular", parameters=SHORT) == ()


def _arc(orf: PredictedOrf, length: int) -> frozenset[int]:
    return frozenset(position % length for position in range(orf.start, orf.end))


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_the_sweep_removes_exactly_the_orfs_inside_another(topology: str) -> None:
    rng = random.Random(20261114)
    for _ in range(40):
        sequence = "".join(rng.choice("ACGT") for _ in range(rng.randrange(150, 400)))
        everything = predict_orfs(
            sequence,
            topology=topology,
            parameters=OrfPredictionParameters(min_protein_length=5, keep_nested=True),
        )
        kept = predict_orfs(
            sequence,
            topology=topology,
            parameters=OrfPredictionParameters(min_protein_length=5),
        )
        arcs = {orf: _arc(orf, len(sequence)) for orf in everything}
        expected = tuple(
            orf
            for orf in everything
            if not any(arcs[orf] < other for other in arcs.values())
        )
        assert kept == expected
        assert all(orf.end - orf.start <= len(sequence) for orf in everything)
        features = predicted_orf_features(
            everything, record_length=len(sequence), genetic_code=1
        )
        for orf, feature in zip(everything, features):
            coding = feature.extract(Seq(sequence))
            # A 3'-partial ORF runs to the record end, over the last incomplete codon (orfind).
            extra = len(coding) - len(orf.protein) * 3
            assert extra in ((0, 1, 2) if orf.partial_3prime else (3,))
            assert str(coding[:3].translate()) == "M"
            assert (
                str(coding[: len(coding) - len(coding) % 3].translate(to_stop=True))
                == orf.protein
            )


def _genbank(sequence: str, orfs: tuple[PredictedOrf, ...], topology: str) -> str:
    record = SeqRecord(
        Seq(sequence),
        id="rec1",
        name="rec1",
        annotations={"molecule_type": "DNA", "topology": topology},
    )
    record.features = predicted_orf_features(
        orfs, record_length=len(sequence), genetic_code=11
    )
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


def test_features_are_cds_with_fuzzy_partial_ends_and_origin_joins() -> None:
    reverse = str(Seq(ORIGIN_SEQUENCE).reverse_complement())
    linear = _genbank(
        ORIGIN_SEQUENCE,
        predict_orfs(ORIGIN_SEQUENCE, topology="linear", parameters=SHORT),
        "linear",
    )
    assert "     CDS             49..>81\n" in linear
    linear_minus = _genbank(
        reverse, predict_orfs(reverse, topology="linear", parameters=SHORT), "linear"
    )
    assert "     CDS             complement(<1..33)\n" in linear_minus
    circular = _genbank(
        ORIGIN_SEQUENCE,
        predict_orfs(ORIGIN_SEQUENCE, topology="circular", parameters=SHORT),
        "circular",
    )
    assert "     CDS             join(49..81,1..18)\n" in circular
    circular_minus = _genbank(
        reverse,
        predict_orfs(reverse, topology="circular", parameters=SHORT),
        "circular",
    )
    assert "     CDS             complement(join(64..81,1..33))\n" in circular_minus


def test_feature_qualifiers_say_the_cds_is_predicted_and_carry_no_label() -> None:
    orfs = predict_orfs(ORIGIN_SEQUENCE, topology="circular", parameters=SHORT)
    (feature,) = predicted_orf_features(
        orfs, record_length=len(ORIGIN_SEQUENCE), genetic_code=11
    )
    assert feature.type == "CDS"
    assert feature.qualifiers == {
        "transl_table": ["11"],
        "translation": ["M" + "P" * 15],
        "inference": [f"ab initio prediction:gbdraw:{__version__}"],
    }
    assert (
        str(feature.extract(Seq(ORIGIN_SEQUENCE)).translate(table=11, to_stop=True))
        == "M" + "P" * 15
    )


@pytest.mark.parametrize(
    "arguments",
    [{"genetic_code": 7}, {"start_codons": "all"}, {"min_protein_length": 0}],
)
def test_parameters_are_checked(arguments: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        OrfPredictionParameters(**arguments)  # type: ignore[arg-type]


def test_topology_is_checked() -> None:
    with pytest.raises(ValidationError):
        predict_orfs("ATG", topology="ring", parameters=SHORT)  # type: ignore[arg-type]
