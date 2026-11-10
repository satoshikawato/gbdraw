"""Circular similarity rings search with one E-value scale per genome (PD-OI-018 revision 6, OV-3319).

The reference genome is every record of the files of the drawn records, the ring genome every record of the ring's
files. A ring search passes their lengths multiplied as ``-searchsp`` (TLOSATX in amino acids), except when both
genomes are one file of one record, which LOSAT's own search space already covers. The Web asks the same Python
owner (``plan_ring_scales_json``), so the CLI and the Web pass one value.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import pytest

from gbdraw.api import CircularDiagramOptions, CircularDiagramRequest, LosatSearchOptions, RecordInput
from gbdraw.api.request_render import plan_circular_request
from gbdraw.api.requests import GenBankInputSource, RecordCardinality
from gbdraw.comparisons.losat_jobs import LOSAT_SEARCH_SPACE_LIMIT, nucleotide_losat_cache_key, sha256_text

pytestmark = pytest.mark.circular

# Sequences without a repeat, long enough that a TLOSATX length counts whole codons and remainders.
_SEQUENCES = {
    "REF1": "ACGTTGCAAGTCCATGGATCCTTAGCGGCTAAGTCGATCGGAT" * 3,  # 129 bp
    "REF2": "TTGACCGGTAACGTAGCTAGGCTAACGTTAGC" * 2,  # 64 bp
    "RINGA": "GGATCCTTAGCGGCTAAGTCGATCGGATTT",  # 30 bp
    "RINGB": "CCATGGATCCTTAGCGGCTAAGTCGATCGGATACG",  # 35 bp
    "RINGC": "AAGTCCATGGATCCTTAGC",  # 19 bp
}


def _genbank(path: Path, *record_ids: str) -> Path:
    from io import StringIO

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    handle = StringIO()
    SeqIO.write(
        [
            SeqRecord(
                Seq(_SEQUENCES[record_id]), id=record_id, name=record_id, description=f"{record_id} record",
                annotations={"molecule_type": "DNA", "topology": "circular"},
            )
            for record_id in record_ids
        ],
        handle,
        "genbank",
    )
    path.write_text(handle.getvalue(), encoding="ascii")
    return path


def _fasta(path: Path, *record_ids: str) -> Path:
    path.write_text("".join(f">{record_id}\n{_SEQUENCES[record_id]}\n" for record_id in record_ids), encoding="ascii")
    return path


def _length(*record_ids: str) -> int:
    return sum(len(_SEQUENCES[record_id]) for record_id in record_ids)


@pytest.fixture
def ring_searches(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, object]]:
    """Each ring search's options, with a runtime that takes ``-searchsp``; the search returns no rows."""

    import gbdraw.comparisons.circular_losat as circular_losat_module

    searches: list[dict[str, object]] = []

    def fake_search(_spec, query_fasta, subject_fasta, *, options, **_kwargs):
        searches.append({"options": options, "query": query_fasta, "subject": subject_fasta})
        return ""

    monkeypatch.setattr(circular_losat_module, "run_losat_search", fake_search)
    monkeypatch.setattr(circular_losat_module, "losat_takes_genome_scale", lambda *_args, **_kwargs: True)
    return searches


def _plan(reference: list[Path], rings: list[list[Path]], program: str = "losatn"):
    return plan_circular_request(CircularDiagramRequest(
        records=tuple(
            RecordInput(source=GenBankInputSource(path=str(path)), cardinality=RecordCardinality.ALL)
            for path in reference
        ),
        options=CircularDiagramOptions(
            losat_search=LosatSearchOptions(program=program),
            conservation_sequence_files=[[str(path) for path in ring] for ring in rings],
        ),
    ))


def test_a_ring_of_two_files_passes_the_search_space_of_the_genome_pair(
    tmp_path: Path, ring_searches: list[dict[str, object]]
) -> None:
    reference = _genbank(tmp_path / "reference.gb", "REF1")
    ring = [_fasta(tmp_path / "chromosome.fasta", "RINGA"), _fasta(tmp_path / "plasmid.fasta", "RINGB")]
    plan = _plan([reference], [ring])

    (search,) = ring_searches
    space = _length("REF1") * _length("RINGA", "RINGB")
    assert search["options"].searchsp == space
    (entry,) = plan.losat_cache_entries
    assert entry["args"] == ["--task", "megablast", "--searchsp", str(space)]
    assert entry["key"] == nucleotide_losat_cache_key(
        program="blastn", args=("--task", "megablast", "--searchsp", str(space)),
        query_hash=sha256_text(str(search["query"])), subject_hash=sha256_text(str(search["subject"])),
        flow="circular-conservation",
    )


def test_a_one_record_ring_around_a_one_record_reference_keeps_the_args_and_key(
    tmp_path: Path, ring_searches: list[dict[str, object]]
) -> None:
    reference = _genbank(tmp_path / "reference.gb", "REF1")
    plan = _plan([reference], [[_fasta(tmp_path / "ring.fasta", "RINGA")]])

    (search,) = ring_searches
    assert search["options"].searchsp is None
    (entry,) = plan.losat_cache_entries
    assert entry["args"] == ["--task", "megablast"]


def test_a_reference_of_several_records_counts_every_record_of_its_files(
    tmp_path: Path, ring_searches: list[dict[str, object]]
) -> None:
    """One record drawn of a two-record file: the reference genome is both records; the ring is one record."""

    from gbdraw.io.record_select import parse_record_selector

    reference = _genbank(tmp_path / "reference.gb", "REF1", "REF2")
    plan_circular_request(CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(path=str(reference)), selector=parse_record_selector("#2")),),
        options=CircularDiagramOptions(
            losat_search=LosatSearchOptions(program="losatn"),
            conservation_sequence_files=[str(_fasta(tmp_path / "ring.fasta", "RINGC"))],
        ),
    ))

    (search,) = ring_searches
    assert search["options"].searchsp == _length("REF1", "REF2") * _length("RINGC")


def test_a_reference_of_two_files_is_one_genome(tmp_path: Path, ring_searches: list[dict[str, object]]) -> None:
    references = [_genbank(tmp_path / "chromosome.gb", "REF1"), _genbank(tmp_path / "plasmid.gb", "REF2")]
    _plan(references, [[_fasta(tmp_path / "ring.fasta", "RINGA")]])

    (search,) = ring_searches
    assert search["options"].searchsp == _length("REF1", "REF2") * _length("RINGA")


def test_tlosatx_counts_each_genome_in_amino_acids(tmp_path: Path, ring_searches: list[dict[str, object]]) -> None:
    reference = _genbank(tmp_path / "reference.gb", "REF1", "REF2")
    plan = _plan([reference], [[_fasta(tmp_path / "ring.fasta", "RINGB")]], program="tlosatx")

    (search,) = ring_searches
    space = (_length("REF1", "REF2") // 3) * (_length("RINGB") // 3)
    assert search["options"].searchsp == space
    (entry,) = plan.losat_cache_entries
    assert entry["args"] == ["--query-gencode", "1", "--db-gencode", "1", "--searchsp", str(space)]


def test_each_ring_has_its_own_genome(tmp_path: Path, ring_searches: list[dict[str, object]]) -> None:
    """Two rings of one file each: neither ring's file joins the other's genome or the reference's."""

    reference = _genbank(tmp_path / "reference.gb", "REF1")
    single = _fasta(tmp_path / "single.fasta", "RINGA")
    double = _fasta(tmp_path / "double.fasta", "RINGB", "RINGC")
    _plan([reference], [[single], [double], [single]])

    assert [search["options"].searchsp for search in ring_searches] == [
        None, _length("REF1") * _length("RINGB", "RINGC"),
    ]  # the third ring repeats the first search, which the raw cache serves


def test_a_runtime_without_searchsp_searches_as_before_with_a_warning(
    tmp_path: Path, ring_searches: list[dict[str, object]], monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    import gbdraw.comparisons.circular_losat as circular_losat_module

    probes: list[str] = []

    def no_genome_scale(program, **_kwargs):
        probes.append(str(program.name))
        return False

    monkeypatch.setattr(circular_losat_module, "losat_takes_genome_scale", no_genome_scale)
    reference = _genbank(tmp_path / "reference.gb", "REF1", "REF2")
    plan = _plan([reference], [[_fasta(tmp_path / "ring.fasta", "RINGA")]])

    assert probes == ["losatn"]
    (search,) = ring_searches
    assert search["options"].searchsp is None
    (entry,) = plan.losat_cache_entries
    assert entry["args"] == ["--task", "megablast"]


def test_a_genome_pair_past_what_losat_reads_is_searched_per_file_with_a_warning(
    caplog: pytest.LogCaptureFixture,
) -> None:
    from gbdraw.comparisons.circular_losat import RingGenomeFile, ring_search_spaces

    big = 3_100_000_000
    reference = [RingGenomeFile(name="human.fa", length=big, record_lengths=(big // 2, big - big // 2))]
    rings = [[RingGenomeFile(name="chimp.fa", length=big, record_lengths=(big,))],
             [RingGenomeFile(name="small.fa", length=100, record_lengths=(100,))]]
    scales = ring_search_spaces("losatn", reference, rings)

    assert big * big > LOSAT_SEARCH_SPACE_LIMIT
    assert scales.spaces == (None, big * 100)
    first, second = scales.warnings
    assert second is None
    assert first is not None and (first.first, first.second) == ("chimp.fa", "human.fa")
    assert "computed per file" in first.message()


def test_the_web_helper_plans_the_values_the_command_line_passes(
    tmp_path: Path, ring_searches: list[dict[str, object]]
) -> None:
    """``plan_ring_scales_json`` reads the files the Web stages and answers each ring's ``--searchsp``, as a decimal
    string, with the warnings: the values of the CLI's searches of the same files."""

    from gbdraw.web_support.comparison_sequences import plan_ring_scales_json

    reference = _genbank(tmp_path / "reference.gb", "REF1", "REF2")
    rings = [[_fasta(tmp_path / "a.fasta", "RINGA"), _fasta(tmp_path / "b.fasta", "RINGB")],
             [_fasta(tmp_path / "c.fasta", "RINGC")]]
    for program in ("losatn", "tlosatx"):
        ring_searches.clear()
        _plan([reference], rings, program=program)
        reply = json.loads(plan_ring_scales_json(
            program, "genbank", json.dumps([str(reference)]), json.dumps(["reference.gb"]),
            json.dumps([[str(path) for path in ring] for ring in rings]),
            json.dumps([[path.name for path in ring] for ring in rings]),
        ))
        assert reply == {
            "scales": [str(search["options"].searchsp) for search in ring_searches],
            "warnings": [],
        }


def test_the_web_helper_reads_a_gff3_reference_as_its_fasta_file(tmp_path: Path) -> None:
    """A GFF3 and FASTA reference counts every FASTA entry, as the CLI ring database does."""

    from gbdraw.web_support.comparison_sequences import plan_ring_scales_json

    fasta = _fasta(tmp_path / "reference.fasta", "REF1", "REF2")
    reply = json.loads(plan_ring_scales_json(
        "losatn", "fasta", json.dumps([str(fasta)]), json.dumps(["reference.fasta"]),
        json.dumps([[str(_fasta(tmp_path / "ring.fasta", "RINGA"))]]), json.dumps([["ring.fasta"]]),
    ))
    assert reply["scales"] == [str(_length("REF1", "REF2") * _length("RINGA"))]


def test_the_command_line_logs_the_per_file_warning_once(
    tmp_path: Path, ring_searches: list[dict[str, object]], monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    import gbdraw.comparisons.losat_jobs as jobs_module

    monkeypatch.setattr(jobs_module, "LOSAT_SEARCH_SPACE_LIMIT", 10)
    reference = _genbank(tmp_path / "reference.gb", "REF1", "REF2")
    ring = _fasta(tmp_path / "ring.fasta", "RINGA")
    with caplog.at_level(logging.WARNING):
        _plan([reference], [[ring], [ring]])
    warnings = [record.getMessage() for record in caplog.records if "computed per file" in record.getMessage()]
    assert warnings == [
        f"The search space of genomes reference.gb and ring.fasta ({_length('REF1', 'REF2') * _length('RINGA'):,}) "
        "is past the largest LOSAT reads (10), so their E-values are computed per file, not per genome."
    ]
    assert ring_searches[0]["options"].searchsp is None
