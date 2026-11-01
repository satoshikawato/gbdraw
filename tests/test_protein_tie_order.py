"""Score ties between protein hits are broken on the CDS, never on its runtime handle (OV-2370, OV-2380).

A runtime handle (``h_...``) hashes the record instance key, which is the
record key on the CLI and a per-upload sequence uid on the Web. Two equal
inputs therefore carry different handles on the two surfaces, and a tie broken
on the handle text keeps a different subject on each.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.protein_colinearity import (
    CdsProtein,
    build_pairwise_protein_blastp_comparisons,
    build_rbh_orthogroup_protein_blastp_comparisons,
    extract_web_stable_cds_proteins,
    select_top_hits_per_query,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS

pytestmark = pytest.mark.linear

# Instance-key sets as the CLI ("record-N") and Web pages (sequence uids) give them.
KEY_SETS = (
    ("record-1", "record-2"),
    ("seq-3f9c2a", "seq-0b71de"),
    ("seq-e4d0aa", "seq-91c3b2"),
    ("seq-77aa01", "seq-2c5e9f"),
)
# One query CDS; three identical subject CDSs, as a repeated transposase gives.
SUBJECT_STARTS = (0, 30, 60)


def _records() -> tuple[SeqRecord, ...]:
    def cds(start: int) -> SeqFeature:
        return SeqFeature(FeatureLocation(start, start + 9, strand=1), type="CDS", qualifiers={"translation": ["MK"]})

    query = SeqRecord(Seq("ATGAAATAG" * 20), id="query_record")
    query.features = [cds(0)]
    subject = SeqRecord(Seq("ATGAAATAG" * 20), id="subject_record")
    subject.features = [cds(start) for start in SUBJECT_STARTS]
    return (query, subject)


def _hit(query: str, subject: str) -> dict[str, object]:
    return {
        "query": query, "subject": subject, "identity": 90.0, "alignment_length": 2,
        "mismatches": 0, "gap_opens": 0, "qstart": 1, "qend": 2, "sstart": 1, "send": 2,
        "evalue": 1e-30, "bitscore": 200.0,
    }


def _fasta_ids(fasta: str) -> list[str]:
    return [line[1:].split()[0] for line in fasta.splitlines() if line.startswith(">")]


def _tied_runner(query_fasta: str, subject_fasta: str) -> pd.DataFrame:
    """Every query protein hits every subject protein with the same scores."""

    rows = [_hit(query, subject) for query in _fasta_ids(query_fasta) for subject in _fasta_ids(subject_fasta)]
    return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)


def _extraction(keys: tuple[str, str]):
    records = _records()
    return records, extract_web_stable_cds_proteins(
        records, record_instance_keys=keys, record_source_ids=tuple(record.id for record in records)
    )


def _link_spans(frame: pd.DataFrame) -> list[tuple[int, int, int, int]]:
    return sorted(
        (int(row.qstart), int(row.qend), int(row.sstart), int(row.send))
        for row in frame.itertuples(index=False)
    )


def _pairwise_spans(keys: tuple[str, str]) -> list[tuple[int, int, int, int]]:
    records, extraction = _extraction(keys)
    result = build_pairwise_protein_blastp_comparisons(
        records, runner=_tied_runner, protein_extraction=extraction,
        max_hits=1, evalue=1.0, bitscore=0.0, identity=0.0,
    )
    return _link_spans(result.comparisons[0])


def _similarity_group_spans(keys: tuple[str, str]) -> list[tuple[int, int, int, int]]:
    records, extraction = _extraction(keys)
    result = build_rbh_orthogroup_protein_blastp_comparisons(
        records, runner=_tied_runner, protein_extraction=extraction,
        max_related_edges_per_orthogroup=1, evalue=1.0, bitscore=0.0, identity=0.0,
    )
    return _link_spans(result.comparisons[0])


def test_top_hits_cut_keeps_the_first_cds_on_a_score_tie() -> None:
    def protein(protein_id: str, start: int) -> CdsProtein:
        return CdsProtein(
            protein_id=protein_id, record_index=1, feature_index=start, record_id="subject_record",
            label=protein_id, start=start, end=start + 9, strand=1, sequence="MK", protein_length=2, feature_svg_id=None,
        )

    # The handle text sorts the other way round from the CDS location.
    protein_map = {"h_z": protein("h_z", 0), "h_a": protein("h_a", 30)}
    hits = pd.DataFrame.from_records([_hit("q", "h_a"), _hit("q", "h_z")], columns=COMPARISON_COLUMNS)

    selected = select_top_hits_per_query(hits, max_hits=1, protein_map=protein_map)

    assert selected["subject"].tolist() == ["h_z"]


@pytest.mark.parametrize("keys", KEY_SETS[1:])
def test_pairwise_tie_keeps_the_same_cds_for_every_record_key(keys: tuple[str, str]) -> None:
    assert _pairwise_spans(keys) == _pairwise_spans(KEY_SETS[0])
    # The first subject CDS wins the tie.
    assert [span[2] for span in _pairwise_spans(keys)] == [1]


@pytest.mark.parametrize("keys", KEY_SETS[1:])
def test_similarity_group_tie_keeps_the_same_links_for_every_record_key(keys: tuple[str, str]) -> None:
    assert _similarity_group_spans(keys) == _similarity_group_spans(KEY_SETS[0])


def _web_protein_entry(protein: CdsProtein) -> dict[str, object]:
    entry = {
        field.name: getattr(protein, field.name)
        for field in dataclasses.fields(protein)
        if isinstance(getattr(protein, field.name), (str, int, float, type(None)))
    }
    entry["protein_length"] = len(protein.sequence)
    return entry


def test_web_pairwise_tie_equals_the_cli(tmp_path: Path, stage_web_losatp_transport) -> None:
    """The Web helper (another record key, as a Web page has) keeps the CLI's CDS."""

    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(helpers_js.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    _records_, extraction = _extraction(KEY_SETS[3])
    query_ids = [protein.protein_id for protein in extraction.proteins_by_record[0]]
    subject_ids = [protein.protein_id for protein in extraction.proteins_by_record[1]]
    hits = pd.DataFrame.from_records(
        [_hit(query, subject) for query in query_ids for subject in subject_ids], columns=COMPARISON_COLUMNS
    )
    payload = {
        "records": [
            {
                "recordIndex": index,
                "recordId": proteins[0].record_id,
                "proteinMap": {protein.protein_id: _web_protein_entry(protein) for protein in proteins},
                "runtimeBindingHash": f"runtime-{index}",
                "displayBindingHash": f"display-{index}",
            }
            for index, proteins in enumerate(extraction.proteins_by_record)
        ],
        "pairs": [
            {
                "pairIndex": 0, "queryIndex": 0, "subjectIndex": 1, "cacheKey": "pair-0-1",
                "blastText": hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            }
        ],
    }
    pairs_path, raw_tsv_path = stage_web_losatp_transport(tmp_path, payload)
    result = json.loads(str(namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path), str(raw_tsv_path), "pairwise", 1, 0, "1", 0, 0,
    )))
    assert "error" not in result, result

    from tests.utils.comparison_tables import converted_pair_rows

    web_spans = sorted(
        (int(row["qstart"]), int(row["qend"]), int(row["sstart"]), int(row["send"]))
        for row in converted_pair_rows(result["pairs"][0])
    )
    assert web_spans == _pairwise_spans(KEY_SETS[0])
