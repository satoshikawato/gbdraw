"""Shared LOSAT job-plan and FASTA vectors (design 3.6, D15).

tests/web/linear-sources.test.mjs runs the same cases against the Web code, so
the CLI and the Web search the same databases under the same raw keys.
"""

from __future__ import annotations

from io import StringIO
import json
from pathlib import Path

from Bio import SeqIO
import pytest

from gbdraw.comparisons.losat_jobs import (
    LosatBatch,
    losat_search_frame_fasta,
    nucleotide_losat_cache_key,
    plan_losat_jobs,
    prepare_losat_batches,
    sha256_text,
    split_losat_batch_result,
)
from gbdraw.comparisons.losat_runtime import LosatSearchArgs, losat_cache_args
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector, select_record
from gbdraw.io.regions import apply_region_specs, parse_region_specs

FIXTURES = Path(__file__).parent / "fixtures"
_PROGRAMS = {"blastn": "losatn", "tblastx": "tlosatx"}


def _cases(name: str) -> list[dict]:
    return json.loads((FIXTURES / name).read_text(encoding="utf-8"))["cases"]


def _case_id(case: dict) -> str:
    return case["name"]


@pytest.mark.parametrize("case", _cases("losat_fasta_extraction_cases.json"), ids=_case_id)
def test_fasta_extraction_matches_web(case: dict) -> None:
    records = list(SeqIO.parse(StringIO(case["text"]), case["fmt"]))
    selector = parse_record_selector(case["recordSelector"])
    records = [records[0]] if selector is None else select_record(records, selector)
    if case["regionSpec"]:
        records = apply_region_specs(records, parse_region_specs([case["regionSpec"]]))
    fasta = losat_search_frame_fasta(records[0])
    expected = case["expected"]
    assert fasta == expected["fasta"]
    assert sha256_text(fasta) == expected["hash"]
    assert records[0].id == expected["recordId"]


def test_search_frame_undoes_a_display_reverse_complement() -> None:
    case = _cases("losat_fasta_extraction_cases.json")[0]
    record = next(SeqIO.parse(StringIO(case["text"]), "genbank"))
    reversed_record = record.reverse_complement(id=record.id)
    reversed_record.annotations["gbdraw_coord_step"] = -1
    assert losat_search_frame_fasta(reversed_record) == case["expected"]["fasta"]


def _build_args(case: dict):
    program = _PROGRAMS[case["program"]]
    records = case["records"]

    def build(query: int, subject: int) -> list[str]:
        if program == "losatn":
            options = LosatSearchArgs(task=case["task"])
        else:
            options = LosatSearchArgs(
                query_gencode=records[query]["gencode"],
                db_gencode=records[subject]["gencode"],
            )
        return losat_cache_args(program, options)

    return build


def _plan(case: dict) -> tuple[LosatBatch, ...]:
    records = case["records"]
    jobs = plan_losat_jobs(
        source_ids=[record["source"] for record in records],
        specs=[tuple(spec) for spec in case["specs"]],
        build_args=_build_args(case),
    )
    return prepare_losat_batches(
        jobs,
        uids=[record["uid"] for record in records],
        record_fasta=lambda index: records[index]["fasta"],
    )


@pytest.mark.parametrize("case", _cases("losat_job_plan_cases.json"), ids=_case_id)
def test_job_plan_matches_web(case: dict) -> None:
    batches = _plan(case)
    actual = [
        {
            "scope": batch.job.scope,
            "args": list(batch.job.args),
            "specs": [list(spec) for spec in batch.job.specs],
            "queryIndexes": list(batch.query.indexes),
            "subjectIndexes": list(batch.subject.indexes),
            "queryIds": {key: list(value) for key, value in batch.query.ids.items()},
            "subjectIds": {key: list(value) for key, value in batch.subject.ids.items()},
            "queryHash": batch.query.hash,
            "subjectHash": batch.subject.hash,
            "searchContext": batch.search_context,
        }
        for batch in batches
    ]
    assert actual == case["expected"]["jobs"]


@pytest.mark.parametrize("case", _cases("losat_job_plan_cases.json"), ids=_case_id)
def test_raw_keys_match_web(case: dict) -> None:
    records = case["records"]
    build_args = _build_args(case)
    batch_by_spec = {spec: batch for batch in _plan(case) for spec in batch.job.specs}
    keys = [
        [
            query,
            subject,
            nucleotide_losat_cache_key(
                program=case["program"],
                args=build_args(query, subject),
                query_hash=sha256_text(records[query]["fasta"]),
                subject_hash=sha256_text(records[subject]["fasta"]),
                search_context=batch_by_spec[(query, subject)].search_context,
            ),
        ]
        for query, subject in (tuple(spec) for spec in case["specs"])
    ]
    assert keys == case["expected"]["rawKeys"]


def test_split_restores_record_ids_and_rejects_unknown_endpoints() -> None:
    case = next(
        item for item in _cases("losat_job_plan_cases.json") if len(item["specs"]) == 2
        and item["expected"]["jobs"][0]["scope"] == "between-sources"
    )
    batch = _plan(case)[0]
    query_ids = list(batch.query.ids)
    subject_id = next(iter(batch.subject.ids))
    text = "".join(
        "\t".join([query_id, subject_id, "99", "10", "0", "0", "1", "10", "1", "10", "1e-5", "40"]) + "\n"
        for query_id in query_ids
    )
    split = split_losat_batch_result(text, batch, batch.job.specs)
    for (query, _subject), tsv in split.items():
        assert tsv.split("\t")[0] == batch.query.ids[query_ids[batch.query.indexes.index(query)]][1]
    with pytest.raises(ValidationError) as excinfo:
        split_losat_batch_result("unknown\tunknown\n", batch, batch.job.specs)
    assert excinfo.value.diagnostic == {"code": "LOSAT_RUNTIME", "reason": "OUTPUT", "row": 1}


# LOSATP uses the same plan (design D7, PR-5).


def _protein_cases() -> list[dict]:
    return json.loads((FIXTURES / "losat_job_plan_cases.json").read_text(encoding="utf-8"))[
        "proteinCases"
    ]


@pytest.mark.parametrize("case", _protein_cases(), ids=_case_id)
def test_losatp_specs_match_web(case: dict) -> None:
    from gbdraw.comparisons.losat_jobs import losatp_job_specs

    specs = losatp_job_specs(
        case["mode"],
        record_count=len(case["records"]),
        pairs=[tuple(pair) for pair in case["edges"]],
        infer_orthogroups=case["inferOrthogroups"],
        search_scope=case["searchScope"],
    )
    assert [list(spec) for spec in specs] == case["specs"]


@pytest.mark.parametrize("case", _protein_cases(), ids=_case_id)
def test_losatp_jobs_and_raw_keys_match_web(case: dict) -> None:
    from gbdraw.analysis.protein_colinearity import (
        ProteinLosatPairIdentity,
        build_protein_losat_cache_key,
    )
    from gbdraw.comparisons.losat_jobs import losatp_job_specs

    records = case["records"]
    specs = losatp_job_specs(
        case["mode"],
        record_count=len(records),
        pairs=[tuple(pair) for pair in case["edges"]],
        infer_orthogroups=case["inferOrthogroups"],
        search_scope=case["searchScope"],
    )
    batches = prepare_losat_batches(
        plan_losat_jobs(
            source_ids=[record["source"] for record in records],
            specs=specs,
            build_args=lambda _query, _subject: case["args"],
        ),
        uids=[record["uid"] for record in records],
        record_fasta=lambda index: records[index]["fasta"],
        protein=True,
    )
    assert [
        {
            "scope": batch.job.scope,
            "args": list(batch.job.args),
            "specs": [list(spec) for spec in batch.job.specs],
            "queryIndexes": list(batch.query.indexes),
            "subjectIndexes": list(batch.subject.indexes),
            "queryIds": {key: list(value) for key, value in batch.query.ids.items()},
            "subjectIds": {key: list(value) for key, value in batch.subject.ids.items()},
            "queryHash": batch.query.hash,
            "subjectHash": batch.subject.hash,
            "searchContext": batch.search_context,
        }
        for batch in batches
    ] == case["expected"]["jobs"]

    def identity(query: int, subject: int) -> ProteinLosatPairIdentity:
        return ProteinLosatPairIdentity(
            query_protein_set_hash=records[query]["identity"]["proteinSetHash"],
            subject_protein_set_hash=records[subject]["identity"]["proteinSetHash"],
            query_runtime_binding_hash=records[query]["identity"]["runtimeBindingHash"],
            subject_runtime_binding_hash=records[subject]["identity"]["runtimeBindingHash"],
            query_record_instance_key=records[query]["uid"],
            subject_record_instance_key=records[subject]["uid"],
        )

    context = {spec: batch.search_context for batch in batches for spec in batch.job.specs}
    assert [
        [query, subject, build_protein_losat_cache_key(
            identity(query, subject), args=case["args"], search_context=context[(query, subject)]
        )]
        for query, subject in specs
    ] == case["expected"]["rawKeys"]
