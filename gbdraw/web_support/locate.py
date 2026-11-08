"""Worker helpers of Find in other files: the locate search plan and the plan of its rows, as JSON.

The page sends the reference File, record and windows, the target Files and
the program; :func:`locate_search_plan_json` replies with one LOSAT job per
target File, which the page runs. :func:`locate_windows_json` takes the same
input with each job's raw rows (or its error) and the dialog's thresholds and
window rules, and replies with the plan and the window pairs the new drawing
compares (the pairs ``gbdraw locate`` writes as its comparisons table). Both call the locate owner
(:mod:`gbdraw.comparisons.locate`) as ``gbdraw locate`` does, so the Web has no
copy of the search plan or of the window rules.
"""

from __future__ import annotations

from dataclasses import replace
from typing import TypedDict

from gbdraw.comparisons.locate import (
    LocateDrawing,
    LocateInputFile,
    LocateOptions,
    LocatePlanJson,
    LocateProgramName,
    LocateSearch,
    LocateSearchResult,
    LocateWindowPlace,
    finish_locate,
    locate_drawing,
    locate_results,
    prepare_locate_search,
)
from gbdraw.io.regions import RecordWindow
from gbdraw.mode_profiles import ComparisonThresholds


class LocateFileInput(TypedDict):
    """One File: its key in the plan and its input paths (GenBank, GFF3 with FASTA, or FASTA alone)."""

    key: str
    gbk: str | None
    gff: str | None
    fasta: str | None


class LocateWindowInput(TypedDict):
    """One reference window, 1-based and inclusive; ``start > end`` crosses the origin of a circular record."""

    start: int
    end: int
    reverseComplement: bool


class LocateSearchInput(TypedDict):
    """The reference File, its record (``None`` for its only record) and that record's key in the drawing it was
    chosen from, its windows in column order (none: the whole record), the target Files and the program."""

    reference: LocateFileInput
    recordId: str | None
    recordKey: str
    windows: list[LocateWindowInput]
    targets: list[LocateFileInput]
    program: LocateProgramName


class LocateJobReply(TypedDict):
    fileKey: str
    queryFasta: str
    subjectFasta: str
    rawKey: str


class LocateSearchPlanReply(TypedDict):
    """The jobs to run with ``program`` and the search options in the raw-cache form."""

    program: LocateProgramName
    cacheArgs: list[str]
    jobs: list[LocateJobReply]


class LocateOutcomeInput(TypedDict):
    """One job's raw outfmt 6 text as the runtime wrote it, or the error that stopped it."""

    fileKey: str
    text: str
    error: str | None


class LocateOptionsInput(TypedDict):
    evalue: float
    bitscore: float
    identity: float
    alignmentLength: int
    flank: int
    mergeDistance: int
    maxWindows: int


class LocateWindowsInput(LocateSearchInput):
    outcomes: list[LocateOutcomeInput]
    options: LocateOptionsInput


class LocateWindowName(TypedDict):
    """One window of the drawing: its File, its record by key, and its span (``start > end`` across the origin)."""

    fileKey: str
    recordKey: str
    start: int
    end: int


class LocateWindowPairReply(TypedDict):
    """Two windows of adjacent rows that share a reference gene, the upper row's first."""

    query: LocateWindowName
    subject: LocateWindowName


class LocateWindowsReply(TypedDict):
    """The plan as ``LocatePlan.to_json`` writes it (``gbdraw locate`` writes the same as ``PREFIX.plan.json``),
    and the window pairs of its drawing (``locate_drawing``)."""

    plan: LocatePlanJson
    pairs: list[LocateWindowPairReply]


def _window_name(drawing: LocateDrawing, place: LocateWindowPlace) -> LocateWindowName:
    row, index, window_index = place
    record = drawing.rows[row].records[index]
    window = record.windows[window_index]
    return {"fileKey": record.file_key, "recordKey": record.record_key, "start": window.start, "end": window.end}


def _file(value: LocateFileInput) -> LocateInputFile:
    return LocateInputFile(str(value["key"]), value["gbk"], value["gff"], value["fasta"])


def _search(payload: LocateSearchInput) -> LocateSearch:
    return prepare_locate_search(
        _file(payload["reference"]),
        [_file(target) for target in payload["targets"]],
        record=payload["recordId"],
        record_key=payload["recordKey"],
        windows=[
            RecordWindow(int(window["start"]), int(window["end"]), bool(window["reverseComplement"]), column)
            for column, window in enumerate(payload["windows"], start=1)
        ],
        program=payload["program"],
    )


def locate_search_plan_json(payload: LocateSearchInput) -> LocateSearchPlanReply:
    """One LOSAT job per target File for the page to run."""

    plan = _search(payload).plan
    return {
        "program": plan.program,
        "cacheArgs": list(plan.cache_args),
        "jobs": [
            {"fileKey": job.file_key, "queryFasta": job.query_fasta, "subjectFasta": job.subject_fasta, "rawKey": job.raw_key}
            for job in plan.jobs
        ],
    }


def locate_windows_json(payload: LocateWindowsInput) -> LocateWindowsReply:
    """The plan of the jobs' rows with the dialog's thresholds and window rules."""

    search = _search(payload)
    values = payload["options"]
    options = replace(
        LocateOptions.for_program(payload["program"]),
        thresholds=ComparisonThresholds(
            float(values["evalue"]), float(values["bitscore"]), float(values["identity"]), int(values["alignmentLength"])
        ),
        flank=int(values["flank"]),
        merge_distance=int(values["mergeDistance"]),
        max_windows=int(values["maxWindows"]),
    )
    outcomes = [
        LocateSearchResult(str(outcome["fileKey"]), str(outcome["text"]), outcome["error"])
        for outcome in payload["outcomes"]
    ]
    plan = finish_locate(search, locate_results(search, outcomes), options)
    drawing = locate_drawing(plan)
    return {
        "plan": plan.to_json(),
        "pairs": [
            {"query": _window_name(drawing, upper), "subject": _window_name(drawing, lower)}
            for upper, lower in drawing.pairs
        ],
    }


__all__ = [
    "LocateFileInput",
    "LocateOptionsInput",
    "LocateOutcomeInput",
    "LocateSearchInput",
    "LocateSearchPlanReply",
    "LocateWindowInput",
    "LocateWindowName",
    "LocateWindowPairReply",
    "LocateWindowsInput",
    "LocateWindowsReply",
    "locate_search_plan_json",
    "locate_windows_json",
]
