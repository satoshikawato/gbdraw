"""Worker helpers of Find in other files (PD-OI-123): the reference, the search plan, the plan of its rows, and the
new drawing, as JSON.

The page sends the shown Linear drawing's draft request (its canonical
request, as the next Generate would send it, with its staged resources), the
reference record by its key in that request with the reference windows, and
the Files to search by their keys (each record of the request is one File of
the drawing, a Linear File card). Python reads the Files the request names:

- :func:`locate_web_reference` checks the reference windows and counts their
  genes (the CDS that lie wholly inside a window) for the dialog;
- :func:`locate_web_search_plan` gives one LOSAT job per File, which the page
  runs;
- :func:`locate_web_plan` turns each job's rows (or its error) into the plan
  with the dialog's thresholds;
- :func:`locate_web_drawing` takes the plan with the windows the user keeps
  and their edits, and derives the new Linear drawing from the drawing it was
  found from, as New drawing from a region does.

Each calls the locate owner (:mod:`gbdraw.comparisons.locate`) as ``gbdraw
locate`` does, so the Web has no copy of the search plan, of the window rules,
or of the drawing's rows and pairs.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Literal, TypedDict

from gbdraw.api.locate_drawing import LocateFile, derive_locate_request
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import (
    FastaInputSource,
    GenBankInputSource,
    GffFastaInputSource,
    LinearDiagramRequest,
    RecordInputSource,
)
from gbdraw.comparisons.locate import (
    LocateDrawing,
    LocateInputFile,
    LocateOptions,
    LocatePlan,
    LocatePlanJson,
    LocateProgramName,
    LocateSearch,
    LocateSearchResult,
    LocateWindowChoice,
    LocateWindowPlace,
    LocateWindowRefusal,
    LocateWindowRefusalReason,
    choose_locate_windows,
    finish_locate,
    locate_drawing,
    locate_results,
    prepare_locate_search,
)
from gbdraw.contract_shapes import Shape, ref, typed_dict_shape
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import RecordWindow
from gbdraw.mode_profiles import ComparisonThresholds
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.web_support.region_drawing import DerivedRegionDrawing, derived_drawing_reply


class LocateWindowInput(TypedDict):
    """One reference window, 1-based and inclusive; ``start > end`` crosses the origin of a circular record."""

    start: int
    end: int
    reverseComplement: bool


class LocateReferenceInput(TypedDict):
    """The reference record by its key in the drawing's request, and its windows in column order."""

    recordKey: str
    windows: list[LocateWindowInput]


class LocateReferenceReady(TypedDict):
    """The reference windows can be searched: the File and record they are of, and their genes (``genes``, the
    CDS that lie wholly inside a window) and the CDS a window edge cuts (``cutGenes``, by name), which are not
    searched."""

    status: Literal["ready"]
    fileKey: str
    recordId: str
    genes: int
    cutGenes: list[str]


class LocateWindowRefused(TypedDict):
    """A window the plan cannot take: the rule it breaks, and its span."""

    status: Literal["refused"]
    reason: LocateWindowRefusalReason
    start: int
    end: int


LocateReferenceReply = LocateReferenceReady | LocateWindowRefused


class LocateSearchInput(TypedDict):
    """The reference, the Files to search by their keys (the request's record keys, in the drawing's order),
    and the program."""

    reference: LocateReferenceInput
    targets: list[str]
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
    """The dialog's thresholds and window rules."""

    evalue: float
    bitscore: float
    identity: float
    alignmentLength: int
    flank: int
    mergeDistance: int
    maxWindows: int


class LocatePlanInput(LocateSearchInput):
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


class LocatePlanReply(TypedDict):
    """The plan as ``LocatePlan.to_json`` writes it (``gbdraw locate`` writes the same as ``PREFIX.plan.json``),
    the window pairs of its drawing (``locate_drawing``), and the name the dialog shows for each reference gene
    by protein ID (its gene name, else its locus tag, else its protein ID)."""

    plan: LocatePlanJson
    pairs: list[LocateWindowPairReply]
    geneNames: dict[str, str]


class LocateWindowChoiceInput(TypedDict):
    """A window the user keeps: its File, its column among the File's windows (its identity in the plan), and
    its span and orientation (the search's, or the user's edit; ``start > end`` across the origin)."""

    fileKey: str
    column: int
    start: int
    end: int
    reverseComplement: bool


class LocateDrawingInput(TypedDict):
    """The plan the search gave, the windows the user keeps, the name of each File's row by File key, and
    whether hand-set sizes whose Auto differs return to Auto."""

    plan: LocatePlanJson
    choices: list[LocateWindowChoiceInput]
    labels: dict[str, str]
    adaptSizes: bool


class LocateDrawingReady(TypedDict):
    """The plan with the windows kept, its window pairs, and the new drawing."""

    status: Literal["derived"]
    plan: LocatePlanJson
    pairs: list[LocateWindowPairReply]
    drawing: DerivedRegionDrawing


class LocateEditRefused(LocateWindowRefused):
    """A kept window whose edit the plan cannot take, by its File and column."""

    fileKey: str
    column: int


LocateDrawingReply = LocateDrawingReady | LocateEditRefused


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL"})


def _input_file(key: str, source: RecordInputSource) -> LocateInputFile:
    if isinstance(source, GenBankInputSource):
        return LocateInputFile(key, str(source.path), None, None)
    if isinstance(source, GffFastaInputSource):
        return LocateInputFile(key, None, str(source.gff_path), str(source.fasta_path))
    if isinstance(source, FastaInputSource):
        return LocateInputFile(key, None, None, str(source.fasta_path))
    raise _protocol_error("Find in other files reads files only.")


@dataclass(frozen=True)
class _Drawing:
    """The source drawing's request and its Files: each request record is one File, by its record key."""

    request: LinearDiagramRequest
    files: dict[str, RecordInputSource]
    # The reference record's File, its ``#index`` selector and its key, by any key the drawing gives it.
    records: dict[str, tuple[str, str, str]]
    # The keys the drawn records of each File have, by File key and 1-based record index.
    record_keys: dict[str, dict[int, str]]


def _drawing(request_payload: Mapping[str, object], resource_paths: Mapping[str, str], output_directory: str) -> _Drawing:
    request = decode_canonical_request(request_payload, resource_paths=resource_paths, output_directory=output_directory)
    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError("Find in other files searches from a Linear drawing.", diagnostic={"code": "INPUT_INVALID"})
    files: dict[str, RecordInputSource] = {}
    for record in request.records:
        if record.record_key is None or record.record_key in files:
            raise _protocol_error("Each File of the drawing needs its own record key.")
        files[record.record_key] = record.source
    keys = list(files)
    records: dict[str, tuple[str, str, str]] = {}
    record_keys: dict[str, dict[int, str]] = {}
    for unit in resolve_request_records(request).provenance:
        file_key = keys[unit.input_index]
        index = unit.source_record_index + 1
        found = (file_key, f"#{index}", unit.record_identity_key)
        records.setdefault(unit.record_identity_key, found)
        records.setdefault(unit.record_key, found)
        record_keys.setdefault(file_key, {}).setdefault(index, unit.record_identity_key)
    return _Drawing(request, files, records, record_keys)


def _search(drawing: _Drawing, payload: LocateSearchInput) -> LocateSearch:
    reference = payload["reference"]
    found = drawing.records.get(str(reference["recordKey"]))
    if found is None:
        raise ValidationError(
            f"The drawing has no record {reference['recordKey']!r}.",
            diagnostic={"code": "INPUT_INVALID", "field": "record_selector"},
        )
    file_key, selector, record_key = found
    targets = [str(key) for key in payload["targets"]]
    if any(key not in drawing.files for key in targets) or file_key in targets or len(set(targets)) != len(targets):
        raise _protocol_error("Find in other files searches other Files of the drawing, each once.")
    return prepare_locate_search(
        _input_file(file_key, drawing.files[file_key]),
        [_input_file(key, drawing.files[key]) for key in targets],
        record=selector,
        record_key=record_key,
        windows=[
            RecordWindow(int(window["start"]), int(window["end"]), bool(window["reverseComplement"]), column)
            for column, window in enumerate(reference["windows"], start=1)
        ],
        program=payload["program"],
        target_record_keys={key: drawing.record_keys.get(key, {}) for key in targets},
    )


def _refused(refusal: LocateWindowRefusal) -> LocateWindowRefused:
    return {"status": "refused", "reason": refusal.reason, "start": refusal.start, "end": refusal.end}


def locate_web_reference(
    request_payload: Mapping[str, object], resource_paths: Mapping[str, str], output_directory: str,
    reference: LocateReferenceInput,
) -> LocateReferenceReply:
    """The genes of the reference windows, or the rule a window breaks; no File is searched."""

    if not reference["windows"]:
        raise _protocol_error("The reference needs at least one window.")
    try:
        search = _search(_drawing(request_payload, resource_paths, output_directory), {
            "reference": reference, "targets": [], "program": "tlosatn",
        })
    except LocateWindowRefusal as refusal:
        return _refused(refusal)
    located = search.query.reference
    return {
        "status": "ready",
        "fileKey": located.file_key,
        "recordId": located.record_id,
        "genes": len(located.genes),
        "cutGenes": [cds.name for cds in located.cut_genes],
    }


def locate_web_search_plan(
    request_payload: Mapping[str, object], resource_paths: Mapping[str, str], output_directory: str,
    payload: LocateSearchInput,
) -> LocateSearchPlanReply:
    """One LOSAT job per File for the page to run."""

    plan = _search(_drawing(request_payload, resource_paths, output_directory), payload).plan
    return {
        "program": plan.program,
        "cacheArgs": list(plan.cache_args),
        "jobs": [
            {"fileKey": job.file_key, "queryFasta": job.query_fasta, "subjectFasta": job.subject_fasta, "rawKey": job.raw_key}
            for job in plan.jobs
        ],
    }


def _window_name(drawing: LocateDrawing, place: LocateWindowPlace) -> LocateWindowName:
    row, index, window_index = place
    record = drawing.rows[row].records[index]
    window = record.windows[window_index]
    return {"fileKey": record.file_key, "recordKey": record.record_key, "start": window.start, "end": window.end}


def _pairs(plan: LocatePlan) -> list[LocateWindowPairReply]:
    drawing = locate_drawing(plan)
    return [{"query": _window_name(drawing, upper), "subject": _window_name(drawing, lower)} for upper, lower in drawing.pairs]


def locate_web_plan(
    request_payload: Mapping[str, object], resource_paths: Mapping[str, str], output_directory: str,
    payload: LocatePlanInput,
) -> LocatePlanReply:
    """The plan of the jobs' rows with the dialog's thresholds and window rules."""

    search = _search(_drawing(request_payload, resource_paths, output_directory), payload)
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
    return {
        "plan": plan.to_json(),
        "pairs": _pairs(plan),
        "geneNames": {gene.protein_id: gene.cds.name for gene in plan.reference.genes},
    }


def _choices(values: Sequence[LocateWindowChoiceInput]) -> list[LocateWindowChoice]:
    try:
        return [
            LocateWindowChoice(
                str(value["fileKey"]), int(value["column"]), int(value["start"]), int(value["end"]),
                bool(value["reverseComplement"]),
            )
            for value in values
        ]
    except (KeyError, TypeError, ValueError) as exc:
        raise _protocol_error("A kept window needs its File, column, start, end and orientation.") from exc


def locate_web_drawing(
    request_payload: Mapping[str, object], resource_paths: Mapping[str, str], output_directory: str,
    payload: LocateDrawingInput,
) -> LocateDrawingReply:
    """The new Linear drawing of the plan with the windows the user keeps, derived from the drawing it was found
    from; the derived request names the source drawing's resources (its Files are the new drawing's inputs)."""

    if not isinstance(payload["adaptSizes"], bool):
        raise _protocol_error("adaptSizes must be true or false.")
    choices = _choices(payload["choices"])
    try:
        plan = choose_locate_windows(LocatePlan.from_json(payload["plan"]), choices)
    except LocateWindowRefusal as refusal:
        assert refusal.place is not None
        file_key, column = refusal.place
        return {**_refused(refusal), "fileKey": file_key, "column": column}
    drawing = _drawing(request_payload, resource_paths, output_directory)
    labels = {str(key): str(value) for key, value in payload["labels"].items()}
    files = {
        key: LocateFile(source=source, label=labels.get(key, key))
        for key, source in drawing.files.items()
    }
    derived = derive_locate_request(drawing.request, plan, files, adapt_sizes=payload["adaptSizes"])
    return {
        "status": "derived",
        "plan": plan.to_json(),
        "pairs": _pairs(plan),
        "drawing": derived_drawing_reply(
            derived.request, "linear", derived.adaptation, derived.not_carried, resource_paths
        ),
    }


def describe_locate() -> dict[str, dict[str, Shape]]:
    """The inputs and replies of the helpers above as named shapes of :mod:`gbdraw.contract_shapes`
    (tools/generate_web_contracts.py writes them as the Web's typedefs)."""

    from typing import is_typeddict

    types_: dict[str, Shape] = {}

    def special(hint: object) -> Shape | None:
        if isinstance(hint, type) and is_typeddict(hint):
            if hint.__name__ not in types_:
                types_[hint.__name__] = {"kind": "unknown"}  # reserves its place before its members'
                types_[hint.__name__] = typed_dict_shape(hint, special)
            return ref(hint.__name__)
        return None

    for root in (
        LocateReferenceInput, LocateReferenceReady, LocateWindowRefused, LocateSearchInput, LocateSearchPlanReply,
        LocatePlanInput, LocatePlanReply, LocateDrawingInput, LocateDrawingReady, LocateEditRefused,
    ):
        special(root)
    return {"types": types_}


def _json_helper(function, request_json: str, resource_paths_json: str, output_directory: str, payload_json: str) -> str:  # type: ignore[no-untyped-def]
    return json.dumps(function(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory),
        json.loads(str(payload_json)),
    ))


def locate_web_reference_json(request_json: str, resource_paths_json: str, output_directory: str, payload_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return _json_helper(locate_web_reference, request_json, resource_paths_json, output_directory, payload_json)


def locate_web_search_plan_json(request_json: str, resource_paths_json: str, output_directory: str, payload_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return _json_helper(locate_web_search_plan, request_json, resource_paths_json, output_directory, payload_json)


def locate_web_plan_json(request_json: str, resource_paths_json: str, output_directory: str, payload_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return _json_helper(locate_web_plan, request_json, resource_paths_json, output_directory, payload_json)


def locate_web_drawing_json(request_json: str, resource_paths_json: str, output_directory: str, payload_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return _json_helper(locate_web_drawing, request_json, resource_paths_json, output_directory, payload_json)


__all__ = [
    "LocateDrawingInput",
    "LocateDrawingReady",
    "LocateDrawingReply",
    "LocateEditRefused",
    "LocateJobReply",
    "LocateOptionsInput",
    "LocateOutcomeInput",
    "LocatePlanInput",
    "LocatePlanReply",
    "LocateReferenceInput",
    "LocateReferenceReady",
    "LocateReferenceReply",
    "LocateSearchInput",
    "LocateSearchPlanReply",
    "LocateWindowChoiceInput",
    "LocateWindowInput",
    "LocateWindowName",
    "LocateWindowPairReply",
    "LocateWindowRefused",
    "describe_locate",
    "locate_web_drawing",
    "locate_web_drawing_json",
    "locate_web_plan",
    "locate_web_plan_json",
    "locate_web_reference",
    "locate_web_reference_json",
    "locate_web_search_plan",
    "locate_web_search_plan_json",
]
