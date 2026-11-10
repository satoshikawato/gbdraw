"""Find in other files: the new Linear drawing of a locate plan.

The drawing's first row holds the reference windows, labelled with the
reference File's name; each target File with a selected window has one row
with its windows in the plan's order, labelled with the File's name
(:func:`gbdraw.comparisons.locate.locate_drawing`). Windows of adjacent rows
that share a reference gene are compared by LOSAT. The source drawing's
settings and edits that apply to the new records are carried as New drawing
from a region carries them (:func:`gbdraw.api.region_drawing.derive_from_source`).
The records table and comparisons table that ``gbdraw locate`` writes read
back to the same records and pairs.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

from gbdraw.api.options import (
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LosatRuntimeOptions,
    LosatSearchOptions,
)
from gbdraw.api.region_drawing import (
    DerivedRecords,
    NotCarried,
    SourceRegion,
    derive_from_source,
)
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import (
    DiagramRequest,
    LinearDiagramRequest,
    RecordInput,
    RecordInputSource,
    RecordPresentation,
)
from gbdraw.auto_sizes import SettingAdaptation
from gbdraw.comparisons.locate import LocateDrawing, LocateDrawnRecord, LocatePlan, locate_drawing
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import RecordWindow
from gbdraw.session import MaterializedSession, SessionDrawingSpec, session_to_request

_COMPARISON_FILES = ("blast_files", "linear_comparisons", "protein_comparisons", "orthogroups", "collinearity_blocks")


@dataclass(frozen=True)
class LocateFile:
    """One File of a locate plan, by its ``file_key``: where its records are read, and the name of its row."""

    source: RecordInputSource
    label: str


@dataclass(frozen=True)
class LocateRequest:
    """The request of the new drawing, the sizes that returned to Auto, and what was not carried."""

    request: LinearDiagramRequest
    adaptation: SettingAdaptation
    not_carried: tuple[NotCarried, ...]


@dataclass(frozen=True)
class LocateNewDrawing:
    """The new drawing, the sizes that returned to Auto, and what was not carried."""

    drawing: SessionDrawingSpec
    adaptation: SettingAdaptation
    not_carried: tuple[NotCarried, ...]


def _record_input(record: LocateDrawnRecord, row: int, label: str | None, files: Mapping[str, LocateFile]) -> RecordInput:
    """One record of the drawing as the records table gives it: a whole record placed and oriented by its
    presentation, one window placed by the presentation and oriented on the window (``RecordInput.oriented``, the
    form the Web's File cards write back), several windows by their own orientation and column."""

    file = files.get(record.file_key)
    if file is None:
        raise ValidationError(f"Locate has no File {record.file_key!r}.", diagnostic={"code": "INPUT_INVALID"})
    first = record.windows[0]
    several = len(record.windows) > 1
    drawn = RecordInput(
        source=file.source,
        selector=parse_record_selector(record.record_id),
        windows=(
            record.windows if several
            else () if record.whole
            else (RecordWindow(first.start, first.end),)
        ),
        presentation=RecordPresentation(label=label, grid_row=row, grid_column=None if several else first.column),
        record_key=record.record_key,
    )
    return drawn if several else drawn.oriented(first.reverse_complement)


def locate_record_inputs(drawing: LocateDrawing, files: Mapping[str, LocateFile]) -> tuple[RecordInput, ...]:
    """The record inputs of ``drawing`` in row order; the first record of a row carries the File's name."""

    return tuple(
        _record_input(record, row_index, files[row.file_key].label if index == 0 else None, files)
        for row_index, row in enumerate(drawing.rows, start=1)
        for index, record in enumerate(row.records)
    )


def locate_unit_pairs(drawing: LocateDrawing) -> tuple[tuple[int, int], ...]:
    """The drawn-unit pairs LOSAT compares: the window pairs of ``drawing``, each window one drawn unit in row,
    record and window order (``locate_record_inputs``)."""

    units = {
        (row_index, record_index, window_index): unit
        for unit, (row_index, record_index, window_index) in enumerate(
            (row_index, record_index, window_index)
            for row_index, row in enumerate(drawing.rows)
            for record_index, record in enumerate(row.records)
            for window_index in range(len(record.windows))
        )
    }
    return tuple((units[upper], units[lower]) for upper, lower in drawing.pairs)


def _locate_comparisons(
    options: LinearDiagramOptions, dropped: list[NotCarried], pairs: Sequence[tuple[int, int]]
) -> dict[str, object]:
    """LOSATP links between the windows of adjacent rows that share a reference gene; the source's comparisons
    are not carried."""

    search = options.losat_search
    if any(getattr(options, name) for name in _COMPARISON_FILES) or options.comparison_table_file is not None:
        dropped.append(NotCarried("COMPARISON_TABLES", "linear"))
    if search is not None and not (search.program == "losatp" and search.losatp_mode == "pairwise"):
        dropped.append(NotCarried("LINEAR_COMPARISONS", "linear"))
    base = (
        search
        if search is not None and search.program == "losatp"
        else LosatSearchOptions(
            "losatp", runtime=search.runtime if search is not None else LosatRuntimeOptions(), losatp_mode="pairwise"
        )
    )
    return {
        **dict.fromkeys((*_COMPARISON_FILES, "comparison_table_file")),
        "losat_search": replace(base, losatp_mode="pairwise", pairs=tuple(pairs)),
    }


def derive_locate_request(
    source: DiagramRequest,
    plan: LocatePlan,
    files: Mapping[str, LocateFile],
    *,
    adapt_sizes: bool = True,
) -> LocateRequest:
    """The request of the new Linear drawing of ``plan`` from the request of its source Linear drawing.

    ``files`` gives each File of the plan by its ``file_key``. A record whose
    ``record_key`` is a record of the source drawing carries that record's
    edits inside its windows; the plan's other records are new.
    """

    if not isinstance(source, LinearDiagramRequest):
        raise ValidationError("Find in other files derives from a Linear drawing.", diagnostic={"code": "INPUT_INVALID"})
    drawing = locate_drawing(plan)
    records = locate_record_inputs(drawing, files)
    drawn = [record for row in drawing.rows for record in row.records]
    collection = resolve_request_records(source)
    sources: dict[str, int] = {}
    for unit in collection.provenance:
        sources.setdefault(unit.record_identity_key, unit.resolved_index)
    regions: list[SourceRegion] = []
    first_unit = 0
    for index, record in enumerate(drawn):
        resolved = sources.get(record.record_key)
        if resolved is not None:
            item = collection.provenance[resolved]
            regions.append(SourceRegion(
                item=item,
                index=index,
                spans=tuple((window.start, window.end) for window in record.windows),
                whole=record.whole,
                first_unit=first_unit,
                length=record.record_length,
            ))
        first_unit += len(record.windows)
    pairs = locate_unit_pairs(drawing)
    derived = derive_from_source(
        source,
        collection,
        DerivedRecords(
            mode="linear",
            records=records,
            record_ids=tuple(record.record_id for record in drawn),
            lengths=tuple(
                window.end - window.start + 1 if window.start <= window.end
                else record.record_length - window.start + 1 + window.end
                for record in drawn
                for window in record.windows
            ),
            regions=tuple(regions),
            # The rows are grid rows; the source's record gap spaces the records of a row.
            layout=LinearMultiRecordOptions(
                record_gap_px=source.layout.record_gap_px if source.layout is not None else LinearMultiRecordOptions().record_gap_px
            ),
        ),
        adapt_sizes=adapt_sizes,
        linear_comparisons=lambda options, dropped: _locate_comparisons(options, dropped, pairs),
    )
    assert isinstance(derived.request, LinearDiagramRequest)
    return LocateRequest(request=derived.request, adaptation=derived.adaptation, not_carried=derived.not_carried)


def derive_locate_drawing(
    materialized: MaterializedSession,
    plan: LocatePlan,
    files: Mapping[str, LocateFile],
    *,
    drawing: str | None = None,
    name: str | None = None,
    adapt_sizes: bool = True,
) -> LocateNewDrawing:
    """The new drawing of ``plan`` from a Linear drawing of ``materialized`` (by ID or name; without it, the only
    drawing). Use the result while ``materialized`` is active: its request reads the materialized files."""

    derived = derive_locate_request(
        session_to_request(materialized, drawing=drawing), plan, files, adapt_sizes=adapt_sizes
    )
    return LocateNewDrawing(
        drawing=SessionDrawingSpec(request=derived.request, mode="linear", name=name),
        adaptation=derived.adaptation,
        not_carried=derived.not_carried,
    )


__all__ = [
    "LocateFile",
    "LocateNewDrawing",
    "LocateRequest",
    "derive_locate_drawing",
    "derive_locate_request",
    "locate_record_inputs",
    "locate_unit_pairs",
]
