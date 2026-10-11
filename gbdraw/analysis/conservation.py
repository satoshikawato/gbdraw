"""Circular conservation-ring input loading and normalization."""

from __future__ import annotations

import io
import logging
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Literal, Sequence, cast

import pandas as pd
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.core.color import normalize_hex_color, tint_color
from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import resolve_color_to_hex
from gbdraw.io.comparison_sequences import comparison_file_stem
from gbdraw.io.comparisons import (
    COMPARISON_COLUMNS,
    filter_comparison_dataframe,
    normalize_comparison_dataframe,
    read_comparison_table,
)

if TYPE_CHECKING:
    from gbdraw.configurators import BlastMatchConfigurator

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ConservationSearchResult:
    """Raw outfmt 6 rows of one similarity ring that a LOSAT search produced.

    The planner sets these in place of the search intent. ``name`` is the raw
    TSV filename (Session resource and ``--losat_output_dir`` file). A Session
    stores them as ``conservationBlastFiles``, so replay reads them as files.
    """

    name: str
    text: str

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValidationError(
                "ConservationSearchResult.name must be a non-empty string.",
                diagnostic={"code": "INPUT_INVALID", "field": "conservation_search_results"},
            )
        if not isinstance(self.text, str):
            raise ValidationError(
                "ConservationSearchResult.text must be a string.",
                diagnostic={"code": "INPUT_INVALID", "field": "conservation_search_results"},
            )


ConservationReferenceSide = Literal["query", "subject"]
ConservationReferenceMode = Literal["query", "subject", "auto"]

NORMALIZED_CONSERVATION_COLUMNS = (
    "source_index",
    "source_hit_index",
    "track_index",
    "track_label",
    "reference_side",
    "reference_match_key",
    "reference_record_id",
    "track_color",
    "query",
    "subject",
    "qstart",
    "qend",
    "sstart",
    "send",
    "start",
    "end",
    "draw_start",
    "draw_end",
    "identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "evalue",
    "bitscore",
    "orientation",
    "full_reference",
)


@dataclass(frozen=True)
class ConservationSource:
    source_index: int
    label: str
    color: str | None
    dataframe: DataFrame | None
    path: str | None = None
    skipped: bool = False
    skip_reason: str | None = None


@dataclass(frozen=True)
class ConservationTrack:
    source_index: int
    track_index: int
    track_label: str
    track_color: str | None
    reference_side: ConservationReferenceSide | None
    hits: DataFrame


@dataclass(frozen=True)
class ConservationLoadResult:
    sources: tuple[ConservationSource, ...]
    skipped_sources: tuple[ConservationSource, ...]


def normalize_conservation_reference(value: object | None) -> ConservationReferenceMode:
    normalized = str(value or "auto").strip().lower()
    if normalized not in {"query", "subject", "auto"}:
        raise ValidationError("conservation_reference must be one of: query, subject, auto")
    return cast(ConservationReferenceMode, normalized)


def empty_normalized_conservation_hits() -> DataFrame:
    return DataFrame(columns=NORMALIZED_CONSERVATION_COLUMNS)


def _default_label(source_index: int, path: "str | ConservationSearchResult | None") -> str:
    """The file name without its last extension (D-03), as for a FASTA ring and in the Web."""

    if isinstance(path, ConservationSearchResult):
        path = path.name
    if path:
        basename = os.path.basename(str(path))
        if basename:
            return comparison_file_stem(basename).strip() or basename
    return f"Conservation {int(source_index) + 1}"


def normalize_conservation_color(value: object | None) -> str | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return normalize_hex_color(resolve_color_to_hex(text))
    except Exception as exc:
        raise ValidationError(f"Invalid conservation color {text!r}. Use an SVG color name or #RRGGBB.") from exc


def _default_conservation_color(value: object | None, leaf: str) -> str:
    """A configured default ring color: a config-accepted value this domain rejects names its leaf."""

    path = f"objects.conservation.{leaf}"
    try:
        color = normalize_conservation_color(value)
    except ValidationError as exc:
        raise ValidationError(
            f"{path} {str(value).strip()!r} cannot color a conservation ring. Use an SVG color name or #RRGGBB.",
            diagnostic={"code": "INPUT_INVALID", "reason": "COLOR", "configPath": path},
        ) from exc
    if color is None:
        raise ValidationError(
            f"{path} must be a color.",
            diagnostic={"code": "INPUT_INVALID", "reason": "COLOR", "configPath": path},
        )
    return color


def conservation_track_gradient_colors(
    track_color: object | None,
    *,
    default_min_color: str,
    default_max_color: str,
) -> tuple[str, str]:
    normalized_track_color = normalize_conservation_color(track_color)
    if normalized_track_color is None:
        return (
            _default_conservation_color(default_min_color, "min_color"),
            _default_conservation_color(default_max_color, "max_color"),
        )
    return tint_color(normalized_track_color), normalized_track_color


def _filter_normalized_dataframe(dataframe: DataFrame, blast_config: BlastMatchConfigurator) -> DataFrame:
    df = dataframe.copy()
    # Keep the original source-row identity through filtering and paint-order sorting.
    df["source_hit_index"] = range(len(df))
    filtered = filter_comparison_dataframe(df, blast_config)
    return filtered.loc[:, [*COMPARISON_COLUMNS, "source_hit_index"]].reset_index(drop=True)


def _load_conservation_file(
    path: "str | ConservationSearchResult",
    blast_config: BlastMatchConfigurator,
) -> tuple[DataFrame | None, str | None]:
    try:
        if isinstance(path, ConservationSearchResult):
            # Resolved LOSAT rows take the file path's reader (byte-identical rings).
            table = read_comparison_table(io.StringIO(path.text), label=path.name)
        else:
            table = read_comparison_table(path)
        return _filter_normalized_dataframe(table, blast_config), None
    except ValidationError as exc:
        return None, str(exc)


def load_conservation_sources(
    *,
    blast_config: object,
    conservation_files: "Sequence[str | ConservationSearchResult] | None" = None,
    conservation_dataframes: Sequence[DataFrame] | None = None,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
) -> ConservationLoadResult:
    """Load conservation sources while preserving user-visible source indexes."""

    files = list(conservation_files or [])
    dataframes = list(conservation_dataframes or [])
    logical_count = max(len(files), len(dataframes))
    if logical_count == 0:
        return ConservationLoadResult(sources=(), skipped_sources=())
    if labels is not None and len(labels) < logical_count:
        raise ValidationError(
            f"Expected at least {logical_count} conservation label(s); got {len(labels)}."
        )
    if colors is not None and len(colors) < logical_count:
        raise ValidationError(
            f"Expected at least {logical_count} conservation color(s); got {len(colors)}."
        )

    sources: list[ConservationSource] = []
    skipped: list[ConservationSource] = []

    for source_index in range(logical_count):
        path = files[source_index] if source_index < len(files) else None
        label = (
            str(labels[source_index]).strip()
            if labels is not None
            else _default_label(source_index, path)
        )
        if not label:
            label = _default_label(source_index, path)
        color = normalize_conservation_color(colors[source_index]) if colors is not None else None

        valid_frames: list[DataFrame] = []
        skip_reasons: list[str] = []

        if path is not None:
            frame, reason = _load_conservation_file(
                path if isinstance(path, ConservationSearchResult) else str(path),
                cast("BlastMatchConfigurator", blast_config),
            )
            if frame is not None:
                valid_frames.append(frame)
            elif reason:
                logger.warning("WARNING: Skipping conservation source %s file: %s", source_index, reason)
                skip_reasons.append(reason)

        if source_index < len(dataframes):
            try:
                valid_frames.append(
                    _filter_normalized_dataframe(
                        normalize_comparison_dataframe(dataframes[source_index]),
                        cast("BlastMatchConfigurator", blast_config),
                    )
                )
            except ValidationError as exc:
                reason = f"error parsing conservation dataframe {source_index}: {exc}"
                logger.warning("WARNING: %s", reason)
                skip_reasons.append(reason)

        if valid_frames:
            merged = (
                pd.concat(valid_frames, ignore_index=True)
                if len(valid_frames) > 1
                else valid_frames[0].copy()
            )
            # DataFrame and file inputs can be merged for one logical source. A
            # single monotonically increasing index keeps IDs unique in that case.
            merged["source_hit_index"] = range(len(merged))
            sources.append(
                ConservationSource(
                    source_index=source_index,
                    label=label,
                    color=color,
                    dataframe=merged.reset_index(drop=True),
                    path=str(path) if path is not None else None,
                )
            )
            continue

        reason = "; ".join(skip_reasons) if skip_reasons else "no valid conservation data"
        skipped_source = ConservationSource(
            source_index=source_index,
            label=label,
            color=color,
            dataframe=None,
            path=str(path) if path is not None else None,
            skipped=True,
            skip_reason=reason,
        )
        sources.append(skipped_source)
        skipped.append(skipped_source)

    return ConservationLoadResult(
        sources=tuple(sources),
        skipped_sources=tuple(skipped),
    )


def _record_match_keys(record: SeqRecord) -> set[str]:
    keys: set[str] = set()

    def add(value: object) -> None:
        text = str(value or "").strip()
        if not text:
            return
        keys.add(text)
        first_token = text.split()[0]
        if first_token:
            keys.add(first_token)

    add(record.id)
    add(record.name)
    annotations = getattr(record, "annotations", {}) or {}
    accessions = annotations.get("accessions")
    if isinstance(accessions, Iterable) and not isinstance(accessions, (str, bytes)):
        for accession in accessions:
            add(accession)
    else:
        add(accessions)
    add(annotations.get("accession"))
    sequence_version = annotations.get("sequence_version")
    if sequence_version is not None and accessions:
        first_accession = next(iter(accessions)) if not isinstance(accessions, str) else accessions
        add(f"{first_accession}.{sequence_version}")
    return keys


def _reference_key_map(records: Sequence[SeqRecord]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for record in records:
        for key in _record_match_keys(record):
            mapping.setdefault(key, str(record.id))
    return mapping


def resolve_conservation_reference_side(
    source: ConservationSource,
    displayed_records: Sequence[SeqRecord],
    reference: ConservationReferenceMode | str = "auto",
) -> ConservationReferenceSide | None:
    mode = normalize_conservation_reference(reference)
    if mode in {"query", "subject"}:
        return mode

    df = source.dataframe
    if df is None or df.empty:
        return None

    reference_keys = set(_reference_key_map(displayed_records))
    query_matches = df["query"].astype(str).isin(reference_keys)
    subject_matches = df["subject"].astype(str).isin(reference_keys)
    mentioned = query_matches | subject_matches
    if not bool(mentioned.any()):
        raise ValidationError(
            f"Could not resolve conservation source {source.source_index} reference side. "
            "Set --conservation_reference to query or subject."
        )

    query_any = bool(query_matches[mentioned].any())
    subject_any = bool(subject_matches[mentioned].any())
    if query_any and not subject_any:
        return "query"
    if subject_any and not query_any:
        return "subject"
    raise ValidationError(
        f"Conservation source {source.source_index} matches displayed records on both BLAST sides. "
        "Set --conservation_reference to query or subject."
    )


def _row_float(row: object, name: str) -> float:
    value = getattr(row, name)
    return float(value)


def _row_text(row: object, name: str) -> str:
    value = getattr(row, name)
    return str(value)


def _normalize_source_hits_for_record(
    source: ConservationSource,
    *,
    track_index: int,
    reference_side: ConservationReferenceSide | None,
    record: SeqRecord,
) -> DataFrame:
    df = source.dataframe
    if df is None or df.empty or reference_side is None:
        return empty_normalized_conservation_hits()

    record_len = len(record.seq)
    if record_len <= 0:
        return empty_normalized_conservation_hits()

    start_column, end_column, id_column = (
        ("qstart", "qend", "query")
        if reference_side == "query"
        else ("sstart", "send", "subject")
    )
    record_keys = _record_match_keys(record)
    matched = df[df[id_column].astype(str).isin(record_keys)]
    rows: list[dict[str, object]] = []
    dropped = 0

    for row in matched.itertuples(index=False):
        try:
            raw_start = _row_float(row, start_column)
            raw_end = _row_float(row, end_column)
            if pd.isna(raw_start) or pd.isna(raw_end):
                dropped += 1
                continue
            start = min(raw_start, raw_end)
            end = max(raw_start, raw_end)
            draw_start = max(0.0, start - 1.0)
            draw_end = min(float(record_len), end)
            if draw_start >= draw_end:
                dropped += 1
                continue
            orientation = "forward" if raw_start <= raw_end else "reverse"
            rows.append(
                {
                    "source_index": int(source.source_index),
                    "source_hit_index": int(getattr(row, "source_hit_index")),
                    "track_index": int(track_index),
                    "track_label": source.label,
                    "reference_side": reference_side,
                    "reference_match_key": _row_text(row, id_column),
                    "reference_record_id": str(record.id),
                    "track_color": source.color or "",
                    "query": _row_text(row, "query"),
                    "subject": _row_text(row, "subject"),
                    "qstart": int(_row_float(row, "qstart")),
                    "qend": int(_row_float(row, "qend")),
                    "sstart": int(_row_float(row, "sstart")),
                    "send": int(_row_float(row, "send")),
                    "start": int(start) if float(start).is_integer() else float(start),
                    "end": int(end) if float(end).is_integer() else float(end),
                    "draw_start": float(draw_start),
                    "draw_end": float(draw_end),
                    "identity": _row_float(row, "identity"),
                    "alignment_length": _row_float(row, "alignment_length"),
                    "mismatches": _row_float(row, "mismatches"),
                    "gap_opens": _row_float(row, "gap_opens"),
                    "evalue": _row_float(row, "evalue"),
                    "bitscore": _row_float(row, "bitscore"),
                    "orientation": orientation,
                    "full_reference": bool(
                        draw_start == 0.0 and draw_end == float(record_len)
                    ),
                }
            )
        except Exception:
            dropped += 1

    if dropped:
        logger.warning(
            "WARNING: Dropped %s invalid conservation hit(s) for source %s and record %s.",
            dropped,
            source.source_index,
            record.id,
        )
    if not rows:
        return empty_normalized_conservation_hits()
    return DataFrame(rows, columns=NORMALIZED_CONSERVATION_COLUMNS)


def normalize_conservation_tracks_for_record(
    load_result: ConservationLoadResult,
    *,
    displayed_records: Sequence[SeqRecord],
    record: SeqRecord,
    conservation_reference: ConservationReferenceMode | str = "auto",
) -> tuple[ConservationTrack, ...]:
    """Create compact render tracks for one displayed circular record."""

    tracks: list[ConservationTrack] = []
    track_index = 0
    for source in load_result.sources:
        if source.skipped:
            continue
        track_index += 1
        reference_side = resolve_conservation_reference_side(
            source,
            displayed_records,
            conservation_reference,
        )
        hits = _normalize_source_hits_for_record(
            source,
            track_index=track_index,
            reference_side=reference_side,
            record=record,
        )
        tracks.append(
            ConservationTrack(
                source_index=int(source.source_index),
                track_index=track_index,
                track_label=source.label,
                track_color=source.color,
                reference_side=reference_side,
                hits=hits,
            )
        )
    return tuple(tracks)


__all__ = [
    "ConservationLoadResult",
    "ConservationSearchResult",
    "ConservationReferenceMode",
    "ConservationReferenceSide",
    "ConservationSource",
    "ConservationTrack",
    "NORMALIZED_CONSERVATION_COLUMNS",
    "conservation_track_gradient_colors",
    "empty_normalized_conservation_hits",
    "load_conservation_sources",
    "normalize_conservation_color",
    "normalize_conservation_reference",
    "normalize_conservation_tracks_for_record",
    "resolve_conservation_reference_side",
]
