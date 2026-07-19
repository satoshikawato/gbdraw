"""Circular conservation-ring input loading and normalization."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Iterable, Literal, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.core.color import normalize_hex_color, tint_color  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.io.colors import resolve_color_to_hex  # type: ignore[reportMissingImports]
from gbdraw.io.comparisons import (  # type: ignore[reportMissingImports]
    COMPARISON_COLUMNS,
    filter_comparison_dataframe,
)


logger = logging.getLogger(__name__)

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

_NUMERIC_COMPARISON_COLUMNS = (
    "identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
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
    return normalized  # type: ignore[return-value]


def empty_normalized_conservation_hits() -> DataFrame:
    return DataFrame(columns=NORMALIZED_CONSERVATION_COLUMNS)


def _default_label(source_index: int, path: str | None) -> str:
    if path:
        basename = os.path.basename(str(path))
        if basename:
            return basename
    return f"Conservation {int(source_index) + 1}"


def normalize_conservation_color(value: object | None) -> str | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return normalize_hex_color(resolve_color_to_hex(text))
    except Exception as exc:
        raise ValidationError(f"Invalid conservation color {text!r}. Use an SVG color name or #RRGGBB.") from exc


def conservation_track_gradient_colors(
    track_color: object | None,
    *,
    default_min_color: str,
    default_max_color: str,
) -> tuple[str, str]:
    normalized_track_color = normalize_conservation_color(track_color)
    if normalized_track_color is None:
        return normalize_hex_color(default_min_color), normalize_hex_color(default_max_color)
    return tint_color(normalized_track_color), normalized_track_color


def _coerce_comparison_dataframe(dataframe: DataFrame) -> DataFrame:
    """Return a BLAST outfmt 6/7 shaped dataframe or raise ValueError."""

    if dataframe is None:
        raise ValueError("dataframe is None")

    if set(COMPARISON_COLUMNS).issubset(set(dataframe.columns)):
        df = dataframe.loc[:, list(COMPARISON_COLUMNS)].copy()
    elif len(dataframe.columns) >= len(COMPARISON_COLUMNS):
        df = dataframe.iloc[:, : len(COMPARISON_COLUMNS)].copy()
        df.columns = list(COMPARISON_COLUMNS)
    else:
        raise ValueError(
            "comparison dataframe must contain BLAST outfmt 6 columns "
            f"({', '.join(COMPARISON_COLUMNS)})"
        )

    for column in ("query", "subject"):
        df[column] = df[column].astype(str)
    for column in _NUMERIC_COMPARISON_COLUMNS:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def _filter_valid_dataframe(dataframe: DataFrame, blast_config: object) -> DataFrame:
    df = _coerce_comparison_dataframe(dataframe)
    # Keep the original source-row identity through filtering and paint-order sorting.
    df["source_hit_index"] = range(len(df))
    filtered = filter_comparison_dataframe(df, blast_config)  # type: ignore[arg-type]
    return filtered.loc[:, [*COMPARISON_COLUMNS, "source_hit_index"]].reset_index(drop=True)


def _load_conservation_file(path: str, blast_config: object) -> tuple[DataFrame | None, str | None]:
    if not os.path.isfile(path):
        return None, f"file does not exist or is not accessible: {path}"
    try:
        raw_df = pd.read_csv(
            path,
            sep="\t",
            comment="#",
            names=COMPARISON_COLUMNS,
        )
        return _filter_valid_dataframe(raw_df, blast_config), None
    except pd.errors.EmptyDataError:
        raw_df = pd.DataFrame(columns=COMPARISON_COLUMNS)
        return _filter_valid_dataframe(raw_df, blast_config), None
    except Exception as exc:
        return None, f"error parsing BLAST file for similarity ring {path}: {exc}"


def load_conservation_sources(
    *,
    blast_config: object,
    conservation_files: Sequence[str] | None = None,
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
            frame, reason = _load_conservation_file(str(path), blast_config)
            if frame is not None:
                valid_frames.append(frame)
            elif reason:
                logger.warning("WARNING: Skipping conservation source %s file: %s", source_index, reason)
                skip_reasons.append(reason)

        if source_index < len(dataframes):
            try:
                valid_frames.append(_filter_valid_dataframe(dataframes[source_index], blast_config))
            except Exception as exc:
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
        return mode  # type: ignore[return-value]

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
