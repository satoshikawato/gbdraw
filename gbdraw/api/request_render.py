"""Normalize and render typed diagram requests without CLI/session orchestration."""

from __future__ import annotations

import copy
import hashlib
import json
import re
from dataclasses import dataclass, fields, is_dataclass, replace
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    ProteinExtractionResult,
    extract_web_stable_cds_proteins,
    is_legacy_protein_losat_cache_entry,
    is_protein_losat_cache_entry,
    make_legacy_protein_raw_candidate,
    promote_legacy_protein_raw_cache_entries,
    proteins_to_fasta,
    validate_protein_raw_entry_references,
)
from gbdraw.features.visibility import (
    compile_feature_visibility_rules,
    read_feature_visibility_file,
    resolve_candidate_feature_types,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.io.record_select import reverse_records, select_record
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload

from .diagram import (
    DEFAULT_SELECTED_FEATURES,
    build_circular_diagram,
    build_circular_multi_diagram,
    build_linear_diagram,
)
from .io import apply_region_specs, load_gbks, load_gff_fasta
from .options import CircularMultiRecordOptions, DiagramOptions, LinearMultiRecordOptions
from .render import save_figure_to
from .requests import (
    CircularDiagramRequest,
    DiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)


@dataclass(frozen=True)
class PreparedDiagramRequest:
    """A validated request with normalized records and its SVG drawing."""

    mode: Literal["circular", "linear"]
    request: DiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    legacy_protein_raw_candidates: tuple[Mapping[str, Any], ...] = ()
    legacy_protein_derived_evidence: tuple[Mapping[str, Any], ...] = ()
    protein_id_map: Mapping[str, str] | None = None
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class RequestRenderResult:
    """Files and normalized inputs produced by one request render."""

    mode: Literal["circular", "linear"]
    request: DiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    output_paths: tuple[Path, ...]
    interactive_context: InteractiveSvgContext | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    legacy_protein_raw_candidates: tuple[Mapping[str, Any], ...] = ()
    legacy_protein_derived_evidence: tuple[Mapping[str, Any], ...] = ()
    protein_id_map: Mapping[str, str] | None = None
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class _PreparedLinearArtifacts:
    request: LinearDiagramRequest
    cache: LosatpCacheManager | None
    extraction: ProteinExtractionResult | None
    nucleotide_entries: tuple[Mapping[str, Any], ...]
    passthrough_derived_entries: tuple[Mapping[str, Any], ...]
    legacy_candidates: tuple[Mapping[str, Any], ...]
    protein_id_map: Mapping[str, str]
    source_mode: str
    warnings: tuple[str, ...]


_LEGACY_PROTEIN_REFERENCE_RE = re.compile(
    r"p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}"
    r"(?:_[2-9][0-9]*)?"
)


def _mode(request: DiagramRequest) -> Literal["circular", "linear"]:
    return "circular" if isinstance(request, CircularDiagramRequest) else "linear"


def _load_record_input(
    record_input: RecordInput,
    *,
    mode: Literal["circular", "linear"],
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    color_table: DataFrame | None,
    feature_visibility_table: DataFrame | None,
) -> SeqRecord:
    selector = record_input.selector
    selector_values = [selector.raw] if selector is not None else None
    reverse = record_input.presentation.reverse_complement
    source = record_input.source

    if isinstance(source, GenBankInputSource):
        records = load_gbks(
            [str(source.path)],
            mode=mode,
            record_selectors=selector_values,
            reverse_flags=[reverse],
        )
    elif isinstance(source, GffFastaInputSource):
        records = load_gff_fasta(
            [str(source.gff_path)],
            [str(source.fasta_path)],
            mode=mode,
            selected_features_set=gff_candidate_features,
            keep_all_features=gff_keep_all_features,
            record_selectors=selector_values,
            reverse_flags=[reverse],
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
    elif isinstance(source, InMemoryRecordSource):
        records = [copy.deepcopy(source.record)]
        try:
            records = select_record(records, selector)
            records = reverse_records(records, reverse)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
    else:  # pragma: no cover - RecordInput validates this union.
        raise ValidationError("Unsupported record input source.")

    if record_input.region is not None:
        records = apply_region_specs(records, [record_input.region])
    if len(records) != 1:
        raise ValidationError(
            "Each RecordInput must resolve to exactly one record; add a selector or region."
        )

    record = records[0]
    presentation = record_input.presentation
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    if presentation.label:
        record.annotations["gbdraw_record_label"] = presentation.label
    if presentation.subtitle:
        record.annotations["gbdraw_record_subtitle"] = presentation.subtitle
    if record_input.record_key:
        record.annotations["gbdraw_record_key"] = record_input.record_key
    return record


def normalize_request_records(request: DiagramRequest) -> tuple[SeqRecord, ...]:
    """Load/copy every RecordInput and return exactly one record per input."""

    if not isinstance(request, (CircularDiagramRequest, LinearDiagramRequest)):
        raise ValidationError("Unsupported diagram request type.")
    mode = _mode(request)
    has_gff_source = any(
        isinstance(record_input.source, GffFastaInputSource)
        for record_input in request.records
    )
    color_table = _color_table(request.options) if has_gff_source else None
    feature_visibility_table = (
        _visibility_table(request.options) if has_gff_source else None
    )
    candidate_features, keep_all_features = (
        resolve_candidate_feature_types(
            request.options.selected_features_set or DEFAULT_SELECTED_FEATURES,
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
        if has_gff_source
        else (set(request.options.selected_features_set or DEFAULT_SELECTED_FEATURES), False)
    )
    return tuple(
        _load_record_input(
            record_input,
            mode=mode,
            gff_candidate_features=tuple(sorted(candidate_features)),
            gff_keep_all_features=keep_all_features,
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
        for record_input in request.records
    )


def _layout_with_record_placements(
    request: CircularDiagramRequest,
) -> CircularMultiRecordOptions | None:
    layout = request.layout
    if layout is None or layout.multi_record_positions:
        return layout
    positioned = [
        (index, record_input.presentation)
        for index, record_input in enumerate(request.records)
        if record_input.presentation.grid_row is not None
    ]
    if not positioned:
        return layout
    positioned.sort(
        key=lambda item: (
            int(item[1].grid_row or 0),
            int(item[1].grid_column)
            if item[1].grid_column is not None
            else item[0],
            item[0],
        )
    )
    positions = tuple(
        f"#{index + 1}@{presentation.grid_row}"
        for index, presentation in positioned
    )
    return replace(layout, multi_record_positions=positions)


def _linear_layout_with_record_placements(
    request: LinearDiagramRequest,
) -> LinearMultiRecordOptions | None:
    layout = request.layout
    if layout is None or layout.multi_record_positions:
        return layout
    positioned = [
        (index, record_input.presentation)
        for index, record_input in enumerate(request.records)
        if record_input.presentation.grid_row is not None
    ]
    if not positioned:
        return layout
    positioned.sort(
        key=lambda item: (
            int(item[1].grid_row or 0),
            int(item[1].grid_column)
            if item[1].grid_column is not None
            else item[0],
            item[0],
        )
    )
    positions = tuple(
        f"#{index + 1}@{presentation.grid_row}"
        for index, presentation in positioned
    )
    return replace(layout, multi_record_positions=positions)


def _empty_protein_identity_manifest() -> dict[str, Any]:
    return {
        "schema": 1,
        "proteinSets": {},
        "recordAnalyses": {},
        "recordInstances": {},
    }


def _artifact_entries(
    artifacts: Mapping[str, Any] | None,
    field: str,
) -> tuple[Mapping[str, Any], ...]:
    if artifacts is None:
        return ()
    container = artifacts.get(field)
    if container is None:
        return ()
    if not isinstance(container, Mapping) or not isinstance(
        container.get("entries", []), list
    ):
        raise ValidationError(f"Session artifact {field}.entries must be an array.")
    entries = container.get("entries", [])
    assert isinstance(entries, list)
    if not all(isinstance(entry, Mapping) for entry in entries):
        raise ValidationError(f"Session artifact {field}.entries must contain objects.")
    return tuple(entries)


def _is_nucleotide_losat_entry(entry: Mapping[str, Any]) -> bool:
    return (
        entry.get("schema") == 2
        and entry.get("kind") == "raw-losat"
        and isinstance(entry.get("text"), str)
        and not is_legacy_protein_losat_cache_entry(entry)
    )


def _legacy_candidate_entries(
    artifacts: Mapping[str, Any] | None,
    direct_legacy_entries: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    candidates: list[Mapping[str, Any]] = []
    if artifacts is not None:
        legacy_artifacts = artifacts.get("legacyArtifacts")
        envelope = (
            legacy_artifacts.get("proteinRawCandidates")
            if isinstance(legacy_artifacts, Mapping)
            else None
        )
        if envelope is not None:
            if not isinstance(envelope, Mapping) or not isinstance(
                envelope.get("entries"), list
            ):
                raise ValidationError(
                    "legacyArtifacts.proteinRawCandidates must use an entries array."
                )
            candidates.extend(
                candidate
                for candidate in envelope["entries"]
                if isinstance(candidate, Mapping)
            )
    candidates.extend(
        make_legacy_protein_raw_candidate(entry) for entry in direct_legacy_entries
    )

    result: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    for candidate in candidates:
        original = candidate.get("originalEntry")
        if not isinstance(original, Mapping):
            raise ValidationError("Legacy protein candidate has no original entry.")
        fingerprint = json.dumps(
            original,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        result.append(candidate)
    return tuple(result)


def _merge_promoted_tsv_id_map(
    target: dict[str, str],
    *,
    legacy_text: str,
    current_text: str,
) -> None:
    def rows(text: str) -> list[list[str]]:
        return [
            line.split("\t")
            for line in str(text).splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]

    legacy_rows = rows(legacy_text)
    current_rows = rows(current_text)
    if len(legacy_rows) != len(current_rows):
        raise ValidationError("Promoted LOSATP TSV changed its row count.")
    for old_row, new_row in zip(legacy_rows, current_rows, strict=True):
        if len(old_row) < 2 or len(new_row) < 2:
            raise ValidationError("Promoted LOSATP TSV has an incomplete row.")
        for old_id, new_id in zip(old_row[:2], new_row[:2], strict=True):
            previous = target.get(old_id)
            if previous is not None and previous != new_id:
                raise ValidationError(
                    f"Legacy protein ID {old_id!r} resolves to multiple current IDs."
                )
            target[old_id] = new_id


def _promote_legacy_candidates(
    candidates: Sequence[Mapping[str, Any]],
    extraction: ProteinExtractionResult,
) -> tuple[
    tuple[Mapping[str, Any], ...],
    tuple[Mapping[str, Any], ...],
    Mapping[str, str],
]:
    manifest = extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Protein cache promotion requires an identity manifest.")
    fastas = tuple(proteins_to_fasta(items) for items in extraction.proteins_by_record)
    promoted_entries: list[Mapping[str, Any]] = []
    unresolved_candidates: list[Mapping[str, Any]] = []
    id_map: dict[str, str] = {}

    for candidate in candidates:
        original = candidate.get("originalEntry")
        if not isinstance(original, Mapping) or not is_legacy_protein_losat_cache_entry(
            original
        ):
            raise ValidationError("Legacy protein candidate is not a schema-2 blastp entry.")
        raw_args = original.get("args")
        expected_args = tuple(str(arg) for arg in raw_args) if isinstance(raw_args, list) else ()
        promotion = None
        rejection_reasons: list[str] = []
        for query_index, query_proteins in enumerate(extraction.proteins_by_record):
            for subject_index, subject_proteins in enumerate(
                extraction.proteins_by_record
            ):
                scan = promote_legacy_protein_raw_cache_entries(
                    (candidate,),
                    query_proteins=query_proteins,
                    subject_proteins=subject_proteins,
                    query_fasta=fastas[query_index],
                    subject_fasta=fastas[subject_index],
                    identity_manifest=manifest,
                    expected_args=expected_args,
                    expected_program=str(original.get("program") or "blastp"),
                    expected_outfmt=str(original.get("outfmt") or "6"),
                )
                if scan.promotion is not None:
                    promotion = scan.promotion
                    break
                rejection_reasons.extend(item.reason for item in scan.rejections)
            if promotion is not None:
                break

        if promotion is None:
            reason = next(
                (
                    item
                    for item in reversed(rejection_reasons)
                    if item and "hash does not match" not in item
                ),
                "No current record pair matched the legacy cache identity.",
            )
            unresolved_candidates.append(
                make_legacy_protein_raw_candidate(
                    original,
                    state="rejected",
                    rejection_reason=reason,
                )
            )
            continue
        if not validate_protein_raw_entry_references(promotion.entry, manifest):
            raise ValidationError("Promoted protein raw entry failed manifest validation.")
        promoted_entries.append(promotion.entry)
        _merge_promoted_tsv_id_map(
            id_map,
            legacy_text=str(original.get("text") or ""),
            current_text=promotion.rewritten_tsv,
        )

    return tuple(promoted_entries), tuple(unresolved_candidates), id_map


def _rewrite_protein_reference_string(value: str, id_map: Mapping[str, str]) -> str:
    if "p_r_" not in value:
        return value
    return _LEGACY_PROTEIN_REFERENCE_RE.sub(
        lambda match: id_map.get(match.group(0), match.group(0)),
        value,
    )


def rewrite_protein_artifact_references(
    value: Any,
    id_map: Mapping[str, str],
) -> Any:
    """Return a detached value with verified legacy protein IDs rewritten."""

    if isinstance(value, str):
        return _rewrite_protein_reference_string(value, id_map)
    if isinstance(value, DataFrame):
        rewritten = value.copy(deep=True)
        for column in rewritten.columns:
            rewritten[column] = rewritten[column].map(
                lambda item: (
                    _rewrite_protein_reference_string(item, id_map)
                    if isinstance(item, str)
                    else item
                )
            )
        return rewritten
    if is_dataclass(value) and not isinstance(value, type):
        return replace(
            value,
            **{
                field.name: rewrite_protein_artifact_references(
                    getattr(value, field.name), id_map
                )
                for field in fields(value)
            },
        )
    if isinstance(value, Mapping):
        return {
            rewrite_protein_artifact_references(key, id_map):
            rewrite_protein_artifact_references(item, id_map)
            for key, item in value.items()
        }
    if isinstance(value, tuple):
        return tuple(rewrite_protein_artifact_references(item, id_map) for item in value)
    if isinstance(value, list):
        return [rewrite_protein_artifact_references(item, id_map) for item in value]
    if isinstance(value, frozenset):
        return frozenset(
            rewrite_protein_artifact_references(item, id_map) for item in value
        )
    if isinstance(value, set):
        return {rewrite_protein_artifact_references(item, id_map) for item in value}
    return copy.deepcopy(value)


def _contains_legacy_protein_reference(value: Any) -> bool:
    if isinstance(value, str):
        return _LEGACY_PROTEIN_REFERENCE_RE.search(value) is not None
    if isinstance(value, DataFrame):
        return any(
            _contains_legacy_protein_reference(item)
            for item in value.to_numpy().ravel()
        )
    if is_dataclass(value) and not isinstance(value, type):
        return any(
            _contains_legacy_protein_reference(getattr(value, field.name))
            for field in fields(value)
        )
    if isinstance(value, Mapping):
        return any(
            _contains_legacy_protein_reference(key)
            or _contains_legacy_protein_reference(item)
            for key, item in value.items()
        )
    if isinstance(value, (tuple, list, set, frozenset)):
        return any(_contains_legacy_protein_reference(item) for item in value)
    return False


def _rewrite_linear_request_protein_references(
    request: LinearDiagramRequest,
    id_map: Mapping[str, str],
) -> LinearDiagramRequest:
    if not id_map:
        return request
    options = request.options
    rewritten_options = replace(
        options,
        linear_comparisons=rewrite_protein_artifact_references(
            options.linear_comparisons, id_map
        ),
        protein_comparisons=rewrite_protein_artifact_references(
            options.protein_comparisons, id_map
        ),
        orthogroups=rewrite_protein_artifact_references(options.orthogroups, id_map),
        collinearity_blocks=rewrite_protein_artifact_references(
            options.collinearity_blocks, id_map
        ),
    )
    for value in (
        rewritten_options.linear_comparisons,
        rewritten_options.protein_comparisons,
        rewritten_options.orthogroups,
        rewritten_options.collinearity_blocks,
    ):
        if _contains_legacy_protein_reference(value):
            raise ValidationError(
                "Precomputed protein artifacts contain an ID that no verified raw cache entry resolved."
            )
    return replace(request, options=rewritten_options)


def _source_protein_mode(
    request: LinearDiagramRequest,
    artifacts: Mapping[str, Any] | None,
) -> str:
    requested = str(request.options.protein_blastp_mode or "none")
    if requested != "none":
        return requested
    if artifacts is not None:
        config = artifacts.get("config")
        losat = config.get("losat") if isinstance(config, Mapping) else None
        blastp = losat.get("blastp") if isinstance(losat, Mapping) else None
        configured = blastp.get("mode") if isinstance(blastp, Mapping) else None
        if str(configured or "") in {"pairwise", "orthogroup", "collinear"}:
            return str(configured)
        for entry in _artifact_entries(artifacts, "losatDerivedCache"):
            mode = str(entry.get("mode") or "")
            if mode in {"pairwise", "orthogroup", "collinear"}:
                return mode
    return "orthogroup" if request.options.orthogroups is not None else "pairwise"


def _prepare_linear_artifacts(
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    artifacts: Mapping[str, Any] | None,
) -> _PreparedLinearArtifacts:
    raw_entries = _artifact_entries(artifacts, "losatCache")
    current_protein = tuple(
        entry for entry in raw_entries if is_protein_losat_cache_entry(entry)
    )
    direct_legacy = tuple(
        entry for entry in raw_entries if is_legacy_protein_losat_cache_entry(entry)
    )
    nucleotide_entries = tuple(
        copy.deepcopy(dict(entry))
        for entry in raw_entries
        if _is_nucleotide_losat_entry(entry)
    )
    candidates = _legacy_candidate_entries(artifacts, direct_legacy)
    needs_protein_identity = bool(
        current_protein
        or candidates
        or request.options.protein_blastp_mode != "none"
        or request.options.protein_comparisons is not None
        or request.options.orthogroups is not None
        or request.options.collinearity_blocks is not None
    )
    passthrough_derived = tuple(
        copy.deepcopy(dict(entry))
        for entry in _artifact_entries(artifacts, "losatDerivedCache")
        if entry.get("schema") == 2
        and entry.get("kind") == "derived-losatp-payload"
    )
    if not needs_protein_identity:
        return _PreparedLinearArtifacts(
            request=request,
            cache=None,
            extraction=None,
            nucleotide_entries=nucleotide_entries,
            passthrough_derived_entries=passthrough_derived,
            legacy_candidates=(),
            protein_id_map={},
            source_mode="none",
            warnings=(),
        )

    record_instance_keys = tuple(
        record_input.record_key or f"record-{index + 1}"
        for index, record_input in enumerate(request.records)
    )
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=record_instance_keys,
        record_source_ids=tuple(str(record.id) for record in records),
        record_selectors=tuple(
            record_input.selector.raw if record_input.selector is not None else None
            for record_input in request.records
        ),
        regions=tuple(
            record_input.region.raw if record_input.region is not None else None
            for record_input in request.records
        ),
        feature_visibility_rules=compile_feature_visibility_rules(
            _visibility_table(request.options)
        ),
    )
    manifest = extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Protein extraction did not produce an identity manifest.")
    promoted, unresolved, id_map = _promote_legacy_candidates(candidates, extraction)
    reusable_current = tuple(
        entry
        for entry in current_protein
        if validate_protein_raw_entry_references(entry, manifest)
    )
    cache = LosatpCacheManager(
        (*reusable_current, *promoted),
        identity_manifest=manifest,
        threads_per_job=request.options.losatp_threads or "auto",
    )
    rewritten_request = _rewrite_linear_request_protein_references(request, id_map)
    # A canonical precomputed replay has completed its current conversion once
    # every referenced protein ID was resolved above.  Unmatched legacy entries
    # are stale cache history, not candidates for the rendered request.
    if request.options.protein_comparisons is not None and promoted:
        unresolved = ()
    warnings = (
        (
            f"{len(unresolved)} legacy protein raw cache candidate(s) could not be promoted."
        ),
    ) if unresolved else ()
    return _PreparedLinearArtifacts(
        request=rewritten_request,
        cache=cache,
        extraction=extraction,
        nucleotide_entries=nucleotide_entries,
        passthrough_derived_entries=(),
        legacy_candidates=unresolved,
        protein_id_map=id_map,
        source_mode=_source_protein_mode(request, artifacts),
        warnings=warnings,
    )


def _comparison_frame_payload(
    frame: DataFrame,
    *,
    pair_index: int,
    query_index: int,
    subject_index: int,
) -> dict[str, Any]:
    missing = [column for column in COMPARISON_COLUMNS if column not in frame.columns]
    if missing:
        raise ValidationError(
            "Protein comparison is missing canonical column(s): "
            + ", ".join(missing)
        )
    tsv = frame.loc[:, list(COMPARISON_COLUMNS)].to_csv(
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n",
    )
    rows = json.loads(frame.to_json(orient="records"))
    return {
        "pair_index": int(pair_index),
        "query_index": int(query_index),
        "subject_index": int(subject_index),
        "tsv": tsv,
        "rows": rows,
        "hit_count": len(frame),
    }


def _resolved_protein_pair_payloads(
    drawing: Drawing,
    request: LinearDiagramRequest,
) -> list[dict[str, Any]]:
    direct = getattr(drawing, "_gbdraw_resolved_protein_comparisons", None)
    if direct is None:
        direct = request.options.protein_comparisons
    explicit = getattr(drawing, "_gbdraw_resolved_linear_comparisons", None)

    result: list[dict[str, Any]] = []
    for pair_index, frame in enumerate(direct or ()):
        if not isinstance(frame, DataFrame):
            continue
        result.append(
            _comparison_frame_payload(
                frame,
                pair_index=pair_index,
                query_index=pair_index,
                subject_index=pair_index + 1,
            )
        )
    if result:
        return result

    explicit = explicit if explicit is not None else request.options.linear_comparisons
    for comparison in explicit or ():
        frame = getattr(comparison, "matches", None)
        if not isinstance(frame, DataFrame) or not {
            "query_protein_id",
            "subject_protein_id",
        }.issubset(frame.columns):
            continue
        result.append(
            _comparison_frame_payload(
                frame,
                pair_index=len(result),
                query_index=int(comparison.query_record_index),
                subject_index=int(comparison.subject_record_index),
            )
        )
    return result


def _build_current_derived_entries(
    drawing: Drawing,
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    artifacts: _PreparedLinearArtifacts,
    raw_entries: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    if artifacts.extraction is None:
        return artifacts.passthrough_derived_entries
    mode = artifacts.source_mode
    if mode not in {"orthogroup", "collinear"}:
        return ()
    manifest = artifacts.extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Derived protein payload requires an identity manifest.")
    pair_payloads = _resolved_protein_pair_payloads(drawing, request)
    orthogroups = getattr(drawing, "_gbdraw_orthogroups", request.options.orthogroups)
    orthogroup_payload = serialize_orthogroups_payload(
        orthogroups,
        records=records,
    )
    if not pair_payloads and not orthogroup_payload:
        return ()

    record_instances = manifest.record_instances
    instance_keys = tuple(
        record_input.record_key or f"record-{index + 1}"
        for index, record_input in enumerate(request.records)
    )
    protein_raw_entries = [
        entry for entry in raw_entries if is_protein_losat_cache_entry(entry)
    ]
    identity = {
        "cacheSchema": 2,
        "converter": "convert_losatp_blastp_pairs_to_genomic_payload",
        "mode": mode,
        "maxHits": int(request.options.protein_blastp_max_hits),
        "thresholds": {
            "bitscore": str(request.options.bitscore),
            "evalue": str(request.options.evalue),
            "identity": str(request.options.identity),
            "alignmentLength": str(request.options.alignment_length),
        },
        "orthogroup": {
            "membershipMode": str(request.options.orthogroup_membership_mode),
            "memberMaxHits": int(request.options.orthogroup_member_max_hits),
        },
        "collinear": {
            "unitMode": str(request.options.collinearity_unit_mode),
            "colorMode": str(request.options.collinearity_color_mode),
            "anchorMode": str(request.options.collinearity_anchor_mode),
            "searchScope": str(request.options.collinearity_search_scope),
            "maxParalogLinksPerOrthogroup": int(
                request.options.collinear_max_paralog_links_per_orthogroup
            ),
        },
        "records": [
            {
                "recordIndex": index,
                "recordInstanceKey": instance_key,
                "bindingHash": str(
                    record_instances.get(instance_key, {}).get("bindingHash") or ""
                ),
                "viewTransform": {
                    "length": len(records[index].seq),
                    "reverse": bool(
                        request.records[index].presentation.reverse_complement
                        or (
                            request.records[index].region is not None
                            and request.records[index].region.reverse_complement
                        )
                    ),
                },
            }
            for index, instance_key in enumerate(instance_keys)
        ],
        "pairs": [
            {
                "pairIndex": int(pair["pair_index"]),
                "queryIndex": int(pair["query_index"]),
                "subjectIndex": int(pair["subject_index"]),
                "rawCacheKeys": sorted(
                    str(entry.get("key") or "")
                    for entry in protein_raw_entries
                    if entry.get("queryRecordInstanceKey")
                    == instance_keys[int(pair["query_index"])]
                    and entry.get("subjectRecordInstanceKey")
                    == instance_keys[int(pair["subject_index"])]
                ),
            }
            for pair in pair_payloads
        ],
    }
    key = hashlib.sha256(
        json.dumps(
            identity,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    entry: Mapping[str, Any] = {
        "schema": 2,
        "kind": "derived-losatp-payload",
        "key": key,
        "mode": mode,
        "payload": {
            "identity": identity,
            "pairs": pair_payloads,
            "orthogroups": orthogroup_payload,
        },
    }
    if _contains_legacy_protein_reference(entry):
        raise ValidationError("Derived protein payload retained a legacy protein ID.")
    return (entry,)


def _legacy_derived_evidence_entries(
    artifacts: Mapping[str, Any] | None,
) -> tuple[Mapping[str, Any], ...]:
    if artifacts is None:
        return ()
    legacy = artifacts.get("legacyArtifacts")
    envelope = (
        legacy.get("proteinDerivedEvidence")
        if isinstance(legacy, Mapping)
        else None
    )
    if envelope is None:
        return ()
    if not isinstance(envelope, Mapping) or not isinstance(
        envelope.get("entries"), list
    ):
        raise ValidationError(
            "legacyArtifacts.proteinDerivedEvidence must use an entries array."
        )
    return tuple(
        copy.deepcopy(dict(entry))
        for entry in envelope["entries"]
        if isinstance(entry, Mapping)
    )


def build_request_diagram(
    request: DiagramRequest,
    *,
    session_artifacts: Mapping[str, Any] | None = None,
) -> PreparedDiagramRequest:
    """Normalize inputs and build a drawing through the high-level API owners."""

    records = normalize_request_records(request)
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    legacy_candidates: tuple[Mapping[str, Any], ...] = ()
    legacy_derived_evidence: tuple[Mapping[str, Any], ...] = ()
    protein_id_map: Mapping[str, str] = {}
    warnings: tuple[str, ...] = ()
    if isinstance(request, CircularDiagramRequest):
        layout = _layout_with_record_placements(request)
        if layout is None:
            drawing = build_circular_diagram(records[0], options=request.options)
        else:
            drawing = build_circular_multi_diagram(
                records,
                options=request.options,
                layout=layout,
            )
        mode: Literal["circular", "linear"] = "circular"
        if session_artifacts is not None:
            raw_entries = _artifact_entries(session_artifacts, "losatCache")
            losat_cache_entries = tuple(
                copy.deepcopy(dict(entry))
                for entry in raw_entries
                if is_protein_losat_cache_entry(entry)
                or _is_nucleotide_losat_entry(entry)
            )
            direct_legacy = tuple(
                entry
                for entry in raw_entries
                if is_legacy_protein_losat_cache_entry(entry)
            )
            legacy_candidates = _legacy_candidate_entries(
                session_artifacts,
                direct_legacy,
            )
            losat_derived_cache_entries = tuple(
                copy.deepcopy(dict(entry))
                for entry in _artifact_entries(
                    session_artifacts,
                    "losatDerivedCache",
                )
                if entry.get("schema") == 2
                and entry.get("kind") == "derived-losatp-payload"
            )
            source_manifest = session_artifacts.get("proteinIdentityManifest")
            protein_identity_manifest = (
                copy.deepcopy(dict(source_manifest))
                if isinstance(source_manifest, Mapping)
                and source_manifest.get("schema") == 1
                else _empty_protein_identity_manifest()
            )
            legacy_derived_evidence = _legacy_derived_evidence_entries(
                session_artifacts
            )
    elif isinstance(request, LinearDiagramRequest):
        artifacts = _prepare_linear_artifacts(request, records, session_artifacts)
        request = artifacts.request
        linear_layout = _linear_layout_with_record_placements(request)
        build_kwargs: dict[str, Any] = {"options": request.options}
        if artifacts.cache is not None:
            build_kwargs["losatp_cache"] = artifacts.cache
        if artifacts.extraction is not None:
            build_kwargs["protein_extraction"] = artifacts.extraction
        if linear_layout is None:
            drawing = build_linear_diagram(records, **build_kwargs)
        else:
            drawing = build_linear_diagram(
                records,
                layout=linear_layout,
                **build_kwargs,
            )
        mode = "linear"
        protein_entries = (
            tuple(artifacts.cache.session_entries())
            if artifacts.cache is not None
            else ()
        )
        losat_cache_entries = (
            *protein_entries,
            *artifacts.nucleotide_entries,
        )
        losat_derived_cache_entries = _build_current_derived_entries(
            drawing,
            request,
            records,
            artifacts,
            losat_cache_entries,
        )
        protein_identity_manifest = (
            artifacts.extraction.identity_manifest.to_dict()
            if artifacts.extraction is not None
            and artifacts.extraction.identity_manifest is not None
            else (
                _empty_protein_identity_manifest()
                if session_artifacts is not None
                else None
            )
        )
        legacy_candidates = artifacts.legacy_candidates
        protein_id_map = artifacts.protein_id_map
        warnings = artifacts.warnings
    else:  # pragma: no cover - normalize_request_records rejects this first.
        raise ValidationError("Unsupported diagram request type.")
    return PreparedDiagramRequest(
        mode=mode,
        request=request,
        records=records,
        drawing=drawing,
        losat_cache_entries=losat_cache_entries,
        losat_derived_cache_entries=losat_derived_cache_entries,
        protein_identity_manifest=protein_identity_manifest,
        legacy_protein_raw_candidates=legacy_candidates,
        legacy_protein_derived_evidence=legacy_derived_evidence,
        protein_id_map=protein_id_map,
        warnings=warnings,
    )


def _visibility_table(options: DiagramOptions) -> DataFrame | None:
    table = (
        options.feature_visibility_table
        if options.feature_visibility_table is not None
        else options.feature_table
    )
    file_path = (
        options.feature_visibility_table_file
        if options.feature_visibility_table_file is not None
        else options.feature_table_file
    )
    if table is None and file_path is not None:
        return read_feature_visibility_file(file_path)
    return table


def _color_table(options: DiagramOptions) -> DataFrame | None:
    colors = options.colors
    if colors is None or colors.color_table is not None:
        return colors.color_table if colors is not None else None
    if colors.color_table_file is not None:
        return read_color_table(colors.color_table_file)
    return None


def _interactive_context(
    prepared: PreparedDiagramRequest,
) -> InteractiveSvgContext | None:
    output = prepared.request.output
    if "interactive_svg" not in output.formats or output.interactive_metadata_policy == "omit":
        return None

    options = prepared.request.options
    colors = options.colors
    color_table = _color_table(options)
    default_colors = colors.default_colors if colors is not None else None
    if default_colors is None:
        default_colors = load_default_colors(
            colors.default_colors_file if colors is not None and colors.default_colors_file else "",
            colors.default_colors_palette if colors is not None else "default",
        )
    return build_interactive_svg_context(
        prepared.records,
        selected_features_set=options.selected_features_set,
        feature_table=_visibility_table(options),
        color_table=color_table,
        default_colors=default_colors,
        orthogroups=options.orthogroups,
        linear_rendered_feature_ids=prepared.mode == "linear",
        annotations=options.annotations,
        mode=prepared.mode,
    )


def render_request(
    request: DiagramRequest,
    *,
    session_artifacts: Mapping[str, Any] | None = None,
) -> RequestRenderResult:
    """Build and save one typed request, returning only paths that were created."""

    prepared = (
        build_request_diagram(request)
        if session_artifacts is None
        else build_request_diagram(request, session_artifacts=session_artifacts)
    )
    output = request.output
    interactive_context = _interactive_context(prepared)
    paths = save_figure_to(
        prepared.drawing,
        output.formats,
        output_dir=(
            str(output.output_directory)
            if output.output_directory is not None
            else None
        ),
        output_prefix=output.output_prefix,
        overwrite=output.overwrite,
        interactive_context=interactive_context,
    )
    return RequestRenderResult(
        mode=prepared.mode,
        request=prepared.request,
        records=prepared.records,
        drawing=prepared.drawing,
        output_paths=tuple(Path(path) for path in paths),
        interactive_context=interactive_context,
        losat_cache_entries=prepared.losat_cache_entries,
        losat_derived_cache_entries=prepared.losat_derived_cache_entries,
        protein_identity_manifest=prepared.protein_identity_manifest,
        legacy_protein_raw_candidates=prepared.legacy_protein_raw_candidates,
        legacy_protein_derived_evidence=prepared.legacy_protein_derived_evidence,
        protein_id_map=prepared.protein_id_map,
        warnings=prepared.warnings,
    )


__all__ = [
    "PreparedDiagramRequest",
    "RequestRenderResult",
    "build_request_diagram",
    "normalize_request_records",
    "render_request",
    "rewrite_protein_artifact_references",
]
