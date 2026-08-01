#!/usr/bin/env python
# coding: utf-8

"""Shared CLI helpers for GUI session JSON input and sidecar output."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence

from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    read_circular_track_table,
    read_conservation_table,
    read_comparisons_table,
    read_records_table,
)
from gbdraw.render.formats import (
    SVG_FORMAT,
    resolve_format_output_path,
    resolve_output_paths,
)
from gbdraw.render.output_paths import preflight_output_paths
from gbdraw.render.track_slot_metadata import (
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
)
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS,
    SessionBuildContext,
    SessionFileBinding,
    build_session_json,
    migrate_legacy_linear_comparison_draft_for_current_writer,
    migrate_persisted_web_state_field_names,
    safe_embedded_filename,
    serialize_file_entry,
    validate_current_web_state_field_names,
    write_session_json,
)

if TYPE_CHECKING:
    from gbdraw.api.requests import DiagramRequest
    from gbdraw.render.interactive_svg import InteractiveSvgContext


@dataclass(frozen=True)
class RenderedSvg:
    output_prefix: str
    svg_path: Path
    result_name: str


@dataclass(frozen=True)
class DiagramRunResult:
    mode: Literal["circular", "linear"]
    render_formats: tuple[str, ...]
    outputs: tuple[RenderedSvg, ...]
    feature_metadata: tuple[Mapping[str, Any], ...] = ()
    orthogroup_metadata: tuple[Mapping[str, Any], ...] | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] | None = None
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] | None = None
    protein_identity_manifest: Mapping[str, Any] | None = None
    legacy_protein_raw_candidates: tuple[Mapping[str, Any], ...] | None = None
    legacy_protein_derived_evidence: tuple[Mapping[str, Any], ...] | None = None
    linear_record_metadata: tuple[Mapping[str, Any], ...] = ()
    run_metadata: Mapping[str, Any] = field(default_factory=dict)
    canonical_request: DiagramRequest | None = None
    biological_feature_metadata: tuple[Mapping[str, Any], ...] = ()
    interactive_contexts: tuple[InteractiveSvgContext | None, ...] = ()


@dataclass(frozen=True)
class SessionCliRequest:
    session_path: str
    output: str | None
    format: str | None
    overwrite: bool
    save_session: bool
    session_output: str | None


def add_session_args(parser: argparse.ArgumentParser) -> None:
    """Add session input/output options to a diagram parser."""

    parser.add_argument(
        "--session",
        help=(
            "Regenerate a diagram from a plain or gzip-compressed gbdraw GUI "
            "session JSON file."
        ),
        type=str,
    )
    parser.add_argument(
        "--save_session",
        help="Write one GUI-loadable .gbdraw-session.json sidecar for this run.",
        action="store_true",
    )
    parser.add_argument(
        "--session_output",
        metavar="PATH",
        help=(
            "Write the session sidecar to PATH; use a .gz suffix for gzip "
            "compression; implies --save_session."
        ),
        type=str,
    )


def parse_session_pre_args(
    cmd_args: Sequence[str],
    *,
    mode: Literal["circular", "linear"],
) -> SessionCliRequest | None:
    """Pre-parse --session invocations and reject unsupported override options."""

    if "-h" in cmd_args or "--help" in cmd_args:
        return None
    if "--session" not in cmd_args:
        return None

    parser = argparse.ArgumentParser(
        prog=f"gbdraw {mode}",
        add_help=False,
    )
    parser.add_argument("--session", required=True)
    parser.add_argument("-o", "--output")
    parser.add_argument("-f", "--format")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--save_session", action="store_true")
    parser.add_argument("--session_output")
    namespace, unknown = parser.parse_known_args(list(cmd_args))
    if unknown:
        parser.error(
            "--session cannot be combined with unsupported option(s): "
            + " ".join(unknown)
        )
    return SessionCliRequest(
        session_path=str(namespace.session),
        output=namespace.output,
        format=namespace.format,
        overwrite=bool(namespace.overwrite),
        save_session=bool(namespace.save_session or namespace.session_output),
        session_output=namespace.session_output,
    )


def resolve_session_sidecar_path(
    *,
    explicit_path: str | None,
    output_prefix: str | None,
    outputs: Sequence[RenderedSvg],
) -> Path:
    """Resolve the run-level session sidecar path."""

    if explicit_path:
        return Path(explicit_path)
    if output_prefix:
        return Path(f"{output_prefix}.gbdraw-session.json")
    if len(outputs) == 1:
        return outputs[0].svg_path.with_suffix(".gbdraw-session.json")
    return Path("gbdraw.gbdraw-session.json")


def preflight_session_sidecar_if_requested(
    *,
    save_session: bool,
    session_output: str | None,
    output_prefix: str | None,
    outputs: Sequence[RenderedSvg] = (),
    diagram_output_paths: Sequence[str | Path] = (),
    overwrite: bool = False,
) -> Path | None:
    """Reject sidecar collisions before rendering when its path is known."""

    if not save_session and not session_output:
        return None
    if not session_output and not output_prefix and not outputs:
        return None
    sidecar_path = resolve_session_sidecar_path(
        explicit_path=session_output,
        output_prefix=output_prefix,
        outputs=outputs,
    )
    preflight_output_paths((sidecar_path,), overwrite=True)
    try:
        sidecar_identity = sidecar_path.resolve(strict=False)
        diagram_identities = tuple(
            (Path(path), Path(path).resolve(strict=False))
            for path in diagram_output_paths
        )
    except (OSError, ValueError) as exc:
        raise ValidationError(
            f"Could not resolve output path: {sidecar_path}."
        ) from exc
    colliding_output = next(
        (
            Path(path)
            for path, identity in diagram_identities
            if identity == sidecar_identity
        ),
        None,
    )
    if colliding_output is not None:
        raise ValidationError(
            f"Session output path collides with diagram output: {sidecar_path}. "
            "Choose a distinct --session_output path."
        )
    if sidecar_path.exists() and not overwrite:
        raise ValidationError(
            f"Session output already exists: {sidecar_path}. "
            "Use --overwrite to replace it."
        )
    return sidecar_path


def diagram_request_output_paths(request: DiagramRequest) -> tuple[Path, ...]:
    """Return every file path one resolved typed request will write."""

    from gbdraw.api.requests import CircularBatchRequest

    outputs = (
        request.outputs
        if isinstance(request, CircularBatchRequest)
        else (request.output,)
    )
    return tuple(
        Path(path)
        for output in outputs
        for path in resolve_output_paths(
            str(
                Path(output.output_directory or ".")
                / output.output_prefix
            ),
            output.formats,
            include_base_svg=True,
        )
    )


def diagram_request_rendered_svgs(
    request: DiagramRequest,
) -> tuple[RenderedSvg, ...]:
    """Project resolved typed outputs to the run-level SVG result contract."""

    from gbdraw.api.requests import CircularBatchRequest

    outputs = (
        request.outputs
        if isinstance(request, CircularBatchRequest)
        else (request.output,)
    )
    return tuple(
        make_rendered_svg(
            str(Path(output.output_directory or ".") / output.output_prefix),
            output.output_prefix,
        )
        for output in outputs
    )


def make_rendered_svg(output_prefix: str, result_name: str | None = None) -> RenderedSvg:
    """Create a RenderedSvg result for a static SVG export."""

    svg_path = Path(resolve_format_output_path(output_prefix, SVG_FORMAT))
    return RenderedSvg(
        output_prefix=str(output_prefix),
        svg_path=svg_path,
        result_name=result_name or svg_path.stem,
    )


def _feature_catalog_for_svg_results(
    svg_results: Sequence[tuple[str, str]],
    contexts: Sequence[InteractiveSvgContext | None],
) -> dict[str, object]:
    from gbdraw.render.interactive_svg import InteractiveSvgContext
    from gbdraw.web_support.feature_catalog import (
        build_feature_catalog,
        build_feature_catalog_item,
    )

    if not svg_results:
        return build_feature_catalog([])
    if contexts and len(contexts) != len(svg_results):
        raise ValidationError(
            "Session feature metadata must contain one context per Result."
        )
    aligned_contexts = (
        tuple(contexts)
        if contexts
        else tuple(None for _ in svg_results)
    )
    items = []
    for result_index, ((result_name, svg_source), context) in enumerate(
        zip(svg_results, aligned_contexts, strict=True)
    ):
        if context is None:
            items.append(
                {
                    "resultIndex": result_index,
                    "resultName": result_name,
                    "recordKeys": [],
                    "features": [],
                    "biologicalFeatures": [],
                    "orthogroups": [],
                    "annotations": [],
                    "comparisonMatches": [],
                }
            )
            continue
        if not isinstance(context, InteractiveSvgContext):
            raise ValidationError(
                "Session feature metadata contains an invalid render context."
            )
        items.append(
            build_feature_catalog_item(
                svg_source,
                context,
                result_index=result_index,
                result_name=result_name,
            )
        )
    return build_feature_catalog(items)


def _replace_current_derived_feature_state(
    payload: dict[str, Any],
    feature_catalog: Mapping[str, object],
) -> None:
    editor_state = payload.get("editorState")
    editor_state = (
        dict(editor_state) if isinstance(editor_state, Mapping) else {}
    )
    editor_state["featureCatalog"] = dict(feature_catalog)
    payload["editorState"] = editor_state

    features = payload.get("features")
    features = dict(features) if isinstance(features, Mapping) else {}
    for key in CURRENT_WRITER_FORBIDDEN_FEATURE_FIELDS:
        features.pop(key, None)
    payload["features"] = features

    orthogroup_state = payload.get("orthogroupState")
    orthogroup_state = (
        dict(orthogroup_state)
        if isinstance(orthogroup_state, Mapping)
        else {}
    )
    orthogroup_state.pop("groups", None)
    payload["orthogroupState"] = orthogroup_state


def save_session_sidecar_if_requested(
    *,
    save_session: bool,
    session_output: str | None,
    output_prefix: str | None,
    run_result: DiagramRunResult,
    cmd_args: Sequence[str] | None = None,
    source_session: Mapping[str, Any] | None = None,
    cli_invocation_args: Sequence[str] = (),
    file_bindings: Sequence[SessionFileBinding] = (),
    overwrite: bool = False,
) -> Path | None:
    """Build and write a GUI session sidecar when requested."""

    if not save_session and not session_output:
        return None
    sidecar_path = preflight_session_sidecar_if_requested(
        save_session=save_session,
        session_output=session_output,
        output_prefix=output_prefix,
        outputs=run_result.outputs,
        diagram_output_paths=(
            diagram_request_output_paths(run_result.canonical_request)
            if run_result.canonical_request is not None
            else tuple(
                Path(path)
                for output in run_result.outputs
                for path in resolve_output_paths(
                    output.output_prefix,
                    run_result.render_formats,
                    include_base_svg=True,
                )
            )
        ),
        overwrite=overwrite,
    )
    assert sidecar_path is not None

    if source_session is not None:
        embedded_files = source_session.get("files")
        if not isinstance(embedded_files, Mapping):
            raise ValidationError("Source session has no files to preserve.")
        session_files: Mapping[str, Any] = embedded_files
        invocation_args = tuple(str(arg) for arg in cli_invocation_args)
        bindings = tuple(file_bindings)
    else:
        invocation_args = tuple(strip_session_output_args(cmd_args or ()))
        session_files, bindings = collect_embedded_files_from_cli_args(
            run_result.mode,
            invocation_args,
        )

    svg_results = _read_svg_results(run_result.outputs)
    interactive_contexts = run_result.interactive_contexts
    if not interactive_contexts and (
        run_result.feature_metadata
        or run_result.biological_feature_metadata
        or run_result.orthogroup_metadata
    ):
        from gbdraw.render.interactive_svg import InteractiveSvgContext

        interactive_contexts = (
            InteractiveSvgContext(
                features=run_result.feature_metadata,
                biological_features=run_result.biological_feature_metadata,
                orthogroups=run_result.orthogroup_metadata or (),
            ),
        )
    feature_catalog = _feature_catalog_for_svg_results(
        svg_results,
        interactive_contexts,
    )
    context_output_prefix = output_prefix
    if context_output_prefix is None and len(run_result.outputs) == 1:
        context_output_prefix = run_result.outputs[0].output_prefix
    payload = build_session_json(
        SessionBuildContext(
            mode=run_result.mode,
            output_prefix=context_output_prefix,
            render_formats=run_result.render_formats,
            source_session=source_session,
            cli_invocation_args=invocation_args,
            file_bindings=tuple(bindings),
            linear_record_metadata=run_result.linear_record_metadata,
        ),
        svg_results=svg_results,
        embedded_files=session_files,
        generated_at=datetime.now(timezone.utc),
        feature_catalog=feature_catalog,
        losat_cache_entries=run_result.losat_cache_entries,
        losat_derived_cache_entries=(),
        protein_identity_manifest=run_result.protein_identity_manifest,
        legacy_protein_raw_candidates=run_result.legacy_protein_raw_candidates,
        legacy_protein_derived_evidence=run_result.legacy_protein_derived_evidence,
        canonical_request=run_result.canonical_request,
    )
    payload.pop("files", None)
    write_session_json(sidecar_path, payload, overwrite=overwrite)
    return sidecar_path


def render_canonical_session_if_present(
    session: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    output_override: str | None,
    format_override: str | None,
    save_session: bool,
    session_output: str | None,
    overwrite: bool = False,
) -> bool:
    """Render an authoritative canonical request and bypass legacy CLI replay."""

    from gbdraw.session import (
        load_session_document,
        materialize_session,
        save_session_document,
        session_to_request,
        with_request_output,
    )
    from gbdraw.session_io import CANONICAL_SESSION_MIN_VERSION

    if int(session.get("version", 0)) < CANONICAL_SESSION_MIN_VERSION:
        return False
    document = load_session_document(session)
    if document.mode != mode:
        raise ValidationError(
            f"Session renderRequest mode is {document.mode!r}; expected {mode!r}."
        )

    output_path = Path(output_override) if output_override else None
    output_directory = (
        output_path.parent if output_path is not None and output_path.parent != Path("") else Path.cwd()
    )
    with materialize_session(
        document,
        output_directory=output_directory,
    ) as materialized:
        request = session_to_request(materialized)
        request = with_request_output(
            request,
            output_prefix=output_path.name if output_path is not None else None,
            output_directory=output_directory,
            formats=format_override,
            overwrite=overwrite,
        )
        replay_prefix: str | None = None
        sidecar_path: Path | None = None
        adjunct: dict[str, Any] | None = None
        web_file_inventory: dict[str, Any] | None = None
        if save_session or session_output:
            from gbdraw.api.requests import CircularBatchRequest

            if isinstance(request, CircularBatchRequest):
                replay_prefix = (
                    output_path.name
                    if output_path is not None
                    else (
                        request.outputs[0].output_prefix
                        if len(request.outputs) == 1
                        else "gbdraw"
                    )
                )
            else:
                replay_prefix = request.output.output_prefix
            sidecar_path = (
                Path(session_output)
                if session_output
                else output_directory / f"{replay_prefix}.gbdraw-session.json"
            )
            preflight_session_sidecar_if_requested(
                save_session=True,
                session_output=str(sidecar_path),
                output_prefix=None,
                diagram_output_paths=diagram_request_output_paths(request),
                overwrite=overwrite,
            )
            adjunct, web_file_inventory = _project_session_adjunct_for_current_write(
                document.to_dict(),
                source_version=document.version,
            )
            validate_current_web_state_field_names(adjunct.get("config"))

        rendered = _render_request(
            request,
            session_document=document.to_dict(),
            include_feature_catalog=sidecar_path is not None,
        )

        if sidecar_path is not None:
            assert adjunct is not None
            protein_id_map = getattr(rendered, "protein_id_map", None) or {}
            losat_entries = getattr(rendered, "losat_cache_entries", ())
            identity_manifest = getattr(rendered, "protein_identity_manifest", None)
            legacy_raw = getattr(rendered, "legacy_protein_raw_candidates", ())
            legacy_derived = getattr(rendered, "legacy_protein_derived_evidence", ())
            if protein_id_map:
                from gbdraw.api.session_compat import (
                    rewrite_protein_artifact_references,
                )

                adjunct = rewrite_protein_artifact_references(
                    adjunct,
                    protein_id_map,
                )
            for artifact_key in (
                "losatCache",
                "losatDerivedCache",
                "proteinIdentityManifest",
                "legacyArtifacts",
            ):
                adjunct.pop(artifact_key, None)
            adjunct["losatCache"] = {
                "entries": [dict(entry) for entry in losat_entries]
            }
            adjunct["losatDerivedCache"] = {
                "entries": []
            }
            adjunct["proteinIdentityManifest"] = dict(
                identity_manifest
                or {
                    "schema": 2,
                    "proteinSets": {},
                    "recordAnalyses": {},
                    "recordInstances": {},
                }
            )
            legacy_artifacts: dict[str, Any] = {}
            if legacy_raw:
                legacy_artifacts["proteinRawCandidates"] = {
                    "schema": 1,
                    "entries": [dict(entry) for entry in legacy_raw],
                }
            if legacy_derived:
                legacy_artifacts["proteinDerivedEvidence"] = {
                    "schema": 1,
                    "entries": [dict(entry) for entry in legacy_derived],
                }
            if legacy_artifacts:
                adjunct["legacyArtifacts"] = legacy_artifacts
            svg_results = []
            for output in rendered.output_paths:
                if output.suffix.lower() != ".svg" or not output.is_file():
                    continue
                if output.name.lower().endswith(".interactive.svg"):
                    continue
                svg_results.append(
                    {"name": output.stem, "content": output.read_text(encoding="utf-8")}
                )
            if svg_results:
                adjunct["results"] = svg_results
            interactive_contexts = (
                rendered.interactive_contexts
                if hasattr(rendered, "interactive_contexts")
                else (rendered.interactive_context,)
            )
            catalog_contexts = tuple(interactive_contexts)
            catalog_results = [
                (str(result["name"]), str(result["content"]))
                for result in svg_results
            ]
            feature_catalog = _feature_catalog_for_svg_results(
                catalog_results,
                catalog_contexts,
            )
            _replace_current_derived_feature_state(
                adjunct,
                feature_catalog,
            )
            drawings = (
                rendered.drawings
                if hasattr(rendered, "drawings")
                else (rendered.drawing,)
            )
            rendered_request = rendered.request
            result_names = (
                tuple(output.output_prefix for output in rendered_request.outputs)
                if isinstance(rendered_request, CircularBatchRequest)
                else (rendered_request.output.output_prefix,)
            )
            geometry_records = [
                record
                for index, (drawing, result_name) in enumerate(
                    zip(drawings, result_names, strict=True)
                )
                for record in collect_track_slot_geometry_records(
                    drawing,
                    result_index=index,
                    result_name=str(result_name),
                )
            ]
            run_metadata = build_track_slot_geometry_run_metadata(
                mode=mode,
                records=geometry_records,
            )
            if run_metadata:
                adjunct["runMetadata"] = run_metadata
            else:
                adjunct.pop("runMetadata", None)
            save_session_document(
                sidecar_path,
                rendered_request,
                title=str(document.to_dict().get("title") or replay_prefix),
                adjunct=adjunct,
                web_file_inventory=web_file_inventory,
                overwrite=overwrite,
            )
    return True


def _web_binding_as_embedded_file(
    resources: Mapping[str, Any],
    binding: Any,
) -> dict[str, Any] | None:
    if not isinstance(binding, Mapping):
        return None
    resource_id = str(binding.get("resourceId") or "")
    resource = resources.get(resource_id)
    if not isinstance(resource, Mapping):
        if isinstance(binding.get("data"), (str, Mapping)):
            return dict(binding)
        return None
    embedded = dict(resource)
    embedded["name"] = str(binding.get("name") or resource.get("name") or "file")
    embedded["type"] = str(binding.get("type") or resource.get("type") or "")
    embedded["lastModified"] = int(
        binding.get("lastModified") or resource.get("lastModified") or 0
    )
    return embedded


def _project_web_file_inventory(
    session: Mapping[str, Any],
) -> dict[str, Any] | None:
    web_files = session.get("webFiles")
    resources = session.get("resources")
    if not isinstance(web_files, Mapping) or not isinstance(resources, Mapping):
        return None
    bindings_value = web_files.get("bindings")
    has_current_bindings = (
        isinstance(bindings_value, Mapping) and bindings_value.get("schema") == 1
    )
    bindings = bindings_value if has_current_bindings else {}
    direct_source_fields = {
        "conservationLosatFastaSources": "c_conservation_fastas",
        "conservationSequenceSources": "c_conservation_sequence_sources",
    }
    has_direct_sources = any(
        isinstance(web_files.get(field), list) for field in direct_source_fields
    )
    if not has_current_bindings and not has_direct_sources:
        return None

    original_names_value = web_files.get("resourceOriginalNames")
    original_names = (
        original_names_value if isinstance(original_names_value, Mapping) else {}
    )

    def restore(value: Any) -> Any:
        if isinstance(value, list):
            return [restore(item) for item in value]
        return _web_binding_as_embedded_file(resources, value)

    def restore_resource_id(value: Any) -> Any:
        if isinstance(value, list):
            return [restore_resource_id(item) for item in value]
        resource_id = str(value or "").strip()
        if not resource_id:
            return None
        return _web_binding_as_embedded_file(
            resources,
            {
                "resourceId": resource_id,
                "name": original_names.get(resource_id),
            },
        )

    files: dict[str, Any] = {}
    for slot in (
        "c_gb",
        "c_gff",
        "c_fasta",
        "c_depth",
        "c_conservation_blasts",
        "c_conservation_fastas",
        "c_conservation_sequence_sources",
        "d_color",
        "t_color",
        "blacklist",
        "whitelist",
        "qualifier_priority",
    ):
        if slot in bindings:
            files[slot] = restore(bindings[slot])
    files["c_conservation_blasts_source"] = (
        "losat-cache"
        if bindings.get("c_conservation_blasts_source") == "losat-cache"
        else None
    )

    linear_sequences = bindings.get("linearSeqs")
    if isinstance(linear_sequences, list):
        files["linearSeqs"] = [
            {
                **dict(sequence),
                "gb": restore(sequence.get("gb")),
                "gff": restore(sequence.get("gff")),
                "fasta": restore(sequence.get("fasta")),
                "depth": restore(sequence.get("depth")),
                "blast": restore(sequence.get("blast")),
            }
            for sequence in linear_sequences
            if isinstance(sequence, Mapping)
        ]
    linear_comparisons = bindings.get("linearComparisons")
    if isinstance(linear_comparisons, list):
        files["linearComparisons"] = [
            {**dict(comparison), "file": restore(comparison.get("file"))}
            for comparison in linear_comparisons
            if isinstance(comparison, Mapping)
        ]
    for source_field, slot in direct_source_fields.items():
        source_ids = web_files.get(source_field)
        if isinstance(source_ids, list) and slot not in files:
            files[slot] = restore_resource_id(source_ids)
    return files


def _project_session_adjunct_for_current_write(
    session: Mapping[str, Any],
    *,
    source_version: int,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """Detach non-canonical state and migrate released Web-owned field names."""

    adjunct = {
        key: value
        for key, value in session.items()
        if key
        not in {
            "format",
            "version",
            "createdAt",
            "renderRequest",
            "resources",
            "files",
        }
    }
    web_file_inventory = _project_web_file_inventory(session)
    if source_version >= CURRENT_SESSION_VERSION:
        return adjunct, web_file_inventory

    config = adjunct.get("config")
    if isinstance(config, Mapping):
        migrated_config = migrate_persisted_web_state_field_names(config)
        assert isinstance(migrated_config, Mapping)
        source_files = session.get("files")
        has_source_file_inventory = (
            isinstance(source_files, Mapping) and bool(source_files)
        ) or web_file_inventory is not None
        migrated_config, migrated_files = (
            migrate_legacy_linear_comparison_draft_for_current_writer(
                migrated_config,
                source_files
                if isinstance(source_files, Mapping)
                else (web_file_inventory or {}),
                force_web_draft=(
                    isinstance(config.get("linearRecordLayout"), Mapping)
                    or not isinstance(config.get("cliOptions"), Mapping)
                ),
            )
        )
        adjunct["config"] = migrated_config
        web_file_inventory = migrated_files if has_source_file_inventory else None
    else:
        adjunct.pop("config", None)
    if isinstance(adjunct.get("ui"), Mapping):
        adjunct["ui"] = dict(adjunct["ui"])
        adjunct["ui"].pop("blastSource", None)
    else:
        adjunct.pop("ui", None)
    web_files_value = adjunct.get("webFiles")
    if isinstance(web_files_value, Mapping):
        web_files = dict(web_files_value)
        bindings_value = web_files.get("bindings")
        if isinstance(bindings_value, Mapping):
            bindings = dict(bindings_value)
            bindings.pop("linearCanonicalComparisons", None)
            if web_file_inventory is not None:
                bindings = {}
            else:
                linear_sequences = bindings.get("linearSeqs")
                if isinstance(linear_sequences, list):
                    bindings["linearSeqs"] = [
                        {
                            key: value
                            for key, value in sequence.items()
                            if key not in {"blast", "losat_filename"}
                        }
                        if isinstance(sequence, Mapping)
                        else sequence
                        for sequence in linear_sequences
                    ]
                bindings["linearComparisons"] = []
            web_files["bindings"] = bindings
        metadata_value = web_files.get("linearRecordMetadata")
        if isinstance(metadata_value, list):
            web_files["linearRecordMetadata"] = [
                {
                    key: value
                    for key, value in metadata.items()
                    if key not in {"losatFilename", "losat_filename"}
                }
                if isinstance(metadata, Mapping)
                else metadata
                for metadata in metadata_value
            ]
        adjunct["webFiles"] = web_files
    return adjunct, web_file_inventory


def _render_request(
    request,
    *,
    session_document=None,
    include_feature_catalog: bool = False,
):
    """Import the request renderer lazily to keep CLI session imports lightweight."""

    if session_document is None:
        from gbdraw.api.request_render import render_request

        return render_request(
            request,
            include_feature_catalog=include_feature_catalog,
        )
    from gbdraw.api.session_compat import render_session_compatible_request

    return render_session_compatible_request(
        request,
        session_document,
        include_feature_catalog=include_feature_catalog,
    )


def strip_session_output_args(cmd_args: Sequence[str]) -> list[str]:
    """Remove sidecar controls and overwrite permission from saved CLI arguments."""

    result: list[str] = []
    index = 0
    while index < len(cmd_args):
        token = str(cmd_args[index])
        if token == "--save_session":
            index += 1
            continue
        if token == "--overwrite":
            index += 1
            continue
        if token.startswith("--session_output="):
            index += 1
            continue
        if token == "--session_output":
            index += 2
            continue
        result.append(token)
        index += 1
    return result


def collect_embedded_files_from_cli_args(
    mode: Literal["circular", "linear"],
    cli_args: Sequence[str],
) -> tuple[dict[str, Any], tuple[SessionFileBinding, ...]]:
    """Embed local CLI input files and build cliInvocation file bindings."""

    files = _empty_files_payload()
    bindings: list[SessionFileBinding] = []
    circular_counts: dict[str, int] = {}
    linear_depth_track_index = 0
    circular_depth_index = 0

    index = 0
    while index < len(cli_args):
        token = str(cli_args[index])
        if _is_cli_table_option(mode, token):
            value_index = index + 1
            if value_index < len(cli_args) and _is_embeddable_path(cli_args[value_index]):
                table_slot = _append_cli_input(files, cli_args[value_index], depth=False)
                bindings.append(_binding(value_index, table_slot, cli_args[value_index]))
                table_entry = {
                    "argIndex": value_index,
                    "kind": _cli_table_kind(token),
                    "slot": table_slot,
                    "dependencies": [],
                }
                for dependency in _read_cli_table_dependencies(token, cli_args[value_index]):
                    if not _is_embeddable_path(dependency.path):
                        continue
                    dependency_slot = _append_cli_input(files, dependency.path, depth=False)
                    table_entry["dependencies"].append(
                        {
                            "rowIndex": dependency.row_index,
                            "rowNumber": dependency.row_number,
                            "column": dependency.column,
                            "slot": dependency_slot,
                        }
                    )
                files.setdefault("cliTables", []).append(table_entry)
            index += 2
            continue
        if mode == "circular" and token in {"--gbk", "--gff", "--fasta", "--conservation_blast", "--depth_track"}:
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                if token == "--gbk":
                    ordinal = circular_counts.get("gbk", 0)
                    slot = "files.c_gb" if ordinal == 0 else _append_cli_input(files, value, depth=False)
                    circular_counts["gbk"] = ordinal + 1
                elif token == "--gff":
                    ordinal = circular_counts.get("gff", 0)
                    slot = "files.c_gff" if ordinal == 0 else _append_cli_input(files, value, depth=False)
                    circular_counts["gff"] = ordinal + 1
                elif token == "--fasta":
                    ordinal = circular_counts.get("fasta", 0)
                    slot = "files.c_fasta" if ordinal == 0 else _append_cli_input(files, value, depth=False)
                    circular_counts["fasta"] = ordinal + 1
                elif token == "--conservation_blast":
                    slot = f"files.c_conservation_blasts[{len(files['c_conservation_blasts'])}]"
                else:
                    slot = f"files.c_depth[{circular_depth_index}]"
                    circular_depth_index += 1
                _set_file_slot(files, slot, value, depth=token == "--depth_track")
                bindings.append(_binding(arg_index, slot, value))
            index = next_index
            continue
        if mode == "linear" and token in {"--gbk", "--gff", "--fasta", "-b", "--blast"}:
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                seq_index = offset
                if token == "--gbk":
                    slot = f"files.linearSeqs[{seq_index}].gb"
                    depth = False
                elif token == "--gff":
                    slot = f"files.linearSeqs[{seq_index}].gff"
                    depth = False
                elif token == "--fasta":
                    slot = f"files.linearSeqs[{seq_index}].fasta"
                    depth = False
                else:
                    slot = f"files.linearSeqs[{seq_index}].blast"
                    depth = False
                _set_file_slot(files, slot, value, depth=depth)
                bindings.append(_binding(arg_index, slot, value))
            index = next_index
            continue
        if mode == "linear" and token == "--depth_track":
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                slot = f"files.linearSeqs[{offset}].depth[{linear_depth_track_index}]"
                _set_file_slot(files, slot, value, depth=True)
                bindings.append(_binding(arg_index, slot, value))
            linear_depth_track_index += 1
            index = next_index
            continue
        if token in _COMMON_SINGLE_FILE_OPTIONS:
            value_index = index + 1
            if value_index < len(cli_args) and _is_embeddable_path(cli_args[value_index]):
                slot = _COMMON_SINGLE_FILE_OPTIONS[token]
                if slot == "files.cliInputs[]":
                    slot = _append_cli_input(files, cli_args[value_index], depth=False)
                else:
                    _set_file_slot(files, slot, cli_args[value_index], depth=False)
                bindings.append(_binding(value_index, slot, cli_args[value_index]))
            index += 2
            continue
        index += 1

    return files, tuple(bindings)


def _is_cli_table_option(mode: Literal["circular", "linear"], token: str) -> bool:
    if token == "--records_table":
        return True
    if mode == "circular" and token in {"--conservation_table", "--circular_track_table"}:
        return True
    if mode == "linear" and token == "--comparisons_table":
        return True
    return False


def _cli_table_kind(token: str) -> str:
    if token == "--records_table":
        return "records"
    if token == "--conservation_table":
        return "conservation"
    if token == "--circular_track_table":
        return "circular_track"
    if token == "--comparisons_table":
        return "comparisons"
    return "unknown"


def _read_cli_table_dependencies(token: str, path: object):
    if token == "--records_table":
        return read_records_table(str(path)).path_dependencies
    if token == "--conservation_table":
        return read_conservation_table(str(path)).path_dependencies
    if token == "--circular_track_table":
        return read_circular_track_table(str(path)).path_dependencies
    if token == "--comparisons_table":
        return read_comparisons_table(str(path)).path_dependencies
    return ()


_COMMON_SINGLE_FILE_OPTIONS = {
    "-d": "files.d_color",
    "--default_colors": "files.d_color",
    "-t": "files.t_color",
    "--table": "files.t_color",
    "--label_whitelist": "files.whitelist",
    "--label_blacklist": "files.blacklist",
    "--qualifier_priority": "files.qualifier_priority",
    "--label_table": "files.cliInputs[]",
    "--feature_visibility_table": "files.cliInputs[]",
}


def _empty_files_payload() -> dict[str, Any]:
    return {
        "c_gb": None,
        "c_gff": None,
        "c_fasta": None,
        "c_depth": None,
        "c_conservation_blasts": [],
        "c_conservation_fastas": [],
        "d_color": None,
        "t_color": None,
        "blacklist": None,
        "whitelist": None,
        "qualifier_priority": None,
        "linearSeqs": [],
        "cliInputs": [],
        "cliTables": [],
    }


def _read_svg_results(outputs: Sequence[RenderedSvg]) -> list[tuple[str, str]]:
    results: list[tuple[str, str]] = []
    for output in outputs:
        try:
            content = output.svg_path.read_text(encoding="utf-8")
        except OSError as exc:
            raise ValidationError(
                f"Session sidecar output cannot read generated static SVG: {output.svg_path}"
            ) from exc
        results.append((output.result_name, content))
    return results


def _collect_option_values(args: Sequence[str], start_index: int) -> tuple[list[str], int]:
    values: list[str] = []
    index = start_index
    while index < len(args):
        token = str(args[index])
        if token.startswith("-") and token.lower() not in {"-", "none", "null"}:
            break
        values.append(token)
        index += 1
    return values, index


def _is_embeddable_path(value: object) -> bool:
    text = str(value or "").strip()
    if not text or text.lower() in {"-", "none", "null"}:
        return False
    return Path(text).is_file()


def _binding(arg_index: int, slot: str, path: object) -> SessionFileBinding:
    return SessionFileBinding(
        argIndex=arg_index,
        slot=slot,
        name=safe_embedded_filename(Path(str(path)).name, fallback="file"),
    )


def _append_cli_input(files: dict[str, Any], path: object, *, depth: bool) -> str:
    slot = f"files.cliInputs[{len(files['cliInputs'])}]"
    _set_file_slot(files, slot, path, depth=depth)
    return slot


def _set_file_slot(files: dict[str, Any], slot: str, path: object, *, depth: bool) -> None:
    entry = serialize_file_entry(str(path), depth=depth)
    if slot == "files.c_gb":
        files["c_gb"] = entry
        return
    if slot == "files.c_gff":
        files["c_gff"] = entry
        return
    if slot == "files.c_fasta":
        files["c_fasta"] = entry
        return
    if slot == "files.c_depth":
        files["c_depth"] = entry
        return
    if slot == "files.d_color":
        files["d_color"] = entry
        return
    if slot == "files.t_color":
        files["t_color"] = entry
        return
    if slot == "files.blacklist":
        files["blacklist"] = entry
        return
    if slot == "files.whitelist":
        files["whitelist"] = entry
        return
    if slot == "files.qualifier_priority":
        files["qualifier_priority"] = entry
        return
    if slot.startswith("files.c_conservation_blasts["):
        index = _slot_index(slot)
        _ensure_list_size(files["c_conservation_blasts"], index + 1, None)
        files["c_conservation_blasts"][index] = entry
        return
    if slot.startswith("files.c_depth["):
        index = _slot_index(slot)
        if not isinstance(files.get("c_depth"), list):
            files["c_depth"] = []
        _ensure_list_size(files["c_depth"], index + 1, None)
        files["c_depth"][index] = entry
        return
    if slot.startswith("files.cliInputs["):
        index = _slot_index(slot)
        _ensure_list_size(files["cliInputs"], index + 1, None)
        files["cliInputs"][index] = entry
        return
    if slot.startswith("files.linearSeqs["):
        _set_linear_seq_slot(files, slot, entry)
        return
    raise ValidationError(f"Unsupported session file slot: {slot}")


def _set_linear_seq_slot(files: dict[str, Any], slot: str, entry: dict[str, Any]) -> None:
    prefix = "files.linearSeqs["
    rest = slot[len(prefix):]
    index_text, suffix = rest.split("]", 1)
    seq_index = int(index_text)
    _ensure_linear_seq(files, seq_index)
    seq = files["linearSeqs"][seq_index]
    if suffix == ".gb":
        seq["gb"] = entry
    elif suffix == ".gff":
        seq["gff"] = entry
    elif suffix == ".fasta":
        seq["fasta"] = entry
    elif suffix == ".blast":
        seq["blast"] = entry
    elif suffix == ".depth":
        seq["depth"] = entry
    elif suffix.startswith(".depth["):
        depth_index = int(suffix[len(".depth["):-1])
        if not isinstance(seq.get("depth"), list):
            seq["depth"] = []
        _ensure_list_size(seq["depth"], depth_index + 1, None)
        seq["depth"][depth_index] = entry
    else:
        raise ValidationError(f"Unsupported linear sequence file slot: {slot}")


def _ensure_linear_seq(files: dict[str, Any], index: int) -> None:
    linear_seqs = files["linearSeqs"]
    while len(linear_seqs) <= index:
        ordinal = len(linear_seqs) + 1
        linear_seqs.append(
            {
                "uid": f"cli-seq-{ordinal}",
                "gb": None,
                "gff": None,
                "fasta": None,
                "depth": None,
                "blast": None,
                "losat_gencode": 1,
                "losat_filename": "",
                "definition": "",
                "record_subtitle": "",
                "region_record_id": "",
                "region_start": None,
                "region_end": None,
                "region_reverse": False,
            }
        )


def _slot_index(slot: str) -> int:
    left = slot.rfind("[")
    right = slot.rfind("]")
    if left < 0 or right < left:
        raise ValidationError(f"Invalid session slot index: {slot}")
    return int(slot[left + 1:right])


def _ensure_list_size(values: list[Any], size: int, fill: Any) -> None:
    while len(values) < size:
        values.append(fill)


__all__ = [
    "DiagramRunResult",
    "RenderedSvg",
    "SessionCliRequest",
    "add_session_args",
    "build_track_slot_geometry_run_metadata",
    "collect_embedded_files_from_cli_args",
    "collect_track_slot_geometry_records",
    "diagram_request_output_paths",
    "diagram_request_rendered_svgs",
    "make_rendered_svg",
    "parse_session_pre_args",
    "preflight_session_sidecar_if_requested",
    "resolve_session_sidecar_path",
    "render_canonical_session_if_present",
    "save_session_sidecar_if_requested",
    "strip_session_output_args",
]
