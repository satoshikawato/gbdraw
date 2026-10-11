#!/usr/bin/env python
# coding: utf-8

"""Shared CLI helpers for GUI session JSON input and sidecar output."""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field, replace
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
    SessionBuildContext,
    SessionFileBinding,
    _write_validated_session_json,
    build_session_json,
    get_session_slot,
    safe_embedded_filename,
    serialize_file_entry,
)
from gbdraw.session_migration import (
    project_session_adjunct_for_current_write,
    replace_current_derived_feature_state,
    with_current_artifacts,
)

if TYPE_CHECKING:
    from gbdraw.api.requests import DiagramRequest
    from gbdraw.render.interactive_svg import InteractiveSvgContext
    from gbdraw.session import SessionDocument


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
    drawings: tuple[str, ...] = ()
    list_drawings: bool = False


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
        "--drawing",
        metavar="ID",
        help=(
            "With --session, the drawing to render, by ID or name (default: "
            "the Session's only drawing of this mode)."
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
        if any(
            str(token) == "--drawing" or str(token).startswith("--drawing=")
            for token in cmd_args
        ):
            raise ValidationError(
                "--drawing selects a drawing of a --session file.",
                diagnostic={"code": "INPUT_INVALID", "field": "input", "reason": "REQUIRED"},
            )
        return None

    parser = argparse.ArgumentParser(
        prog=f"gbdraw {mode}",
        add_help=False,
    )
    parser.add_argument("--session", required=True)
    parser.add_argument("--drawing", action="append")
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
    if namespace.drawing and len(namespace.drawing) > 1:
        parser.error(
            f"gbdraw {mode} renders one drawing; use gbdraw render for several."
        )
    return SessionCliRequest(
        session_path=str(namespace.session),
        output=namespace.output,
        format=namespace.format,
        overwrite=bool(namespace.overwrite),
        save_session=bool(namespace.save_session or namespace.session_output),
        session_output=namespace.session_output,
        drawings=tuple(namespace.drawing or ()),
    )


def parse_render_args(cmd_args: Sequence[str]) -> SessionCliRequest:
    """Parse ``gbdraw render``: render the drawings of a saved Session."""

    parser = argparse.ArgumentParser(
        prog="gbdraw render",
        description=(
            "Render the drawings of a plain or gzip-compressed gbdraw Session "
            "file. Without --drawing, every drawing with a committed render "
            "is rendered and the others are skipped with a notice."
        ),
    )
    parser.add_argument(
        "--session",
        required=True,
        metavar="FILE",
        help="The gbdraw Session file (.gbdraw-session.json or .json.gz).",
    )
    parser.add_argument(
        "--drawing",
        action="extend",
        nargs="+",
        metavar="ID",
        help=(
            "Render only these drawings, by ID or name; repeatable. Naming a "
            "drawing without a committed render is an error."
        ),
    )
    parser.add_argument(
        "--list_drawings",
        action="store_true",
        help=(
            "Print one tab-separated line per drawing (ID, mode, name, and yes "
            "or no for a committed render) and exit."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output path prefix (default: each drawing's saved prefix). With "
            "several drawings, each name ends in _<ID>."
        ),
    )
    parser.add_argument(
        "-f",
        "--format",
        help=(
            "Comma-separated list of output file formats (svg, interactive_svg, "
            "png, pdf, eps, ps; default: the saved formats; png/pdf/eps/ps "
            "require CairoSVG)."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output files (default: refuse to overwrite).",
    )
    sidecar = parser.add_mutually_exclusive_group()
    sidecar.add_argument(
        "--save_session",
        action="store_true",
        help=(
            "Write the Session again with the rendered drawings replaced, next "
            "to the diagrams; several drawings need -o."
        ),
    )
    sidecar.add_argument(
        "--session_output",
        metavar="PATH",
        help=(
            "Write the Session again with the rendered drawings replaced to "
            "PATH; use a .gz suffix for gzip compression."
        ),
    )
    namespace = parser.parse_args(list(cmd_args))
    return SessionCliRequest(
        session_path=str(namespace.session),
        output=namespace.output,
        format=namespace.format,
        overwrite=bool(namespace.overwrite),
        save_session=bool(namespace.save_session or namespace.session_output),
        session_output=namespace.session_output,
        drawings=tuple(namespace.drawing or ()),
        list_drawings=bool(namespace.list_drawings),
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
    from gbdraw.api.request_render import diagram_request_output_paths

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
        _canonical_request_is_resolved=True,
    )
    # build_session_json validated the payload it returned (without "files").
    _write_validated_session_json(sidecar_path, payload, overwrite=overwrite)
    return sidecar_path


def _require_committed_render(document: SessionDocument) -> None:
    if not any(drawing.has_canonical_request for drawing in document.drawings):
        raise ValidationError("Settings-only Session has no biological render request; load a source in Web before generating.")


def render_canonical_session_if_present(
    session: SessionDocument | Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    output_override: str | None,
    format_override: str | None,
    save_session: bool,
    session_output: str | None,
    overwrite: bool = False,
    drawing: str | None = None,
) -> bool:
    """Render one drawing of ``mode`` from a canonical Session.

    ``drawing`` names it by ID or name; without it the Session must have one
    drawing of ``mode``. Sessions 27-30 return ``False`` for legacy replay.
    """

    from gbdraw.session import load_session_document
    from gbdraw.session_io import CANONICAL_SESSION_MIN_VERSION

    document = load_session_document(session)
    if document.version < CANONICAL_SESSION_MIN_VERSION:
        return False
    _require_committed_render(document)
    selected = document.drawing(drawing, mode=mode)
    render_session_drawings_cli(
        document,
        drawings=(selected.id,),
        output_override=output_override,
        format_override=format_override,
        save_session=save_session,
        session_output=session_output,
        overwrite=overwrite,
    )
    return True


def _session_sidecar_path(
    plans: Sequence[Any],
    *,
    session_output: str | None,
    output_path: Path | None,
    output_directory: Path,
) -> tuple[Path, str]:
    """The sidecar path and the title fallback of a Session re-save."""

    from gbdraw.api.requests import CircularBatchRequest

    if len(plans) > 1:
        if output_path is None and not session_output:
            raise ValidationError(
                "--save_session with several drawings needs -o or --session_output.",
                diagnostic={"code": "INPUT_INVALID", "field": "output_prefix", "reason": "REQUIRED"},
            )
        prefix = output_path.name if output_path is not None else Path(str(session_output)).name
    else:
        request = plans[0].request
        if isinstance(request, CircularBatchRequest):
            prefix = (
                output_path.name
                if output_path is not None
                else (
                    request.outputs[0].output_prefix
                    if len(request.outputs) == 1
                    else "gbdraw"
                )
            )
        else:
            prefix = request.output.output_prefix
    path = (
        Path(session_output)
        if session_output
        else output_directory / f"{prefix}.gbdraw-session.json"
    )
    return path, prefix


def _rendered_drawing_build(
    plan: Any,
    rendered: Any,
    *,
    source_version: int,
) -> tuple[Any, dict[str, Any] | None]:
    """The re-save of one rendered drawing: its request, Results and artifacts."""

    from gbdraw.api.requests import CircularBatchRequest
    from gbdraw.session import _DrawingBuild

    state, web_file_inventory = project_session_adjunct_for_current_write(
        copy.deepcopy(dict(plan.drawing.artifacts.fields)),
        source_version=source_version,
    )
    state = with_current_artifacts(
        state,
        losat_cache_entries=getattr(rendered, "losat_cache_entries", ()),
        protein_identity_manifest=getattr(rendered, "protein_identity_manifest", None),
        legacy_protein_raw_candidates=getattr(rendered, "legacy_protein_raw_candidates", ()),
        legacy_protein_derived_evidence=getattr(rendered, "legacy_protein_derived_evidence", ()),
        protein_id_map=getattr(rendered, "protein_id_map", None),
    )
    svg_results = [
        {"name": output.stem, "content": output.read_text(encoding="utf-8")}
        for output in rendered.output_paths
        if output.suffix.lower() == ".svg"
        and output.is_file()
        and not output.name.lower().endswith(".interactive.svg")
    ]
    if svg_results:
        state["results"] = svg_results
    interactive_contexts = (
        rendered.interactive_contexts
        if hasattr(rendered, "interactive_contexts")
        else (rendered.interactive_context,)
    )
    replace_current_derived_feature_state(
        state,
        _feature_catalog_for_svg_results(
            [(str(result["name"]), str(result["content"])) for result in svg_results],
            tuple(interactive_contexts),
        ),
    )
    request = rendered.request
    svg_drawings = (
        rendered.drawings if hasattr(rendered, "drawings") else (rendered.drawing,)
    )
    result_names = (
        tuple(output.output_prefix for output in request.outputs)
        if isinstance(request, CircularBatchRequest)
        else (request.output.output_prefix,)
    )
    run_metadata = build_track_slot_geometry_run_metadata(
        mode=plan.drawing.mode,
        records=[
            record
            for index, (svg_drawing, result_name) in enumerate(
                zip(svg_drawings, result_names, strict=True)
            )
            for record in collect_track_slot_geometry_records(
                svg_drawing,
                result_index=index,
                result_name=str(result_name),
            )
        ],
    )
    if run_metadata:
        state["runMetadata"] = run_metadata
    else:
        state.pop("runMetadata", None)
    build = _DrawingBuild(
        mode=plan.drawing.mode,
        request=request,
        state=state,
        id=plan.drawing.id,
    )
    return build, web_file_inventory


def render_session_drawings_cli(
    document: SessionDocument,
    *,
    drawings: Sequence[str] | None,
    output_override: str | None,
    format_override: str | None,
    save_session: bool,
    session_output: str | None,
    overwrite: bool = False,
) -> None:
    """Render drawings of a canonical Session and optionally save it again.

    ``drawings`` selects drawings by ID or name; ``None`` renders every drawing
    with a committed render. ``-o`` splits into the output directory and
    prefix; several drawings get ``<prefix>_<ID>``. With a re-save, the
    document is brought to the current version and validated before any
    render; every diagram path and the sidecar are checked together before the
    first write; and the re-save replaces only the rendered drawings.
    """

    from gbdraw.api.request_render import (
        CircularBatchRenderResult,
        diagram_request_output_paths,
        preflight_diagram_request_outputs,
    )
    from gbdraw.features.overrides import log_feature_identity_notices
    from gbdraw.session import (
        _build_session_document_from_drawings,
        _plan_session_drawings,
        _render_session_drawing_plans,
        _write_session_document,
        materialize_session,
        upgrade_session_document,
    )

    save = bool(save_session or session_output)
    # The document the re-save is built from is validated before any render.
    # The upgrade logs a warning for each drawing whose Results it drops; the
    # render below writes new Results for the rendered drawings.
    base = upgrade_session_document(document).document if save else document
    output_path = Path(output_override) if output_override else None
    output_directory = (
        output_path.parent if output_path is not None and output_path.parent != Path("") else Path.cwd()
    )
    with materialize_session(base, output_directory=output_directory) as materialized:
        plans = _plan_session_drawings(
            materialized,
            drawings,
            output_prefix=output_path.name if output_path is not None else None,
            formats=format_override,
            overwrite=overwrite,
        )
        sidecar: tuple[Path, str] | None = None
        if save:
            sidecar = _session_sidecar_path(
                plans,
                session_output=session_output,
                output_path=output_path,
                output_directory=output_directory,
            )
            preflight_session_sidecar_if_requested(
                save_session=True,
                session_output=str(sidecar[0]),
                output_prefix=None,
                diagram_output_paths=tuple(
                    path
                    for plan in plans
                    for path in diagram_request_output_paths(plan.request)
                ),
                overwrite=overwrite,
            )
        preflight_diagram_request_outputs(tuple(plan.request for plan in plans))
        rendered = _render_session_drawing_plans(
            materialized,
            plans,
            include_feature_catalog=save,
        )
        for result in rendered.values():
            items = (
                result.items
                if isinstance(result, CircularBatchRenderResult)
                else (result,)
            )
            for item in items:
                log_feature_identity_notices(item.feature_identity_notices)
        if sidecar is None:
            return
        sidecar_path, title_fallback = sidecar
        builds = []
        web_file_inventory: dict[str, Any] | None = None
        for plan in plans:
            build, inventory = _rendered_drawing_build(
                plan,
                rendered[plan.drawing.id],
                source_version=base.version,
            )
            builds.append(build)
            # The drawings of a current Session share their Web file inventory.
            web_file_inventory = web_file_inventory or inventory
        _write_session_document(
            sidecar_path,
            _build_session_document_from_drawings(
                builds,
                base=base,
                title=str(plans[0].drawing.artifacts.fields.get("title") or title_fallback),
                web_file_inventory=web_file_inventory,
                resources=plans[0].drawing.artifacts.fields["resources"],
            ),
            overwrite=overwrite,
        )


def _legacy_session_replay(mode: str) -> Any:
    """The mode command's replay of a Session 27-30 (CLI arguments)."""

    if mode == "linear":
        from gbdraw.linear import replay_legacy_session

        return replay_legacy_session
    from gbdraw.circular import replay_legacy_session as replay_circular

    return replay_circular


def render_main(cmd_args: Sequence[str]) -> None:
    """``gbdraw render``: render the drawings of a saved Session."""

    from gbdraw.session import SessionDrawingSelectionError, load_session_document
    from gbdraw.session_io import CANONICAL_SESSION_MIN_VERSION

    request = parse_render_args(cmd_args)
    document = load_session_document(request.session_path)
    if request.list_drawings:
        for drawing in document.drawings:
            print(
                "\t".join(
                    (
                        drawing.id,
                        drawing.mode,
                        drawing.name,
                        "yes" if drawing.has_canonical_request else "no",
                    )
                )
            )
        return
    if document.version < CANONICAL_SESSION_MIN_VERSION:
        # A Session 27-30 is one drawing that its mode command replays.
        if len(request.drawings) > 1:
            raise SessionDrawingSelectionError(
                f"Session version {document.version} has one drawing; name it once."
            )
        selected = document.drawing(request.drawings[0] if request.drawings else None)
        _legacy_session_replay(selected.mode)(document, request)
        return
    _require_committed_render(document)
    render_session_drawings_cli(
        document,
        drawings=request.drawings or None,
        output_override=request.output,
        format_override=request.format,
        save_session=request.save_session,
        session_output=request.session_output,
        overwrite=request.overwrite,
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
    circular_genbank_bindings: list[int] = []
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
                table_entry: dict[str, Any] = {
                    "argIndex": value_index,
                    "kind": _cli_table_kind(token),
                    "slot": table_slot,
                    "dependencies": [],
                }
                for dependency in _read_cli_table_dependencies(
                    token,
                    cli_args[value_index],
                    ring_losat=mode == "circular" and _has_cli_losat(cli_args),
                ):
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
                if token == "--gbk":
                    circular_genbank_bindings.append(len(bindings))
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

    if len(circular_genbank_bindings) > 1:
        components = [
            get_session_slot({"files": files}, bindings[index].slot)
            for index in circular_genbank_bindings
        ]
        files["c_gb"] = {
            "kind": "composite",
            "components": components,
            **{key: components[0][key] for key in ("name", "type", "lastModified")},
        }
        for component_index, binding_index in enumerate(circular_genbank_bindings):
            bindings[binding_index] = replace(
                bindings[binding_index], slot=f"files.c_gb.components[{component_index}]"
            )
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


def _has_cli_losat(cli_args) -> bool:
    return any(
        str(token) == "--losat" or str(token).startswith("--losat=") for token in cli_args
    )


def _read_cli_table_dependencies(token: str, path: object, *, ring_losat: bool = False):
    if token == "--records_table":
        return read_records_table(str(path)).path_dependencies
    if token == "--conservation_table":
        return read_conservation_table(str(path), losat=ring_losat).path_dependencies
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
    "--feature_placement_table": "files.cliInputs[]",
    "--feature_override_table": "files.cliInputs[]",
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
    "diagram_request_rendered_svgs",
    "make_rendered_svg",
    "parse_session_pre_args",
    "preflight_session_sidecar_if_requested",
    "resolve_session_sidecar_path",
    "render_canonical_session_if_present",
    "save_session_sidecar_if_requested",
    "strip_session_output_args",
]
