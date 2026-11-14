#!/usr/bin/env python
# coding: utf-8

"""Shared CLI helpers for GUI session JSON input and sidecar output."""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence, TypedDict, TypeVar

from gbdraw.api.record_planning import CLI_RECORD_KEY_PREFIX
from gbdraw.cli_utils.common import add_dpi_arg, add_font_file_args, png_dpi, register_font_files, setup_logging
from gbdraw.config.models.root import read_stored_font_size
from gbdraw.core.fonts import custom_font_families
from gbdraw.render.font_embedding import remove_embedded_font_css
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    read_circular_track_table,
    read_conservation_table,
    read_comparisons_table,
    read_records_table,
)
from gbdraw.render.formats import (
    SVG_FORMAT,
    parse_output_formats,
    resolve_format_output_path,
    resolve_output_paths,
)
from gbdraw.render.output_paths import OutputExistsError, preflight_output_paths
from gbdraw.api.run_metadata import RenderedResult, build_run_metadata
from gbdraw.session_container import SessionResourceSource
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    CircularSourceCard,
    SessionBuildContext,
    SessionFileBinding,
    SessionResourceFile,
    SessionResources,
    _write_validated_session_file,
    refuse_json_session_output,
    build_session_json,
    safe_embedded_filename,
    file_resource_entry,
)
from gbdraw.session_contract import SessionRunMetadata
from gbdraw.session_migration import (
    project_session_adjunct_for_current_write,
    replace_current_derived_feature_state,
    with_current_artifacts,
)

if TYPE_CHECKING:
    from gbdraw.api.request_render import (
        CircularBatchRenderResult,
        RequestRenderResult,
    )
    from gbdraw.api.requests import DiagramRequest
    from gbdraw.render.interactive_svg import InteractiveSvgContext
    from gbdraw.session import SessionDocument, SessionPngDpi, _DrawingBuild, _SessionDrawingRender


logger = logging.getLogger(__name__)

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
    run_metadata: SessionRunMetadata | None = None
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
    font_files: tuple[str, ...] = ()
    embed_fonts: bool = True
    dpi: int | None = None


def read_stored_font_size_options(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Read each replayed ``--*font_size`` option <= 0 as Auto (OV-1282, D-PYB-7).

    The released CLI stored such a value in ``cliInvocation.args`` and drew it;
    replay keeps the option's default instead (``auto`` for one value of a
    per-track list), with the notice of the other Session readers. The options
    are the parsed ones whose destination names a font size.
    """

    for dest, value in vars(args).items():
        if not dest.endswith("font_size") or value is None:
            continue
        path = f"cliInvocation.args --{dest}"
        if isinstance(value, list):
            setattr(args, dest, [
                "auto" if read_stored_font_size(item, path=path) is None and item is not None else item
                for item in value
            ])
        elif read_stored_font_size(value, path=path) is None:
            setattr(args, dest, parser.get_default(dest))


def add_session_args(parser: argparse.ArgumentParser) -> None:
    """Add session input/output options to a diagram parser."""

    parser.add_argument(
        "--session",
        help=(
            "Regenerate a diagram from a gbdraw Session file (.gbdraw, or an "
            "older .json or .json.gz)."
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
        help="Write one GUI-loadable .gbdraw Session file for this run.",
        action="store_true",
    )
    parser.add_argument(
        "--session_output",
        metavar="PATH",
        help=(
            "Write the Session file to PATH (.gbdraw; not .json or .json.gz); "
            "implies --save_session."
        ),
        type=str,
    )


def _format_override(value: str | None) -> str | None:
    """Read a ``-f/--format`` override typed with ``--session`` (OV-1244)."""

    return None if value is None else ",".join(parse_output_formats(value))


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
    add_dpi_arg(parser)
    add_font_file_args(parser)
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
        format=_format_override(namespace.format),
        overwrite=bool(namespace.overwrite),
        save_session=bool(namespace.save_session or namespace.session_output),
        session_output=namespace.session_output,
        drawings=tuple(namespace.drawing or ()),
        font_files=tuple(namespace.font_file),
        embed_fonts=bool(namespace.embed_fonts),
        dpi=namespace.dpi,
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
        help="The gbdraw Session file (.gbdraw, or an older .json or .json.gz).",
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
            "PATH (.gbdraw; not .json or .json.gz)."
        ),
    )
    add_dpi_arg(parser)
    add_font_file_args(parser)
    namespace = parser.parse_args(list(cmd_args))
    return SessionCliRequest(
        session_path=str(namespace.session),
        output=namespace.output,
        format=_format_override(namespace.format),
        overwrite=bool(namespace.overwrite),
        save_session=bool(namespace.save_session or namespace.session_output),
        session_output=namespace.session_output,
        drawings=tuple(namespace.drawing or ()),
        list_drawings=bool(namespace.list_drawings),
        font_files=tuple(namespace.font_file),
        embed_fonts=bool(namespace.embed_fonts),
        dpi=namespace.dpi,
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
        return Path(f"{output_prefix}.gbdraw")
    if len(outputs) == 1:
        return outputs[0].svg_path.with_suffix(".gbdraw")
    return Path("gbdraw.gbdraw")


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
    refuse_json_session_output(sidecar_path)
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
        raise OutputExistsError((sidecar_path,), subject="Session output")
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
    source_session: Mapping[str, object] | None = None,
    source_resources: Mapping[str, SessionResourceSource] | None = None,
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
    payload, sources = build_session_json(
        SessionBuildContext(
            mode=run_result.mode,
            output_prefix=context_output_prefix,
            render_formats=run_result.render_formats,
            source_session=source_session,
            source_resources=source_resources,
            cli_invocation_args=invocation_args,
            file_bindings=tuple(bindings),
            run_metadata=run_result.run_metadata,
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
    _write_validated_session_file(sidecar_path, payload, sources, overwrite=overwrite)
    _warn_custom_fonts_saved_by_name(run_result.canonical_request)
    return sidecar_path


def _warn_custom_fonts_saved_by_name(request: DiagramRequest | None) -> None:
    """D-02: a Session stores a font's name, never its file."""

    if request is None:
        return
    from gbdraw.api.config import apply_config_overrides

    options = request.options
    font_family = apply_config_overrides(options.config, options.config_overrides).objects.text.font_family
    custom = custom_font_families(font_family)
    if custom:
        names = ", ".join(f"'{name}'" for name in custom)
        logger.warning(
            "The Session stores the font name %s but not the font file; a computer without "
            "that font draws the text with the nearest available font.",
            names,
        )


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
    embed_fonts: bool = True,
    dpi: int | None = None,
) -> bool:
    """Render one drawing of ``mode`` from a canonical Session.

    ``drawing`` names it by ID or name; without it the Session must have one
    drawing of ``mode``. Sessions 27-30 return ``False`` for legacy replay.
    The caller registers the font files first (``register_font_files``), so a
    Session that names its family draws with it instead of a substitute.
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
        embed_fonts=embed_fonts,
        dpi=dpi,
    )
    return True


def _session_sidecar_path(
    plans: Sequence[_SessionDrawingRender],
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
        else output_directory / f"{prefix}.gbdraw"
    )
    refuse_json_session_output(path)
    return path, prefix


def _cache_entries(cache: object) -> Sequence[Mapping[str, object]]:
    """The entries of a stored LOSAT cache (a checked Session's list), or none."""

    entries = cache.get("entries", ()) if isinstance(cache, Mapping) else ()
    return entries if isinstance(entries, (list, tuple)) else ()


def _rendered_drawing_build(
    plan: _SessionDrawingRender,
    rendered: RequestRenderResult | CircularBatchRenderResult,
    resources: SessionResources,
) -> tuple[_DrawingBuild, dict[str, object] | None]:
    """The re-save of one rendered drawing: its request, Results and artifacts.

    ``plan`` is a drawing of the current Session the re-save was upgraded to,
    so its artifacts are that drawing's own view (OV-1910).
    """

    from gbdraw.api.request_render import CircularBatchRenderResult, stored_protein_result_request
    from gbdraw.session import _DrawingBuild, require_search_frame_protein_results

    require_search_frame_protein_results(rendered.migration_report)

    # These steps replace values and never change one in place; the document
    # built from the drawing copies its data once (SessionDocument).
    state, web_file_inventory = project_session_adjunct_for_current_write(
        plan.drawing.artifacts.fields,
        resources,
        source_version=CURRENT_SESSION_VERSION,
    )
    state = with_current_artifacts(
        state,
        losat_cache_entries=rendered.losat_cache_entries,
        # A replay builds none; the Session keeps the entries it was saved with.
        losat_derived_cache_entries=_cache_entries(state.get("losatDerivedCache")),
        protein_identity_manifest=rendered.protein_identity_manifest,
        legacy_protein_raw_candidates=rendered.migration_report.protein_raw_candidates,
        legacy_protein_derived_evidence=rendered.migration_report.protein_derived_evidence,
        protein_id_map=rendered.migration_report.protein_id_map,
    )
    svg_results = [
        {"name": output.stem, "content": remove_embedded_font_css(output.read_text(encoding="utf-8"))}
        for output in rendered.output_paths
        if output.suffix.lower() == ".svg"
        and output.is_file()
        and not output.name.lower().endswith(".interactive.svg")
    ]
    if svg_results:
        state["results"] = svg_results
    if isinstance(rendered, CircularBatchRenderResult):
        interactive_contexts = rendered.interactive_contexts
        items = rendered.items
        result_names = tuple(output.output_prefix for output in rendered.request.outputs)
    else:
        interactive_contexts = (rendered.interactive_context,)
        items = (rendered,)
        result_names = (rendered.request.output.output_prefix,)
    replace_current_derived_feature_state(
        state,
        _feature_catalog_for_svg_results(
            [(str(result["name"]), str(result["content"])) for result in svg_results],
            interactive_contexts,
        ),
    )
    state["runMetadata"] = build_run_metadata(
        plan.drawing.mode,
        [RenderedResult(item, str(name)) for item, name in zip(items, result_names, strict=True)],
    )
    build = _DrawingBuild(
        mode=plan.drawing.mode,
        # The Session stores the protein results this render computed (OV-1068).
        request=stored_protein_result_request(
            rendered.request,
            None if isinstance(rendered, CircularBatchRenderResult) else rendered.linear_metadata,
        ),
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
    embed_fonts: bool = True,
    dpi: int | None = None,
) -> None:
    """Render drawings of a canonical Session and optionally save it again.

    ``drawings`` selects drawings by ID or name; ``None`` renders every drawing
    with a committed render. ``-o`` splits into the output directory and
    prefix; several drawings get ``<prefix>_<ID>``. ``dpi`` (``--dpi``) sets
    the PNG resolution of every drawing; without it a PNG is drawn at the PNG
    DPI the Session saved (``ui.downloadDpi``, D-93), else at the default
    (300), and only ``--dpi`` is noted when no PNG is written. With a re-save, the
    document, validated when it was loaded, is brought to the current version
    before any render; every diagram path and the sidecar are checked together before the
    first write; and the re-save replaces only the rendered drawings. Without
    a re-save, the resources' embedded data is dropped from ``document`` once
    its files are written: afterwards ``document`` and the drawing artifacts
    that share its resources are not valid for saving or reuse.
    """

    from gbdraw.api.request_render import (
        CircularBatchRenderResult,
        diagram_request_output_paths,
        preflight_diagram_request_outputs,
    )
    from gbdraw.cli_utils.diagram_run import log_render_notices
    from gbdraw.session import (
        _build_session_document_from_drawings,
        _plan_session_drawings,
        _release_materialized_resource_data,
        _render_session_drawing_plans,
        _session_parse_cache,
        _upgraded_session,
        _write_session_document,
        materialize_session,
    )

    save = bool(save_session or session_output)
    # The document the re-save is built from was validated when it was loaded;
    # its move to drawings is checked with the document written (OV-1815).
    # The upgrade logs a warning for each drawing whose Results it drops; the
    # render below writes new Results for the rendered drawings.
    base = _upgraded_session(document, None, checked=False).document if save else document
    output_path = Path(output_override) if output_override else None
    output_directory = (
        output_path.parent if output_path is not None and output_path.parent != Path("") else Path.cwd()
    )
    with materialize_session(base, output_directory=output_directory) as materialized:
        if not save:
            # Nothing reads the embedded data again: the files hold the bytes.
            _release_materialized_resource_data(materialized)
        # The drawings are decoded and rendered in one cache transaction, so a
        # decoded resource is keyed by its resource; the logs and the re-save
        # follow it.
        with _session_parse_cache(materialized):
            plans = _plan_session_drawings(
                materialized,
                drawings,
                output_prefix=output_path.name if output_path is not None else None,
                formats=format_override,
                overwrite=overwrite,
                embed_fonts=embed_fonts,
                dpi=dpi,
            )
            png_dpi(
                dpi,
                tuple(
                    path.suffix.lstrip(".")
                    for plan in plans
                    for path in diagram_request_output_paths(plan.request)
                ),
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
            # Each drawing is rendered and then logs its notices, as its mode
            # command logs its run: after its own figures, and before a later
            # drawing can fail.
            rendered: dict[str, RequestRenderResult | CircularBatchRenderResult] = {}
            for plan in plans:
                result = _render_session_drawing_plans(
                    materialized,
                    (plan,),
                    include_feature_catalog=save,
                )[plan.drawing.id]
                rendered[plan.drawing.id] = result
                log_render_notices(
                    result.items if isinstance(result, CircularBatchRenderResult) else (result,)
                )
        if sidecar is None:
            return
        sidecar_path, title_fallback = sidecar
        builds = []
        web_file_inventory: dict[str, object] | None = None
        for plan in plans:
            build, inventory = _rendered_drawing_build(plan, rendered[plan.drawing.id], base.resources)
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
                resources=base.resources,
            ),
            overwrite=overwrite,
        )
        for plan in plans:
            _warn_custom_fonts_saved_by_name(plan.request)


def legacy_replay_png_dpi(
    document: SessionDocument,
    session_request: SessionCliRequest,
    formats: str,
) -> SessionPngDpi:
    """The PNG DPI of a Session 27-30 replay: ``--dpi``, else the PNG DPI the Session saved.

    The saved one is passed on only when the replay writes a PNG
    (``formats``): the replay notes a DPI without a PNG (``png_dpi``), and
    the saved one is not the caller's choice.
    """

    from gbdraw.session import SessionPngDpi, session_png_dpi

    png = session_png_dpi(document._data, session_request.dpi)
    if png.from_session and "png" not in parse_output_formats(formats):
        return SessionPngDpi(None, from_session=False)
    return png


def _legacy_session_replay(mode: str) -> Any:
    """The mode command's replay of a Session 27-30 (CLI arguments)."""

    if mode == "linear":
        from gbdraw.linear import replay_legacy_session

        return replay_legacy_session
    from gbdraw.circular import replay_legacy_session as replay_circular

    return replay_circular


def render_main(cmd_args: Sequence[str]) -> None:
    """``gbdraw render``: render the drawings of a saved Session."""

    setup_logging()

    from gbdraw.session import SessionDrawingSelectionError, load_session_document
    from gbdraw.session_io import CANONICAL_SESSION_MIN_VERSION

    request = parse_render_args(cmd_args)
    # Before the Session is read, so its font names resolve to these files.
    register_font_files(request.font_files)
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
        embed_fonts=request.embed_fonts,
        dpi=request.dpi,
    )


def strip_session_output_args(cmd_args: Sequence[str]) -> list[str]:
    """Remove sidecar controls, overwrite permission, the font options and
    ``--dpi`` from saved CLI arguments.

    A Session stores neither a font file nor its path (D-02); font embedding
    and the PNG DPI are the caller's export choices, like ``--overwrite``.
    """

    result: list[str] = []
    index = 0
    while index < len(cmd_args):
        token = str(cmd_args[index])
        if token in {"--save_session", "--overwrite", "--no_font_embedding"}:
            index += 1
            continue
        if token.startswith(("--session_output=", "--font_file=", "--dpi=")):
            index += 1
            continue
        if token in {"--session_output", "--font_file", "--dpi"}:
            index += 2
            continue
        result.append(token)
        index += 1
    return result



class CliLinearSeqFiles(TypedDict):
    """One Linear record's files (``files.linearSeqs[i]``) and the fields the Web keeps beside them."""

    uid: str
    gb: SessionResourceFile | None
    gff: SessionResourceFile | None
    fasta: SessionResourceFile | None
    # One file, or one per depth track (``depth[j]``).
    depth: SessionResourceFile | list[SessionResourceFile | None] | None
    blast: SessionResourceFile | None
    losat_gencode: int
    losat_filename: str
    definition: str
    record_subtitle: str
    region_record_id: str
    region_start: int | None
    region_end: int | None
    region_reverse: bool


class CliTableDependency(TypedDict):
    """A file a CLI table row names, and the slot it is kept in."""

    rowIndex: int
    rowNumber: int
    column: str
    slot: str


class CliTable(TypedDict):
    """A CLI table option's file (``files.cliTables[i]``): where it is in the command and the files it names."""

    argIndex: int
    kind: str
    slot: str
    dependencies: list[CliTableDependency]


class CliSessionFiles(TypedDict):
    """The ``files`` of a CLI Session sidecar: each input file of the command in the
    slot the Web keeps it in, bound by ``cliInvocation`` (``SessionFileBinding.slot``)."""

    # One File card per Circular record input (``circularSources[i]``).
    circularSources: list[CircularSourceCard]
    # One file, or one per Circular depth track (``c_depth[j]``).
    c_depth: SessionResourceFile | list[SessionResourceFile | None] | None
    c_conservation_blasts: list[SessionResourceFile | None]
    c_conservation_fastas: list[SessionResourceFile | None]
    d_color: SessionResourceFile | None
    t_color: SessionResourceFile | None
    blacklist: SessionResourceFile | None
    whitelist: SessionResourceFile | None
    qualifier_priority: SessionResourceFile | None
    linearSeqs: list[CliLinearSeqFiles]
    cliInputs: list[SessionResourceFile | None]
    cliTables: list[CliTable]


def collect_embedded_files_from_cli_args(
    mode: Literal["circular", "linear"],
    cli_args: Sequence[str],
) -> tuple[CliSessionFiles, tuple[SessionFileBinding, ...]]:
    """Embed local CLI input files and build cliInvocation file bindings."""

    files = _empty_files_payload()
    bindings: list[SessionFileBinding] = []
    linear_depth_track_index = 0
    circular_depth_index = 0

    index = 0
    while index < len(cli_args):
        token = str(cli_args[index])
        if _is_cli_table_option(mode, token):
            value_index = index + 1
            if value_index < len(cli_args) and _is_embeddable_path(cli_args[value_index]):
                table_slot = _append_cli_input(files, cli_args[value_index])
                bindings.append(_binding(value_index, table_slot, cli_args[value_index]))
                table_entry: CliTable = {
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
                    dependency_slot = _append_cli_input(files, dependency.path)
                    table_entry["dependencies"].append(
                        {
                            "rowIndex": dependency.row_index,
                            "rowNumber": dependency.row_number,
                            "column": dependency.column,
                            "slot": dependency_slot,
                        }
                    )
                files["cliTables"].append(table_entry)
            index += 2
            continue
        if mode == "circular" and token in {"--gbk", "--gff", "--fasta", "--conservation_blast", "--depth_track"}:
            values, next_index = _collect_option_values(cli_args, index + 1)
            for offset, value in enumerate(values):
                arg_index = index + 1 + offset
                if not _is_embeddable_path(value):
                    continue
                if token in _CIRCULAR_SOURCE_FILE_KEYS:
                    # One File card per input, as the request has one record
                    # input per --gbk value or per --gff/--fasta pair.
                    slot = f"files.circularSources[{offset}].{_CIRCULAR_SOURCE_FILE_KEYS[token]}"
                elif token == "--conservation_blast":
                    slot = f"files.c_conservation_blasts[{len(files['c_conservation_blasts'])}]"
                else:
                    slot = f"files.c_depth[{circular_depth_index}]"
                    circular_depth_index += 1
                _set_file_slot(files, slot, value)
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
                elif token == "--gff":
                    slot = f"files.linearSeqs[{seq_index}].gff"
                elif token == "--fasta":
                    slot = f"files.linearSeqs[{seq_index}].fasta"
                else:
                    slot = f"files.linearSeqs[{seq_index}].blast"
                _set_file_slot(files, slot, value)
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
                _set_file_slot(files, slot, value)
                bindings.append(_binding(arg_index, slot, value))
            linear_depth_track_index += 1
            index = next_index
            continue
        if token in _COMMON_SINGLE_FILE_OPTIONS:
            value_index = index + 1
            if value_index < len(cli_args) and _is_embeddable_path(cli_args[value_index]):
                slot = _COMMON_SINGLE_FILE_OPTIONS[token]
                if slot == "files.cliInputs[]":
                    slot = _append_cli_input(files, cli_args[value_index])
                else:
                    _set_file_slot(files, slot, cli_args[value_index])
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
    "--qualifier_priority": "files.qualifier_priority",
    "--label_table": "files.cliInputs[]",
    "--feature_visibility_table": "files.cliInputs[]",
    "--feature_placement_table": "files.cliInputs[]",
    "--feature_override_table": "files.cliInputs[]",
}


# The file of a Circular File card that each input option fills.
_CIRCULAR_SOURCE_FILE_KEYS = {"--gbk": "gb", "--gff": "gff", "--fasta": "fasta"}


def _empty_files_payload() -> CliSessionFiles:
    return {
        "circularSources": [],
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
            # D-02: a Session's preview never carries font data.
            content = remove_embedded_font_css(output.svg_path.read_text(encoding="utf-8"))
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


def _append_cli_input(files: CliSessionFiles, path: object) -> str:
    slot = f"files.cliInputs[{len(files['cliInputs'])}]"
    _set_file_slot(files, slot, path)
    return slot


def _set_file_slot(files: CliSessionFiles, slot: str, path: object) -> SessionResourceFile:
    """Keep the file ``path`` in ``slot`` of ``files``; the kept entry."""

    entry = file_resource_entry(str(path))
    _place_file(files, slot, entry)
    return entry


def _place_file(files: CliSessionFiles, slot: str, entry: SessionResourceFile) -> None:
    if slot.startswith("files.circularSources["):
        _set_circular_source_slot(files, slot, entry)
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
    if slot == "files.whitelist":
        files["whitelist"] = entry
        return
    if slot == "files.qualifier_priority":
        files["qualifier_priority"] = entry
        return
    if slot.startswith("files.c_conservation_blasts["):
        _place_at(files["c_conservation_blasts"], _slot_index(slot), entry)
        return
    if slot.startswith("files.c_depth["):
        depths = files["c_depth"]
        if not isinstance(depths, list):
            depths = []
            files["c_depth"] = depths
        _place_at(depths, _slot_index(slot), entry)
        return
    if slot.startswith("files.cliInputs["):
        _place_at(files["cliInputs"], _slot_index(slot), entry)
        return
    if slot.startswith("files.linearSeqs["):
        _set_linear_seq_slot(files, slot, entry)
        return
    raise ValidationError(f"Unsupported session file slot: {slot}")


def _set_linear_seq_slot(files: CliSessionFiles, slot: str, entry: SessionResourceFile) -> None:
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
        depths = seq["depth"]
        if not isinstance(depths, list):
            depths = []
            seq["depth"] = depths
        _place_at(depths, int(suffix[len(".depth["):-1]), entry)
    else:
        raise ValidationError(f"Unsupported linear sequence file slot: {slot}")


def _set_circular_source_slot(files: CliSessionFiles, slot: str, entry: SessionResourceFile) -> None:
    """Fill ``files.circularSources[<i>].<gb|gff|fasta>``; card ``i`` is the
    request's record input ``i``, so its uid is that input's CLI record key."""

    index_text, suffix = slot[len("files.circularSources["):].split("]", 1)
    index = int(index_text)
    cards = files["circularSources"]
    while len(cards) <= index:
        cards.append(CircularSourceCard(uid=f"{CLI_RECORD_KEY_PREFIX}{len(cards) + 1}", gb=None, gff=None, fasta=None))
    card = cards[index]
    if suffix == ".gb":
        card["gb"] = entry
    elif suffix == ".gff":
        card["gff"] = entry
    elif suffix == ".fasta":
        card["fasta"] = entry
    else:
        raise ValidationError(
            f"Unsupported session file slot: {slot}",
            diagnostic={"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"},
        )


def _ensure_linear_seq(files: CliSessionFiles, index: int) -> None:
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


_Item = TypeVar("_Item")


def _place_at(values: list[_Item | None], index: int, item: _Item) -> None:
    """Put ``item`` at ``index`` of ``values``, filling the slots before it with None."""

    while len(values) <= index:
        values.append(None)
    values[index] = item


__all__ = [
    "DiagramRunResult",
    "RenderedSvg",
    "SessionCliRequest",
    "add_session_args",
    "read_stored_font_size_options",
    "collect_embedded_files_from_cli_args",
    "diagram_request_rendered_svgs",
    "make_rendered_svg",
    "parse_session_pre_args",
    "preflight_session_sidecar_if_requested",
    "resolve_session_sidecar_path",
    "render_canonical_session_if_present",
    "save_session_sidecar_if_requested",
    "strip_session_output_args",
]
