"""Canonical typed-request bridge used by the browser render worker."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import re
import shutil
from time import perf_counter
from typing import AbstractSet, Any, Iterator, Literal, Mapping, MutableMapping, TypedDict, cast

from gbdraw.api.request_render import (
    CircularBatchRenderResult,
    RequestRenderResult,
    _capture_request_render_diagnostics,
    _without_derived_losat_entries,
    render_request,
)
from gbdraw.api.requests import DiagramRequest, LinearDiagramRequest
from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
)
from gbdraw.exceptions import ExportError, LabelsDoNotFitError, ValidationError
from gbdraw.io.colors import rule_captions_as_written
from gbdraw.web_support.error_adapter import WebErrorAttributes, web_error_stage
from gbdraw.legend.row_facts import collect_legend_row_facts
from gbdraw.render.document import ContentPiece, DocumentParts, document_parts
from gbdraw.render.formats import SVG_FORMAT, resolve_format_output_path
from gbdraw.render.track_slot_metadata import (
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
)
from gbdraw.render_stages import render_stage_revision
from gbdraw.session import with_request_output
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.session_io import materialize_embedded_file, safe_embedded_filename
from gbdraw.web_support.feature_catalog import (
    build_feature_catalog,
    build_feature_catalog_item,
)
from gbdraw.web_support.label_fit import LabelFitMode, LabelFitRowSearches, label_fit_report_payload

# OD-F14: Generate on the Web draws crossing label leaders with their counts
# only; the dialog computes the rows when the user asks for them.
WEB_LABEL_FIT_WARN_SEARCH_SECONDS = 0.0


class _WebContentPieceFields(TypedDict):
    kind: Literal["content"]
    start: str
    revision: str
    bytes: int


class WebContentPiece(_WebContentPieceFields, total=False):
    """A placed content group: ``children`` only when the page does not hold ``revision``."""

    children: str


class WebRawPiece(TypedDict):
    kind: Literal["raw"]
    text: str


class WebDocument(TypedDict):
    """The root start tag, the top-level pieces in order, and the root end tag."""

    head: str
    tail: str
    pieces: list[WebContentPiece | WebRawPiece]


class WebRenderResult(TypedDict):
    """One SVG Result: the written file text, and the same text as pieces for the Worker."""

    name: str
    content: str
    document: WebDocument


class WebRenderReply(TypedDict):
    """One Web render: the SVG Results, the run metadata, and the feature catalog beside them.

    ``featureCatalog`` is ``None`` only when ``featureCatalogRevision`` equals
    the revision the page holds: the page keeps the catalog it admitted for it
    (gbdraw/web/CLAUDE.md, Runtime data flow). A render without a prepared
    input cache has no revision and always returns the catalog.
    """

    results: list[WebRenderResult]
    metadata: dict[str, Any]
    featureCatalog: dict[str, object] | None
    featureCatalogRevision: str | None


def _label_fit_mode_of(request: DiagramRequest) -> LabelFitMode:
    return "linear" if isinstance(request, LinearDiagramRequest) else "circular"


_RESOURCE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")
_RESOURCE_CACHE_TOKEN_RE = re.compile(r"^render-resource-[1-9][0-9]*$")
_STAGED_WORKSPACE_RE = re.compile(r"^gbdraw-web-render-[1-9][0-9]*$")
_STAGED_WORKSPACE_MARKER = ".gbdraw-worker-render-workspace"


@contextmanager
def _web_render_diagnostic_phase(
    diagnostics: MutableMapping[str, Any] | None,
    name: str,
) -> Iterator[None]:
    stage = {"decode": "request-validation", "renderRequest": "render"}.get(name, "result-admission")
    started_at = perf_counter()
    try:
        with web_error_stage(stage):
            yield
    finally:
        if diagnostics is not None:
            timings = diagnostics.setdefault("timingsMs", {})
            timings[name] = float(timings.get(name, 0.0)) + (
                perf_counter() - started_at
            ) * 1000.0


def _attach_exception_note(error: BaseException, note: str) -> None:
    """Attach a diagnostic on Python versions before BaseException.add_note."""

    cast(WebErrorAttributes, error)._web_error_secondary = [{"code": "CLEANUP_FAILED", "stage": "cleanup"}]
    add_note = getattr(error, "add_note", None)
    if callable(add_note):
        add_note(note)
        return
    notes = getattr(error, "__notes__", None)
    if isinstance(notes, list):
        notes.append(note)
    else:
        cast(WebErrorAttributes, error).__notes__ = [note]


def _base_svg_path(item: RequestRenderResult) -> Path:
    output = item.request.output
    base_prefix = Path(output.output_directory or ".") / output.output_prefix
    expected = Path(resolve_format_output_path(str(base_prefix), SVG_FORMAT))
    expected_identity = expected.resolve(strict=False)
    for output_path in item.output_paths:
        path = Path(output_path)
        if path.resolve(strict=False) == expected_identity:
            return path
    raise ValidationError(
        "The canonical Web request did not produce its expected base SVG: "
        f"{expected.name}."
    )


def _web_document(
    parts: DocumentParts,
    file_text: str,
    held_segments: AbstractSet[str],
) -> WebDocument:
    """The reply form of ``parts``: a held revision is named without its text."""

    pieces: list[WebContentPiece | WebRawPiece] = []
    length = len(parts.prolog) + len(parts.head) + len(parts.tail)
    for piece in parts.pieces:
        if isinstance(piece, ContentPiece):
            content: WebContentPiece = {
                "kind": "content",
                "start": piece.start,
                "revision": piece.revision,
                "bytes": piece.byte_length,
            }
            if piece.revision not in held_segments:
                content["children"] = piece.children
            pieces.append(content)
            length += len(piece.start) + len(piece.children) + len("</g>")
        else:
            pieces.append({"kind": "raw", "text": piece.text})
            length += len(piece.text)
    if length != len(file_text):
        raise ExportError("The Web document pieces do not match the written SVG.")
    return {"head": parts.head, "tail": parts.tail, "pieces": pieces}


def render_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
    _diagnostics: MutableMapping[str, Any] | None = None,
    _label_fit_searches: LabelFitRowSearches | None = None,
) -> WebRenderReply:
    """Decode and render one browser request without a CLI translation layer."""

    return _render_canonical_web_request_with_prepared_inputs(
        payload,
        resource_paths=resource_paths,
        output_directory=output_directory,
        _diagnostics=_diagnostics,
        _prepared_input_cache=None,
        _resource_identities=None,
        _label_fit_searches=_label_fit_searches,
        _held_feature_catalog_revision=None,
        _held_segments=frozenset(),
    )


def _render_canonical_web_request_with_prepared_inputs(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
    _diagnostics: MutableMapping[str, Any] | None,
    _prepared_input_cache: PreparedBiologicalInputCache | None,
    _resource_identities: Mapping[str, PreparedResourceIdentity] | None,
    _label_fit_searches: LabelFitRowSearches | None,
    _held_feature_catalog_revision: str | None,
    _held_segments: AbstractSet[str],
) -> WebRenderReply:
    def render() -> WebRenderReply:
        # The request's color table is the app's rule list, not a saved file.
        # No reply field reads the derived LOSAT entries; the page builds its own.
        with (
            _capture_request_render_diagnostics(_diagnostics),
            rule_captions_as_written(),
            _without_derived_losat_entries(),
        ):
            return _render_canonical_web_request(
                payload,
                resource_paths=resource_paths,
                output_directory=output_directory,
                diagnostics=_diagnostics,
                label_fit_searches=_label_fit_searches,
                held_feature_catalog_revision=_held_feature_catalog_revision,
                held_segments=_held_segments,
            )

    if _prepared_input_cache is None:
        if _resource_identities is not None:
            raise ValidationError(
                "Prepared resource identities require a prepared input cache."
            )
        return render()
    if _resource_identities is None:
        raise ValidationError(
            "A prepared input cache requires validated resource identities."
        )
    cache_paths = {
        resource_paths[resource_id]: identity
        for resource_id, identity in _resource_identities.items()
    }
    with _prepared_input_cache.transaction(
        resource_paths=cache_paths,
        diagnostics=_diagnostics,
    ):
        return render()


def _render_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
    diagnostics: MutableMapping[str, Any] | None,
    label_fit_searches: LabelFitRowSearches | None,
    held_feature_catalog_revision: str | None,
    held_segments: AbstractSet[str],
) -> WebRenderReply:

    if not isinstance(payload, Mapping):
        raise ValidationError("The Web render request must be an object.")
    if not isinstance(resource_paths, Mapping):
        raise ValidationError("The Web render resource map must be an object.")
    output_root = Path(output_directory)
    output_root.mkdir(parents=True, exist_ok=True)
    with _web_render_diagnostic_phase(diagnostics, "decode"):
        # The preview is not an export: the page shows the fonts it loaded, and
        # a Session keeps the preview without font data (D-02). The Web export
        # adds the fonts itself (gbdraw.render.font_embedding).
        request = with_request_output(
            decode_canonical_request(
                payload,
                resource_paths=resource_paths,
                output_directory=output_root,
            ),
            embed_fonts=False,
        )
    with _web_render_diagnostic_phase(diagnostics, "renderRequest"):
        try:
            rendered = render_request(
                request,
                include_feature_catalog=True,
                label_fit_warn_search_seconds=WEB_LABEL_FIT_WARN_SEARCH_SECONDS,
            )
        except LabelsDoNotFitError as error:
            # Nothing was drawn: keep the failure for the rows its bounded search did not reach.
            if label_fit_searches is not None:
                label_fit_searches.retain_failure(_label_fit_mode_of(request), error)
            raise
    items: tuple[RequestRenderResult, ...]
    if isinstance(rendered, CircularBatchRenderResult):
        items = rendered.items
    else:
        items = (rendered,)

    results: list[WebRenderResult] = []
    geometry_records: list[dict[str, Any]] = []
    feature_catalog_items: list[dict[str, Any]] = []
    legend_rows: list[dict[str, Any]] = []
    for result_index, item in enumerate(items):
        path = _base_svg_path(item)
        result_name = path.name
        with _web_render_diagnostic_phase(diagnostics, "svgReadback"):
            svg_content = path.read_text(encoding="utf-8")
            document = _web_document(document_parts(item.drawing), svg_content, held_segments)
        results.append(
            {
                "name": result_name,
                "content": svg_content,
                "document": document,
            }
        )
        with _web_render_diagnostic_phase(diagnostics, "featureCatalog"):
            feature_catalog_items.append(
                build_feature_catalog_item(
                    svg_content,
                    item.interactive_context,
                    result_index=result_index,
                    result_name=result_name,
                    drawing=item.drawing,
                    _diagnostics=diagnostics,
                )
            )
        with _web_render_diagnostic_phase(diagnostics, "geometryMetadata"):
            legend_rows.append(
                collect_legend_row_facts(
                    item.drawing,
                    result_index=result_index,
                    result_name=result_name,
                )
            )
            geometry_records.extend(
                collect_track_slot_geometry_records(
                    item.drawing,
                    result_index=result_index,
                    result_name=result_name,
                )
            )
    if not results:
        with web_error_stage("result-admission"):
            raise ValidationError("The canonical Web request did not produce an SVG.")
    with _web_render_diagnostic_phase(diagnostics, "geometryMetadata"):
        metadata = build_track_slot_geometry_run_metadata(
            mode=rendered.mode,
            records=geometry_records,
        )
    with _web_render_diagnostic_phase(diagnostics, "featureCatalog"):
        # The page keeps the catalog it holds when this render reused the same
        # item builds in this cache (gbdraw/web/CLAUDE.md, Runtime data flow).
        catalog_revision = render_stage_revision("catalog_item", feature_catalog_items)
        feature_catalog = (
            None
            if catalog_revision is not None and catalog_revision == held_feature_catalog_revision
            else build_feature_catalog(feature_catalog_items)
        )
    # Per Result: the Legend rows drawn and the feature rows the draft removed (OV-63).
    metadata["legendRows"] = legend_rows
    metadata["annotationWarnings"] = [
        {
            "code": warning.code,
            "setId": warning.set_id,
            "annotationId": warning.annotation_id,
            "recordId": warning.record_id,
            "recordIndex": warning.record_index,
            "missingCount": warning.missing_count,
            "message": warning.message,
            "resultIndex": result_index,
            "resultName": results[result_index]["name"],
        }
        for result_index, item in enumerate(items)
        for warning in item.annotation_warnings
    ]
    metadata["featureIdentityNotices"] = [
        {
            "recordKey": notice.record_key,
            "biologicalFeatureId": notice.biological_feature_id,
            "status": notice.status,
            "kinds": list(notice.kinds),
            "resultIndex": result_index,
        }
        for result_index, item in enumerate(items)
        for notice in item.feature_identity_notices
    ]
    # Browser execution discards logging, so table-ID warnings travel as metadata.
    metadata["comparisonWarnings"] = [
        {
            "code": warning.code,
            "queryRecordIndex": warning.query_record_index,
            "subjectRecordIndex": warning.subject_record_index,
            "queryRecordId": warning.query_record_id,
            "subjectRecordId": warning.subject_record_id,
            "rowCount": warning.row_count,
            "exampleIds": list(warning.example_ids),
            "message": warning.message,
            "resultIndex": result_index,
            "resultName": results[result_index]["name"],
        }
        for result_index, item in enumerate(items)
        if item.linear_metadata is not None
        for warning in item.linear_metadata.comparison_warnings
    ]
    # Crossing label leaders drawn anyway: the counts, and a token that names
    # the figure whose rows the dialog computes later (OD-F14).
    row_tokens = (
        label_fit_searches.retain(rendered.mode, [item.drawing for item in items])
        if label_fit_searches is not None
        else [None] * len(items)
    )
    metadata["labelFitReports"] = [
        {
            "resultIndex": result_index,
            "resultName": results[result_index]["name"],
            "report": label_fit_report_payload(report, rendered.mode),
            "rowToken": row_tokens[result_index],
        }
        for result_index, item in enumerate(items)
        if (report := item.label_fit_report) is not None
    ]
    return {
        "results": results,
        "metadata": metadata,
        "featureCatalog": feature_catalog,
        "featureCatalogRevision": catalog_revision,
    }


def render_embedded_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resources: Mapping[str, Mapping[str, Any]],
    workspace: str | Path,
) -> WebRenderReply:
    """Materialize browser-owned resources, then render their typed request."""

    if not isinstance(resources, Mapping):
        raise ValidationError("The Web render resources must be an object.")
    workspace_path = Path(workspace)
    owns_workspace = False
    render_error: BaseException | None = None
    try:
        workspace_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            workspace_path.mkdir()
        except FileExistsError as exc:
            raise ValidationError(
                "The Web render workspace must not already exist."
            ) from exc
        owns_workspace = True
        resource_directory = workspace_path / "resources"
        output_directory = workspace_path / "outputs"
        resource_directory.mkdir()
        output_directory.mkdir()
        resource_paths: dict[str, Path] = {}
        materialized_names: set[str] = set()
        for resource_id, entry in resources.items():
            if (
                not isinstance(resource_id, str)
                or not _RESOURCE_ID_RE.fullmatch(resource_id)
            ):
                raise ValidationError(
                    f"Invalid Web render resource ID: {resource_id!r}."
                )
            if not isinstance(entry, Mapping):
                raise ValidationError(
                    f"Web render resource {resource_id!r} must be an object."
                )
            name = safe_embedded_filename(entry.get("name"), fallback="")
            if not name:
                raise ValidationError(
                    f"Web render resource {resource_id!r} requires a safe filename."
                )
            if name in materialized_names:
                raise ValidationError(
                    f"Duplicate Web render resource filename: {name!r}."
                )
            materialized_names.add(name)
            resource_paths[resource_id] = materialize_embedded_file(
                entry,
                temp_dir=resource_directory,
                role=resource_id,
                prefix_role=False,
            )
        return render_canonical_web_request(
            payload,
            resource_paths=resource_paths,
            output_directory=output_directory,
        )
    except BaseException as exc:
        render_error = exc
        raise
    finally:
        if owns_workspace:
            try:
                shutil.rmtree(workspace_path)
            except FileNotFoundError:
                pass
            except OSError as cleanup_error:
                if render_error is not None:
                    _attach_exception_note(
                        render_error,
                        "Temporary Web render workspace cleanup also failed: "
                        f"{cleanup_error}",
                    )
                else:
                    raise ValidationError(
                        "The temporary Web render workspace could not be cleaned up."
                    ) from cleanup_error


def render_staged_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    workspace: str | Path,
    _diagnostics: MutableMapping[str, Any] | None = None,
) -> WebRenderReply:
    """Render a browser request from paths staged by the diagram Worker."""

    return _render_staged_canonical_web_request_with_prepared_inputs(
        payload,
        resource_paths=resource_paths,
        workspace=workspace,
        _diagnostics=_diagnostics,
        _prepared_input_cache=None,
        _resource_identities=None,
    )


@web_error_stage("resource-staging")
def _render_staged_canonical_web_request_with_prepared_inputs(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    workspace: str | Path,
    _diagnostics: MutableMapping[str, Any] | None = None,
    _prepared_input_cache: PreparedBiologicalInputCache | None = None,
    _resource_identities: Mapping[str, Mapping[str, Any]] | None = None,
    _label_fit_searches: LabelFitRowSearches | None = None,
    _held_feature_catalog_revision: str | None = None,
    _held_segments: AbstractSet[str] = frozenset(),
) -> WebRenderReply:
    if not isinstance(resource_paths, Mapping):
        raise ValidationError("The staged Web render resource map must be an object.")
    workspace_path = Path(workspace)
    resource_directory = workspace_path / "resources"
    cache_directory = workspace_path.parent / "gbdraw-web-render-resource-cache"
    output_directory = workspace_path / "outputs"
    owns_workspace = False
    render_error: BaseException | None = None
    try:
        marker_path = workspace_path / _STAGED_WORKSPACE_MARKER
        if (
            not _STAGED_WORKSPACE_RE.fullmatch(workspace_path.name)
            or not workspace_path.is_dir()
            or not resource_directory.is_dir()
            or not marker_path.is_file()
            or marker_path.is_symlink()
        ):
            raise ValidationError("The staged Web render workspace is incomplete.")
        owns_workspace = True
        workspace_identity = workspace_path.resolve(strict=True)
        resource_identity = resource_directory.resolve(strict=True)
        allowed_resource_parents = {resource_identity}
        if cache_directory.is_dir():
            allowed_resource_parents.add(cache_directory.resolve(strict=True))
        validated_paths: dict[str, Path] = {}
        for resource_id, raw_path in resource_paths.items():
            if (
                not isinstance(resource_id, str)
                or not _RESOURCE_ID_RE.fullmatch(resource_id)
            ):
                raise ValidationError(
                    f"Invalid staged Web render resource ID: {resource_id!r}."
                )
            path = Path(raw_path)
            if path.parent.resolve(strict=True) != resource_identity:
                raise ValidationError(
                    f"Staged Web render resource {resource_id!r} is outside its workspace."
                )
            resolved = path.resolve(strict=True)
            if resolved.parent not in allowed_resource_parents or not resolved.is_file():
                raise ValidationError(
                    f"Staged Web render resource {resource_id!r} is outside its workspace."
                )
            validated_paths[resource_id] = path
        validated_identities = _validate_prepared_resource_identities(
            _resource_identities,
            resource_paths=validated_paths,
            cache_requested=_prepared_input_cache is not None,
        )
        if output_directory.exists():
            raise ValidationError("The staged Web render output directory already exists.")
        output_directory.mkdir()
        if output_directory.resolve(strict=True).parent != workspace_identity:
            raise ValidationError("The staged Web render output directory is invalid.")
        return _render_canonical_web_request_with_prepared_inputs(
            payload,
            resource_paths=validated_paths,
            output_directory=output_directory,
            _diagnostics=_diagnostics,
            _prepared_input_cache=_prepared_input_cache,
            _resource_identities=validated_identities,
            _label_fit_searches=_label_fit_searches,
            _held_feature_catalog_revision=_held_feature_catalog_revision,
            _held_segments=_held_segments,
        )
    except BaseException as exc:
        render_error = exc
        raise
    finally:
        if owns_workspace:
            try:
                shutil.rmtree(workspace_path)
            except FileNotFoundError:
                pass
            except OSError as cleanup_error:
                if render_error is not None:
                    _attach_exception_note(
                        render_error,
                        "Temporary staged Web render workspace cleanup also failed: "
                        f"{cleanup_error}",
                    )
                else:
                    raise ValidationError(
                        "The temporary staged Web render workspace could not be cleaned up."
                    ) from cleanup_error


def _validate_prepared_resource_identities(
    value: Mapping[str, Mapping[str, Any]] | None,
    *,
    resource_paths: Mapping[str, Path],
    cache_requested: bool,
) -> dict[str, PreparedResourceIdentity] | None:
    if value is None:
        if cache_requested:
            raise ValidationError(
                "The staged Web render prepared-resource map is missing."
            )
        return None
    if not cache_requested:
        raise ValidationError(
            "Prepared resource identities require a prepared input cache."
        )
    if not isinstance(value, Mapping) or set(value) != set(resource_paths):
        raise ValidationError(
            "The staged Web render prepared-resource map does not match its manifest."
        )
    identities: dict[str, PreparedResourceIdentity] = {}
    for resource_id, raw_identity in value.items():
        if (
            not isinstance(resource_id, str)
            or not _RESOURCE_ID_RE.fullmatch(resource_id)
            or not isinstance(raw_identity, Mapping)
            or set(raw_identity) != {"cacheToken", "size"}
        ):
            raise ValidationError(
                f"Invalid prepared resource identity for {resource_id!r}."
            )
        cache_token = raw_identity.get("cacheToken")
        size = raw_identity.get("size")
        if (
            not isinstance(cache_token, str)
            or not _RESOURCE_CACHE_TOKEN_RE.fullmatch(cache_token)
            or not isinstance(size, int)
            or isinstance(size, bool)
            or size < 0
        ):
            raise ValidationError(
                f"Invalid prepared resource identity for {resource_id!r}."
            )
        try:
            actual_size = resource_paths[resource_id].resolve(strict=True).stat().st_size
        except OSError as exc:
            raise ValidationError(
                f"Prepared resource {resource_id!r} could not be inspected."
            ) from exc
        if actual_size != size:
            raise ValidationError(
                f"Prepared resource {resource_id!r} does not match its byte size."
            )
        identities[resource_id] = PreparedResourceIdentity(
            resource_id=resource_id,
            cache_token=cache_token,
            size=size,
        )
    return identities


__all__ = [
    "render_canonical_web_request",
    "render_embedded_canonical_web_request",
    "render_staged_canonical_web_request",
]
