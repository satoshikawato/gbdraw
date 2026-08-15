"""Canonical typed-request bridge used by the browser render worker."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import re
import shutil
from time import perf_counter
from typing import Any, Iterator, Mapping, MutableMapping

from gbdraw.api.request_render import (
    CircularBatchRenderResult,
    RequestRenderResult,
    _capture_request_render_diagnostics,
    render_request,
)
from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
)
from gbdraw.exceptions import ValidationError
from gbdraw.render.formats import SVG_FORMAT, resolve_format_output_path
from gbdraw.render.track_slot_metadata import (
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
)
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.session_io import materialize_embedded_file, safe_embedded_filename
from gbdraw.web_support.feature_catalog import (
    build_feature_catalog,
    build_feature_catalog_item,
)


_RESOURCE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")
_RESOURCE_CACHE_TOKEN_RE = re.compile(r"^render-resource-[1-9][0-9]*$")
_STAGED_WORKSPACE_RE = re.compile(r"^gbdraw-web-render-[1-9][0-9]*$")
_STAGED_WORKSPACE_MARKER = ".gbdraw-worker-render-workspace"


@contextmanager
def _web_render_diagnostic_phase(
    diagnostics: MutableMapping[str, Any] | None,
    name: str,
) -> Iterator[None]:
    if diagnostics is None:
        yield
        return
    started_at = perf_counter()
    try:
        yield
    finally:
        timings = diagnostics.setdefault("timingsMs", {})
        timings[name] = float(timings.get(name, 0.0)) + (
            perf_counter() - started_at
        ) * 1000.0


def _attach_exception_note(error: BaseException, note: str) -> None:
    """Attach a diagnostic on Python versions before BaseException.add_note."""

    add_note = getattr(error, "add_note", None)
    if callable(add_note):
        add_note(note)
        return
    notes = getattr(error, "__notes__", None)
    if isinstance(notes, list):
        notes.append(note)
    else:
        error.__notes__ = [note]  # type: ignore[attr-defined]


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


def render_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
    _diagnostics: MutableMapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Decode and render one browser request without a CLI translation layer."""

    return _render_canonical_web_request_with_prepared_inputs(
        payload,
        resource_paths=resource_paths,
        output_directory=output_directory,
        _diagnostics=_diagnostics,
        _prepared_input_cache=None,
        _resource_identities=None,
    )


def _render_canonical_web_request_with_prepared_inputs(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
    _diagnostics: MutableMapping[str, Any] | None,
    _prepared_input_cache: PreparedBiologicalInputCache | None,
    _resource_identities: Mapping[str, PreparedResourceIdentity] | None,
) -> dict[str, Any]:
    def render() -> dict[str, Any]:
        with _capture_request_render_diagnostics(_diagnostics):
            return _render_canonical_web_request(
                payload,
                resource_paths=resource_paths,
                output_directory=output_directory,
                diagnostics=_diagnostics,
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
) -> dict[str, Any]:

    if not isinstance(payload, Mapping):
        raise ValidationError("The Web render request must be an object.")
    if not isinstance(resource_paths, Mapping):
        raise ValidationError("The Web render resource map must be an object.")
    output_root = Path(output_directory)
    output_root.mkdir(parents=True, exist_ok=True)
    with _web_render_diagnostic_phase(diagnostics, "decode"):
        request = decode_canonical_request(
            payload,
            resource_paths=resource_paths,
            output_directory=output_root,
        )
    with _web_render_diagnostic_phase(diagnostics, "renderRequest"):
        rendered = render_request(request, include_feature_catalog=True)
    items: tuple[RequestRenderResult, ...]
    if isinstance(rendered, CircularBatchRenderResult):
        items = rendered.items
    else:
        items = (rendered,)

    results: list[dict[str, str]] = []
    geometry_records: list[dict[str, Any]] = []
    feature_catalog_items: list[dict[str, Any]] = []
    for result_index, item in enumerate(items):
        path = _base_svg_path(item)
        result_name = path.name
        with _web_render_diagnostic_phase(diagnostics, "svgReadback"):
            svg_content = path.read_text(encoding="utf-8")
        results.append(
            {
                "name": result_name,
                "content": svg_content,
            }
        )
        with _web_render_diagnostic_phase(diagnostics, "featureCatalog"):
            feature_catalog_items.append(
                build_feature_catalog_item(
                    svg_content,
                    item.interactive_context,
                    result_index=result_index,
                    result_name=result_name,
                    _diagnostics=diagnostics,
                )
            )
        with _web_render_diagnostic_phase(diagnostics, "geometryMetadata"):
            geometry_records.extend(
                collect_track_slot_geometry_records(
                    item.drawing,
                    result_index=result_index,
                    result_name=result_name,
                )
            )
    if not results:
        raise ValidationError("The canonical Web request did not produce an SVG.")
    with _web_render_diagnostic_phase(diagnostics, "geometryMetadata"):
        metadata = build_track_slot_geometry_run_metadata(
            mode=rendered.mode,
            records=geometry_records,
        )
    with _web_render_diagnostic_phase(diagnostics, "featureCatalog"):
        metadata["featureCatalog"] = build_feature_catalog(feature_catalog_items)
    return {
        "results": results,
        "metadata": metadata,
    }


def render_embedded_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resources: Mapping[str, Mapping[str, Any]],
    workspace: str | Path,
) -> dict[str, Any]:
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
) -> dict[str, Any]:
    """Render a browser request from paths staged by the diagram Worker."""

    return _render_staged_canonical_web_request_with_prepared_inputs(
        payload,
        resource_paths=resource_paths,
        workspace=workspace,
        _diagnostics=_diagnostics,
        _prepared_input_cache=None,
        _resource_identities=None,
    )


def _render_staged_canonical_web_request_with_prepared_inputs(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    workspace: str | Path,
    _diagnostics: MutableMapping[str, Any] | None = None,
    _prepared_input_cache: PreparedBiologicalInputCache | None = None,
    _resource_identities: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
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
