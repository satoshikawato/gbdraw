"""Canonical typed-request bridge used by the browser render worker."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
from typing import Any, Mapping

from gbdraw.api.request_render import (
    CircularBatchRenderResult,
    RequestRenderResult,
    render_request,
)
from gbdraw.exceptions import ValidationError
from gbdraw.render.track_slot_metadata import (
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
)
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.session_io import materialize_embedded_file, safe_embedded_filename


_RESOURCE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")


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


def render_canonical_web_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> dict[str, Any]:
    """Decode and render one browser request without a CLI translation layer."""

    if not isinstance(payload, Mapping):
        raise ValidationError("The Web render request must be an object.")
    if not isinstance(resource_paths, Mapping):
        raise ValidationError("The Web render resource map must be an object.")
    output_root = Path(output_directory)
    output_root.mkdir(parents=True, exist_ok=True)
    request = decode_canonical_request(
        payload,
        resource_paths=resource_paths,
        output_directory=output_root,
    )
    rendered = render_request(request)
    items: tuple[RequestRenderResult, ...]
    if isinstance(rendered, CircularBatchRenderResult):
        items = rendered.items
    else:
        items = (rendered,)

    results: list[dict[str, str]] = []
    geometry_records: list[dict[str, Any]] = []
    for result_index, item in enumerate(items):
        result_name = item.request.output.output_prefix
        for output_path in item.output_paths:
            path = Path(output_path)
            if path.suffix.lower() != ".svg":
                continue
            result_name = path.name
            results.append(
                {
                    "name": path.name,
                    "content": path.read_text(encoding="utf-8"),
                }
            )
        geometry_records.extend(
            collect_track_slot_geometry_records(
                item.drawing,
                result_index=result_index,
                result_name=result_name,
            )
        )
    if not results:
        raise ValidationError("The canonical Web request did not produce an SVG.")
    return {
        "results": results,
        "metadata": build_track_slot_geometry_run_metadata(
            mode=rendered.mode,
            records=geometry_records,
        ),
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


__all__ = [
    "render_canonical_web_request",
    "render_embedded_canonical_web_request",
]
