"""Public canonical session document lifecycle.

Canonical sessions carry a CLI-independent ``renderRequest`` and a mapping of
embedded resources. Resource paths exposed by :class:`MaterializedSession` are
temporary and are valid only while the materialization context is active.
"""

from __future__ import annotations

import copy
import logging
import tempfile
from contextlib import AbstractContextManager, contextmanager
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterator, Literal, Mapping, NamedTuple, Sequence, cast

from gbdraw.exceptions import ValidationError
from gbdraw.session_contract import SESSION_ID_RE, GbdrawSession49, SessionResult
from gbdraw.io.table_text import legacy_table_rows
from gbdraw.render.output_paths import OutputExistsError, preflight_output_paths
from gbdraw.session_io import (
    CURRENT_AUTHORITY_SESSION_MIN_VERSION,
    CURRENT_SESSION_VERSION,
    CANONICAL_SESSION_MIN_VERSION,
    DEPTH_FILE_ENCODING,
    DRAWINGS_SESSION_MIN_VERSION,
    MODE_SCOPED_SESSION_MAX_VERSION,
    MODE_SCOPED_SESSION_MIN_VERSION,
    ROW_KEYED_LEGEND_SESSION_MIN_VERSION,
    SESSION_FORMAT,
    _read_session_json,
    _attach_current_web_file_bindings,
    _embedded_resource_bytes,
    CLI_SESSION_ADMISSION,
    WEB_READ_SESSION_ADMISSION,
    SessionAdmission,
    _write_validated_session_json,
    expand_session_feature_catalog,
    materialize_embedded_file,
    normalize_current_session_artifacts,
    safe_embedded_filename,
    seal_protein_identity_manifest,
    session_legend_edits_by_row_key,
    thaw_protein_identity_manifest,
    validate_session,
)
from gbdraw.session_drawings import (
    DRAWING_NAMES,
    DrawingMode,
    SessionDrawingArtifacts,
    SessionDrawingParts,
    active_session_drawing,
    allocate_session_drawing,
    drawings_document_from_mode_scoped,
    other_mode_fields,
    read_session_drawings,
    replace_session_drawings,
    write_session_drawings,
)

if TYPE_CHECKING:
    from gbdraw.api.migration_report import SessionMigrationReport
    from gbdraw.api.request_render import (
        CircularBatchRenderResult,
        RequestRenderResult,
    )
    from gbdraw.api.requests import DiagramRequest
    from gbdraw.session_migration import DraftMigrationRepair


logger = logging.getLogger(__name__)

_RESOURCE_REQUIRED_FIELDS = frozenset(
    {"kind", "name", "type", "size", "encoding", "data"}
)
_RESOURCE_OPTIONAL_FIELDS = frozenset({"lastModified", "checksum"})
class SessionError(ValidationError):
    """Base error for the public session bridge."""


class SessionFormatError(SessionError):
    """Raised when a session document or canonical envelope is malformed."""


class SessionVersionError(SessionError):
    """Raised when a session version has no requested conversion path."""


class SessionResourceError(SessionError):
    """Raised when an embedded canonical resource cannot be materialized."""


class SessionConversionError(SessionError):
    """Raised when a canonical payload cannot be converted to a typed request."""


class SessionRenderError(SessionError):
    """Raised when canonical session rendering fails."""


def require_search_frame_protein_results(report: SessionMigrationReport) -> None:
    """Refuse to write a Session whose replay kept protein results in the drawn frame.

    A replay draws an item the Session 47 migrator could not place as stored;
    a current Session must hold every protein result in the search frame
    (OV-1145).
    """

    if report.unconverted_protein_results:
        raise SessionConversionError(
            f"Cannot write a version {CURRENT_SESSION_VERSION} Session: "
            f"{report.unconverted_protein_results[0]}"
            + (
                f" (and {len(report.unconverted_protein_results) - 1} more)"
                if len(report.unconverted_protein_results) > 1
                else ""
            )
            + ". Generate the drawing again to compute its protein results."
        )


class SessionDrawingSelectionError(SessionError):
    """Raised when a drawing selector is missing, unknown, or ambiguous.

    The message lists the Session's drawings: ID, mode and name.
    """


@dataclass(frozen=True)
class SessionDrawing:
    """One drawing of a Session (project): an ID, a name and a mode.

    ``has_canonical_request`` is false for a drawing without a committed
    render, such as a settings-only Session. This is not an SVG drawing; see
    ``RequestRenderResult.drawing`` for that.
    """

    id: str
    name: str
    mode: DrawingMode
    has_canonical_request: bool


@dataclass(frozen=True)
class SessionDrawingSpec:
    """One drawing for :func:`build_session_document`.

    ``mode`` is required without a ``request`` and must match it otherwise.
    ``id`` and ``name`` default to the layout's names (``"circular"``,
    ``"Circular"``). ``state`` holds Web-owned drawing fields such as
    ``config``, ``ui``, ``features``, ``results`` or ``editorState``.
    """

    request: DiagramRequest | None = None
    mode: DrawingMode | None = None
    id: str | None = None
    name: str | None = None
    state: Mapping[str, Any] | None = None


def _describe_drawings(drawings: Sequence[SessionDrawingParts]) -> str:
    if not drawings:
        return "no drawing"
    return ", ".join(
        f"{drawing.id} ({drawing.mode}, {drawing.name!r}"
        + ("" if drawing.request is not None else ", no committed render")
        + ")"
        for drawing in drawings
    )


def _select_drawing(
    drawings: Sequence[SessionDrawingParts],
    selector: str | None,
    *,
    mode: str | None = None,
) -> SessionDrawingParts:
    """One drawing by exact ID, else by a unique exact name, else the only one."""

    if selector is not None:
        matches = [drawing for drawing in drawings if drawing.id == selector] or [
            drawing for drawing in drawings if drawing.name == selector
        ]
        if len(matches) != 1:
            raise SessionDrawingSelectionError(
                (
                    f"Session has no drawing {selector!r}"
                    if not matches
                    else f"Session has several drawings named {selector!r}"
                )
                + f"; it has: {_describe_drawings(drawings)}."
            )
        if mode is not None and matches[0].mode != mode:
            raise SessionDrawingSelectionError(
                f"Drawing {matches[0].id!r} is a {matches[0].mode} drawing, not {mode}; "
                f"the Session has: {_describe_drawings(drawings)}."
            )
        return matches[0]
    candidates = [drawing for drawing in drawings if mode is None or drawing.mode == mode]
    if len(candidates) == 1:
        return candidates[0]
    kind = f"{mode} drawing" if mode is not None else "drawing"
    if not candidates:
        raise SessionDrawingSelectionError(
            f"Session has no {kind}; it has: {_describe_drawings(drawings)}."
        )
    raise SessionDrawingSelectionError(
        f"Session has {len(candidates)} {kind}s; select one: {_describe_drawings(candidates)}."
    )


@dataclass(frozen=True)
class SessionDocument:
    """Validated session envelope detached from caller-owned mutable mappings.

    ``render_session_drawings_cli`` without a re-save drops the resources'
    embedded data (``_release_materialized_resource_data``). After that the
    document, and any drawing artifacts that share its resources, are no
    longer valid for saving or reuse.
    """

    _data: Mapping[str, Any]
    source_path: Path | None = None
    _drawings: tuple[SessionDrawingParts, ...] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        self._adopt(copy.deepcopy(dict(self._data)), CLI_SESSION_ADMISSION)

    @classmethod
    def _from_parsed(cls, data: dict[str, Any], source_path: Path) -> SessionDocument:
        """Keep a payload just parsed from ``source_path``: no caller holds it."""

        document = object.__new__(cls)
        object.__setattr__(document, "source_path", source_path)
        document._adopt(data, CLI_SESSION_ADMISSION)
        return document

    @classmethod
    def _admitted(cls, data: Mapping[str, Any], admission: SessionAdmission) -> SessionDocument:
        """A copy of ``data`` checked with the checks of ``admission``."""

        document = object.__new__(cls)
        object.__setattr__(document, "source_path", None)
        document._adopt(copy.deepcopy(dict(data)), admission)
        return document

    def _adopt(self, data: dict[str, Any], admission: SessionAdmission) -> None:
        expanded = expand_session_feature_catalog(data)
        # Checked once here; every later reader of the document trusts the proof.
        seal_protein_identity_manifest(expanded)
        _validate_document(expanded, admission)
        object.__setattr__(self, "_data", expanded)
        object.__setattr__(self, "_drawings", read_session_drawings(expanded))
        if self.source_path is not None:
            object.__setattr__(self, "source_path", Path(self.source_path))

    @classmethod
    def _from_checked_move(cls, data: Mapping[str, Any]) -> SessionDocument:
        """Keep the Session 49 that a checked Session 46-48 document moved to, unchecked.

        Only the CLI re-save keeps one (``_upgraded_session``): the document it
        writes holds every drawing of this one and is checked (OV-1815).
        """

        document = object.__new__(cls)
        object.__setattr__(document, "source_path", None)
        moved = dict(data)
        seal_protein_identity_manifest(moved)
        object.__setattr__(document, "_data", moved)
        object.__setattr__(document, "_drawings", read_session_drawings(moved))
        return document

    @property
    def version(self) -> int:
        """Session envelope version."""

        return int(self._data["version"])

    @property
    def drawings(self) -> tuple[SessionDrawing, ...]:
        """The Session's drawings in document order."""

        return tuple(_public_drawing(drawing) for drawing in self._drawings)

    def drawing(
        self,
        selector: str | None = None,
        *,
        mode: DrawingMode | None = None,
    ) -> SessionDrawing:
        """One drawing by ID or unique name; without one, the only drawing.

        ``mode`` restricts the choice to drawings of that mode. An unknown,
        ambiguous or missing selection raises
        :class:`SessionDrawingSelectionError`, which lists the drawings.
        """

        return _public_drawing(self._drawing_parts(selector, mode=mode))

    @property
    def active_drawing_id(self) -> str:
        """The ID of the drawing the Web app opens."""

        active = active_session_drawing(self._data, self._drawings)
        if active is None:
            raise SessionDrawingSelectionError("Session has no drawing.")
        return active.id

    @property
    def mode(self) -> DrawingMode | None:
        """The mode of the only drawing, or ``None`` without a drawing.

        A Session with several drawings raises
        :class:`SessionDrawingSelectionError`; use :meth:`drawing`.
        """

        if not self._drawings:
            return None
        return self._drawing_parts().mode

    @property
    def has_canonical_request(self) -> bool:
        """Whether the only drawing can enter the canonical typed bridge.

        A Session with several drawings raises
        :class:`SessionDrawingSelectionError`; use :meth:`drawing`.
        """

        if not self._drawings:
            return False
        return self._drawing_parts().request is not None

    def to_dict(self) -> dict[str, Any]:
        """Return a detached JSON-compatible copy of the document."""

        data = copy.deepcopy(dict(self._data))
        thaw_protein_identity_manifest(data)
        return data

    def _drawing_parts(
        self,
        selector: str | None = None,
        *,
        mode: str | None = None,
    ) -> SessionDrawingParts:
        return _select_drawing(self._drawings, selector, mode=mode)


def _public_drawing(drawing: SessionDrawingParts) -> SessionDrawing:
    return SessionDrawing(
        id=drawing.id,
        name=drawing.name,
        mode=drawing.mode,
        has_canonical_request=drawing.request is not None,
    )


def session_drawing_artifacts(
    document: SessionDocument,
    drawing: str | None = None,
) -> SessionDrawingArtifacts:
    """The validated view of one drawing that the compatibility adapter reads."""

    return document._drawing_parts(drawing).artifacts


@dataclass
class _MaterializationLifetime:
    active: bool = True


@dataclass(frozen=True)
class MaterializedSession:
    """Materialized canonical resources with context-bounded path lifetime.

    ``resource_paths`` and any paths stored in a decoded request become invalid
    as soon as the owning ``materialize_session(...)`` context exits.
    """

    document: SessionDocument
    temp_directory: Path
    resource_paths: Mapping[str, Path]
    output_directory: Path
    _lifetime: _MaterializationLifetime

    @property
    def active(self) -> bool:
        """Whether temporary resources are still owned by an active context."""

        return self._lifetime.active


class _SessionMaterializationContext(AbstractContextManager[MaterializedSession]):
    def __init__(
        self,
        document: SessionDocument,
        *,
        output_directory: str | Path,
        temporary_directory: str | Path | None,
    ) -> None:
        raw_output_directory = str(output_directory).strip()
        if not raw_output_directory:
            raise SessionResourceError("A replay output directory is required.")
        self._document = document
        self._output_directory = Path(raw_output_directory)
        self._temporary_directory = (
            Path(temporary_directory) if temporary_directory is not None else None
        )
        self._owner: tempfile.TemporaryDirectory[str] | None = None
        self._materialized: MaterializedSession | None = None
        self._table_rows: AbstractContextManager[None] | None = None

    def __enter__(self) -> MaterializedSession:
        if self._materialized is not None:
            raise SessionResourceError("A session materialization context is single-use.")
        try:
            if self._temporary_directory is not None:
                self._temporary_directory.mkdir(parents=True, exist_ok=True)
            self._owner = tempfile.TemporaryDirectory(
                prefix="gbdraw-session-v31-",
                dir=(
                    str(self._temporary_directory)
                    if self._temporary_directory is not None
                    else None
                ),
            )
            temp_directory = Path(self._owner.name)
            resource_paths = _materialize_resources(
                self._document,
                temp_directory=temp_directory,
            )
            self._materialized = MaterializedSession(
                document=self._document,
                temp_directory=temp_directory,
                resource_paths=resource_paths,
                output_directory=self._output_directory,
                _lifetime=_MaterializationLifetime(),
            )
            # Sessions 31-39 read their Default colors, Label whitelist and
            # Qualifier priority rows as the current writer writes them, as
            # Web Load does (OV-148); the saved resource is not rewritten.
            version = self._document.version
            self._table_rows = legacy_table_rows(
                CANONICAL_SESSION_MIN_VERSION <= version < CURRENT_AUTHORITY_SESSION_MIN_VERSION
            )
            self._table_rows.__enter__()
            return self._materialized
        except SessionError:
            self._cleanup_after_failed_enter()
            raise
        except (OSError, ValidationError) as exc:
            self._cleanup_after_failed_enter()
            raise SessionResourceError(
                f"Canonical session resources could not be materialized: {exc}"
            ) from exc

    def __exit__(self, exc_type, exc_value, traceback) -> Literal[False]:
        if self._materialized is not None:
            self._materialized._lifetime.active = False
        if self._table_rows is not None:
            self._table_rows.__exit__(None, None, None)
        if self._owner is None:
            return False
        try:
            self._owner.cleanup()
        except OSError as cleanup_error:
            if exc_value is not None:
                if hasattr(exc_value, "add_note"):
                    exc_value.add_note(
                        f"Temporary session cleanup also failed: {cleanup_error}"
                    )
                return False
            raise SessionResourceError(
                "Temporary canonical session resources could not be cleaned up."
            ) from cleanup_error
        return False

    def _cleanup_after_failed_enter(self) -> None:
        if self._owner is None:
            return
        try:
            self._owner.cleanup()
        except OSError as cleanup_error:
            raise SessionResourceError(
                "Session materialization failed and temporary resources could not be cleaned up."
            ) from cleanup_error


def load_session_document(
    source: str | Path | Mapping[str, Any] | SessionDocument,
) -> SessionDocument:
    """Load and validate a session document without materializing resources."""

    if isinstance(source, SessionDocument):
        return source
    if isinstance(source, Mapping):
        try:
            return SessionDocument(source)
        except SessionError:
            raise
        except ValidationError as exc:
            raise _classify_validation_error(exc) from exc

    path = Path(source)
    try:
        payload = _read_session_json(path)
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    if not isinstance(payload, dict):
        raise SessionFormatError("Session JSON must be an object.")
    try:
        return SessionDocument._from_parsed(payload, path)
    except SessionError:
        raise
    except ValidationError as exc:
        raise _classify_validation_error(exc) from exc


def materialize_session(
    document: SessionDocument | Mapping[str, Any] | str | Path,
    *,
    output_directory: str | Path,
    temporary_directory: str | Path | None = None,
) -> AbstractContextManager[MaterializedSession]:
    """Create the sole owner of temporary canonical session resources."""

    return _SessionMaterializationContext(
        load_session_document(document),
        output_directory=output_directory,
        temporary_directory=temporary_directory,
    )


def _release_materialized_resource_data(materialized: MaterializedSession) -> None:
    """Drop each resource's embedded data once its file is materialized.

    Only for a document its caller owns and neither reads nor saves again:
    the CLI render without a re-save. The resource descriptors keep every
    other field, the materialized files hold the bytes, and the Session file
    is untouched. A re-save and the Python API keep the data.
    """

    _active_materialization(materialized)
    resources = materialized.document._data.get("resources", {})
    for resource_id in materialized.resource_paths:
        resources[resource_id].pop("data", None)


def _active_materialization(materialized: MaterializedSession) -> None:
    if not isinstance(materialized, MaterializedSession):
        raise SessionConversionError("A MaterializedSession is required.")
    if not materialized.active:
        raise SessionResourceError(
            "Materialized session resources are no longer active; decode inside the context."
        )


def _with_saved_rule_captions(request: DiagramRequest) -> DiagramRequest:
    """The request of a saved Session with the captions of its color table read as saved.

    A color table file of the request is read so by ``read_color_table``.
    """

    from gbdraw.legend.table import saved_rule_captions

    colors = request.options.colors
    if colors is None or colors.color_table is None:
        return request
    color_table = saved_rule_captions(colors.color_table)
    if color_table is colors.color_table:
        return request
    options = replace(request.options, colors=replace(colors, color_table=color_table))
    # Each request type takes its own options type, which mypy cannot pair.
    return cast("DiagramRequest", replace(cast(Any, request), options=options))


def _decode_drawing(
    materialized: MaterializedSession,
    drawing: SessionDrawingParts,
) -> DiagramRequest:
    """Decode one drawing's committed request within its resource lifetime."""

    version = drawing.artifacts.version
    if drawing.request is None:
        if version >= CANONICAL_SESSION_MIN_VERSION:
            raise SessionConversionError("Settings-only Session has no biological render request; load a source in Web before generating.")
        raise SessionVersionError(
            "Sessions version 27 through 30 support internal CLI replay only and do not "
            "have a public typed-request conversion."
        )
    from gbdraw.session_request_codec import CanonicalRequestCodecError
    from gbdraw.api.session_compat import (
        decode_session_request,
        promote_legacy_session_similarity_alignment_request,
    )

    try:
        request = decode_session_request(
            version,
            drawing.request,
            resource_paths=materialized.resource_paths,
            output_directory=materialized.output_directory,
        )
        return promote_legacy_session_similarity_alignment_request(
            _with_saved_rule_captions(request),
            drawing.artifacts,
        )
    except (CanonicalRequestCodecError, ValidationError) as exc:
        raise SessionConversionError(str(exc)) from exc


def session_to_request(
    materialized: MaterializedSession,
    *,
    drawing: str | None = None,
) -> DiagramRequest:
    """Decode a canonical request within its resource lifetime.

    ``drawing`` selects a drawing by ID or name; without it the Session must
    have one drawing.
    """

    _active_materialization(materialized)
    return _decode_drawing(
        materialized,
        materialized.document._drawing_parts(drawing),
    )


@dataclass(frozen=True)
class _SessionDrawingRender:
    """One selected drawing and its decoded request with output overrides."""

    drawing: SessionDrawingParts
    request: DiagramRequest


def _drawing_output_request(
    request: DiagramRequest,
    drawing_id: str,
    *,
    several: bool,
    output_prefix: str | None,
    formats: str | tuple[str, ...] | None,
    overwrite: bool | None,
    embed_fonts: bool | None = None,
) -> DiagramRequest:
    """Apply the output overrides; several drawings get ``<base>_<id>``."""

    from gbdraw.api.requests import CircularBatchRequest

    if several:
        if output_prefix is not None:
            output_prefix = f"{output_prefix}_{drawing_id}"
        elif isinstance(request, CircularBatchRequest):
            request = replace(
                request,
                outputs=tuple(
                    replace(output, output_prefix=f"{output.output_prefix}_{drawing_id}")
                    for output in request.outputs
                ),
            )
        else:
            output_prefix = f"{request.output.output_prefix}_{drawing_id}"
    return with_request_output(
        request,
        output_prefix=output_prefix,
        formats=formats,
        overwrite=overwrite,
        embed_fonts=embed_fonts,
    )


def _plan_session_drawings(
    materialized: MaterializedSession,
    selectors: Sequence[str] | None,
    *,
    output_prefix: str | None = None,
    formats: str | tuple[str, ...] | None = None,
    overwrite: bool | None = None,
    embed_fonts: bool | None = None,
) -> tuple[_SessionDrawingRender, ...]:
    """Select drawings, decode their requests and name their outputs.

    Without selectors every drawing with a committed render is selected and
    the others are skipped with a notice; a named drawing without one is an
    error. The plan keeps document order.
    """

    _active_materialization(materialized)
    document = materialized.document
    if selectors is None:
        selected = [drawing for drawing in document._drawings if drawing.request is not None]
        for drawing in document._drawings:
            if drawing.request is None:
                logger.warning(
                    "Drawing %r (%s) has no committed render and is skipped; "
                    "Generate it in the Web app first.",
                    drawing.id,
                    drawing.mode,
                )
        if not selected:
            raise SessionDrawingSelectionError(
                "Session has no drawing with a committed render; it has: "
                f"{_describe_drawings(document._drawings)}."
            )
    else:
        if isinstance(selectors, str):
            raise SessionDrawingSelectionError(
                "Select drawings with a sequence of IDs or names, not one string."
            )
        if not selectors:
            raise SessionDrawingSelectionError("Select at least one drawing.")
        chosen = [document._drawing_parts(selector) for selector in selectors]
        chosen_ids = [drawing.id for drawing in chosen]
        repeated = sorted({drawing_id for drawing_id in chosen_ids if chosen_ids.count(drawing_id) > 1})
        if repeated:
            raise SessionDrawingSelectionError(
                "Drawing(s) selected more than once: " + ", ".join(repeated) + "."
            )
        for drawing in chosen:
            if drawing.request is None:
                raise SessionDrawingSelectionError(
                    f"Drawing {drawing.id!r} has no committed render; Generate it in the Web app first."
                )
        selected = [drawing for drawing in document._drawings if drawing.id in chosen_ids]
    several = len(selected) > 1
    return tuple(
        _SessionDrawingRender(
            drawing=drawing,
            request=_drawing_output_request(
                _decode_drawing(materialized, drawing),
                drawing.id,
                several=several,
                output_prefix=output_prefix,
                formats=formats,
                overwrite=overwrite,
                embed_fonts=embed_fonts,
            ),
        )
        for drawing in selected
    )


@contextmanager
def _session_parse_cache(materialized: MaterializedSession) -> Iterator[None]:
    """Parse each materialized resource once for every drawing of one call.

    One ``PreparedBiologicalInputCache`` transaction spans the call, so its
    retention is the selected drawings' resources; the cache ends with the
    call. The callers decode the drawings inside it, as a Worker Generate
    decodes inside its transaction, so a decoded resource is keyed by its
    resource identity instead of by its value. An error of the body keeps its
    own type; opening or committing the cache raises
    :class:`SessionRenderError`.
    """

    from gbdraw.api.prepared import (
        PreparedBiologicalInputCache,
        PreparedResourceIdentity,
    )

    # Before the files are read: a closed materialization is a resource error.
    _active_materialization(materialized)
    body_failed = False
    try:
        identities: dict[str | Path, PreparedResourceIdentity] = {
            path: PreparedResourceIdentity(
                resource_id=resource_id,
                cache_token=f"session-{resource_id}",
                size=path.stat().st_size,
            )
            for resource_id, path in materialized.resource_paths.items()
        }
        with PreparedBiologicalInputCache().transaction(
            resource_paths=identities,
            diagnostics=None,
        ):
            try:
                yield
            except BaseException:
                body_failed = True
                raise
    except SessionError:
        raise
    except Exception as exc:
        if body_failed:
            raise
        raise SessionRenderError(f"Canonical session rendering failed: {exc}") from exc


def _render_session_drawing_plans(
    materialized: MaterializedSession,
    plans: Sequence[_SessionDrawingRender],
    *,
    include_feature_catalog: bool = False,
) -> dict[str, RequestRenderResult | CircularBatchRenderResult]:
    """Render planned drawings inside the call's :func:`_session_parse_cache`.

    Its callers (``render_session``, ``render_session_drawings`` and the CLI's
    ``render_session_drawings_cli``) open that cache around decoding and
    rendering; without it a decoded resource is keyed by its value and the
    drawings parse their records again.
    """

    from gbdraw.api.session_compat import render_session_compatible_request

    _active_materialization(materialized)
    results: dict[str, RequestRenderResult | CircularBatchRenderResult] = {}
    try:
        for plan in plans:
            results[plan.drawing.id] = render_session_compatible_request(
                plan.request,
                plan.drawing.artifacts,
                include_feature_catalog=include_feature_catalog,
            )
    except SessionError:
        raise
    except Exception as exc:
        raise SessionRenderError(f"Canonical session rendering failed: {exc}") from exc
    return results


def render_session(
    materialized: MaterializedSession,
    *,
    drawing: str | None = None,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Decode and render a canonical session while its resources are active.

    ``drawing`` selects a drawing by ID or name; without it the Session must
    have one drawing.
    """

    with _session_parse_cache(materialized):
        selected = materialized.document._drawing_parts(drawing)
        plan = _SessionDrawingRender(selected, _decode_drawing(materialized, selected))
        return _render_session_drawing_plans(materialized, (plan,))[selected.id]


def render_session_drawings(
    materialized: MaterializedSession,
    *,
    drawings: Sequence[str] | None = None,
    output_prefix: str | None = None,
    formats: str | tuple[str, ...] | None = None,
    overwrite: bool | None = None,
) -> dict[str, RequestRenderResult | CircularBatchRenderResult]:
    """Render several drawings of a materialized Session together.

    ``drawings`` selects drawings by ID or name. By default every drawing
    with a committed render is rendered and the others are skipped with a
    logged notice; naming a drawing without one raises
    :class:`SessionDrawingSelectionError`.

    One drawing keeps its output names. Several drawings write
    ``<base>_<id>``, where ``<base>`` is ``output_prefix`` or the drawing's
    own prefix; a Circular batch inside still appends ``_<n>``. Every output
    path of every selected drawing is checked before the first file is
    written, and each resource is parsed once for all drawings. The results
    are keyed by drawing ID in document order.
    """

    from gbdraw.api.request_render import preflight_diagram_request_outputs

    with _session_parse_cache(materialized):
        plans = _plan_session_drawings(
            materialized,
            drawings,
            output_prefix=output_prefix,
            formats=formats,
            overwrite=overwrite,
        )
        try:
            preflight_diagram_request_outputs(tuple(plan.request for plan in plans))
        except ValidationError as exc:
            raise SessionRenderError(f"Canonical session rendering failed: {exc}") from exc
        # A copy typed as the public union (``dict`` is invariant in its values).
        return dict(_render_session_drawing_plans(materialized, plans))


@dataclass(frozen=True)
class _DrawingBuild:
    """One drawing to write: its mode, resolved request and Web-owned fields."""

    mode: str
    request: DiagramRequest | None
    state: Mapping[str, Any]
    id: str | None = None
    name: str | None = None


_RESERVED_DRAWING_FIELDS = frozenset(
    {"format", "version", "createdAt", "renderRequest", "resources"}
)


def _build_session_document_from_drawings(
    drawings: Sequence[_DrawingBuild],
    *,
    base: SessionDocument | None = None,
    title: str | None = None,
    created_at: datetime | None = None,
    active_drawing: str | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
    resources: Mapping[str, Mapping[str, Any]] | None = None,
) -> SessionDocument:
    """Build a current (Session 49) document from already-resolved drawings.

    With a ``base`` (a Session 49) the drawings replace its drawings of their
    IDs and the other drawings stay; otherwise the drawings are the
    Session's, in the given order (D-15-125). ``resources`` are the
    resources of the Session being saved again: equal bytes keep their IDs,
    and those that no request or Web file binding names are dropped. A
    drawing's ``state`` holds its own mode's slice only: another mode's
    (the retired ``modes`` key of Sessions 46-48) is refused.
    """

    for drawing in drawings:
        other = other_mode_fields(drawing.state, drawing.mode)
        if other:
            raise SessionFormatError(
                f"A {drawing.mode} drawing's state cannot hold "
                + ", ".join(other)
                + "; give that mode a drawing of its own."
            )
    encoded = _encoded_drawings(drawings, allocate=base is None, resources=resources)
    active = _active_drawing_index(encoded, active_drawing)
    if base is not None:
        if base.version < DRAWINGS_SESSION_MIN_VERSION:
            raise SessionFormatError(
                f"A Session {base.version} moves to drawings before it is saved again."
            )
        return _drawings_session_document(
            base,
            dict(zip(encoded.ids, encoded.views, strict=True)),
            resources=encoded.resources,
            title=title,
            created_at=created_at,
            web_file_inventory=web_file_inventory,
        )
    timestamp = created_at or datetime.now(timezone.utc)
    views = [
        _enveloped(view, DRAWINGS_SESSION_MIN_VERSION, timestamp, encoded.resources)
        for view in encoded.views
    ]
    _finish_views(views, title=title, web_file_inventory=web_file_inventory)
    if resources:
        from gbdraw.session_resources import drawing_resource_ids

        referenced = set().union(*(drawing_resource_ids(view) for view in views))
        for view in views:
            _keep_resources(view, referenced | (set(encoded.resources) - set(resources)))
    try:
        fields = write_session_drawings(views, ids=encoded.ids, names=encoded.names, active=active)
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    return SessionDocument(fields)


def _mode_scoped_document(
    drawing: _DrawingBuild,
    *,
    created_at: datetime | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
    resources: Mapping[str, Mapping[str, Any]] | None = None,
) -> SessionDocument:
    """The checked mode-scoped (Session 48) document of one drawing upgraded
    from a Session 31-45, before its move to drawings."""

    return SessionDocument(
        _mode_scoped_data(
            drawing, created_at=created_at, web_file_inventory=web_file_inventory, resources=resources,
        )
    )


def _mode_scoped_document_of_request(
    request: DiagramRequest,
    *,
    resolved: bool,
    created_at: datetime | None,
) -> SessionDocument:
    """:func:`_mode_scoped_document` of one request, resolved here unless ``resolved``
    (:func:`gbdraw.session_io.build_session_json` builds on it)."""

    (build,) = (
        (_DrawingBuild(mode=_diagram_request_mode(request), request=request, state={}),)
        if resolved
        else _drawing_builds(request, None)
    )
    return _mode_scoped_document(build, created_at=created_at)


def _mode_scoped_data(
    drawing: _DrawingBuild,
    *,
    created_at: datetime | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
    resources: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """The mode-scoped (Session 48) data of one drawing of an older Session,
    not yet checked; the move to drawings
    (:func:`gbdraw.session_drawings.drawings_document_from_mode_scoped`)
    reads its ``modes`` and ``ui.mode``. ``resources`` are the source's."""

    encoded = _encoded_drawings((drawing,), allocate=True, resources=resources)
    (view,) = encoded.views
    data = _enveloped(
        view, MODE_SCOPED_SESSION_MAX_VERSION, created_at or datetime.now(timezone.utc), encoded.resources,
    )
    if web_file_inventory is not None:
        _attach_current_web_file_bindings(data, web_file_inventory)
    if resources:
        _drop_unreferenced_resources(data, resources)
    normalize_current_session_artifacts(data)
    return data


@dataclass(frozen=True)
class _EncodedDrawings:
    """Drawings with their requests encoded through one resource table."""

    views: list[dict[str, Any]]
    ids: list[str]
    names: list[str]
    resources: dict[str, Any]


def _encoded_drawings(
    drawings: Sequence[_DrawingBuild],
    *,
    allocate: bool,
    resources: Mapping[str, Mapping[str, Any]] | None,
) -> _EncodedDrawings:
    """The one-drawing view of each drawing, its request encoded through one
    resource table that starts from ``resources``. New drawings (``allocate``)
    get their IDs and names from their modes; drawings saved again keep theirs."""

    from gbdraw.session_request_codec import (
        CanonicalRequestCodecError,
        encode_canonical_request,
    )
    from gbdraw.session_resources import SessionResourceTable
    from gbdraw.api.session_compat import (
        project_legacy_similarity_alignment_for_current_write,
    )

    table = SessionResourceTable(resources)
    views: list[dict[str, Any]] = []
    taken_ids: list[str] = []
    names: list[str] = []
    for drawing in drawings:
        if allocate:
            try:
                drawing_id, name = allocate_session_drawing(
                    drawing.mode,
                    requested_id=drawing.id,
                    requested_name=drawing.name,
                    taken_ids=taken_ids,
                )
            except ValidationError as exc:
                raise SessionFormatError(str(exc)) from exc
        else:
            drawing_id, name = drawing.id or drawing.mode, drawing.name or ""
        taken_ids.append(drawing_id)
        names.append(name)
        # The views are assembled from new mappings; SessionDocument(data)
        # detaches the document from ``drawing.state`` with one deep copy.
        state = dict(drawing.state)
        conflicting = _RESERVED_DRAWING_FIELDS & set(state)
        if conflicting:
            raise SessionFormatError(
                "Session drawing state cannot replace canonical field(s): "
                + ", ".join(sorted(conflicting))
                + "."
            )
        try:
            payload = (
                encode_canonical_request(
                    project_legacy_similarity_alignment_for_current_write(drawing.request),
                    table=table,
                ).payload
                if drawing.request is not None
                else None
            )
        except (CanonicalRequestCodecError, ValidationError) as exc:
            raise SessionConversionError(str(exc)) from exc
        view: dict[str, Any] = {
            "renderRequest": payload,
            "results": [],
            "editorState": {"featureCatalog": None},
        }
        view.update(state)
        editor_state = view.get("editorState")
        if isinstance(editor_state, Mapping):
            normalized_editor_state = dict(editor_state)
            normalized_editor_state.setdefault("featureCatalog", None)
            view["editorState"] = normalized_editor_state
        if payload is None:
            view["ui"] = {**dict(view.get("ui") or {}), "mode": drawing.mode}
        views.append(view)
    try:
        descriptors = table.descriptors()
    except (CanonicalRequestCodecError, ValidationError) as exc:
        raise SessionConversionError(str(exc)) from exc
    return _EncodedDrawings(views, taken_ids, names, descriptors)


def _active_drawing_index(encoded: _EncodedDrawings, active_drawing: str | None) -> int | None:
    """The index of the drawing whose ID or name is ``active_drawing``."""

    if active_drawing is None:
        return None
    for index, (drawing_id, name) in enumerate(zip(encoded.ids, encoded.names, strict=True)):
        if active_drawing in (drawing_id, name):
            return index
    raise SessionDrawingSelectionError(
        f"active_drawing {active_drawing!r} names none of the drawings: "
        + ", ".join(encoded.ids)
        + "."
    )


def _enveloped(
    view: Mapping[str, Any],
    version: int,
    created_at: datetime,
    resources: Mapping[str, Any],
) -> dict[str, Any]:
    """``view`` behind the Session envelope, in the written key order."""

    return {
        "format": SESSION_FORMAT,
        "version": version,
        "createdAt": created_at.isoformat(),
        "renderRequest": view["renderRequest"],
        "resources": resources,
        **{key: value for key, value in view.items() if key != "renderRequest"},
    }


def _finish_views(
    views: Sequence[dict[str, Any]],
    *,
    title: str | None,
    web_file_inventory: Mapping[str, Any] | None,
) -> None:
    """Bind each view's Web files, set the Session title and normalize its artifacts."""

    for view in views:
        if web_file_inventory is not None:
            _attach_current_web_file_bindings(view, web_file_inventory)
        if title is not None:
            view["title"] = str(title)
        normalize_current_session_artifacts(view)


def _build_session_document_from_resolved_request(
    request: DiagramRequest,
    *,
    title: str | None = None,
    created_at: datetime | None = None,
    adjunct: Mapping[str, Any] | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
    resources: Mapping[str, Mapping[str, Any]] | None = None,
) -> SessionDocument:
    """Build a one-drawing current document from an already-resolved request."""

    return _build_session_document_from_drawings(
        (_DrawingBuild(mode=_diagram_request_mode(request), request=request, state=adjunct or {}),),
        title=title,
        created_at=created_at,
        web_file_inventory=web_file_inventory,
        resources=resources,
    )


def _drawings_session_document(
    base: SessionDocument,
    views: dict[str, dict[str, Any]],
    *,
    resources: Mapping[str, Mapping[str, Any]],
    title: str | None,
    created_at: datetime | None,
    web_file_inventory: Mapping[str, Any] | None,
) -> SessionDocument:
    """A Session 49 saved again: each view replaces the base drawing of its ID.

    The resources that no drawing names any more are dropped before the
    writer, which refuses them.
    """

    from gbdraw.session_resources import drawing_resource_ids

    shared = dict(resources)
    for view in views.values():
        view["resources"] = shared
    _finish_views(list(views.values()), title=title, web_file_inventory=web_file_inventory)
    kept = [drawing.artifacts.fields for drawing in base._drawings if drawing.id not in views]
    referenced = set().union(*(drawing_resource_ids(fields) for fields in (*views.values(), *kept)))
    for view in views.values():
        _keep_resources(view, referenced)
    try:
        fields = replace_session_drawings(base._data, views)
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    return SessionDocument({
        "format": SESSION_FORMAT,
        "createdAt": (created_at or datetime.now(timezone.utc)).isoformat(),
        **fields,
    })


def _drop_unreferenced_resources(
    data: dict[str, Any],
    previous: Mapping[str, Any],
) -> None:
    """Drop the previous resources that the requests and Web files no longer name."""

    from gbdraw.session_resources import drawing_resource_ids

    referenced = set().union(
        *(drawing_resource_ids(drawing.artifacts.fields) for drawing in read_session_drawings(data))
    )
    _keep_resources(data, referenced | (set(data["resources"]) - set(previous)))


def _keep_resources(data: dict[str, Any], kept: set[str]) -> None:
    """Keep only the resources ``kept`` names in ``data``, and their original names."""

    web_files = data.get("webFiles")
    data["resources"] = {
        resource_id: descriptor
        for resource_id, descriptor in data["resources"].items()
        if resource_id in kept
    }
    original_names = web_files.get("resourceOriginalNames") if isinstance(web_files, Mapping) else None
    if isinstance(web_files, Mapping) and isinstance(original_names, Mapping):
        data["webFiles"] = {
            **web_files,
            "resourceOriginalNames": {
                resource_id: name
                for resource_id, name in original_names.items()
                if resource_id in data["resources"]
            },
        }


def _diagram_request_mode(request: DiagramRequest) -> str:
    from gbdraw.api.requests import LinearDiagramRequest

    return "linear" if isinstance(request, LinearDiagramRequest) else "circular"


def _drawing_builds(
    request: DiagramRequest | None,
    drawings: Sequence[DiagramRequest | SessionDrawingSpec] | None,
) -> tuple[_DrawingBuild, ...]:
    """Resolve the requests of ``request`` or ``drawings`` (exactly one)."""

    from gbdraw.api.requests import CircularBatchRequest, CircularDiagramRequest, LinearDiagramRequest
    from gbdraw.api.session_compat import resolve_session_write_request

    if (request is None) == (drawings is None):
        raise SessionFormatError("Pass either request or drawings.")
    items: Sequence[DiagramRequest | SessionDrawingSpec] = (
        (request,) if request is not None else tuple(drawings or ())
    )
    builds: list[_DrawingBuild] = []
    for item in items:
        spec = (
            item
            if isinstance(item, SessionDrawingSpec)
            else SessionDrawingSpec(request=item)
        )
        if spec.request is not None and not isinstance(
            spec.request,
            (CircularDiagramRequest, CircularBatchRequest, LinearDiagramRequest),
        ):
            raise SessionFormatError(
                "A Session drawing takes a typed diagram request or a SessionDrawingSpec."
            )
        mode = _diagram_request_mode(spec.request) if spec.request is not None else spec.mode
        if mode is None:
            raise SessionFormatError("A Session drawing without a request needs a mode.")
        if spec.mode is not None and spec.mode != mode:
            raise SessionFormatError(
                f"Session drawing mode {spec.mode!r} does not match its {mode} request."
            )
        try:
            resolved = resolve_session_write_request(spec.request) if spec.request is not None else None
        except ValidationError as exc:
            raise SessionConversionError(str(exc)) from exc
        builds.append(
            _DrawingBuild(
                mode=mode,
                request=resolved,
                state=spec.state or {},
                id=spec.id,
                name=spec.name,
            )
        )
    return tuple(builds)


def build_session_document(
    request: DiagramRequest | None = None,
    *,
    drawings: Sequence[DiagramRequest | SessionDrawingSpec] | None = None,
    title: str | None = None,
    created_at: datetime | None = None,
    active_drawing: str | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
) -> SessionDocument:
    """Build a current-version document from one request or several drawings.

    Pass ``request`` for a one-drawing Session, or ``drawings`` (typed
    requests or :class:`SessionDrawingSpec` values) in drawing order.
    ``active_drawing`` names the drawing the Web app opens. A drawing's
    ``state`` may contain Web/editor fields such as ``ui`` or ``results``; it
    cannot replace canonical envelope fields, and it holds its own mode's
    draft only: another mode's slice (``state["modes"]["linear"]`` on a
    Circular drawing) raises :class:`SessionFormatError`. The current
    Session version is 49: it keeps the drawings in the given order, any
    number of them, several of one mode, and drawings without a request.
    Each drawing gets its ID and name from its mode (``circular``,
    ``circular-2``; ``Circular``, ``Circular 2``).
    """

    return _build_session_document_from_drawings(
        _drawing_builds(request, drawings),
        title=title,
        created_at=created_at,
        active_drawing=active_drawing,
        web_file_inventory=web_file_inventory,
    )


def _write_session_document(
    path: str | Path,
    document: SessionDocument,
    *,
    overwrite: bool,
) -> SessionDocument:
    """Write one built session document through a staged commit."""

    try:
        _preflight_session_output(path, overwrite=overwrite)
        # The document was validated when it was built.
        _write_validated_session_json(path, document._data, overwrite=overwrite)
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    return document


def _preflight_session_output(path: str | Path, *, overwrite: bool) -> None:
    """Refuse a Session output path before anything is built or written."""

    output_path = Path(path)
    preflight_output_paths((output_path,), overwrite=True)
    if (output_path.exists() or output_path.is_symlink()) and not overwrite:
        raise OutputExistsError((output_path,), subject="Session output")


def save_session_document(
    path: str | Path,
    request: DiagramRequest | None = None,
    *,
    drawings: Sequence[DiagramRequest | SessionDrawingSpec] | None = None,
    title: str | None = None,
    created_at: datetime | None = None,
    active_drawing: str | None = None,
    web_file_inventory: Mapping[str, Any] | None = None,
    overwrite: bool = False,
) -> SessionDocument:
    """Build and write a current canonical session through a staged commit."""

    try:
        # A refused output is reported before the drawings are rendered.
        _preflight_session_output(path, overwrite=overwrite)
    except ValidationError as exc:
        raise SessionFormatError(str(exc)) from exc
    document = build_session_document(
        request,
        drawings=drawings,
        title=title,
        created_at=created_at,
        active_drawing=active_drawing,
        web_file_inventory=web_file_inventory,
    )
    return _write_session_document(path, document, overwrite=overwrite)


@dataclass(frozen=True)
class SessionUpgrade:
    """A Session in the current version, from :func:`upgrade_session_document`.

    ``warnings`` has one line per drawing whose Results the upgrade dropped,
    naming each dropped Result; the same lines are logged as warnings.
    Rendering the drawing and saving the Session writes new Results.
    """

    document: SessionDocument
    warnings: tuple[str, ...] = ()


def _dropped_results_warning(
    source_version: int,
    drawing_id: str,
    results: Any,
) -> str | None:
    if not isinstance(results, list) or not results:
        return None
    names = ", ".join(
        repr(result.get("name"))
        if isinstance(result, Mapping) and isinstance(result.get("name"), str)
        else f"#{index + 1}"
        for index, result in enumerate(results)
    )
    return (
        f"Upgrading Session {source_version} to {CURRENT_SESSION_VERSION} dropped "
        f"the {drawing_id} drawing's {len(results)} Result(s) {names}: Session "
        f"{source_version} saved no feature catalog for them. Render the drawing "
        "and save the Session to write new Results."
    )


# The stored Results of a Session 27-39 that a Session 49 cannot hold, by
# the ID of the Session 49 drawing they belong to (OV-1659).
StoredResults = dict[str, list[SessionResult]]


# Where a Session 27-44 holds its Result set: the top level, its one set (a
# Session before 46 has no ``otherModeResult``).
LegacyResultSetPath = Literal[""]
# The file Result set each Session 49 drawing of a Session 27-44 took, by the
# drawing's ID: the drawing that holds the file's request or its set-aside
# Results (D-15-130). A drawing that took none is not named.
LegacyResultSetSources = dict[str, LegacyResultSetPath]


class LegacyDrawings(NamedTuple):
    """:func:`drawings_document_from_legacy` of a Session 27-44."""

    document: GbdrawSession49
    stored_results: StoredResults
    result_set_sources: LegacyResultSetSources
    repairs: tuple[DraftMigrationRepair, ...]


def _set_aside_results(stored_results: StoredResults, drawing_id: str, state: dict[str, Any]) -> None:
    """Move the Results a Session 49 cannot hold out of ``state`` (OV-1659).

    Those of a Session 27-30 have no request and those of a Session 31-39 no
    feature catalog; a current Result needs both. Web Load shows them from
    ``stored_results`` and Save refuses them until one Generate.
    """

    results = state.get("results")
    if isinstance(results, list) and results:
        stored_results[drawing_id] = cast(list[SessionResult], results)
    state["results"] = []
    state.pop("runMetadata", None)


def _upgrade_mode_scoped_session(
    loaded: SessionDocument,
    temporary_directory: str | Path | None,
) -> SessionDocument:
    """A Session 46 or 47 as the current version, every drawing kept.

    Session 46 gets the Legend edits of each slice keyed by row key (47).
    The stored protein results of each Linear drawing go to each record's
    search frame (48, OV-1145): only their resources' bytes change, under
    their stored IDs.
    """

    from gbdraw.api.session_compat import (
        describe_unmatched_frame_item,
        drawn_frame_protein_results_to_search_frame,
        replace_resource_bytes,
        search_frame_protein_resources,
        stores_protein_results,
    )

    data = dict(loaded._data)
    if loaded.version < ROW_KEYED_LEGEND_SESSION_MIN_VERSION:
        data, _outcomes = session_legend_edits_by_row_key(data)
    drawings = [
        drawing for drawing in loaded._drawings
        if drawing.mode == "linear" and stores_protein_results(drawing.request)
    ]
    if drawings:
        from gbdraw.api.request_render import plan_request

        replaced: dict[str, bytes] = {}
        with materialize_session(
            loaded, output_directory=Path("."), temporary_directory=temporary_directory,
        ) as materialized:
            for drawing in drawings:
                assert drawing.request is not None  # stores_protein_results
                request = _decode_drawing(materialized, drawing)
                try:
                    converted, unmatched = drawn_frame_protein_results_to_search_frame(
                        plan_request(request), request,
                    )
                    if unmatched:
                        raise SessionConversionError(
                            f"Cannot upgrade drawing {drawing.id!r} to Session {CURRENT_SESSION_VERSION}: "
                            f"{describe_unmatched_frame_item(unmatched[0])}."
                        )
                    for resource_id, content in search_frame_protein_resources(
                        drawing.request, request, converted,
                    ).items():
                        if replaced.get(resource_id, content) != content:
                            raise SessionConversionError(
                                f"Session resource {resource_id!r} holds two different protein results."
                            )
                        replaced[resource_id] = content
                except SessionError:
                    raise
                except ValidationError as exc:
                    raise SessionConversionError(str(exc)) from exc
        resources = loaded._data["resources"]
        data["resources"] = {
            **resources,
            **{key: replace_resource_bytes(resources[key], content) for key, content in replaced.items()},
        }
    return SessionDocument({**data, "version": MODE_SCOPED_SESSION_MAX_VERSION})


def upgrade_session_document(
    document: SessionDocument | Mapping[str, Any] | str | Path,
    *,
    temporary_directory: str | Path | None = None,
) -> SessionUpgrade:
    """Return a Session in the current version (49) without rendering it.

    A Session 27-47 is upgraded to the mode-scoped Session 48 as below, and a
    Session 46-48 then moves to drawings
    (:func:`gbdraw.session_drawings.drawings_document_from_mode_scoped`).

    A current document is returned unchanged. A Session 46 or 47 keeps its
    drawings; a Session 46 gets its Legend color and stroke edits keyed by
    row key (:func:`gbdraw.session_io.legend_edits_by_row_key`), and the
    stored protein results of a reverse-complemented Linear record are
    written in its search frame under their resource IDs (OV-1145); a result
    matching neither strand of its record raises
    :class:`SessionConversionError`. A Session 31-44 gets the
    migrations a CLI re-save applies: its request is decoded, adapted to
    current typed state (comparison frames, similarity alignment, LOSAT
    artifacts) and encoded again with the same resource IDs, and its
    Web-owned fields are migrated; the rendered-ID feature edits of a Session
    31-39 are named through its GenBank sources read again with their crops
    and orientations, as Web Load names them, or else through its saved
    feature metadata. Results with a feature catalog (Sessions
    40-44) are kept. Sessions 31-39 saved no catalog, so their Results are
    dropped and the drawing waits for its next render; each drawing with
    dropped Results gets a warning in :attr:`SessionUpgrade.warnings` that
    names them, and the warning is logged. The resources are materialized
    under ``temporary_directory`` while the request is adapted. Sessions
    27-30 have no canonical request and raise :class:`SessionVersionError`.
    """

    return _upgraded_session(load_session_document(document), temporary_directory, checked=True)


def _upgraded_session(
    document: SessionDocument,
    temporary_directory: str | Path | None,
    *,
    checked: bool,
) -> SessionUpgrade:
    """:func:`upgrade_session_document`. ``checked=False`` (the CLI re-save,
    which checks the document it writes) keeps the move to drawings unchecked."""

    upgrade = _upgrade_to_current(document, temporary_directory)
    if upgrade.document.version >= DRAWINGS_SESSION_MIN_VERSION:
        return upgrade
    # The move reads the document's own data, whose manifest stays sealed
    # (``to_dict`` thaws it), and builds a new document that SessionDocument
    # copies: the manifest is not checked again.
    drawings = cast(Mapping[str, Any], _drawings_of_current(upgrade.document._data))
    moved = SessionDocument(drawings) if checked else SessionDocument._from_checked_move(drawings)
    return SessionUpgrade(moved, upgrade.warnings)


def _drawings_of_current(
    data: Mapping[str, Any],
    *,
    stored_results: StoredResults | None = None,
) -> GbdrawSession49:
    """A mode-scoped (46-48) document's data moved to drawings, the drawings of
    ``stored_results`` counting as drawings with a Result; a refusal is a
    :class:`SessionConversionError`."""

    try:
        drawings, _source_facts = drawings_document_from_mode_scoped(
            data, stored_result_modes=frozenset(stored_results or ())
        )
    except SessionError:
        raise
    except ValidationError as exc:
        raise SessionConversionError(str(exc)) from exc
    return drawings


def drawings_document_from_legacy(
    data: SessionDocument | Mapping[str, Any],
    *,
    temporary_directory: str | Path | None = None,
) -> LegacyDrawings:
    """A Session 27-44 as a Session 49, as Web Load reads it, with its stored
    Results, where each drawing's Result set was in the file, and Load repairs.

    Web Load's readers read the file first and send what they kept; it is
    admitted with :data:`gbdraw.session_io.WEB_READ_SESSION_ADMISSION`, here
    and in the check of the Session 49 (D-15-131).

    The migrations are those of :func:`upgrade_session_document` (the mode
    split, feature edits, annotation targets, hash rules, placements, legacy
    Linear comparison drafts, the sources of a Session 31-39 read again, and
    the request adaptation with its comparison and protein frames), then the
    move of :func:`gbdraw.session_drawings.drawings_document_from_mode_scoped`.
    Unlike the CLI re-save, no stored Result is
    dropped: the Results of a Session 40-44 stay in the document with their
    catalog, and those a Session 49 cannot hold (27-30 have no request,
    31-39 no feature catalog) are returned beside it, by the ID of their
    drawing, which keeps ``results: []`` (OV-1659). A Session 27-30
    is a drawing without a request; its document goes to drawings without
    the Session 48 check, which refuses a settings-only file with biological
    sources while Session 49 admits an undrawn drawing with its bindings.
    The returned document is checked as a Session 49. The repairs are what
    the draft migrations could not keep, as Web Load reports them. Other
    versions raise :class:`SessionVersionError`.
    """

    loaded = data if isinstance(data, SessionDocument) else _web_read_document(data)
    if not loaded.version < MODE_SCOPED_SESSION_MIN_VERSION:
        raise SessionVersionError(
            f"Only a Session 27-{MODE_SCOPED_SESSION_MIN_VERSION - 2} moves to drawings "
            f"through this step; this one is {loaded.version}.",
            diagnostic={"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"},
        )
    repairs: list[DraftMigrationRepair] = []
    stored_results: StoredResults = {}
    current = (
        _request_less_data(loaded, repairs=repairs, stored_results=stored_results)
        if loaded._drawing_parts().request is None
        else _upgrade_to_current(
            loaded, temporary_directory, stored_results=stored_results, repairs=repairs,
        ).document.to_dict()
    )
    drawings = _drawings_of_current(current, stored_results=stored_results)
    drawing_ids = {drawing["id"] for drawing in drawings["drawings"]}
    if not stored_results.keys() <= drawing_ids:
        raise SessionConversionError(
            "The drawings " + ", ".join(sorted(stored_results.keys() - drawing_ids))
            + " of the stored Results are not in the Session 49."
        )
    _web_read_document(cast(Mapping[str, Any], drawings))
    # The file's one set went to the drawing of its request, or was set aside
    # for the drawing of its Results.
    set_sources: LegacyResultSetSources = {
        drawing["id"]: ""
        for drawing in drawings["drawings"]
        if drawing.get("renderRequest") is not None or drawing["id"] in stored_results
    }
    if len(set_sources) > 1:
        raise SessionConversionError(
            "The Result set of a Session 27-44 went to more than one drawing: "
            + ", ".join(sorted(set_sources)) + "."
        )
    return LegacyDrawings(drawings, stored_results, set_sources, tuple(repairs))


def _web_read_document(data: Mapping[str, Any]) -> SessionDocument:
    """``data`` as a document checked as Web Load admits it after its own
    readers (:data:`WEB_READ_SESSION_ADMISSION`)."""

    try:
        return SessionDocument._admitted(data, WEB_READ_SESSION_ADMISSION)
    except SessionError:
        raise
    except ValidationError as exc:
        raise _classify_validation_error(exc) from exc


def _project_adjunct(
    loaded: SessionDocument,
    *,
    repairs: list[DraftMigrationRepair] | None,
    source_features: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """The current drawing state and Web file inventory of a Session 27-44."""

    from gbdraw.session_migration import project_session_adjunct_for_current_write

    try:
        return project_session_adjunct_for_current_write(
            loaded._data,
            source_version=loaded.version,
            source_features=source_features,
            repairs=repairs,
        )
    except SessionError:
        raise
    except ValidationError as exc:
        raise SessionConversionError(str(exc)) from exc


def _request_less_data(
    loaded: SessionDocument,
    *,
    repairs: list[DraftMigrationRepair] | None = None,
    stored_results: StoredResults | None = None,
) -> dict[str, Any]:
    """The current data of a Session 27-44 whose drawing has no request, not yet checked.

    A drawing without a request has no Results' catalog (Session 48 and 49);
    its Results go to ``stored_results`` when given (Web Load), by the ID of
    the drawing, else stay for the check to refuse. The document keeps the
    mode shown at Save (``ui.mode``), which is a use of that mode,
    so the drawing of a set-aside Result may be another mode.
    """

    from gbdraw.session_migration import replace_current_derived_feature_state

    drawing = loaded._drawing_parts()
    state, web_file_inventory = _project_adjunct(loaded, repairs=repairs)
    if stored_results is not None:
        _set_aside_results(stored_results, drawing.id, state)
    replace_current_derived_feature_state(state, None)
    ui = state.get("ui")
    shown = ui.get("mode") if isinstance(ui, Mapping) else None
    mode = shown if isinstance(shown, str) and shown in DRAWING_NAMES else drawing.mode
    return _mode_scoped_data(
        _DrawingBuild(mode=mode, request=None, state=state),
        web_file_inventory=web_file_inventory,
        resources=loaded._data.get("resources"),
    )


def _upgrade_to_current(
    loaded: SessionDocument,
    temporary_directory: str | Path | None,
    *,
    stored_results: StoredResults | None = None,
    repairs: list[DraftMigrationRepair] | None = None,
) -> SessionUpgrade:
    """``loaded`` in the newest mode-scoped version (48), or a Session 48 or
    49 unchanged. With ``stored_results`` (Web Load's reading), Results
    without a feature catalog go there instead of being dropped (OV-1659)."""

    if loaded.version >= MODE_SCOPED_SESSION_MAX_VERSION:
        return SessionUpgrade(loaded)
    if loaded.version >= MODE_SCOPED_SESSION_MIN_VERSION:
        return SessionUpgrade(_upgrade_mode_scoped_session(loaded, temporary_directory))
    if loaded.version < CANONICAL_SESSION_MIN_VERSION:
        raise SessionVersionError(
            "Sessions version 27 through 30 have no canonical request; replay one with "
            "gbdraw circular or gbdraw linear --session and --session_output to write a "
            "current Session."
        )
    from gbdraw.session_migration import (
        read_legacy_source_features,
        replace_current_derived_feature_state,
        with_current_artifacts,
    )

    source = loaded._data
    # A Session before the current version holds one drawing.
    drawing = loaded._drawing_parts()
    if drawing.request is None:
        return SessionUpgrade(SessionDocument(_request_less_data(loaded, repairs=repairs)))
    from gbdraw.api.session_compat import (
        adapt_session_request,
        project_legacy_similarity_alignment_for_current_write,
    )
    from gbdraw.web_support.feature_catalog import promote_legacy_feature_catalog

    with materialize_session(
        loaded,
        output_directory=Path("."),
        temporary_directory=temporary_directory,
    ) as materialized:
        # Rendered-ID edits without a saved catalog are named through the
        # sources read again, as Web Load reads them.
        state, web_file_inventory = _project_adjunct(
            loaded,
            repairs=repairs,
            source_features=read_legacy_source_features(source, materialized.resource_paths),
        )
        request = _decode_drawing(materialized, drawing)
        try:
            adapted = adapt_session_request(request, drawing.artifacts)
        except SessionError:
            raise
        except ValidationError as exc:
            raise SessionConversionError(str(exc)) from exc
        require_search_frame_protein_results(adapted.migration_report)
        state = with_current_artifacts(
            state,
            losat_cache_entries=adapted.artifacts.losat_cache_entries,
            protein_identity_manifest=adapted.artifacts.protein_identity_manifest,
            legacy_protein_raw_candidates=adapted.migration_report.protein_raw_candidates,
            legacy_protein_derived_evidence=adapted.migration_report.protein_derived_evidence,
            protein_id_map=adapted.migration_report.protein_id_map,
        )
        editor_state = state.get("editorState")
        catalog = (
            editor_state.get("featureCatalog")
            if isinstance(editor_state, Mapping)
            else None
        )
        if isinstance(catalog, Mapping) and catalog.get("schema") in (3, 4):
            catalog = promote_legacy_feature_catalog(catalog)
        warnings: list[str] = []
        if not isinstance(catalog, Mapping):
            catalog = None
        if catalog is None and stored_results is not None:
            _set_aside_results(stored_results, drawing.id, state)
        if catalog is None:
            dropped = _dropped_results_warning(
                loaded.version, drawing.id, state.get("results")
            )
            if dropped is not None:
                warnings.append(dropped)
            state["results"] = []
            state.pop("runMetadata", None)
        replace_current_derived_feature_state(state, catalog)
        upgraded = _mode_scoped_document(
            _DrawingBuild(
                mode=drawing.mode,
                request=project_legacy_similarity_alignment_for_current_write(
                    adapted.request,
                    legacy_source=request,
                ),
                state=state,
            ),
            web_file_inventory=web_file_inventory,
            resources=source["resources"],
        )
    for warning in warnings:
        logger.warning("WARNING: %s", warning)
    return SessionUpgrade(upgraded, tuple(warnings))


def with_request_output(
    request: DiagramRequest,
    *,
    output_prefix: str | None = None,
    output_directory: str | Path | None = None,
    formats: str | tuple[str, ...] | None = None,
    overwrite: bool | None = None,
    embed_fonts: bool | None = None,
) -> DiagramRequest:
    """Return a request with caller-owned replay output overrides."""

    from gbdraw.api.requests import CircularBatchRequest, RenderOutputRequest

    def updated(current: RenderOutputRequest, prefix: str) -> RenderOutputRequest:
        return replace(
            current,
            output_prefix=prefix,
            resolve_prefix_from_first_record=False,
            output_directory=(
                output_directory
                if output_directory is not None
                else current.output_directory
            ),
            formats=formats if formats is not None else current.formats,
            overwrite=current.overwrite if overwrite is None else overwrite,
            embed_fonts=current.embed_fonts if embed_fonts is None else embed_fonts,
        )

    if isinstance(request, CircularBatchRequest):
        item_count = len(request.outputs)
        if output_prefix is None:
            prefixes = tuple(output.output_prefix for output in request.outputs)
        elif item_count == 1:
            prefixes = (output_prefix,)
        else:
            prefixes = tuple(
                f"{output_prefix}_{index}"
                for index in range(1, item_count + 1)
            )
        outputs = tuple(
            updated(current, prefix)
            for prefix, current in zip(
                prefixes,
                request.outputs,
                strict=True,
            )
        )
        return replace(request, outputs=outputs)

    return replace(
        request,
        output=updated(request.output, output_prefix or request.output.output_prefix),
    )


def _validate_document(data: Mapping[str, Any], admission: SessionAdmission) -> None:
    try:
        validate_session(data, admission)
    except ValidationError as exc:
        raise _classify_validation_error(exc) from exc
    if int(data.get("version", 0)) < CANONICAL_SESSION_MIN_VERSION:
        return
    resources = data.get("resources")
    assert isinstance(resources, Mapping)
    sanitized_names: set[str] = set()
    for resource_id, entry in resources.items():
        if not isinstance(resource_id, str) or not SESSION_ID_RE.fullmatch(resource_id):
            raise SessionResourceError(
                f"Invalid canonical resource ID: {resource_id!r}."
            )
        if not isinstance(entry, Mapping):
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} must be an object."
            )
        fields = set(entry)
        missing = _RESOURCE_REQUIRED_FIELDS - fields
        unknown = fields - _RESOURCE_REQUIRED_FIELDS - _RESOURCE_OPTIONAL_FIELDS
        if missing:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} is missing field(s): "
                + ", ".join(sorted(missing))
                + "."
            )
        if unknown:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} has unknown field(s): "
                + ", ".join(sorted(unknown))
                + "."
            )
        if not isinstance(entry.get("kind"), str) or not str(entry["kind"]).strip():
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} requires a kind."
            )
        name = safe_embedded_filename(entry.get("name"), fallback="")
        if not name:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} requires a safe basename."
            )
        if name in sanitized_names:
            raise SessionResourceError(
                f"Duplicate canonical resource filename after sanitization: {name!r}."
            )
        sanitized_names.add(name)
        encoding = entry.get("encoding")
        if encoding not in {"base64", DEPTH_FILE_ENCODING}:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} has unsupported encoding {encoding!r}."
            )
        if encoding == "base64" and not isinstance(entry.get("data"), str):
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} requires base64 string data."
            )
        if encoding == DEPTH_FILE_ENCODING and not isinstance(entry.get("data"), Mapping):
            raise SessionResourceError(
                f"Canonical depth resource {resource_id!r} requires an object payload."
            )
        size = entry.get("size")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} has invalid size metadata."
            )
        if "checksum" in entry:
            try:
                _embedded_resource_bytes(entry)
            except ValidationError as exc:
                raise SessionResourceError(
                    f"Canonical resource {resource_id!r} is invalid: {exc}"
                ) from exc
    for drawing in read_session_drawings(data):
        unresolved = drawing.resource_ids - set(resources)
        if unresolved:
            raise SessionResourceError(
                f"The {drawing.id} drawing's renderRequest references missing canonical resource(s): "
                + ", ".join(sorted(unresolved))
                + "."
            )


def _materialize_resources(
    document: SessionDocument,
    *,
    temp_directory: Path,
) -> dict[str, Path]:
    if document.version < CANONICAL_SESSION_MIN_VERSION:
        return {}
    raw_resources = document._data.get("resources")
    assert isinstance(raw_resources, Mapping)
    result: dict[str, Path] = {}
    for resource_id, entry in raw_resources.items():
        assert isinstance(resource_id, str)
        assert isinstance(entry, Mapping)
        try:
            result[resource_id] = materialize_embedded_file(
                entry,
                temp_dir=temp_directory,
                role=resource_id,
                prefix_role=False,
            )
        except ValidationError as exc:
            raise SessionResourceError(
                f"Canonical resource {resource_id!r} could not be materialized: {exc}"
            ) from exc
    return result


def _classify_validation_error(exc: ValidationError) -> SessionError:
    message = str(exc)
    if "version" in message.lower():
        return SessionVersionError(message)
    return SessionFormatError(message)


__all__ = [
    "MaterializedSession",
    "SessionConversionError",
    "SessionDocument",
    "SessionDrawing",
    "SessionDrawingSelectionError",
    "SessionDrawingSpec",
    "SessionError",
    "SessionFormatError",
    "SessionRenderError",
    "SessionResourceError",
    "SessionUpgrade",
    "SessionVersionError",
    "build_session_document",
    "load_session_document",
    "materialize_session",
    "render_session",
    "render_session_drawings",
    "save_session_document",
    "session_to_request",
    "upgrade_session_document",
    "with_request_output",
]
