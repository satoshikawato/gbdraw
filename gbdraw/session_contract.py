"""The Session document of the current version, declared once.

Each class names the keys of one JSON object of a Session as the current
writers write it: the Web's Save, the CLI sidecar and
:func:`gbdraw.session.build_session_document`. The keys are the wire names.
The validators of :mod:`gbdraw.session_io` read their field sets from these
classes (:func:`keys_of`), and ``tools/generate_web_contracts.py`` writes the
Web's typedefs from :func:`describe_session_document`.

A value type is the type a writer that builds the value writes.
``validate_session`` checks keys and the parts it reads, not every leaf, and
the CLI copies the Web-owned parts of a source Session without reading them,
so a reader takes a Session as unvalidated data. ``object`` marks a value
whose owner has no description yet; the comment names that owner.

``get_type_hints`` reads these class bodies at run time, so every name in
them is importable at run time.
"""

from __future__ import annotations

from typing import Any, Literal, NewType, TypedDict, get_type_hints

from .features.overrides import FeatureVisibilityMode, LabelVisibilityMode
from .web_support.mode_scoped_settings import (
    DIAGRAM_MODES,
    LOSAT_EXECUTION_FIELDS,
    SLICE_CONTAINERS,
    DiagramMode,
)

# A value another owner describes; describe_session_document() names its shape.
# The canonical render request of the current schema (session_request_codec).
CanonicalRequestJson = NewType("CanonicalRequestJson", dict[str, object])
# ``modes``: one slice per diagram mode, each holding MODE_SCOPED_SETTINGS rows only.
ModeSlicesJson = NewType("ModeSlicesJson", dict[str, object])
# ``ui.losatExecution``: LOSAT_EXECUTION_FIELDS only.
LosatExecutionJson = NewType("LosatExecutionJson", dict[str, object])


class SessionResult(TypedDict):
    """One Result: a named plain SVG."""

    name: str
    content: str


class SessionRunMetadata(TypedDict, total=False):
    """What a Generate reported beside its Results.

    The rows are checked by ``_validate_other_mode_run_metadata`` and by the
    Web's warning and notice validators.
    """

    trackSlotGeometry: object
    annotationWarnings: list[object]
    featureIdentityNotices: list[object]
    comparisonWarnings: list[object]


class _CliInvocationRequired(TypedDict):
    schema: int
    mode: str
    args: list[str]
    renderFormats: list[str]
    fileBindings: list[object]
    generatedBy: str


class SessionCliInvocation(_CliInvocationRequired, total=False):
    """The command line that reproduces the Results: the CLI's own
    (``build_session_json``), or the one the Web shows for its run."""

    # Written by the Web: whether every file of the command is in the Session.
    sessionExportable: bool


class _ResultSetRequired(TypedDict):
    results: list[SessionResult]


class SessionResultSet(_ResultSetRequired, total=False):
    """The fields that one committed Result set names the same way at the top
    level and in ``otherModeResult``; each adds its own ``renderRequest``."""

    runMetadata: SessionRunMetadata
    cliInvocation: SessionCliInvocation


class SessionOtherModeLegend(TypedDict, total=False):
    """The Legend of a Result set as Generate drew it."""

    originalOrder: list[str]
    originalColors: dict[str, str]


class SessionSvgStroke(TypedDict):
    color: str | None
    width: float | None


class _OtherModeEditorStateRequired(TypedDict):
    # The admitted feature catalog (gbdraw.web_support.feature_catalog), or
    # null without Results; stored in its compact form.
    featureCatalog: object


class SessionOtherModeEditorState(_OtherModeEditorStateRequired, total=False):
    """The part of ``editorState`` that belongs to one Result set."""

    # Checked against the request by _validate_alignment_reset_receipt.
    alignmentResetReceipt: dict[str, object] | None
    legend: SessionOtherModeLegend
    originalSvgStroke: SessionSvgStroke


class SessionOtherModeUi(TypedDict, total=False):
    """The part of ``ui`` that belongs to one Result set: the shown Result and
    the layout and palette it was generated with. The Web owns the values."""

    selectedResultIndex: object
    generatedLegendPosition: object
    generatedMultiRecordCanvas: object
    generatedCircularPlotTitlePosition: object
    appliedPaletteName: object
    appliedPaletteColors: object


class _OtherModeResultRequired(SessionResultSet):
    renderRequest: CanonicalRequestJson
    editorState: SessionOtherModeEditorState


class SessionOtherModeResult(_OtherModeResultRequired, total=False):
    """The other diagram mode's committed Result set (E1), written only when
    both modes keep a Result. Its request's resources are in the top-level
    ``resources`` table."""

    ui: SessionOtherModeUi


class SessionEditorState(SessionOtherModeEditorState):
    """``editorState``: the shown Result set's part. A Session may hold other
    keys here; the Legend and stroke edits are each mode's (``modes``)."""


class SessionUi(SessionOtherModeUi, total=False):
    """``ui``: the shown Result set's part, the shown mode, and the Web's
    app-level view state, which only the Web reads."""

    mode: DiagramMode
    losatExecution: LosatExecutionJson


class SessionOrthogroupState(TypedDict, total=False):
    """``orthogroupState``: the Web's orthogroup selection; the group names
    are the Linear drawing's (``modes.linear.config.webEdits``)."""

    selectedOrthogroupId: str
    selectedOrthogroupAlignmentFeature: str


class SessionWebFiles(TypedDict, total=False):
    """``webFiles``: the Web's file inventory and metadata of the resources."""

    # Checked by _validate_web_file_bindings (schema 1 or 2).
    bindings: dict[str, object]


class SessionArtifactCache(TypedDict):
    """``losatCache`` and ``losatDerivedCache``."""

    # Checked by validate_current_session_artifacts.
    entries: list[object]


class SessionLegacyArtifacts(TypedDict, total=False):
    """Protein artifacts of older Sessions, kept until a Generate replaces them."""

    proteinRawCandidates: dict[str, object]
    proteinDerivedEvidence: dict[str, object]


class FeatureOverrideDraftRow(TypedDict):
    """A per-feature edit draft row of a mode slice (``features.featureOverrides``):
    a request ``featureOverrides`` row (``FeatureOverride``) plus the Web-only
    ``labelSourceText``, the label's text before any edit."""

    recordKey: str
    biologicalFeatureId: str
    featureVisibility: FeatureVisibilityMode | None
    labelVisibility: LabelVisibilityMode | None
    labelText: str | None
    labelSourceText: str | None


class _GbdrawSessionRequired(SessionResultSet):
    format: Literal["gbdraw-session"]
    version: int
    createdAt: str
    # null in a settings-only Session.
    renderRequest: CanonicalRequestJson | None
    # Resource descriptors by id (gbdraw.session_resources).
    resources: dict[str, object]
    editorState: SessionEditorState
    losatCache: SessionArtifactCache
    losatDerivedCache: SessionArtifactCache
    # gbdraw.analysis.protein_colinearity (PROTEIN_IDENTITY_MANIFEST_SCHEMA).
    proteinIdentityManifest: dict[str, object]


class GbdrawSession(_GbdrawSessionRequired, total=False):
    """A gbdraw Session of the current version (``CURRENT_SESSION_VERSION``)."""

    title: str
    webFiles: SessionWebFiles
    ui: SessionUi
    orthogroupState: SessionOrthogroupState
    # Each diagram mode's Web draft (PD-OI-086).
    modes: ModeSlicesJson
    # The CLI options a Session 27-44 draft kept in ``config``.
    cliOptions: dict[str, object]
    legacyArtifacts: SessionLegacyArtifacts
    otherModeResult: SessionOtherModeResult


class FlatDraftFields(TypedDict, total=False):
    """The flat Web draft of Sessions 44 and older; the current version keeps
    it per mode in ``modes``."""

    config: object
    features: object


class LegacySessionFields(FlatDraftFields, total=False):
    """The top-level keys that only readers of older versions take."""

    # Embedded input files of Sessions older than version 40.
    files: object


# The objects that hold keys their class does not name.
_OPEN_OBJECTS: tuple[type, ...] = (SessionEditorState, SessionUi, SessionOrthogroupState, SessionWebFiles)


def keys_of(cls: type) -> frozenset[str]:
    """The keys a declared object may hold."""

    return frozenset(get_type_hints(cls))


def describe_session_document() -> dict[str, Any]:
    """Describe the JSON of a current Session as neutral data.

    Returns ``{"root": name, "types": {name: shape}}`` with the shapes of
    :mod:`gbdraw.contract_shapes`. ``CanonicalRenderRequest`` is a reference
    into ``describe_canonical_request()``. ``SessionInput`` is what a reader
    takes before any check: every key a Session of any supported version may
    hold at the top level, with unvalidated values.
    """

    from typing import is_typeddict

    from .contract_shapes import STRING, UNKNOWN, Shape, array, obj, record, ref, typed_dict_shape

    types_: dict[str, Shape] = {}

    def named(name: str, shape: Shape) -> Shape:
        types_[name] = shape
        return ref(name)

    # A slice container holds registry keys only, each optional; the registry
    # names no value type, so a leaf is unknown unless a validator checks it.
    slice_leaves: dict[str, Shape] = {
        "config.recordsOff": array(STRING),  # _validate_records_off
    }

    def slice_container(container: str) -> Shape:
        fields: dict[str, Shape] = {}
        for key in sorted(SLICE_CONTAINERS[container]):
            child = f"{container}.{key}" if container else key
            fields[key] = slice_leaves.get(child) or (
                slice_container(child) if child in SLICE_CONTAINERS else UNKNOWN
            )
        name = "SessionModeSlice" + "".join(part[0].upper() + part[1:] for part in container.split(".") if part)
        return named(name, obj(fields, optional=list(fields)))

    def special(hint: object) -> Shape | None:
        if hint is CanonicalRequestJson:
            return ref("CanonicalRenderRequest")
        if hint is ModeSlicesJson:
            if "SessionModes" not in types_:
                # _validate_feature_override_drafts
                slice_leaves["features.featureOverrides"] = record(special(FeatureOverrideDraftRow) or UNKNOWN)
                mode_slice = slice_container("")
                named("SessionModes", obj({mode: mode_slice for mode in DIAGRAM_MODES}, optional=DIAGRAM_MODES))
            return ref("SessionModes")
        if hint is LosatExecutionJson:
            return named("SessionLosatExecution", obj(
                {field: UNKNOWN for field in LOSAT_EXECUTION_FIELDS}, optional=LOSAT_EXECUTION_FIELDS
            ))
        if isinstance(hint, type) and is_typeddict(hint):
            if hint.__name__ not in types_:
                types_[hint.__name__] = UNKNOWN  # reserves the typedef's place before its members'
                types_[hint.__name__] = typed_dict_shape(hint, special, open=hint in _OPEN_OBJECTS)
            return ref(hint.__name__)
        return None

    root = special(GbdrawSession)
    assert root is not None
    top_level = [*types_[root["name"]]["fields"], *get_type_hints(LegacySessionFields)]
    named("SessionInput", obj({key: UNKNOWN for key in top_level}, optional=top_level))
    return {"root": root["name"], "types": types_}


__all__ = [
    "FeatureOverrideDraftRow",
    "FlatDraftFields",
    "GbdrawSession",
    "LegacySessionFields",
    "SessionOtherModeEditorState",
    "SessionOtherModeLegend",
    "SessionOtherModeResult",
    "SessionOtherModeUi",
    "SessionRunMetadata",
    "describe_session_document",
    "keys_of",
]
