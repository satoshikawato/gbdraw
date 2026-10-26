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

import re
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


# The ID of a resource (``resources``) and of a drawing (``drawings[].id``):
# lowercase words joined by hyphens.
SESSION_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")
DRAWING_ID_MAX_LENGTH = 40
# A Gallery Session shows one diagram: this many of its drawings have a
# committed render (``renderRequest`` not null). tools/refresh_gallery_sessions.py
# checks the Sessions it reads and writes; the Web publication reads it generated.
GALLERY_SESSION_DRAWN_DRAWINGS = 1


class SessionProjectUi(TypedDict, total=False):
    """A Session 49 project's ``ui``: the Web's app-level view state, and the
    drawing the Web opens. A drawing's ``ui`` holds the rest."""

    activeDrawingId: str
    losatExecution: LosatExecutionJson
    downloadDpi: object
    autoLabelReflow: object
    paletteInstantPreviewEnabled: object
    richFeaturePopup: object
    featurePanelTab: object


class SessionProjectWebFiles(TypedDict, total=False):
    """A Session 49 project's ``webFiles``: what the resources were named."""

    resourceOriginalNames: dict[str, object]


class _SessionDrawingRequired(TypedDict):
    id: str
    name: str
    mode: DiagramMode
    # null for a drawing not generated yet.
    renderRequest: CanonicalRequestJson | None
    results: list[SessionResult]


class SessionDrawing(_SessionDrawingRequired, total=False):
    """One drawing of a Session 49. :mod:`gbdraw.session_drawings` maps it to
    the one-drawing Session that Sessions 46-48 give, which ``validate_session``
    checks; the parts are checked there."""

    runMetadata: SessionRunMetadata
    cliInvocation: SessionCliInvocation
    cliOptions: dict[str, object]
    orthogroupState: SessionOrthogroupState
    losatCache: SessionArtifactCache
    losatDerivedCache: SessionArtifactCache
    proteinIdentityManifest: dict[str, object]
    # The drawing's mode slice (``modes[mode]`` of Sessions 46-48).
    config: dict[str, object]
    features: dict[str, object]
    # The Result set's and the slice's editor state, and the drawing's ``ui``.
    editorState: dict[str, object]
    ui: dict[str, object]
    # Its own mode's bindings and request metadata (D-15-79, D-15-83).
    webFiles: dict[str, object]


class _GbdrawSession49Required(TypedDict):
    format: Literal["gbdraw-session"]
    version: int
    createdAt: str
    resources: dict[str, object]
    ui: SessionProjectUi
    drawings: list[SessionDrawing]


class GbdrawSession49(_GbdrawSession49Required, total=False):
    """A Session 49: the project and its drawings, in tab order (DESIGN §1)."""

    title: str
    legacyArtifacts: SessionLegacyArtifacts
    webFiles: SessionProjectWebFiles


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


def describe_session_admission() -> dict[str, Any]:
    """Describe what ``validate_session`` proves of a stored Session, per version band.

    Returns ``{"bands": [{"minVersion", "versions", "root"}], "types"}`` with
    the shapes of :mod:`gbdraw.contract_shapes`. A band starts at a threshold
    of :mod:`gbdraw.session_io` where a check of ``validate_session`` starts;
    its root is what every Session of its ``versions`` is once Python has
    accepted it. A shape keeps the keys of the writer's object
    (:func:`describe_session_document`) and declares only what Python refuses
    otherwise: a field is ``unknown`` unless a validator refuses its other
    kinds, required only where Python requires it, and closed only where
    Python refuses another key. Values, references and cross-field conditions
    stay in the validators, so a document outside its band's shape is one that
    ``validate_session`` refuses (``tests/test_session_contract.py``).
    ``defaults`` gives a required part's empty form as the writer writes it.
    """

    from . import session_io as sio
    from .contract_shapes import (
        BOOLEAN, INTEGER, STRING, UNKNOWN, Shape, array, literal, nullable, obj, record, ref,
    )

    writer = describe_session_document()["types"]
    types_: dict[str, Shape] = {}

    def named(name: str, shape: Shape) -> Shape:
        if types_.setdefault(name, shape) != shape:
            raise RuntimeError(f"Two admission shapes are named {name}.")
        return ref(name)

    def relaxed(
        writer_name: str,
        checked: dict[str, Shape] | None = None,
        *,
        required: tuple[str, ...] = (),
        open: bool = False,
        defaults: dict[str, object] | None = None,
        name: str | None = None,
    ) -> Shape:
        """The writer's object ``writer_name``: its keys, each ``unknown``
        unless ``checked`` names the shape Python proves."""
        keys = list(writer[writer_name]["fields"])
        undeclared = (set(checked or ()) | set(required)) - set(keys)
        if undeclared:
            raise RuntimeError(f"{writer_name} declares no field(s) {sorted(undeclared)}.")
        fields = {key: (checked or {}).get(key, UNKNOWN) for key in keys}
        return named(name or f"Checked{writer_name}", obj(
            fields, optional=[key for key in keys if key not in required], open=open, defaults=defaults
        ))

    # A Result: a named plain SVG (_validate_current_feature_catalog_authority).
    result = relaxed("SessionResult", {"name": STRING, "content": STRING}, required=("name", "content"), open=True)

    # The caches of Session 39+ (validate_current_session_artifacts); an entry
    # is refused unless it holds these.
    def cache(name: str, entry: dict[str, Shape]) -> Shape:
        entry_shape = named(f"Checked{name}Entry", obj(entry, open=True))
        return nullable(relaxed("SessionArtifactCache", {"entries": array(entry_shape)}, open=True, name=f"Checked{name}"))

    losat_cache = cache("LosatCache", {"kind": literal("raw-losat"), "text": STRING, "key": STRING})
    derived_cache = cache("LosatDerivedCache", {
        "kind": literal("derived-losatp-payload"), "key": STRING, "payload": record(UNKNOWN),
    })

    # The drafts a Session 46 slice holds (_validate_display_placement_draft_config,
    # _validate_feature_override_drafts); the other leaves are the Web's.
    anchor_fields: dict[str, Shape] = {
        "schema": literal(1), "recordKey": STRING, "biologicalFeatureId": STRING,
        "placement": literal("anchor", "feature-end"),
        "anchor": nullable(literal("five-prime", "midpoint", "three-prime")),
        "offsetBp": INTEGER, "orientForward": BOOLEAN,
    }
    display_fields: dict[str, Shape] = {
        "sourceUid": STRING, "selector": STRING, "recordId": STRING,
        "topologyOverride": UNKNOWN, "startCoordinate": UNKNOWN,
        "reverseComplementOverride": nullable(BOOLEAN),
        "anchorIntent": nullable(named("CheckedRecordAnchorIntent", obj(anchor_fields))),
    }
    if (set(anchor_fields) != sio.RECORD_ANCHOR_INTENT_FIELDS
            or set(display_fields) != sio.RECORD_DISPLAY_DRAFT_FIELDS | sio.TYPED_RECORD_DISPLAY_DRAFT_FIELDS):
        raise RuntimeError("The record display draft shapes name other fields than the validator.")
    slice_leaves: dict[str, Shape] = {
        "config.recordsOff": array(STRING),  # _validate_records_off
        "config.recordDisplayDrafts": array(named("CheckedRecordDisplayDraft", obj(display_fields))),
        "config.featurePlacementOverrides": record(UNKNOWN),
        "config.adv.feature_overlap_tolerance_bp": INTEGER,
        "features.featureOverrides": record(relaxed(
            "FeatureOverrideDraftRow",
            {"recordKey": STRING, "biologicalFeatureId": STRING, "labelSourceText": nullable(STRING)},
            required=tuple(writer["FeatureOverrideDraftRow"]["fields"]),
        )),
    }

    def mode_slice(writer_name: str, container: str = "") -> Shape:
        """A slice container: the registry's keys only, each optional (_validate_mode_slice_fields)."""
        checked: dict[str, Shape] = {}
        for key, shape in writer[writer_name]["fields"].items():
            path = f"{container}.{key}" if container else key
            if shape["kind"] == "ref" and shape["name"].startswith("SessionModeSlice"):
                checked[key] = mode_slice(shape["name"], path)
            elif path in slice_leaves:
                checked[key] = slice_leaves[path]
        return relaxed(writer_name, checked)

    other_mode_result = relaxed("SessionOtherModeResult", {  # _validate_other_mode_result
        "results": array(result),
        "editorState": relaxed(
            "SessionOtherModeEditorState", {"legend": relaxed("SessionOtherModeLegend")}, required=("featureCatalog",)
        ),
        "ui": relaxed("SessionOtherModeUi"),
        "runMetadata": relaxed("SessionRunMetadata"),
    }, required=("renderRequest", "results", "editorState"))

    supported = sorted(sio.SUPPORTED_SESSION_VERSIONS)
    minimums = [
        supported[0],
        sio.CANONICAL_SESSION_MIN_VERSION,
        sio.CURRENT_ARTIFACT_SESSION_MIN_VERSION,
        sio.CURRENT_AUTHORITY_SESSION_MIN_VERSION,
        sio.MODE_SCOPED_SESSION_MIN_VERSION,
        sio.DRAWINGS_SESSION_MIN_VERSION,
    ]
    top_level = list(writer["SessionInput"]["fields"])
    bands: list[dict[str, Any]] = []
    for index, minimum in enumerate(minimums):
        below = minimums[index + 1] if index + 1 < len(minimums) else supported[-1] + 1
        versions = [version for version in supported if minimum <= version < below]
        fields: dict[str, Shape] = {"format": literal(sio.SESSION_FORMAT), "version": literal(*versions)}
        root = f"SessionChecked{minimum}"
        if minimum >= sio.DRAWINGS_SESSION_MIN_VERSION:
            # The closed project and its drawings (_validate_session_drawings);
            # each drawing's parts are checked on its one-drawing view.
            drawing = obj(
                # Sorted: a set's order changes between runs, and the generated Web modules list fields in order.
                {key: UNKNOWN for key in sorted(keys_of(SessionDrawing))}
                | {"id": STRING, "name": STRING, "mode": literal(*DIAGRAM_MODES)},
                optional=sorted(keys_of(SessionDrawing) - {"id", "name", "mode"}),
            )
            project = {key: UNKNOWN for key in sorted(keys_of(GbdrawSession49))} | fields | {
                "resources": record(UNKNOWN),
                "drawings": array(named("SessionCheckedDrawing", drawing)),
            }
            named(root, obj(project, optional=sorted(set(project) - {"format", "version", "resources", "drawings"})))
            bands.append({"minVersion": minimum, "versions": versions, "root": root})
            continue
        required = ["format", "version"]
        defaults: dict[str, object] = {}
        if minimum < sio.MODE_SCOPED_SESSION_MIN_VERSION:
            fields["otherModeResult"] = literal()  # absent below 46
        # Below 40 `files` stays unknown: Python requires an object below 31 and
        # refuses another kind below 40, but the Web loads a Session without it
        # and one whose `files` is not an object (tests/fixtures/session_admission_vectors.json).
        if minimum >= sio.CANONICAL_SESSION_MIN_VERSION:
            fields["renderRequest"] = UNKNOWN  # the request's own reader checks it
            fields["resources"] = record(UNKNOWN)
            required += ["renderRequest", "resources"]
            defaults["resources"] = {}
        if minimum >= sio.CURRENT_ARTIFACT_SESSION_MIN_VERSION:
            fields["losatCache"] = losat_cache
            fields["losatDerivedCache"] = derived_cache
            fields["proteinIdentityManifest"] = nullable(record(UNKNOWN))
            # `legacyArtifacts` stays unknown: the Web loads one that is not an object.
        if minimum < sio.CURRENT_AUTHORITY_SESSION_MIN_VERSION:
            # Below 40 no top-level key is checked.
            named(root, obj(fields, optional=[key for key in fields if key not in required], open=True, defaults=defaults))
        else:
            # From 40 the top level is closed (_validate_current_top_level_fields),
            # without the keys of the other homes (_validate_mode_scoped_fields).
            absent = {"files"} | (
                {"modes", "cliOptions", "otherModeResult"}
                if minimum < sio.MODE_SCOPED_SESSION_MIN_VERSION
                else set(sio.FLAT_DRAFT_TOP_LEVEL_FIELDS)
            )
            fields = {key: fields.get(key, UNKNOWN) for key in top_level if key not in absent}
            fields["results"] = array(result)
            fields["editorState"] = relaxed("SessionEditorState", required=("featureCatalog",), open=True)
            fields["orthogroupState"] = relaxed("SessionOrthogroupState", open=True)
            required += ["results", "editorState"]
            defaults |= {"results": [], "editorState": {"featureCatalog": None}}
            if minimum < sio.MODE_SCOPED_SESSION_MIN_VERSION:
                fields["features"] = record(UNKNOWN)
            else:
                fields["modes"] = nullable(relaxed("SessionModes", {
                    mode: mode_slice("SessionModeSlice") for mode in DIAGRAM_MODES
                }))
                fields["cliOptions"] = record(UNKNOWN)
                fields["otherModeResult"] = other_mode_result
            named(root, obj(fields, optional=[key for key in fields if key not in required], defaults=defaults))
        bands.append({"minVersion": minimum, "versions": versions, "root": root})
    return {"bands": bands, "types": types_}


__all__ = [
    "FeatureOverrideDraftRow",
    "DRAWING_ID_MAX_LENGTH",
    "FlatDraftFields",
    "GbdrawSession",
    "GbdrawSession49",
    "SESSION_ID_RE",
    "SessionDrawing",
    "SessionProjectUi",
    "SessionProjectWebFiles",
    "LegacySessionFields",
    "SessionOtherModeEditorState",
    "SessionOtherModeLegend",
    "SessionOtherModeResult",
    "SessionOtherModeUi",
    "SessionRunMetadata",
    "describe_session_admission",
    "describe_session_document",
    "keys_of",
]
