"""Render stages: function calls whose results a Worker transaction may reuse.

This module is the one owner of render-stage keys (design STEP2 2.1-2.3).
``render_stage(name, function, *args, **kwargs)`` calls ``function`` directly
when no prepared-input cache transaction is active (CLI, Python API). In a
transaction it keys the call by its arguments and reuses a result stored by an
earlier successful transaction (the fifth layer of
``gbdraw.api.prepared.PreparedBiologicalInputCache``).

Keys follow a closed table (``_key``): plain values (with tuples, lists,
sets, dicts and read-only dict views of them) by value; a frozen dataclass by
its fields; a module-level class by its qualified name; a ``GbdrawConfig`` by the leaves of the sections (and
the single leaves) the stage declares in ``RENDER_STAGE_SETTINGS``; a configurator by its attributes;
a record, a decoded table or another stage's result by the identity registered
for it in the transaction. Anything else makes the call uncacheable for that
request.

Settings view: in a transaction the function receives a copy of each config
whose undeclared sections and leaves flag a read. A view that outlives its
build (the result, or anything else, holds a part of it) also flags the build,
and a copy of a view section reads the section: a stored result never hands a
later reader the sections its key left out. A flagged result is returned but
not stored, gets no identity, and is counted. Views made in a nested stage count for every
enclosing stage too. A view kept alive only by a reference cycle until the collector runs also makes
the build unstored (a lost store, never a stale hit).
"""

from __future__ import annotations

import copy
import dataclasses
import inspect
import logging
import re
import typing
import weakref
from contextvars import ContextVar
from collections import defaultdict
from dataclasses import dataclass, field
from types import MappingProxyType
from enum import Enum
from typing import Any, Callable, Hashable, Iterable, Literal, Mapping, ParamSpec, Sequence, TypeVar, cast

from Bio.SeqRecord import SeqRecord

from gbdraw.config.models.root import GbdrawConfig

logger = logging.getLogger(__name__)

P = ParamSpec("P")
R = TypeVar("R")

RenderStageName = Literal[
    "color_rules",
    "default_colors",
    "visibility_rules",
    "label_filter_maps",
    "gc_skew_data",
    "depth_frame",
    "depth_track_config",
    "legend_table",
    "legend_measurement",
    "legend_content",
    "plot_title_content",
    "circular_plot_title_content",
    "composition_plan",
    "gc_content_content",
    "gc_skew_content",
    "depth_content",
    "length_bar_content",
    "definition_content",
    "circular_definition_content",
    "annotation_content",
    "definition_track_overlaps",
    "record_content",
    "circular_record_content",
    "comparison_content",
    "interactive_context",
    "catalog_normalization",
    "catalog_snapshot",
    "catalog_item",
    "comparison_rows",
    "protein_extraction",
    "protein_analysis",
    "collinearity_projection",
    "derived_losat_entries",
    "feature_layers",
    "feature_dom_index",
    "circular_labels",
    "circular_radial_labels",
    "axis_content",
    "ticks_content",
    "conservation_sources",
    "conservation_tracks",
    "conservation_content",
    "label_leaders_content",
    "label_text_content",
    "record_part",
    "record_merge",
    "record_gc_percent",
]

# The sections of GbdrawConfig each stage reads (design 2.2). A section is a
# top node ("canvas", "labels", "objects": its direct values) or one of the
# dataclasses below it ("labels.filtering", "objects.features", ...). A stage
# that reads one direct value of a top node declares that leaf alone
# ("canvas.strandedness"), so a change of the node's other values (the track
# switches in ``canvas``) does not rebuild it.
_CIRCULAR_LABEL_SECTIONS = frozenset(
    {
        "canvas",
        "canvas.circular",
        "labels",
        "labels.font_size",
        "labels.length_threshold",
        "labels.spacing",
        "labels.unified_adjustment",
        "objects.features",
        "objects.text",
    }
)
# The Circular tick group sizes itself from the length threshold and the
# preset track table, and places its ring by the strand setting.
_CIRCULAR_TRACK_SECTIONS = frozenset({"canvas.circular", "canvas.strandedness", "labels.length_threshold"})
RENDER_STAGE_SETTINGS: Mapping[RenderStageName, frozenset[str]] = {
    "color_rules": frozenset(),
    "default_colors": frozenset(),
    "visibility_rules": frozenset(),
    "label_filter_maps": frozenset(),
    "gc_skew_data": frozenset(),
    "depth_frame": frozenset(),
    "depth_track_config": frozenset(),
    "legend_table": frozenset(),
    "legend_measurement": frozenset(),
    "legend_content": frozenset(),
    "plot_title_content": frozenset(),
    # The Circular title and definition take their font and interval as a
    # ``CircularDefinitionStyle`` argument, not the config (OV-1115).
    "circular_plot_title_content": frozenset(),
    "composition_plan": frozenset(),
    # The Circular track groups take their norm factor from the radial plan.
    "gc_content_content": frozenset(),
    "gc_skew_content": frozenset(),
    "depth_content": frozenset(),
    "length_bar_content": frozenset(),
    "definition_content": frozenset(),
    "circular_definition_content": frozenset(),
    "annotation_content": frozenset(),
    # The centre definition lines over the Circular track paint: the text font, the conservation background and
    # the tick ring (the tick group's sections) are the settings that paint reads (OV-2088).
    "definition_track_overlaps": _CIRCULAR_TRACK_SECTIONS
    | {"objects.conservation", "objects.scale", "objects.text", "objects.ticks"},
    # Linear: the record group draws the axis, the ruler, features (on one or two strands) and labels (with
    # their style).
    "record_content": frozenset(
        {
            "canvas.strandedness",
            "labels.linear",
            "labels.stroke_color",
            "objects.axis",
            "objects.scale",
            "objects.text",
        }
    ),
    # Circular: the record group sizes features (on one or two strands) and
    # embedded labels (the label drawer reads the preset track factors).
    "circular_record_content": frozenset(
        {
            "canvas.circular",
            "canvas.strandedness",
            "labels.font_size",
            "labels.length_threshold",
            "objects.text",
        }
    ),
    "comparison_content": frozenset(),
    "interactive_context": frozenset(),
    "catalog_normalization": frozenset(),
    "catalog_snapshot": frozenset(),
    "catalog_item": frozenset(),
    "comparison_rows": frozenset(),
    "protein_extraction": frozenset(),
    "protein_analysis": frozenset(),
    "collinearity_projection": frozenset(),
    "derived_losat_entries": frozenset(),
    "feature_layers": frozenset(),
    "feature_dom_index": frozenset(),
    # One function (CA ``circular_labels``) under two names: radial labels also
    # replan the radial layout as they grow the radius, which reads the axis
    # and tick sections; horizontal labels never do.
    "circular_labels": _CIRCULAR_LABEL_SECTIONS,
    "circular_radial_labels": _CIRCULAR_LABEL_SECTIONS | {"objects.axis", "objects.scale", "objects.ticks"},
    "axis_content": frozenset({"objects.axis"}),
    "ticks_content": _CIRCULAR_TRACK_SECTIONS | {"objects.scale", "objects.text", "objects.ticks"},
    "conservation_sources": frozenset(),
    "conservation_tracks": frozenset(),
    "conservation_content": frozenset({"objects.conservation"}),
    # The leaders end on their underlines, which the text layer draws with the leader's stroke.
    "label_leaders_content": frozenset({"labels.circular", "labels.stroke_color"}),
    # The text layer places embedded labels on their strand's lane.
    "label_text_content": frozenset(
        {
            "canvas.strandedness",
            "labels.circular",
            "labels.font_size",
            "labels.length_threshold",
            "labels.stroke_color",
            "objects.text",
        }
    ),
    "record_part": frozenset(),
    "record_merge": frozenset(),
    "record_gc_percent": frozenset(),
}

_TOP_SECTIONS = ("canvas", "labels", "objects")


def _validate_stage_settings(table: Mapping[RenderStageName, frozenset[str]]) -> None:
    """Every stage is a ``RenderStageName`` and every section names a config node or a top node's leaf."""

    stage_names = frozenset(typing.get_args(RenderStageName))
    top_hints = typing.get_type_hints(GbdrawConfig)
    for stage, sections in table.items():
        if stage not in stage_names:
            raise ValueError(f"Unknown render stage {stage!r} in RENDER_STAGE_SETTINGS.")
        for section in sections:
            top, _dot, child = section.partition(".")
            top_type = top_hints.get(top) if top in _TOP_SECTIONS else None
            known = top_type is not None and (not child or child in typing.get_type_hints(top_type))
            if not known or "." in child:
                raise ValueError(f"Render stage {stage!r} declares unknown section {section!r}.")


_validate_stage_settings(RENDER_STAGE_SETTINGS)
_PLAIN_TYPES = (bool, int, float, str, bytes)
# ``gbdraw.api.prepared`` is imported inside the functions that use it (about
# 0.2 us once loaded): feature, label, depth and legend modules call
# ``render_stage`` and are imported while ``gbdraw.api`` initializes, so this
# module cannot import the ``gbdraw.api`` package when it is loaded.


class _Unkeyable(Exception):
    """An argument has no key; the call is not cached in this request."""


_T = TypeVar("_T")


class _PreKeyed(tuple[_T, ...]):
    """A tuple keyed once, for an argument that several stages take unchanged.

    Built by :func:`_pre_keyed`; ``key`` is the key of the equal plain tuple,
    or ``None`` outside a cache transaction or when an item has no key.
    """

    key: Hashable | None

    def __new__(cls, items: Iterable[_T]) -> _PreKeyed[_T]:
        keyed = super().__new__(cls, items)
        keyed.key = None
        return keyed


def _pre_keyed(items: Iterable[_T]) -> _PreKeyed[_T]:
    """Key ``items`` once, so each stage that takes them reads the key.

    The items must be values nobody else holds (copies), so they cannot change
    after the key is computed. Outside a cache transaction nothing is keyed.
    """

    from gbdraw.api import prepared

    keyed = _PreKeyed(items)
    if prepared.active_render_stage_transaction() is not None:
        try:
            keyed.key = _key(tuple(keyed), config_path=None)
        except _Unkeyable:
            # The stage meets the same item and takes its uncacheable path.
            pass
    return keyed


@dataclass(eq=False)
class _StageFrame:
    """One stage build in progress: its declarations, its views, the inner stage keys and the read flag."""

    name: str
    declared: frozenset[str]
    active: bool = True
    undeclared_sections: set[str] = field(default_factory=set)
    inner_keys: set[Hashable] = field(default_factory=set)
    # The views built for this stage (each config and section copy), held weakly.
    views: list[weakref.ref[object]] = field(default_factory=list)


_FRAMES: ContextVar[tuple[_StageFrame, ...]] = ContextVar(
    "gbdraw_render_stage_frames",
    default=(),
)


# --- settings view ----------------------------------------------------------


_FLAG_STATE = "_gbdraw_render_stage_flag"
_FLAGGING_CLASSES: dict[tuple[type, bool], type] = {}


@dataclass(frozen=True)
class _FlagState:
    section: str
    source: object
    frames: tuple[_StageFrame, ...]
    top: bool
    field_names: frozenset[str]

    # A copy of a view (copy.copy, copy.deepcopy) keeps flagging the live
    # frames of the stages that built the view.
    def __copy__(self) -> _FlagState:
        return self

    def __deepcopy__(self, memo: dict[int, object]) -> _FlagState:
        return self


def _flag_state(value: object) -> _FlagState | None:
    try:
        state = object.__getattribute__(value, "__dict__").get(_FLAG_STATE)
    except AttributeError:
        return None
    return state if isinstance(state, _FlagState) else None


def _mark(frames: tuple[_StageFrame, ...], section: str) -> None:
    """Flag a read of ``section`` (a section or a top node's leaf) in each build that did not declare it."""

    for frame in frames:
        if frame.active and section not in frame.declared:
            frame.undeclared_sections.add(section)


def _flagging_class(cls: type, top: bool) -> type:
    sub = _FLAGGING_CLASSES.get((cls, top))
    if sub is not None:
        return sub

    def __getattribute__(self: object, name: str) -> object:
        value = object.__getattribute__(self, name)
        if name.startswith("__") and name != "__dict__":
            return value
        state = object.__getattribute__(self, "__dict__").get(_FLAG_STATE)
        if not isinstance(state, _FlagState):
            return value
        if not state.top or name == "__dict__":
            # ``vars()`` hands over every value of the node.
            _mark(state.frames, state.section)
        elif name in state.field_names and not _is_dataclass_value(value):
            # A top node flags its direct values, each as its leaf; a child section flags itself.
            _mark(state.frames, f"{state.section}.{name}")
        return value

    def __eq__(self: object, other: object) -> bool:
        # A comparison reads the section; it compares the sections themselves.
        for value in (self, other):
            state = _flag_state(value)
            if state is not None:
                _mark(state.frames, state.section)
        return bool(_section_source(self) == _section_source(other))

    def __ne__(self: object, other: object) -> bool:
        return not __eq__(self, other)

    def __hash__(self: object) -> int:
        state = _flag_state(self)
        if state is not None:
            _mark(state.frames, state.section)
        return hash(_section_source(self))

    def __reduce_ex__(self: object, protocol: int) -> object:
        # A copy (copy.copy, copy.deepcopy) holds every value of the section
        # without keeping the view alive: making it reads the section.
        state = _flag_state(self)
        if state is not None:
            _mark(state.frames, state.section)
        return object.__reduce_ex__(self, protocol)

    sub = type(
        cls.__name__,
        (cls,),
        {
            "__getattribute__": __getattribute__,
            "__eq__": __eq__,
            "__ne__": __ne__,
            "__hash__": __hash__,
            "__reduce_ex__": __reduce_ex__,
        },
    )
    _FLAGGING_CLASSES[(cls, top)] = sub
    return sub


def _is_dataclass_value(value: object) -> bool:
    return dataclasses.is_dataclass(value) and not isinstance(value, type)


def _raw_fields(value: object) -> dict[str, object]:
    """Read every field without flagging (building a view is not a read)."""

    return {
        config_field.name: object.__getattribute__(value, config_field.name)
        for config_field in dataclasses.fields(cast(Any, value))
    }


def _section_source(value: object) -> object:
    state = _flag_state(value)
    return value if state is None else state.source


def _flagging_copy(
    value: object,
    section: str,
    frames: tuple[_StageFrame, ...],
    *,
    top: bool,
    overrides: Mapping[str, object] | None = None,
) -> object:
    state = _flag_state(value)
    source = _section_source(value)
    enclosing = () if state is None else state.frames
    values = _raw_fields(value)
    values.update(overrides or {})
    new: object = object.__new__(_flagging_class(type(source), top))
    for name, child in values.items():
        object.__setattr__(new, name, child)
    object.__getattribute__(new, "__dict__")[_FLAG_STATE] = _FlagState(
        section=section,
        source=source,
        frames=enclosing + frames,
        top=top,
        field_names=frozenset(values),
    )
    return new


def _plain_copy(value: object, overrides: Mapping[str, object]) -> object:
    """Copy a node with replaced children; keep the flags it already has."""

    state = _flag_state(value)
    new = copy.copy(value)
    for name, child in overrides.items():
        object.__setattr__(new, name, child)
    if state is not None:
        object.__getattribute__(new, "__dict__")[_FLAG_STATE] = state
    return new


def _settings_view(
    config: GbdrawConfig,
    declared: frozenset[str],
    frame: _StageFrame,
) -> GbdrawConfig:
    """Copy ``config`` so that every section not in ``declared`` flags ``frame``."""

    frames = (frame,)
    top_values: dict[str, object] = {}
    for top_name, node in _raw_fields(config).items():
        children: dict[str, object] = {}
        for child_name, child in _raw_fields(node).items():
            if not _is_dataclass_value(child):
                continue
            section = f"{top_name}.{child_name}"
            if section in declared:
                continue
            children[child_name] = _flagging_copy(child, section, frames, top=False)
        if top_name in declared:
            top_values[top_name] = _plain_copy(node, children)
        else:
            top_values[top_name] = _flagging_copy(
                node,
                top_name,
                frames,
                top=True,
                overrides=children,
            )
        frame.views.extend(weakref.ref(copied) for copied in (top_values[top_name], *children.values()))
    view = copy.copy(config)
    for top_name, node in top_values.items():
        object.__setattr__(view, top_name, node)
    frame.views.append(weakref.ref(view))
    return view


def _section_leaves(config: GbdrawConfig, section: str) -> Hashable:
    """Key the leaves of one section, or one leaf; reading through a view flags its frames."""

    top_name, _dot, child_name = section.partition(".")
    node = object.__getattribute__(config, top_name)
    if child_name:
        # Through a top node's view, the attribute read flags a leaf.
        node = getattr(node, child_name)
        if not _is_dataclass_value(node):
            return _key(node, config_path=None)
    state = _flag_state(node)
    if state is not None:
        _mark(state.frames, state.section)
    values = _raw_fields(_section_source(node))
    return tuple(
        (name, _key(value, config_path=None))
        for name, value in values.items()
        # A top node's own leaves; its child sections are keyed when declared.
        if child_name or not _is_dataclass_value(value)
    )


def _config_key(config: GbdrawConfig, declared: frozenset[str]) -> Hashable:
    return ("config", tuple((section, _section_leaves(config, section)) for section in sorted(declared)))


# --- keys -------------------------------------------------------------------


_ArgumentPath = tuple[int | str, ...]
_ConfigPaths = list[tuple[_ArgumentPath, GbdrawConfig]]
_CONFIGURATOR_TYPES: tuple[type, ...] | None = None


def _configurator_types() -> tuple[type, ...]:
    global _CONFIGURATOR_TYPES
    if _CONFIGURATOR_TYPES is None:
        from gbdraw.configurators import (
            BlastMatchConfigurator,
            DepthConfigurator,
            FeatureDrawingConfigurator,
            GcContentConfigurator,
            GcSkewConfigurator,
            LegendDrawingConfigurator,
        )

        _CONFIGURATOR_TYPES = (
            FeatureDrawingConfigurator,
            GcContentConfigurator,
            GcSkewConfigurator,
            DepthConfigurator,
            LegendDrawingConfigurator,
            BlastMatchConfigurator,
        )
    return _CONFIGURATOR_TYPES


def _type_name(value: object) -> str:
    cls = type(value)
    return f"{cls.__module__}.{cls.__qualname__}"


def _key(
    value: object,
    *,
    config_path: _ArgumentPath | None,
    configs: _ConfigPaths | None = None,
    declared: frozenset[str] = frozenset(),
) -> Hashable:
    """Key one value by the closed table of design 2.1.

    ``config_path`` is the field path from the stage argument while only
    frozen dataclasses were crossed; a config is keyable only there.
    """

    if value is None:
        return None
    if type(value) is _PreKeyed:
        if value.key is not None:
            return value.key
        value = tuple(value)
    if isinstance(value, Enum):
        return ("enum", _type_name(value), value.value)
    if isinstance(value, float):
        # hex() tells -0.0 from 0.0, which compare equal.
        return (type(value).__name__, value.hex())
    if isinstance(value, _PLAIN_TYPES):
        return (type(value).__name__, value)
    if isinstance(value, type) and _is_keyable_function(value):
        # A module-level class a stage builds with, keyed like the stage function.
        return ("class", f"{value.__module__}.{value.__qualname__}")
    state = _flag_state(value)
    if state is not None:
        # A section of a view passed on by value: keying it reads it, so the
        # stages whose view it is are marked, and the section itself is keyed.
        _mark(state.frames, state.section)
        value = state.source
    if isinstance(value, GbdrawConfig):
        # Before the identity lookup: a config is keyed by its declared
        # sections only, whatever identity a stage result gave it.
        if config_path is None or configs is None:
            raise _Unkeyable("a GbdrawConfig outside an argument or a frozen dataclass field")
        configs.append((config_path, value))
        return _config_key(value, declared)
    from gbdraw.api import prepared

    identity = prepared.prepared_value_identity(value)
    if identity is prepared.NO_IDENTITY:
        raise _Unkeyable(f"unstored stage result {_type_name(value)}")
    if identity is not None:
        return ("identity", identity)
    if isinstance(value, SeqRecord):
        raise _Unkeyable("a record without a registered identity")
    if isinstance(value, defaultdict):
        raise _Unkeyable("a defaultdict (its default factory is state)")
    if isinstance(value, (tuple, list)):
        return (
            _type_name(value),
            tuple(_key(item, config_path=None) for item in value),
        )
    if isinstance(value, (set, frozenset)):
        return (
            _type_name(value),
            tuple(sorted((_key(item, config_path=None) for item in value), key=repr)),
        )
    if isinstance(value, (dict, MappingProxyType)):
        return (
            _type_name(value),
            tuple(
                sorted(
                    (
                        (_key(name, config_path=None), _key(item, config_path=None))
                        for name, item in value.items()
                    ),
                    key=repr,
                )
            ),
        )
    if isinstance(value, re.Pattern):
        return ("pattern", value.pattern, value.flags)
    if _is_dataclass_value(value) and getattr(type(value), "__dataclass_params__").frozen:
        return (
            "dataclass",
            _type_name(value),
            tuple(
                (
                    name,
                    _key(
                        item,
                        config_path=None if config_path is None else (*config_path, name),
                        configs=configs,
                        declared=declared,
                    ),
                )
                for name, item in _raw_fields(value).items()
            ),
        )
    if isinstance(value, _configurator_types()):
        return (
            "configurator",
            _type_name(value),
            tuple(
                sorted(
                    (name, _key(item, config_path=None))
                    for name, item in vars(value).items()
                )
            ),
        )
    raise _Unkeyable(f"no key for {_type_name(value)}")


def _with_views(
    value: object,
    path: _ArgumentPath,
    views: Mapping[_ArgumentPath, GbdrawConfig],
) -> object:
    """Rebuild one argument with each config at ``views`` replaced by its view."""

    if path in views:
        return views[path]
    if not any(view_path[: len(path)] == path for view_path in views):
        return value
    new = copy.copy(value)
    for name, child in _raw_fields(value).items():
        replaced = _with_views(child, (*path, name), views)
        if replaced is not child:
            object.__setattr__(new, name, replaced)
    return new


def _is_keyable_function(function: Callable[..., object]) -> bool:
    if isinstance(function, type):
        return function.__qualname__ == function.__name__
    return (
        inspect.isfunction(function)
        and function.__closure__ is None
        and function.__name__ != "<lambda>"
        and function.__qualname__ == function.__name__
    )


def _direct_children(value: object) -> Iterable[tuple[Hashable, object]]:
    if _is_dataclass_value(value):
        return _raw_fields(value).items()
    if isinstance(value, tuple) and hasattr(value, "_fields"):
        return zip(cast(Any, value)._fields, value)
    if isinstance(value, (tuple, list)):
        return enumerate(value)
    if isinstance(value, dict):
        return ((name, child) for name, child in value.items() if isinstance(name, _PLAIN_TYPES))
    return ()


# Where a stored result's value came from among the stage's arguments: the
# argument (position or keyword) and, for a direct child of it, the child name.
_ArgumentLocator = tuple[int | str, Hashable | None]


def _argument_locators(
    args: tuple[object, ...],
    kwargs: Mapping[str, object],
) -> dict[int, _ArgumentLocator]:
    """The arguments and their direct children: a result may pass these through."""

    locators: dict[int, _ArgumentLocator] = {}
    positions: list[tuple[int | str, object]] = [*enumerate(args), *kwargs.items()]
    for position, arg in positions:
        locators.setdefault(id(arg), (position, None))
        for name, child in _direct_children(arg):
            locators.setdefault(id(child), (position, name))
    return locators


def _locate(
    locator: _ArgumentLocator,
    args: tuple[object, ...],
    kwargs: Mapping[str, object],
) -> object:
    position, name = locator
    value = args[position] if isinstance(position, int) else kwargs[position]
    if name is None:
        return value
    return dict(_direct_children(value))[name]


@dataclass(frozen=True)
class _StoredStage:
    """A stored result and the parts of it that were arguments of its build.

    ``passed_through`` maps ``None`` (the result itself) or a direct child name
    to the argument it was. On a hit those parts take the identity of the
    current request's argument, as they kept their argument's identity on the
    build, so a consumer gets the same key on a build and on a hit.
    """

    result: object
    passed_through: tuple[tuple[Hashable | None, _ArgumentLocator], ...]


def _passed_through(
    result: object,
    locators: Mapping[int, _ArgumentLocator],
) -> tuple[tuple[Hashable | None, _ArgumentLocator], ...]:
    parts: list[tuple[Hashable | None, _ArgumentLocator]] = []
    if id(result) in locators:
        parts.append((None, locators[id(result)]))
    for name, child in _direct_children(result):
        if id(child) in locators:
            parts.append((name, locators[id(child)]))
    return tuple(parts)


def _plain(value: object) -> bool:
    """A value that needs no identity; empty tuples and frozensets are shared singletons."""

    return (
        value is None
        or isinstance(value, (*_PLAIN_TYPES, Enum))
        or (type(value) in (tuple, frozenset) and not value)
    )


def _register_result(value: object, identity: Hashable, passed: Mapping[int, object]) -> None:
    """Give a result and its direct fields the stage key, or NO_IDENTITY.

    A value the stage received (an argument or a direct child of one) keeps the
    key it was given as an argument, so its key does not depend on whether a
    stage has returned it.
    """

    if _plain(value) or id(value) in passed:
        return
    from gbdraw.api import prepared

    prepared.register_prepared_value_identity(value, identity)
    for name, child in _direct_children(value):
        if _plain(child) or id(child) in passed:
            continue
        prepared.register_prepared_value_identity(
            child,
            prepared.NO_IDENTITY if identity is prepared.NO_IDENTITY else ("field", identity, name),
        )


def _register_hit(
    stored: _StoredStage,
    key: Hashable,
    args: tuple[object, ...],
    kwargs: Mapping[str, object],
) -> None:
    """Register a hit's result as its build registered it (see ``_StoredStage``)."""

    from gbdraw.api import prepared

    passed = dict(stored.passed_through)

    def register(value: object, part: Hashable | None, identity: Hashable) -> None:
        if _plain(value):
            return
        if part in passed:
            current = prepared.prepared_value_identity(_locate(passed[part], args, kwargs))
            if current is not None:
                prepared.register_prepared_value_identity(value, current)
            return
        prepared.register_prepared_value_identity(value, identity)

    result = stored.result
    register(result, None, key)
    for name, child in _direct_children(result):
        if None in passed and name not in passed:
            continue
        register(child, name, ("field", key, name))


def _metric(kind: str, name: str) -> None:
    from gbdraw.api import prepared

    prepared.record_prepared_input_metric(f"renderStage{kind}Count")
    prepared.record_prepared_input_metric(f"renderStage{kind}:{name}")


def render_stage(
    name: RenderStageName,
    function: Callable[P, R],
    /,
    *args: P.args,
    **kwargs: P.kwargs,
) -> R:
    """Call ``function``; in a cache transaction, reuse a stored result by key."""

    if name not in RENDER_STAGE_SETTINGS:
        raise ValueError(f"Unknown render stage {name!r}.")
    from gbdraw.api import prepared

    transaction = prepared.active_render_stage_transaction()
    if transaction is None:
        return function(*args, **kwargs)
    declared = RENDER_STAGE_SETTINGS[name]
    configs: _ConfigPaths = []
    try:
        if not _is_keyable_function(function):
            raise _Unkeyable(f"stage function {function!r} is not module-level")
        key: Hashable = (
            "render-stage",
            name,
            f"{function.__module__}.{function.__qualname__}",
            transaction.font_environment,
            tuple(
                _key(arg, config_path=(index,), configs=configs, declared=declared)
                for index, arg in enumerate(args)
            ),
            tuple(
                sorted(
                    (
                        argument,
                        _key(value, config_path=(argument,), configs=configs, declared=declared),
                    )
                    for argument, value in kwargs.items()
                )
            ),
        )
    except _Unkeyable as reason:
        logger.debug("render stage %s is not cached: %s", name, reason)
        _metric("Uncacheable", name)
        # It may take a record without an identity (O-2).
        prepared.mark_record_touched(transaction)
        prepared.cancel_checkpoint()
        result = function(*args, **kwargs)
        _register_result(result, prepared.NO_IDENTITY, _argument_locators(args, kwargs))
        return result

    enclosing = _FRAMES.get()
    entry = prepared.lookup_render_stage(transaction, key)
    if entry is not None:
        _metric("Hit", name)
        for frame in enclosing:
            frame.inner_keys.add(key)
            frame.inner_keys.update(entry.inner_keys)
        stored = cast(_StoredStage, entry.value)
        _register_hit(stored, key, args, kwargs)
        return cast(R, stored.result)

    # A canceled Worker request stops before the build, never inside it (D7).
    prepared.cancel_checkpoint()
    if prepared.key_names_record(key):
        prepared.mark_record_touched(transaction)
    frame = _StageFrame(name, declared)
    views = {path: _settings_view(config, declared, frame) for path, config in configs}
    call_args = tuple(_with_views(arg, (index,), views) for index, arg in enumerate(args))
    call_kwargs = {
        argument: _with_views(value, (argument,), views) for argument, value in kwargs.items()
    }
    token = _FRAMES.set((*enclosing, frame))
    try:
        result = cast(Callable[..., R], function)(*call_args, **call_kwargs)
    finally:
        frame.active = False
        # An enclosing stage may keep what this build made, by its return value or an error.
        for outer in enclosing:
            outer.views.extend(frame.views)
        _FRAMES.reset(token)
    _metric("Build", name)
    # The views and profile copies the function received count as arguments.
    locators = {**_argument_locators(call_args, call_kwargs), **_argument_locators(args, kwargs)}
    # Released, the views live on only where the build left them (in the result, at any depth, OV-2105).
    del views, call_args, call_kwargs
    if any(view() is not None for view in frame.views):
        frame.undeclared_sections.add("(a settings view held after the build)")
    if frame.undeclared_sections:
        logger.debug(
            "render stage %s read undeclared sections %s",
            name,
            sorted(frame.undeclared_sections),
        )
        _metric("UndeclaredRead", name)
        _register_result(result, prepared.NO_IDENTITY, locators)
        return result
    for outer in enclosing:
        outer.inner_keys.add(key)
        outer.inner_keys.update(frame.inner_keys)
    prepared.store_render_stage(
        transaction,
        key,
        _StoredStage(result, _passed_through(result, _argument_locators(args, kwargs))),
        frozenset(frame.inner_keys),
    )
    _register_result(result, key, locators)
    return result


def render_stage_revision(name: RenderStageName, results: Sequence[object]) -> str | None:
    """Name the builds of ``results``, the stage ``name`` results of this transaction, in order.

    The revision is the cache epoch and the build number of each result
    (gbdraw/web/CLAUDE.md, Runtime data flow): equal revisions name the same
    builds in one cache. ``None`` when a result was not stored under a
    ``name`` key (uncacheable, an undeclared read, or not built by the stage),
    outside a transaction, or when the cache has no epoch.
    """

    from gbdraw.api import prepared

    transaction = prepared.active_render_stage_transaction()
    if transaction is None or not results:
        return None
    epoch: str | None = None
    build_ids: list[str] = []
    for result in results:
        key = prepared.prepared_value_identity(result)
        if not (isinstance(key, tuple) and key[:2] == ("render-stage", name)):
            return None
        build = prepared.render_stage_build(transaction, key)
        if build is None:
            return None
        epoch = build[0]
        build_ids.append(str(build[1]))
    return f"{epoch}:{'.'.join(build_ids)}"


__all__ = ["RENDER_STAGE_SETTINGS", "RenderStageName", "render_stage", "render_stage_revision"]
