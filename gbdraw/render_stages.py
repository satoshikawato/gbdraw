"""Render stages: function calls whose results a Worker transaction may reuse.

This module is the one owner of render-stage keys (design STEP2 2.1-2.3).
``render_stage(name, function, *args, **kwargs)`` calls ``function`` directly
when no prepared-input cache transaction is active (CLI, Python API). In a
transaction it keys the call by its arguments and reuses a result stored by an
earlier successful transaction (the fifth layer of
``gbdraw.api.prepared.PreparedBiologicalInputCache``).

Keys follow a closed table (``_key``): plain values by value; a frozen
dataclass by its fields; a ``GbdrawConfig`` by the leaves of the sections the
stage declares in ``RENDER_STAGE_SETTINGS``; a configurator by its attributes;
a record, a decoded table or another stage's result by the identity registered
for it in the transaction. Anything else makes the call uncacheable for that
request.

Settings view: in a transaction the function receives a copy of each config
whose undeclared sections flag a read. A flagged result is returned but not
stored, gets no identity, and is counted.
"""

from __future__ import annotations

import copy
import dataclasses
import inspect
import logging
import re
from contextvars import ContextVar
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Hashable, Iterable, Literal, Mapping, ParamSpec, TypeVar, cast

from Bio.SeqRecord import SeqRecord

from gbdraw.api.prepared import (
    NO_IDENTITY,
    active_render_stage_transaction,
    lookup_render_stage,
    prepared_value_identity,
    record_prepared_input_metric,
    register_prepared_value_identity,
    store_render_stage,
)
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
    "depth_track_data",
    "legend_table",
    "legend_measurement",
    "interactive_context",
    "catalog_normalization",
    "comparison_rows",
    "protein_extraction",
    "derived_losat_entries",
    "protein_comparison_analysis",
    "feature_layers",
    "feature_dom_index",
    "linear_labels",
    "definition_metrics",
]

# The sections of GbdrawConfig each stage reads (design 2.2). A section is a
# top node ("canvas", "labels", "objects": its direct values) or one of the
# dataclasses below it ("labels.filtering", "objects.features", ...).
RENDER_STAGE_SETTINGS: Mapping[RenderStageName, frozenset[str]] = {
    "color_rules": frozenset(),
    "default_colors": frozenset(),
    "visibility_rules": frozenset(),
    "label_filter_maps": frozenset(),
    "gc_skew_data": frozenset(),
    "depth_frame": frozenset(),
    "depth_track_data": frozenset(),
    "legend_table": frozenset(),
    "legend_measurement": frozenset(),
    "interactive_context": frozenset(),
    "catalog_normalization": frozenset(),
    "comparison_rows": frozenset(),
    "protein_extraction": frozenset(),
    "derived_losat_entries": frozenset(),
    "protein_comparison_analysis": frozenset(),
    "feature_layers": frozenset(),
    "feature_dom_index": frozenset(),
    "linear_labels": frozenset(),
    "definition_metrics": frozenset(),
}

_TOP_SECTIONS = ("canvas", "labels", "objects")
_PLAIN_TYPES = (bool, int, float, str, bytes)


class _Unkeyable(Exception):
    """An argument has no key; the call is not cached in this request."""


@dataclass(eq=False)
class _StageFrame:
    """One stage build in progress: the inner stage keys and the read flag."""

    name: str
    active: bool = True
    undeclared_sections: set[str] = field(default_factory=set)
    inner_keys: set[Hashable] = field(default_factory=set)


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


def _flag_state(value: object) -> _FlagState | None:
    try:
        state = object.__getattribute__(value, "__dict__").get(_FLAG_STATE)
    except AttributeError:
        return None
    return state if isinstance(state, _FlagState) else None


def _mark(frames: tuple[_StageFrame, ...], section: str) -> None:
    for frame in frames:
        if frame.active:
            frame.undeclared_sections.add(section)


def _flagging_class(cls: type, top: bool) -> type:
    sub = _FLAGGING_CLASSES.get((cls, top))
    if sub is not None:
        return sub

    def __getattribute__(self: object, name: str) -> object:
        value = object.__getattribute__(self, name)
        if name.startswith("__"):
            return value
        state = object.__getattribute__(self, "__dict__").get(_FLAG_STATE)
        if not isinstance(state, _FlagState):
            return value
        if state.top and (name not in state.field_names or _is_dataclass_value(value)):
            # A top node flags its direct values; a child section flags itself.
            return value
        _mark(state.frames, state.section)
        return value

    sub = type(cls.__name__, (cls,), {"__getattribute__": __getattribute__})
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
    view = copy.copy(config)
    for top_name, node in top_values.items():
        object.__setattr__(view, top_name, node)
    return view


def _section_leaves(config: GbdrawConfig, section: str) -> Hashable:
    """Key the leaves of one section; reading through a view flags its frames."""

    path = section.split(".")
    node: object = config
    for part in path:
        node = object.__getattribute__(node, part)
        state = _flag_state(node)
        if state is not None and not state.top:
            _mark(state.frames, state.section)
    state = _flag_state(node)
    source = _section_source(node)
    if len(path) == 1:
        if state is not None:
            _mark(state.frames, state.section)
        values = {
            name: value
            for name, value in _raw_fields(source).items()
            if not _is_dataclass_value(value)
        }
    else:
        values = _raw_fields(source)
    return tuple((name, _key(value, config_path=None)) for name, value in values.items())


def _config_key(config: GbdrawConfig, declared: frozenset[str]) -> Hashable:
    return ("config", tuple((section, _section_leaves(config, section)) for section in sorted(declared)))


# --- keys -------------------------------------------------------------------


_ConfigPaths = list[tuple[tuple[Any, ...], GbdrawConfig]]
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
    config_path: tuple[Any, ...] | None,
    configs: _ConfigPaths | None = None,
    declared: frozenset[str] = frozenset(),
) -> Hashable:
    """Key one value by the closed table of design 2.1.

    ``config_path`` is the field path from the stage argument while only
    frozen dataclasses were crossed; a config is keyable only there.
    """

    if value is None:
        return None
    if isinstance(value, Enum):
        return ("enum", _type_name(value), value.value)
    if isinstance(value, _PLAIN_TYPES):
        return (type(value).__name__, value)
    identity = prepared_value_identity(value)
    if identity is NO_IDENTITY:
        raise _Unkeyable(f"unstored stage result {_type_name(value)}")
    if identity is not None:
        return ("identity", identity)
    if isinstance(value, GbdrawConfig):
        if config_path is None or configs is None:
            raise _Unkeyable("a GbdrawConfig outside an argument or a frozen dataclass field")
        configs.append((config_path, value))
        return _config_key(value, declared)
    if isinstance(value, SeqRecord):
        raise _Unkeyable("a record without a registered identity")
    if isinstance(value, (tuple, list)):
        return (
            type(value).__name__,
            tuple(_key(item, config_path=None) for item in value),
        )
    if isinstance(value, (set, frozenset)):
        return (
            type(value).__name__,
            tuple(sorted((_key(item, config_path=None) for item in value), key=repr)),
        )
    if isinstance(value, dict):
        return (
            "dict",
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
    path: tuple[Any, ...],
    views: Mapping[tuple[Any, ...], GbdrawConfig],
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
        return True
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


def _argument_ids(args: tuple[object, ...], kwargs: Mapping[str, object]) -> set[int]:
    """The arguments and their direct children: a result may pass these through."""

    ids: set[int] = set()
    for arg in (*args, *kwargs.values()):
        ids.add(id(arg))
        ids.update(id(child) for _name, child in _direct_children(arg))
    return ids


def _register_result(value: object, identity: Hashable, passed: set[int]) -> None:
    """Give a result and its direct fields the stage key, or NO_IDENTITY.

    A value the stage received (an argument or a direct child of one) keeps the
    key it was given as an argument, so its key does not depend on whether a
    stage has returned it.
    """

    if value is None or isinstance(value, (*_PLAIN_TYPES, Enum)) or id(value) in passed:
        return
    register_prepared_value_identity(value, identity)
    for name, child in _direct_children(value):
        if child is None or isinstance(child, (*_PLAIN_TYPES, Enum)) or id(child) in passed:
            continue
        register_prepared_value_identity(
            child,
            NO_IDENTITY if identity is NO_IDENTITY else ("field", identity, name),
        )


def _metric(kind: str, name: str) -> None:
    record_prepared_input_metric(f"renderStage{kind}Count")
    record_prepared_input_metric(f"renderStage{kind}:{name}")


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
    transaction = active_render_stage_transaction()
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
        result = function(*args, **kwargs)
        _register_result(result, NO_IDENTITY, _argument_ids(args, kwargs))
        return result

    enclosing = _FRAMES.get()
    entry = lookup_render_stage(transaction, key)
    if entry is not None:
        _metric("Hit", name)
        for frame in enclosing:
            frame.inner_keys.add(key)
            frame.inner_keys.update(entry.inner_keys)
        result = cast(R, entry.value)
        _register_result(result, key, _argument_ids(args, kwargs))
        return result

    frame = _StageFrame(name)
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
        _FRAMES.reset(token)
    _metric("Build", name)
    if frame.undeclared_sections:
        logger.debug(
            "render stage %s read undeclared sections %s",
            name,
            sorted(frame.undeclared_sections),
        )
        _metric("UndeclaredRead", name)
        _register_result(result, NO_IDENTITY, _argument_ids(args, kwargs))
        return result
    for outer in enclosing:
        outer.inner_keys.add(key)
        outer.inner_keys.update(frame.inner_keys)
    store_render_stage(transaction, key, result, frozenset(frame.inner_keys))
    _register_result(result, key, _argument_ids(args, kwargs))
    return result


__all__ = ["RENDER_STAGE_SETTINGS", "RenderStageName", "render_stage"]
