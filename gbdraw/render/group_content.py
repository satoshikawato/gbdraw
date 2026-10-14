"""Group contents: what a group class built, kept by a render stage (design 2.4).

A ``GroupContent`` is a stage result: the attributes and the child elements a
group class produced, and the definitions the group needs in ``defs``. It is
never positioned, composed or changed after its stage. For each request a
builder places it with ``place_content``: the definitions go to the drawing's
``defs`` by today's rules, and a new ``ContentGroup`` over the same child
objects is returned for the builder to position and add. The document writer
then writes the stored text of the children instead of serializing them again.

``GroupKey`` is the identity of a placed group and ``group_manifest`` lists
the placed groups of a drawing in document order.
"""

from __future__ import annotations

import copy
import hashlib
from dataclasses import dataclass, field
from enum import Enum
from functools import cached_property
from types import MappingProxyType
from typing import Callable, Mapping, ParamSpec, Protocol

from svgwrite import Drawing
from svgwrite.base import BaseElement
from svgwrite.container import Group
from svgwrite.etree import etree
from svgwrite.params import Parameter

from gbdraw.layout.spatial import Aabb

P = ParamSpec("P")


class GroupKind(str, Enum):
    """The closed set of group kinds a content can have."""

    # Linear: the record group that draws the axis and the ruler (it carries the record span).
    AXIS = "axis"
    # Circular: the axis circle alone.
    AXIS_CIRCLE = "axis_circle"
    FEATURES = "features"
    LABEL_LEADERS = "label_leaders"
    LABEL_TEXT = "label_text"
    TICKS = "ticks"
    GC_CONTENT = "gc_content"
    GC_SKEW = "gc_skew"
    DEPTH = "depth"
    CONSERVATION = "conservation"
    ANNOTATIONS = "annotations"
    DEFINITION = "definition"
    COMPARISON = "comparison"
    LENGTH_BAR = "length_bar"
    LEGEND = "legend"
    PLOT_TITLE = "plot_title"
    RECORD = "record"


CANVAS_SCOPE = "canvas"


@dataclass(frozen=True)
class GroupKey:
    """The identity of a placed group.

    ``scope`` is ``"canvas"``, a record key, or a pair of record keys for a
    comparison. ``slot`` is the track slot id, or the ordinal among the groups
    of ``kind`` in the scope when the request has no slot ids.
    """

    scope: str | tuple[str, str]
    kind: GroupKind
    slot: str | int = 0


def element_text(element: BaseElement) -> str:
    """svgwrite's serialization of one element, as ``Drawing.tostring`` writes it."""

    return etree.tostring(element.get_xml(), encoding="utf-8").decode("utf-8")


def start_tag(element: BaseElement) -> str:
    """The start tag of ``element`` by svgwrite's attribute rules (a child-less copy)."""

    bare = _childless(element)
    text = element_text(bare)
    if not text.endswith(" />"):  # pragma: no cover - ElementTree writes empty elements so
        raise ValueError(f"unexpected empty-element text for <{element.elementname}>")
    return text[:-3] + ">"


def _childless(element: BaseElement) -> BaseElement:
    bare = object.__new__(type(element))
    bare.__dict__.update(element.__dict__)
    bare.elements = []
    return bare


@dataclass(frozen=True)
class DefinitionMeasurement:
    """What a record definition measured besides its group (design 2.4, M6)."""

    width: float
    height: float
    text_anchor: str
    line_count: int


@dataclass(frozen=True, eq=False)
class GroupContent:
    """The attributes, children and definitions one group class built.

    ``local_bounds`` is the extent the group class measured around its own
    origin, when the caller places or composes by it; ``definition`` is the
    measurement of a record definition, which the caller places by.
    """

    attributes: Mapping[str, object]
    children: tuple[BaseElement, ...]
    definitions: tuple[BaseElement, ...] = ()
    local_bounds: Aabb | None = None
    definition: DefinitionMeasurement | None = None
    # svgwrite's per-element parameter (debug flag, profile) of the built group.
    parameter: Parameter | None = field(default=None, repr=False)

    @classmethod
    def of(
        cls,
        group: Group,
        definitions: tuple[BaseElement, ...] = (),
        *,
        local_bounds: Aabb | None = None,
        definition: DefinitionMeasurement | None = None,
    ) -> GroupContent:
        """The content of a group a group class has just built."""

        if group.elementname != "g":
            raise TypeError(f"a group content needs a <g>, not <{group.elementname}>")
        return cls(
            attributes=MappingProxyType(dict(group.attribs)),
            children=tuple(group.elements),
            definitions=tuple(definitions),
            local_bounds=local_bounds,
            definition=definition,
            parameter=getattr(group, "_parameter", None),
        )

    @cached_property
    def children_text(self) -> str:
        """The serialized children, computed when first asked for."""

        return "".join(element_text(child) for child in self.children)

    @cached_property
    def text(self) -> str:
        """The group as its group class built it: its own attributes and children, unplaced."""

        group = Group(debug=False)
        if self.parameter is not None:
            group._parameter = self.parameter
        group.attribs = dict(self.attributes)
        return start_tag(group) + self.children_text + "</g>"

    @cached_property
    def revision(self) -> str:
        """A digest of the serialized children."""

        return hashlib.sha256(self.children_text.encode("utf-8")).hexdigest()

    @property
    def definition_ids(self) -> tuple[str, ...]:
        return tuple(str(item.attribs.get("id", "")) for item in self.definitions)


class ContentGroup(Group):
    """A new group of one request over a stored content's child objects."""

    def __init__(self, content: GroupContent, key: GroupKey) -> None:
        super().__init__(debug=False)
        if content.parameter is not None:
            self._parameter = content.parameter
        self.attribs = dict(content.attributes)
        self.elements = list(content.children)
        self.content = content
        self.key = key

    def __deepcopy__(self, memo: dict[int, object]) -> Group:
        """A deep copy is a plain group: its children are copies, not the content's."""

        copied = Group.__new__(Group)
        memo[id(self)] = copied
        for name, value in self.__dict__.items():
            if name not in ("content", "key"):
                copied.__dict__[name] = copy.deepcopy(value, memo)
        return copied

    def holds_content_children(self) -> bool:
        """Whether the children are still exactly the content's children."""

        children = self.content.children
        return len(self.elements) == len(children) and all(
            element is child for element, child in zip(self.elements, children)
        )


class BuildsGroup(Protocol):
    """A group class instance: it built its group in ``__init__``."""

    def get_group(self) -> Group: ...


def build_group_content(
    group_class: Callable[P, BuildsGroup],
    /,
    *args: P.args,
    **kwargs: P.kwargs,
) -> GroupContent:
    """A render-stage function: build one group with ``group_class`` and keep its content.

    The group class is keyed by its qualified name, its arguments as the
    stage's arguments. A ``local_bounds`` or a ``definition_measurement`` the
    instance measured is kept.
    """

    built = group_class(*args, **kwargs)
    bounds = getattr(built, "local_bounds", None)
    definition = getattr(built, "definition_measurement", None)
    return GroupContent.of(
        built.get_group(),
        local_bounds=bounds if isinstance(bounds, Aabb) else None,
        definition=definition if isinstance(definition, DefinitionMeasurement) else None,
    )


def add_content_definitions(drawing: Drawing, content: GroupContent) -> None:
    """Add the definitions of ``content`` to ``drawing.defs`` by today's rules.

    A pattern (hatch) is added unless an element with its id is already there;
    any other definition (a clip path) is always added.
    """

    present = {
        str(getattr(item, "attribs", {}).get("id", ""))
        for item in getattr(drawing.defs, "elements", ())
    }
    for definition in content.definitions:
        definition_id = str(definition.attribs.get("id", ""))
        if definition.elementname == "pattern" and definition_id in present:
            continue
        drawing.defs.add(definition)
        present.add(definition_id)


def place_content(drawing: Drawing, content: GroupContent, key: GroupKey) -> ContentGroup:
    """Add the definitions of ``content`` and return a new group over its children."""

    add_content_definitions(drawing, content)
    return ContentGroup(content, key)


@dataclass(frozen=True)
class RenderedGroup:
    """One placed group of a drawing: identity, DOM id, position and content."""

    key: GroupKey
    dom_id: str | None
    parent: GroupKey | None
    position: str
    content: GroupContent
    definition_ids: tuple[str, ...]


def group_manifest(drawing: Drawing) -> tuple[RenderedGroup, ...]:
    """The placed top-level groups of ``drawing`` in document order."""

    return tuple(
        RenderedGroup(
            key=element.key,
            dom_id=(
                str(element.attribs["id"]) if element.attribs.get("id") is not None else None
            ),
            parent=None,
            position=str(element.attribs.get("transform", "")),
            content=element.content,
            definition_ids=element.content.definition_ids,
        )
        for element in drawing.elements
        if isinstance(element, ContentGroup)
    )


__all__ = [
    "BuildsGroup",
    "CANVAS_SCOPE",
    "ContentGroup",
    "DefinitionMeasurement",
    "GroupContent",
    "GroupKey",
    "GroupKind",
    "RenderedGroup",
    "add_content_definitions",
    "build_group_content",
    "element_text",
    "group_manifest",
    "place_content",
    "start_tag",
]
