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
from typing import Callable, Mapping, ParamSpec, Protocol, cast
from xml.etree.ElementTree import Element

from svgwrite import Drawing
from svgwrite.base import BaseElement
from svgwrite.container import Group
from svgwrite.etree import CDATA_TAG, etree
from svgwrite.params import Parameter

from gbdraw.exceptions import ExportError
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
    """svgwrite's serialization of one element, as ``Drawing.tostring`` writes it.

    An empty element (most of a drawing: paths, lines, circles) is written by
    ``empty_element_markup``; any other through ElementTree.
    """

    xml = element.get_xml()
    markup = empty_element_markup(xml)
    if markup is not None:
        return markup
    return etree.tostring(xml, encoding="utf-8").decode("utf-8")


def _attribute_value(value: str) -> str:
    # ElementTree's attribute escaping (``_escape_attrib``), in its order.
    if "&" in value:
        value = value.replace("&", "&amp;")
    if "<" in value:
        value = value.replace("<", "&lt;")
    if ">" in value:
        value = value.replace(">", "&gt;")
    if '"' in value:
        value = value.replace('"', "&quot;")
    if "\r" in value:
        value = value.replace("\r", "&#13;")
    if "\n" in value:
        value = value.replace("\n", "&#10;")
    if "\t" in value:
        value = value.replace("\t", "&#09;")
    return value


def empty_element_markup(xml: Element) -> str | None:
    """``etree.tostring(xml, encoding="utf-8")`` as text, for an element with
    no children, text or tail, plain string attributes and a plain tag; None
    for any other element (the caller serializes it with ElementTree).

    ElementTree writes the attributes in their order, escapes the values with
    ``_escape_attrib`` and, for UTF-8, writes a character UTF-8 cannot encode
    (a lone surrogate) as a character reference.
    """

    tag = xml.tag
    if (
        len(xml)
        or xml.text
        or xml.tail
        or not isinstance(tag, str)
        or tag == CDATA_TAG
        or tag[:1] == "{"
    ):
        return None
    parts = ["<", tag]
    for name, value in xml.items():
        if not isinstance(name, str) or not isinstance(value, str) or name[:1] == "{":
            return None
        parts.append(f' {name}="{_attribute_value(value)}"')
    parts.append(" />")
    markup = "".join(parts)
    if markup.isascii():
        return markup
    return markup.encode("utf-8", "xmlcharrefreplace").decode("utf-8")


_NOT_AN_ELEMENT = (
    "An SVG document piece holds markup that is not an element "
    "(a comment, CDATA or a processing instruction)."
)


def require_element_markup(text: str) -> str:
    """Return ``text`` if it holds no ``<!`` or ``<?`` (precondition P, gbdraw/web/CLAUDE.md).

    svgwrite escapes ``<`` in text and attribute values, so ``<!`` and ``<?``
    appear only as a comment, CDATA or a processing instruction, which the
    page could not sanitize piece by piece as it sanitizes the whole document.
    """

    if "<!" in text or "<?" in text:
        raise ExportError(_NOT_AN_ELEMENT)
    return text


def start_tag(element: BaseElement) -> str:
    """The start tag of ``element`` by svgwrite's attribute rules (a child-less copy)."""

    bare = _childless(element)
    text = element_text(bare)
    if not text.endswith(" />"):  # pragma: no cover - ElementTree writes empty elements so
        raise ValueError(f"unexpected empty-element text for <{element.elementname}>")
    return text[:-3] + ">"


# Values an element keeps that nothing changes in place: shared by a copy.
_SHARED_VALUES = (str, int, float, bool, bytes, type(None), Parameter)


def copy_element_tree(element: BaseElement) -> BaseElement:
    """A copy of ``element`` and its subtree that shares nothing a later step changes.

    Every element is a new object with a new ``attribs`` dict, a new
    ``elements`` list and new copies of the other lists, dicts and objects it
    keeps (``Path.commands``; a ``TextPath``'s ``href`` element). As with
    ``copy.deepcopy`` of the element, an element referenced twice in the copy
    (a ``TextPath``'s path that is also a child) is copied once, and an element
    referenced from outside the subtree is copied apart from it. Strings,
    numbers and svgwrite's per-element ``Parameter`` (debug flag and profile,
    which elements a ``Drawing`` creates share and nothing changes) are
    shared. A value of any other type is copied with ``copy.deepcopy``.
    """

    return _copy_element(element, {})


def _copy_element(element: BaseElement, memo: dict[int, object]) -> BaseElement:
    found = memo.get(id(element))
    if found is not None:
        return cast(BaseElement, found)
    duplicate = object.__new__(type(element))
    memo[id(element)] = duplicate
    state = duplicate.__dict__
    for name, value in element.__dict__.items():
        state[name] = _copy_value(value, memo)
    return duplicate


def _copy_value(value: object, memo: dict[int, object]) -> object:
    if isinstance(value, _SHARED_VALUES):
        return value
    if isinstance(value, BaseElement):
        return _copy_element(value, memo)
    if type(value) is list:
        return [_copy_value(item, memo) for item in value]
    if type(value) is dict:
        return {key: _copy_value(item, memo) for key, item in value.items()}
    return copy.deepcopy(value, memo)


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
        """The serialized children, computed when first asked for.

        Each child is one whole element and the text holds no ``<!`` or ``<?``
        (precondition P, ``require_element_markup``); otherwise this raises.
        """

        texts = []
        for child in self.children:
            text = element_text(child)
            name = child.elementname
            if not text.startswith(f"<{name}") or not (
                text.endswith("/>") or text.endswith(f"</{name}>")
            ):
                raise ExportError(_NOT_AN_ELEMENT)
            texts.append(text)
        return require_element_markup("".join(texts))

    @cached_property
    def text(self) -> str:
        """The group as its group class built it: its own attributes and children, unplaced."""

        group = Group(debug=False)
        if self.parameter is not None:
            group._parameter = self.parameter
        group.attribs = dict(self.attributes)
        return start_tag(group) + self.children_text + "</g>"

    @cached_property
    def _children_digest(self) -> tuple[str, int]:
        data = self.children_text.encode("utf-8")
        return hashlib.sha256(data).hexdigest(), len(data)

    @property
    def revision(self) -> str:
        """The sha256 of the serialized children's UTF-8 bytes."""

        return self._children_digest[0]

    @property
    def children_byte_length(self) -> int:
        """The UTF-8 length of the serialized children."""

        return self._children_digest[1]

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
    "copy_element_tree",
    "element_text",
    "empty_element_markup",
    "group_manifest",
    "place_content",
    "require_element_markup",
    "start_tag",
]
