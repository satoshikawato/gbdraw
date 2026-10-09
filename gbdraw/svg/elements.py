"""svgwrite element classes that gbdraw builds without svgwrite's validation.

svgwrite checks every attribute name and value of an element built in its
default debug mode, when the attribute is set and again when the element is
serialized. gbdraw checks each user color and SVG keyword once when it reads
it (``gbdraw.io.colors.check_user_color`` and
``gbdraw.config.models.root.validate_style_leaf``) and builds every other
attribute value itself, so its elements skip that check. The classes keep svgwrite's
names and constructors; only the default of ``debug`` differs.
"""

from __future__ import annotations

from typing import Any

from svgwrite import container, gradients, masking, path, shapes, text


class Group(container.Group):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class Path(path.Path):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class Text(text.Text):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class TSpan(text.TSpan):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class TextPath(text.TextPath):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class Line(shapes.Line):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class Rect(shapes.Rect):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class Circle(shapes.Circle):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class ClipPath(masking.ClipPath):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


class LinearGradient(gradients.LinearGradient):
    def __init__(self, *args: Any, **extra: Any) -> None:
        extra.setdefault("debug", False)
        super().__init__(*args, **extra)


__all__ = [
    "Circle",
    "ClipPath",
    "Group",
    "Line",
    "LinearGradient",
    "Path",
    "Rect",
    "TSpan",
    "Text",
    "TextPath",
]
