#!/usr/bin/env python
# coding: utf-8

"""
Core feature data models.

`FeatureObject.location` is a list of 6-tuples:
  (kind, segment_id, strand, start, end, is_last)

- kind: "block" (exon / feature block) or "line" (intron / connector)
- segment_id: zero-padded id string (e.g. "001")
- strand: "positive" / "negative" / "undefined"
- start/end: 0/1-based-ish coordinates depending on upstream parsing (kept as-is for compatibility)
- is_last: True only for the last exon block in a multi-part feature
"""

from __future__ import annotations

from typing import Any, List, Literal, NamedTuple, Optional, cast

from .shapes import FeatureGlyph


FeatureSegmentKind = Literal["block", "line"]
Strand = Literal["positive", "negative", "undefined"]


def _resolve_glyph_kind(
    *,
    is_directional: bool | None,
    glyph_kind: FeatureGlyph | None,
) -> FeatureGlyph:
    if is_directional is not None and not isinstance(is_directional, bool):
        raise ValueError("is_directional must be a Boolean or None")

    if glyph_kind is None:
        if is_directional is None:
            raise ValueError("either glyph_kind or is_directional must be provided")
        return "arrow" if is_directional else "rectangle"

    normalized = str(glyph_kind).strip().lower()
    if normalized not in {"arrow", "arrowhead", "rectangle"}:
        raise ValueError(
            "glyph_kind must be 'arrow', 'arrowhead', or 'rectangle'"
        )
    resolved = cast(FeatureGlyph, normalized)
    if (
        is_directional is not None
        and is_directional != (resolved in {"arrow", "arrowhead"})
    ):
        raise ValueError("glyph_kind contradicts is_directional")
    return resolved


class FeatureLocationPart(NamedTuple):
    """
    A single location segment.

    This is intentionally a NamedTuple so it supports BOTH:
      - tuple-style access: `part[0]`, `part[3]`, ...
      - attribute access: `part.kind`, `part.start`, ...
    """

    kind: FeatureSegmentKind
    segment_id: str
    strand: Strand
    start: int
    end: int
    is_last: bool
FeatureLocation = list[FeatureLocationPart]


class FeatureObject:
    def __init__(
        self,
        feature_id: str,
        location: FeatureLocation,
        is_directional: bool | None,
        color: str,
        note: str,
        label_text: str,
        coordinates: Any,
        type: str,
        qualifiers: dict,
        record_id: Optional[str] = None,
        glyph_kind: FeatureGlyph | None = None,
    ) -> None:
        """
        Represents a general genomic feature.
        """
        self.feature_id: str = feature_id
        # Accept both raw 6-tuples and FeatureLocationPart instances for compatibility.
        self.location: FeatureLocation = [
            part if isinstance(part, FeatureLocationPart) else FeatureLocationPart(*part)  # type: ignore[arg-type]
            for part in (location or [])
        ]
        self.glyph_kind: FeatureGlyph = _resolve_glyph_kind(
            is_directional=is_directional,
            glyph_kind=glyph_kind,
        )
        self.coordinates = coordinates
        self.color: str = color
        self.note: str = note
        self.label_text: str = label_text
        self.qualifiers: dict = qualifiers
        self.record_id: Optional[str] = record_id
        self.feature_track_id: int = 0
        self.source_feature_index: int | None = None
        self._feature_type: str = type
        # Frequently used derived attribute (historically added dynamically elsewhere)
        self.strand: Strand = self.location[0].strand if self.location else "undefined"

    @property
    def feature_type(self) -> str:
        """Preferred name for the feature kind (e.g. CDS, tRNA, repeat_region)."""
        return self._feature_type

    @feature_type.setter
    def feature_type(self, value: str) -> None:
        self._feature_type = value

    @property
    def is_directional(self) -> bool:
        """Backward-compatible directionality derived from the stored glyph."""

        return self.glyph_kind in {"arrow", "arrowhead"}

    @property
    def type(self) -> str:
        """Backward-compatible alias for `feature_type`."""
        return self._feature_type

    @type.setter
    def type(self, value: str) -> None:
        self._feature_type = value

    def __str__(self) -> str:
        lines: List[str] = []
        if self.feature_id:
            lines.append(f"ID: {self.feature_id}")
        if self.feature_type:
            lines.append(f"Type: {self.feature_type}")
        if self.location:
            lines.append(f"Location: {self.location}")
        if self.coordinates:
            lines.append(f"Coordinates: {self.coordinates}")
        if self.is_directional is not None:
            lines.append(f"Is directional: {self.is_directional}")
        if self.color:
            lines.append(f"Color: {self.color}")
        if self.note:
            lines.append(f"Note: {self.note}")
        if self.label_text:
            lines.append(f"Label text: {self.label_text}")
        return "\n".join(lines)

class GeneObject(FeatureObject):
    def __init__(
        self,
        feature_id: str,
        location: FeatureLocation,
        is_directional: bool | None,
        color: str,
        note: str,
        product: str,
        gene_biotype: str,
        gene: str,
        label_text: str,
        coordinates: Any,
        type: str,
        qualifiers: dict,
        record_id: Optional[str] = None,
        glyph_kind: FeatureGlyph | None = None,
    ) -> None:
        """
        Represents a gene, inheriting from FeatureObject.
        """
        super().__init__(
            feature_id,
            location,
            is_directional,
            color,
            note,
            label_text,
            coordinates,
            type,
            qualifiers,
            record_id=record_id,
            glyph_kind=glyph_kind,
        )
        self.gene_biotype: str = gene_biotype
        self.product: str = product
        self.gene: str = gene

    def __str__(self) -> str:
        base_str = super().__str__()
        additional_info = []
        if self.gene_biotype:
            additional_info.append(f"Gene Biotype: {self.gene_biotype}")
        if self.gene:
            additional_info.append(f"Gene: {self.gene}")
        if self.product:
            additional_info.append(f"Product: {self.product}")
        return base_str + "\n" + "\n".join(additional_info) if additional_info else base_str

class RepeatObject(FeatureObject):
    def __init__(
        self,
        feature_id: str,
        location: FeatureLocation,
        is_directional: bool | None,
        color: str,
        note: str,
        rpt_family: str,
        rpt_type: str,
        label_text: str,
        coordinates: Any,
        type: str,
        qualifiers: dict,
        record_id: Optional[str] = None,
        glyph_kind: FeatureGlyph | None = None,
    ) -> None:
        """
        Represents a repeat region, inheriting from FeatureObject.
        """
        super().__init__(
            feature_id,
            location,
            is_directional,
            color,
            note,
            label_text,
            coordinates,
            type,
            qualifiers,
            record_id=record_id,
            glyph_kind=glyph_kind,
        )
        self.rpt_family: str = rpt_family
        self.rpt_type: str = rpt_type

    def __str__(self) -> str:
        base_str = super().__str__()
        additional_info = []
        if self.rpt_family:
            additional_info.append(f"Repeat Family: {self.rpt_family}")
        if self.rpt_type:
            additional_info.append(f"Repeat Type: {self.rpt_type}")
        return base_str + "\n" + "\n".join(additional_info) if additional_info else base_str


__all__ = [
    "FeatureLocation",
    "FeatureLocationPart",
    "FeatureObject",
    "FeatureSegmentKind",
    "GeneObject",
    "RepeatObject",
    "Strand",
]
