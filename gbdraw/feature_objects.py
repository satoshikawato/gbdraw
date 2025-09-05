#!/usr/bin/env python
# coding: utf-8
from typing import Tuple, List


class FeatureObject:
    def __init__(self, feature_id: str, location: list[Tuple[str, str, str, int, int, bool]], is_directional: bool, color: str, note: str, label_text: str, coordinates, type, qualifiers: dict) -> None:
        """
        Represents a general genomic feature.
        """
        self.feature_id: str = feature_id
        self.location: list[Tuple[str, str, str, int, int, bool]] = location
        self.is_directional: bool = is_directional
        self.coordinates = coordinates
        self.color: str = color
        self.note: str = note
        self.label_text: str = label_text
        self.qualifiers: dict = qualifiers
        self.feature_track_id = 0
        self.type = type

    def __str__(self) -> str:
        lines: List[str] = []
        if self.feature_id:
            lines.append(f"ID: {self.feature_id}")
        if self.type:
            lines.append(f"Type: {self.type}")
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
    def __init__(self, feature_id: str, location: list[Tuple[str, str, str, int, int, bool]], is_directional: bool, color: str, note: str, product: str, gene_biotype: str, gene: str, label_text: str, coordinates, type, qualifiers: dict) -> None:
        """
        Represents a gene, inheriting from FeatureObject.
        """
        super().__init__(feature_id, location, is_directional, color, note, label_text, coordinates, type, qualifiers)
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
    def __init__(self, feature_id: str, location: list[Tuple[str, str, str, int, int, bool]], is_directional: bool, color: str, note: str, rpt_family: str, rpt_type: str, label_text: str, coordinates, type, qualifiers: dict) -> None:
        """
        Represents a repeat region, inheriting from FeatureObject.
        """
        super().__init__(feature_id, location, is_directional, color, note, label_text, coordinates, type, qualifiers)
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