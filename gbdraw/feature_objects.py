#!/usr/bin/env python
# coding: utf-8
from typing import Tuple, List


class FeatureObject:
    def __init__(self, feature_id: str, location: list[Tuple[str, str, str, int, int, bool]], is_directional: bool, color: str, note: str) -> None:
        """
        Represents a general genomic feature.

        Args:
            feature_id (str): Unique identifier for the feature.
            location (list[Tuple[str, str, str, int, int, bool]]): List of tuples representing the feature's location details.
            is_directional (bool): Indicates if the feature has a direction.
            color (str): Color code (e.g., HEX) for the feature's visualization.
            note (str): Additional notes or description of the feature.

        Each tuple in the location list contains information about the feature type, its strand, start and end positions, and a boolean indicating if it's the last part.
        """
        self.feature_id: str = feature_id
        self.location: list[Tuple[str, str, str, int, int, bool]] = location
        self.is_directional: bool = is_directional
        self.color: str = color
        self.note: str = note
    def __str__(self) -> str:
        lines: List[str] = []
        if self.feature_id:
            lines.append(f"ID: {self.feature_id}")
        if self.location:
            lines.append(f"Location: {self.location}")
        if self.is_directional is not None:
            lines.append(f"Is directional: {self.is_directional}")
        if self.color:
            lines.append(f"Color: {self.color}")
        if self.note:
            lines.append(f"Note: {self.note}")
        return "\n".join(lines)

class GeneObject(FeatureObject):
    def __init__(self, feature_id: str, location: list[Tuple[str, str, str, int, int, bool]], is_directional: bool, color: str, note: str, product: str, gene_biotype: str) -> None:
        """
        Represents a gene, inheriting from FeatureObject.

        Args:
            feature_id (str): Unique identifier for the gene.
            location (list[Tuple[str, str, str, int, int, bool]]): List of tuples representing the gene's location details.
            is_directional (bool): Indicates if the gene has a direction.
            color (str): Color code for the gene's visualization.
            note (str): Additional notes or description of the gene.
            product (str): Product of the gene.
            gene_biotype (str): Biotype of the gene.

        This class extends FeatureObject with specific properties for genes, such as product and gene_biotype.
        """
        super().__init__(feature_id, location, is_directional, color, note)
        self.gene_biotype: str = gene_biotype
        self.product: str = product
    def __str__(self) -> str:
        base_str = super().__str__()
        additional_info = []
        if self.gene_biotype:
            additional_info.append(f"Gene Biotype: {self.gene_biotype}")
        if self.product:
            additional_info.append(f"Product: {self.product}")
        return base_str + "\n" + "\n".join(additional_info) if additional_info else base_str

class RepeatObject(FeatureObject):
    def __init__(self, feature_id: str, location: list[Tuple[str, str, str, int, int, bool]], is_directional: bool, color: str, note: str, rpt_family: str, rpt_type: str) -> None:
        """
        Represents a repeat region, inheriting from FeatureObject.

        Args:
            feature_id (str): Unique identifier for the repeat region.
            location (list[Tuple[str, str, str, int, int, bool]]): List of tuples representing the repeat region's location details.
            is_directional (bool): Indicates if the repeat region has a direction.
            color (str): Color code for the repeat region's visualization.
            note (str): Additional notes or description of the repeat region.
            rpt_family (str): Family of the repeat region.
            rpt_type (str): Type of the repeat region.

        This class extends FeatureObject with specific properties for repeat regions, such as rpt_family and rpt_type.
        """
        super().__init__(feature_id, location, is_directional, color, note)
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