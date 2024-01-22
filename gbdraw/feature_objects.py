#!/usr/bin/env python
# coding: utf-8
from typing import Tuple


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
        self.product: str = product
        self.gene_biotype: str = gene_biotype


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
