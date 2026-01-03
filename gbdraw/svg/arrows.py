#!/usr/bin/env python
# coding: utf-8


def set_arrow_shoulder(feat_strand: str, arrow_end: float, cds_arrow_length: float) -> float:
    """
    Calculates the shoulder position for an arrow representing a genomic feature.
    """
    if feat_strand == "positive":
        shoulder = float(arrow_end - cds_arrow_length)
    else:
        shoulder = float(arrow_end + cds_arrow_length)
    return shoulder


__all__ = ["set_arrow_shoulder"]


