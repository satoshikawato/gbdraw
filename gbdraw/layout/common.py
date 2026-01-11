#!/usr/bin/env python
# coding: utf-8

def calculate_cds_ratio(track_ratio, length_param, track_ratio_factor):
    if length_param == "short":
        cds_ratio = float(track_ratio * track_ratio_factor)
        offset = float(0.01)
    else:
        cds_ratio = float(track_ratio * track_ratio_factor)
        offset = float(0.005)
    return cds_ratio, offset


__all__ = ["calculate_cds_ratio"]


