#!/usr/bin/env python
# coding: utf-8

def calculate_gc_percent(sequence: str) -> float:
    """
    Calculate the GC content of a DNA sequence.
    """
    g_count: int = sequence.count("G")
    c_count: int = sequence.count("C")
    gc_percent: float = round(100 * ((g_count + c_count) / len(sequence)), 2)
    return gc_percent


__all__ = ["calculate_gc_percent"]


