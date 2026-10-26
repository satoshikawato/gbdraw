#!/usr/bin/env python
# coding: utf-8

from typing import overload

import numpy as np
from numpy.typing import NDArray


@overload
def normalize_position_to_linear_track(
    position: float, genome_length: float, alignment_width: float
) -> float: ...


@overload
def normalize_position_to_linear_track(
    position: NDArray[np.float64], genome_length: float, alignment_width: float
) -> NDArray[np.float64]: ...


def normalize_position_to_linear_track(
    position: float | NDArray[np.float64],
    genome_length: float,
    alignment_width: float,
) -> float | NDArray[np.float64]:
    """A position (or an array of them) as an x on a Linear track of ``alignment_width``."""

    return alignment_width * (position / genome_length)


__all__ = ["normalize_position_to_linear_track"]


