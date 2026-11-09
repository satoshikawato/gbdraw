"""A Linear GC content, GC skew or Depth path reaches the end of what it draws (OV-2871)."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.gc import gc_content_percent_df
from gbdraw.analysis.skew import skew_df
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.scalar_axis import ScalarRow, ScalarRows
from gbdraw.svg.linear_tracks import (
    calculate_depth_path_desc,
    calculate_gc_content_path_desc,
    calculate_gc_skew_path_desc,
    calculate_linear_scalar_area_path_desc,
)

# 671 bp sampled every 100 bp: the last window starts at 600, 71 bp short of the end.
RECORD_LENGTH = 671
WIDTH = 1000.0
HEIGHT = 20.0


def _gc_df() -> pd.DataFrame:
    rng = np.random.default_rng(2865)
    sequence = "".join(rng.choice(list("ACGT"), RECORD_LENGTH))
    return skew_df(SeqRecord(Seq(sequence)), window=1000, step=100, nt="GC", circular=True)


def _depth_df() -> pd.DataFrame:
    positions = np.arange(0, RECORD_LENGTH, 100)
    return pd.DataFrame({"position": positions, "depth_normalized": np.linspace(0.2, 0.8, len(positions))})


def _path_points(desc: str) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in re.findall(r"[ML](-?[\d.]+) (-?[\d.]+)", desc)]


def _paths(record_len: int, width: float, row: ScalarRow | None) -> dict[str, str]:
    gc = _gc_df()
    return {
        "gc content": calculate_gc_content_path_desc(0.0, 0.0, gc, record_len, width, HEIGHT, "GC", row=row),
        "gc skew": calculate_gc_skew_path_desc(0.0, 0.0, gc, record_len, width, HEIGHT, row=row),
        "gc percent": calculate_linear_scalar_area_path_desc(
            0.0, 0.0, gc_content_percent_df(gc, dinucleotide="GC"), record_len, width, HEIGHT, row=row
        ),
        "depth": calculate_depth_path_desc(0.0, 0.0, _depth_df(), record_len, width, HEIGHT, row=row),
    }


@pytest.mark.parametrize("track", ["gc content", "gc skew", "gc percent", "depth"])
def test_a_linear_track_line_reaches_the_record_end(track: str) -> None:
    points = _path_points(_paths(RECORD_LENGTH, WIDTH, None)[track])
    line = points[1:-2]  # Between the baseline start and the two baseline points that close the area.
    assert line[-1][0] == pytest.approx(WIDTH, abs=1e-3)
    # The line holds the last window's value from its sample to the end (no value is made up).
    last_sample_x = WIDTH * 600 / RECORD_LENGTH
    held = [y for x, y in line if x >= last_sample_x - 1e-6]
    assert len(set(held)) == 1


@pytest.mark.parametrize("track", ["gc content", "gc skew", "gc percent", "depth"])
def test_every_wrapped_row_of_a_linear_track_reaches_its_row_end(track: str) -> None:
    transform = RecordDisplayTransform(
        length=RECORD_LENGTH, source_base=1, source_step=1, start_coordinate=None, is_circular=False
    ).with_wrap(300, RECORD_LENGTH)
    samples = ScalarRows(transform)
    for row in range(transform.row_count):
        row_length = transform.row_length(row)
        row_width = WIDTH * row_length / transform.axis_length
        points = _path_points(_paths(row_length, row_width, ScalarRow(samples, row))[track])
        assert max(x for x, _ in points) == pytest.approx(row_width, abs=1e-3), (track, row)
