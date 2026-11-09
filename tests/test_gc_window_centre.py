"""Every GC content and GC skew track draws a window's value at the window's centre (OV-2870).

The record is all A except a run of G bases. A window centred on a sample counts the run when the sample lies
within half a window of it, so the samples that see the run are centred on the run, in both modes, for the
content and the skew track, and at the record's ends (the circular windows wrap) and in a crop.
"""

from __future__ import annotations

import math
import re
from xml.etree import ElementTree as ET

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.skew import skew_df
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.io.regions import RegionSpec

LENGTH = 10_000
WINDOW = 1_000
STEP = 100
_POINT = re.compile(r"[ML](-?[\d.e+-]+)[ ,](-?[\d.e+-]+)")


def _record(run_start: int, run_end: int) -> SeqRecord:
    bases = ["A"] * LENGTH
    bases[run_start:run_end] = ["G"] * (run_end - run_start)
    return SeqRecord(Seq("".join(bases)), id="run", annotations={"topology": "circular", "molecule_type": "DNA"})


def _circular_mean(positions: list[float], length: int) -> float:
    angles = [2.0 * math.pi * p / length for p in positions]
    mean = math.atan2(sum(math.sin(a) for a in angles), sum(math.cos(a) for a in angles))
    return (mean / (2.0 * math.pi) * length) % length


def _circular_distance(a: float, b: float, length: int) -> float:
    return abs((a - b + length / 2.0) % length - length / 2.0)


@pytest.mark.parametrize("column", ["GC content", "GC skew"])
@pytest.mark.parametrize("run", [(4000, 4100), (0, 100), (9900, 10_000)], ids=["middle", "record-start", "record-end"])
def test_the_samples_that_see_a_run_are_centred_on_it(column: str, run: tuple[int, int]) -> None:
    frame = skew_df(_record(*run), WINDOW, STEP, "GC", circular=True)
    seen = [float(p) for p in frame.index[frame[column].to_numpy() > 0]]
    # The window centred on p, [p - 500, p + 500), meets the run [s, e) when s - 500 < p < e + 500; a circular
    # record's windows wrap, so p is read on the circle.
    expected = sorted({p % LENGTH for p in range(-WINDOW, LENGTH + WINDOW, STEP)
                       if run[0] - WINDOW // 2 < p < run[1] + WINDOW // 2})
    assert sorted(int(p) for p in seen) == expected
    assert _circular_distance(_circular_mean(seen, LENGTH), sum(run) / 2.0, LENGTH) < 1.0


def test_a_window_value_counts_the_bases_centred_on_its_sample() -> None:
    record = _record(4000, 4100)
    frame = skew_df(record, WINDOW, STEP, "GC", circular=True)
    sequence = str(record.seq)
    for position in (0, 3500, 4000, 4500, 9900):
        start = position - WINDOW // 2
        window = "".join(sequence[(start + i) % LENGTH] for i in range(WINDOW))
        assert frame.loc[position, "GC content"] == pytest.approx(window.count("G") / WINDOW)


def _track_paths(mode: str, record: SeqRecord, region: RegionSpec | None) -> dict[str, str]:
    config = apply_config_overrides(None, {"canvas.show_gc": True, "canvas.show_skew": True})
    request_type, options_type = (
        (CircularDiagramRequest, CircularDiagramOptions) if mode == "circular"
        else (LinearDiagramRequest, LinearDiagramOptions)
    )
    request = request_type(
        records=(RecordInput(InMemoryRecordSource(record), region=region),),
        options=options_type(config=config, window=WINDOW, step=STEP),
    )
    result = plan_request(request).build()
    root = ET.fromstring((result.drawing if hasattr(result, "drawing") else result).tostring())
    paths: dict[str, str] = {}
    for group in root.iter():
        # Linear tags its slot groups; Circular's default groups keep their preset ids.
        track = {"dinucleotide_content": "GC content", "dinucleotide_skew": "GC skew"}.get(
            group.get("data-gbdraw-slot-renderer") or ""
        ) or {"gc_content": "GC content", "skew": "GC skew"}.get(group.get("id") or "")
        if track is None:
            continue
        # The data line is the longest path of the group (a clip or a background rectangle is a few points).
        descs = [node.get("d", "") for node in group.iter() if node.tag.endswith("path")]
        paths[track] = max(descs, key=lambda d: len(_POINT.findall(d)))
    return paths


def _seen_positions(mode: str, desc: str, length: int) -> list[float]:
    """The genome positions of the drawn points that stand out from the track (the samples that see the run)."""

    points = np.array([(float(x), float(y)) for x, y in _POINT.findall(desc)])
    if mode == "circular":
        heights = np.hypot(points[:, 0], points[:, 1])
        positions = (np.arctan2(points[:, 0], -points[:, 1]) % (2.0 * math.pi)) / (2.0 * math.pi) * length
    else:
        line = points[1:]  # After the move to the baseline start; the line runs from base 0 to the drawn end.
        heights = -line[:, 1]
        positions = (line[:, 0] - line[0, 0]) / (line[:, 0].max() - line[0, 0]) * length
    background = float(np.median(heights))
    peak = float(heights.max())
    return [float(p) for p in positions[heights > background + 0.5 * (peak - background)]]


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize(
    "run,region,length,centre",
    [
        ((4000, 4100), None, LENGTH, 4050.0),
        ((0, 100), None, LENGTH, 50.0),
        # The crop 3001..6000 draws bases 3000..5999; the run sits at 1000..1100 of it.
        ((4000, 4100), RegionSpec("run:3001-6000", None, "run", None, 3001, 6000, False), 3000, 1050.0),
    ],
    ids=["middle", "record-start", "crop"],
)
def test_both_modes_draw_content_and_skew_at_the_window_centre(
    mode: str, run: tuple[int, int], region: RegionSpec | None, length: int, centre: float
) -> None:
    paths = _track_paths(mode, _record(*run), region)
    assert set(paths) == {"GC content", "GC skew"}
    for track, desc in paths.items():
        seen = _seen_positions(mode, desc, length)
        assert seen, track
        # Start-anchored windows put the run half a window (500 bp) before its place.
        assert _circular_distance(_circular_mean(seen, length), centre, length) < STEP / 2, (mode, track)
