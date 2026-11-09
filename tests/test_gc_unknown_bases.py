"""GC content and GC skew count only the known bases A, C, G and T (OV-2799).

N and the other IUPAC codes are unknown bases. A window's content is the share of G and C among its known bases and
its skew their difference over their sum; a window with no known base has no value, so its tracks break there (they
lie on the baseline), the cumulative skew adds nothing for it, and the deviation baseline and scale read only the
windows with a value. A record of unknown bases only draws no GC track.
"""

from __future__ import annotations

import math
import random
import re
import warnings
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.gc import calculate_gc_percent, gc_content_percent_df
from gbdraw.analysis.skew import skew_df
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.svg.circular_tracks import (
    circular_gc_content_line,
    circular_gc_skew_line,
    generate_circular_gc_content_path_desc,
    generate_circular_gc_skew_path_desc,
    generate_circular_scalar_area_path_desc,
)
from gbdraw.svg.linear_tracks import calculate_gc_content_path_desc, calculate_gc_skew_path_desc

WINDOW = 1_000
STEP = 100
# Bases 3,000 .. 5,999 are unknown: the windows centred on 3,500 .. 5,500 read no known base.
GAP = (3_000, 6_000)
_POINT = re.compile(r"[ML](-?[\d.e+-]+)[ ,](-?[\d.e+-]+)")
_ARC_RADIUS = re.compile(r"A ([\d.]+),")


def _known(length: int, seed: int) -> str:
    rng = random.Random(seed)
    # GC-rich and AT-rich stretches, so the content and the skew vary along the record.
    return "".join(rng.choice("GGGCCA" if (i // 700) % 2 else "AATTGC") for i in range(length))


def _gapped_sequence(unknown: str = "N") -> str:
    return _known(GAP[0], 1) + unknown * (GAP[1] - GAP[0]) + _known(10_000 - GAP[1], 2)


def _record(sequence: str, topology: str = "circular") -> SeqRecord:
    return SeqRecord(Seq(sequence), id="gap", annotations={"topology": topology, "molecule_type": "DNA"})


def _expected(sequence: str, *, circular: bool) -> pd.DataFrame:
    """Content and skew of each sample's window, counted base by base over the known bases."""

    rows = {}
    length = len(sequence)
    for p in range(0, length, STEP):
        indices = range(p - WINDOW // 2, p - WINDOW // 2 + WINDOW)
        bases = "".join(sequence[i % length] for i in indices if circular or 0 <= i < length).upper()
        g, c = bases.count("G"), bases.count("C")
        known = sum(bases.count(base) for base in "ACGTU")
        content = math.nan if known == 0 else (g + c) / known
        skew = math.nan if known == 0 else 0.0 if g + c == 0 else (g - c) / (g + c)
        rows[p] = (content, skew)
    return pd.DataFrame.from_dict(rows, orient="index", columns=["GC content", "GC skew"])


@pytest.mark.parametrize("topology", ["circular", "linear"])
@pytest.mark.parametrize("unknown", ["N", "n", "RYKMSWBDHV", "-"])
def test_a_window_counts_only_its_known_bases(topology: str, unknown: str) -> None:
    sequence = _known(GAP[0], 1) + (unknown * 3_000)[: GAP[1] - GAP[0]] + _known(10_000 - GAP[1], 2).lower()
    frame = skew_df(_record(sequence, topology), WINDOW, STEP, "GC", circular=topology == "circular")
    expected = _expected(sequence, circular=topology == "circular")

    pd.testing.assert_frame_equal(frame[["GC content", "GC skew"]], expected, check_exact=False, rtol=1e-12)
    no_value = frame.index[frame["GC content"].isna()]
    assert list(no_value) == list(range(3_500, 5_501, STEP))
    assert frame.loc[no_value].isna().all().all()


def test_a_window_of_known_bases_without_g_or_c_has_skew_zero() -> None:
    frame = skew_df(_record("AT" * 2_000 + "N" * 2_000), WINDOW, STEP, "GC", circular=False)

    assert frame.loc[1_000, "GC content"] == 0.0
    assert frame.loc[1_000, "GC skew"] == 0.0
    assert math.isnan(frame.loc[5_000, "GC skew"])


def test_the_cumulative_skew_adds_nothing_for_a_window_without_value() -> None:
    frame = skew_df(_record(_gapped_sequence()), WINDOW, STEP, "GC", circular=True)
    valued = frame["GC skew"].notna()
    cumulative = frame["GC skew"][valued].cumsum()
    factor = frame["GC skew"][valued].abs().max() / cumulative.abs().max()

    np.testing.assert_allclose(frame["Cumulative GC skew, normalized"][valued], cumulative * factor, rtol=1e-12)
    assert frame["Cumulative GC skew, normalized"][~valued].isna().all()


def _transform(length: int) -> RecordDisplayTransform:
    return RecordDisplayTransform(length, 1, 1, None, True)


@pytest.mark.parametrize(
    ("line_of", "column"),
    [(circular_gc_content_line, "GC content"), (circular_gc_skew_line, "GC skew")],
)
def test_the_circular_deviation_scale_reads_only_windows_with_a_value(line_of, column: str) -> None:
    # The windows without a value are set here, so the drawing is checked apart from the counting.
    frame = skew_df(_record(_known(10_000, 3)), WINDOW, STEP, "GC", circular=True)
    frame.loc[3_500:5_500] = np.nan
    valued = frame[column].notna()
    if line_of is circular_gc_content_line:
        def line(df):
            return line_of(100.0, 10_000, df, 20.0, 1.0, "GC", _transform(10_000))
    else:
        def line(df):
            return line_of(100.0, df, 10_000, 20.0, 1.0, "GC", _transform(10_000))

    gapped = line(frame)
    alone = line(frame[valued])
    radii = dict(zip(gapped.positions, gapped.radii))
    # The windows with a value draw as if the windows without one were not there; those lie on the baseline.
    np.testing.assert_allclose([radii[p] for p in alone.positions], alone.radii, rtol=1e-12)
    assert all(radii[float(p)] == 100.0 for p in frame.index[~valued])
    assert max(abs(r - 100.0) for r in gapped.radii) == pytest.approx(10.0)


def _line_points(desc: str) -> np.ndarray:
    return np.array([(float(x), float(y)) for x, y in _POINT.findall(desc)])


@pytest.mark.parametrize("track", ["content", "skew"])
def test_the_linear_deviation_path_lies_on_its_baseline_where_windows_have_no_value(track: str) -> None:
    frame = skew_df(_record(_gapped_sequence(), "linear"), WINDOW, STEP, "GC", circular=False)
    if track == "content":
        desc = calculate_gc_content_path_desc(0.0, 0.0, frame, 10_000, 1_000.0, 100.0, "GC")
    else:
        desc = calculate_gc_skew_path_desc(0.0, 0.0, frame, 10_000, 1_000.0, 100.0)
    points = _line_points(desc)[1:]  # After the move to the baseline start.

    # x is a tenth of the position: the windows centred on 3,500 .. 5,500 have no value.
    inside = (points[:, 0] >= 350.0) & (points[:, 0] <= 550.0)
    assert inside.any()
    assert np.all(points[inside, 1] == 0.0)
    assert np.abs(points[:, 1]).max() == pytest.approx(50.0, abs=1e-3)


def test_a_record_of_unknown_bases_only_draws_no_gc_track() -> None:
    record = _record("N" * 5_000)
    transform = _transform(5_000)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        frame = skew_df(record, WINDOW, STEP, "GC", circular=True)
        assert frame.isna().all().all() and len(frame) == 50
        assert generate_circular_gc_content_path_desc(100.0, 5_000, frame, 20.0, 1.0, "GC", transform) == ""
        assert generate_circular_gc_skew_path_desc(100.0, frame, 5_000, 20.0, 1.0, "GC", transform) == ""
        assert calculate_gc_content_path_desc(0.0, 0.0, frame, 5_000, 1_000.0, 100.0, "GC") == ""
        assert calculate_gc_skew_path_desc(0.0, 0.0, frame, 5_000, 1_000.0, 100.0) == ""
        percent = gc_content_percent_df(frame, dinucleotide="GC")
        assert percent.empty
        assert generate_circular_scalar_area_path_desc(100.0, 5_000, percent, 20.0, 1.0, record_transform=transform) == ""


def test_the_percent_track_draws_a_window_without_value_on_its_axis_minimum() -> None:
    frame = skew_df(_record(_gapped_sequence()), WINDOW, STEP, "GC", circular=True)
    percent = gc_content_percent_df(frame, dinucleotide="GC", min_percent=20.0, max_percent=80.0)
    no_value = frame["GC content"].isna()

    assert percent["percent"][no_value].isna().all()
    assert (percent["value_normalized"][no_value] == 0.0).all()
    assert percent["value_normalized"][~no_value].between(0.0, 1.0).all()


@pytest.mark.parametrize(("sequence", "expected"), [("GCNN", 100.0), ("gcatRYNN", 50.0), ("NNNN", 0.0)])
def test_the_record_gc_percent_reads_its_known_bases(sequence: str, expected: float) -> None:
    assert calculate_gc_percent(Seq(sequence)) == expected


def _drawn_tracks(mode: str, record: SeqRecord) -> dict[str, list[str]]:
    config = apply_config_overrides(None, {"canvas.show_gc": True, "canvas.show_skew": True})
    request_type, options_type = (
        (CircularDiagramRequest, CircularDiagramOptions) if mode == "circular"
        else (LinearDiagramRequest, LinearDiagramOptions)
    )
    request = request_type(
        records=(RecordInput(InMemoryRecordSource(record)),),
        options=options_type(config=config, window=WINDOW, step=STEP),
    )
    result = plan_request(request).build()
    root = ET.fromstring((result.drawing if hasattr(result, "drawing") else result).tostring())
    tracks: dict[str, list[str]] = {}
    for group in root.iter():
        track = {"gc_content": "GC content", "skew": "GC skew", "gc_skew": "GC skew"}.get(group.get("id") or "")
        if track is not None:
            tracks[track] = [node.get("d", "") for node in group.iter() if node.tag.endswith("path")]
    return tracks


def test_a_circular_diagram_draws_the_gap_on_the_baseline_circle() -> None:
    tracks = _drawn_tracks("circular", _record(_gapped_sequence()))
    assert set(tracks) == {"GC content", "GC skew"}
    for track, descs in tracks.items():
        desc = max(descs, key=lambda d: len(_POINT.findall(d)))
        baseline = float(_ARC_RADIUS.findall(desc)[-1])
        points = _line_points(desc)
        heights = np.hypot(points[:, 0], points[:, 1])
        positions = (np.arctan2(points[:, 0], -points[:, 1]) % (2.0 * math.pi)) / (2.0 * math.pi) * 10_000
        inside = (positions > 3_500) & (positions < 5_500)
        assert inside.any(), track
        np.testing.assert_allclose(heights[inside], baseline, atol=0.01, err_msg=track)


def test_a_linear_diagram_of_unknown_bases_only_draws_no_gc_track() -> None:
    tracks = _drawn_tracks("linear", _record("N" * 5_000, "linear"))
    assert all(not any(_POINT.findall(d) for d in descs) for descs in tracks.values()), tracks
