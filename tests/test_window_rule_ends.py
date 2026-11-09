"""Sliding windows read the bases next to them on the molecule, also at a record's ends (OV-2890, OV-2919).

GC content, GC skew and Depth take their windows from one rule (``gbdraw.analysis.windows``): the window of the
sample at ``p`` is ``[p - w // 2, p - w // 2 + w)``, centred on it. Past an end of the drawn record the window
reads what lies there on the molecule: a circular record's windows wrap round its origin, a crop (a region) reads
its flanking bases from the record it was cut from, and a linear record's windows stop at its ends, where a value
is the mean of the bases inside.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.depth import depth_df
from gbdraw.analysis.skew import skew_df
from gbdraw.api.diagram import build_linear_diagram
from gbdraw.api.options import DepthTrackInput, LinearDiagramOptions
from gbdraw.io.regions import RegionSpec, crop_record_to_region
from tests.test_gc_window_centre import STEP, WINDOW, _seen_positions, _track_paths

LENGTH = 5_000
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _sequence() -> str:
    rng = random.Random(2890)
    return "".join(rng.choice("ACGT") for _ in range(LENGTH))


def _record(sequence: str, topology: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id="src", annotations={"topology": topology, "molecule_type": "DNA"})


def _region(start: int, end: int, reverse: bool = False) -> RegionSpec:
    return RegionSpec(f"src:{start}-{end}", None, "src", None, start, end, reverse)


def _window_indices(first: int, length: int, reverse: bool, circular: bool, n: int, window: int, step: int):
    """The source indices each sample's window reads: local ``x`` is source index ``first + x`` (``first - x``
    for a reverse complement), read round the circle or kept inside ``0 .. n - 1``."""

    for p in range(0, length, step):
        indices = []
        for x in range(p - window // 2, p - window // 2 + window):
            i = first - x if reverse else first + x
            if circular:
                indices.append(i % n)
            elif 0 <= i < n:
                indices.append(i)
        yield p, indices


def _expected_gc(source: str, first: int, length: int, *, reverse: bool, circular: bool, window: int, step: int):
    rows = {}
    for p, indices in _window_indices(first, length, reverse, circular, len(source), window, step):
        bases = "".join(source[i] for i in indices)
        if reverse:
            bases = bases.translate(_COMPLEMENT)
        g, c = bases.count("G"), bases.count("C")
        rows[p] = ((g + c) / len(indices), 0.0 if g + c == 0 else (g - c) / (g + c))
    return rows


# id, topology, the regions cropped one after the other, the first source index of the drawn record, whether it is
# reverse-complemented, and the window.
_CASES = [
    ("circular", "circular", (), 0, False, 600),
    ("linear", "linear", (), 0, False, 600),
    ("crop-circular", "circular", ((1001, 3500, False),), 1000, False, 600),
    ("crop-linear", "linear", ((1001, 3500, False),), 1000, False, 600),
    ("crop-linear-rc", "linear", ((1001, 3500, True),), 3499, True, 600),
    ("crop-circular-rc", "circular", ((1001, 3500, True),), 3499, True, 600),
    # A crop whose start lies less than half a window after the start of a linear record: its first windows read
    # the 200 bases before it and stop at the record's start.
    ("crop-near-start-linear", "linear", ((201, 2000, False),), 200, False, 600),
    ("crop-near-start-circular", "circular", ((201, 2000, False),), 200, False, 600),
    # A crop across the origin of a circular record reads on round the origin, also reverse-complemented.
    ("crop-across-origin", "circular", ((4501, 500, False),), 4500, False, 600),
    ("crop-across-origin-rc", "circular", ((4501, 500, True),), 499, True, 600),
    # A crop of a crop reads the flanks of the first record: local 501..1500 of the crop 1001..3500.
    ("crop-of-crop-linear", "linear", ((1001, 3500, False), (501, 1500, False)), 1500, False, 600),
    ("crop-of-crop-circular", "circular", ((1001, 3500, False), (501, 1500, False)), 1500, False, 600),
    # A window longer than the drawn record: a crop of 400 bases, and a whole record of 5,000 with 7,001-base
    # windows (round a circular record more than once).
    ("long-window-crop", "linear", ((2001, 2400, False),), 2000, False, 1_001),
    ("long-window-linear", "linear", (), 0, False, 7_001),
    ("long-window-circular", "circular", (), 0, False, 7_001),
]
_IDS = [case[0] for case in _CASES]


def _drawn(topology: str, regions: tuple[tuple[int, int, bool], ...]) -> tuple[str, SeqRecord]:
    source = _sequence()
    record = _record(source, topology)
    for region in regions:
        record = crop_record_to_region(record, _region(*region), circular=topology == "circular")
    return source, record


def _step(kind: str, window: int) -> int:
    """A fine step, or one wider than half the window, so the samples leave bases of the record unread."""
    return 50 if kind == "fine" else window + 37


@pytest.mark.parametrize("step_kind", ["fine", "wide"])
@pytest.mark.parametrize("case", _CASES, ids=_IDS)
def test_gc_windows_read_the_bases_next_to_them_on_the_molecule(case, step_kind: str) -> None:
    _id, topology, regions, first, reverse, window = case
    step = _step(step_kind, window)
    source, record = _drawn(topology, regions)
    frame = skew_df(record, window, step, "GC", circular=topology == "circular")
    expected = _expected_gc(
        source, first, len(record.seq), reverse=reverse, circular=topology == "circular", window=window, step=step
    )
    assert frame.index.tolist() == list(expected)
    assert frame["GC content"].tolist() == pytest.approx([content for content, _ in expected.values()], abs=1e-12)
    assert frame["GC skew"].tolist() == pytest.approx([skew for _, skew in expected.values()], abs=1e-12)


def test_the_test_sequence_is_mixed() -> None:
    """The brute-force cases above see a wrong window only on a sequence whose windows differ."""

    source = _sequence()
    assert {source.count(base) for base in "ACGT"} != {0} and min(source.count(base) for base in "ACGT") > 1_000


def test_a_crops_windows_equal_the_same_windows_of_its_source() -> None:
    """The crop 1001..3500 draws source bases 1000..3499; its samples lie on the source's sample grid."""

    source = _sequence()
    for topology in ("circular", "linear"):
        record = _record(source, topology)
        crop = crop_record_to_region(record, _region(1001, 3500), circular=topology == "circular")
        whole = skew_df(record, 600, 50, "GC", circular=topology == "circular")
        cropped = skew_df(crop, 600, 50, "GC", circular=topology == "circular")
        for column in ("GC content", "GC skew"):
            assert cropped[column].to_numpy() == pytest.approx(whole.loc[1000:3499, column].to_numpy(), abs=0)


def _depths() -> np.ndarray:
    rng = random.Random(2919)
    return np.array([rng.randint(0, 60) for _ in range(LENGTH)], dtype=np.int64)


def _depth_table(depths: np.ndarray, rows: range | None = None) -> pd.DataFrame:
    positions = np.arange(1, len(depths) + 1) if rows is None else np.array(rows) + 1
    return pd.DataFrame(
        {"reference_name": ["src"] * len(positions), "position": positions, "depth": depths[positions - 1]}
    )


def _expected_depths(depths, rows: range, drawn: set[int], case, step: int) -> dict[int, float]:
    """Each window's mean over the bases it reads: a drawn base, or a flank base within ``rows`` (the source
    indexes from the first to the last row); a base inside both with no row adds 0 (D-64)."""

    _id, topology, _regions, first, reverse, window = case
    held = np.zeros(LENGTH)
    held[list(rows)] = depths[list(rows)]
    expected = {}
    length = len(drawn)
    for p, indices in _window_indices(first, length, reverse, topology == "circular", LENGTH, window, step):
        read = [i for i in indices if i in drawn or rows.start <= i < rows.stop]
        expected[p] = float(held[read].sum()) / len(read)
    return expected


@pytest.mark.parametrize("step_kind", ["fine", "wide"])
@pytest.mark.parametrize("rows", [range(0, LENGTH), range(0, 3_200), range(1_200, 3_300)], ids=["all", "to-3200", "1200-3300"])
@pytest.mark.parametrize("case", _CASES, ids=_IDS)
def test_depth_windows_follow_the_same_rule(case, rows: range, step_kind: str) -> None:
    """Depth reads the same windows; past the drawn record it reads only where the rows of the file reach, so the
    rows of the drawn region alone draw no fall at its ends (D-64)."""

    _id, topology, regions, first, reverse, window = case
    depths = _depths()
    _source, record = _drawn(topology, regions)
    drawn = {(first - x if reverse else first + x) % LENGTH for x in range(len(record.seq))}
    step = _step(step_kind, window)
    frame = depth_df(record, _depth_table(depths, rows), window, step, circular=topology == "circular")
    expected = _expected_depths(depths, rows, drawn, case, step)
    assert frame["position"].tolist() == list(expected)
    assert frame["depth"].tolist() == pytest.approx(list(expected.values()), rel=1e-12)


def test_a_wide_step_still_counts_the_end_of_a_circular_record() -> None:
    """A step wider than half the window leaves the end of the record unread by any sample; the first window still
    reads it across the origin, and its positions without a row count 0."""

    record = SeqRecord(Seq("ACGT" * 2_500), id="c", annotations={"topology": "circular", "molecule_type": "DNA"})
    table = pd.DataFrame({"reference_name": ["c"] * 9_000, "position": range(1, 9_001), "depth": [50.0] * 9_000})
    frame = depth_df(record, table, 100, 5_000, circular=True)
    assert frame["depth"].tolist()[0] == 25.0


def test_rows_of_the_drawn_region_alone_draw_a_flat_line() -> None:
    """A flat region-only depth file (``samtools depth -r``) draws flat, as before the flanks were read (D-64)."""

    record = crop_record_to_region(_record("A" * 100_000, "linear"), _region(40_001, 60_000))
    rows = pd.DataFrame({"reference_name": "src", "position": np.arange(40_001, 60_001), "depth": 50})
    assert set(depth_df(record, rows, 2_000, 200, circular=False)["depth"]) == {50.0}


# Through the drawing: the samples whose windows see a G run stand out from an all-A record.


def _run_record(run: tuple[int, int], topology: str, length: int = 10_000) -> SeqRecord:
    bases = ["A"] * length
    bases[run[0]:run[1]] = ["G"] * (run[1] - run[0])
    return SeqRecord(Seq("".join(bases)), id="run", annotations={"topology": topology, "molecule_type": "DNA"})


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_linear_records_end_windows_do_not_read_its_other_end(mode: str) -> None:
    """A run at the start of a linear record is seen within half a window of the start only (circular: and
    within half a window before the end). Past its last sample a Circular line of a linear record holds that
    sample's value to 12 o'clock, as a Linear track holds it to the end; it does not ramp to the first sample's
    value, which would draw bases that are not there (OV-2920)."""

    for topology in ("linear", "circular"):
        paths = _track_paths(mode, _run_record((0, 100), topology), None)
        for track, desc in paths.items():
            seen = _seen_positions(mode, desc, 10_000)
            near_end = [p for p in seen if p > 10_000 - WINDOW]
            assert bool(near_end) == (topology == "circular"), (mode, topology, track, near_end)
            assert all(p < WINDOW // 2 + 100 + STEP or p > 10_000 - WINDOW for p in seen), (mode, topology, track)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_a_crop_reads_its_flanks_from_the_source(mode: str, topology: str) -> None:
    """The crop 3001..6000 draws bases 3000..5999; the run 2800..2900 lies before it, within half a window of
    its first base, so the crop's first samples see it, and its last samples do not."""

    region = RegionSpec("run:3001-6000", None, "run", None, 3001, 6000, False)
    paths = _track_paths(mode, _run_record((2800, 2900), topology), region)
    for track, desc in paths.items():
        seen = _seen_positions(mode, desc, 3000)
        assert seen, (mode, topology, track)
        # Local 0 is base 3000: the windows of samples below 2900 - 3000 + 500 = 400 reach the run.
        assert max(seen) < 400 + STEP / 2 and min(seen) > -STEP / 2 or mode == "circular" and all(
            p < 400 + STEP / 2 or p > 3000 - STEP / 2 for p in seen
        ), (mode, topology, track, seen)


def test_depth_is_drawn_at_the_window_centre() -> None:
    """One deep position at base 3000 of a 10 kb record: the samples whose windows hold it are centred on it."""

    depths = np.zeros(10_000, dtype=np.int64)
    depths[3000] = 1_000
    record = _run_record((0, 0), "linear")
    record.id = "run"
    svg = build_linear_diagram(
        [record],
        options=LinearDiagramOptions(
            config_overrides={"canvas.show_gc": False, "canvas.show_skew": False},
            depth_tracks=(DepthTrackInput(source=_depth_table(depths).assign(reference_name="run")),),
            depth_window=2_000,
            depth_step=100,
        ),
    ).tostring()
    from xml.etree import ElementTree as ET

    root = ET.fromstring(svg)
    group = next(node for node in root.iter() if node.get("id") == "depth")
    desc = max((node.get("d", "") for node in group.iter() if node.tag.endswith("path")), key=len)
    seen = _seen_positions("linear", desc, 10_000)
    assert seen
    centre = sum(seen) / len(seen)
    # Start-anchored windows put the deep base half a window (1,000 bp) before its place.
    assert math.isclose(centre, 3_000, abs_tol=STEP), seen


@pytest.mark.parametrize("close_at_record_len", [False, True], ids=["depth", "percent"])
def test_a_circular_area_of_a_linear_record_holds_its_last_value_to_the_end(close_at_record_len: bool) -> None:
    """Depth and Percent GC areas: a record drawn as circular closes on its first sample (Percent GC) or ends at
    its last sample (Depth) as before; a linear record's area holds its last value on to 12 o'clock (OV-2920)."""

    from gbdraw.layout.record_transform import unset_record_transform
    from gbdraw.svg.circular_tracks import circular_scalar_area_line

    data = pd.DataFrame({"position": [0, 400, 800], "value_normalized": [0.2, 0.5, 0.9]})
    lines = {}
    for topology in ("linear", "circular"):
        transform = unset_record_transform(_run_record((0, 0), topology, length=1_000))
        line = circular_scalar_area_line(
            100.0, 1_000, data, 10.0, 1.0, close_at_record_len=close_at_record_len, record_transform=transform
        )
        assert line is not None
        lines[topology] = line
    linear = lines["linear"]
    assert list(linear.positions) == [0, 400, 800, 1_000]
    assert linear.radii[-1] == linear.radii[-2]
    circular = lines["circular"]
    if close_at_record_len:
        assert list(circular.positions) == [0, 400, 800, 1_000] and circular.radii[-1] == circular.radii[0]
    else:
        assert list(circular.positions) == [0, 400, 800]


@pytest.mark.parametrize("track", ["gc_content", "gc_skew", "depth", "percent"])
def test_a_circular_crop_of_a_circular_record_holds_its_last_value_to_12_oclock(track: str) -> None:
    """A crop of a circular record has no origin in it: its ends are two places on the molecule, so past its last
    sample its Circular line holds that value to 12 o'clock and steps there, as a linear record's does; it does
    not ramp to its first sample's value (D-64). The whole circular record still closes on its first sample."""

    from gbdraw.layout.record_transform import unset_record_transform
    from gbdraw.svg.circular_tracks import circular_gc_content_line, circular_gc_skew_line, circular_scalar_area_line

    source = _run_record((0, 0), "circular", length=10_000)
    crop = crop_record_to_region(source, _region(3_001, 4_000), circular=True)
    data = pd.DataFrame({"position": [0, 400, 800], "value_normalized": [0.2, 0.5, 0.9]})
    frame = pd.DataFrame({"GC content": [0.2, 0.5, 0.9], "GC skew": [-0.2, 0.1, 0.4]}, index=[0, 400, 800])
    lines = {}
    for name, record in (("crop", crop), ("whole", _run_record((0, 0), "circular", length=1_000))):
        transform = unset_record_transform(record)
        if track == "gc_content":
            line = circular_gc_content_line(100.0, 1_000, frame, 10.0, 1.0, "GC", transform)
        elif track == "gc_skew":
            line = circular_gc_skew_line(100.0, frame, 1_000, 10.0, 1.0, "GC", transform)
        else:
            line = circular_scalar_area_line(
                100.0, 1_000, data, 10.0, 1.0, close_at_record_len=track == "percent", record_transform=transform
            )
        assert line is not None
        lines[name] = (list(line.positions), list(line.radii))
    positions, radii = lines["crop"]
    assert positions == [0, 400, 800, 1_000] and radii[-1] == radii[-2], lines["crop"]
    whole_positions, whole_radii = lines["whole"]
    if track in ("gc_content", "percent"):
        assert whole_positions == [0, 400, 800, 1_000] and whole_radii[-1] == whole_radii[0], lines["whole"]
    else:
        # The skew line closes along the circle in its path (``_closed_line_desc``); Depth ends at its last sample.
        assert whole_positions == [0, 400, 800], lines["whole"]
