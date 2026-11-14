"""A feature that starts at the record's first base is drawn from that base.

Feature parts are half-open, 0-based bounds (Biopython's ``SimpleLocation``):
GenBank ``1..1000`` is ``[0, 1000)``. The block of a part therefore starts at
its Biopython start in both modes, whether the record is drawn from its
location or from a display start, and lane placement reads the same bounds.
The connector line between two parts spans the half-open gap between them,
``[previous end, next start)``, as the display-start path draws it.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw import RecordDisplayOptions, draw_circular, draw_linear, read_genbank
from gbdraw.features.coordinates import get_exon_and_intron_coordinates
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.features.tracks import arrange_feature_tracks, get_feature_ends
from gbdraw.layout.spatial import Interval, split_circular_interval

_LENGTH = 10_000
_NUMBER = re.compile(r"-?\d+(?:\.\d+)?")


def _record_file(tmp_path: Path, *locations: object) -> Path:
    record = SeqRecord(
        Seq("ACGT" * (_LENGTH // 4)), id="r1", name="r1", description="first base",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    record.features = [
        SeqFeature(location, type="CDS", qualifiers={"product": [f"p{index}"]})
        for index, location in enumerate(locations)
    ]
    path = tmp_path / "first_base.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _block_paths(svg: str) -> list[str]:
    """The path data of each feature block, in drawing order."""
    return [
        match.group(1)
        for match in re.finditer(r'<path d="([^"]+)"[^>]*data-gbdraw-feature-part="block"', svg)
    ]


def _points(d: str) -> list[tuple[float, float]]:
    """The end points of each path command (an arc's radii and flags are skipped)."""
    points: list[tuple[float, float]] = []
    for command in re.findall(r"[MLA][^MLAz]*", d):
        numbers = [float(n) for n in _NUMBER.findall(command)]
        points.append((numbers[-2], numbers[-1]))
    return points


def _linear_extent_bp(d: str, px_per_bp: float) -> tuple[float, float]:
    xs = [x for x, _ in _points(d)]
    return round(min(xs) / px_per_bp, 2), round(max(xs) / px_per_bp, 2)


def _circular_extent_bp(d: str) -> tuple[float, float]:
    """The record positions the block covers, from the angle of each point (12 o'clock is 0, clockwise)."""
    positions = []
    for x, y in _points(d):
        position = (math.atan2(x, -y) / (2 * math.pi) * _LENGTH) % _LENGTH
        positions.append(0.0 if position > _LENGTH - 0.05 else position)
    return round(min(positions), 1), round(max(positions), 1)


_FIRST = SimpleLocation(0, 1000, 1)  # 1..1000
_SECOND = SimpleLocation(1, 1000, -1)  # complement(2..1000)


def test_linear_whole_record_draws_a_block_at_the_first_base_from_that_base(tmp_path: Path) -> None:
    records = read_genbank(_record_file(tmp_path, _FIRST, _SECOND))
    first, second = _block_paths(draw_linear(records).to_svg())
    px_per_bp = 2000 / _LENGTH  # The default alignment width over the record.

    assert _linear_extent_bp(first, px_per_bp) == (0.0, 1000.0)
    assert _linear_extent_bp(second, px_per_bp) == (1.0, 1000.0)


def test_linear_location_path_and_display_start_path_draw_the_same_block(tmp_path: Path) -> None:
    records = read_genbank(_record_file(tmp_path, _FIRST, _SECOND))
    whole = _block_paths(draw_linear(records).to_svg())
    shown = _block_paths(
        draw_linear(records, record_displays=[RecordDisplayOptions(start_coordinate=4001)]).to_svg()
    )
    px_per_bp = 2000 / _LENGTH
    for drawn_whole, drawn_shown in zip(whole, shown, strict=True):
        start, end = _linear_extent_bp(drawn_whole, px_per_bp)
        shown_start, shown_end = _linear_extent_bp(drawn_shown, px_per_bp)
        # Base 4001 is drawn at position 0: each block moves by 6,000 bp and keeps its length.
        assert (shown_start, shown_end) == (start + 6000, end + 6000)


def test_circular_draws_a_block_at_the_first_base_from_that_base(tmp_path: Path) -> None:
    records = read_genbank(_record_file(tmp_path, _FIRST, _SECOND))
    first, second = _block_paths(draw_circular(records).to_svg())

    assert _circular_extent_bp(first) == (0.0, 1000.0)
    assert _circular_extent_bp(second) == (1.0, 1000.0)


def test_circular_location_path_and_display_start_path_draw_the_same_block(tmp_path: Path) -> None:
    records = read_genbank(_record_file(tmp_path, _FIRST, _SECOND))
    whole = _block_paths(draw_circular(records).to_svg())
    shown = _block_paths(
        draw_circular(records, record_displays=[RecordDisplayOptions(start_coordinate=4001)]).to_svg()
    )
    for drawn_whole, drawn_shown in zip(whole, shown, strict=True):
        start, end = _circular_extent_bp(drawn_whole)
        assert _circular_extent_bp(drawn_shown) == (start + 6000, end + 6000)


def _lines(location: list[SimpleLocation]) -> list[tuple[int, int]]:
    return [
        (part.start, part.end)
        for part in get_exon_and_intron_coordinates(location)
        if part.kind == "line"
    ]


def test_a_connector_line_spans_the_gap_between_two_parts() -> None:
    # join(1..100,201..300) and complement(join(1..100,201..300)): the gap is bases 101..200, [100, 200).
    assert _lines([SimpleLocation(0, 100, 1), SimpleLocation(200, 300, 1)]) == [(100, 200)]
    assert _lines([SimpleLocation(200, 300, -1), SimpleLocation(0, 100, -1)]) == [(100, 200)]


def test_a_connector_line_from_the_first_base_starts_where_the_part_ends() -> None:
    # join(1..1,3..10): the gap is base 2, [1, 2).
    assert _lines([SimpleLocation(0, 1, 1), SimpleLocation(2, 10, 1)]) == [(1, 2)]
    # Parts over the origin of an unrotated record have no connection (as the display path draws them).
    assert _lines([SimpleLocation(9000, 9500, 1), SimpleLocation(0, 500, 1)]) == []


def _feature(feature_id: str, *parts: tuple[int, int]) -> FeatureObject:
    location = [
        FeatureLocationPart("block", str(i).zfill(3), "positive", start, end, i == len(parts))
        for i, (start, end) in enumerate(parts, start=1)
    ]
    return FeatureObject(
        feature_id=feature_id, location=location, is_directional=True, color="#54bcf8",
        note="", label_text="", coordinates=location, type="CDS", qualifiers={},
    )


def test_lane_placement_reads_a_part_at_the_first_base_from_that_base() -> None:
    assert get_feature_ends(_feature("a", (0, 1000)), _LENGTH) == (0, 1000, "positive")
    # Only a part at the first base touches the origin; one at base 2 does not.
    assert get_feature_ends(_feature("o", (9989, _LENGTH), (0, 20)), _LENGTH) == (9989, 20, "positive")
    assert get_feature_ends(_feature("n", (9989, _LENGTH), (1, 20)), _LENGTH) == (1, _LENGTH, "positive")
    assert split_circular_interval(9989, 20, _LENGTH) == [Interval(9989, _LENGTH), Interval(0, 20)]


@pytest.mark.parametrize("genome_length", [_LENGTH])
def test_two_features_that_share_only_the_first_base_take_two_lanes(genome_length: int) -> None:
    # 1..1 and join(9990..10000,1..1) both cover base 1.
    features = {"a": _feature("a", (0, 1)), "o": _feature("o", (9989, _LENGTH), (0, 1))}
    arranged = arrange_feature_tracks(features, False, True, genome_length=genome_length)
    assert arranged["a"].feature_track_id != arranged["o"].feature_track_id
