"""A feature that starts at the record's first base is drawn from that base.

Feature parts are half-open, 0-based bounds (Biopython's ``SimpleLocation``):
GenBank ``1..1000`` is ``[0, 1000)``. The block of a part therefore starts at
its Biopython start in both modes, whether the record is drawn from its
location or from a display start, and lane placement reads the same bounds.
The connector line between two parts spans the half-open gap between them,
``[previous end, next start)``, as the display-start path draws it; on a
circular record two parts of one strand that cross the origin connect over it
in every view. A block over the whole record is a full ring.
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
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.features.coordinates import get_exon_and_intron_coordinates
from gbdraw.features.objects import FeatureObject
from gbdraw.features.tracks import arrange_feature_tracks, get_feature_ends
from gbdraw.io.regions import parse_region_spec
from gbdraw.layout.spatial import Interval, split_circular_interval
from gbdraw.svg.circular_features import (
    generate_circular_arrow_path_with_radii,
    generate_circular_narrow_arrow_path_with_radii,
)

_LENGTH = 10_000
_NUMBER = re.compile(r"-?\d+(?:\.\d+)?")


def _record_file(tmp_path: Path, *locations: object, topology: str = "circular") -> Path:
    record = SeqRecord(
        Seq("ACGT" * (_LENGTH // 4)), id="r1", name="r1", description="first base",
        annotations={"molecule_type": "DNA", "topology": topology},
    )
    record.features = [
        SeqFeature(location, type="CDS", qualifiers={"product": [f"p{index}"]})
        for index, location in enumerate(locations)
    ]
    path = tmp_path / "first_base.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _block_paths(svg: str, part: str = "block") -> list[str]:
    """The path data of each feature block (or ``connector`` line), in drawing order."""
    return [
        match.group(1)
        for match in re.finditer(rf'<path d="([^"]+)"[^>]*data-gbdraw-feature-part="{part}"', svg)
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


def test_a_region_that_starts_inside_a_feature_draws_its_block_from_the_first_drawn_base(tmp_path: Path) -> None:
    # GenBank 191..736 cut by the region 701..760 (OV-3652): the crop shifts the
    # part to [0, 36), so its block starts at the region's first base.
    path = _record_file(tmp_path, SimpleLocation(190, 736, 1), topology="linear")
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(read_genbank(path)[0]), region=parse_region_spec("r1:701-760")),),
    )
    drawing = plan_request(request).build().drawing
    block = _block_paths(drawing.tostring())[0]

    assert _linear_extent_bp(block, 2000 / 60) == (0.0, 36.0)


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


def _lines(location: list[SimpleLocation], origin_length: int | None = None) -> list[tuple[int, int]]:
    return [
        (part.start, part.end)
        for part in get_exon_and_intron_coordinates(location, origin_length=origin_length)
        if part.kind == "line"
    ]


def test_a_connector_line_spans_the_gap_between_two_parts() -> None:
    # join(1..100,201..300) and complement(join(1..100,201..300)): the gap is bases 101..200, [100, 200).
    assert _lines([SimpleLocation(0, 100, 1), SimpleLocation(200, 300, 1)]) == [(100, 200)]
    assert _lines([SimpleLocation(200, 300, -1), SimpleLocation(0, 100, -1)]) == [(100, 200)]


def test_a_connector_line_from_the_first_base_starts_where_the_part_ends() -> None:
    # join(1..1,3..10): the gap is base 2, [1, 2).
    assert _lines([SimpleLocation(0, 1, 1), SimpleLocation(2, 10, 1)]) == [(1, 2)]
    # Parts over the origin connect over it only on a record whose spans may cross its origin.
    assert _lines([SimpleLocation(9000, 9500, 1), SimpleLocation(0, 500, 1)]) == []
    assert _lines([SimpleLocation(9000, 9500, 1), SimpleLocation(0, 500, 1)], _LENGTH) == [(9500, _LENGTH)]


def test_parts_of_one_strand_over_the_origin_connect_over_it_on_a_circular_record() -> None:
    # join(9501..9990,11..500): the gaps are 9991..10000 and 1..10.
    plus = [SimpleLocation(9500, 9990, 1), SimpleLocation(10, 500, 1)]
    assert _lines(plus, _LENGTH) == [(9990, _LENGTH), (0, 10)]
    # complement(join(9501..9990,11..500)): Biopython lists the parts 11..500 first.
    minus = [SimpleLocation(10, 500, -1), SimpleLocation(9500, 9990, -1)]
    assert _lines(minus, _LENGTH) == [(0, 10), (9990, _LENGTH)]
    # Parts that meet at the origin leave no gap; parts of two strands connect by the first part's strand only.
    assert _lines([SimpleLocation(9500, _LENGTH, 1), SimpleLocation(0, 500, 1)], _LENGTH) == []
    assert _lines([SimpleLocation(9500, 9990, 1), SimpleLocation(10, 500, -1)], _LENGTH) == []


def _feature(feature_id: str, *parts: tuple[int, int], trans_spliced: bool = False) -> FeatureObject:
    """A plus-strand feature on a circular record of ``_LENGTH`` bp, its blocks and lines built as the factory builds them."""
    location = get_exon_and_intron_coordinates(
        [SimpleLocation(start, end, 1) for start, end in parts], trans_spliced, origin_length=_LENGTH
    )
    return FeatureObject(
        feature_id=feature_id, location=location, is_directional=True, color="#54bcf8",
        note="", label_text="", coordinates=location, type="CDS",
        qualifiers={"trans_splicing": [""]} if trans_spliced else {},
    )


def test_lane_placement_reads_a_part_at_the_first_base_from_that_base() -> None:
    assert get_feature_ends(_feature("a", (0, 1000)), _LENGTH) == (0, 1000, "positive")
    # Only a part at the first base touches the origin; one at base 2 does not.
    assert get_feature_ends(_feature("o", (9989, _LENGTH), (0, 20)), _LENGTH) == (9989, 20, "positive")
    # A part at base 2 after one that runs to the end still crosses the origin: the line over base 1 reaches it.
    assert get_feature_ends(_feature("n", (9989, _LENGTH), (1, 20)), _LENGTH) == (9989, 20, "positive")
    assert split_circular_interval(9989, 20, _LENGTH) == [Interval(9989, _LENGTH), Interval(0, 20)]


@pytest.mark.parametrize("genome_length", [_LENGTH])
def test_two_features_that_share_only_the_first_base_take_two_lanes(genome_length: int) -> None:
    # 1..1 and join(9990..10000,1..1) both cover base 1.
    features = {"a": _feature("a", (0, 1)), "o": _feature("o", (9989, _LENGTH), (0, 1))}
    arranged = arrange_feature_tracks(features, False, True, genome_length=genome_length)
    assert arranged["a"].feature_track_id != arranged["o"].feature_track_id


def test_a_join_over_the_origin_with_a_gap_at_the_first_base_keeps_its_lane() -> None:
    # join(9901..10000,2..100) next to 1001..1500 and 1601..2200: nothing overlaps, so all take one lane.
    features = {
        "o": _feature("o", (9900, _LENGTH), (1, 100)),
        "a": _feature("a", (1000, 1500)),
        "b": _feature("b", (1600, 2200)),
    }
    arranged = arrange_feature_tracks(features, False, True, genome_length=_LENGTH)
    assert {feature.feature_track_id for feature in arranged.values()} == {0}
    # A join over the origin of three parts reserves all three.
    three = _feature("t", (9000, 9500), (9600, _LENGTH), (0, 20))
    assert get_feature_ends(three, _LENGTH) == (9000, 20, "positive")


def test_the_part_order_of_a_trans_spliced_feature_does_not_cross_the_origin() -> None:
    # A trans-spliced rps12: its 5' part lies after its 3' parts, which is not a path over the origin.
    rps12 = _feature("r", (7000, 7100), (3000, 3200), (3500, 3600), trans_spliced=True)
    assert get_feature_ends(rps12, _LENGTH) == (3000, 7100, "positive")
    # Its parts still cross the origin where they reach both ends of the record.
    rps12_origin = _feature("o", (9900, _LENGTH), (0, 100), trans_spliced=True)
    assert get_feature_ends(rps12_origin, _LENGTH) == (9900, 100, "positive")


def _arc_sweeps(d: str) -> list[float]:
    """The signed sweep in degrees of each arc command (positive clockwise on screen)."""
    sweeps: list[float] = []
    current = (0.0, 0.0)
    for command in re.findall(r"[MLAz][^MLAz]*", d):
        numbers = [float(n) for n in _NUMBER.findall(command)]
        if command[0] == "A":
            sweep_flag, end = int(numbers[4]), (numbers[5], numbers[6])
            start_deg = math.degrees(math.atan2(current[1], current[0]))
            end_deg = math.degrees(math.atan2(end[1], end[0]))
            delta = (end_deg - start_deg) % 360.0
            sweeps.append(delta if sweep_flag == 1 else delta - 360.0)
            current = end
        elif numbers:
            current = (numbers[-2], numbers[-1])
    return sweeps


@pytest.mark.parametrize("strand", ["positive", "negative"])
@pytest.mark.parametrize("narrow", [False, True], ids=["arrow", "narrow_arrow"])
def test_a_circular_arrow_over_the_whole_record_is_a_full_ring(strand: str, narrow: bool) -> None:
    # 1..10000 and complement(1..10000): [0, 10000). The shaft runs round the ring up to the head, which ends at
    # the shaft's tail.
    coord = {"coord_type": "block", "coord_strand": strand, "coord_start": 0, "coord_end": _LENGTH}
    head_bp = 300
    if narrow:
        path = generate_circular_narrow_arrow_path_with_radii(coord, _LENGTH, head_bp, 90.0, 100.0, 110.0, 0.5)[1]
    else:
        path = generate_circular_arrow_path_with_radii(coord, _LENGTH, head_bp, 90.0, 100.0, 110.0)[1]
    sweeps = _arc_sweeps(path)
    shaft_deg = 360.0 * (_LENGTH - head_bp) / _LENGTH
    direction = 1.0 if strand == "positive" else -1.0
    assert sum(s for s in sweeps if s * direction > 0) == pytest.approx(direction * shaft_deg, abs=0.01)
    assert sum(s for s in sweeps if s * direction < 0) == pytest.approx(-direction * shaft_deg, abs=0.01)


def test_a_circular_feature_over_the_whole_record_is_drawn_round_the_ring(tmp_path: Path) -> None:
    whole = SimpleLocation(0, _LENGTH, 1)
    svg = draw_circular(read_genbank(_record_file(tmp_path, whole, SimpleLocation(0, _LENGTH, -1)))).to_svg()
    for d in _block_paths(svg):
        sweeps = _arc_sweeps(d)
        assert sum(abs(s) for s in sweeps) > 2 * 0.9 * 360.0


def test_a_circular_record_draws_its_connector_over_the_origin_in_every_view(tmp_path: Path) -> None:
    # join(9701..9990,1..300): the gap 9991..10000 is drawn with and without a display start.
    records = read_genbank(_record_file(tmp_path, SimpleLocation(9700, 9990, 1) + SimpleLocation(0, 300, 1)))
    whole = [_circular_extent_bp(d) for d in _block_paths(draw_circular(records).to_svg(), "connector")]
    shown_svg = draw_circular(records, record_displays=[RecordDisplayOptions(start_coordinate=1001)]).to_svg()
    shown = [_circular_extent_bp(d) for d in _block_paths(shown_svg, "connector")]
    # Its ends are base 9990 and the origin, which the extent reads as 0.
    assert whole == [(0.0, 9990.0)]
    # Base 1001 is drawn at position 0: the line moves back by 1,000 bp.
    assert shown == [(8990.0, 9000.0)]
    px_per_bp = 2000 / _LENGTH
    linear = [_linear_extent_bp(d, px_per_bp) for d in _block_paths(draw_linear(records).to_svg(), "connector")]
    assert linear == [(9990.0, 10000.0)]


def test_a_linear_record_draws_no_connector_over_the_origin(tmp_path: Path) -> None:
    records = read_genbank(
        _record_file(tmp_path, SimpleLocation(9700, 9990, 1) + SimpleLocation(0, 300, 1), topology="linear")
    )
    assert _block_paths(draw_linear(records).to_svg(), "connector") == []


# DnaFeaturesViewer#85: CJ_0003 join(complement(1501..1700),complement(1801..2000)) lists its parts low to high on
# the minus strand. It overlaps nothing and covers 1501..2000; it does not run round the ring.
_DFV85 = Path(__file__).parent / "fixtures" / "input_variants" / "dfv85_complement_join.gb"


def _dfv85(topology: str) -> RecordInput:
    record = SeqIO.read(_DFV85, "genbank")
    record.annotations["topology"] = topology
    return RecordInput(InMemoryRecordSource(record))


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_a_join_listed_low_to_high_on_the_minus_strand_keeps_one_lane_in_linear(topology: str) -> None:
    options = LinearDiagramOptions(config=apply_config_overrides(None, {"canvas.resolve_overlaps": True}))
    svg = plan_request(LinearDiagramRequest(records=(_dfv85(topology),), options=options)).build().drawing.tostring()
    assert {min(y for _, y in _points(d)) for d in _block_paths(svg)} == {-20.0}
    assert len(_block_paths(svg, "connector")) == 3


@pytest.mark.parametrize("resolve_overlaps", [False, True])
@pytest.mark.parametrize("start", [None, 1001])
def test_a_join_listed_low_to_high_on_the_minus_strand_does_not_cross_the_origin(
    resolve_overlaps: bool, start: int | None
) -> None:
    options = CircularDiagramOptions(config=apply_config_overrides(None, {"canvas.resolve_overlaps": resolve_overlaps}))
    record = _dfv85("circular")
    if start is not None:
        record = RecordInput(record.source, display=RecordDisplayOptions(start_coordinate=start))
    svg = plan_request(CircularDiagramRequest(records=(record,), options=options)).build().tostring()
    # One lane: every block's inner edge lies on one circle.
    assert len({round(min(math.hypot(x, y) for x, y in _points(d)), 1) for d in _block_paths(svg)}) == 1
    # The connector lines span the three gaps, 400 of 3,000 bp (48 degrees), and none runs round the ring.
    sweep = sum(abs(s) for d in _block_paths(svg, "connector") for s in _arc_sweeps(d))
    assert sweep == pytest.approx(48.0, abs=0.5)


def _ends_record(topology: str, *, trans_spliced: bool) -> RecordInput:
    """3,000 bp: a join over 201..500 and 2501..2800 that draws nothing over the record's ends, and CDSs at both ends."""
    record = SeqRecord(
        Seq("ACGT" * 750), id="ends", name="ends", description="ends",
        annotations={"molecule_type": "DNA", "topology": topology},
    )
    parts = [SimpleLocation(200, 500, 1), SimpleLocation(2500, 2800, 1)]
    join = SeqFeature(parts[0] + parts[1] if trans_spliced else parts[1] + parts[0], type="CDS")
    if trans_spliced:
        join.qualifiers["trans_splicing"] = [""]
    record.features = [join, SeqFeature(SimpleLocation(20, 150, 1), type="CDS"),
                       SeqFeature(SimpleLocation(2850, 2980, 1), type="CDS")]
    return RecordInput(InMemoryRecordSource(record))


@pytest.mark.parametrize("topology", ["linear", "circular"])
@pytest.mark.parametrize("trans_spliced", [True, False], ids=["trans_spliced", "other_order"])
def test_a_join_that_draws_nothing_over_the_record_ends_leaves_them_to_other_features(
    topology: str, trans_spliced: bool
) -> None:
    options = LinearDiagramOptions(config=apply_config_overrides(None, {"canvas.resolve_overlaps": True}))
    request = LinearDiagramRequest(records=(_ends_record(topology, trans_spliced=trans_spliced),), options=options)
    svg = plan_request(request).build().drawing.tostring()
    lanes = {min(y for _, y in _points(d)) for d in _block_paths(svg)}
    # A trans-spliced join draws no line; parts in the other order connect over the origin only on a circular record.
    assert lanes == ({-20.0, -72.0} if topology == "circular" and not trans_spliced else {-20.0})


def test_a_trans_spliced_join_on_a_circular_record_leaves_the_origin_to_other_features() -> None:
    options = CircularDiagramOptions(config=apply_config_overrides(None, {"canvas.resolve_overlaps": True}))
    request = CircularDiagramRequest(records=(_ends_record("circular", trans_spliced=True),), options=options)
    svg = plan_request(request).build().tostring()
    assert len({round(min(math.hypot(x, y) for x, y in _points(d)), 1) for d in _block_paths(svg)}) == 1


def _part_order(svg: str) -> list[str]:
    """The feature parts (``block`` or ``connector``) of the drawing, in document order."""
    return re.findall(r'data-gbdraw-feature-part="(block|connector)"', svg)


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("display_start", [None, 4001])
def test_a_feature_draws_its_connector_lines_under_its_blocks(tmp_path: Path, mode: str, display_start: int | None) -> None:
    # join(101..200,301..400,501..600): the lines end on the blocks' edges; each block's outline must cover their caps.
    location = SimpleLocation(100, 200, 1) + SimpleLocation(300, 400, 1) + SimpleLocation(500, 600, 1)
    records = read_genbank(_record_file(tmp_path, location))
    draw = draw_linear if mode == "linear" else draw_circular
    displays = None if display_start is None else [RecordDisplayOptions(start_coordinate=display_start)]
    order = _part_order(draw(records, record_displays=displays).to_svg())
    assert order == ["connector", "connector", "block", "block", "block"]
