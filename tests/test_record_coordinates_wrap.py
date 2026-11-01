"""One record drawn over several rows: the display transform owns source -> rows."""

from __future__ import annotations

from dataclasses import replace
from math import ceil
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import pytest

from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    GenBankInputSource,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
)
from gbdraw.io.regions import RegionSpec, parse_region_spec
from gbdraw.exceptions import ValidationError
from gbdraw.features.coordinates import get_exon_and_intron_coordinates, project_feature_parts
from gbdraw.layout.record_coordinates import (
    DisplayFragment,
    RecordDisplayTransform,
    SeriesPoint,
    SourceInterval,
)
from gbdraw.layout.record_labels import RulerTick, source_ruler_ticks


LAMBDA_LENGTH = 48_502
LENGTH = 1000
ROW = 300
# name: (transform, displayed bases). A crop keeps the complete source length; a rotation needs a circular record.
AXES = {
    "plain": (RecordDisplayTransform(LENGTH, 1, 1, None, False), LENGTH),
    "rotated": (RecordDisplayTransform(LENGTH, 1, 1, 401, True), LENGTH),
    "reversed": (RecordDisplayTransform(LENGTH, LENGTH, -1, None, False), LENGTH),
    "reversed_rotated": (RecordDisplayTransform(LENGTH, LENGTH, -1, 401, True), LENGTH),
    "cropped": (RecordDisplayTransform(5000, 1001, 1, None, False), LENGTH),
    "cropped_reversed": (RecordDisplayTransform(5000, 2000, -1, None, False), LENGTH),
}


def wrapped(name: str, row: int = ROW) -> RecordDisplayTransform:
    transform, displayed = AXES[name]
    return transform.with_wrap(row, displayed)


@pytest.mark.parametrize("name", AXES)
def test_axis_length_and_row_count(name: str) -> None:
    transform = wrapped(name)
    assert transform.projects is True
    assert transform.axis_length == ROW
    assert transform.row_count == ceil(LENGTH / ROW) == 4
    assert wrapped(name, 2000).axis_length == LENGTH
    assert wrapped(name, 2000).row_count == 1


def test_an_unwrapped_transform_projects_only_when_rotated() -> None:
    assert AXES["plain"][0].projects is False
    assert AXES["cropped"][0].projects is False
    assert AXES["rotated"][0].projects is True
    assert AXES["plain"][0].wrap_length is None
    assert AXES["plain"][0].axis_length == LENGTH
    assert AXES["plain"][0].row_count == 1


@pytest.mark.parametrize("bad", [0, -1, True, 1.5, "10"])
def test_wrap_length_is_validated(bad: object) -> None:
    with pytest.raises(ValidationError, match="Wrap length"):
        RecordDisplayTransform(LENGTH, 1, 1, None, False, wrap_length=bad)  # type: ignore[arg-type]


@pytest.mark.parametrize("bad", [0, LENGTH + 1])
def test_display_length_is_validated(bad: int) -> None:
    with pytest.raises(ValidationError, match="Display length"):
        RecordDisplayTransform(LENGTH, 1, 1, None, False, display_length=bad)


def _whole_axis(transform: RecordDisplayTransform, displayed: int) -> tuple[DisplayFragment, ...]:
    # A feature on every displayed base, in the record's local frame.
    fragments = transform.project_local_parts((SourceInterval(0, displayed),))
    assert isinstance(fragments, tuple)
    return fragments  # type: ignore[return-value]


@pytest.mark.parametrize("name", AXES)
def test_rows_tile_the_axis_exactly_once(name: str) -> None:
    transform, displayed = AXES[name]
    fragments = _whole_axis(wrapped(name), displayed)
    covered: dict[tuple[int, int], int] = {}
    for fragment in fragments:
        assert 0 <= fragment.display_start < fragment.display_end <= ROW
        for offset in range(fragment.display_end - fragment.display_start):
            cell = (fragment.row, fragment.display_start + offset)
            assert cell not in covered
            # Going right along the row reads the source as the unwrapped axis reads it.
            base = (
                fragment.source_start + offset + 1
                if transform.source_step == 1
                else fragment.source_end - offset
            )
            covered[cell] = base
    assert len(covered) == displayed
    for (row, x), base in covered.items():
        assert transform.source_base_to_display_index(base) == row * ROW + x
    assert {row for row, _x in covered} == {0, 1, 2, 3}


@pytest.mark.parametrize("name", AXES)
def test_local_bounds_follow_the_display_bounds(name: str) -> None:
    for fragment in _whole_axis(wrapped(name), AXES[name][1]):
        assert fragment.local_end - fragment.local_start == fragment.display_end - fragment.display_start
        assert fragment.source_end - fragment.source_start == fragment.display_end - fragment.display_start


@pytest.mark.parametrize("name", AXES)
def test_a_cut_edge_is_artificial_and_open_and_biological_ends_stay_true(name: str) -> None:
    transform, displayed = AXES[name]
    fragments = _whole_axis(wrapped(name), displayed)
    # A row edge that is also the local origin of a rotated axis is an end of this interval, not a cut.
    for fragment in fragments:
        if fragment.display_start == 0 and fragment.row > 0 and fragment.local_start != 0:
            assert fragment.artificial_start and fragment.open_start
        if fragment.display_end == ROW and fragment.local_end != displayed:
            assert fragment.artificial_end and fragment.open_end
    assert sum(fragment.biological_start for fragment in fragments) == 1
    assert sum(fragment.biological_end for fragment in fragments) == 1
    assert fragments[0].biological_start
    assert fragments[-1].biological_end


@pytest.mark.parametrize(
    ("step", "strand", "low_row_is_first"),
    [(1, 1, True), (1, -1, False), (-1, 1, False), (-1, -1, True)],
)
def test_traversal_order_and_terminals_follow_the_strand(step, strand, low_row_is_first) -> None:
    # Source 250..350 of a 1000 bp record crosses a row boundary on either orientation.
    transform = RecordDisplayTransform(LENGTH, 1 if step == 1 else LENGTH, step, None, False).with_wrap(ROW, LENGTH)
    start, end = (250, 350) if step == 1 else (650, 750)
    first, second = transform.project_interval(SourceInterval(start, end, strand))
    assert (first.row < second.row) is low_row_is_first
    assert first.biological_start and not first.biological_end
    assert second.biological_end and not second.biological_start
    assert first.orientation == second.orientation == strand * step
    low, high = sorted((first, second), key=lambda fragment: fragment.row)
    assert low.artificial_end and low.open_end and not low.artificial_start and not low.open_start
    assert high.artificial_start and high.open_start and not high.artificial_end and not high.open_end
    assert (low.display_start, low.display_end) == (250, 300)
    assert (high.display_start, high.display_end) == (0, 50)


def test_an_interval_ending_on_a_row_boundary_is_not_cut() -> None:
    fragment, = wrapped("plain").project_interval(SourceInterval(100, 300))
    assert fragment.row == 0
    assert (fragment.artificial_start, fragment.artificial_end) == (False, False)
    assert (fragment.open_start, fragment.open_end) == (False, False)
    fragment, = wrapped("plain").project_interval(SourceInterval(300, 400))
    assert fragment.row == 1
    assert (fragment.display_start, fragment.display_end) == (0, 100)


def test_an_interval_over_the_origin_and_a_row_boundary() -> None:
    # Display 0 is base 401; source 900..1000 + 0..100 lies at display 500..700, one row boundary at 600.
    transform = wrapped("rotated")
    parts = transform.project_parts((SourceInterval(900, 1000, 1, 0), SourceInterval(0, 100, 1, 1)))
    assert [(f.row, f.display_start, f.display_end, f.part_index) for f in parts] == [
        (1, 200, 300, 0),
        (2, 0, 100, 1),
    ]
    # The row boundary falls exactly on the origin: an existing fragment edge, not a new cut.
    assert not parts[0].artificial_end and not parts[1].artificial_start
    # An anchor seam stays artificial on a wrapped rotated record.
    seam = transform.project_interval(SourceInterval(350, 450))
    assert [(f.row, f.display_start, f.display_end, f.artificial_start, f.artificial_end) for f in seam] == [
        (3, 50, 100, False, True),
        (0, 0, 50, True, False),
    ]


def test_partial_ends_stay_on_the_true_outer_ends() -> None:
    fragments = wrapped("plain").project_interval(
        SourceInterval(250, 650, 1, 0, biological_start=False, biological_end=True, partial_start=True)
    )
    assert [f.row for f in fragments] == [0, 1, 2]
    assert fragments[0].open_start and not fragments[0].biological_start
    assert fragments[1].open_start and fragments[1].open_end
    assert fragments[2].open_end is False and fragments[2].biological_end


@pytest.mark.parametrize("name", AXES)
def test_one_row_wrapped_fragments_equal_the_unwrapped_geometry(name: str) -> None:
    transform, displayed = AXES[name]
    one_row = transform.with_wrap(displayed, displayed)
    intervals = (SourceInterval(100, 250), SourceInterval(300, 340, -1))
    for interval in intervals:
        got = one_row.project_interval(interval) if name not in ("cropped", "cropped_reversed") else None
        if got is not None:
            assert got == transform.project_interval(interval)
            assert all(fragment.row == 0 for fragment in got)
    fragments = one_row.project_local_parts(intervals)
    assert [(f.display_start, f.display_end) for f in fragments] == [
        (
            transform.local_boundary_to_display_offset(interval.start),
            transform.local_boundary_to_display_offset(interval.start) + interval.end - interval.start,
        )
        for interval in intervals
    ] or transform.start_coordinate is not None


@pytest.mark.parametrize("strand", [1, -1])
def test_a_join_over_a_row_boundary_has_one_arrow_block_and_a_cut_connector(strand: int) -> None:
    exons = [SimpleLocation(100, 200, strand), SimpleLocation(400, 500, strand)]
    if strand == -1:
        exons.reverse()
    parts = project_feature_parts(exons, wrapped("plain"))
    lines = [part for part in parts if part.kind == "line"]
    blocks = [part for part in parts if part.kind == "block"]
    assert [(p.fragment.row, p.fragment.display_start, p.fragment.display_end) for p in sorted(lines, key=lambda p: p.fragment.row)] == [
        (0, 200, 300),
        (1, 0, 100),
    ]
    assert sum(block.location.is_last for block in blocks) == 1
    assert all(line.location.is_last is False for line in lines)


def _connectors(exons: list[SimpleLocation], transform: RecordDisplayTransform) -> list[tuple[int, int, int]]:
    return [
        (part.fragment.row, part.fragment.display_start, part.fragment.display_end)
        for part in project_feature_parts(exons, transform)
        if part.kind == "line"
    ]


_WRAPPED = RecordDisplayTransform(30_000, 1, 1, None, True).with_wrap(10_000, 30_000)


@pytest.mark.parametrize(
    "exons",
    [
        # Over the origin of a circular record: the unwrapped figure draws two blocks and no connection.
        [(28_000, 28_500, 1), (1_500, 2_000, 1)],
        # Out of order without trans_splicing: the same.
        [(12_000, 13_000, 1), (5_000, 6_000, 1)],
        [(5_000, 6_000, -1), (12_000, 13_000, -1)],
    ],
    ids=["origin", "out_of_order", "out_of_order_minus"],
)
def test_a_wrapped_record_adds_no_connection_its_unwrapped_figure_does_not_draw(exons) -> None:
    locations = [SimpleLocation(start, end, strand=strand) for start, end, strand in exons]
    assert not [part for part in get_exon_and_intron_coordinates(locations, 30_000) if part.kind == "line"]
    assert _connectors(locations, _WRAPPED) == []


@pytest.mark.parametrize(
    ("exons", "expected"),
    [
        ([(5_000, 6_000, 1), (12_000, 13_000, 1)], [(0, 6_000, 10_000), (1, 0, 2_000)]),
        ([(12_000, 13_000, -1), (5_000, 6_000, -1)], [(1, 0, 2_000), (0, 6_000, 10_000)]),
        # Strands that differ: the unwrapped figure connects them by the first part's strand.
        ([(5_000, 6_000, 1), (7_000, 8_000, -1)], [(0, 6_000, 7_000)]),
        ([(7_000, 8_000, -1), (5_000, 6_000, 1)], [(0, 6_000, 7_000)]),
    ],
    ids=["plus", "minus", "plus_then_minus", "minus_then_plus"],
)
def test_a_wrapped_record_drops_no_connection_its_unwrapped_figure_draws(exons, expected) -> None:
    locations = [SimpleLocation(start, end, strand=strand) for start, end, strand in exons]
    assert [part for part in get_exon_and_intron_coordinates(locations, 30_000) if part.kind == "line"]
    assert sorted(_connectors(locations, _WRAPPED)) == sorted(expected)


def test_a_rotated_wrapped_record_keeps_its_connection_over_the_origin() -> None:
    rotated = RecordDisplayTransform(30_000, 1, 1, 20_000, True).with_wrap(10_000, 30_000)
    locations = [SimpleLocation(28_000, 28_500, strand=1), SimpleLocation(1_500, 2_000, strand=1)]
    assert _connectors(locations, rotated)


def test_a_block_cut_at_a_row_has_one_arrowhead_on_its_true_end() -> None:
    parts = project_feature_parts([SimpleLocation(250, 350, 1)], wrapped("plain"))
    assert [part.location.is_last for part in parts] == [False, True]
    assert parts[0].location.open_end and parts[1].location.open_start


def test_a_partial_source_feature_keeps_its_open_end_and_gets_no_arrowhead() -> None:
    parts = project_feature_parts(
        [SimpleLocation(250, AfterPosition(350), 1)], wrapped("plain")
    )
    assert [part.location.is_last for part in parts] == [False, False]
    assert parts[1].location.open_end
    parts = project_feature_parts(
        [SimpleLocation(BeforePosition(250), 350, 1)], wrapped("plain")
    )
    assert parts[0].location.open_start and parts[0].location.open_end
    assert parts[1].location.is_last is True


def _whole(points) -> list[tuple[float, float]]:
    return [(point.position, point.value) for point in points]


def _concatenate(segments, row: int = ROW) -> list[tuple[float, float, int | None]]:
    rows = []
    for segment in segments:
        rows.extend(
            (segment.row * row + point.position, point.value, index)
            for point, index in zip(segment.points, segment.source_indices, strict=True)
        )
    return rows


SAMPLES = tuple(SeriesPoint(float(position), position / 10.0) for position in range(0, LENGTH, 70))


@pytest.mark.parametrize("name", AXES)
def test_series_per_row_concatenate_to_the_whole_series(name: str) -> None:
    transform, displayed = AXES[name]
    whole = (
        transform.project_local_series(SAMPLES)
        if transform.projects
        else ()
    )
    segments = wrapped(name).project_local_series(SAMPLES)
    assert [segment.row for segment in segments] == [0, 1, 2, 3]
    joined = _concatenate(segments)
    # Only the row edges (index None) are new; every sample stays, once, in display order.
    originals = [item for item in joined if item[2] is not None]
    if whole:
        assert [(p.position, p.value) for p in whole[0].points] == [
            (position, value) for position, value, _ in originals
        ] or any(index is None for index in whole[0].source_indices)
    else:
        assert [(position, value) for position, value, _ in originals] == _whole(SAMPLES)
    positions = [position for position, _value, _index in joined]
    assert positions == sorted(positions)
    for segment in segments:
        assert all(0 <= point.position <= ROW for point in segment.points)


@pytest.mark.parametrize("name", ["plain", "reversed", "cropped", "cropped_reversed"])
def test_a_row_seam_interpolates_one_shared_edge(name: str) -> None:
    segments = wrapped(name).project_local_series(SAMPLES)
    for left, right in zip(segments, segments[1:], strict=False):
        assert left.points[-1].position == ROW
        assert right.points[0].position == 0
        assert left.points[-1].value == right.points[0].value
        assert left.source_indices[-1] is None and right.source_indices[0] is None
    # Row 0 ends at 300 between the samples at 280 and 350, on the line between them.
    assert segments[0].points[-1].value == pytest.approx(30.0)


def test_a_sample_on_a_row_seam_is_duplicated_not_interpolated() -> None:
    samples = (SeriesPoint(0, 1.0), SeriesPoint(300, 2.0), SeriesPoint(600, 3.0))
    first, second = wrapped("plain").project_local_series(samples)
    assert _whole(first.points) == [(0, 1.0), (ROW, 2.0)]
    assert first.source_indices == (0, 1)
    assert _whole(second.points) == [(0, 2.0), (ROW, 3.0)]
    # The last sample closes its row: no row follows with one point.
    assert second.source_indices == (1, 2)


def test_rows_after_the_last_sample_have_no_series() -> None:
    samples = (SeriesPoint(0, 1.0), SeriesPoint(100, 2.0))
    segments = wrapped("plain").project_local_series(samples)
    assert [segment.row for segment in segments] == [0]


def test_a_rotated_series_opens_at_the_display_start_then_at_the_rows() -> None:
    segments = wrapped("rotated").project_series(SAMPLES)
    assert [segment.row for segment in segments] == [0, 1, 2, 3]
    assert segments[0].points[0].position == 0
    # The display ends at 1000, row-local 100 of the last row, on the same value it started with.
    assert segments[-1].points[-1].position == 100
    assert segments[-1].points[-1].value == segments[0].points[0].value


def test_a_series_ending_exactly_on_the_last_row_edge_opens_no_empty_row() -> None:
    transform = RecordDisplayTransform(900, 1, 1, 401, True).with_wrap(ROW, 900)
    segments = transform.project_series(tuple(SeriesPoint(float(p), 1.0) for p in range(0, 900, 70)))
    assert [segment.row for segment in segments] == [0, 1, 2]
    assert segments[-1].points[-1].position == ROW


def test_one_row_wrapped_series_equals_the_rotated_series() -> None:
    transform, displayed = AXES["rotated"]
    assert transform.with_wrap(displayed, displayed).project_series(SAMPLES) == transform.project_series(SAMPLES)


def test_unwrapped_series_keep_their_paths() -> None:
    assert AXES["plain"][0].project_local_series(SAMPLES) is SAMPLES
    single = AXES["plain"][0].project_series(SAMPLES)
    assert len(single) == 1 and single[0].row == 0 and single[0].seam_sampled is None


def _old_affine_ticks(start: int, end: int, interval: int) -> list[int]:
    # The unrotated axis ruler of the coord map: multiples strictly inside the displayed range.
    low, high = min(start, end), max(start, end)
    return [tick for tick in range(interval, high + interval, interval) if low < tick < high]


@pytest.mark.parametrize(
    ("name", "start", "end"),
    [("plain", 1, 1000), ("reversed", 1000, 1), ("cropped", 1001, 2000), ("cropped_reversed", 2000, 1001)],
)
def test_ticks_of_an_affine_axis_equal_the_coord_map_ticks(name: str, start: int, end: int) -> None:
    transform, displayed = AXES[name]
    expected = _old_affine_ticks(start, end, 100)
    assert sorted(tick.base for tick in source_ruler_ticks(transform.with_wrap(displayed, displayed), 100)) == expected
    assert sorted(tick.base for tick in source_ruler_ticks(wrapped(name), 100)) == expected


@pytest.mark.parametrize("name", AXES)
def test_ticks_land_on_their_row(name: str) -> None:
    transform, _displayed = AXES[name]
    ticks = source_ruler_ticks(wrapped(name), 100)
    assert ticks
    for tick in ticks:
        assert isinstance(tick, RulerTick)
        assert 0 <= tick.offset < ROW
        assert tick.row * ROW + tick.offset == transform.source_base_to_display_index(tick.base)


def test_a_rotated_ruler_keeps_every_base() -> None:
    transform = AXES["rotated"][0]
    assert [(tick.base, tick.row, tick.offset) for tick in source_ruler_ticks(transform, 400)] == [
        (400, 0, 999),
        (800, 0, 399),
    ]
    assert [(tick.base, tick.row, tick.offset) for tick in source_ruler_ticks(wrapped("rotated"), 400)] == [
        (400, 3, 99),
        (800, 1, 99),
    ]


LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
CROP = parse_region_spec("JOIN1:2000-9000")
_NUMBER = re.compile(r'\b(x|x1|x2)="([^"]+)"')


def _join_record(path: Path) -> Path:
    # A spliced feature on each strand, and one plain block.
    record = SeqRecord(Seq("ACGTAACCGG" * 1200), id="JOIN1", name="JOIN1", description="spliced")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"
    for parts, name in (
        ([SimpleLocation(500, 900, 1), SimpleLocation(1500, 1800, 1)], "plus"),
        ([SimpleLocation(5000, 5300, -1), SimpleLocation(4000, 4400, -1)], "minus"),
        ([SimpleLocation(8000, 8400, 1)], "single"),
    ):
        location = CompoundLocation(parts) if len(parts) > 1 else parts[0]
        record.features.append(SeqFeature(location, type="CDS", qualifiers={"gene": [name]}))
    SeqIO.write(record, path, "genbank")
    return path


def _without_row_groups(svg: str) -> str:
    """The SVG with a row group's children in its place, so the one row compares with the unwrapped drawing."""

    ET.register_namespace("", "http://www.w3.org/2000/svg")
    root = ET.fromstring(svg)

    def hoist(parent: ET.Element) -> None:
        for index, child in reversed(list(enumerate(parent))):
            hoist(child)
            if "data-gbdraw-wrap-row" in child.attrib:
                assert "transform" not in child.attrib
                parent.remove(child)
                for offset, grandchild in enumerate(child):
                    parent.insert(index + offset, grandchild)

    hoist(root)
    return ET.tostring(root, encoding="unicode")


def _svg(
    *,
    source: Path = LAMBDA,
    display: RecordDisplayOptions | None = None,
    presentation: RecordPresentation | None = None,
    region: RegionSpec | None = None,
    one_row: bool,
    overrides: dict[str, object],
) -> str:
    plan = plan_request(
        LinearDiagramRequest(
            records=(
                RecordInput(
                    GenBankInputSource(str(source)),
                    region=region,
                    display=display or RecordDisplayOptions(),
                    presentation=presentation or RecordPresentation(),
                ),
            ),
            options=LinearDiagramOptions(config_overrides=overrides),
        )
    )
    if one_row:
        length = len(plan.records[0].seq)
        plan = replace(plan, transforms=(plan.transforms[0].with_wrap(length, length),))
    return _without_row_groups(plan.build().drawing.tostring())


@pytest.mark.parametrize(
    "display",
    [RecordDisplayOptions(), RecordDisplayOptions(is_circular=True, start_coordinate=1000)],
    ids=["plain", "rotated"],
)
@pytest.mark.parametrize(
    "overrides",
    [
        {},
        {"canvas.show_gc": True, "canvas.show_skew": True, "labels.linear.scope": "all"},
    ],
    ids=["features", "gc_skew_labels"],
)
def test_one_row_wrapped_svg_equals_the_unwrapped_svg(display, overrides) -> None:
    assert _svg(display=display, one_row=True, overrides=overrides) == _svg(
        display=display, one_row=False, overrides=overrides
    )


SPLICED_DISPLAYS = {
    "plain": {},
    "rotated": {"display": RecordDisplayOptions(is_circular=True, start_coordinate=1000)},
    "reversed": {"presentation": RecordPresentation(reverse_complement=True)},
    "cropped": {"region": CROP},
    "cropped_reversed": {"region": CROP, "presentation": RecordPresentation(reverse_complement=True)},
}


@pytest.mark.parametrize("fields", SPLICED_DISPLAYS.values(), ids=SPLICED_DISPLAYS)
def test_one_row_wrapped_svg_of_spliced_features_differs_only_in_connector_ends(tmp_path: Path, fields) -> None:
    # A wrapped record draws connectors as a rotated record does: over the whole gap, with butt
    # ends, since a row boundary may cut one. Blocks, labels, and tracks are byte-equal.
    source = _join_record(tmp_path / "join.gb")
    overrides = {"canvas.show_gc": True, "labels.linear.scope": "all"}
    wrapped_elements = _svg(source=source, one_row=True, overrides=overrides, **fields).split("<")
    plain_elements = _svg(source=source, one_row=False, overrides=overrides, **fields).split("<")
    assert len(wrapped_elements) == len(plain_elements)
    rotated = "display" in fields
    for got, plain in zip(wrapped_elements, plain_elements, strict=True):
        if got == plain:
            continue
        assert not rotated
        assert 'data-gbdraw-feature-part="connector"' in got and 'data-gbdraw-feature-part="connector"' in plain
        assert 'stroke-linecap="butt"' in got and 'stroke-linecap="round"' in plain
        same_but_ends = re.compile(r' d="[^"]*"| stroke-linecap="[a-z]+"')
        assert same_but_ends.sub("", got) == same_but_ends.sub("", plain)


def test_one_row_wrapped_ruler_ticks_sit_on_their_base_within_a_base() -> None:
    # An unprojected ruler spreads ticks over L - 1 intervals; a row places them as features, by index / L.
    overrides = {
        "objects.scale.style": "ruler",
        "canvas.linear.track_layout": "above",
        "canvas.linear.ruler_on_axis": True,
        "objects.scale.interval": 10_000,
    }
    got_elements = _svg(one_row=True, overrides=overrides).split("<")
    plain_elements = _svg(one_row=False, overrides=overrides).split("<")
    assert len(got_elements) == len(plain_elements)
    differing = [pair for pair in zip(got_elements, plain_elements, strict=True) if pair[0] != pair[1]]
    assert differing
    for got, plain in differing:
        assert got.startswith(("line ", "text "))
        assert _NUMBER.sub("", got) == _NUMBER.sub("", plain)
        for (_name, got_x), (_other, plain_x) in zip(_NUMBER.findall(got), _NUMBER.findall(plain), strict=True):
            assert float(got_x) == pytest.approx(float(plain_x), abs=2000 / LAMBDA_LENGTH)
