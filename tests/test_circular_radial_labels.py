from __future__ import annotations

import math
import re
from copy import deepcopy
from pathlib import Path

import pytest

import gbdraw.labels.circular as circular_labels_module
from gbdraw.api.config import apply_config_overrides, load_default_config
from gbdraw.circular import _get_args, run_circular_from_namespace
from gbdraw.exceptions import ValidationError
from gbdraw.labels.circular_leaders import leader_polyline
from gbdraw.labels.fit import crossing_pairs
from gbdraw.labels.circular_radial import place_radial_labels, radial_text_orientation


def _label(index: int, angle_deg: float, *, inner: bool = False, width: float = 80.0) -> dict:
    middle = (angle_deg % 360.0) / 360.0 * 10_000.0
    radius = 380.0 if not inner else 270.0
    return {
        "stable_id": f"feature-{index}",
        "input_order": index,
        "label_text": f"label {index}",
        "middle": middle,
        "preferred_middle": middle,
        "width_px": width,
        "height_px": 12.0,
        "middle_x": radius,
        "middle_y": 0.0,
        "feature_anchor_x": 340.0 if not inner else 300.0,
        "feature_anchor_y": 0.0,
        "feature_inner_radius_px": 290.0,
        "feature_outer_radius_px": 350.0,
        "is_inner": inner,
    }


@pytest.mark.parametrize("angle", [0.0, 90.0, 180.0, 270.0, 359.999999])
@pytest.mark.parametrize("inner", [False, True])
def test_radial_orientation_is_readable_and_stable(angle: float, inner: bool) -> None:
    rotation, anchor, baseline = radial_text_orientation(angle, is_inner=inner)
    assert -90.0 <= rotation <= 90.0
    assert anchor in {"start", "end"}
    assert baseline == "middle"


def test_radial_dense_seam_cluster_is_collision_free_and_deterministic() -> None:
    labels = [_label(index, angle) for index, angle in enumerate((358.0, 359.0, 0.0, 0.0, 1.0, 2.0))]
    first = place_radial_labels(
        labels,
        (),
        total_length=10_000,
        spacing_px=3.0,
        outer_reserved_radius_px=360.0,
    )
    second = place_radial_labels(
        deepcopy(labels),
        (),
        total_length=10_000,
        spacing_px=3.0,
        outer_reserved_radius_px=360.0,
    )

    assert first.collision_free
    assert len(first.labels) == len(labels)
    assert [label["placed_angle_deg"] for label in first.labels] == [
        label["placed_angle_deg"] for label in second.labels
    ]
    assert all(label["placement"] == "radial" for label in first.labels)
    assert all(len(label["oriented_corners"]) == 4 for label in first.labels)


def _label_on_feature(index: int, angle_deg: float, *, inner: bool) -> dict:
    """A label whose leader starts on its own feature, at the feature's angle (no shared anchor)."""

    label = _label(index, angle_deg, inner=inner)
    direction = math.radians(angle_deg) - 0.5 * math.pi
    anchor_radius, elbow_radius = (300.0, 280.0) if inner else (340.0, 380.0)
    label["feature_anchor_x"] = anchor_radius * math.cos(direction)
    label["feature_anchor_y"] = anchor_radius * math.sin(direction)
    label["middle_x"] = elbow_radius * math.cos(direction)
    label["middle_y"] = elbow_radius * math.sin(direction)
    return label


@pytest.mark.parametrize("inner", [False, True])
def test_radial_crowded_cluster_spreads_to_both_sides_of_its_features(inner: bool) -> None:
    """A cluster too crowded for its features' angles spreads around them, not one way only (OV-3648).

    Clockwise packing alone keeps the first label on its feature and pushes every other one clockwise, which slants
    the leaders of the whole cluster; whether it ran alone used to depend on whether its leaders crossed, so a small
    change of leader ends turned a cluster by several degrees.
    """

    angles = [90.0 + 0.25 * index for index in range(9)]
    labels = [_label_on_feature(index, angle, inner=inner) for index, angle in enumerate(angles)]
    result = place_radial_labels(
        () if inner else labels,
        labels if inner else (),
        total_length=10_000,
        spacing_px=3.0,
        outer_reserved_radius_px=360.0,
        inner_reserved_outer_radius_px=100.0,
    )
    assert result.collision_free
    offsets = {
        label["stable_id"]: ((label["placed_angle_deg"] - label["preferred_angle_deg"] + 180.0) % 360.0) - 180.0
        for label in result.labels
    }
    spread = max(offsets.values()) - min(offsets.values())
    assert spread > 1.0, offsets  # the cluster had to spread
    assert min(offsets.values()) < -0.25 * spread and max(offsets.values()) > 0.25 * spread, offsets
    assert abs(sum(offsets.values()) / len(offsets)) < 0.1 * spread, offsets


def test_long_inner_label_reports_required_radius_growth_without_truncation() -> None:
    inner = _label(0, 45.0, inner=True, width=500.0)
    result = place_radial_labels(
        (),
        (inner,),
        total_length=10_000,
        spacing_px=3.0,
        inner_reserved_outer_radius_px=100.0,
    )
    assert result.required_radius_growth_px > 0.0
    assert result.labels[0]["label_text"] == inner["label_text"]
    assert result.labels[0]["width_px"] == 500.0


def test_mode_label_config_defaults_and_circular_placement_validation() -> None:
    old_config = load_default_config()
    old_config["labels"].pop("circular")
    old_config["labels"].pop("linear")
    labels = apply_config_overrides(old_config, None).labels
    assert (labels.circular.scope, labels.circular.placement) == ("none", "horizontal")
    assert (labels.linear.scope, labels.linear.placement, labels.linear.rotation) == (
        "none",
        "auto",
        0.0,
    )
    assert apply_config_overrides(
        None,
        {"labels.circular.placement": "radial"},
    ).labels.circular.placement == "radial"
    with pytest.raises(ValidationError, match="horizontal.*radial"):
        apply_config_overrides(
            None,
            {"labels.circular.placement": "diagonal"},
        )


def test_circular_cli_placement_default_and_radial_parse() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gbk"
    base_args = ["--gbk", str(input_path)]
    assert _get_args(base_args).label_placement == "horizontal"
    assert _get_args([*base_args, "--label_placement", "radial"]).label_placement == "radial"


def test_radial_solver_keeps_500_labels_without_collisions_or_cutoffs() -> None:
    label_count = 500
    total_length = 100_000
    labels = []
    for index in range(label_count):
        angle_deg = 360.0 * index / label_count
        angle_rad = math.radians(angle_deg - 90.0)
        labels.append(
            {
                "stable_id": f"feature-{index}",
                "input_order": index,
                "label_text": f"label {index}",
                "middle": angle_deg / 360.0 * total_length,
                "preferred_middle": angle_deg / 360.0 * total_length,
                "width_px": 22.0,
                "height_px": 8.0,
                "middle_x": 360.0 * math.cos(angle_rad),
                "middle_y": 360.0 * math.sin(angle_rad),
                "feature_anchor_x": 340.0 * math.cos(angle_rad),
                "feature_anchor_y": 340.0 * math.sin(angle_rad),
                "feature_inner_radius_px": 330.0,
                "feature_outer_radius_px": 350.0,
                "is_inner": False,
            }
        )

    result = place_radial_labels(
        labels,
        (),
        total_length=total_length,
        spacing_px=2.0,
        outer_reserved_radius_px=360.0,
    )

    assert len(result.labels) == label_count
    assert result.collision_free


def test_candidate_extraction_runs_once_across_inner_radius_preflight(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    calls = 0
    original = circular_labels_module.build_circular_label_candidates

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(circular_labels_module, "build_circular_label_candidates", counted)
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gbk"
    args = _get_args(
        [
            "--gbk",
            str(input_path),
            "--labels",
            "both",
            "--label_placement",
            "radial",
            "--label_rendering",
            "external_only",
            "--legend",
            "none",
            # The full name is drawn over the GC rows on one line (DF-D25), and the definition rows of its warning
            # draw the figure again; the abbreviation keeps this to the one drawing whose preflight is counted.
            "--species",
            "MjeNMV",
            "--format",
            "svg",
            "--output",
            str(tmp_path / "radial-candidate-cache"),
        ]
    )

    run_circular_from_namespace(args)

    assert calls == 1


def _rotated_tobacco_argv(tmp_path: Path, output: Path) -> list[str]:
    """Docs recipe H-CLI-14: the tobacco plastome drawn from bp 5500, one feature placed outward, radial labels."""

    repository = Path(__file__).parents[1]
    placements = tmp_path / "placements.tsv"
    placements.write_text(
        "record\tfeature_selector\tplacement\tlevel\nNC_001879.2\tprotein_id=NP_054479.1\toutward\t1\n",
        encoding="utf-8",
    )
    priority = tmp_path / "qualifier_priority.tsv"
    priority.write_text("CDS\tgene,old_locus_tag\n", encoding="utf-8")
    return [
        "--gbk", str(repository / "tests/test_inputs/NC_001879.gbk"),
        "-t", "builtin:plastid",
        "-k", "CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin",
        "--track_type", "tuckin",
        "--display_start_coordinate", "5500",
        "--feature_placement_table", str(placements),
        "--feature_overlap_tolerance_bp", "1",
        "--resolve_overlaps",
        "--gc",
        "--no-skew",
        "--labels", "both",
        "--label_placement", "radial",
        "--outer_label_x_radius_offset", "0.9",
        "--outer_label_y_radius_offset", "0.9",
        "--inner_label_x_radius_offset", "0.975",
        "--inner_label_y_radius_offset", "0.975",
        "--qualifier_priority", str(priority),
        "--annotation_table", str(repository / "examples/nicotiana-tabacum-regions.tsv"),
        "--circular_track_slot", "features:features@side=overlay,lane_direction=split",
        "--circular_track_slot",
        "plastome_regions:annotations@set_id=plastome_regions,side=inside,r=0.65,w=20px,inner_gap_px=1,"
        "outer_gap_px=1,show_labels=true,padding_px=1,overflow=compress",
        "--circular_track_slot",
        "gc_content:dinucleotide_content@side=inside,r=0.56,w=0.08,nt=GC,legend_label=GC content",
        "--legend", "none",
        "--format", "svg",
        "--output", str(output),
    ]


# The ribosomal protein cluster under ycf2 in the rotated tobacco plastome (docs recipe H-CLI-14, OV-3648): rpoA,
# rps11, rpl36, rps8, rpl14, rpl16, rps3, rpl22, rps19, rpl2, rpl23 and tRNA-Ile, by feature id (the inverted-repeat
# copies of rpl2, rpl23, rps19 and tRNA-Ile carry other ids).
_TOBACCO_SOUTH_FEATURE_IDS = {
    "fab2c2b1d", "f46c497f8", "f7cae9638", "fbaf469e4", "ff70975fb", "f2126f55f",
    "f7371cc8d", "f33c446c4", "f3d0f7434", "f2c0d3c71", "f4e9aae7d", "fc77059df",
}


def test_rotated_tobacco_plastome_keeps_its_south_inner_labels_over_their_features(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """H-CLI-14's inner labels below ycf2 spread around their features (OV-3648).

    The recipe places one feature outward, so inner leaders start on each feature's drawn outline; that removed a
    leader crossing from the clockwise packing, and the packing alone then turned the whole cluster up to 12 degrees
    one way, slanting every leader. The cluster spreads to both sides of its features instead.
    """

    results = []
    original = circular_labels_module.place_radial_labels

    def captured(*args, **kwargs):
        result = original(*args, **kwargs)
        results.append(result)
        return result

    monkeypatch.setattr(circular_labels_module, "place_radial_labels", captured)
    args = _get_args(_rotated_tobacco_argv(tmp_path, tmp_path / "rotated_placed_chloroplast"))
    run_circular_from_namespace(args)

    cluster = [
        label for label in results[-1].labels
        if label["feature_id"] in _TOBACCO_SOUTH_FEATURE_IDS
    ]
    assert len(cluster) == len(_TOBACCO_SOUTH_FEATURE_IDS)
    offsets = [
        ((float(label["placed_angle_deg"]) - float(label["preferred_angle_deg"]) + 180.0) % 360.0) - 180.0
        for label in cluster
    ]
    # Before the fix: every offset 0 to +12 degrees, mean +5.6. Spread around the features: -6 to +6, mean 0.
    assert min(offsets) < -2.0 and max(offsets) < 8.0, offsets
    assert abs(sum(offsets) / len(offsets)) < 2.0, offsets


def _block_outlines(svg_text: str) -> list[list[tuple[float, float]]]:
    """The drawn feature blocks of a Circular SVG as polygons: its arcs (about the circle centre) in 0.25 degree steps."""

    outlines = []
    for element in re.findall(r'<path [^>]*data-gbdraw-feature-part="block"[^>]*/>', svg_text):
        tokens = re.findall(r"[MLAz]|-?[0-9.]+(?:e-?[0-9]+)?", re.search(r' d="([^"]+)"', element).group(1))
        points: list[tuple[float, float]] = []
        index = 0
        while index < len(tokens):
            command = tokens[index]
            if command in {"M", "L"}:
                points.append((float(tokens[index + 1]), float(tokens[index + 2])))
                index += 3
            elif command == "A":
                sweep = tokens[index + 5] == "1"
                end = (float(tokens[index + 6]), float(tokens[index + 7]))
                start_angle = math.atan2(points[-1][1], points[-1][0])
                turn = (math.atan2(end[1], end[0]) - start_angle) % (2 * math.pi)
                if not sweep:
                    turn -= 2 * math.pi
                radius = float(tokens[index + 1])
                steps = max(1, int(abs(math.degrees(turn)) / 0.25))
                points.extend(
                    (radius * math.cos(start_angle + turn * step / steps), radius * math.sin(start_angle + turn * step / steps))
                    for step in range(1, steps)
                )
                points.append(end)
                index += 8
            else:
                index += 1
        outlines.append(points)
    return outlines


def _inside(point: tuple[float, float], polygon: list[tuple[float, float]]) -> bool:
    x, y = point
    inside = False
    for (x1, y1), (x2, y2) in zip(polygon, polygon[1:] + polygon[:1]):
        if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def _lambda_two_slot_argv(font_size: str, output: Path) -> list[str]:
    return [
        "--gbk", str(Path(__file__).parents[1] / "tests/test_inputs/NC_001416.gb"),
        "--labels", "both",
        "--label_placement", "radial",
        "--label_font_size", font_size,
        "--circular_track_slot", "outer:features@side=outside,lane_direction=outside",
        "--circular_track_slot", "ticks:ticks",
        "--circular_track_slot", "features:features",
        "--format", "svg",
        "--output", str(output),
    ]


def _radial_case_argv(case: str, tmp_path: Path, output: Path) -> list[str]:
    if case == "rotated-tobacco":
        return _rotated_tobacco_argv(tmp_path, output)
    return _lambda_two_slot_argv(case.removeprefix("lambda-").removesuffix("pt"), output)


@pytest.mark.parametrize("case", ["lambda-8.5pt", "lambda-8pt", "rotated-tobacco"])
def test_radial_leaders_run_outside_the_feature_band(
    case: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A radial leader enters no feature block before its end on its feature's outline (OV-3841).

    Leaders are drawn below the features, so a leader that runs through a block is hidden there. In lambda with
    labels on both sides, the inner leaders ran from an elbow on the band's inner edge, at the label's angle, to the
    middle of the lane at the feature's angle, through up to 14 degrees of neighbouring features; with the elbow off
    the band they ran side by side along it, and 8.5 pt no longer fit until growth widened the inner room. In the
    rotated tobacco plastome (docs recipe H-CLI-14) the CDS of the lane above covers the whole outer edge of the
    replication origin NitaCp079: a leader at the label's own angle ran under that CDS.
    """

    results = []
    original = circular_labels_module.place_radial_labels

    def captured(*args, **kwargs):
        result = original(*args, **kwargs)
        results.append(result)
        return result

    monkeypatch.setattr(circular_labels_module, "place_radial_labels", captured)
    output = tmp_path / "radial"
    run_circular_from_namespace(_get_args(_radial_case_argv(case, tmp_path, output)))

    labels = [label for label in results[-1].labels if label.get("placement") == "radial"]
    assert len(labels) > 50
    assert not _hidden_leaders(labels, output.with_suffix(".svg").read_text(encoding="utf-8"))


def _hidden_leaders(labels, svg_text: str) -> dict[tuple[str, str], float]:
    """The radial leaders that run inside a drawn feature block before their end, with the length hidden (px)."""

    outlines = _block_outlines(svg_text)
    assert outlines
    buried = {}
    for label in labels:
        anchor, elbow, contact = leader_polyline(label)
        hidden = 0.0
        for start, end in ((contact, elbow), (elbow, anchor)):
            length = math.dist(start, end)
            steps = max(1, int(length / 0.25))
            for step in range(steps + 1):
                share = step / steps
                point = (start[0] + (end[0] - start[0]) * share, start[1] + (end[1] - start[1]) * share)
                # The end on the feature's outline is the leader's last point.
                if math.dist(point, anchor) > 0.5 and any(_inside(point, outline) for outline in outlines):
                    hidden += length / steps
        if hidden > 0.0:
            buried[(label["feature_id"], label["label_text"])] = round(hidden, 1)
    return buried


# Two overlapping hypothetical proteins of the tobacco plastome (NC_001879.gbk) drawn in one lane as arrowheads that
# point at each other.
TOBACCO_OVERLAPPING_HEADS = ("f8ec4b35c", "fa1b8e657")


def test_leaders_of_overlapping_arrowheads_end_in_the_order_of_their_labels(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """OV-3841: the nearest points of two overlapping arrowheads' slanted edges are their tips, in the opposite order
    of the labels; the leaders crossed at every radius and packing never converged, so the tobacco plastome with
    radial labels outside, drawn before, was refused. Their leaders now end in the order of their labels."""

    results = []
    original = circular_labels_module.place_radial_labels

    def captured(*args, **kwargs):
        result = original(*args, **kwargs)
        results.append(result)
        return result

    monkeypatch.setattr(circular_labels_module, "place_radial_labels", captured)
    run_circular_from_namespace(
        _get_args(
            [
                "--gbk", str(Path(__file__).parents[1] / "tests/test_inputs/NC_001879.gbk"),
                "--labels", "out", "--label_placement", "radial", "--label_font_size", "10",
                "--format", "svg", "--output", str(tmp_path / "tobacco"),
            ]
        )
    )

    by_id = {label["feature_id"]: label for label in results[-1].labels}
    first, second = (by_id[feature_id] for feature_id in TOBACCO_OVERLAPPING_HEADS)
    assert crossing_pairs([leader_polyline(first), leader_polyline(second)]) == []


# Two neighbouring features of one lane in the M16-5 genome (M16-5.gb): "hypothetical protein" and the tip of
# "Large ribosomal subunit protein uL10", which points at it.
M16_5_NEIGHBOURS = ("fd3da9502", "fd1c87198")


@pytest.mark.parametrize("track_type", ["default", "tuckin"])
def test_a_radial_leader_ends_where_it_passes_the_neighbouring_features_of_its_lane(
    track_type: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """OV-3841, M16-5 with radial labels on both sides: a crowded inner label moved far along the circle ended its
    leader at the nearest point of its feature, the tip of an arrowhead in the middle of the lane. That leader came in
    along the band, ran under the neighbouring feature of the same lane, and crossed the neighbour's leader at every
    radius, so growth stalled and the figure, drawn before, was refused. A leader now ends at the nearest point from
    which it passes the neighbouring features."""

    results = []
    original = circular_labels_module.place_radial_labels

    def captured(*args, **kwargs):
        result = original(*args, **kwargs)
        results.append(result)
        return result

    monkeypatch.setattr(circular_labels_module, "place_radial_labels", captured)
    output = tmp_path / "m16-5"
    layout = ["--separate_strands", "--track_type", "tuckin", "--legend", "right"] if track_type == "tuckin" else []
    run_circular_from_namespace(
        _get_args(
            [
                "--gbk", str(Path(__file__).parents[1] / "tests/test_inputs/M16-5.gb"),
                *layout,
                "--labels", "both", "--label_placement", "radial", "--label_font_size", "10",
                "--format", "svg", "--output", str(output),
            ]
        )
    )

    labels = [label for label in results[-1].labels if label.get("placement") == "radial"]
    by_id = {label["feature_id"]: label for label in labels}
    first, second = (by_id[feature_id] for feature_id in M16_5_NEIGHBOURS)
    assert crossing_pairs([leader_polyline(first), leader_polyline(second)]) == []
    assert results[-1].leader_crossing_count == 0
    assert not _hidden_leaders(labels, output.with_suffix(".svg").read_text(encoding="utf-8"))


# Feature id of "AAA family ATPase" in lambda (NC_001416.gb).
LAMBDA_AAA_ATPASE = "f53e02f57"


def test_a_label_too_long_for_its_feature_stays_outside_in_the_two_slot_radial_figure(tmp_path: Path) -> None:
    """OV-3646 and OV-3649: in the Owner's lambda figure (two features slots, radial labels on both sides, 8.5 pt),
    "AAA family ATPase" is longer than its feature can hold, at any position along it: it stays outside."""

    output = tmp_path / "radial"
    run_circular_from_namespace(_get_args(_lambda_two_slot_argv("8.5", output)))

    svg_text = output.with_suffix(".svg").read_text(encoding="utf-8")
    texts = re.findall(
        rf'<text [^>]*data-label-feature-id="{LAMBDA_AAA_ATPASE}"[^>]*>(<textPath[^>]*>)?([^<]*)', svg_text
    )
    assert texts == [("", "AAA family ATPase")]


def _slot_markup(svg_text: str, slot_id: str) -> str:
    """The markup of one track slot's group, up to the next top-level group."""

    start = svg_text.index(f'data-gbdraw-slot-id="{slot_id}"')
    end = svg_text.find('<g data-gbdraw-composition-role="primary"', start)
    return svg_text[start : end if end >= 0 else len(svg_text)]


def _path_point_radii(markup: str) -> list[float]:
    """The radii of the points the paths in ``markup`` move, draw lines and arcs to (about the circle centre)."""

    radii = []
    for path_data in re.findall(r'<path d="([^"]+)"', markup):
        command = ""
        numbers: list[float] = []
        for token in re.findall(r"[MLAZmlaz]|-?\d+(?:\.\d+)?(?:e-?\d+)?", path_data) + ["Z"]:
            if not token[0].isalpha():
                numbers.append(float(token))
                continue
            size = 7 if command == "A" else 2
            radii.extend(
                math.hypot(numbers[index + size - 2], numbers[index + size - 1])
                for index in range(0, len(numbers) - size + 1, size)
            )
            command, numbers = token.upper(), []
    return radii


def test_radial_growth_holds_only_the_tracks_inside_the_inner_label_lanes(tmp_path: Path) -> None:
    """Growth for radial inner labels keeps the tracks inside the inner labels' lanes and moves the rest (OV-3841).

    In the lambda repro a scale stands between the two features slots. Held at its radius, it held the inner lanes
    in place too, so growing the circle never widened the inner labels' room and the figure was refused. Now the
    figure draws: the scale stays outside the inner lanes, and the inner labels lie between those lanes and the
    dinucleotide track inside them.
    """

    repository = Path(__file__).parents[1]
    output = tmp_path / "lambda_growth"
    args = _get_args(
        [
            "--gbk", str(repository / "tests/test_inputs/NC_001416.gb"),
            "--labels", "both",
            "--label_placement", "radial",
            "--label_font_size", "8.5",
            "--circular_track_slot", "outer:features@side=outside,lane_direction=outside",
            "--circular_track_slot", "ticks:ticks",
            "--circular_track_slot", "features:features",
            "--circular_track_slot", "gc_content:dinucleotide_content@side=inside,r=0.42,w=0.04",
            "--legend", "none",
            "--format", "svg",
            "--output", str(output),
        ]
    )
    run_circular_from_namespace(args)  # growth widened the inner labels' room: no LabelsDoNotFitError

    svg_text = output.with_suffix(".svg").read_text(encoding="utf-8")
    inner_lanes = [math.hypot(*point) for outline in _block_outlines(_slot_markup(svg_text, "features")) for point in outline]
    ticks = _path_point_radii(_slot_markup(svg_text, "ticks"))
    dinucleotide = _path_point_radii(_slot_markup(svg_text, "gc_content"))
    labels_text = svg_text[svg_text.index('id="label_text"') :]
    inner_label_radii = [
        radius
        for x_coord, y_coord in re.findall(r'<text [^>]*data-label-feature-id="[^"]+"[^>]*x="([^"]+)" y="([^"]+)"', labels_text)
        if (radius := math.hypot(float(x_coord), float(y_coord))) < min(inner_lanes)
    ]
    assert inner_lanes and ticks and dinucleotide and len(inner_label_radii) > 20
    assert min(ticks) >= max(inner_lanes) - 1e-6  # outside the inner lanes: it moves with the radius
    assert max(dinucleotide) < min(inner_label_radii)  # inside the inner labels: it keeps its radius


def test_label_placement_reads_feature_blocks_where_they_are_drawn(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The intervals labels keep clear of lie where the blocks are drawn, also from a display start (OV-3861).

    ``feature_block_intervals`` read the record's location; the rotated tobacco plastome draws its blocks 5,500 bp
    on, so the clearance of labels and the ends of radial leaders looked at other blocks.
    """

    import gbdraw.diagrams.circular.assemble as assemble_module
    from gbdraw.layout.circular import feature_block_intervals
    from gbdraw.svg.circular_arcs import circle_point, position_angle

    captured = []
    original = assemble_module.prepare_label_list

    def capture(*args, **kwargs):
        captured.append(args[0])
        return original(*args, **kwargs)

    monkeypatch.setattr(assemble_module, "prepare_label_list", capture)
    output = tmp_path / "rotated"
    run_circular_from_namespace(_get_args(_rotated_tobacco_argv(tmp_path, output)))

    svg_text = output.with_suffix(".svg").read_text(encoding="utf-8")
    drawn: dict[str, list[list[tuple[float, float]]]] = {}
    for element in re.findall(r'<path [^>]*data-gbdraw-feature-part="block"[^>]*/>', svg_text):
        feature_id = re.search(r'data-gbdraw-feature-id="([^"]+)"', element).group(1)
        drawn.setdefault(feature_id, []).extend(_block_outlines(element))
    total_length = 155943
    misplaced = []
    for feature in captured[-1].values():
        for _block, start, end in feature_block_intervals(feature, total_length):
            outlines = drawn.get(str(feature.feature_hash), [])
            middle_angle = position_angle((start + end) / 2.0, total_length)
            if not any(
                _inside(circle_point(0.5 * (min(radii) + max(radii)), middle_angle), outline)
                for outline in outlines
                for radii in ([math.hypot(x, y) for x, y in outline],)
            ):
                misplaced.append((feature.feature_hash, start, end))
    assert captured and drawn
    assert not misplaced, misplaced[:5]


def test_radial_leader_ends_leave_out_positions_under_another_lane() -> None:
    """A radial leader ends where no feature of another lane lies between its feature and the labels (OV-3841).

    An outer leader that ended under the feature of the lane above ran under it. When such features cover the whole
    labelled segment, every point stays and the outlines of the other features go to the placer, grown by the
    leader's clearance.
    """

    from gbdraw.features.objects import FeatureLocationPart, FeatureObject
    from gbdraw.labels.circular import _radial_leader_ends
    from gbdraw.layout.circular import FeatureRadiusBlock

    def feature(name: str) -> FeatureObject:
        return FeatureObject(
            name, [FeatureLocationPart("block", "1", "positive", 1000, 2000, True)], False, "#000000", "", "", None, "CDS", {}
        )

    own = feature("own")
    radii = (100.0, 110.0, 120.0)
    own_lane = FeatureRadiusBlock(own, 1000.0, 2000.0, 100.0, 120.0)
    lower_lane = FeatureRadiusBlock(feature("lower"), 0.0, 10_000.0, 70.0, 90.0)

    def ends(blocks, *, outer=True, clearance=0.0):
        return _radial_leader_ends(
            own, (1000.0, 2000.0), 10_000, 0.0, radii, 1.0, outer=outer, blocks=blocks, clearance_px=clearance,
            block_outline=lambda block, _position: (block.inner_px, block.outer_px),
        )

    def angles(blocks, **kwargs):
        points = ends(blocks, **kwargs).points
        outer = kwargs.get("outer", True)
        assert all(math.hypot(*point) == pytest.approx(120.0 if outer else 100.0) for point in points)
        return [(math.degrees(math.atan2(y, x)) + 90.0) % 360.0 for x, y in points]

    free = angles([own_lane, lower_lane])
    assert len(free) == 33 and min(free) == pytest.approx(36.0) and max(free) == pytest.approx(72.0)
    upper_half = FeatureRadiusBlock(feature("upper"), 1000.0, 1500.0, 125.0, 145.0)
    half = angles([own_lane, upper_half, lower_lane])
    assert min(half) >= 54.0 - 1e-9 and max(half) == pytest.approx(72.0)  # 1500 bp is 54 degrees
    # Covered throughout: nothing to choose from, and the covering outline goes to the placer.
    covering = FeatureRadiusBlock(feature("cover"), 0.0, 10_000.0, 125.0, 145.0)
    covered = ends([own_lane, covering, lower_lane], clearance=2.0)
    assert len(covered.points) == 33 and len(covered.obstacles) == 1
    assert min(math.hypot(*point) for point in covered.obstacles[0]) == pytest.approx(123.0)
    # An inner label's leader comes from inside: the lane above never blocks it, a lane below over half does.
    assert len(angles([own_lane, upper_half], outer=False)) == 33
    lower_half = FeatureRadiusBlock(feature("lower-half"), 1000.0, 1500.0, 70.0, 90.0)
    assert min(angles([own_lane, lower_half], outer=False)) >= 54.0 - 1e-9
