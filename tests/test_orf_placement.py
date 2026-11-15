"""Predicted ORFs beside annotated features (Owner D-131, 2026-11-15).

H1=A: on a record with features of its own, predicted ORFs take their own lanes, outside the circle (Circular) and
below the line (Linear); the annotated features stay where they are drawn without prediction, and predicted ORFs
that overlap each other take further lanes. The Legend keeps the annotated row "CDS" and adds "Predicted ORF".
H2=A: LOSATP reads the proteins of predicted ORFs, also when it counts a subject genome's proteins.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
from Bio import SeqIO

from gbdraw.analysis.orf_prediction import (
    PREDICTED_ORF_CAPTION,
    PREDICTED_ORF_COLOR,
    OrfPredictionParameters,
    predict_record_orf_features,
)
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    FastaInputSource,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.diagrams.circular.radial_layout import build_circular_feature_layout
from gbdraw.features.factory import create_feature_layers
from gbdraw.features.placement import FeaturePlacementSlot
from gbdraw.features.tracks import feature_overlap_bp, get_feature_ends
from tests.test_genome_search_space import DE3, LAMBDA, _fake_losat, _search_spaces
from tests.utils.cli_in_process import run_gbdraw

INPUTS = Path(__file__).resolve().parent / "test_inputs"
PLASMID = INPUTS / "NZ_CP046414.1.gbk"
DRAWN = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region", "regulatory", "rep_origin", "misc_feature"]


def _fasta(tmp_path: Path, path: Path) -> Path:
    record = SeqIO.read(path, "genbank")
    fasta = tmp_path / f"{path.stem}.fasta"
    fasta.write_text(f">{record.id}\n{record.seq}\n")
    return fasta


def _plan(source, mode: str, prediction: OrfPredictionParameters | None):
    request, options = (
        (CircularDiagramRequest, CircularDiagramOptions) if mode == "circular"
        else (LinearDiagramRequest, LinearDiagramOptions)
    )
    return plan_request(request(
        records=(RecordInput(source, record_key="one"),), options=options(orf_prediction=prediction),
    ))


def _features(plan, mode: str, direction: str, separate: bool, resolve: bool):
    layers = create_feature_layers(
        plan.records[0], {}, DRAWN, {"default": "#999999", **{kind: "#999999" for kind in DRAWN}},
        separate, resolve, {}, compute_label_text=False, record_transform=plan.transforms[0],
        placement_slot=FeaturePlacementSlot(mode, direction, separate),
    )
    return list(layers.foreground_features.values())


SLOTS = [
    ("circular", "inside", False, False),
    ("circular", "inside", True, False),
    ("circular", "split", False, True),
    ("circular", "outside", False, False),
    ("linear", "overlay", False, False),
    ("linear", "overlay", True, False),
    ("linear", "overlay", False, True),
    ("linear", "below", False, False),
    ("linear", "above", False, False),
]


@pytest.mark.parametrize("path", [LAMBDA, PLASMID], ids=["lambda", "plasmid"])
@pytest.mark.parametrize(("mode", "direction", "separate", "resolve"), SLOTS)
def test_predicted_orfs_take_own_lanes_beside_the_annotated_features(
    path: Path, mode: str, direction: str, separate: bool, resolve: bool
) -> None:
    before = _features(_plan(GenBankInputSource(path), mode, None), mode, direction, separate, resolve)
    after = _features(_plan(GenBankInputSource(path), mode, OrfPredictionParameters()), mode, direction, separate, resolve)
    annotated = [feature for feature in after if not feature.predicted]
    predicted = [feature for feature in after if feature.predicted]
    assert predicted and len(annotated) == len(before)

    # The annotated features stay where they are drawn without prediction.
    def place(feature):
        return feature.feature_track_id, feature.placement.side, feature.placement.level

    assert {f.source_feature_index: place(f) for f in annotated} == {f.source_feature_index: place(f) for f in before}
    # No predicted ORF shares a lane with an annotated feature, so none covers one.
    assert not {f.feature_track_id for f in predicted} & {f.feature_track_id for f in annotated}
    # Outside the circle, below the line (a one-sided Linear slot keeps its side: the lanes follow the others).
    sides = {f.placement.side for f in predicted}
    assert sides == ({"outward"} if mode == "circular" else {"above"} if direction == "above" else {"below"})
    if direction in ("outside", "above", "below"):
        # A one-sided slot: the predicted lanes come after every lane of the annotated features.
        assert min(f.placement.level for f in predicted) > max(f.placement.level for f in annotated)
    # Predicted ORFs that overlap each other take further lanes.
    length = len(_plan(GenBankInputSource(path), mode, None).records[0].seq)
    lanes: dict[int, list] = {}
    for feature in predicted:
        lanes.setdefault(feature.feature_track_id, []).append(feature)
    for members in lanes.values():
        ends = [dict(zip(("start", "end", "strand"), get_feature_ends(f, length))) for f in members]
        assert not any(
            feature_overlap_bp(a, b, genome_length=length) for i, a in enumerate(ends) for b in ends[i + 1:]
        )


@pytest.mark.parametrize(("mode", "direction"), [("circular", "inside"), ("linear", "overlay")])
def test_a_record_without_features_draws_its_predicted_orfs_as_before(tmp_path: Path, mode: str, direction: str) -> None:
    features = _features(
        _plan(FastaInputSource(_fasta(tmp_path, LAMBDA)), mode, OrfPredictionParameters()), mode, direction, False, False,
    )
    assert len(features) == 110
    # Every ORF on the main lane, as every feature of a slot without --resolve_overlaps.
    assert {feature.feature_track_id for feature in features} == {0}


def test_the_circular_lanes_of_predicted_orfs_lie_outside_the_circle() -> None:
    plan = _plan(GenBankInputSource(LAMBDA), "circular", OrfPredictionParameters())
    features = {f.feature_id: f for f in _features(plan, "circular", "inside", False, False)}
    before = {
        f.feature_id: f for f in _features(_plan(GenBankInputSource(LAMBDA), "circular", None), "circular", "inside", False, False)
    }

    def layout(drawn):
        return build_circular_feature_layout(
            drawn, axis_radius_px=1000.0, width_px=20.0, lane_direction="inside", strandedness=False,
            anchor_radius_px=980.0, lane_spacing_px=4.0,
        )

    with_prediction, without = layout(features), layout(before)
    predicted = {f.feature_track_id for f in features.values() if f.predicted}
    for track_id, lane in with_prediction.lanes_by_track_id.items():
        if track_id in predicted:
            assert lane.inner_px >= 1000.0
        else:
            assert lane == without.lanes_by_track_id[track_id]
    assert with_prediction.all_band_px.outer_px > 1000.0
    assert with_prediction.packing_band_px == without.all_band_px


def test_the_figure_and_the_labels_hold_the_lanes_outside_the_circle() -> None:
    from gbdraw.diagrams.circular.radial_layout import resolve_circular_radial_layout
    from gbdraw.labels.circular import circular_label_radii
    from gbdraw.tracks import CircularTrackSlot
    from tests.test_circular_radial_layout import _small_radial_canvas

    canvas_config, _cfg = _small_radial_canvas()

    def radial(path_plan):
        features = {f.feature_id: f for f in _features(path_plan, "circular", "inside", False, False)}
        return resolve_circular_radial_layout(
            total_length=48502, geometry=canvas_config.geometry(), cfg=canvas_config.profile.config,
            slots=[CircularTrackSlot(id="features", renderer="features"), CircularTrackSlot(id="ticks", renderer="ticks")],
            feature_dict=features,
        )

    layout = radial(_plan(GenBankInputSource(LAMBDA), "circular", OrfPredictionParameters()))
    before = radial(_plan(GenBankInputSource(LAMBDA), "circular", None))
    assert layout.features is not None and before.features is not None
    outside = layout.features.all_band_px.outer_px
    axis = layout.axis.radius_px
    assert outside > axis
    # The slot packs where it packed without prediction; the lanes outside are room no other row takes.
    assert layout.features.packing_band_px == before.features.all_band_px
    assert [slot.packing_band_px for slot in layout.slots] == [slot.packing_band_px for slot in before.slots]
    # The figure holds them and the outer labels keep past them.
    assert layout.outer_content_radius_px >= outside and layout.ring_outer_radius_px >= outside
    assert circular_label_radii(layout).outer_floor_px >= outside
    assert circular_label_radii(before).outer_floor_px == 0.0


def _svg(tmp_path: Path, mode: str, *args: str) -> str:
    result = run_gbdraw([mode, *args, "-o", str(tmp_path / "out"), "-f", "svg", "--overwrite"], cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    return (tmp_path / "out.svg").read_text()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_legend_keeps_the_annotated_cds_row_and_adds_the_predicted_row(tmp_path: Path, mode: str) -> None:
    svg = _svg(tmp_path, mode, "--gbk", str(LAMBDA), "--predict_orfs")
    assert f">{PREDICTED_ORF_CAPTION}<" in svg
    assert ">CDS<" in svg
    assert "other proteins" not in svg
    fasta_only = _svg(tmp_path, mode, "--fasta", str(_fasta(tmp_path, LAMBDA)), "--predict_orfs")
    assert f">{PREDICTED_ORF_CAPTION}<" in fasta_only
    assert ">CDS<" not in fasta_only


def test_a_record_whose_cds_are_all_predicted_has_no_cds_row_beside_an_unused_cds_rule(tmp_path: Path) -> None:
    """OV-3830: a FASTA record draws every CDS as a predicted ORF on its own colour. A CDS rule that colors none of
    them drew a "CDS" row (in the CDS colour, which no feature has); now only "Predicted ORF" names them."""

    table = tmp_path / "colors.tsv"
    table.write_text("CDS\tproduct\tno such product\t#123456\tNever\n", encoding="utf-8")
    svg = _svg(tmp_path, "circular", "--fasta", str(_fasta(tmp_path, LAMBDA)), "--predict_orfs", "-t", str(table))
    assert f'data-legend-key="{PREDICTED_ORF_CAPTION}"' in svg
    assert 'data-legend-key="CDS"' not in svg and 'data-legend-key="Never"' not in svg


def test_an_inside_stack_draws_outside_only_the_lanes_of_their_own() -> None:
    """OV-3829: the lanes an inside stack draws outside the circle line are the placement's own lanes
    (``FeaturePlacementAssignment.own_lane``), not every lane placed outward: with two features slots the features
    carry the first slot's sides, and without prediction no lane moves outside."""

    from gbdraw.diagrams.circular.radial_layout import measure_circular_feature_stack
    from gbdraw.features.placement import FeaturePlacementAssignment

    def measure(**assignment):
        return measure_circular_feature_stack(
            axis_radius_px=400.0, lane_width_px=10.0, lane_spacing_px=2.0, preset="tuckin", strandedness=False,
            track_ids=(0, 1), placements_by_track_id={
                0: FeaturePlacementAssignment("outward", "combined", 0),
                1: FeaturePlacementAssignment("outward", "combined", 1, **assignment),
            },
        )

    placed_outward, own = measure(), measure(own_lane=True)
    assert placed_outward.outside_track_ids == frozenset()
    assert max(placed_outward.lane_centers_by_track_id.values()) < 400.0
    assert own.outside_track_ids == frozenset({1})
    assert own.lane_centers_by_track_id[1] > 400.0 > own.lane_centers_by_track_id[0]


def _predicted_letters(path: Path) -> int:
    from gbdraw.analysis.protein_colinearity import source_protein_letters

    record = copy.copy(SeqIO.read(path, "genbank"))
    record.features = predict_record_orf_features(record, topology="linear", parameters=OrfPredictionParameters())
    return source_protein_letters([record])


@pytest.mark.linear
def test_a_fasta_pair_compares_through_its_predicted_proteins(tmp_path: Path) -> None:
    lam, de3 = _fasta(tmp_path, LAMBDA), _fasta(tmp_path, DE3)
    table = tmp_path / "records.tsv"
    table.write_text(
        "fasta\trecord_id\trow\tcolumn\tregion\n"
        f"{lam}\tNC_001416.1\t1\t1\t\n"
        f"{de3}\tNC_042057.1\t2\t1\t1-20000\n",
        encoding="utf-8",
    )
    losat, log = _fake_losat(tmp_path)
    result = run_gbdraw([
        "linear", "--records_table", str(table), "--predict_orfs", "--losat", "losatp", "--losat_bin", losat,
        "--bitscore", "0", "--evalue", "10", "--identity", "0", "--alignment_length", "0",
        "-o", str(tmp_path / "out"), "-f", "svg", "--overwrite",
    ], cwd=tmp_path)
    assert result.returncode == 0, result.stdout + result.stderr
    # A subject's database is the proteins of every ORF predicted on its whole genome (a search of a whole record
    # against a whole record passes none), never the annotated proteins a FASTA record lacks.
    assert set(_search_spaces(log, "-dbsize")) == {None, _predicted_letters(DE3), _predicted_letters(LAMBDA)}


def _draw(tmp_path: Path, record, *args: str) -> str:
    SeqIO.write(record, tmp_path / "record.gb", "genbank")
    return _svg(tmp_path, "circular", "--gbk", str(tmp_path / "record.gb"), *args)


def test_a_user_rule_captioned_like_the_predicted_row_is_an_ordinary_rule_without_prediction(tmp_path: Path) -> None:
    """Prediction off: a user's CDS rule captioned "Predicted ORF" in #e58fa6 renames the default CDS row as any rule
    does ("other proteins"), as on dev; the predicted row is not a caption and colour pair."""

    table = tmp_path / "colors.tsv"
    table.write_text(f"CDS\tproduct\tterminase\t#e58fa6\t{PREDICTED_ORF_CAPTION}\n", encoding="utf-8")
    svg = _draw(tmp_path, SeqIO.read(LAMBDA, "genbank"), "-t", str(table))
    assert f">{PREDICTED_ORF_CAPTION}<" in svg
    assert ">other proteins<" in svg
    assert ">CDS<" not in svg


def test_an_annotated_cds_with_a_gbdraw_inference_stays_an_annotated_cds(tmp_path: Path) -> None:
    """Prediction on: only the typed mark makes a predicted ORF. An annotated CDS whose /inference reads like
    gbdraw's keeps the default colour, the annotated lanes inside the circle and the "CDS" row."""

    import math
    import re

    from tests.test_circular_leaders_over_lanes import _path_points

    record = SeqIO.read(LAMBDA, "genbank")[0:9000]
    record.annotations = {"molecule_type": "DNA", "topology": "circular"}
    record.features = [feature for feature in record.features if feature.type == "CDS"]
    record.features[0].qualifiers["inference"] = ["ab initio prediction:gbdraw:0.14.5"]
    svg = _draw(tmp_path, record, "--predict_orfs")

    axis = float(re.search(r'<circle cx="0" cy="0" fill="none" r="([0-9.]+)"', svg).group(1))
    paths = re.findall(r'<path d="([^"]+)"[^>]*data-gbdraw-source-feature-index="(\d+)" fill="([^"]+)"', svg)
    marked = [(d, fill) for d, index, fill in paths if index == "0"]
    assert marked
    for d, fill in marked:
        assert fill != PREDICTED_ORF_COLOR
        radii = [math.hypot(float(x), float(y)) for x, y in _path_points(d)]
        assert max(radii) <= axis
    assert any(fill == PREDICTED_ORF_COLOR for _d, _index, fill in paths)
    assert ">CDS<" in svg and f">{PREDICTED_ORF_CAPTION}<" in svg
