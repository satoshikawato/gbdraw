from __future__ import annotations

from dataclasses import replace
import re
from types import SimpleNamespace
from xml.etree import ElementTree

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.diagram as api_diagram_module
import gbdraw.api.request_render as request_render_module
import gbdraw.linear as linear_cli_module
from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    LinearComparison,
    RegionAnnotation,
    parse_record_selector,
)
from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.api.options import (
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearOutputOptions,
    LinearRecordTranslation,
)
from gbdraw.api.record_planning import (
    ResolvedRecordCollection,
    resolve_cli_similarity_alignment_plan,
)
from gbdraw.api.diagram import LinearDiagramMetadata
from gbdraw.api.request_render import build_request_diagram, plan_linear_request
from gbdraw.api.session_compat import _replace_plan_request
from gbdraw.api.requests import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordDisplayOptions,
    RecordPresentation,
)
from gbdraw.diagrams.linear.assemble import _final_record_translations
from gbdraw.exceptions import ValidationError
from gbdraw.features.source import build_source_feature_catalog
from gbdraw.layout.linear_multi_record import LinearRecordPlacement, RecordKey
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentDecisionStatus,
    AlignmentRecordDecision,
    AlignmentResolutionRationale,
    SimilarityAlignmentMode,
    SimilarityAlignmentPlan,
)
from gbdraw.io.regions import parse_region_spec
from gbdraw.io.comparisons import COMPARISON_COLUMNS


def _record(
    record_id: str,
    length: int,
    start: int,
    end: int,
    *,
    strand: int = 1,
    protein_id: str | None = None,
) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, description=record_id)
    record.annotations["topology"] = "linear"
    record.features = [
        SeqFeature(
            SimpleLocation(start, end, strand=strand),
            type="CDS",
            qualifiers={
                "protein_id": [protein_id or f"{record_id}-protein"],
                "translation": ["M" * max(1, (end - start) // 3)],
            },
        )
    ]
    return record


def _anchor(record: SeqRecord, record_key: str) -> AlignmentAnchorIdentity:
    feature = build_source_feature_catalog(record)[0]
    return AlignmentAnchorIdentity(
        record_key=record_key,
        biological_feature_id=feature.biological_feature_id,
        source_feature_index=feature.source_feature_index,
        stable_feature_svg_id=feature.stable_feature_id,
    )


def _plan(
    reference: AlignmentAnchorIdentity,
    target: AlignmentAnchorIdentity,
    *,
    orient_target: bool | None = None,
    mode: SimilarityAlignmentMode = SimilarityAlignmentMode.POSITION,
) -> SimilarityAlignmentPlan:
    return SimilarityAlignmentPlan(
        mode=mode,
        group_id="og-1",
        reference=reference,
        records=(
            AlignmentRecordDecision(
                reference.record_key,
                AlignmentDecisionStatus.REFERENCE,
                AlignmentResolutionRationale.REFERENCE,
                reference,
            ),
            AlignmentRecordDecision(
                target.record_key,
                AlignmentDecisionStatus.ALIGNED,
                AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
                target,
                orient_target,
            ),
        ),
    )


def _request(
    first: SeqRecord,
    second: SeqRecord,
    plan: SimilarityAlignmentPlan,
    *,
    translations: tuple[LinearRecordTranslation, ...] | None = None,
    positions: tuple[str, ...] | None = None,
    overrides: dict[str, object] | None = None,
) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=(
            RecordInput(
                InMemoryRecordSource(first),
                record_key="first",
                presentation=RecordPresentation(reverse_complement=False),
            ),
            RecordInput(
                InMemoryRecordSource(second),
                record_key="second",
                presentation=RecordPresentation(reverse_complement=False),
            ),
        ),
        options=LinearDiagramOptions(
            selected_features_set=("CDS",),
            output=LinearOutputOptions(legend="none"),
            config_overrides={
                "canvas.show_gc": False,
                "canvas.show_skew": False,
                "objects.scale.show": False,
                **(overrides or {}),
            },
        ),
        layout=LinearMultiRecordOptions(
            multi_record_positions=positions,
            record_translations=translations
            or (
                LinearRecordTranslation("first"),
                LinearRecordTranslation("second"),
            ),
        ),
        similarity_alignment=plan,
    )


def _placement(index: int, key: str, *, x: float, width: float) -> LinearRecordPlacement:
    return LinearRecordPlacement(
        record_index=index,
        record_key=RecordKey(key),
        row=index,
        column=0,
        x=x,
        axis_y=10.0 * index,
        sequence_width=width,
        left_inset=0.0,
        right_inset=0.0,
        top_extent=1.0,
        bottom_extent=1.0,
        comparison_top_y=10.0 * index - 1.0,
        comparison_bottom_y=10.0 * index + 1.0,
        px_per_bp=width / 100.0,
    )


def test_absolute_translation_formula_preserves_reference_skipped_and_every_y() -> None:
    reference = AlignmentAnchorIdentity("reference", "ref")
    target = AlignmentAnchorIdentity("target", "target")
    skipped_key = "missing"
    plan = SimilarityAlignmentPlan(
        SimilarityAlignmentMode.POSITION,
        "og-1",
        reference,
        (
            AlignmentRecordDecision(
                "reference",
                AlignmentDecisionStatus.REFERENCE,
                AlignmentResolutionRationale.REFERENCE,
                reference,
            ),
            AlignmentRecordDecision(
                "target",
                AlignmentDecisionStatus.ALIGNED,
                AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
                target,
            ),
            AlignmentRecordDecision(
                skipped_key,
                AlignmentDecisionStatus.SKIPPED,
                AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE,
            ),
        ),
    )
    layout = LinearMultiRecordOptions(
        record_translations=(
            LinearRecordTranslation("reference", 5.0, 3.0),
            LinearRecordTranslation("target", 99.0, -4.0),
            LinearRecordTranslation(skipped_key, -7.0, 8.0),
        )
    )
    placements = {
        0: _placement(0, "reference", x=10.0, width=100.0),
        1: _placement(1, "target", x=20.0, width=200.0),
        2: _placement(2, skipped_key, x=30.0, width=100.0),
    }

    translations = _final_record_translations(
        record_keys=("reference", "target", skipped_key),
        placements=placements,
        layout=layout,
        similarity_alignment=plan,
        anchor_centers=(20.0, 40.0, None),
    )

    # R = 10 + 20 + 5 = 35; T = 20 + (40 * 2) = 100; X = R - T.
    assert translations == ((5.0, 3.0), (-65.0, -4.0), (-7.0, 8.0))

    positive = _final_record_translations(
        record_keys=("reference", "target", skipped_key),
        placements={
            **placements,
            1: _placement(1, "target", x=0.0, width=100.0),
        },
        layout=layout,
        similarity_alignment=plan,
        anchor_centers=(20.0, 5.0, None),
    )
    assert positive[1][0] == pytest.approx(30.0)


def test_planning_derives_effective_orientation_once_then_projects_centers() -> None:
    first = _record("first", 100, 20, 30)
    second = _record("second", 120, 60, 70, strand=-1)
    plan = _plan(
        _anchor(first, "first"),
        _anchor(second, "second"),
        orient_target=True,
        mode=SimilarityAlignmentMode.POSITION_AND_ORIENTATION,
    )
    request = _request(first, second, plan)

    planned = plan_linear_request(request)

    assert [item.presentation.reverse_complement for item in planned.provenance] == [
        False,
        False,
    ]
    assert [transform.source_step for transform in planned.transforms] == [1, -1]
    assert planned.alignment_anchor_centers == pytest.approx((25.0, 55.0))
    assert int(planned.records[1].features[0].location.start) == 50
    assert int(planned.records[1].features[0].location.end) == 60
    assert planned.records[1].features[0].location.strand == 1
    assert int(second.features[0].location.strand) == -1


def test_session_request_replacement_does_not_reapply_effective_orientation() -> None:
    first = _record("first", 100, 20, 30)
    second = _record("second", 120, 60, 70, strand=-1)
    request = _request(
        first,
        second,
        _plan(
            _anchor(first, "first"),
            _anchor(second, "second"),
            orient_target=True,
            mode=SimilarityAlignmentMode.POSITION_AND_ORIENTATION,
        ),
    )
    planned = plan_linear_request(request)

    replaced = _replace_plan_request(
        planned,
        replace(planned.request, options=replace(planned.request.options, evalue=1e-6)),
    )

    assert replaced.records == planned.records
    assert replaced.transforms == planned.transforms
    assert replaced.alignment_anchor_centers == planned.alignment_anchor_centers


@pytest.mark.parametrize(
    ("base_orientation", "override", "expected_step"),
    (
        (False, None, 1),
        (True, None, -1),
        (False, False, 1),
        (False, True, -1),
        (True, False, 1),
        (True, True, -1),
    ),
)
def test_effective_orientation_is_absolute_against_record_presentation(
    base_orientation: bool,
    override: bool | None,
    expected_step: int,
) -> None:
    first = _record("first", 100, 20, 30)
    second = _record("second", 120, 60, 70)
    mode = (
        SimilarityAlignmentMode.POSITION
        if override is None
        else SimilarityAlignmentMode.POSITION_AND_ORIENTATION
    )
    plan = _plan(
        _anchor(first, "first"),
        _anchor(second, "second"),
        orient_target=override,
        mode=mode,
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(first), record_key="first"),
            RecordInput(
                InMemoryRecordSource(second),
                record_key="second",
                presentation=RecordPresentation(
                    reverse_complement=base_orientation
                ),
            ),
        ),
        layout=LinearMultiRecordOptions(
            record_translations=(
                LinearRecordTranslation("first"),
                LinearRecordTranslation("second"),
            )
        ),
        similarity_alignment=plan,
    )

    planned = plan_linear_request(request)

    assert planned.provenance[1].presentation.reverse_complement is base_orientation
    assert planned.transforms[1].source_step == expected_step


def test_crop_and_circular_display_offset_are_applied_before_anchor_translation() -> None:
    first = _record("first", 100, 20, 30)
    cropped = _record("second", 120, 60, 70)
    plan = _plan(_anchor(first, "first"), _anchor(cropped, "second"))
    crop_request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(first), record_key="first"),
            RecordInput(
                InMemoryRecordSource(cropped),
                record_key="second",
                region=parse_region_spec("second:41-100"),
            ),
        ),
        layout=LinearMultiRecordOptions(
            record_translations=(
                LinearRecordTranslation("first"),
                LinearRecordTranslation("second"),
            )
        ),
        similarity_alignment=plan,
    )
    cropped_plan = plan_linear_request(crop_request)
    assert cropped_plan.alignment_anchor_centers[1] == pytest.approx(
        cropped_plan.transforms[1].source_position_to_display_offset(65.0)
    )

    circular = _record("second", 120, 60, 70)
    circular.annotations["topology"] = "circular"
    circular_plan = _plan(_anchor(first, "first"), _anchor(circular, "second"))
    display_request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(first), record_key="first"),
            RecordInput(
                InMemoryRecordSource(circular),
                record_key="second",
                display=RecordDisplayOptions(
                    is_circular=True,
                    start_coordinate=51,
                ),
            ),
        ),
        layout=LinearMultiRecordOptions(
            record_translations=(
                LinearRecordTranslation("first"),
                LinearRecordTranslation("second"),
            )
        ),
        similarity_alignment=circular_plan,
    )
    displayed = plan_linear_request(display_request)
    assert displayed.alignment_anchor_centers[1] == pytest.approx(15.0)


@pytest.mark.parametrize(
    ("overrides", "positions"),
    (
        ({"canvas.linear.normalize_length": True}, None),
        ({"canvas.linear.align_center": True}, None),
        (
            {
                "canvas.linear.normalize_length": True,
                "canvas.linear.align_center": True,
            },
            None,
        ),
        ({}, ("#1@1", "#2@1")),
    ),
    ids=(
        "normalize-length",
        "center-alignment",
        "normalize-and-center",
        "same-row",
    ),
)
def test_final_placements_align_anchors_and_bounds_do_not_clip(
    overrides: dict[str, object],
    positions: tuple[str, ...] | None,
) -> None:
    first = _record("first", 100, 10, 20)
    second = _record("second", 200, 150, 170)
    plan = _plan(_anchor(first, "first"), _anchor(second, "second"))
    request = _request(
        first,
        second,
        plan,
        translations=(
            LinearRecordTranslation("first", -30.0, 4.0),
            LinearRecordTranslation("second", 80.0, -3.0),
        ),
        positions=positions,
        overrides=overrides,
    )

    prepared = build_request_diagram(request)
    geometry = prepared.drawing._gbdraw_track_slot_geometry["records"]
    centers = plan_linear_request(request).alignment_anchor_centers
    world_centers = [
        item["axisXpx"]
        + centers[index] * item["sequenceWidthPx"] / len(prepared.records[index])
        for index, item in enumerate(geometry)
    ]
    composition = prepared.drawing._gbdraw_linear_composition_plan
    primary = composition.placement_for("primary").final_bounds

    assert world_centers[0] == pytest.approx(world_centers[1])
    assert primary.min_x >= composition.canvas_bounds.min_x
    assert primary.max_x <= composition.canvas_bounds.max_x
    assert primary.min_y >= composition.canvas_bounds.min_y
    assert primary.max_y <= composition.canvas_bounds.max_y
    if positions is not None:
        assert geometry[1]["axisYpx"] - geometry[0]["axisYpx"] == pytest.approx(-7.0)


def test_final_translation_is_shared_by_every_record_geometry_consumer() -> None:
    first = _record("first", 100, 10, 20)
    second = _record("second", 200, 150, 170)
    first.features[0].qualifiers["gene"] = ["first-gene"]
    second.features[0].qualifiers["gene"] = ["second-gene"]
    request = _request(
        first,
        second,
        _plan(_anchor(first, "first"), _anchor(second, "second")),
        translations=(
            LinearRecordTranslation("first", -30.0, 4.0),
            LinearRecordTranslation("second", 80.0, -3.0),
        ),
        overrides={
            "objects.scale.show": True,
            "canvas.linear.ruler_on_axis": True,
            "labels.linear.scope": "all",
        },
    )
    comparison = LinearComparison(
        0,
        1,
        pd.DataFrame(
            [["q", "s", 90.0, 10, 0, 0, 11, 20, 151, 170, 1e-5, 100]],
            columns=COMPARISON_COLUMNS,
        ),
    )
    annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "target-region",
                        CoordinateSpan(
                            parse_record_selector("#2"),
                            140,
                            180,
                        ),
                        label="target",
                        mark="band",
                    ),
                ),
            ),
        )
    )
    request = replace(
        request,
        options=replace(
            request.options,
            annotations=annotations,
            linear_comparisons=(comparison,),
        ),
    )

    base = build_request_diagram(replace(request, similarity_alignment=None)).drawing
    aligned = build_request_diagram(request).drawing

    def parsed(drawing):
        return ElementTree.fromstring(drawing.tostring())

    def find(root, predicate):
        return next(element for element in root.iter() if predicate(element))

    def local_translation(element):
        translations = re.findall(
            r"translate\(\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)\s*\)",
            element.attrib["transform"],
        )
        return tuple(float(value) for value in translations[-1])

    base_root = parsed(base)
    aligned_root = parsed(aligned)
    targets = []
    for root in (base_root, aligned_root):
        record = find(
            root,
            lambda element: element.attrib.get("data-gbdraw-record-id") == "second"
            and element.attrib.get("id", "").startswith("record_group_"),
        )
        definition = find(
            root,
            lambda element: element.attrib.get("data-gbdraw-record-id") == "second"
            and element.attrib.get("data-gbdraw-role") == "record-definition",
        )
        annotation = find(
            root,
            lambda element: element.attrib.get("id", "").startswith(
                "gbdraw-annotation-track-annotations_1-2"
            ),
        )
        comparison_group = find(
            root,
            lambda element: element.attrib.get("id") == "comparison1",
        )
        comparison_path = next(
            child for child in comparison_group if child.tag.endswith("path")
        )
        path_x = [
            float(value)
            for value in re.findall(r"[-+0-9.eE]+", comparison_path.attrib["d"])
        ][::2]
        targets.append(
            (
                local_translation(record),
                local_translation(definition),
                local_translation(annotation),
                path_x[-2:],
            )
        )
        assert any(child.tag.endswith("line") for child in record)
        assert any(child.text == "second-gene" for child in record)

    base_target, aligned_target = targets
    target_delta_x = aligned_target[0][0] - base_target[0][0]
    assert target_delta_x < 0.0
    assert aligned_target[0][1] == pytest.approx(base_target[0][1])
    assert aligned_target[1][0] - base_target[1][0] == pytest.approx(target_delta_x)
    assert aligned_target[1][1] == pytest.approx(base_target[1][1])
    assert aligned_target[2][0] - base_target[2][0] == pytest.approx(target_delta_x)
    assert aligned_target[2][1] == pytest.approx(base_target[2][1])
    assert [
        after - before
        for before, after in zip(base_target[3], aligned_target[3], strict=True)
    ] == pytest.approx([target_delta_x, target_delta_x])

    for drawing, target in zip((base, aligned), targets, strict=True):
        record_geometry = drawing._gbdraw_track_slot_geometry["records"][1]
        primary = drawing._gbdraw_linear_composition_plan.placement_for("primary")
        assert primary is not None
        assert record_geometry["axisYpx"] - primary.dy == pytest.approx(
            target[0][1]
        )
    base_axis = base._gbdraw_track_slot_geometry["records"][1]["axisXpx"]
    aligned_axis = aligned._gbdraw_track_slot_geometry["records"][1]["axisXpx"]
    base_primary = base._gbdraw_linear_composition_plan.placement_for("primary")
    aligned_primary = aligned._gbdraw_linear_composition_plan.placement_for("primary")
    assert base_primary is not None and aligned_primary is not None
    assert (aligned_axis - aligned_primary.dx) - (
        base_axis - base_primary.dx
    ) == pytest.approx(target_delta_x)


def test_repeated_planning_and_rendering_are_geometry_idempotent() -> None:
    first = _record("first", 100, 15, 25)
    second = _record("second", 160, 100, 120)
    request = _request(
        first,
        second,
        _plan(_anchor(first, "first"), _anchor(second, "second")),
    )

    first_plan = plan_linear_request(request)
    second_plan = plan_linear_request(request)
    first_render = build_request_diagram(request)
    second_render = build_request_diagram(request)

    assert first_plan.alignment_anchor_centers == second_plan.alignment_anchor_centers
    assert first_render.drawing.tostring() == second_render.drawing.tostring()


def test_typed_request_rejects_legacy_string_before_rendering() -> None:
    from gbdraw.api import SimilarityAlignmentPlan as PublicSimilarityAlignmentPlan

    assert PublicSimilarityAlignmentPlan is SimilarityAlignmentPlan
    with pytest.raises(ValidationError, match="similarity_alignment"):
        LinearDiagramRequest(
            records=(
                RecordInput(
                    InMemoryRecordSource(_record("record", 100, 10, 20)),
                    record_key="record",
                ),
            ),
            similarity_alignment="og-legacy",  # type: ignore[arg-type]
        )


def _orthogroups_for_collection(
    collection: ResolvedRecordCollection,
    *,
    ambiguous_second: bool = False,
) -> OrthogroupResult:
    members: list[OrthogroupMember] = []
    for index, (record, provenance) in enumerate(
        zip(collection.records, collection.provenance, strict=True)
    ):
        catalog = provenance.source_feature_catalog
        assert catalog is not None
        aliases = ("ref-protein",) if index == 0 else ("target-a", "target-b")
        selected = catalog[:2] if ambiguous_second and index == 1 else catalog[:1]
        for feature, protein_id in zip(
            selected, aliases[: len(selected)], strict=True
        ):
            members.append(
                OrthogroupMember(
                    orthogroup_id="og-1",
                    protein_id=protein_id,
                    record_index=index,
                    feature_index=feature.source_feature_index,
                    record_id=record.id,
                    label=protein_id,
                    start=min(part[0] for part in feature.location_parts),
                    end=max(part[1] for part in feature.location_parts),
                    strand=1,
                    feature_svg_id=feature.stable_feature_id,
                    source_protein_id=protein_id,
                )
            )
    return OrthogroupResult(
        orthogroups={"og-1": members},
        member_by_protein_id={member.protein_id: member for member in members},
    )


def test_cli_adapter_is_exact_noninteractive_and_reports_record_candidates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = _record("first", 100, 10, 20, protein_id="ref-protein")
    second = _record("second", 100, 30, 40, protein_id="target-a")
    second.features.append(
        SeqFeature(
            SimpleLocation(60, 70, strand=1),
            type="CDS",
            qualifiers={"protein_id": ["target-b"], "translation": ["MMM"]},
        )
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(first), record_key="first"),
            RecordInput(InMemoryRecordSource(second), record_key="second"),
        )
    )
    planned = plan_linear_request(request)
    collection = ResolvedRecordCollection(planned.records, planned.provenance)
    orthogroups = _orthogroups_for_collection(collection, ambiguous_second=True)
    monkeypatch.setattr(
        "builtins.input",
        lambda *_args, **_kwargs: pytest.fail("CLI alignment must not prompt"),
    )

    with pytest.raises(ValidationError, match="not Similarity Group ID"):
        resolve_cli_similarity_alignment_plan(
            collection,
            orthogroups,
            exact_reference="og-1",
        )
    with pytest.raises(ValidationError, match="record 'second'.*target-[ab].*target-[ab]"):
        resolve_cli_similarity_alignment_plan(
            collection,
            orthogroups,
            exact_reference="ref-protein",
        )

    unique = _orthogroups_for_collection(collection)
    resolved = resolve_cli_similarity_alignment_plan(
        collection,
        unique,
        exact_reference="ref-protein",
    )
    assert resolved.reference.record_key == "first"
    assert [decision.status for decision in resolved.records] == [
        AlignmentDecisionStatus.REFERENCE,
        AlignmentDecisionStatus.ALIGNED,
    ]

    reference_only = OrthogroupResult(
        orthogroups={"og-1": unique.orthogroups["og-1"][:1]},
        member_by_protein_id={"ref-protein": unique.orthogroups["og-1"][0]},
    )
    missing = resolve_cli_similarity_alignment_plan(
        collection,
        reference_only,
        exact_reference="ref-protein",
    )
    assert missing.records[1].status is AlignmentDecisionStatus.SKIPPED
    assert missing.records[1].rationale is (
        AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE
    )


def test_supplied_plan_does_not_invoke_protein_analysis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = _record("first", 100, 10, 20)
    second = _record("second", 100, 30, 40)
    request = _request(
        first,
        second,
        _plan(_anchor(first, "first"), _anchor(second, "second")),
    )
    monkeypatch.setattr(
        api_diagram_module,
        "_invoke_protein_analysis_helper",
        lambda *_args, **_kwargs: pytest.fail(
            "A supplied alignment plan must not start LOSATP"
        ),
    )

    build_request_diagram(request)


def test_cli_materializes_typed_plan_and_reuses_completed_analysis(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = (
        _record("first", 100, 10, 20, protein_id="ref-protein"),
        _record("second", 100, 30, 40, protein_id="target-a"),
    )
    captured: dict[str, object] = {"analysis_builds": 0}
    monkeypatch.setattr(
        request_render_module,
        "load_gbks",
        lambda *args, **_kwargs: [
            records[1] if "b.gb" in str(args[0]) else records[0]
        ],
    )
    monkeypatch.setattr(
        request_render_module,
        "read_color_table",
        lambda _path: None,
    )
    monkeypatch.setattr(
        request_render_module,
        "read_feature_visibility_file",
        lambda _path: None,
    )
    monkeypatch.setattr(
        "builtins.input",
        lambda *_args, **_kwargs: pytest.fail("CLI alignment must not prompt"),
    )

    def fake_analysis_build(plan, *, artifacts):
        captured["analysis_builds"] = int(captured["analysis_builds"]) + 1
        collection = ResolvedRecordCollection(plan.records, plan.provenance)
        return SimpleNamespace(
            linear_metadata=LinearDiagramMetadata(
                protein_comparisons=(),
                linear_comparisons=(),
                orthogroups=_orthogroups_for_collection(collection),
            ),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
        )

    def fake_render(request, *, artifacts, **_kwargs):
        captured["request"] = request
        captured["artifacts"] = artifacts
        planned = plan_linear_request(request)
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=planned.records,
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=planned.request,
        )

    monkeypatch.setattr(linear_cli_module, "build_request_plan_diagram", fake_analysis_build)
    monkeypatch.setattr(linear_cli_module, "render_request", fake_render)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "orthogroup",
            "--align_orthogroup_feature",
            "ref-protein",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "aligned"),
        ]
    )

    request = captured["request"]
    assert isinstance(request, LinearDiagramRequest)
    assert request.similarity_alignment is not None
    assert request.options.protein_blastp_mode == "none"
    assert not hasattr(request.options, "align_orthogroup_feature")
    assert captured["analysis_builds"] == 1
