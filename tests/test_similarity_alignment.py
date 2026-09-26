from __future__ import annotations

import ast
import inspect
from itertools import permutations

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.layout import similarity_alignment as alignment_module
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentChoiceKind,
    AlignmentDecisionStatus,
    AlignmentEvidenceEdge,
    AlignmentRecordChoice,
    AlignmentRecordDecision,
    AlignmentResolutionRationale,
    AmbiguousAlignmentRecord,
    SimilarityAlignmentCandidate,
    AlignmentStrandRelation,
    AlignmentRecommendationReason,
    SimilarityAlignmentPlan,
    resolve_similarity_alignment,
)


GROUP = "group-1"


def _anchor(
    record_key: str,
    feature_id: str,
    *,
    source_index: int | None = None,
    stable_id: str | None = None,
) -> AlignmentAnchorIdentity:
    return AlignmentAnchorIdentity(
        record_key=record_key,
        biological_feature_id=feature_id,
        source_feature_index=source_index,
        stable_feature_svg_id=stable_id,
    )


def _candidate(
    record_key: str,
    feature_id: str,
    *,
    source_index: int | None = None,
    stable_id: str | None = None,
    group_id: str = GROUP,
    strand: int | None = 1,
    center: float | None = 10.0,
    unique: bool = True,
    hidden: bool = False,
    representative: bool = False,
    role: str = "anchor",
) -> SimilarityAlignmentCandidate:
    return SimilarityAlignmentCandidate(
        group_id=group_id,
        anchor=_anchor(
            record_key,
            feature_id,
            source_index=source_index,
            stable_id=stable_id,
        ),
        displayed_strand=strand,
        center_mappable=center is not None,
        display_center=center,
        identity_is_unique=unique,
        hidden=hidden,
        representative=representative,
        role=role,
    )


def _decision(resolution, record_key: str) -> AlignmentRecordDecision:
    outcome = next(
        record for record in resolution.records if record.record_key == record_key
    )
    assert isinstance(outcome, AlignmentRecordDecision)
    return outcome


def test_exact_non_representative_reference_and_zero_or_one_candidate() -> None:
    reference = _candidate(
        "record-b",
        "clicked-inparalog",
        source_index=4,
        stable_id="stable-shared",
        representative=False,
    )
    same_record_representative = _candidate(
        "record-b",
        "representative",
        source_index=8,
        representative=True,
    )
    only_target = _candidate(
        "record-a",
        "only-target",
        source_index=2,
        representative=False,
    )

    resolution = resolve_similarity_alignment(
        record_keys=("record-a", "record-b", "record-c"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(same_record_representative, only_target, reference),
    )

    assert [record.record_key for record in resolution.records] == [
        "record-a",
        "record-b",
        "record-c",
    ]
    assert resolution.reference == reference.anchor
    assert _decision(resolution, "record-a") == AlignmentRecordDecision(
        record_key="record-a",
        status="aligned",
        rationale="only_usable_candidate",
        anchor=only_target.anchor,
    )
    assert _decision(resolution, "record-b").anchor == reference.anchor
    assert _decision(resolution, "record-c").rationale is (
        AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE
    )
    plan = resolution.require_plan()
    plan.validate_record_coverage(("record-c", "record-b", "record-a"))


def test_several_candidates_remain_ambiguous_without_ranking_fallbacks() -> None:
    reference = _candidate("reference", "clicked", center=30)
    representative = _candidate(
        "target",
        "feature-z",
        center=900,
        representative=True,
        role="coortholog",
    )
    non_representative = _candidate(
        "target",
        "feature-a",
        center=1,
        representative=False,
        role="inparalog",
    )

    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(representative, reference, non_representative),
    )

    assert resolution.plan is None
    assert len(resolution.ambiguities) == 1
    ambiguity = resolution.ambiguities[0]
    assert isinstance(ambiguity, AmbiguousAlignmentRecord)
    assert [item.anchor for item in ambiguity.candidates] == [
        non_representative.anchor,
        representative.anchor,
    ]
    assert ambiguity.recommended_anchor == representative.anchor
    assert ambiguity.recommendation_reason is (
        AlignmentRecommendationReason.UNIQUE_REPRESENTATIVE
    )
    with pytest.raises(ValidationError, match="Select or Skip.*target"):
        resolution.require_plan()


def test_explicit_select_and_skip_have_first_priority() -> None:
    reference = _candidate("reference", "clicked")
    candidate_a = _candidate("selected", "feature-a")
    candidate_b = _candidate("selected", "feature-b")
    skipped = _candidate("skipped", "only")
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "selected", "skipped"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, candidate_a, candidate_b, skipped),
        edges=(
            AlignmentEvidenceEdge(
                GROUP, reference.anchor, candidate_a.anchor, "rbh"
            ),
        ),
        choices=(
            AlignmentRecordChoice(
                "selected", AlignmentChoiceKind.SELECT, candidate_b.anchor
            ),
            AlignmentRecordChoice("skipped", AlignmentChoiceKind.SKIP),
        ),
    )

    selected = _decision(resolution, "selected")
    assert selected.anchor == candidate_b.anchor
    assert selected.rationale is AlignmentResolutionRationale.USER_SELECTED
    skipped_decision = _decision(resolution, "skipped")
    assert skipped_decision.status is AlignmentDecisionStatus.SKIPPED
    assert skipped_decision.rationale is (
        AlignmentResolutionRationale.SKIPPED_BY_USER
    )


@pytest.mark.parametrize("reverse_edge", [False, True])
def test_one_direct_rbh_resolves_independent_of_edge_direction(
    reverse_edge: bool,
) -> None:
    reference = _candidate("reference", "clicked")
    connected = _candidate("target", "connected")
    other = _candidate("target", "other")
    query, subject = (
        (connected.anchor, reference.anchor)
        if reverse_edge
        else (reference.anchor, connected.anchor)
    )
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(other, reference, connected),
        edges=(AlignmentEvidenceEdge(GROUP, query, subject, "rbh"),),
    )

    decision = _decision(resolution, "target")
    assert decision.anchor == connected.anchor
    assert decision.rationale is AlignmentResolutionRationale.UNIQUE_DIRECT_RBH


def test_duplicate_edge_evidence_does_not_create_multiple_rbh_candidates() -> None:
    reference = _candidate("reference", "clicked")
    connected = _candidate("target", "connected")
    other = _candidate("target", "other")
    edge = AlignmentEvidenceEdge(GROUP, reference.anchor, connected.anchor, "rbh")

    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, connected, other),
        edges=(edge, edge),
    )

    assert _decision(resolution, "target").anchor == connected.anchor


def test_multiple_direct_rbhs_remain_ambiguous() -> None:
    reference = _candidate("reference", "clicked")
    first = _candidate("target", "first")
    second = _candidate("target", "second")
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, first, second),
        edges=(
            AlignmentEvidenceEdge(GROUP, reference.anchor, first.anchor, "rbh"),
            AlignmentEvidenceEdge(GROUP, second.anchor, reference.anchor, "rbh"),
        ),
    )

    ambiguity = resolution.ambiguities[0]
    assert [anchor.biological_feature_id for anchor in ambiguity.direct_rbh_candidates] == [
        "first",
        "second",
    ]


def test_non_rbh_and_multihop_evidence_do_not_resolve_ambiguity() -> None:
    reference = _candidate("reference", "clicked")
    intermediate = _candidate("middle", "middle")
    first = _candidate("target", "first")
    second = _candidate("target", "second")
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "middle", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, intermediate, first, second),
        edges=(
            AlignmentEvidenceEdge(
                GROUP, reference.anchor, first.anchor, "coortholog"
            ),
            AlignmentEvidenceEdge(
                GROUP, reference.anchor, intermediate.anchor, "rbh"
            ),
            AlignmentEvidenceEdge(GROUP, intermediate.anchor, first.anchor, "rbh"),
        ),
    )

    assert resolution.ambiguities[0].record_key == "target"
    assert resolution.ambiguities[0].direct_rbh_candidates == ()


def test_hidden_mappable_candidate_is_usable_and_crop_excluded_center_is_not() -> None:
    reference = _candidate("reference", "clicked")
    hidden = _candidate("hidden", "hidden-feature", hidden=True)
    partially_overlapping = _candidate("cropped", "partial", center=None)
    non_unique = _candidate("non-unique", "duplicate", unique=False)
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "hidden", "cropped", "non-unique"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, hidden, partially_overlapping, non_unique),
    )

    assert _decision(resolution, "hidden").anchor == hidden.anchor
    assert _decision(resolution, "cropped").rationale is (
        AlignmentResolutionRationale.SKIPPED_UNMAPPABLE
    )
    assert _decision(resolution, "non-unique").rationale is (
        AlignmentResolutionRationale.SKIPPED_UNMAPPABLE
    )


def test_other_group_member_is_not_a_candidate() -> None:
    reference = _candidate("reference", "clicked")
    other_group = _candidate("target", "other", group_id="group-2")
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, other_group),
    )
    assert _decision(resolution, "target").rationale is (
        AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE
    )


@pytest.mark.parametrize("reference_strand", (None, -1, 1))
@pytest.mark.parametrize("target_strand", (None, -1, 1))
def test_review_candidate_strand_relation(
    reference_strand: int | None,
    target_strand: int | None,
) -> None:
    reference = _candidate("reference", "clicked", strand=reference_strand)
    target = _candidate("target", "target", strand=target_strand)
    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, target),
    )
    expected = (
        AlignmentStrandRelation.UNKNOWN
        if reference_strand is None or target_strand is None
        else AlignmentStrandRelation.SAME
        if reference_strand == target_strand
        else AlignmentStrandRelation.OPPOSITE
    )
    assert resolution.review_rows[1].candidates[0].strand_relation is expected
    assert _decision(resolution, "target").anchor == target.anchor


def test_choice_and_decision_reject_removed_orientation_fields() -> None:
    anchor = _anchor("target", "feature")
    with pytest.raises(TypeError, match="orientation_policy"):
        AlignmentRecordChoice("target", "select", anchor, orientation_policy="preserve")
    with pytest.raises(TypeError, match="effective_reverse_complement"):
        AlignmentRecordDecision(
            "target", "aligned", "user_selected", anchor,
            effective_reverse_complement=False,
        )


def test_input_permutations_produce_identical_resolution() -> None:
    reference = _candidate("reference", "clicked")
    rbh = _candidate("target", "rbh")
    other = _candidate("target", "other")
    edges = (
        AlignmentEvidenceEdge(GROUP, reference.anchor, rbh.anchor, "rbh"),
        AlignmentEvidenceEdge(GROUP, reference.anchor, other.anchor, "coortholog"),
    )
    candidates = (reference, rbh, other)
    expected = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=candidates,
        edges=edges,
    )

    for candidate_order in permutations(candidates):
        for edge_order in permutations(edges):
            assert resolve_similarity_alignment(
                record_keys=("reference", "target"),
                group_id=GROUP,
                reference=reference.anchor,
                candidates=candidate_order,
                edges=edge_order,
            ) == expected


@pytest.mark.parametrize(
    ("record_keys", "match"),
    [
        (("record", "record"), "record keys must be unique"),
        ((), "at least one record"),
    ],
)
def test_invalid_record_key_coverage_is_rejected(record_keys, match: str) -> None:
    reference = _candidate("record", "clicked")
    with pytest.raises(ValidationError, match=match):
        resolve_similarity_alignment(
            record_keys=record_keys,
            group_id=GROUP,
            reference=reference.anchor,
            candidates=(reference,),
        )


def test_conflicting_identity_aliases_are_rejected() -> None:
    reference = _candidate("reference", "clicked")
    first = _candidate("target", "first", source_index=4)
    second = _candidate("target", "second", source_index=4)
    with pytest.raises(ValidationError, match="conflicting source-feature aliases"):
        resolve_similarity_alignment(
            record_keys=("reference", "target"),
            group_id=GROUP,
            reference=reference.anchor,
            candidates=(reference, first, second),
        )

    with pytest.raises(ValidationError, match="Reference identity aliases conflict"):
        resolve_similarity_alignment(
            record_keys=("reference",),
            group_id=GROUP,
            reference=_anchor("reference", "clicked", source_index=9),
            candidates=(
                _candidate("reference", "clicked", source_index=3),
            ),
        )


def test_duplicate_candidate_identity_is_rejected() -> None:
    reference = _candidate("reference", "clicked")
    target = _candidate("target", "target")
    with pytest.raises(ValidationError, match="Duplicate alignment candidate"):
        resolve_similarity_alignment(
            record_keys=("reference", "target"),
            group_id=GROUP,
            reference=reference.anchor,
            candidates=(reference, target, target),
        )


def test_stale_explicit_selection_is_rejected() -> None:
    reference = _candidate("reference", "clicked")
    current = _candidate("target", "current")
    stale = _anchor("target", "removed")
    with pytest.raises(ValidationError, match="stale or unusable"):
        resolve_similarity_alignment(
            record_keys=("reference", "target"),
            group_id=GROUP,
            reference=reference.anchor,
            candidates=(reference, current),
            choices=(AlignmentRecordChoice("target", "select", stale),),
        )


def test_missing_or_unmappable_reference_is_rejected() -> None:
    reference = _candidate("reference", "clicked")
    with pytest.raises(ValidationError, match="not a member"):
        resolve_similarity_alignment(
            record_keys=("reference",),
            group_id=GROUP,
            reference=reference.anchor,
            candidates=(),
        )
    with pytest.raises(ValidationError, match="does not resolve uniquely"):
        resolve_similarity_alignment(
            record_keys=("reference",),
            group_id=GROUP,
            reference=reference.anchor,
            candidates=(_candidate("reference", "clicked", center=None),),
        )


def test_duplicate_decisions_and_invalid_plan_are_rejected() -> None:
    anchor = _anchor("record", "feature")
    reference_decision = AlignmentRecordDecision(
        "record", "reference", "reference", anchor
    )
    with pytest.raises(ValidationError, match="duplicate record decisions"):
        SimilarityAlignmentPlan(
            group_id=GROUP,
            reference=anchor,
            records=(reference_decision, reference_decision),
        )
    with pytest.raises(ValidationError, match="exactly one matching"):
        SimilarityAlignmentPlan(
            group_id=GROUP,
            reference=anchor,
            records=(
                AlignmentRecordDecision(
                    "record", "skipped", "skipped_no_candidate"
                ),
            ),
        )


def test_plan_rejects_unknown_or_missing_record_coverage() -> None:
    anchor = _anchor("reference", "feature")
    plan = SimilarityAlignmentPlan(
        group_id=GROUP,
        reference=anchor,
        records=(
            AlignmentRecordDecision(
                "reference", "reference", "reference", anchor
            ),
        ),
    )
    with pytest.raises(ValidationError, match="missing=.*target"):
        plan.validate_record_coverage(("reference", "target"))


def test_ambiguous_outcome_rejects_conflicting_direct_evidence_alias() -> None:
    first = _candidate("target", "first", source_index=1)
    second = _candidate("target", "second", source_index=2)
    with pytest.raises(ValidationError, match="aliases conflict"):
        AmbiguousAlignmentRecord(
            record_key="target",
            candidates=(first, second),
            direct_rbh_candidates=(
                _anchor("target", "first", source_index=9),
            ),
        )


@pytest.mark.parametrize("center", [float("nan"), float("inf"), -float("inf")])
def test_candidate_center_must_be_finite(center: float) -> None:
    with pytest.raises(ValidationError, match="finite"):
        _candidate("record", "feature", center=center)


@pytest.mark.parametrize("strand", [True, False, 1.0, "1"])
def test_candidate_strand_requires_exact_supported_integer(strand: object) -> None:
    with pytest.raises(ValidationError, match="displayed_strand"):
        _candidate("record", "feature", strand=strand)  # type: ignore[arg-type]


def test_plan_schema_requires_exact_integer_two() -> None:
    anchor = _anchor("record", "feature")
    decision = AlignmentRecordDecision("record", "reference", "reference", anchor)
    with pytest.raises(ValidationError, match="schema must be 2"):
        SimilarityAlignmentPlan(
            group_id=GROUP,
            reference=anchor,
            records=(decision,),
            schema=1.0,  # type: ignore[arg-type]
        )


def test_candidate_one_recommendation_is_canonical_and_ignores_non_direct_facts() -> None:
    reference = _candidate("reference", "clicked")
    first = _candidate("target", "alpha", center=900, hidden=True, role="inparalog")
    second = _candidate("target", "zeta", center=1, role="coortholog")
    middle = _candidate("middle", "bridge")
    supportive_edges = (
        AlignmentEvidenceEdge(GROUP, reference.anchor, middle.anchor, "rbh"),
        AlignmentEvidenceEdge(GROUP, middle.anchor, second.anchor, "rbh"),
        AlignmentEvidenceEdge(GROUP, reference.anchor, second.anchor, "coortholog"),
    )
    for order in (("reference", "middle", "target"), ("target", "reference", "middle")):
        for candidates in permutations((reference, first, second, middle)):
            result = resolve_similarity_alignment(
                record_keys=order,
                group_id=GROUP,
                reference=reference.anchor,
                candidates=candidates,
                edges=supportive_edges * 3,
            )
            ambiguity = result.ambiguities[0]
            assert ambiguity.record_key == "target"
            assert ambiguity.recommended_anchor == first.anchor
            assert ambiguity.recommendation_reason is (
                AlignmentRecommendationReason.DETERMINISTIC_CANDIDATE_1
            )
            assert result.plan is None

    candidate_fields = set(SimilarityAlignmentCandidate.__dataclass_fields__)
    resolver_inputs = set(inspect.signature(resolve_similarity_alignment).parameters)
    assert not {"viewport", "scroll", "ribbon_geometry", "confidence_score", "supporting_edges"} & (
        candidate_fields | resolver_inputs
    )


def test_recommendation_can_be_replaced_or_skipped_only_by_explicit_choice() -> None:
    reference = _candidate("reference", "clicked")
    first = _candidate("target", "alpha", representative=True)
    second = _candidate("target", "zeta")
    facts = dict(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, first, second),
    )
    unresolved = resolve_similarity_alignment(**facts)
    assert unresolved.ambiguities[0].recommended_anchor == first.anchor
    assert unresolved.plan is None

    replaced = resolve_similarity_alignment(
        **facts,
        choices=(AlignmentRecordChoice("target", "select", second.anchor),),
    ).require_plan()
    assert replaced.records[1].anchor == second.anchor
    assert replaced.records[1].rationale is AlignmentResolutionRationale.USER_SELECTED

    skipped = resolve_similarity_alignment(
        **facts,
        choices=(AlignmentRecordChoice("target", "skip"),),
    ).require_plan()
    assert skipped.records[1].status is AlignmentDecisionStatus.SKIPPED


def test_plan_decision_consistency_and_schema() -> None:
    reference = _anchor("reference", "reference")
    target = _anchor("target", "target")
    assert AlignmentRecordDecision("target", "aligned", "user_selected", target).anchor == target
    with pytest.raises(ValidationError, match="inconsistent"):
        AlignmentRecordDecision("target", "skipped", "skipped_by_user", target)
    with pytest.raises(ValidationError, match="schema must be 2"):
        SimilarityAlignmentPlan(
            group_id=GROUP,
            reference=reference,
            records=(AlignmentRecordDecision("reference", "reference", "reference", reference),),
            schema=1,
        )


def test_domain_module_has_no_analysis_ui_or_render_runtime_dependency() -> None:
    tree = ast.parse(inspect.getsource(alignment_module))
    imports = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imports.update(
        node.module or ""
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    )
    forbidden = (
        "gbdraw.analysis",
        "gbdraw.diagrams",
        "gbdraw.render",
        "gbdraw.web",
        "Bio",
        "pandas",
    )
    assert not {
        module
        for module in imports
        if any(module == prefix or module.startswith(f"{prefix}.") for prefix in forbidden)
    }
