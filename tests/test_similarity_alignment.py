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
    SimilarityAlignmentMode,
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
    reversed_: bool = False,
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
        effective_reverse_complement=reversed_,
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


def test_opposite_known_strands_record_absolute_orientation_override() -> None:
    reference = _candidate("reference", "clicked", strand=1)
    target = _candidate("target", "target", strand=-1, reversed_=True)
    oriented = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, target),
        mode=SimilarityAlignmentMode.POSITION_AND_ORIENTATION,
    )
    position_only = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, target),
        mode=SimilarityAlignmentMode.POSITION,
    )
    unknown = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, _candidate("target", "target", strand=None)),
        mode=SimilarityAlignmentMode.POSITION_AND_ORIENTATION,
    )

    assert _decision(oriented, "target").effective_reverse_complement is False
    assert _decision(position_only, "target").effective_reverse_complement is None
    assert _decision(unknown, "target").effective_reverse_complement is None
    assert _decision(oriented, "reference").effective_reverse_complement is None


@pytest.mark.parametrize(
    ("reference_strand", "target_strand", "target_reversed", "expected"),
    (
        (1, 1, False, None),
        (-1, -1, True, None),
        (1, -1, False, True),
        (-1, 1, True, False),
        (None, 1, False, None),
        (1, None, False, None),
    ),
    ids=(
        "same-forward",
        "same-reverse",
        "opposite-to-reverse",
        "opposite-to-forward",
        "unknown-reference",
        "unknown-target",
    ),
)
def test_position_and_orientation_strand_matrix(
    reference_strand: int | None,
    target_strand: int | None,
    target_reversed: bool,
    expected: bool | None,
) -> None:
    reference = _candidate("reference", "clicked", strand=reference_strand)
    target = _candidate(
        "target",
        "target",
        strand=target_strand,
        reversed_=target_reversed,
    )

    resolution = resolve_similarity_alignment(
        record_keys=("reference", "target"),
        group_id=GROUP,
        reference=reference.anchor,
        candidates=(reference, target),
        mode=SimilarityAlignmentMode.POSITION_AND_ORIENTATION,
    )

    assert _decision(resolution, "target").effective_reverse_complement is expected


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


def test_duplicate_decisions_and_invalid_plan_enums_are_rejected() -> None:
    anchor = _anchor("record", "feature")
    reference_decision = AlignmentRecordDecision(
        "record", "reference", "reference", anchor
    )
    with pytest.raises(ValidationError, match="duplicate record decisions"):
        SimilarityAlignmentPlan(
            mode="position",
            group_id=GROUP,
            reference=anchor,
            records=(reference_decision, reference_decision),
        )
    with pytest.raises(ValidationError, match="alignment mode must be one of"):
        SimilarityAlignmentPlan(
            mode="smart",
            group_id=GROUP,
            reference=anchor,
            records=(reference_decision,),
        )
    with pytest.raises(ValidationError, match="exactly one matching"):
        SimilarityAlignmentPlan(
            mode="position",
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
        mode="position",
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


def test_plan_schema_requires_exact_integer_one() -> None:
    anchor = _anchor("record", "feature")
    decision = AlignmentRecordDecision("record", "reference", "reference", anchor)
    with pytest.raises(ValidationError, match="schema must be 1"):
        SimilarityAlignmentPlan(
            mode="position",
            group_id=GROUP,
            reference=anchor,
            records=(decision,),
            schema=1.0,  # type: ignore[arg-type]
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
