"""Typed, side-effect-free Similarity Group alignment resolution.

This module owns semantic anchor selection only.  Callers are responsible for
building validated candidate facts from their surface-specific inputs and for
turning a completed plan into render geometry.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from enum import Enum
from math import isfinite
from numbers import Real

from gbdraw.exceptions import ValidationError


SIMILARITY_ALIGNMENT_PLAN_SCHEMA = 2


class AlignmentOrientationPolicy(str, Enum):
    PRESERVE = "preserve"
    MATCH_REFERENCE = "match_reference"


class AlignmentOrientationEffect(str, Enum):
    PRESERVE = "preserve"
    REVERSE_WHOLE_RECORD = "reverse_whole_record"
    PRESERVE_UNKNOWN_STRAND = "preserve_unknown_strand"


class AlignmentRecommendationReason(str, Enum):
    UNIQUE_REPRESENTATIVE = "unique_representative"
    DETERMINISTIC_CANDIDATE_1 = "deterministic_candidate_1"


class AlignmentDecisionStatus(str, Enum):
    REFERENCE = "reference"
    ALIGNED = "aligned"
    SKIPPED = "skipped"


class AlignmentResolutionRationale(str, Enum):
    REFERENCE = "reference"
    USER_SELECTED = "user_selected"
    ONLY_USABLE_CANDIDATE = "only_usable_candidate"
    UNIQUE_DIRECT_RBH = "unique_direct_rbh"
    SKIPPED_BY_USER = "skipped_by_user"
    SKIPPED_NO_CANDIDATE = "skipped_no_candidate"
    SKIPPED_UNMAPPABLE = "skipped_unmappable"


class AlignmentChoiceKind(str, Enum):
    SELECT = "select"
    SKIP = "skip"


def _required_text(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\0" in value:
        raise ValidationError(f"{name} must be a non-empty string without NUL.")
    return value.strip()


def _optional_text(value: object, name: str) -> str | None:
    if value is None:
        return None
    return _required_text(value, name)


def _enum_value(value: object, enum_type: type[Enum], name: str) -> Enum:
    try:
        return enum_type(value)
    except (TypeError, ValueError) as exc:
        allowed = ", ".join(str(item.value) for item in enum_type)
        raise ValidationError(f"{name} must be one of: {allowed}.") from exc


def _record_order(record_keys: Sequence[str]) -> tuple[str, ...]:
    if isinstance(record_keys, (str, bytes)) or not isinstance(
        record_keys, Sequence
    ):
        raise ValidationError("record_keys must be a sequence of record identities.")
    normalized = tuple(
        _required_text(record_key, "record_key") for record_key in record_keys
    )
    if not normalized:
        raise ValidationError("Similarity alignment requires at least one record.")
    if len(set(normalized)) != len(normalized):
        raise ValidationError("Similarity alignment record keys must be unique.")
    return normalized


@dataclass(frozen=True)
class AlignmentAnchorIdentity:
    """Canonical biological feature identity plus optional consistency evidence."""

    record_key: str
    biological_feature_id: str
    source_feature_index: int | None = None
    stable_feature_svg_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "record_key", _required_text(self.record_key, "record_key")
        )
        object.__setattr__(
            self,
            "biological_feature_id",
            _required_text(self.biological_feature_id, "biological_feature_id"),
        )
        if self.source_feature_index is not None and (
            isinstance(self.source_feature_index, bool)
            or not isinstance(self.source_feature_index, int)
            or self.source_feature_index < 0
        ):
            raise ValidationError(
                "source_feature_index must be a non-negative integer or None."
            )
        object.__setattr__(
            self,
            "stable_feature_svg_id",
            _optional_text(self.stable_feature_svg_id, "stable_feature_svg_id"),
        )

    @property
    def canonical_key(self) -> tuple[str, str]:
        return self.record_key, self.biological_feature_id

    @property
    def sort_key(self) -> tuple[str, int, str]:
        return (
            self.biological_feature_id,
            -1 if self.source_feature_index is None else self.source_feature_index,
            self.stable_feature_svg_id or "",
        )


def _identities_agree(
    expected: AlignmentAnchorIdentity,
    actual: AlignmentAnchorIdentity,
) -> bool:
    if expected.canonical_key != actual.canonical_key:
        return False
    return not (
        expected.source_feature_index is not None
        and actual.source_feature_index is not None
        and expected.source_feature_index != actual.source_feature_index
        or expected.stable_feature_svg_id is not None
        and actual.stable_feature_svg_id is not None
        and expected.stable_feature_svg_id != actual.stable_feature_svg_id
    )


@dataclass(frozen=True)
class SimilarityAlignmentCandidate:
    """One group member projected into the current displayed-record domain."""

    group_id: str
    anchor: AlignmentAnchorIdentity
    displayed_strand: int | None
    center_mappable: bool
    display_center: float | None = None
    identity_is_unique: bool = True
    hidden: bool = False
    effective_reverse_complement: bool = False
    representative: bool = False
    role: str = ""
    source_start: int | None = None
    source_end: int | None = None
    display_name: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "group_id", _required_text(self.group_id, "group_id")
        )
        if not isinstance(self.anchor, AlignmentAnchorIdentity):
            raise ValidationError("Candidate anchor must be AlignmentAnchorIdentity.")
        if self.displayed_strand is not None and (
            isinstance(self.displayed_strand, bool)
            or not isinstance(self.displayed_strand, int)
            or self.displayed_strand not in (-1, 1)
        ):
            raise ValidationError("displayed_strand must be -1, 1, or None.")
        for name in (
            "center_mappable",
            "identity_is_unique",
            "hidden",
            "effective_reverse_complement",
            "representative",
        ):
            if not isinstance(getattr(self, name), bool):
                raise ValidationError(f"{name} must be a boolean.")
        if self.display_center is not None:
            if (
                isinstance(self.display_center, bool)
                or not isinstance(self.display_center, Real)
                or not isfinite(float(self.display_center))
            ):
                raise ValidationError("display_center must be a finite number or None.")
            object.__setattr__(self, "display_center", float(self.display_center))
        if self.center_mappable != (self.display_center is not None):
            raise ValidationError(
                "A mappable candidate requires one finite display center; an "
                "unmappable candidate must omit it."
            )
        if not isinstance(self.role, str) or "\0" in self.role:
            raise ValidationError("Candidate role must be text without NUL.")
        object.__setattr__(self, "role", self.role.strip())
        if (self.source_start is None) != (self.source_end is None):
            raise ValidationError("Candidate source coordinates must both be supplied or omitted.")
        if self.source_start is not None and (
            isinstance(self.source_start, bool)
            or not isinstance(self.source_start, int)
            or self.source_start < 0
            or isinstance(self.source_end, bool)
            or not isinstance(self.source_end, int)
            or self.source_end < self.source_start
        ):
            raise ValidationError("Candidate source coordinates must be ordered non-negative integers.")
        object.__setattr__(
            self, "display_name", _optional_text(self.display_name, "display_name")
        )

    def is_usable_for(self, group_id: str) -> bool:
        return (
            self.group_id == group_id
            and self.identity_is_unique
            and self.center_mappable
        )


@dataclass(frozen=True)
class AlignmentEvidenceEdge:
    """One direct evidence edge; only exact-reference RBH affects selection."""

    group_id: str
    query: AlignmentAnchorIdentity
    subject: AlignmentAnchorIdentity
    edge_kind: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "group_id", _required_text(self.group_id, "group_id")
        )
        if not isinstance(self.query, AlignmentAnchorIdentity) or not isinstance(
            self.subject, AlignmentAnchorIdentity
        ):
            raise ValidationError(
                "Evidence endpoints must be AlignmentAnchorIdentity values."
            )
        object.__setattr__(
            self, "edge_kind", _required_text(self.edge_kind, "edge_kind").lower()
        )


@dataclass(frozen=True)
class AlignmentRecordChoice:
    """An explicit Select or Skip answer for one non-reference record."""

    record_key: str
    kind: AlignmentChoiceKind
    anchor: AlignmentAnchorIdentity | None = None
    orientation_policy: AlignmentOrientationPolicy = AlignmentOrientationPolicy.PRESERVE

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "record_key", _required_text(self.record_key, "record_key")
        )
        object.__setattr__(
            self,
            "orientation_policy",
            _enum_value(
                self.orientation_policy, AlignmentOrientationPolicy, "orientation policy"
            ),
        )
        object.__setattr__(
            self,
            "kind",
            _enum_value(self.kind, AlignmentChoiceKind, "choice kind"),
        )
        if self.kind is AlignmentChoiceKind.SELECT:
            if not isinstance(self.anchor, AlignmentAnchorIdentity):
                raise ValidationError("A Select choice requires an anchor.")
            if self.anchor.record_key != self.record_key:
                raise ValidationError(
                    "A selected anchor must belong to the choice record."
                )
        elif self.anchor is not None:
            raise ValidationError("A Skip choice must not contain an anchor.")
        elif self.orientation_policy is not AlignmentOrientationPolicy.PRESERVE:
            raise ValidationError("A Skip choice must preserve orientation.")


@dataclass(frozen=True)
class AlignmentRecordDecision:
    """A completed reference, aligned, or unchanged/skipped record decision."""

    record_key: str
    status: AlignmentDecisionStatus
    rationale: AlignmentResolutionRationale
    anchor: AlignmentAnchorIdentity | None = None
    orientation_policy: AlignmentOrientationPolicy = AlignmentOrientationPolicy.PRESERVE
    effective_reverse_complement: bool | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "record_key", _required_text(self.record_key, "record_key")
        )
        object.__setattr__(
            self,
            "status",
            _enum_value(self.status, AlignmentDecisionStatus, "decision status"),
        )
        object.__setattr__(
            self,
            "rationale",
            _enum_value(
                self.rationale, AlignmentResolutionRationale, "decision rationale"
            ),
        )
        object.__setattr__(
            self,
            "orientation_policy",
            _enum_value(
                self.orientation_policy, AlignmentOrientationPolicy, "orientation policy"
            ),
        )
        if self.effective_reverse_complement is not None and not isinstance(
            self.effective_reverse_complement, bool
        ):
            raise ValidationError(
                "effective_reverse_complement must be a boolean or None."
            )

        aligned_rationales = {
            AlignmentResolutionRationale.USER_SELECTED,
            AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
            AlignmentResolutionRationale.UNIQUE_DIRECT_RBH,
        }
        skipped_rationales = {
            AlignmentResolutionRationale.SKIPPED_BY_USER,
            AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE,
            AlignmentResolutionRationale.SKIPPED_UNMAPPABLE,
        }
        if self.status is AlignmentDecisionStatus.REFERENCE:
            valid = (
                isinstance(self.anchor, AlignmentAnchorIdentity)
                and self.anchor.record_key == self.record_key
                and self.rationale is AlignmentResolutionRationale.REFERENCE
                and self.effective_reverse_complement is None
                and self.orientation_policy is AlignmentOrientationPolicy.PRESERVE
            )
        elif self.status is AlignmentDecisionStatus.ALIGNED:
            valid = (
                isinstance(self.anchor, AlignmentAnchorIdentity)
                and self.anchor.record_key == self.record_key
                and self.rationale in aligned_rationales
                and isinstance(self.effective_reverse_complement, bool)
            )
        else:
            valid = (
                self.anchor is None
                and self.rationale in skipped_rationales
                and self.effective_reverse_complement is None
                and self.orientation_policy is AlignmentOrientationPolicy.PRESERVE
            )
        if not valid:
            raise ValidationError(
                "Alignment decision status, rationale, anchor, and orientation "
                "policy/effect are inconsistent."
            )

    @property
    def review_reason(self) -> AlignmentResolutionRationale | None:
        if self.rationale in (
            AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
            AlignmentResolutionRationale.UNIQUE_DIRECT_RBH,
        ):
            return self.rationale
        return None


@dataclass(frozen=True)
class AmbiguousAlignmentRecord:
    """A record whose remaining usable members require Select or Skip."""

    record_key: str
    candidates: tuple[SimilarityAlignmentCandidate, ...]
    direct_rbh_candidates: tuple[AlignmentAnchorIdentity, ...] = ()
    recommended_anchor: AlignmentAnchorIdentity = field(init=False)
    recommendation_reason: AlignmentRecommendationReason = field(init=False)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "record_key", _required_text(self.record_key, "record_key")
        )
        candidates = tuple(self.candidates)
        if len(candidates) < 2 or any(
            not isinstance(candidate, SimilarityAlignmentCandidate)
            or candidate.anchor.record_key != self.record_key
            or not candidate.identity_is_unique
            or not candidate.center_mappable
            for candidate in candidates
        ):
            raise ValidationError(
                "An ambiguous record requires at least two usable candidates "
                "from that record."
            )
        if len({candidate.anchor.canonical_key for candidate in candidates}) != len(
            candidates
        ):
            raise ValidationError("Ambiguous candidates must have unique identities.")
        candidates = tuple(sorted(candidates, key=lambda item: item.anchor.sort_key))
        object.__setattr__(self, "candidates", candidates)

        direct = tuple(self.direct_rbh_candidates)
        if any(not isinstance(anchor, AlignmentAnchorIdentity) for anchor in direct):
            raise ValidationError(
                "Direct-RBH candidate identities must be AlignmentAnchorIdentity values."
            )
        candidate_keys = {candidate.anchor.canonical_key for candidate in candidates}
        if (
            len({anchor.canonical_key for anchor in direct}) != len(direct)
            or any(anchor.canonical_key not in candidate_keys for anchor in direct)
        ):
            raise ValidationError(
                "Direct-RBH identities must be distinct ambiguous candidates."
            )
        candidate_by_key = {
            candidate.anchor.canonical_key: candidate for candidate in candidates
        }
        if any(
            not _identities_agree(anchor, candidate_by_key[anchor.canonical_key].anchor)
            for anchor in direct
        ):
            raise ValidationError(
                "Direct-RBH identity aliases conflict with ambiguous candidates."
            )
        object.__setattr__(
            self,
            "direct_rbh_candidates",
            tuple(sorted(direct, key=lambda item: item.sort_key)),
        )
        representatives = tuple(item for item in candidates if item.representative)
        if len(representatives) == 1:
            recommended = representatives[0]
            reason = AlignmentRecommendationReason.UNIQUE_REPRESENTATIVE
        else:
            recommended = candidates[0]
            reason = AlignmentRecommendationReason.DETERMINISTIC_CANDIDATE_1
        object.__setattr__(self, "recommended_anchor", recommended.anchor)
        object.__setattr__(self, "recommendation_reason", reason)


@dataclass(frozen=True)
class AlignmentReviewCandidate:
    """Python-owned facts for one candidate in a local review row."""

    candidate: SimilarityAlignmentCandidate
    usable: bool
    direct_evidence: tuple[str, ...]
    match_reference_effective_reverse_complement: bool
    match_reference_effect: AlignmentOrientationEffect


@dataclass(frozen=True)
class AlignmentReviewRow:
    record_key: str
    candidates: tuple[AlignmentReviewCandidate, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "record_key", _required_text(self.record_key, "record_key")
        )
        rows = tuple(self.candidates)
        if any(
            not isinstance(row, AlignmentReviewCandidate)
            or row.candidate.anchor.record_key != self.record_key
            for row in rows
        ):
            raise ValidationError("Review candidates must belong to their record.")
        object.__setattr__(self, "candidates", rows)


@dataclass(frozen=True)
class SimilarityAlignmentPlan:
    """Immutable, fully resolved Similarity Group alignment intent."""

    group_id: str
    reference: AlignmentAnchorIdentity
    records: tuple[AlignmentRecordDecision, ...]
    schema: int = SIMILARITY_ALIGNMENT_PLAN_SCHEMA

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema, bool)
            or not isinstance(self.schema, int)
            or self.schema != SIMILARITY_ALIGNMENT_PLAN_SCHEMA
        ):
            raise ValidationError(
                f"Similarity alignment plan schema must be {SIMILARITY_ALIGNMENT_PLAN_SCHEMA}."
            )
        object.__setattr__(self, "group_id", _required_text(self.group_id, "group_id"))
        if not isinstance(self.reference, AlignmentAnchorIdentity):
            raise ValidationError("Plan reference must be AlignmentAnchorIdentity.")
        records = tuple(self.records)
        if not records or any(
            not isinstance(record, AlignmentRecordDecision) for record in records
        ):
            raise ValidationError(
                "Plan records must contain AlignmentRecordDecision values."
            )
        keys = [record.record_key for record in records]
        if len(set(keys)) != len(keys):
            raise ValidationError("Plan contains duplicate record decisions.")
        references = [
            record
            for record in records
            if record.status is AlignmentDecisionStatus.REFERENCE
        ]
        if (
            len(references) != 1
            or references[0].anchor != self.reference
            or references[0].record_key != self.reference.record_key
        ):
            raise ValidationError(
                "Plan reference must have exactly one matching reference decision."
            )
        object.__setattr__(self, "records", records)

    def validate_record_coverage(self, record_keys: Sequence[str]) -> None:
        expected = set(_record_order(record_keys))
        actual = {record.record_key for record in self.records}
        if actual != expected:
            missing = sorted(expected - actual)
            unknown = sorted(actual - expected)
            raise ValidationError(
                "Similarity alignment plan record coverage differs from the "
                f"displayed records (missing={missing}, unknown={unknown})."
            )


AlignmentResolutionRecord = AlignmentRecordDecision | AmbiguousAlignmentRecord


@dataclass(frozen=True)
class SimilarityAlignmentResolution:
    """Per-record outcome and review facts, possibly awaiting explicit choices."""

    group_id: str
    reference: AlignmentAnchorIdentity
    reference_candidate: SimilarityAlignmentCandidate
    records: tuple[AlignmentResolutionRecord, ...]
    review_rows: tuple[AlignmentReviewRow, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "group_id", _required_text(self.group_id, "group_id"))
        if not isinstance(self.reference, AlignmentAnchorIdentity):
            raise ValidationError("Resolution reference must be AlignmentAnchorIdentity.")
        if (
            not isinstance(self.reference_candidate, SimilarityAlignmentCandidate)
            or self.reference_candidate.anchor != self.reference
            or not self.reference_candidate.is_usable_for(self.group_id)
        ):
            raise ValidationError("Resolution reference facts are inconsistent.")
        records = tuple(self.records)
        if not records or any(
            not isinstance(record, (AlignmentRecordDecision, AmbiguousAlignmentRecord))
            for record in records
        ):
            raise ValidationError("Resolution records must contain decisions or ambiguities.")
        keys = [record.record_key for record in records]
        if len(set(keys)) != len(keys):
            raise ValidationError("Resolution contains duplicate record outcomes.")
        references = [
            record for record in records
            if isinstance(record, AlignmentRecordDecision)
            and record.status is AlignmentDecisionStatus.REFERENCE
        ]
        if len(references) != 1 or references[0].anchor != self.reference:
            raise ValidationError("Resolution requires exactly one matching reference outcome.")
        rows = tuple(self.review_rows)
        if (
            len(rows) != len(records)
            or any(not isinstance(row, AlignmentReviewRow) for row in rows)
            or tuple(row.record_key for row in rows) != tuple(keys)
            or any(row.candidates for row in rows if row.record_key == self.reference.record_key)
        ):
            raise ValidationError("Resolution review rows must cover each record in order.")
        if any(
            isinstance(record, AmbiguousAlignmentRecord)
            and any(candidate.group_id != self.group_id for candidate in record.candidates)
            for record in records
        ):
            raise ValidationError("Ambiguous candidates must belong to the selected group.")
        object.__setattr__(self, "records", records)
        object.__setattr__(self, "review_rows", rows)

    @property
    def ambiguities(self) -> tuple[AmbiguousAlignmentRecord, ...]:
        return tuple(
            record for record in self.records
            if isinstance(record, AmbiguousAlignmentRecord)
        )

    @property
    def decisions(self) -> tuple[AlignmentRecordDecision, ...]:
        return tuple(
            record for record in self.records
            if isinstance(record, AlignmentRecordDecision)
        )

    @property
    def plan(self) -> SimilarityAlignmentPlan | None:
        if self.ambiguities:
            return None
        return SimilarityAlignmentPlan(
            group_id=self.group_id,
            reference=self.reference,
            records=self.decisions,
        )

    def require_plan(self) -> SimilarityAlignmentPlan:
        plan = self.plan
        if plan is None:
            unresolved = ", ".join(record.record_key for record in self.ambiguities)
            raise ValidationError(
                "Similarity alignment requires Select or Skip for records: "
                f"{unresolved}."
            )
        return plan


def _orientation_result(
    policy: AlignmentOrientationPolicy,
    reference: SimilarityAlignmentCandidate,
    target: SimilarityAlignmentCandidate,
) -> tuple[bool, AlignmentOrientationEffect]:
    base = target.effective_reverse_complement
    if policy is AlignmentOrientationPolicy.PRESERVE:
        return base, AlignmentOrientationEffect.PRESERVE
    if reference.displayed_strand is None or target.displayed_strand is None:
        return base, AlignmentOrientationEffect.PRESERVE_UNKNOWN_STRAND
    if reference.displayed_strand != target.displayed_strand:
        return not base, AlignmentOrientationEffect.REVERSE_WHOLE_RECORD
    return base, AlignmentOrientationEffect.PRESERVE


def _aligned_decision(
    candidate: SimilarityAlignmentCandidate,
    rationale: AlignmentResolutionRationale,
    *,
    policy: AlignmentOrientationPolicy,
    reference: SimilarityAlignmentCandidate,
) -> AlignmentRecordDecision:
    effective, _effect = _orientation_result(policy, reference, candidate)
    return AlignmentRecordDecision(
        record_key=candidate.anchor.record_key,
        status=AlignmentDecisionStatus.ALIGNED,
        rationale=rationale,
        anchor=candidate.anchor,
        orientation_policy=policy,
        effective_reverse_complement=effective,
    )


def _validate_candidate_aliases(
    candidates: tuple[SimilarityAlignmentCandidate, ...],
) -> dict[tuple[str, str], SimilarityAlignmentCandidate]:
    by_identity: dict[tuple[str, str], SimilarityAlignmentCandidate] = {}
    by_source_index: dict[tuple[str, int], tuple[str, str]] = {}
    for candidate in candidates:
        key = candidate.anchor.canonical_key
        if key in by_identity:
            raise ValidationError(
                f"Duplicate alignment candidate identity {key!r}."
            )
        by_identity[key] = candidate
        if candidate.anchor.source_feature_index is None:
            continue
        source_key = (
            candidate.anchor.record_key,
            candidate.anchor.source_feature_index,
        )
        existing = by_source_index.get(source_key)
        if existing is not None and existing != key:
            raise ValidationError(
                "Alignment candidates contain conflicting source-feature aliases."
            )
        by_source_index[source_key] = key
    return by_identity


def _validate_known_alias(
    identity: AlignmentAnchorIdentity,
    candidates: dict[tuple[str, str], SimilarityAlignmentCandidate],
    *,
    description: str,
) -> None:
    candidate = candidates.get(identity.canonical_key)
    if candidate is not None and not _identities_agree(identity, candidate.anchor):
        raise ValidationError(f"{description} conflicts with candidate identity aliases.")


def _direct_rbh_candidates(
    *,
    group_id: str,
    reference: AlignmentAnchorIdentity,
    candidates: tuple[SimilarityAlignmentCandidate, ...],
    edges: tuple[AlignmentEvidenceEdge, ...],
) -> tuple[SimilarityAlignmentCandidate, ...]:
    candidates_by_key = {
        candidate.anchor.canonical_key: candidate for candidate in candidates
    }
    connected: set[tuple[str, str]] = set()
    for edge in edges:
        if edge.group_id != group_id or edge.edge_kind != "rbh":
            continue
        if edge.query.canonical_key == reference.canonical_key:
            other = edge.subject.canonical_key
        elif edge.subject.canonical_key == reference.canonical_key:
            other = edge.query.canonical_key
        else:
            continue
        if other in candidates_by_key:
            connected.add(other)
    return tuple(
        sorted(
            (candidates_by_key[key] for key in connected),
            key=lambda item: item.anchor.sort_key,
        )
    )


def _direct_evidence_kinds(
    *,
    group_id: str,
    reference: AlignmentAnchorIdentity,
    candidate: AlignmentAnchorIdentity,
    edges: tuple[AlignmentEvidenceEdge, ...],
) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                edge.edge_kind
                for edge in edges
                if edge.group_id == group_id
                and (
                    (
                        edge.query.canonical_key == reference.canonical_key
                        and edge.subject.canonical_key == candidate.canonical_key
                    )
                    or (
                        edge.subject.canonical_key == reference.canonical_key
                        and edge.query.canonical_key == candidate.canonical_key
                    )
                )
            }
        )
    )


def _review_row(
    *,
    record_key: str,
    group_id: str,
    reference: SimilarityAlignmentCandidate,
    candidates: tuple[SimilarityAlignmentCandidate, ...],
    edges: tuple[AlignmentEvidenceEdge, ...],
) -> AlignmentReviewRow:
    rows = []
    for candidate in candidates:
        matched, effect = _orientation_result(
            AlignmentOrientationPolicy.MATCH_REFERENCE, reference, candidate
        )
        rows.append(
            AlignmentReviewCandidate(
                candidate=candidate,
                usable=candidate.is_usable_for(group_id),
                direct_evidence=_direct_evidence_kinds(
                    group_id=group_id,
                    reference=reference.anchor,
                    candidate=candidate.anchor,
                    edges=edges,
                ),
                match_reference_effective_reverse_complement=matched,
                match_reference_effect=effect,
            )
        )
    return AlignmentReviewRow(record_key=record_key, candidates=tuple(rows))


def resolve_similarity_alignment(
    *,
    record_keys: Sequence[str],
    group_id: str,
    reference: AlignmentAnchorIdentity,
    candidates: Sequence[SimilarityAlignmentCandidate],
    edges: Sequence[AlignmentEvidenceEdge] = (),
    choices: Sequence[AlignmentRecordChoice] = (),
) -> SimilarityAlignmentResolution:
    """Resolve one exact anchor decision for every displayed record.

    Priority is explicit Select/Skip, only usable candidate, one distinct
    direct RBH candidate, ambiguity, then unchanged when none is usable.
    """

    order = _record_order(record_keys)
    known_records = set(order)
    resolved_group_id = _required_text(group_id, "group_id")
    if not isinstance(reference, AlignmentAnchorIdentity):
        raise ValidationError("reference must be AlignmentAnchorIdentity.")
    if reference.record_key not in known_records:
        raise ValidationError("Reference record is not a displayed record.")

    if isinstance(candidates, (str, bytes)) or not isinstance(candidates, Sequence):
        raise ValidationError("candidates must be a sequence.")
    candidate_values = tuple(candidates)
    if any(
        not isinstance(candidate, SimilarityAlignmentCandidate)
        for candidate in candidate_values
    ):
        raise ValidationError(
            "candidates must contain SimilarityAlignmentCandidate values."
        )
    unknown_candidate_records = sorted(
        {
            candidate.anchor.record_key
            for candidate in candidate_values
            if candidate.anchor.record_key not in known_records
        }
    )
    if unknown_candidate_records:
        raise ValidationError(
            "Alignment candidates reference unknown displayed records: "
            f"{unknown_candidate_records}."
        )
    candidates_by_identity = _validate_candidate_aliases(candidate_values)

    reference_candidate = candidates_by_identity.get(reference.canonical_key)
    if (
        reference_candidate is None
        or reference_candidate.group_id != resolved_group_id
    ):
        raise ValidationError(
            "The exact reference is not a member of the selected Similarity Group."
        )
    if not _identities_agree(reference, reference_candidate.anchor):
        raise ValidationError("Reference identity aliases conflict with group metadata.")
    if not reference_candidate.is_usable_for(resolved_group_id):
        raise ValidationError(
            "The exact reference does not resolve uniquely inside the current crop."
        )

    if isinstance(edges, (str, bytes)) or not isinstance(edges, Sequence):
        raise ValidationError("edges must be a sequence.")
    edge_values = tuple(edges)
    if any(not isinstance(edge, AlignmentEvidenceEdge) for edge in edge_values):
        raise ValidationError("edges must contain AlignmentEvidenceEdge values.")
    for edge in edge_values:
        _validate_known_alias(
            edge.query,
            candidates_by_identity,
            description="Evidence query endpoint",
        )
        _validate_known_alias(
            edge.subject,
            candidates_by_identity,
            description="Evidence subject endpoint",
        )

    if isinstance(choices, (str, bytes)) or not isinstance(choices, Sequence):
        raise ValidationError("choices must be a sequence.")
    choice_values = tuple(choices)
    if any(not isinstance(choice, AlignmentRecordChoice) for choice in choice_values):
        raise ValidationError("choices must contain AlignmentRecordChoice values.")
    choice_by_record: dict[str, AlignmentRecordChoice] = {}
    for choice in choice_values:
        if choice.record_key not in known_records:
            raise ValidationError(
                f"Explicit choice references unknown record {choice.record_key!r}."
            )
        if choice.record_key == reference_candidate.anchor.record_key:
            raise ValidationError("The reference record cannot have a target choice.")
        if choice.record_key in choice_by_record:
            raise ValidationError(
                f"Duplicate explicit choice for record {choice.record_key!r}."
            )
        if choice.anchor is not None:
            _validate_known_alias(
                choice.anchor,
                candidates_by_identity,
                description="Explicit choice",
            )
        choice_by_record[choice.record_key] = choice

    records: list[AlignmentResolutionRecord] = []
    review_rows: list[AlignmentReviewRow] = []
    for record_key in order:
        if record_key == reference_candidate.anchor.record_key:
            review_rows.append(AlignmentReviewRow(record_key, ()))
            records.append(
                AlignmentRecordDecision(
                    record_key=record_key,
                    status=AlignmentDecisionStatus.REFERENCE,
                    rationale=AlignmentResolutionRationale.REFERENCE,
                    anchor=reference_candidate.anchor,
                )
            )
            continue

        group_candidates = tuple(
            sorted(
                (
                    candidate
                    for candidate in candidate_values
                    if candidate.anchor.record_key == record_key
                    and candidate.group_id == resolved_group_id
                ),
                key=lambda item: item.anchor.sort_key,
            )
        )
        review_rows.append(
            _review_row(
                record_key=record_key,
                group_id=resolved_group_id,
                reference=reference_candidate,
                candidates=group_candidates,
                edges=edge_values,
            )
        )
        usable = tuple(
            candidate
            for candidate in group_candidates
            if candidate.is_usable_for(resolved_group_id)
        )
        choice = choice_by_record.get(record_key)
        if choice is not None:
            if choice.kind is AlignmentChoiceKind.SKIP:
                records.append(
                    AlignmentRecordDecision(
                        record_key=record_key,
                        status=AlignmentDecisionStatus.SKIPPED,
                        rationale=AlignmentResolutionRationale.SKIPPED_BY_USER,
                    )
                )
                continue
            selected = next(
                (
                    candidate
                    for candidate in usable
                    if candidate.anchor.canonical_key == choice.anchor.canonical_key
                    and _identities_agree(choice.anchor, candidate.anchor)
                ),
                None,
            )
            if selected is None:
                raise ValidationError(
                    f"Explicit choice for record {record_key!r} is stale or unusable."
                )
            records.append(
                _aligned_decision(
                    selected,
                    AlignmentResolutionRationale.USER_SELECTED,
                    policy=choice.orientation_policy,
                    reference=reference_candidate,
                )
            )
            continue

        if not usable:
            records.append(
                AlignmentRecordDecision(
                    record_key=record_key,
                    status=AlignmentDecisionStatus.SKIPPED,
                    rationale=(
                        AlignmentResolutionRationale.SKIPPED_UNMAPPABLE
                        if group_candidates
                        else AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE
                    ),
                )
            )
            continue
        if len(usable) == 1:
            records.append(
                _aligned_decision(
                    usable[0],
                    AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE,
                    policy=AlignmentOrientationPolicy.PRESERVE,
                    reference=reference_candidate,
                )
            )
            continue

        direct_rbh = _direct_rbh_candidates(
            group_id=resolved_group_id,
            reference=reference_candidate.anchor,
            candidates=usable,
            edges=edge_values,
        )
        if len(direct_rbh) == 1:
            records.append(
                _aligned_decision(
                    direct_rbh[0],
                    AlignmentResolutionRationale.UNIQUE_DIRECT_RBH,
                    policy=AlignmentOrientationPolicy.PRESERVE,
                    reference=reference_candidate,
                )
            )
            continue
        records.append(
            AmbiguousAlignmentRecord(
                record_key=record_key,
                candidates=usable,
                direct_rbh_candidates=tuple(
                    candidate.anchor for candidate in direct_rbh
                ),
            )
        )

    resolution = SimilarityAlignmentResolution(
        group_id=resolved_group_id,
        reference=reference_candidate.anchor,
        reference_candidate=reference_candidate,
        records=tuple(records),
        review_rows=tuple(review_rows),
    )
    actual_order = tuple(record.record_key for record in resolution.records)
    if actual_order != order:
        raise AssertionError("Resolver output did not preserve canonical record order.")
    return resolution


__all__ = [
    "AlignmentAnchorIdentity",
    "AlignmentChoiceKind",
    "AlignmentDecisionStatus",
    "AlignmentEvidenceEdge",
    "AlignmentRecordChoice",
    "AlignmentRecordDecision",
    "AlignmentOrientationEffect",
    "AlignmentOrientationPolicy",
    "AlignmentRecommendationReason",
    "AlignmentResolutionRationale",
    "AlignmentReviewCandidate",
    "AlignmentReviewRow",
    "AmbiguousAlignmentRecord",
    "SIMILARITY_ALIGNMENT_PLAN_SCHEMA",
    "SimilarityAlignmentCandidate",
    "SimilarityAlignmentPlan",
    "SimilarityAlignmentResolution",
    "resolve_similarity_alignment",
]
