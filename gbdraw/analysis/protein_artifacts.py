"""Current derived protein artifacts: their identity, key and payload, and
the Session-neutral validation of stored ones."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import re
from typing import TYPE_CHECKING, Literal, Mapping, Sequence, TypeGuard, Union, cast

from gbdraw.exceptions import ValidationError

if TYPE_CHECKING:
    from pandas import DataFrame

    from .collinearity import CollinearityResult, LosslessCollinearityParameters
    from .protein_colinearity import ProteinIdentityManifest

CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA = 3

_RUNTIME_HANDLE_RE = re.compile(r"h_[a-z2-7]{26}")
_FEATURE_ANALYSIS_ID_RE = re.compile(r"(?<![A-Za-z0-9_])f_[0-9a-f]{64}(?![A-Za-z0-9_])")
_LEGACY_PROTEIN_REFERENCE_RE = re.compile(
    r"(?<![A-Za-z0-9._%+-])"
    r"p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}"
    r"(?:_[2-9][0-9]*)?"
    r"(?![A-Za-z0-9._%+-])"
)
_SCALAR_REFERENCE_KEYS = frozenset(
    {
        "proteinId",
        "queryProteinId",
        "subjectProteinId",
        "protein_id",
        "query_protein_id",
        "subject_protein_id",
    }
)
_UNIT_REFERENCE_KEYS = frozenset(
    {
        "queryUnitId",
        "subjectUnitId",
        "query_unit_id",
        "subject_unit_id",
    }
)
_ARRAY_REFERENCE_KEYS = frozenset(
    {
        "proteinIds",
        "sharedProteinIds",
        "protein_ids",
        "shared_protein_ids",
    }
)
_COMPOUND_EDGE_REFERENCE_KEYS = frozenset(
    {
        "supportingEdge",
        "supportingEdges",
        "supporting_edge",
        "supporting_edges",
        "edgeId",
        "edgeIds",
        "edge_id",
        "edge_ids",
    }
)
_COMPOUND_SUPPORTING_EDGE_RE = re.compile(
    rf"^(?P<query>{_RUNTIME_HANDLE_RE.pattern})->"
    rf"(?P<subject>{_RUNTIME_HANDLE_RE.pattern}):"
    r"[A-Za-z][A-Za-z0-9._-]*$"
)
_COMPOUND_PATH_EDGE_RE = re.compile(
    rf"^[^:\s]+:\d+:(?P<query>{_RUNTIME_HANDLE_RE.pattern})->"
    rf"\d+:(?P<subject>{_RUNTIME_HANDLE_RE.pattern}):"
    r"[A-Za-z][A-Za-z0-9._-]*$"
)


def is_current_derived_protein_artifact(entry: object) -> bool:
    """Return whether an object has the current derived artifact envelope."""

    return (
        isinstance(entry, Mapping)
        and entry.get("schema") == CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA
        and entry.get("kind") == "derived-losatp-payload"
        and entry.get("idEncoding") == "runtime-handle-v1"
        and isinstance(entry.get("key"), str)
        and bool(entry.get("key"))
        and isinstance(entry.get("payload"), Mapping)
    )


def _is_nonnegative_integer(value: object) -> TypeGuard[int]:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _has_one_empty_table_form(pair: Mapping[str, object]) -> bool:
    """A zero-hit pair holds exactly one empty table: the header-only
    ``canonicalTsv`` the current writer emits, or the ``rows == []`` that
    derived entries carried from #312 until the canonical table replaced it."""

    has_table = "canonicalTsv" in pair
    has_rows = "rows" in pair
    if has_table == has_rows:
        return False
    if has_rows:
        return pair.get("rows") == []
    table = pair.get("canonicalTsv")
    return isinstance(table, str) and table.endswith("\n") and table.count("\n") == 1


def _is_strict_empty_result(entry: Mapping[str, object]) -> bool:
    mode = entry.get("mode")
    if mode not in {"orthogroup", "collinear"}:
        return False
    payload = entry.get("payload")
    if not isinstance(payload, Mapping):
        return False

    allowed_keys = {"identity", "provenance", "pairs", "orthogroups"}
    if mode == "collinear":
        allowed_keys.update(
            {
                "collinearGroups",
                "collinearGroupScope",
                "collinearityBlocks",
            }
        )
    if not set(payload).issubset(allowed_keys):
        return False

    if "identity" in payload:
        identity = payload["identity"]
        if (
            not isinstance(identity, Mapping)
            or identity.get("cacheSchema") != CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA
            or identity.get("idEncoding") != "runtime-handle-v1"
            or identity.get("mode") != mode
        ):
            return False
        raw_cache_keys = identity.get("rawCacheKeys")
        if not isinstance(raw_cache_keys, list) or not all(
            isinstance(key, str) and bool(key) for key in raw_cache_keys
        ):
            return False

    pairs = payload.get("pairs")
    orthogroups = payload.get("orthogroups")
    if not isinstance(pairs, list) or not pairs or orthogroups != []:
        return False
    pair_indices: set[int] = set()
    for pair in pairs:
        if not isinstance(pair, Mapping) or set(pair).difference(
            {
                "pair_index",
                "query_index",
                "subject_index",
                "tsv",
                "canonicalTsv",
                "rows",
                "hit_count",
            }
        ):
            return False
        pair_index = pair.get("pair_index")
        if (
            not _is_nonnegative_integer(pair_index)
            or pair_index in pair_indices
            or pair.get("tsv") != ""
            or not _has_one_empty_table_form(pair)
            or not _is_nonnegative_integer(pair.get("hit_count"))
            or pair.get("hit_count") != 0
        ):
            return False
        pair_indices.add(pair_index)
        has_query_index = "query_index" in pair
        has_subject_index = "subject_index" in pair
        if has_query_index != has_subject_index:
            return False
        if has_query_index and (
            not _is_nonnegative_integer(pair.get("query_index"))
            or not _is_nonnegative_integer(pair.get("subject_index"))
        ):
            return False

    if mode == "collinear":
        for collection_key in ("collinearGroups", "collinearityBlocks"):
            if collection_key in payload and payload[collection_key] != []:
                return False
        if "collinearGroupScope" in payload and payload["collinearGroupScope"] not in {
            "adjacent_local",
            "global_collinear",
        }:
            return False
    return True


def validate_current_derived_protein_artifacts(
    entries: Sequence[Mapping[str, object]],
    manifest: ProteinIdentityManifest | Mapping[str, object] | None,
) -> None:
    """Require current envelopes and resolve every protein reference."""

    if manifest is None:
        raise ValidationError(
            "Current derived protein artifacts require protein_identity_manifest."
        )
    from .protein_colinearity import (
        ProteinIdentityManifest,
        validate_protein_identity_manifest,
    )

    try:
        authority = (
            manifest
            if isinstance(manifest, ProteinIdentityManifest)
            else validate_protein_identity_manifest(manifest)
        )
    except (TypeError, ValueError, ValidationError) as exc:
        raise ValidationError(
            "Current derived protein artifacts require a valid "
            "protein_identity_manifest."
        ) from exc

    runtime_handles = {
        str(handle)
        for binding in authority.record_instances.values()
        for handle in (
            cast(Mapping[str, object], binding.get("runtimeIds", {})).values()
            if isinstance(binding.get("runtimeIds"), Mapping)
            else ()
        )
    }

    def forbidden_reference(value: str) -> bool:
        return (
            _LEGACY_PROTEIN_REFERENCE_RE.search(value) is not None
            or _FEATURE_ANALYSIS_ID_RE.search(value) is not None
        )

    def compound_edge_references(value: str) -> tuple[str, str] | None:
        match = _COMPOUND_SUPPORTING_EDGE_RE.fullmatch(
            value
        ) or _COMPOUND_PATH_EDGE_RE.fullmatch(value)
        if match is None:
            return None
        return match.group("query"), match.group("subject")

    def visit(value: object, owner_key: str = "") -> tuple[bool, bool]:
        if owner_key in _COMPOUND_EDGE_REFERENCE_KEYS and not isinstance(
            value, (str, list)
        ):
            return False, False
        if isinstance(value, str):
            if forbidden_reference(value):
                return False, False
            if owner_key in _COMPOUND_EDGE_REFERENCE_KEYS:
                references = compound_edge_references(value)
                return (
                    references is not None
                    and all(reference in runtime_handles for reference in references),
                    references is not None,
                )
            if owner_key in _SCALAR_REFERENCE_KEYS and value:
                scalar_references = [reference.strip() for reference in value.split(";")]
                return (
                    bool(scalar_references)
                    and all(scalar_references)
                    and all(reference in runtime_handles for reference in scalar_references),
                    True,
                )
            if owner_key in _UNIT_REFERENCE_KEYS and value:
                unit_references = [reference.strip() for reference in value.split(";")]
                runtime_references = [
                    reference for reference in unit_references if reference.startswith("h_")
                ]
                if not runtime_references:
                    return True, False
                return (
                    all(unit_references)
                    and all(
                        _RUNTIME_HANDLE_RE.fullmatch(reference) is not None
                        and reference in runtime_handles
                        for reference in runtime_references
                    ),
                    True,
                )
            return True, False
        if isinstance(value, list):
            if owner_key in _COMPOUND_EDGE_REFERENCE_KEYS:
                saw_reference = False
                for item in value:
                    if not isinstance(item, str):
                        return False, False
                    valid, saw_item = visit(item, owner_key)
                    if not valid:
                        return False, False
                    saw_reference = saw_reference or saw_item
                return True, saw_reference
            if owner_key in _ARRAY_REFERENCE_KEYS:
                if not all(
                    isinstance(item, str) and item in runtime_handles for item in value
                ):
                    return False, False
                return True, bool(value)
            saw_reference = False
            for item in value:
                valid, saw_item = visit(item)
                if not valid:
                    return False, False
                saw_reference = saw_reference or saw_item
            return True, saw_reference
        if isinstance(value, Mapping):
            saw_reference = False
            for raw_key, item in value.items():
                key = str(raw_key)
                if forbidden_reference(key):
                    return False, False
                if _RUNTIME_HANDLE_RE.fullmatch(key):
                    if key not in runtime_handles:
                        return False, False
                    saw_reference = True
                valid, saw_item = visit(item, key)
                if not valid:
                    return False, False
                saw_reference = saw_reference or saw_item
            return True, saw_reference
        return True, False

    for index, entry in enumerate(entries):
        if not is_current_derived_protein_artifact(entry):
            raise ValidationError(
                f"Unsupported current derived protein artifact at index {index}."
            )
        payload = entry["payload"]
        assert isinstance(payload, Mapping)
        valid, saw_reference = visit(payload)
        if not valid or (
            payload and not saw_reference and not _is_strict_empty_result(entry)
        ):
            raise ValidationError(
                "Current derived protein artifact contains unresolved protein "
                f"references at index {index}."
            )


# The derived LOSATP entry (OV-1063). One owner computes what a converted
# protein payload depends on (its identity), the key that names it, its
# provenance and its pair payloads. The Web converter, the Web lookup (through
# the `buildProteinLosatCacheKeys` reply) and the Python API call it, so equal
# inputs give one key on every surface.

DERIVED_LOSATP_SEMANTICS = "derived-option-conformance-v2"
DERIVED_LOSATP_CONVERTER = "convert_losatp_blastp_pairs_to_genomic_payload"

DerivedLosatpMode = Literal["pairwise", "orthogroup", "collinear"]
# A setting as a surface passes it: the Web sends numbers or text, and a blank
# or JS null/undefined value means the default.
DerivedSettingValue = Union[str, int, float, None]


@dataclass(frozen=True)
class DerivedLosatpSettings:
    """The converter settings of one derived payload, normalized once.

    Settings of another mode keep their defaults, so they cannot change the
    identity of this mode.
    """

    mode: DerivedLosatpMode
    max_hits: int
    bitscore: float
    evalue: float
    identity: float
    alignment_length: int
    membership_mode: str
    member_max_hits: int | None
    min_anchors: int
    max_unit_gap: int
    max_diagonal_drift: int
    max_conflicts: int
    merge_orientation: str
    unit_mode: str
    color_mode: str
    anchor_mode: str
    search_scope: str
    max_paralog_links: int
    infer_orthogroups: bool
    explicit_display_pairs: bool

    def collinearity_params(self) -> LosslessCollinearityParameters:
        from .collinearity import LosslessCollinearityParameters

        return LosslessCollinearityParameters(
            min_anchors=self.min_anchors,
            max_unit_gap=self.max_unit_gap,
            max_diagonal_drift=self.max_diagonal_drift,
            max_conflicts=self.max_conflicts,
            merge_orientation=cast(
                Literal["strand", "order", "either"], self.merge_orientation
            ),
        )


@dataclass(frozen=True)
class DerivedLosatpRecord:
    """One record the converter reads: its protein bindings and display frame."""

    index: int
    runtime_binding_hash: str
    display_binding_hash: str
    length: int
    reverse: bool


@dataclass(frozen=True)
class DerivedLosatpPair:
    """One directional search the converter reads.

    ``pair_index`` numbers a displayed pair (its payload carries it); a hidden
    pair (reverse, self, non-adjacent) has none.
    """

    query_index: int
    subject_index: int
    raw_key: str
    display_pair: bool
    pair_index: int | None = None


def _setting_int(value: DerivedSettingValue, default: int) -> int:
    return default if value is None else int(cast(Union[str, int, float], value))


def _setting_text(value: DerivedSettingValue, default: str) -> str:
    return default if value is None else str(value)


def normalize_derived_losatp_settings(
    *,
    mode: DerivedSettingValue = "pairwise",
    max_hits: DerivedSettingValue = 5,
    bitscore: DerivedSettingValue = 50,
    evalue: DerivedSettingValue = "1e-2",
    identity: DerivedSettingValue = 0,
    alignment_length: DerivedSettingValue = 0,
    membership_mode: DerivedSettingValue = None,
    member_max_hits: DerivedSettingValue = None,
    min_anchors: DerivedSettingValue = None,
    max_unit_gap: DerivedSettingValue = None,
    max_diagonal_drift: DerivedSettingValue = None,
    max_conflicts: DerivedSettingValue = None,
    merge_orientation: DerivedSettingValue = None,
    unit_mode: DerivedSettingValue = None,
    color_mode: DerivedSettingValue = None,
    anchor_mode: DerivedSettingValue = None,
    search_scope: DerivedSettingValue = None,
    max_paralog_links: DerivedSettingValue = None,
    infer_orthogroups: object = True,
    explicit_display_pairs: object = False,
) -> DerivedLosatpSettings:
    """Validate and normalize the settings of one LOSATP conversion.

    Numbers normalize to numbers, so ``0``, ``0.0`` and ``"0"`` give one
    identity; spellings normalize through their mode owners. ``None`` is an
    unset setting, which takes its default.
    """

    from .collinearity import (
        LosslessCollinearityParameters,
        normalize_collinearity_anchor_mode,
        normalize_collinearity_color_mode,
        normalize_collinearity_search_scope,
    )
    from .collinearity_units import normalize_collinearity_unit_mode
    from .protein_colinearity import normalize_orthogroup_membership_mode

    normalized_mode = str(mode or "pairwise").strip().lower()
    if normalized_mode not in {"pairwise", "orthogroup", "collinear"}:
        raise ValueError(f"Unsupported LOSATP blastp mode: {mode!r}")
    normalized_max_hits = 5
    if normalized_mode == "pairwise":
        normalized_max_hits = _setting_int(max_hits, 5)
        if normalized_max_hits <= 0:
            raise ValueError("protein_blastp_max_hits must be > 0")
    normalized_membership_mode = "anchor_core_v1"
    normalized_member_max_hits: int | None = None
    if normalized_mode in {"orthogroup", "collinear"}:
        normalized_membership_mode = normalize_orthogroup_membership_mode(
            _setting_text(membership_mode, "anchor_core_v1")
        )
        if member_max_hits is not None:
            normalized_member_max_hits = _setting_int(member_max_hits, 0)
            if normalized_member_max_hits <= 0:
                raise ValueError("orthogroup_member_max_hits must be > 0 or None")
    paralog_links = 2
    if normalized_mode in {"orthogroup", "collinear"}:
        paralog_links = _setting_int(max_paralog_links, 2)
        if paralog_links <= 0:
            raise ValueError("collinear_max_paralog_links_per_orthogroup must be > 0")
    params = LosslessCollinearityParameters()
    unit = "auto"
    color = "orientation"
    anchor = "rbh"
    scope = "adjacent"
    infer = True
    explicit = False
    if normalized_mode == "collinear":
        if not isinstance(infer_orthogroups, bool):
            raise ValueError("collinear_infer_orthogroups must be a boolean")
        infer = infer_orthogroups
        explicit = explicit_display_pairs is True
        unit = normalize_collinearity_unit_mode(_setting_text(unit_mode, "auto"))
        color = normalize_collinearity_color_mode(_setting_text(color_mode, "orientation"))
        anchor = normalize_collinearity_anchor_mode(_setting_text(anchor_mode, "rbh"))
        scope = normalize_collinearity_search_scope(_setting_text(search_scope, "adjacent"))
        params = LosslessCollinearityParameters(
            min_anchors=_setting_int(min_anchors, 1),
            max_unit_gap=_setting_int(max_unit_gap, 0),
            max_diagonal_drift=_setting_int(max_diagonal_drift, 0),
            max_conflicts=_setting_int(max_conflicts, 1),
            merge_orientation=cast(
                Literal["strand", "order", "either"],
                _setting_text(merge_orientation, "either").strip().lower(),
            ),
        )
        params.validate()
    normalized_bitscore = float(cast(Union[str, int, float], bitscore))
    normalized_evalue = float(cast(Union[str, int, float], evalue))
    normalized_identity = float(cast(Union[str, int, float], identity))
    normalized_alignment_length = int(cast(Union[str, int, float], alignment_length))
    if not math.isfinite(normalized_bitscore) or normalized_bitscore < 0:
        raise ValueError("bitscore must be a finite value >= 0")
    if not math.isfinite(normalized_evalue) or normalized_evalue < 0:
        raise ValueError("evalue must be a finite value >= 0")
    if (
        not math.isfinite(normalized_identity)
        or normalized_identity < 0
        or normalized_identity > 100
    ):
        raise ValueError("identity must be a finite value between 0 and 100")
    if normalized_alignment_length < 0:
        raise ValueError("alignment_length must be >= 0")
    return DerivedLosatpSettings(
        mode=cast(DerivedLosatpMode, normalized_mode),
        max_hits=normalized_max_hits,
        bitscore=normalized_bitscore,
        evalue=normalized_evalue,
        identity=normalized_identity,
        alignment_length=normalized_alignment_length,
        membership_mode=str(normalized_membership_mode),
        member_max_hits=normalized_member_max_hits,
        min_anchors=int(params.min_anchors),
        max_unit_gap=int(params.max_unit_gap),
        max_diagonal_drift=int(params.max_diagonal_drift),
        max_conflicts=int(params.max_conflicts),
        merge_orientation=str(params.merge_orientation),
        unit_mode=str(unit),
        color_mode=str(color),
        anchor_mode=str(anchor),
        search_scope=str(scope),
        max_paralog_links=paralog_links,
        infer_orthogroups=infer,
        explicit_display_pairs=explicit,
    )


def _derived_mode_settings(settings: DerivedLosatpSettings) -> dict[str, object]:
    """The settings blocks of the active mode, shared by identity and provenance."""

    blocks: dict[str, object] = {
        "thresholds": {
            "bitscore": settings.bitscore,
            "evalue": settings.evalue,
            "identity": settings.identity,
            "alignmentLength": settings.alignment_length,
        }
    }
    if settings.mode == "pairwise":
        blocks["pairwise"] = {"maxHits": settings.max_hits}
    else:
        blocks["orthogroup"] = {
            "membershipMode": settings.membership_mode,
            "memberMaxHits": settings.member_max_hits,
            "maxParalogLinksPerOrthogroup": settings.max_paralog_links,
        }
    return blocks


def _derived_collinear_settings(settings: DerivedLosatpSettings) -> dict[str, object]:
    return {
        "anchorMode": settings.anchor_mode,
        "searchScope": settings.search_scope,
        "inferOrthogroups": settings.infer_orthogroups,
        "colorMode": settings.color_mode,
        "parameters": {
            "minAnchors": settings.min_anchors,
            "maxUnitGap": settings.max_unit_gap,
            "maxDiagonalDrift": settings.max_diagonal_drift,
            "maxConflicts": settings.max_conflicts,
            "mergeOrientation": settings.merge_orientation,
        },
    }


def derived_losatp_identity(
    settings: DerivedLosatpSettings,
    records: Sequence[DerivedLosatpRecord],
    pairs: Sequence[DerivedLosatpPair],
) -> dict[str, object]:
    """Every input a converted LOSATP payload depends on, and nothing else.

    ``pairs`` keep the search order the converter reads. Every displayed pair
    has its own number (OV-1069); a hidden pair enters through its raw key only.
    """

    record_indexes = [int(record.index) for record in records]
    if len(set(record_indexes)) != len(record_indexes):
        raise ValidationError("Derived LOSATP identity has a duplicate record index.")
    known_records = set(record_indexes)
    display_numbers: set[int] = set()
    pair_identities: list[dict[str, object]] = []
    for position, pair in enumerate(pairs):
        if not isinstance(pair.raw_key, str) or not pair.raw_key:
            raise ValidationError(
                f"Derived LOSATP pair #{position + 1} has no raw cache key."
            )
        if (
            int(pair.query_index) not in known_records
            or int(pair.subject_index) not in known_records
        ):
            raise ValidationError(
                f"Derived LOSATP pair #{position + 1} references a missing record."
            )
        identity: dict[str, object] = {
            "queryIndex": int(pair.query_index),
            "subjectIndex": int(pair.subject_index),
            "rawKey": pair.raw_key,
            "displayPair": pair.display_pair is True,
        }
        if pair.display_pair is True:
            number = pair.pair_index
            if (
                not isinstance(number, int)
                or isinstance(number, bool)
                or number < 0
            ):
                raise ValidationError(
                    f"Derived LOSATP displayed pair #{position + 1} has no pair number."
                )
            if number in display_numbers:
                raise ValidationError(
                    f"Derived LOSATP displayed pairs share the pair number {number}."
                )
            display_numbers.add(number)
            identity["pairIndex"] = number
        pair_identities.append(identity)
    result: dict[str, object] = {
        "cacheSchema": CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA,
        "semantics": DERIVED_LOSATP_SEMANTICS,
        "idEncoding": "runtime-handle-v1",
        "converter": DERIVED_LOSATP_CONVERTER,
        "mode": settings.mode,
        "rawCacheKeys": sorted({pair.raw_key for pair in pairs}),
        **_derived_mode_settings(settings),
        "records": [
            {
                "recordIndex": int(record.index),
                "runtimeBindingHash": str(record.runtime_binding_hash),
                "displayBindingHash": str(record.display_binding_hash),
                "viewTransform": {
                    "length": int(record.length),
                    "reverse": record.reverse is True,
                },
            }
            for record in sorted(records, key=lambda current: int(current.index))
        ],
        "pairs": pair_identities,
    }
    if settings.mode != "pairwise":
        result["pathRepresentation"] = "lossless-graph-v1"
    if settings.mode == "collinear":
        result["collinear"] = {
            **_derived_collinear_settings(settings),
            "unitMode": settings.unit_mode,
            # With the search scope 'all' a row grid limits the output to its displayed pairs.
            "explicitDisplayPairs": settings.explicit_display_pairs,
        }
    return result


def derived_losatp_cache_key(identity: Mapping[str, object]) -> str:
    """The SHA-256 hex key of one derived identity (its canonical JSON)."""

    return hashlib.sha256(
        json.dumps(
            identity,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def effective_collinearity_unit_kinds(result: CollinearityResult | None) -> list[str]:
    """The unit kinds the anchors of a collinearity result use."""

    if result is None:
        return []
    anchors = [anchor for block in result.blocks for anchor in block.anchors]
    anchors.extend(result.unblocked_anchors)
    return sorted(
        {
            str(unit_kind)
            for anchor in anchors
            for unit_kind in (anchor.query_unit_kind, anchor.subject_unit_kind)
        }
    )


def derived_losatp_provenance(
    identity: Mapping[str, object],
    settings: DerivedLosatpSettings,
    *,
    effective_unit_kinds: Sequence[str] = (),
) -> dict[str, object]:
    """What a derived payload was computed from, for a reader of the entry."""

    provenance: dict[str, object] = {
        "schema": 1,
        "upstreamRawKeys": list(cast(Sequence[str], identity["rawCacheKeys"])),
        **_derived_mode_settings(settings),
    }
    if settings.mode == "collinear":
        provenance["collinear"] = {
            "unitMode": {
                "requested": settings.unit_mode,
                "effectiveKinds": sorted(str(kind) for kind in effective_unit_kinds),
            },
            **_derived_collinear_settings(settings),
        }
    return provenance


def derived_pair_payload(
    pair_index: int,
    query_index: int,
    subject_index: int,
    frame: DataFrame,
) -> dict[str, object]:
    """One displayed pair of a derived payload: the headerless BLAST-style rows
    for Run info and the canonical comparison table a request uses."""

    from gbdraw.io.comparisons import COMPARISON_COLUMNS
    from gbdraw.session_request_codec import encode_canonical_table

    tsv = frame.loc[:, list(COMPARISON_COLUMNS)].to_csv(
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n",
    )
    return {
        "pair_index": int(pair_index),
        "query_index": int(query_index),
        "subject_index": int(subject_index),
        "tsv": tsv,
        "canonicalTsv": encode_canonical_table(frame, name=f"comparison pair {pair_index}"),
        "hit_count": int(frame.shape[0]),
    }


__all__ = [
    "CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA",
    "DERIVED_LOSATP_SEMANTICS",
    "DerivedLosatpPair",
    "DerivedLosatpRecord",
    "DerivedLosatpSettings",
    "derived_losatp_cache_key",
    "derived_losatp_identity",
    "derived_losatp_provenance",
    "derived_pair_payload",
    "effective_collinearity_unit_kinds",
    "is_current_derived_protein_artifact",
    "normalize_derived_losatp_settings",
    "validate_current_derived_protein_artifacts",
]
