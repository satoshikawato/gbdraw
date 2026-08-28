"""Session-neutral validation for current derived protein artifacts."""

from __future__ import annotations

import re
from typing import Mapping, Sequence

from gbdraw.exceptions import ValidationError

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


def _is_nonnegative_integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _is_strict_empty_result(entry: Mapping[str, object]) -> bool:
    mode = entry.get("mode")
    if mode not in {"pairwise", "orthogroup", "collinear"}:
        return False
    payload = entry.get("payload")
    if not isinstance(payload, Mapping):
        return False

    allowed_keys = {
        "identity",
        "provenance",
        "pairs",
        "orthogroups",
        "orthogroupResult",
        "collinearityResult",
    }
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

    for resource_key, resource_kind in (
        ("orthogroupResult", "orthogroupResult"),
        ("collinearityResult", "result"),
    ):
        if resource_key not in payload:
            continue
        resource = payload[resource_key]
        if (
            not isinstance(resource, Mapping)
            or resource.get("schema") != 2
            or resource.get("kind") != resource_kind
            or "value" not in resource
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
            or pair.get("rows") != []
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
    manifest: Mapping[str, object] | None,
) -> None:
    """Require current envelopes and resolve every protein reference."""

    if manifest is None:
        raise ValidationError(
            "Current derived protein artifacts require protein_identity_manifest."
        )
    from .protein_colinearity import validate_protein_identity_manifest

    try:
        authority = validate_protein_identity_manifest(manifest)
    except (TypeError, ValueError, ValidationError) as exc:
        raise ValidationError(
            "Current derived protein artifacts require a valid "
            "protein_identity_manifest."
        ) from exc

    runtime_handles = {
        str(handle)
        for binding in authority.record_instances.values()
        for handle in (
            binding.get("runtimeIds", {}).values()
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
                references = [reference.strip() for reference in value.split(";")]
                return (
                    bool(references)
                    and all(references)
                    and all(reference in runtime_handles for reference in references),
                    True,
                )
            if owner_key in _UNIT_REFERENCE_KEYS and value:
                references = [reference.strip() for reference in value.split(";")]
                runtime_references = [
                    reference for reference in references if reference.startswith("h_")
                ]
                if not runtime_references:
                    return True, False
                return (
                    all(references)
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


__all__ = [
    "CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA",
    "is_current_derived_protein_artifact",
    "validate_current_derived_protein_artifacts",
]
