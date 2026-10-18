"""The identity of a derived LOSATP entry has one owner (OV-1063).

G-F(1) (Web GUI audit 2026-09-30; OIPC-C04, R5): a converted protein payload
is reused only when every input it depends on is unchanged (CO-02, CO-03).
Each row changes one input and must change the key; a display-only or
inactive input must keep it. A builder that drops an input from its identity
fails the table (the last test). The Web lookup, the Web converter and the
Python API all hash this identity, so the table covers every surface.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Callable

import pytest

from gbdraw.analysis.protein_artifacts import (
    DERIVED_LOSATP_SEMANTICS,
    DerivedLosatpPair,
    DerivedLosatpRecord,
    derived_losatp_cache_key,
    derived_losatp_identity,
    derived_losatp_provenance,
    effective_collinearity_unit_kinds,
    normalize_derived_losatp_settings,
)
from gbdraw.exceptions import ValidationError

Inputs = dict[str, object]
Build = Callable[[Inputs], dict[str, object]]

_BASELINE: Inputs = {
    "settings": {
        "mode": "orthogroup",
        "max_hits": 5,
        "bitscore": 50,
        "evalue": "1e-5",
        "identity": 70,
        "alignment_length": 0,
        "membership_mode": "anchor_core_v1",
        "member_max_hits": 5,
    },
    "records": (
        DerivedLosatpRecord(0, "runtime-a", "display-a", 100, False),
        DerivedLosatpRecord(1, "runtime-b", "display-b", 200, True),
    ),
    "pairs": (
        DerivedLosatpPair(0, 1, "raw-a-b", True, 0),
        DerivedLosatpPair(1, 0, "raw-b-a", False),
        DerivedLosatpPair(0, 0, "raw-a-a", False),
    ),
}
_COLLINEAR: Inputs = {
    **_BASELINE,
    "settings": {
        **_BASELINE["settings"],  # type: ignore[dict-item]
        "mode": "collinear",
        "min_anchors": 3,
        "max_unit_gap": 2,
        "unit_mode": "auto",
        "color_mode": "orientation",
        "anchor_mode": "rbh",
        "merge_orientation": "either",
        "max_diagonal_drift": 5,
        "max_conflicts": 1,
        "max_paralog_links": 2,
        "infer_orthogroups": True,
        "search_scope": "adjacent",
    },
}


def _identity(inputs: Inputs) -> dict[str, object]:
    settings = normalize_derived_losatp_settings(**inputs["settings"])  # type: ignore[arg-type]
    return derived_losatp_identity(settings, inputs["records"], inputs["pairs"])  # type: ignore[arg-type]


def _set(**changes: object) -> Callable[[Inputs], Inputs]:
    return lambda inputs: {**inputs, "settings": {**inputs["settings"], **changes}}  # type: ignore[dict-item]


def _record(position: int, **changes: object) -> Callable[[Inputs], Inputs]:
    def change(inputs: Inputs) -> Inputs:
        records = list(inputs["records"])  # type: ignore[call-overload]
        records[position] = replace(records[position], **changes)
        return {**inputs, "records": tuple(records)}

    return change


def _pair(position: int, **changes: object) -> Callable[[Inputs], Inputs]:
    def change(inputs: Inputs) -> Inputs:
        pairs = list(inputs["pairs"])  # type: ignore[call-overload]
        pairs[position] = replace(pairs[position], **changes)
        return {**inputs, "pairs": tuple(pairs)}

    return change


_SCIENTIFIC = [
    ("the record set", lambda inputs: {**inputs, "records": (*inputs["records"], DerivedLosatpRecord(2, "r", "d", 9, False))}),
    ("a record region length", _record(0, length=90)),
    ("a record orientation", _record(0, reverse=True)),
    ("a record runtime binding", _record(1, runtime_binding_hash="runtime-b2")),
    ("a record display binding (protein set, translation table, Feature visibility)", _record(1, display_binding_hash="display-b2")),
    ("a pair endpoint", _pair(0, query_index=1, subject_index=0)),
    ("a displayed pair's raw search", _pair(0, raw_key="raw-a-b-rerun")),
    ("a hidden pair's raw search (reverse)", _pair(1, raw_key="raw-b-a-rerun")),
    ("a hidden pair's raw search (self)", _pair(2, raw_key="raw-a-a-rerun")),
    ("the displayed direction of a pair (displayPair)", _pair(0, display_pair=False)),
    ("the number of a displayed pair (pairIndex)", _pair(0, pair_index=1)),
    ("the search order", lambda inputs: {**inputs, "pairs": tuple(reversed(inputs["pairs"]))}),  # type: ignore[call-overload]
    ("the bit score threshold", _set(bitscore=60)),
    ("the E-value threshold", _set(evalue="1e-6")),
    ("the identity threshold", _set(identity=71)),
    ("the alignment length threshold", _set(alignment_length=30)),
    ("the orthogroup member max hits", _set(member_max_hits=6)),
    ("the maximum paralog links per orthogroup", _set(max_paralog_links=1)),
]
_DISPLAY_ONLY = [
    ("record input order", lambda inputs: {**inputs, "records": tuple(reversed(inputs["records"]))}),  # type: ignore[call-overload]
    ("the number of a hidden pair", _pair(1, pair_index=7)),
    ("a number written as text or float", _set(bitscore="50", identity=70.0, alignment_length="0")),
    ("an equal E-value spelled another way", _set(evalue=1e-5)),
    ("the Pairwise max hits", _set(max_hits=99)),
    ("Collinear-only settings", _set(unit_mode="cds", anchor_mode="all", search_scope="all", explicit_display_pairs=True)),
]
_ORTHOGROUP_TABLE = {
    "scientific": [*_SCIENTIFIC, ("the comparison mode", _set(mode="pairwise"))],
    "displayOnly": _DISPLAY_ONLY,
}
_COLLINEAR_TABLE = {
    "scientific": [
        *_SCIENTIFIC,
        ("the comparison mode", _set(mode="orthogroup")),
        ("the minimum anchors", _set(min_anchors=4)),
        ("the maximum unit gap", _set(max_unit_gap=3)),
        ("the unit mode", _set(unit_mode="cds")),
        ("the color mode", _set(color_mode="orientation_identity")),
        ("the anchor mode", _set(anchor_mode="all")),
        ("the merge orientation", _set(merge_orientation="strand")),
        ("the maximum diagonal drift", _set(max_diagonal_drift=6)),
        ("the maximum conflicts in a merge gap", _set(max_conflicts=2)),
        ("orthogroup inference", _set(infer_orthogroups=False)),
        ("the search scope", _set(search_scope="all")),
        # Only the Collinear converter reads it: with the search scope 'all' it
        # limits the output to the displayed pairs of a CLI row grid.
        ("explicit display pairs (CLI row grid)", _set(explicit_display_pairs=True)),
    ],
    "displayOnly": [row for row in _DISPLAY_ONLY if row[0] != "Collinear-only settings"]
    + [("Pairwise-only settings", _set(max_hits=99))],
}


def _problems(build: Build, baseline: Inputs, table: dict[str, list]) -> list[str]:
    baseline_key = derived_losatp_cache_key(build(baseline))
    problems = []
    for name, change in table["scientific"]:
        if derived_losatp_cache_key(build(change(baseline))) == baseline_key:
            problems.append(f"reused after a change of {name}")
    for name, change in table["displayOnly"]:
        if derived_losatp_cache_key(build(change(baseline))) != baseline_key:
            problems.append(f"not reused after display-only {name}")
    return problems


def test_derived_identity_changes_with_every_input_it_depends_on() -> None:
    assert _problems(_identity, _BASELINE, _ORTHOGROUP_TABLE) == []
    assert _problems(_identity, _COLLINEAR, _COLLINEAR_TABLE) == []
    pairwise = _set(mode="pairwise")(_BASELINE)
    assert derived_losatp_cache_key(_identity(_set(max_hits=6)(pairwise))) != derived_losatp_cache_key(
        _identity(pairwise)
    )


def test_the_sensitivity_table_flags_a_builder_that_drops_an_input() -> None:
    def dropping(*paths: tuple[str, ...]) -> Build:
        def build(inputs: Inputs) -> dict[str, object]:
            identity = _identity(inputs)
            for path in paths:
                *parents, name = path
                for item in _walk(identity, parents):
                    item.pop(name, None)
            return identity

        return build

    def _walk(value: object, parents: list[str]) -> list[dict[str, object]]:
        if not parents:
            return [value] if isinstance(value, dict) else [item for item in value]  # type: ignore[union-attr]
        child = value.get(parents[0])  # type: ignore[union-attr]
        return [] if child is None else _walk(child, parents[1:])

    assert _problems(dropping(("records", "viewTransform")), _BASELINE, _ORTHOGROUP_TABLE) == [
        "reused after a change of a record region length",
        "reused after a change of a record orientation",
    ]
    # A displayed pair also carries its number, so only dropping both loses the direction.
    assert _problems(dropping(("pairs", "displayPair"), ("pairs", "pairIndex")), _BASELINE, _ORTHOGROUP_TABLE) == [
        "reused after a change of the displayed direction of a pair (displayPair)",
        "reused after a change of the number of a displayed pair (pairIndex)",
    ]
    assert _problems(dropping(("pairs", "pairIndex")), _BASELINE, _ORTHOGROUP_TABLE) == [
        "reused after a change of the number of a displayed pair (pairIndex)",
    ]
    assert _problems(dropping(("collinear", "explicitDisplayPairs")), _COLLINEAR, _COLLINEAR_TABLE) == [
        "reused after a change of explicit display pairs (CLI row grid)"
    ]


def test_derived_identity_shape_and_semantics() -> None:
    identity = _identity(_COLLINEAR)
    assert identity["semantics"] == DERIVED_LOSATP_SEMANTICS == "derived-option-conformance-v2"
    assert identity["cacheSchema"] == 3
    assert identity["pathRepresentation"] == "lossless-graph-v1"
    assert identity["rawCacheKeys"] == ["raw-a-a", "raw-a-b", "raw-b-a"]
    assert identity["pairs"] == [
        {"queryIndex": 0, "subjectIndex": 1, "rawKey": "raw-a-b", "displayPair": True, "pairIndex": 0},
        {"queryIndex": 1, "subjectIndex": 0, "rawKey": "raw-b-a", "displayPair": False},
        {"queryIndex": 0, "subjectIndex": 0, "rawKey": "raw-a-a", "displayPair": False},
    ]
    assert [record["recordIndex"] for record in identity["records"]] == [0, 1]  # type: ignore[index]
    assert identity["thresholds"] == {"bitscore": 50.0, "evalue": 1e-5, "identity": 70.0, "alignmentLength": 0}
    # Record instance keys, protein cache keys and feature identity are implied
    # by the binding hashes, so no surface-specific name enters the identity.
    for name in ("featureIdentity", "proteinCacheKey", "recordInstanceKey"):
        assert name not in str(identity)
    assert "collinear" not in _identity(_BASELINE) and "pairwise" not in _identity(_BASELINE)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (_pair(0, raw_key=""), "no raw cache key"),
        (_pair(0, pair_index=None), "has no pair number"),
        (_pair(0, query_index=5), "missing record"),
        # OV-1069: two displayed pairs never share one number.
        (lambda inputs: {**inputs, "pairs": (*inputs["pairs"], DerivedLosatpPair(1, 1, "raw-b-b", True, 0))}, "share the pair number 0"),  # type: ignore[call-overload]
        (_record(1, index=0), "duplicate record index"),
    ],
)
def test_derived_identity_refuses_inputs_no_converter_reads(
    change: Callable[[Inputs], Inputs], message: str
) -> None:
    with pytest.raises(ValidationError, match=message):
        _identity(change(_BASELINE))


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        ({"mode": "other"}, "Unsupported LOSATP blastp mode"),
        ({"mode": "pairwise", "max_hits": 0}, "max_hits must be > 0"),
        ({"mode": "orthogroup", "member_max_hits": 0}, "member_max_hits must be > 0"),
        ({"mode": "collinear", "infer_orthogroups": "yes"}, "must be a boolean"),
        ({"mode": "collinear", "min_anchors": 0}, "collinear_min_anchors must be > 0"),
        ({"mode": "orthogroup", "identity": 101}, "identity must be a finite value"),
        ({"mode": "orthogroup", "evalue": "nan"}, "evalue must be a finite value"),
    ],
)
def test_derived_settings_are_validated_once(settings: dict[str, object], message: str) -> None:
    with pytest.raises((ValueError, ValidationError), match=message):
        normalize_derived_losatp_settings(**settings)  # type: ignore[arg-type]


def test_unset_settings_take_the_defaults() -> None:
    # The Web adapter turns JS null, undefined, and blank text into None
    # (tests/test_api_request_render.py); the owner reads only None as unset.
    defaults = normalize_derived_losatp_settings(mode="collinear")
    assert normalize_derived_losatp_settings(
        mode="collinear", min_anchors=None, max_unit_gap=None, max_conflicts=None,
        unit_mode=None, member_max_hits=None, max_paralog_links=None,
    ) == defaults
    with pytest.raises(ValueError):
        normalize_derived_losatp_settings(mode="collinear", min_anchors="")


def test_derived_provenance_records_requested_and_effective_units() -> None:
    from types import SimpleNamespace

    settings = normalize_derived_losatp_settings(**_COLLINEAR["settings"])  # type: ignore[arg-type]
    identity = _identity(_COLLINEAR)
    anchor = SimpleNamespace(query_unit_kind="locus", subject_unit_kind="cds")
    result = SimpleNamespace(blocks=(SimpleNamespace(anchors=(anchor,)),), unblocked_anchors=())
    provenance = derived_losatp_provenance(
        identity, settings, effective_unit_kinds=effective_collinearity_unit_kinds(result)  # type: ignore[arg-type]
    )
    assert provenance["upstreamRawKeys"] == identity["rawCacheKeys"]
    assert provenance["thresholds"] == identity["thresholds"]
    assert provenance["collinear"]["unitMode"] == {"requested": "auto", "effectiveKinds": ["cds", "locus"]}  # type: ignore[index]
    assert provenance["orthogroup"]["maxParalogLinksPerOrthogroup"] == 2  # type: ignore[index]
    assert effective_collinearity_unit_kinds(None) == []
