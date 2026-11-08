from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.utils.session_drawing_view import drawing_view
from gbdraw.exceptions import ValidationError
from gbdraw.web_support.similarity_alignment import (
    resolve_similarity_alignment_json,
    resolve_similarity_alignment_payload,
)


def _anchor(record_key: str, feature_id: str, index: int) -> dict[str, object]:
    return {
        "recordKey": record_key,
        "biologicalFeatureId": feature_id,
        "sourceFeatureIndex": index,
        "stableFeatureSvgId": f"stable-{feature_id}",
    }


def _member(
    record_key: str,
    feature_id: str,
    index: int,
    *,
    start: int,
    end: int,
    strand: int = 1,
) -> dict[str, object]:
    return {
        "groupId": "og-1",
        "anchor": _anchor(record_key, feature_id, index),
        "sourceStart": start,
        "sourceEnd": end,
        "sourceStrand": strand,
        "identityIsUnique": True,
        "hidden": False,
        "representative": False,
        "role": "inparalog",
    }


def _request() -> dict[str, object]:
    reference = _anchor("record-a", "reference", 0)
    return {
        "schema": 2,
        "groupId": "og-1",
        "records": [
            {
                "recordKey": "record-a",
                "recordLength": 1000,
                "windows": [],
                "presentation": {"reverseComplement": False},
            },
            {
                "recordKey": "record-b",
                "recordLength": 400,
                "windows": [{"start": 101, "end": 300, "reverseComplement": True, "column": None}],
                "presentation": {"reverseComplement": False},
            },
        ],
        "reference": reference,
        "members": [
            _member("record-a", "reference", 0, start=10, end=40),
            _member("record-b", "target", 1, start=140, end=170, strand=-1),
        ],
        "directEdges": [
            {
                "groupId": "og-1",
                "query": reference,
                "subject": _anchor("record-b", "target", 1),
                "edgeKind": "rbh",
            }
        ],
        "choices": [],
    }


def test_web_adapter_projects_crop_facts_and_returns_canonical_plan() -> None:
    result = resolve_similarity_alignment_payload(_request())

    assert result["status"] == "resolved"
    assert result["plan"] == {
        "schema": 2,
        "groupId": "og-1",
        "reference": _anchor("record-a", "reference", 0),
        "records": [
            {
                "recordKey": "record-a",
                "status": "reference",
                "rationale": "reference",
                "anchor": _anchor("record-a", "reference", 0),
            },
            {
                "recordKey": "record-b",
                "status": "aligned",
                "rationale": "only_usable_candidate",
                "anchor": _anchor("record-b", "target", 1),
            },
        ],
    }
    assert json.loads(resolve_similarity_alignment_json(json.dumps(_request()))) == result


def test_web_adapter_returns_ambiguity_then_uses_explicit_skip() -> None:
    request = _request()
    request["directEdges"] = []
    request["members"].append(  # type: ignore[union-attr]
        _member("record-b", "other", 2, start=190, end=220)
    )

    ambiguous = resolve_similarity_alignment_payload(request)
    assert ambiguous["status"] == "ambiguous"
    assert ambiguous["plan"] is None
    assert ambiguous["records"][1]["kind"] == "ambiguous"  # type: ignore[index]

    request["choices"] = [
        {"recordKey": "record-b", "kind": "skip", "anchor": None}
    ]
    resolved = resolve_similarity_alignment_payload(request)
    assert resolved["status"] == "resolved"
    assert resolved["plan"]["records"][1]["rationale"] == "skipped_by_user"  # type: ignore[index]


def test_web_adapter_rejects_unknown_fields_and_unmappable_reference() -> None:
    unknown = _request()
    unknown["unexpected"] = True
    with pytest.raises(ValidationError, match="invalid fields"):
        resolve_similarity_alignment_payload(unknown)

    cropped = _request()
    cropped["records"][0]["windows"] = [{  # type: ignore[index]
        "start": 500,
        "end": 600,
        "reverseComplement": False,
        "column": None,
    }]
    with pytest.raises(ValidationError, match="exact reference.*current crop"):
        resolve_similarity_alignment_payload(cropped)


def test_web_adapter_reads_a_window_across_the_origin_and_places_a_record_of_several_windows() -> None:
    request = _request()
    # record-b (400 bp) drawn from 301 across the origin to 160: the target 140..170 centres at 155,
    # which the window holds at 100 + 155 = 255, reversed at 160 - 155 = 5 (OV-2771).
    window = {"start": 301, "end": 160, "reverseComplement": False, "column": None}
    request["records"][1]["windows"] = [window]  # type: ignore[index]
    assert resolve_similarity_alignment_payload(request)["records"][1]["candidates"][0]["displayCenter"] == 255.0  # type: ignore[index]
    window["reverseComplement"] = True
    assert resolve_similarity_alignment_payload(request)["records"][1]["candidates"][0]["displayCenter"] == 5.0  # type: ignore[index]
    window["end"] = 150
    outside = resolve_similarity_alignment_payload(request)["records"][1]  # type: ignore[index]
    assert outside["candidates"][0]["displayCenter"] is None  # type: ignore[index]
    assert outside["status"] == "skipped"  # type: ignore[index]

    several = _request()
    several["records"][1]["windows"] = [  # type: ignore[index]
        {"start": 1, "end": 100, "reverseComplement": False, "column": None},
        {"start": 131, "end": 200, "reverseComplement": False, "column": None},
    ]
    # The record is anchored in the window that holds the target (131..200 draws 155 at 25), and aligns (D-50).
    row = resolve_similarity_alignment_payload(several)["records"][1]  # type: ignore[index]
    assert (row["candidates"][0]["displayCenter"], row["candidates"][0]["displayedStrand"]) == (25.0, -1)  # type: ignore[index]
    assert row["status"] == "aligned"  # type: ignore[index]


def test_web_projection_discloses_candidate_evidence_and_strand_relation() -> None:
    request = _request()
    request["members"][1]["displayName"] = "target gene"  # type: ignore[index]
    result = resolve_similarity_alignment_payload(request)
    row = result["records"][1]
    assert row["reviewReason"] == "only_usable_candidate"
    assert row["anchor"] == _anchor("record-b", "target", 1)
    candidate = row["candidates"][0]
    assert candidate["displayName"] == "target gene"
    assert (candidate["sourceStart"], candidate["sourceEnd"]) == (140, 170)
    assert candidate["displayCenter"] == 145.0
    assert candidate["displayedStrand"] == 1
    assert candidate["representative"] is False
    assert candidate["directEvidence"] == ["rbh"]
    assert candidate["strandRelation"] == "same"
    assert result["referenceDisplayedStrand"] == 1
    assert result["referenceDisplayCenter"] == 25.0


def test_web_recommendation_is_transient_and_explicit_select_is_validated() -> None:
    request = _request()
    request["directEdges"] = []
    request["members"][1]["sourceStrand"] = 1  # type: ignore[index]
    replacement = _member(
        "record-b", "other", 2, start=190, end=220, strand=1
    )
    replacement["representative"] = True
    request["members"].append(replacement)  # type: ignore[union-attr]
    unresolved = resolve_similarity_alignment_payload(request)
    row = unresolved["records"][1]
    assert row["recommendationReason"] == "unique_representative"
    assert row["recommendedAnchor"] == _anchor("record-b", "other", 2)
    assert row["candidates"][0]["strandRelation"] == "opposite"
    assert unresolved["plan"] is None

    request["choices"] = [{
        "recordKey": "record-b", "kind": "select",
        "anchor": _anchor("record-b", "target", 1),
    }]
    selected = resolve_similarity_alignment_payload(request)
    assert selected["plan"]["schema"] == 2
    assert selected["plan"]["records"][1]["anchor"] == _anchor("record-b", "target", 1)
    assert set(selected["plan"]["records"][1]) == {
        "recordKey", "status", "rationale", "anchor"
    }


def test_web_rejects_removed_choice_field_and_invalid_crop() -> None:
    request = _request()
    request["choices"] = [{
        "recordKey": "record-b", "kind": "skip", "anchor": None,
    }]
    request["choices"][0]["orientationPolicy"] = "match_reference"  # type: ignore[index]
    with pytest.raises(ValidationError, match="invalid fields"):
        resolve_similarity_alignment_payload(request)

    # A window whose End runs past the record end is drawn clamped to the record end (the crop's
    # `region_slice`), so record-b's window 101..300 of a 200-bp record draws 101..200, reversed: the
    # target centre 155 is at 200 - 155 = 45, not at 300 - 155.
    cropped = _request()
    cropped["records"][1]["recordLength"] = 200  # type: ignore[index]
    row = resolve_similarity_alignment_payload(cropped)["records"][1]  # type: ignore[index]
    assert row["candidates"][0]["displayCenter"] == 45.0  # type: ignore[index]
    assert row["status"] == "aligned"  # type: ignore[index]
    # A window that starts past the record end is refused, as the crop refuses it.
    cropped["records"][1]["recordLength"] = 100  # type: ignore[index]
    with pytest.raises(ValidationError, match="starts after the end of record") as refused:
        resolve_similarity_alignment_payload(cropped)
    assert refused.value.diagnostic == {"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"}


def test_web_crop_after_presentation_reverse_uses_source_coordinates() -> None:
    request = _request()
    target_record = request["records"][1]
    target_record["presentation"]["reverseComplement"] = True
    target_record["windows"][0]["reverseComplement"] = False
    request["members"][1]["sourceStrand"] = 1
    result = resolve_similarity_alignment_payload(request)
    candidate = result["records"][1]["candidates"][0]
    # Source center 155 maps to 245 after the 400-bp presentation reversal,
    # then to 145 in the 101..300 crop.
    assert candidate["displayCenter"] == 145.0
    assert candidate["displayedStrand"] == -1
    assert candidate["strandRelation"] == "opposite"

    target_record["windows"][0]["reverseComplement"] = True
    with pytest.raises(ValidationError, match="both window and presentation"):
        resolve_similarity_alignment_payload(request)


BGC_SESSION = (
    Path(__file__).resolve().parents[1]
    / "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json"
)


@pytest.fixture(scope="module")
def bgc_og18(tmp_path_factory):
    """og_18 helper inputs built as the Web UI builds them: catalog members, no edges."""

    import base64

    from gbdraw.analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupResult
    from gbdraw.session_io import materialize_embedded_file
    from gbdraw.session_request_codec import decode_canonical_typed_resource
    from gbdraw.web_support.request_render import render_embedded_canonical_web_request

    session = drawing_view(json.loads(BGC_SESSION.read_text()))
    resource = decode_canonical_typed_resource(
        base64.b64decode(session["resources"]["comparison-orthogroups"]["data"]),
        value_kind="orthogroupResult", expected=OrthogroupResult | OrthogroupGraphResult,
    )
    gene = {m.protein_id: m.gene or m.label for m in resource.orthogroups["og_18"]}
    tmp = tmp_path_factory.mktemp("bgc-og18")
    rendered = render_embedded_canonical_web_request(
        session["renderRequest"], resources=session["resources"], workspace=str(tmp / "render")
    )
    item = rendered["featureCatalog"]["items"][0]
    group = next(entry for entry in item["orthogroups"] if entry["id"] == "og_18")
    features = {(f["recordKey"], f["biologicalFeatureId"]): f for f in item["biologicalFeatures"]}
    lengths = dict(zip(item["recordKeys"], (len(s["sequence"]) for s in item["sequenceSources"])))
    members, anchors = [], {}
    for member in group["members"]:
        feature = features[(member["recordKey"], member["biologicalFeatureId"])]
        anchor = {
            "recordKey": member["recordKey"],
            "biologicalFeatureId": member["biologicalFeatureId"],
            "sourceFeatureIndex": feature.get("sourceFeatureIndex"),
            "stableFeatureSvgId": feature.get("stableFeatureId") or member["biologicalFeatureId"],
        }
        members.append({
            "groupId": "og_18", "anchor": anchor, "sourceStart": feature["start"],
            "sourceEnd": feature["end"], "sourceStrand": {"+": 1, "-": -1}.get(feature["strand"], feature["strand"]),
            "identityIsUnique": True, "hidden": bool(member.get("hidden")),
            "representative": bool(member.get("representative")), "role": str(member.get("role") or ""),
        })
        anchors[gene[feature["protein_id"]]] = anchor
    paths = {
        rid: str(materialize_embedded_file(entry, temp_dir=tmp / "resources", role=rid, prefix_role=False))
        for rid, entry in session["resources"].items()
    }
    records = [{
        "recordKey": record["recordKey"], "recordLength": lengths[record["recordKey"]],
        "windows": record["windows"],
        "presentation": {"reverseComplement": bool(record["presentation"]["reverseComplement"])},
    } for record in session["renderRequest"]["records"]]
    return session, members, anchors, paths, records, tmp


def _resolve_bgc(bgc_og18, reference, *, choices=(), canonical=None, records=None, edges=()):
    session, members, anchors, paths, default_records, tmp = bgc_og18
    request = {"schema": 2, "groupId": "og_18", "records": records or default_records,
               "reference": anchors[reference], "members": members,
               "directEdges": list(edges), "choices": list(choices)}
    projection = {"canonicalRequest": canonical or session["renderRequest"], "orientations": None}
    return json.loads(resolve_similarity_alignment_json(
        json.dumps(request), json.dumps(projection), json.dumps(paths), str(tmp / "helper")))


def _record5(response):
    return next(r for r in response["records"] if r["recordKey"] == "record-5")


def test_committed_orthogroup_resource_supplies_direct_rbh_evidence(bgc_og18) -> None:
    anchors = bgc_og18[2]
    liv_a = _resolve_bgc(bgc_og18, "livA")
    assert liv_a["status"] == "resolved"
    row = _record5(liv_a)
    assert (row["rationale"], row["anchor"]) == ("unique_direct_rbh", anchors["racM"])
    evidence = {tuple(sorted(c["anchor"].items())): c["directEvidence"] for c in row["candidates"]}
    assert evidence[tuple(sorted(anchors["racM"].items()))] == ["rbh"]
    assert evidence[tuple(sorted(anchors["racL"].items()))] == ["coortholog"]

    # A representative-only selector would pick racM here too; parA's RBH is racL.
    par_a = _record5(_resolve_bgc(bgc_og18, "parA"))
    assert (par_a["rationale"], par_a["anchor"]) == ("unique_direct_rbh", anchors["racL"])

    explicit = _record5(_resolve_bgc(bgc_og18, "livA", choices=[
        {"recordKey": "record-5", "kind": "select", "anchor": anchors["racL"]}]))
    assert (explicit["rationale"], explicit["anchor"]) == ("user_selected", anchors["racL"])


def test_resource_edges_are_the_only_authority_and_never_guessed(bgc_og18) -> None:
    session, _members, anchors, _paths, records, _tmp = bgc_og18
    with pytest.raises(ValidationError, match="directEdges must be empty"):
        _resolve_bgc(bgc_og18, "livA", edges=[{
            "groupId": "og_18", "query": anchors["livA"], "subject": anchors["racL"], "edgeKind": "rbh"}])

    # A reorder with index-bound pairwise tables is rejected at decode.
    canonical = json.loads(json.dumps(session["renderRequest"]))
    order = [4, 1, 2, 3, 0]
    canonical["records"] = [canonical["records"][index] for index in order]
    reordered = [records[index] for index in order]
    with pytest.raises(ValidationError):
        _resolve_bgc(bgc_og18, "livA", canonical=canonical, records=reordered)
    # With only the orthogroup resource left, its member record indexes are
    # stale; the endpoints no longer bind, so record-5 needs Review, not a guess.
    canonical["comparisons"] = [
        item for item in canonical["comparisons"] if item["kind"] == "orthogroupResult"]
    stale = _resolve_bgc(bgc_og18, "livA", canonical=canonical, records=reordered)
    assert stale["status"] == "ambiguous"
    assert _record5(stale)["kind"] == "ambiguous"
    assert _record5(stale)["directRbhCandidates"] == []
