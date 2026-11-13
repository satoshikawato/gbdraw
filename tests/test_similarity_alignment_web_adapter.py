from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.utils.session_drawing_view import drawing_view
from tests.utils.session_files import read_session_payload
from gbdraw.exceptions import ValidationError
from gbdraw.web_support.similarity_alignment import (
    resolve_similarity_alignment_json,
    resolve_similarity_alignment_payload,
)


BGC_SESSION = (
    Path(__file__).resolve().parents[1]
    / "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw"
)


@pytest.fixture(scope="module")
def bgc_og18(tmp_path_factory):
    """og_18 helper inputs built as the Web UI builds them: catalog members, no edges."""

    import base64

    from gbdraw.analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupResult
    from gbdraw.session_io import materialize_embedded_file
    from gbdraw.session_request_codec import decode_canonical_typed_resource
    from gbdraw.web_support.request_render import render_embedded_canonical_web_request

    session = drawing_view(read_session_payload(BGC_SESSION))
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
            "groupId": "og_18", "anchor": anchor, "identityIsUnique": True, "hidden": bool(member.get("hidden")),
            "representative": bool(member.get("representative")), "role": str(member.get("role") or ""),
        })
        anchors[gene[feature["protein_id"]]] = anchor
    paths = {
        rid: str(materialize_embedded_file(entry, temp_dir=tmp / "resources", role=rid, prefix_role=False))
        for rid, entry in session["resources"].items()
    }
    records = [{
        "recordKey": record["recordKey"],
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


def test_a_payload_without_a_projection_is_refused(bgc_og18) -> None:
    """Python places every candidate on the committed request it decodes (D5-13): a payload without one is
    malformed and is refused as the helper refuses any malformed payload."""

    from gbdraw.web_support.error_adapter import serialize_web_error

    _session, members, anchors, paths, records, tmp = bgc_og18
    request = {"schema": 2, "groupId": "og_18", "records": records, "reference": anchors["livA"],
               "members": members, "directEdges": [], "choices": []}
    with pytest.raises(ValidationError, match="projection has invalid fields") as refused:
        resolve_similarity_alignment_json(json.dumps(request), "null", json.dumps(paths), str(tmp / "helper"))
    assert serialize_web_error(refused.value, operation="resolveSimilarityAlignment", stage="helper")["code"] \
        == "VALIDATION_UNCLASSIFIED"
    with pytest.raises(ValidationError, match="projection has invalid fields"):
        resolve_similarity_alignment_payload(request, projection=None, resource_paths=paths)


def test_web_adapter_refuses_unknown_and_removed_fields(bgc_og18) -> None:
    _session, members, anchors, paths, records, tmp = bgc_og18

    def resolve(request):
        projection = {"canonicalRequest": bgc_og18[0]["renderRequest"], "orientations": None}
        return resolve_similarity_alignment_payload(
            request, projection=projection, resource_paths=paths, output_directory=str(tmp / "helper"))

    def request():
        return json.loads(json.dumps({"schema": 2, "groupId": "og_18", "records": records,
                                      "reference": anchors["livA"], "members": members,
                                      "directEdges": [], "choices": []}))

    assert resolve(request())["status"] == "resolved"
    unknown = request()
    unknown["unexpected"] = True
    with pytest.raises(ValidationError, match="request has invalid fields"):
        resolve(unknown)
    removed_choice = request()
    removed_choice["choices"] = [{"recordKey": "record-5", "kind": "skip", "anchor": None,
                                  "orientationPolicy": "match_reference"}]
    with pytest.raises(ValidationError, match=r"choices\[0\] has invalid fields"):
        resolve(removed_choice)
    # The page sends no record length: Python reads it from the committed request (D5).
    with_length = request()
    with_length["records"][1]["recordLength"] = 400
    with pytest.raises(ValidationError, match=r"records\[1\] has invalid fields"):
        resolve(with_length)


def _two_record_helper(tmp_path, *, target_strands=(1,), representative=None, reference_windows=()):
    """A reference record ``a`` (one CDS) and a target record ``b`` (one CDS per strand in ``target_strands``) as
    the page sends them, with the committed request that places them."""

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api.options import LinearDiagramOptions, LinearOutputOptions
    from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordInput
    from gbdraw.features.source import capture_source_features
    from gbdraw.session_request_codec import encode_canonical_request

    inputs, members = [], []
    for key, strands in (("a", (1,)), ("b", target_strands)):
        record = SeqRecord(Seq("ACGT" * 100), id=key, description=key)
        record.annotations.update(molecule_type="DNA", topology="linear")
        record.features = [
            SeqFeature(SimpleLocation(100 + 80 * index, 130 + 80 * index, strand=strand), type="CDS",
                       qualifiers={"protein_id": [f"{key}-{index}"]})
            for index, strand in enumerate(strands)
        ]
        path = tmp_path / f"{key}.gbk"
        SeqIO.write(record, path, "genbank")
        for index, identity in enumerate(capture_source_features(SeqIO.read(path, "genbank")).entries):
            members.append({
                "groupId": "og-1",
                "anchor": {"recordKey": key, "biologicalFeatureId": identity.biological_feature_id,
                           "sourceFeatureIndex": identity.source_feature_index,
                           "stableFeatureSvgId": identity.stable_feature_id},
                "identityIsUnique": True, "hidden": False,
                "representative": (key, index) == ("b", representative), "role": "member",
            })
        inputs.append(RecordInput(GenBankInputSource(path), record_key=key,
                                  windows=reference_windows if key == "a" else ()))
    request = LinearDiagramRequest(records=tuple(inputs), options=LinearDiagramOptions(
        selected_features_set=("CDS",), output=LinearOutputOptions(legend="none"),
        config_overrides={"canvas.show_gc": False, "canvas.show_skew": False, "objects.scale.show": False}))
    encoded = encode_canonical_request(request)
    paths = {item.resource_id: str(item.source_path) for item in encoded.resources}
    helper = {"schema": 2, "groupId": "og-1", "reference": members[0]["anchor"], "members": members,
              "directEdges": [], "choices": [],
              "records": [{"recordKey": stored["recordKey"], "windows": stored["windows"],
                           "presentation": {"reverseComplement": stored["presentation"]["reverseComplement"]}}
                          for stored in encoded.payload["records"]]}

    def resolve(payload):
        return resolve_similarity_alignment_payload(
            payload, projection={"canonicalRequest": encoded.payload, "orientations": None},
            resource_paths=paths, output_directory=str(tmp_path / "helper"))

    return helper, resolve


def test_web_adapter_discloses_candidate_facts_and_strand_relation(tmp_path) -> None:
    helper, resolve = _two_record_helper(tmp_path, target_strands=(-1,))
    helper["members"][1]["displayName"] = "target gene"
    result = resolve(helper)
    assert result["status"] == "resolved"
    row = result["records"][1]
    assert (row["status"], row["rationale"], row["reviewReason"]) == (
        "aligned", "only_usable_candidate", "only_usable_candidate")
    assert row["anchor"] == helper["members"][1]["anchor"]
    (candidate,) = row["candidates"]
    assert candidate["displayName"] == "target gene"
    assert (candidate["sourceStart"], candidate["sourceEnd"]) == (100, 130)
    assert (candidate["displayedStrand"], candidate["strandRelation"]) == (-1, "opposite")
    assert candidate["directEvidence"] == []
    assert (result["referenceDisplayedStrand"], result["referenceDisplayCenter"]) == (1, 115.0)
    assert set(result["plan"]["records"][1]) == {"recordKey", "status", "rationale", "anchor"}


def test_web_adapter_returns_ambiguity_then_uses_an_explicit_choice(tmp_path) -> None:
    helper, resolve = _two_record_helper(tmp_path, target_strands=(1, 1))
    ambiguous = resolve(helper)
    assert (ambiguous["status"], ambiguous["plan"]) == ("ambiguous", None)
    assert ambiguous["records"][1]["kind"] == "ambiguous"
    helper["choices"] = [{"recordKey": "b", "kind": "skip", "anchor": None}]
    skipped = resolve(helper)
    assert skipped["status"] == "resolved"
    assert skipped["plan"]["records"][1]["rationale"] == "skipped_by_user"


def test_web_recommendation_is_transient_and_explicit_select_is_validated(tmp_path) -> None:
    helper, resolve = _two_record_helper(tmp_path, target_strands=(-1, 1), representative=1)
    unresolved = resolve(helper)
    row = unresolved["records"][1]
    assert (row["recommendationReason"], row["recommendedAnchor"]) == (
        "unique_representative", helper["members"][2]["anchor"])
    assert [candidate["strandRelation"] for candidate in row["candidates"]] == ["opposite", "same"]
    assert unresolved["plan"] is None
    helper["choices"] = [{"recordKey": "b", "kind": "select", "anchor": helper["members"][1]["anchor"]}]
    selected = resolve(helper)
    assert selected["plan"]["records"][1]["anchor"] == helper["members"][1]["anchor"]


def test_web_adapter_refuses_a_reference_outside_the_committed_crop(tmp_path) -> None:
    from gbdraw.io.regions import RecordWindow

    helper, resolve = _two_record_helper(tmp_path, reference_windows=(RecordWindow(201, 300),))
    with pytest.raises(ValidationError, match="exact reference.*current crop"):
        resolve(helper)
