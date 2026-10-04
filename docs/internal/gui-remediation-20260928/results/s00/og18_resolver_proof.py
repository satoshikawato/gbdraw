#!/usr/bin/env python3
"""Evidence: og_18 (livA anchor) Similarity Group alignment, BGC0000708-BGC0000713 session.

Run:  PYTHONPATH=<gbdraw snapshot> python og18_resolver_proof.py [SESSION_JSON]
(dev snapshot: shared resolver; main snapshot: reports its legacy representative selector)
Read-only: it writes only to a temporary directory that it deletes.
"""
from __future__ import annotations

import base64
import hashlib
import json
import sys
import tempfile
from pathlib import Path

import gbdraw
from gbdraw.analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupResult
from gbdraw.session_io import materialize_embedded_file
from gbdraw.session_request_codec import decode_canonical_typed_resource
from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload
from gbdraw.web_support.request_render import render_embedded_canonical_web_request

DEFAULT = Path(__file__).resolve().parents[5] / "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json"
RID, GROUP, ANCHOR_GENE = "comparison-canonical-orthogroups-1", "og_18", "livA"

path = Path(sys.argv[1] if len(sys.argv) > 1 else DEFAULT)
data = path.read_bytes()
session = json.loads(data)
out = {
    "python": sys.version.split()[0],
    "gbdraw.__file__": gbdraw.__file__,
    "fixture": {"path": str(path), "sha256": hashlib.sha256(data).hexdigest(),
                "gitBlob": hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest(),
                "sessionVersion": session["version"]},
}

# 1. Session resource: base64 envelope around canonical typed JSON bytes.
entry = session["resources"][RID]
raw = base64.b64decode(entry["data"], validate=True)
body = json.loads(raw)
out["resource"] = {"id": RID, "kind": entry["kind"], "envelopeEncoding": entry["encoding"],
                   "declaredSize": entry["size"], "decodedBytes": len(raw),
                   "decodedSha256": hashlib.sha256(raw).hexdigest(),
                   "typedKind": body["kind"], "typedSchema": body["schema"]}
result = decode_canonical_typed_resource(raw, value_kind="orthogroupResult",
                                         expected=OrthogroupResult | OrthogroupGraphResult)
out["resource"]["decodedType"] = type(result).__name__

# 2. og_18 members, memberByProteinId, and direct edges.
record_keys = [r["recordKey"] for r in session["renderRequest"]["records"]]
members = result.orthogroups[GROUP]
gene = {m.protein_id: (m.gene or m.label) for m in members}
out["og18Members"] = [{
    "gene": m.gene, "proteinId": m.protein_id, "sourceProteinId": m.source_protein_id,
    "recordIndex": m.record_index, "recordKey": record_keys[m.record_index],
    "recordId": m.record_id,
    "sourceResource": session["renderRequest"]["records"][m.record_index]["source"]["resourceId"],
    "featureIndex": m.feature_index, "featureSvgId": m.feature_svg_id,
    "start": m.start, "end": m.end, "strand": m.strand, "representative": m.representative,
    "role": m.role, "assignmentReason": m.assignment_reason,
    "memberByProteinIdGroup": result.member_by_protein_id[m.protein_id].orthogroup_id,
} for m in members]
edges = result.ortholog_edges_by_orthogroup_id.get(GROUP, ())
livA = next(m for m in members if m.gene == ANCHOR_GENE)
fmt = lambda e: {"query": gene[e.query_protein_id], "subject": gene[e.subject_protein_id],
                 "edgeKind": e.edge_kind, "renderRole": e.render_role, "identity": e.identity,
                 "bitscore": e.bitscore}
out["og18Edges"] = {
    "count": len(edges),
    "livA": [fmt(e) for e in edges if livA.protein_id in (e.query_protein_id, e.subject_protein_id)],
    "parA_to_record5": [fmt(e) for e in edges if gene[e.query_protein_id] == "parA"
                        and e.subject_record_index == 4],
    "relatedEdges": len(result.related_edges_by_orthogroup_id.get(GROUP, ())),
}

# 3. Catalog: serializer keeps orthologEdges; the Web feature catalog emits only a count.
serialized = next(g for g in serialize_orthogroups_payload(result) if g["id"] == GROUP)
with tempfile.TemporaryDirectory() as tmp:
    rendered = render_embedded_canonical_web_request(
        session["renderRequest"], resources=session["resources"], workspace=f"{tmp}/render")
    item = rendered["metadata"]["featureCatalog"]["items"][0]
    cat_group = next(g for g in item["orthogroups"] if g["id"] == GROUP)
    saved_items = ((session.get("editorState") or {}).get("featureCatalog") or {}).get("items") or [{}]
    saved_group = next((g for g in saved_items[0].get("orthogroups", []) if g.get("id") == GROUP), {})
    out["catalog"] = {
        "serializedPayload.orthologEdges": len(serialized["orthologEdges"]),
        "freshCatalog.hasOrthologEdges": "orthologEdges" in cat_group,
        "freshCatalog.orthologEdgeCount": cat_group.get("orthologEdgeCount"),
        "freshCatalog.orthologPathCount": cat_group.get("orthologPathCount"),
        "freshCatalog.groupKeys": sorted(cat_group),
        "savedCatalog.hasOrthologEdges": "orthologEdges" in saved_group,
        "savedCatalog.orthologEdgeCount": saved_group.get("orthologEdgeCount"),
    }

    # 4. Helper request built as similarity-alignment.js buildHelperRequest does.
    bio = {(b["recordKey"], b["biologicalFeatureId"]): b for b in item["biologicalFeatures"]}
    lengths = {k: len(s["sequence"]) for k, s in zip(item["recordKeys"], item["sequenceSources"])}
    records = [{"recordKey": r["recordKey"], "recordLength": lengths.get(r["recordKey"]),
                "region": r["region"], "presentation": {
                    "reverseComplement": bool(r["presentation"]["reverseComplement"])}}
               for r in session["renderRequest"]["records"]]
    strand = {"+": 1, "-": -1, 1: 1, -1: -1}
    payload_members, endpoint, label = [], {}, {}
    for m in cat_group["members"]:
        b = bio[(m["recordKey"], m["biologicalFeatureId"])]
        anchor = {"recordKey": m["recordKey"], "biologicalFeatureId": m["biologicalFeatureId"],
                  "sourceFeatureIndex": b.get("sourceFeatureIndex"),
                  "stableFeatureSvgId": b.get("stableFeatureId") or m["biologicalFeatureId"]}
        payload_members.append({"groupId": GROUP, "anchor": anchor, "sourceStart": b["start"],
                                "sourceEnd": b["end"], "sourceStrand": strand.get(b["strand"]),
                                "identityIsUnique": True, "hidden": bool(m.get("hidden")),
                                "representative": bool(m.get("representative")),
                                "role": str(m.get("role") or "")})
        record_index = item["recordKeys"].index(m["recordKey"])
        for pid in (b.get("protein_id"), b.get("source_protein_id")):
            endpoint[(record_index, pid)] = anchor
        label[m["biologicalFeatureId"]] = gene.get(b.get("protein_id"))
    keys = [(p["anchor"]["recordKey"], p["anchor"]["biologicalFeatureId"]) for p in payload_members]
    for p, k in zip(payload_members, keys):
        p["identityIsUnique"] = keys.count(k) == 1
    reference = next(p["anchor"] for p in payload_members if label[p["anchor"]["biologicalFeatureId"]] == ANCHOR_GENE)
    saved_edges = [{"groupId": GROUP,
                    "query": endpoint[(e["queryRecordIndex"], e["queryProteinId"])],
                    "subject": endpoint[(e["subjectRecordIndex"], e["subjectProteinId"])],
                    "edgeKind": e["edgeKind"]} for e in serialized["orthologEdges"]]
    request = {"schema": 2, "groupId": GROUP, "records": records, "reference": reference,
               "members": payload_members, "choices": []}

    # The UI sends group.orthologEdges; the catalog group (spread by expandOrthogroup) has none.
    ui_raw = cat_group.get("orthologEdges") if isinstance(cat_group.get("orthologEdges"), list) else []
    ui_edges = [{"groupId": GROUP,
                 "query": endpoint[(e["queryRecordIndex"], e["queryProteinId"])],
                 "subject": endpoint[(e["subjectRecordIndex"], e["subjectProteinId"])],
                 "edgeKind": e["edgeKind"]} for e in ui_raw]
    out["uiInputs"] = {"reference": reference, "directEdgeCountUiWouldSend": len(ui_edges),
                       "directEdgeCountAvailable": len(saved_edges)}
    try:
        from gbdraw.layout import similarity_alignment as sa
        from gbdraw.web_support import similarity_alignment as web
    except ImportError as exc:  # main snapshot: no shared resolver; report its legacy selector
        from gbdraw.diagrams.linear import orthogroup_alignment as legacy
        by_group = legacy._collect_alignment_members_from_orthogroups(result)
        picks = {}
        for target in (GROUP, livA.protein_id):
            gid, anchor_member = legacy._resolve_target_member(by_group, target)
            chosen = {}  # transcription of calculate_orthogroup_alignment_offsets per-record loop
            for m in by_group[gid]:
                cur = chosen.get(m.record_index)
                if m.representative and (cur is None or not cur.representative):
                    chosen[m.record_index] = m
                elif cur is None or (not cur.representative and legacy._is_better_member(m, cur)):
                    chosen[m.record_index] = m
            chosen[anchor_member.record_index] = anchor_member
            picks[target] = {"anchor": gene[anchor_member.protein_id],
                             "perRecord": {record_keys[i]: gene[m.protein_id] for i, m in sorted(chosen.items())}}
        out["resolver"] = {"available": False, "importError": str(exc), "legacySelector": picks}
        print(json.dumps(out, indent=1))
        raise SystemExit(0)

    # 4a. Pure resolver, candidates from the Web adapter's own member projection.
    facts = {r["recordKey"]: web._record_fact(r, "r") for r in records}
    candidates = [web._member_candidate(p, "m", group_id=GROUP, record_facts=facts)
                  for p in payload_members]
    A = lambda a: web._anchor(a, "a")
    evidence = [sa.AlignmentEvidenceEdge(GROUP, A(e["query"]), A(e["subject"]), e["edgeKind"])
                for e in saved_edges]

    def summarize(resolution):
        rows = []
        for rec, row in zip(resolution.records, resolution.review_rows):
            entry = {"recordKey": rec.record_key}
            if isinstance(rec, sa.AmbiguousAlignmentRecord):
                entry.update(kind="ambiguous",
                             candidates=[label[c.anchor.biological_feature_id] for c in rec.candidates],
                             directRbhCandidates=[label[a.biological_feature_id] for a in rec.direct_rbh_candidates],
                             recommended=label[rec.recommended_anchor.biological_feature_id],
                             recommendationReason=rec.recommendation_reason.value)
            else:
                entry.update(kind="decision", status=rec.status.value, rationale=rec.rationale.value,
                             anchor=rec.anchor and label[rec.anchor.biological_feature_id])
            if len(row.candidates) > 1:
                entry["directEvidence"] = {label[c.candidate.anchor.biological_feature_id]:
                                           list(c.direct_evidence) for c in row.candidates}
            rows.append(entry)
        return {"ambiguityCount": len(resolution.ambiguities), "planResolved": resolution.plan is not None,
                "records": rows}

    out["pureResolver"] = {
        "a_noEdges": summarize(sa.resolve_similarity_alignment(
            record_keys=record_keys, group_id=GROUP, reference=A(reference),
            candidates=candidates, edges=())),
        "b_savedOg18Edges": summarize(sa.resolve_similarity_alignment(
            record_keys=record_keys, group_id=GROUP, reference=A(reference),
            candidates=candidates, edges=evidence)),
    }

    # 4b. Same request through the Worker JSON boundary, with the renderer projection.
    res_dir = Path(tmp) / "resources"
    paths = {rid: str(materialize_embedded_file(e, temp_dir=res_dir, role=rid, prefix_role=False))
             for rid, e in session["resources"].items()}
    orient = {r["recordKey"]: bool((r["region"] or {}).get("reverseComplement")
                                   if r["region"] else r["presentation"]["reverseComplement"])
              for r in session["renderRequest"]["records"]}
    projection = {"canonicalRequest": session["renderRequest"], "orientations": orient}

    def boundary(direct_edges):
        response = json.loads(web.resolve_similarity_alignment_json(
            json.dumps({**request, "directEdges": direct_edges}), json.dumps(projection),
            json.dumps(paths), f"{tmp}/helper"))
        rows = []
        for r in response["records"]:
            row = {"recordKey": r["recordKey"], "kind": r["kind"]}
            if r["kind"] == "decision":
                row.update(status=r["status"], rationale=r["rationale"],
                           anchor=r["anchor"] and label[r["anchor"]["biologicalFeatureId"]])
            else:
                row.update(directRbhCandidates=[label[a["biologicalFeatureId"]] for a in r["directRbhCandidates"]],
                           recommended=label[r["recommendedAnchor"]["biologicalFeatureId"]],
                           recommendationReason=r["recommendationReason"])
            if len(r["candidates"]) > 1:
                row["directEvidence"] = {label[c["anchor"]["biologicalFeatureId"]]: c["directEvidence"]
                                         for c in r["candidates"]}
            rows.append(row)
        return {"status": response["status"], "planIsNull": response["plan"] is None, "records": rows}

    out["workerJsonBoundary"] = {"a_directEdgesAsUiSends": boundary(ui_edges),
                                 "b_directEdgesFromSavedOg18": boundary(saved_edges)}

print(json.dumps(out, indent=1))
