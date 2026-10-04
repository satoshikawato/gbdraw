#!/usr/bin/env python3
"""S05 evidence: og_18 direct edges reach the Web resolver from the committed resource.

Run:  PYTHONPATH=<worktree> python og18_resource_edges_proof.py [SESSION_JSON]
Builds the helper request as similarity-alignment.js does (catalog members, no
edges) and calls the Worker JSON boundary with the committed request projection.
Read-only: it writes only to a temporary directory that it deletes.
"""
from __future__ import annotations

import hashlib
import json
import sys
import tempfile
from pathlib import Path

from gbdraw.analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupResult
from gbdraw.session_io import materialize_embedded_file
from gbdraw.session_request_codec import decode_canonical_typed_resource
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from gbdraw.web_support.similarity_alignment import resolve_similarity_alignment_json

import base64

DEFAULT = Path(__file__).resolve().parents[5] / "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json"
RID, GROUP = "comparison-canonical-orthogroups-1", "og_18"

path = Path(sys.argv[1] if len(sys.argv) > 1 else DEFAULT)
data = path.read_bytes()
session = json.loads(data)
out = {"fixture": {"path": path.name, "sha256": hashlib.sha256(data).hexdigest()}}
result = decode_canonical_typed_resource(
    base64.b64decode(session["resources"][RID]["data"], validate=True),
    value_kind="orthogroupResult", expected=OrthogroupResult | OrthogroupGraphResult)
gene = {m.protein_id: m.gene or m.label for m in result.orthogroups[GROUP]}
out["resourceEdges"] = [
    {"query": gene[e.query_protein_id], "subject": gene[e.subject_protein_id], "edgeKind": e.edge_kind}
    for e in result.ortholog_edges_by_orthogroup_id[GROUP]
    if {gene[e.query_protein_id], gene[e.subject_protein_id]} & {"livA", "parA"}
    and "rac" in gene[e.query_protein_id] + gene[e.subject_protein_id]
]

with tempfile.TemporaryDirectory() as tmp:
    rendered = render_embedded_canonical_web_request(
        session["renderRequest"], resources=session["resources"], workspace=f"{tmp}/render")
    item = rendered["metadata"]["featureCatalog"]["items"][0]
    group = next(g for g in item["orthogroups"] if g["id"] == GROUP)
    bio = {(b["recordKey"], b["biologicalFeatureId"]): b for b in item["biologicalFeatures"]}
    lengths = {k: len(s["sequence"]) for k, s in zip(item["recordKeys"], item["sequenceSources"])}
    records = [{"recordKey": r["recordKey"], "recordLength": lengths.get(r["recordKey"]),
                "region": r["region"], "presentation": {
                    "reverseComplement": bool(r["presentation"]["reverseComplement"])}}
               for r in session["renderRequest"]["records"]]
    strand = {"+": 1, "-": -1, 1: 1, -1: -1}
    members, label = [], {}
    for m in group["members"]:
        b = bio[(m["recordKey"], m["biologicalFeatureId"])]
        anchor = {"recordKey": m["recordKey"], "biologicalFeatureId": m["biologicalFeatureId"],
                  "sourceFeatureIndex": b.get("sourceFeatureIndex"),
                  "stableFeatureSvgId": b.get("stableFeatureId") or m["biologicalFeatureId"]}
        members.append({"groupId": GROUP, "anchor": anchor, "sourceStart": b["start"],
                        "sourceEnd": b["end"], "sourceStrand": strand.get(b["strand"]),
                        "identityIsUnique": True, "hidden": bool(m.get("hidden")),
                        "representative": bool(m.get("representative")),
                        "role": str(m.get("role") or "")})
        label[m["biologicalFeatureId"]] = gene.get(b.get("protein_id"))
    paths = {rid: str(materialize_embedded_file(e, temp_dir=Path(tmp) / "resources", role=rid,
                                                prefix_role=False))
             for rid, e in session["resources"].items()}
    projection = {"canonicalRequest": session["renderRequest"], "orientations": None}

    def resolve(reference_gene, choices=()):
        reference = next(p["anchor"] for p in members if label[p["anchor"]["biologicalFeatureId"]] == reference_gene)
        request = {"schema": 2, "groupId": GROUP, "records": records, "reference": reference,
                   "members": members, "directEdges": [], "choices": list(choices)}
        response = json.loads(resolve_similarity_alignment_json(
            json.dumps(request), json.dumps(projection), json.dumps(paths), f"{tmp}/helper"))
        rows = []
        for r in response["records"]:
            row = {"recordKey": r["recordKey"], "kind": r["kind"]}
            if r["kind"] == "decision":
                row.update(status=r["status"], rationale=r["rationale"],
                           anchor=r["anchor"] and label[r["anchor"]["biologicalFeatureId"]])
            else:
                row.update(recommended=label[r["recommendedAnchor"]["biologicalFeatureId"]],
                           recommendationReason=r["recommendationReason"])
            if len(r["candidates"]) > 1:
                row["directEvidence"] = {label[c["anchor"]["biologicalFeatureId"]]: c["directEvidence"]
                                         for c in r["candidates"]}
            rows.append(row)
        return {"status": response["status"], "records": rows}

    out["workerJsonBoundary"] = {"livA": resolve("livA"), "parA": resolve("parA")}
    rac_l = next(p["anchor"] for p in members if label[p["anchor"]["biologicalFeatureId"]] == "racL")
    out["workerJsonBoundary"]["livA_explicit_racL"] = resolve(
        "livA", [{"recordKey": rac_l["recordKey"], "kind": "select", "anchor": rac_l}])

print(json.dumps(out, indent=1))
