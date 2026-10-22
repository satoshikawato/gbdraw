"""Build the main-shaped Session 47 Collinear fixture with a reversed Linear record (OV-1145).

Run with origin/main's code only (no browser), from a detached worktree at
150d26706064e7ca20eef4d3553ea266b51b79f6, from the repository root:

    git worktree add --detach <main worktree> 150d26706064e7ca20eef4d3553ea266b51b79f6
    PYTHONPATH=<main worktree> python3 tools/build_main_collinear_v47_fixture.py \
        <main worktree> <out dir> majanivirus_orthogroup

It writes `majanivirus_collinear-rc.v47.gbdraw-session.json.gz`, the fixture in
`tests/fixtures/sessions/`, byte for byte (`gzip -n -9`).

Steps, mirroring main's Web Generate (app/run-analysis.js) and Save
(services/config.js exportSessionDocument, services/session-request.js
buildComparisons/addJson):
1. Start from main's Gallery Session hepatoplasmataceae_collinear (v47).
2. Keep the first two adjacent Linear records (record-1, record-2) and set
   presentation.reverseComplement on record-2.
3. Exec the Python embedded in main's app/python-helpers.js; call
   extract_cds_protein_fasta (reverse flag "0", as getSeqEntry does) and
   convert_losatp_blastp_pairs_to_genomic_payload (mode collinear, the
   Session's settings, viewTransform {length, reverse} per record), with the
   pair payloads buildLosatJobSpecs gives for collinear/adjacent/infer
   orthogroups: (0,0), (1,1), (0,1 display), (1,0), all pairIndex 0.
4. Store the canonical collinearity-result bytes as
   comparison-canonical-collinearity-1 plus generatedProteinComparison mode
   none, keep the pair's raw LOSATP entries, empty losatDerivedCache.
5. Replay the request with main's CLI (--session_output) to get the Result SVG,
   its feature catalog and runMetadata, and store them as the Web stores a
   Result.  Validate with gbdraw.session.load_session_document.
"""

from __future__ import annotations

import base64
import copy
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

MAIN = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
OUT.mkdir(parents=True, exist_ok=True)
# Optional third argument: the Gallery Session stem (default hepatoplasmataceae_collinear).
SOURCE = sys.argv[3] if len(sys.argv) > 3 else "hepatoplasmataceae_collinear"
NAME = SOURCE.replace("_orthogroup", "_collinear") + "-rc"
KEEP = ("record-1", "record-2")
REVERSED = "record-2"

src = MAIN / f"gbdraw/web/gallery/sessions/{SOURCE}.gbdraw-session.json.gz"
session = json.loads(gzip.decompress(src.read_bytes()))
assert session["version"] == 47, session["version"]
rr = session["renderRequest"]
assert rr["diagramOptions"]["featureOverrides"] == []
# No feature visibility rows or rules: extraction gets no visibility table.
assert not session["modes"]["linear"]["features"].get("featureVisibilityManualRules")

# The Python helpers of main's Web Worker, executed natively.
helpers_js = (MAIN / "gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
match = re.search(r"export const PYTHON_HELPERS = `(.*)`;\s*$", helpers_js, re.S)
assert match, "PYTHON_HELPERS template not found"
helpers_src = match.group(1)
assert "${" not in helpers_src, "template interpolation in PYTHON_HELPERS"
helpers_src = helpers_src.replace("\\`", "`").replace("\\\\", "\\")
H: dict = {"__name__": "gbdraw_web_python_helpers"}
exec(compile(helpers_src, "python-helpers.js:PYTHON_HELPERS", "exec"), H)

# 1-2. Records and resources.
records = [copy.deepcopy(r) for r in rr["records"] if r["recordKey"] in KEEP]
assert [r["recordKey"] for r in records] == list(KEEP)
for r in records:
    r["presentation"]["reverseComplement"] = r["recordKey"] == REVERSED

work = Path(tempfile.mkdtemp(prefix="ov1145-collinear-", dir=os.environ.get("TMPDIR")))
gb_paths = []
for r in records:
    res = session["resources"][r["source"]["resourceId"]]
    path = work / f"{r['recordKey']}.gb"
    path.write_bytes(base64.b64decode(res["data"]))
    gb_paths.append(path)

# 3a. Protein extraction, search frame (reverse flag "0").
manifest = session["proteinIdentityManifest"]
entries = []
for idx, (r, path) in enumerate(zip(records, gb_paths)):
    res = json.loads(H["extract_cds_protein_fasta"](
        str(path), "genbank", None, None, None, "0", idx, r["recordKey"], None, None))
    assert "error" not in res, res.get("error")
    inst = manifest["recordInstances"][r["recordKey"]]
    assert res["runtime_binding_hash"] == inst["runtimeBindingHash"], r["recordKey"]
    assert res["display_binding_hash"] == inst["displayBindingHash"], r["recordKey"]
    entries.append(res)

# 3b. Pair payloads in buildLosatJobSpecs order and the raw TSV buffer.
key_of = {r["recordKey"]: i for i, r in enumerate(records)}
raw_by_dir = {}
for e in session["losatCache"]["entries"]:
    q, s = e["queryRecordInstanceKey"], e["subjectRecordInstanceKey"]
    if q in key_of and s in key_of:
        raw_by_dir[(key_of[q], key_of[s])] = e
assert set(raw_by_dir) == {(0, 0), (1, 1), (0, 1), (1, 0)}, sorted(raw_by_dir)
spec_order = [(0, 0), (1, 1), (0, 1), (1, 0)]
pair_payloads, raw_parts, offset = [], [], 0
for q, s in spec_order:
    e = raw_by_dir[(q, s)]
    data = e["text"].encode("utf-8")
    pair_payloads.append({
        "pairIndex": 0, "ordinal": 0, "edgeKey": "", "queryIndex": q, "subjectIndex": s,
        "displayPair": (q, s) == (0, 1), "cacheKey": e["key"],
        "rawTsvOffset": offset, "rawTsvBytes": len(data),
    })
    raw_parts.append(data)
    offset += len(data)
    ids = {line.split("\t", 1)[0] for line in e["text"].splitlines() if line.strip()}
    assert ids <= set(entries[q]["protein_map"]), "raw LOSAT IDs are not protein map keys"

record_payloads = []
for idx, (r, res) in enumerate(zip(records, entries)):
    record_payloads.append({
        "recordIndex": idx,
        "recordId": res["record_id"],
        "proteinMap": res["protein_map"],
        "proteinCacheKey": res["display_binding_hash"],
        "runtimeBindingHash": res["runtime_binding_hash"],
        "displayBindingHash": res["display_binding_hash"],
        "viewTransform": {"length": res["record_length"], "reverse": r["presentation"]["reverseComplement"]},
    })
pairs_path = work / "pairs.json"
pairs_path.write_text(json.dumps({"records": record_payloads, "pairs": pair_payloads}), encoding="utf-8")
raw_path = work / "raw.tsv"
raw_path.write_bytes(b"".join(raw_parts))

# 3c. Conversion with the Session's settings (generatedProteinComparison.settings, adv thresholds).
generated = next(c for c in rr["comparisons"] if c["kind"] == "generatedProteinComparison")
settings = copy.deepcopy(generated["settings"])
params = settings["collinearityParams"]["parameters"]
opts = rr["diagramOptions"]
adv = session["modes"]["linear"]["config"]["adv"]
canonical_path = work / "canonical-comparison.json"
converted = json.loads(H["convert_losatp_blastp_pairs_to_genomic_payload"](
    str(pairs_path), str(raw_path),
    "collinear", 5,
    opts["bitscore"], str(adv.get("evalue") or opts["evalue"]), opts["identity"], opts["alignmentLength"],
    params["minAnchors"], params["maxUnitGap"], settings["collinearityUnitMode"],
    settings["collinearityColorMode"], settings["collinearityAnchorMode"],
    params["maxDiagonalDrift"], params["maxConflicts"],
    settings["collinearMaxParalogLinksPerOrthogroup"], settings["collinearitySearchScope"],
    settings["orthogroupMembershipMode"], settings["orthogroupMemberMaxHits"],
    params["mergeOrientation"], settings["collinearInferOrthogroups"],
    str(canonical_path), False,
))
assert "error" not in converted, converted.get("error")
canonical = converted["canonicalResource"]
canonical_bytes = canonical_path.read_bytes()
assert canonical["kind"] == "collinearity-result" and canonical["size"] == len(canonical_bytes)

# 4. Request and resources as main's Web Save writes them.
RID = "comparison-canonical-collinearity-1"
new_rr = copy.deepcopy(rr)
new_rr["records"] = records
new_rr["comparisons"] = [
    {"kind": "collinearityResult", "resourceId": RID, "encoding": "canonicalJson", "valueKind": "result"},
    {"kind": "generatedProteinComparison", "mode": "none", "pairs": [], "settings": settings},
]
dropped_sources = {r["source"]["resourceId"] for r in rr["records"] if r["recordKey"] not in KEEP}
resources = {
    rid: res for rid, res in session["resources"].items()
    if rid not in dropped_sources and not rid.startswith("comparison-")
}
resources[RID] = {
    "kind": "collinearity-result",
    "name": f"{RID}-{RID}.json",
    "type": "application/json",
    "size": len(canonical_bytes),
    "lastModified": 0,
    "encoding": "base64",
    "data": base64.b64encode(canonical_bytes).decode("ascii"),
}

kept_instances = {k: manifest["recordInstances"][k] for k in KEEP}
kept_analyses = {i["recordAnalysisId"]: manifest["recordAnalyses"][i["recordAnalysisId"]] for i in kept_instances.values()}
kept_sets = {a["proteinSetHash"]: manifest["proteinSets"][a["proteinSetHash"]] for a in kept_analyses.values()}
new_manifest = {"schema": manifest["schema"], "proteinSets": kept_sets,
                "recordAnalyses": kept_analyses, "recordInstances": kept_instances}

web_files = copy.deepcopy(session["webFiles"])
keep_rids = set(resources)
web_files["resourceOriginalNames"] = {k: v for k, v in web_files["resourceOriginalNames"].items() if k in keep_rids}
web_files["linearRecordMetadata"] = [m for m in web_files["linearRecordMetadata"] if m["recordKey"] in KEEP]

modes = copy.deepcopy(session["modes"])
blastp = modes["linear"]["config"]["losat"]["blastp"]
blastp["mode"] = "collinear"
# The raw entries were searched without a candidate limit (args []); keep that
# limit for the collinear mode so the stored search stays reusable.
blastp["hitLimitsByMode"]["collinear"] = {
    "candidateLimit": settings["proteinBlastpCandidateLimit"],
    "orthogroupMemberMaxHits": settings["orthogroupMemberMaxHits"],
}
blastp["candidateLimit"] = settings["proteinBlastpCandidateLimit"]
blastp["orthogroupMemberMaxHits"] = settings["orthogroupMemberMaxHits"]

base = {
    "format": "gbdraw-session",
    "version": 47,
    "createdAt": "2026-10-18T00:00:00.000Z",
    "title": NAME,
    "ui": session["ui"],
    "renderRequest": new_rr,
    "resources": resources,
    "webFiles": web_files,
    "results": [],
    "editorState": {"featureCatalog": None},
    "orthogroupState": {"selectedOrthogroupId": "", "selectedOrthogroupAlignmentFeature": ""},
    "modes": modes,
    "losatCache": {"entries": [raw_by_dir[p] for p in spec_order]},
    "losatDerivedCache": {"entries": []},
    "proteinIdentityManifest": new_manifest,
}
pre_path = work / f"{NAME}.pre.gbdraw-session.json"
pre_path.write_text(json.dumps(base), encoding="utf-8")

# 5. Main's CLI replay of this request gives the Result, its catalog and runMetadata.
prefix = OUT / f"{NAME}.main-replay"
sidecar = work / "replay-sidecar.gbdraw-session.json"
Path(f"{prefix}.svg").unlink(missing_ok=True)
env = dict(os.environ, PYTHONPATH=str(MAIN))
subprocess.run(
    [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear",
     "--session", str(pre_path), "-o", str(prefix), "-f", "svg", "--session_output", str(sidecar)],
    check=True, env=env, cwd=work,
)
main_svg = Path(f"{prefix}.svg").read_text(encoding="utf-8")
side = json.loads(sidecar.read_text(encoding="utf-8"))
assert len(side["results"]) == 1 and side["results"][0]["content"] == main_svg, "sidecar Result differs from the SVG"
result_name = new_rr["output"]["prefix"]
catalog = copy.deepcopy(side["editorState"]["featureCatalog"])
assert len(catalog["items"]) == 1
catalog["items"][0]["resultName"] = result_name
final = dict(base)
final["results"] = [{"name": result_name, "content": main_svg}]
final["editorState"] = {"featureCatalog": catalog}
if side.get("runMetadata"):
    final["runMetadata"] = side["runMetadata"]
order = ["format", "version", "createdAt", "title", "ui", "renderRequest", "resources", "webFiles",
         "results", "runMetadata", "editorState", "orthogroupState", "modes", "losatCache",
         "losatDerivedCache", "proteinIdentityManifest"]
final = {k: final[k] for k in order if k in final}
raw_json = json.dumps(final).encode("utf-8")
plain = work / f"{NAME}.v47.gbdraw-session.json"
plain.write_bytes(raw_json)
gz = OUT / f"{NAME}.v47.gbdraw-session.json.gz"
with open(gz, "wb") as handle:
    subprocess.run(["gzip", "-n", "-9", "-c", str(plain)], check=True, stdout=handle)

from gbdraw.session import load_session_document  # noqa: E402

document = load_session_document(gz)
print(json.dumps({
    "fixture": str(gz),
    "bytes": gz.stat().st_size,
    "sha256": hashlib.sha256(gz.read_bytes()).hexdigest(),
    "loadedVersion": getattr(document, "version", None),
    "records": [(r["recordKey"], e["record_id"], e["record_length"], r["presentation"]["reverseComplement"])
                for r, e in zip(records, entries)],
    "canonicalBytes": len(canonical_bytes),
    "mainReplaySvg": f"{prefix}.svg",
    "work": str(work),
}, indent=1))
