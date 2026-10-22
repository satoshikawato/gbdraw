"""Build a main-shaped Session 39 or early Session 40 fixture with a stored protein table of a reversed record.

OV-1181 / OV-1183: main's Web writers of Session 39, and of Session 40 before
3bca0e8d9 (#319), stored a reversed record's drawn-frame feature ID in the
table's `*_feature_svg_id` with no `*_view_feature_svg_id`.

Run with origin/main's code only (no browser), from a detached worktree at a
first-parent main commit, from the repository root:

    git worktree add --detach <main worktree> <commit>
    PYTHONPATH=<main worktree> python3 tools/build_main_pairwise_rc_fixture.py \
        <main worktree> <out dir>

With cdd310134 (#317, Session 39) it writes
`BGC0000708-BGC0000709-rc.v39.gbdraw-session.json.gz`; with 7aad9e3e0
(Session 40, before 3bca0e8d9) `BGC0000708-BGC0000709-rc.v40-early.gbdraw-session.json.gz`.
Each is the fixture in `tests/fixtures/sessions/`, byte for byte (`gzip -n -9`).

Steps, mirroring main's Web Generate (app/run-analysis.js, pairwise blastp) and
Save (services/config.js exportSession, services/session-request.js
buildComparisons/addCanonicalTable):
1. Start from main's Gallery Session BGC0000708-BGC0000713 (v39 or v40).
2. Keep the first two Linear records and set presentation.reverseComplement on
   record-2.
3. Exec the Python embedded in main's app/python-helpers.js; call
   extract_cds_protein_fasta (reverse flag "0", as getSeqEntry does) and
   convert_losatp_blastp_pairs_to_genomic_payload (mode pairwise, the
   Session's thresholds, viewTransform {length, reverse} per record) with the
   one job spec pairwise mode gives two records: (0, 1), pairIndex 0 (the
   one adjacent edge, ordinal 0, in Session 40), and the stored raw LOSATP
   entry for that pair.
4. Store the converted rows as main's addCanonicalTable writes them (main's own
   JavaScript: normalizeResourceName, addText and addCanonicalTable, run by
   node with Buffer for textToBase64) under comparison-resolved-1, with
   generatedProteinComparison mode none; keep the pair's raw LOSATP entry and
   an empty losatDerivedCache.
5. Replay the request with main's CLI (--session_output) and keep that
   sidecar's Result, features, editor state and runMetadata, as a Web Save
   after that Generate holds them.  Validate with
   gbdraw.session.load_session_document.
"""

from __future__ import annotations

import base64
import copy
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
NAME = "BGC0000708-BGC0000709-rc"
KEEP = ("record-1", "record-2")
REVERSED = "record-2"
RID = "comparison-resolved-1"

src = MAIN / "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json"
session = json.loads(src.read_bytes())
VERSION = session["version"]
assert VERSION in (39, 40), VERSION
SUFFIX = "v39" if VERSION == 39 else "v40-early"
rr = session["renderRequest"]

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

work = Path(tempfile.mkdtemp(prefix="ov1181-v39-", dir=os.environ.get("TMPDIR")))
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
        str(path), "genbank", None, None, "", "0", idx, r["recordKey"], None))
    assert "error" not in res, res.get("error")
    inst = manifest["recordInstances"][r["recordKey"]]
    assert res["runtime_binding_hash"] == inst["runtimeBindingHash"], r["recordKey"]
    assert res["display_binding_hash"] == inst["displayBindingHash"], r["recordKey"]
    entries.append(res)

# 3b. The pairwise job spec (0, 1) and its stored raw LOSATP entry.
raw = next(
    e for e in session["losatCache"]["entries"]
    if (e["queryRecordInstanceKey"], e["subjectRecordInstanceKey"]) == KEEP
)
ids = {line.split("\t", 1)[0] for line in raw["text"].splitlines() if line.strip()}
assert ids <= set(entries[0]["protein_map"]), "raw LOSAT IDs are not protein map keys"
record_payloads = [{
    "recordIndex": idx,
    "recordId": res["record_id"],
    "proteinMap": res["protein_map"],
    "proteinCacheKey": res["display_binding_hash"],
    "runtimeBindingHash": res["runtime_binding_hash"],
    "displayBindingHash": res["display_binding_hash"],
    "viewTransform": {"length": res["record_length"], "reverse": r["presentation"]["reverseComplement"]},
} for idx, (r, res) in enumerate(zip(records, entries))]
pair_payloads = [{
    "pairIndex": 0, "queryIndex": 0, "subjectIndex": 1, "cacheKey": raw["key"], "blastText": raw["text"],
}]

# 3c. Conversion with the Session's thresholds (run-analysis.js passes adv and losat.blastp).
adv = session["config"]["adv"]
blastp = session["config"]["losat"]["blastp"]
converted = json.loads(H["convert_losatp_blastp_pairs_to_genomic_payload"](
    json.dumps({"records": record_payloads, "pairs": pair_payloads}),
    "pairwise", max(1, blastp["maxHits"]),
    adv["min_bitscore"], adv["evalue"], adv["identity"], adv["alignment_length"],
    blastp["collinearMinAnchors"], blastp["collinearMaxUnitGap"], "cds",
    blastp["collinearColorMode"], "rbh", blastp["collinearMaxDiagonalDrift"],
    blastp["collinearMaxConflictsInMergeGap"], blastp["collinearMaxParalogLinksPerOrthogroup"],
    blastp["collinearSearchScope"], blastp["orthogroupMembershipMode"], blastp["orthogroupMemberMaxHits"],
))
assert "error" not in converted, converted.get("error")
[pair] = converted["pairs"]
assert pair["pair_index"] == 0 and pair["rows"], "no converted rows"

# 4. main's addCanonicalTable / normalizeResourceName, run by node on the JSON rows.
request_js = (MAIN / "gbdraw/web/js/services/session-request.js").read_text(encoding="utf-8")
name_fn = re.search(r"const normalizeResourceName = .*?\n};\n", request_js, re.S).group(0)
text_fn = re.search(r"  const addText = .*?\n  };\n", request_js, re.S).group(0)
table_fn = re.search(r"  const addCanonicalTable = .*?\n  };\n", request_js, re.S).group(0)
rows_path = work / "rows.json"
rows_path.write_text(json.dumps(pair["rows"]), encoding="utf-8")
node_src = "\n".join([
    "const fs = require('fs');",
    name_fn,
    "const textToBase64 = (text) => Buffer.from(text, 'utf8').toString('base64');",
    "const resources = {};",
    text_fn,
    table_fn,
    f"addCanonicalTable({json.dumps(RID)}, JSON.parse(fs.readFileSync(process.argv[1], 'utf8')), null);",
    "process.stdout.write(JSON.stringify(resources));",
])
table_resource = json.loads(subprocess.run(
    ["node", "-e", node_src, str(rows_path)], check=True, capture_output=True, text=True,
).stdout)[RID]

# Request, resources and Web state as main's Web Save writes them after that Generate.
generated = next(c for c in rr["comparisons"] if c["kind"] == "generatedProteinComparison")
new_rr = copy.deepcopy(rr)
new_rr["records"] = records
new_rr["comparisons"] = [
    {"kind": "precomputedProteinComparison", "resourceId": RID, "encoding": "canonicalTsv",
     "queryRecordIndex": 0, "subjectRecordIndex": 1},
    # generatedProteinSettings: no orthogroup is selected after a pairwise Generate.
    {"kind": "generatedProteinComparison", "mode": "none", "pairs": [],
     "settings": {**generated["settings"], "alignOrthogroupFeature": None}},
]
dropped_sources = {r["source"]["resourceId"] for r in rr["records"] if r["recordKey"] not in KEEP}
resources = {
    rid: res for rid, res in session["resources"].items()
    if rid not in dropped_sources and not rid.startswith("comparison-")
}
resources[RID] = table_resource

kept_instances = {k: manifest["recordInstances"][k] for k in KEEP}
kept_analyses = {i["recordAnalysisId"]: manifest["recordAnalyses"][i["recordAnalysisId"]] for i in kept_instances.values()}
kept_sets = {a["proteinSetHash"]: manifest["proteinSets"][a["proteinSetHash"]] for a in kept_analyses.values()}
new_manifest = {"schema": manifest["schema"], "proteinSets": kept_sets,
                "recordAnalyses": kept_analyses, "recordInstances": kept_instances}

web_files = copy.deepcopy(session["webFiles"])
web_files["resourceOriginalNames"] = {k: v for k, v in web_files["resourceOriginalNames"].items() if k in resources}
web_files["linearRecordMetadata"] = [m for m in web_files["linearRecordMetadata"] if m["recordKey"] in KEEP]

config = copy.deepcopy(session["config"])
config["losat"]["blastp"]["mode"] = "pairwise"

base = {
    "format": "gbdraw-session",
    "version": VERSION,
    "createdAt": "2026-10-18T00:00:00.000Z",
    "title": NAME,
    "config": config,
    "ui": session["ui"],
    "renderRequest": new_rr,
    "resources": resources,
    "webFiles": web_files,
    "results": [],
    # Session 40 holds no derived feature payloads (featureRecordIds).
    "features": {"selectedFeatureRecordIdx": 0, **({"featureRecordIds": []} if VERSION == 39 else {})},
    "editorState": {} if VERSION == 39 else {"featureCatalog": None},
    # Session 40 holds no derived orthogroup groups.
    "orthogroupState": {**({"groups": []} if VERSION == 39 else {}), "selectedOrthogroupId": "",
                        "selectedOrthogroupAlignmentFeature": "",
                        "orthogroupNameOverrides": {}, "orthogroupDescriptionOverrides": {}},
    "losatCache": {"entries": [raw]},
    "losatDerivedCache": {"entries": []},
    "proteinIdentityManifest": new_manifest,
}
pre_path = work / f"{NAME}.pre.gbdraw-session.json"
pre_path.write_text(json.dumps(base), encoding="utf-8")

# 5. Main's CLI replay of this request gives the Result, features, editor state and runMetadata.
prefix = OUT / f"{NAME}.{SUFFIX}.main-replay"
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
final = dict(base)
final["results"] = side["results"]
for key in ("features", "editorState", "runMetadata"):
    if side.get(key):
        final[key] = side[key]
order = ["format", "version", "createdAt", "title", "config", "ui", "renderRequest", "resources", "webFiles",
         "results", "features", "editorState", "orthogroupState", "runMetadata", "losatCache",
         "losatDerivedCache", "proteinIdentityManifest"]
final = {k: final[k] for k in order if k in final}
plain = work / f"{NAME}.{SUFFIX}.gbdraw-session.json"
plain.write_bytes(json.dumps(final).encode("utf-8"))
gz = OUT / f"{NAME}.{SUFFIX}.gbdraw-session.json.gz"
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
    "rows": len(pair["rows"]),
    "sidecarKeys": sorted(side),
    "mainReplaySvg": f"{prefix}.svg",
    "work": str(work),
}, indent=1))
