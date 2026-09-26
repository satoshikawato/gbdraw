"""Inspect each persisted namespace on fetched main/tag and current source history."""

import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]


def git(*args):
    return subprocess.check_output(["git", *args], cwd=root, stderr=subprocess.DEVNULL)


paths = {
    "session": "gbdraw/web/js/services/config.js",
    "request": "gbdraw/web/js/services/session-request.js",
    "bindings": "gbdraw/web/js/services/session-resources.js",
    "catalog": "gbdraw/web/js/services/feature-catalog.js",
    "cache": "gbdraw/web/js/app/losat-cache.js",
}
constants = {
    "session": ["SESSION_VERSION"],
    "request": ["CANONICAL_REQUEST_SCHEMA"],
    "catalog": ["FEATURE_CATALOG_SCHEMA"],
    "cache": [
        "PROTEIN_LOSAT_CACHE_SCHEMA",
        "NUCLEOTIDE_LOSAT_CACHE_SCHEMA",
        "LOSAT_DERIVED_CACHE_SCHEMA",
        "PROTEIN_IDENTITY_MANIFEST_SCHEMA",
    ],
}
result = {
    "mainFirstParentHead": git("rev-parse", "origin/main").decode().strip(),
    "refs": {},
    "positiveFixtures": [],
}
for ref in ["origin/main", "0.13.0", "HEAD"]:
    entry = {
        "sha": git("rev-parse", ref + "^{commit}").decode().strip(),
        "namespaces": {},
    }
    for namespace, path in paths.items():
        if not git("ls-tree", ref, "--", path).strip():
            entry["namespaces"][namespace] = {"path": path, "present": False}
            continue
        blob = git("show", ref + ":" + path)
        text = blob.decode()
        values = {
            name: int(re.search(r"\b" + name + r"\s*=\s*(\d+)", text)[1])
            for name in constants.get(namespace, [])
            if re.search(r"\b" + name + r"\s*=\s*(\d+)", text)
        }
        if namespace == "bindings":
            writer = re.search(r"const bindings = \{\s*schema:\s*(\d+)", text)
            assert writer, f"No binding writer schema witness in {ref}:{path}"
            values = {"writerSchema": int(writer[1])}
        entry["namespaces"][namespace] = {
            "path": path,
            "blob": git("rev-parse", ref + ":" + path).decode().strip(),
            "sha256": hashlib.sha256(blob).hexdigest(),
            "values": values,
        }
    result["refs"][ref] = entry
firstparent = set(
    git("rev-list", "--first-parent", "origin/main").decode().splitlines()
)
fixtures = [
    (
        "tests/fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz",
        "17e2c9dee32724219aa7c96a02d183280ffbe438",
    ),
    (
        "tests/fixtures/sessions/single.v41-bindings1.json",
        "4e8c93804186f9c4b163b584bd81d759b1b3522d",
    ),
    (
        "tests/fixtures/sessions/settings-only.v42.json.gz",
        "3548fe14baf1fe3f9bfa83dcef5f34337523e8c8",
    ),
    (
        "tests/fixtures/sessions/HmmtDNA_basic_circular.v44-schema7.json.gz",
        "4d1cf93514d0f75fa7a0ee32c1c4b2f176e03682",
    ),
]
for path, commit in fixtures:
    payload = (root / path).read_bytes()
    expanded = gzip.decompress(payload) if payload[:2] == b"\x1f\x8b" else payload
    d = json.loads(expanded)
    entry = {
        "path": path,
        "compressedSha256": hashlib.sha256(payload).hexdigest(),
        "expandedSha256": hashlib.sha256(expanded).hexdigest(),
        "witnessCommit": commit,
        "onMainFirstParent": commit in firstparent,
        "containingTags": git("tag", "--contains", commit).decode().splitlines(),
        "version": d.get("version"),
        "requestSchema": (d.get("renderRequest") or {}).get("schema"),
        "bindingsSchema": (d.get("webFiles", {}).get("bindings") or {}).get("schema"),
        "catalogSchema": (d.get("editorState", {}).get("featureCatalog") or {}).get(
            "schema"
        ),
        "cacheSchemas": sorted(
            {e.get("schema") for e in d.get("losatCache", {}).get("entries", [])}
        ),
    }
    entry["releasedWitnessVerified"] = entry["onMainFirstParent"] or bool(
        entry["containingTags"]
    )
    if not entry["releasedWitnessVerified"]:
        witnesses = (
            git("log", "--first-parent", "--format=%H", "origin/main", "--", path)
            .decode()
            .splitlines()
        )
        entry["mainFixtureWitnesses"] = witnesses
        entry["releasedFixtureBytesVerified"] = any(
            git("show", c + ":" + path) == payload for c in witnesses
        )
    result["positiveFixtures"].append(entry)
result["sourceFingerprints"] = {
    p: hashlib.sha256((root / p).read_bytes()).hexdigest()
    for p in list(paths.values())
    + [
        "gbdraw/web/js/app/run-analysis.js",
        "gbdraw/web/js/app/watchers.js",
        "gbdraw/web/js/app/record-discovery.js",
        "gbdraw/web/js/services/session-file.js",
        "gbdraw/web/js/services/svg-result-ingestion.js",
        "gbdraw/web/js/services/svg-sanitization.js",
        "tools/web-change-policy.json",
        "tools/check-web-change-budget.mjs",
        "tools/web-architecture-detectors.mjs",
        "docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md",
    ]
}
print(json.dumps(result, indent=2))
