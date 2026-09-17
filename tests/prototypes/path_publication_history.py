"""Read-only S02 namespace/publication evidence from Git objects and fixtures.

Run from any directory: python tests/prototypes/path_publication_history.py --output FILE
No historical source is executed. No user or Gallery artifacts are rewritten.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    "gbdraw/analysis/protein_colinearity.py", "gbdraw/api/__init__.py",
    "gbdraw/session_request_codec.py", "gbdraw/session_io.py",
    "gbdraw/analysis/protein_artifacts.py", "gbdraw/web_support/feature_catalog.py",
    "gbdraw/render/interactive_svg.py", "gbdraw/web/js/app/losat-cache.js",
    "gbdraw/web/js/services/config.js",
)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def source_snapshot(ref):
    result = {}
    paths = set(git("ls-tree", "-r", "--name-only", ref, *SOURCE_PATHS).decode().splitlines())
    for path in SOURCE_PATHS:
        if path not in paths:
            continue
        raw = git("show", f"{ref}:{path}")
        text = raw.decode()
        constants = re.findall(
            r"^(?:(?:export )?const )?[A-Z_]*(?:SCHEMA|VERSION)[A-Z_]*\s*=\s*[^\n]+", text, re.M)
        typed = re.search(r'def (?:encode_canonical_typed_resource|_typed_json_resource).*?"schema": (\d+)', text, re.S)
        result[path] = {
            "sha256": sha(raw), "schemaDeclarations": constants,
            "typedWriterSchema": int(typed[1]) if typed else None,
            "hasPathClassOrExport": '"OrthologPath"' in text or "class OrthologPath:" in text,
            "hasPathArray": "orthologPaths" in text or "ortholog_paths_by_orthogroup_id" in text,
        }
    return result


def artifact(raw):
    encoded_hash = sha(raw)
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    data = json.loads(raw)
    observations = {}

    def visit(value, path):
        if isinstance(value, dict):
            for key, child in value.items():
                if key in {"orthologPaths", "orthologPathsByOrthogroupId", "orthologPathCount"}:
                    observations[path + "." + key] = {
                        "type": type(child).__name__, "size": len(child) if isinstance(child, (list, dict)) else child,
                        "sha256": sha(json.dumps(child, sort_keys=True, separators=(",", ":")).encode()),
                    }
                visit(child, path + "." + key)
        elif isinstance(value, list):
            for index, child in enumerate(value):
                visit(child, f"{path}[{index}]")

    visit(data, "$")
    typed = {}
    for key, value in data.get("resources", {}).items():
        if value.get("kind") not in {"orthogroup-result", "collinearity-result"}:
            continue
        body = base64.b64decode(value["data"], validate=True)
        resource = json.loads(body)
        typed[key] = {"sha256": sha(body), "schema": resource["schema"],
                      "kind": resource["kind"], "type": resource["value"].get("type")}
        visit(resource, f"decoded-resource:{key}")
    return {
        "encodedSha256": encoded_hash, "decodedSha256": sha(raw),
        "sessionVersion": data.get("version"), "requestSchema": (data.get("renderRequest") or {}).get("schema"),
        "typedResources": typed,
        "rawSchemas": sorted({x["schema"] for x in data.get("losatCache", {}).get("entries", [])}),
        "derivedSchemas": sorted({x["schema"] for x in data.get("losatDerivedCache", {}).get("entries", [])}),
        "manifestSchema": data.get("proteinIdentityManifest", {}).get("schema"),
        "catalogSchema": (data.get("editorState", {}).get("featureCatalog") or {}).get("schema"),
        "bindingsSchema": data.get("webFiles", {}).get("bindings", {}).get("schema"),
        "pathFields": observations,
        "svgPathIdAttributes": sum(x.get("content", "").count("data-ortholog-path-id=") for x in data.get("results", [])),
    }


def collect():
    refs = git("log", "--first-parent", "--format=%H", "origin/main", "--", *SOURCE_PATHS).decode().splitlines()
    snapshots = {ref: source_snapshot(ref) for ref in reversed(refs)}
    tags = git("tag", "--list", "--sort=version:refname").decode().splitlines()
    # Latest release containing the class; inspect every tag's class/export presence.
    tag_evidence = {}
    for tag in tags:
        tree = set(git("ls-tree", "-r", "--name-only", tag, SOURCE_PATHS[0]).decode().splitlines())
        text = git("show", f"{tag}:{SOURCE_PATHS[0]}").decode() if tree else ""
        tag_evidence[tag] = {"commit": git("rev-parse", f"{tag}^{{commit}}").decode().strip(),
                             "hasPathClass": "class OrthologPath:" in text}
    session = "gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json"
    historical = {}
    for ref in ("6b89c781", "17e2c9de", "3bca0e8d"):
        full = git("rev-parse", ref).decode().strip()
        assert full in set(git("rev-list", "--first-parent", "origin/main").decode().splitlines())
        historical[f"{full}:{session}"] = artifact(git("show", f"{full}:{session}"))
    historical["0.13.0:gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json"] = artifact(
        git("show", "0.13.0:gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json"))
    fixtures = {}
    for path in (
        "tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz",
        "tests/fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz",
        "tests/fixtures/sessions/settings-only.v42.json.gz",
        "tests/fixtures/sessions/single.v41-bindings1.json",
        "gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
        "gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    ):
        fixtures[path] = artifact((ROOT / path).read_bytes())
    # Current Gallery bytes are themselves present in first-parent main.
    current_public = {}
    for path in fixtures:
        if path.startswith("gbdraw/web/gallery/"):
            public = git("show", f"origin/main:{path}")
            assert sha(public) == fixtures[path]["encodedSha256"]
            current_public[path] = sha(public)
    return {
        "main": git("rev-parse", "origin/main").decode().strip(),
        "dev": git("rev-parse", "origin/dev").decode().strip(),
        "firstParentSourceSnapshots": snapshots, "releaseTags": tag_evidence,
        "historicalPositiveArtifacts": historical, "workingPositiveFixtures": fixtures,
        "currentGalleryBytesEqualMain": current_public,
        "limits": "Source snapshots prove main inclusion, not hosted deployment or a release that was never tagged. Fixture provenance is verified by bytes, not README assertions.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    result = collect()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(gzip.compress(json.dumps(result, sort_keys=True, separators=(",", ":")).encode(), mtime=0))
