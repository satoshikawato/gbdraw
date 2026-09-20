#!/usr/bin/env python3
"""Pin publicly downloadable LOSAT native assets after release approval."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gbdraw.losat_setup import RELEASE_URL, TARGETS, validate_release_lock  # noqa: E402


def download(name: str, limit: int) -> bytes:
    if Path(name).name != name or "\\" in name:
        raise ValueError("invalid release asset name")
    with urlopen(f"{RELEASE_URL}/v0.1.0/{name}", timeout=60) as response:
        if not response.geturl().startswith("https://"):
            raise ValueError("release redirected away from HTTPS")
        data = response.read(limit + 1)
    if len(data) > limit:
        raise ValueError(f"oversized release asset: {name}")
    return data


def public_release_lock(candidate_sha: str) -> dict:
    if not re.fullmatch(r"[0-9a-f]{40}", candidate_sha):
        raise ValueError("candidate SHA must be 40 lowercase hex digits")
    handoff = json.loads(download("RC-HANDOFF.json", 1024 * 1024))
    if (handoff["decision"] != "RC_HANDOFF_READY" or handoff["release"] != "v0.1.0"
            or handoff["candidate_sha"] != candidate_sha):
        raise ValueError("published handoff does not match the approved candidate")
    checksums = {}
    for line in download("SHA256SUMS", 1024 * 1024).decode().splitlines():
        digest, name = line.split("  ", 1)
        if name in checksums or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("invalid or duplicate release checksum")
        checksums[name] = digest
    artifacts = {}
    for asset in handoff["artifacts"]:
        if asset["kind"] != "native":
            continue
        target = asset["target"]
        if target not in TARGETS.values() or target in artifacts:
            raise ValueError("unexpected or duplicate native target")
        filename = asset["filename"]
        binary = "LOSAT.exe" if target.endswith("windows-msvc") else "LOSAT"
        suffix = ".zip" if binary.endswith(".exe") else ".tar.gz"
        if filename != f"LOSAT-0.1.0-{target}{suffix}" or asset["metadata"] != filename + ".metadata.json":
            raise ValueError("unexpected release asset name")
        metadata = json.loads(download(asset["metadata"], 1024 * 1024))
        if (metadata["status"] != "PASS" or metadata["candidate_sha"] != candidate_sha
                or metadata["release"] != "v0.1.0"
                or metadata["artifact"]["target"] != target
                or metadata["artifact"]["filename"] != filename
                or metadata["binary"]["filename"] != binary
                or metadata["binary"]["version"] != "losat 0.1.0"
                or metadata["extracted_artifact"]["status"] != "PASS"
                or metadata["extracted_artifact"]["binary_sha256"] != metadata["binary"]["sha256"]):
            raise ValueError("native metadata differs from the approved handoff")
        payload = download(filename, asset["size"])
        digest = hashlib.sha256(payload).hexdigest()
        if (len(payload) != asset["size"] or digest != asset["sha256"]
                or digest != checksums[filename] or digest != metadata["artifact_sha256"]):
            raise ValueError("public archive differs from the handoff/checksums")
        artifacts[target] = {
            "filename": filename, "sha256": digest, "size": len(payload),
            "binary": binary, "binary_sha256": metadata["binary"]["sha256"],
            "binary_size": metadata["binary"]["size"],
        }
    if set(artifacts) != set(TARGETS.values()):
        raise ValueError("handoff lacks one or more supported native targets")
    lock = {"schema_version": 1, "version": "0.1.0", "candidate_sha": candidate_sha, "artifacts": artifacts}
    validate_release_lock(lock)
    return lock


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-sha", required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "gbdraw/data/losat-release.json")
    args = parser.parse_args()
    # No write until every fixed public URL has been downloaded and verified.
    lock = public_release_lock(args.candidate_sha)
    args.output.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
