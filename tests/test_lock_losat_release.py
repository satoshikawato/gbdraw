from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location(
    "lock_losat_release", Path(__file__).resolve().parents[1] / "tools/lock_losat_release.py"
)
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


@pytest.fixture
def published_release(monkeypatch):
    candidate = "a" * 40
    assets = []
    files = {}
    metadata = []
    for target in release.TARGETS.values():
        binary = "LOSAT.exe" if target.endswith("windows-msvc") else "LOSAT"
        suffix = ".zip" if binary.endswith(".exe") else ".tar.gz"
        filename = f"LOSAT-0.1.0-{target}{suffix}"
        payload = f"archive for {target}".encode()
        digest = hashlib.sha256(payload).hexdigest()
        binary_digest = hashlib.sha256(target.encode()).hexdigest()
        asset = {"kind": "native", "target": target, "filename": filename,
                 "metadata": filename + ".metadata.json", "size": len(payload), "sha256": digest}
        item = {"status": "PASS", "candidate_sha": candidate, "release": "v0.1.0",
                "artifact": {"target": target, "filename": filename},
                "binary": {"filename": binary, "version": "losat 0.1.0", "sha256": binary_digest, "size": 1024},
                "extracted_artifact": {"status": "PASS", "binary_sha256": binary_digest},
                "artifact_sha256": digest}
        assets.append(asset)
        metadata.append(item)
        files[filename] = payload
        files[asset["metadata"]] = item
    handoff = {"decision": "RC_HANDOFF_READY", "release": "v0.1.0", "candidate_sha": candidate, "artifacts": assets}
    files["RC-HANDOFF.json"] = handoff
    files["SHA256SUMS"] = "".join(f"{a['sha256']}  {a['filename']}\n" for a in assets).encode()

    def download(name, limit):
        value = files[name]
        data = json.dumps(value).encode() if isinstance(value, dict) else value
        assert len(data) <= limit
        return data
    monkeypatch.setattr(release, "download", download)
    return candidate, files, handoff, metadata


def test_public_identity_generates_complete_lock(published_release):
    candidate, _, handoff, _ = published_release
    lock = release.public_release_lock(candidate)
    assert lock["candidate_sha"] == candidate
    assert set(lock["artifacts"]) == set(release.TARGETS.values())
    assert {entry["sha256"] for entry in lock["artifacts"].values()} == {a["sha256"] for a in handoff["artifacts"]}


@pytest.mark.parametrize("mutation", ["candidate", "decision", "missing", "duplicate", "checksum", "payload", "metadata", "binary_hash", "binary_size"])
def test_inconsistent_public_release_cannot_generate_lock(published_release, mutation):
    candidate, files, handoff, metadata = published_release
    if mutation == "candidate":
        handoff["candidate_sha"] = "b" * 40
    elif mutation == "decision":
        handoff["decision"] = "FAILED"
    elif mutation == "missing":
        handoff["artifacts"].pop()
    elif mutation == "duplicate":
        handoff["artifacts"].append(handoff["artifacts"][0])
    elif mutation == "checksum":
        files["SHA256SUMS"] = files["SHA256SUMS"].replace(handoff["artifacts"][0]["sha256"].encode(), b"0" * 64)
    elif mutation == "payload":
        name = handoff["artifacts"][0]["filename"]
        files[name] = b"x" * len(files[name])
    elif mutation == "metadata":
        metadata[0]["candidate_sha"] = "b" * 40
    elif mutation == "binary_hash":
        metadata[0]["binary"]["sha256"] = metadata[0]["extracted_artifact"]["binary_sha256"] = "invalid"
    elif mutation == "binary_size":
        metadata[0]["binary"]["size"] = -1
    with pytest.raises(ValueError):
        release.public_release_lock(candidate)


def test_generator_does_not_overwrite_existing_lock_on_failure(published_release, tmp_path, monkeypatch):
    candidate, _, handoff, _ = published_release
    handoff["candidate_sha"] = "b" * 40
    output = tmp_path / "lock.json"
    output.write_text("existing lock\n")
    monkeypatch.setattr(release.sys, "argv", ["lock_losat_release.py", "--candidate-sha", candidate, "--output", str(output)])
    with pytest.raises(ValueError):
        release.main()
    assert output.read_text() == "existing lock\n"
