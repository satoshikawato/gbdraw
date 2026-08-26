#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Mapping
from urllib.parse import urlparse

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
MANIFEST_PATH = GALLERY_ROOT / "artifact-manifest.json"
REMOTE_MANIFEST_PATH = Path("gallery/remote-assets.json")
FULL_SHA_RE = re.compile(r"^[0-9a-f]{40}$")


def _relative(path: Path) -> str:
    return path.relative_to(WEB_ROOT).as_posix()


def _load_examples() -> tuple[Any, ...]:
    from tools.prepare_interactive_gallery_assets import EXAMPLES

    return EXAMPLES


def gallery_artifact_paths(examples: tuple[Any, ...] | None = None) -> tuple[Path, ...]:
    examples = _load_examples() if examples is None else examples
    paths = {GALLERY_ROOT / "examples.json"}
    for example in examples:
        paths.update(
            {
                example.session_path,
                example.source_svg_path,
                example.gallery_svg_path,
                example.thumbnail_path,
            }
        )
    return tuple(sorted(paths, key=_relative))


def _file_record(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def build_gallery_artifact_manifest() -> dict[str, Any]:
    examples = _load_examples()
    return {
        "schema": 1,
        "inventory": "tools.prepare_interactive_gallery_assets.EXAMPLES",
        "exampleIds": [example.id for example in examples],
        "files": {
            _relative(path): _file_record(path)
            for path in gallery_artifact_paths(examples)
        },
    }


def write_gallery_artifact_manifest() -> None:
    MANIFEST_PATH.write_text(
        json.dumps(
            build_gallery_artifact_manifest(),
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _catalog_inventory() -> tuple[list[str], dict[Path, set[str]], set[str]]:
    catalog = _read_json(GALLERY_ROOT / "examples.json")
    if not isinstance(catalog, list) or not all(
        isinstance(item, Mapping) for item in catalog
    ):
        raise ValueError("Gallery examples.json must be an array of objects.")
    ids = [item.get("id") for item in catalog]
    sessions = [Path(str(item.get("session") or "")).name for item in catalog]
    if not all(isinstance(example_id, str) and example_id for example_id in ids):
        raise ValueError("Gallery examples.json contains an invalid example ID.")
    if len(ids) != len(set(ids)) or len(sessions) != len(set(sessions)):
        raise ValueError("Gallery examples.json contains duplicate artifacts.")
    membership = {
        GALLERY_ROOT / "sessions": set(sessions),
        GALLERY_ROOT / "sources": {f"{example_id}.svg" for example_id in ids},
        GALLERY_ROOT / "examples": {f"{example_id}.svg" for example_id in ids},
        GALLERY_ROOT / "thumbnails": {f"{example_id}.webp" for example_id in ids},
    }
    expected_paths = {"gallery/examples.json"}
    for directory, names in membership.items():
        expected_paths.update(f"{_relative(directory)}/{name}" for name in names)
    return ids, membership, expected_paths


def _expected_physical_membership(
    membership: Mapping[Path, set[str]],
) -> None:
    for directory, expected_names in membership.items():
        actual_names = {path.name for path in directory.iterdir() if path.is_file()}
        if actual_names != expected_names:
            raise ValueError(
                f"Gallery physical membership drift at {_relative(directory)}: "
                f"missing={sorted(expected_names - actual_names)}, "
                f"unexpected={sorted(actual_names - expected_names)}"
            )


def _verify_inventory(manifest: Mapping[str, Any]) -> None:
    expected_ids, membership, expected_paths = _catalog_inventory()
    if manifest.get("schema") != 1 or manifest.get("exampleIds") != expected_ids:
        raise ValueError("Gallery artifact manifest inventory does not match examples.json.")
    actual_paths = set(manifest.get("files") or {})
    if actual_paths != expected_paths:
        raise ValueError(
            "Gallery artifact manifest membership drift: "
            f"missing={sorted(expected_paths - actual_paths)}, "
            f"unexpected={sorted(actual_paths - expected_paths)}"
        )
    _expected_physical_membership(membership)


def _remote_path(url: str, commit_sha: str) -> str | None:
    parsed = urlparse(url)
    marker = f"/{commit_sha}/gbdraw/web/"
    if parsed.scheme != "https" or marker not in parsed.path:
        return None
    return parsed.path.split(marker, 1)[1]


def verify_gallery_artifacts(
    *,
    package_root: Path | None = None,
    commit_sha: str | None = None,
) -> None:
    manifest = _read_json(MANIFEST_PATH)
    if not isinstance(manifest, Mapping):
        raise ValueError("Gallery artifact manifest must be an object.")
    _verify_inventory(manifest)
    records = manifest.get("files")
    assert isinstance(records, Mapping)
    for relative_path, expected in records.items():
        checkout_path = WEB_ROOT / relative_path
        if not isinstance(expected, Mapping) or _file_record(checkout_path) != expected:
            raise ValueError(f"Gallery checkout hash mismatch: {relative_path}")
    if package_root is None:
        return
    package_root = package_root.resolve()
    packaged_manifest = package_root / MANIFEST_PATH.relative_to(WEB_ROOT)
    if packaged_manifest.read_bytes() != MANIFEST_PATH.read_bytes():
        raise ValueError("Packaged Gallery artifact manifest differs from checkout.")
    remote_manifest_path = package_root / REMOTE_MANIFEST_PATH
    remote_assets = _read_json(remote_manifest_path) if remote_manifest_path.is_file() else {}
    if not isinstance(remote_assets, Mapping):
        raise ValueError("Cloudflare remote Gallery manifest must be an object.")
    if remote_assets and (commit_sha is None or not FULL_SHA_RE.fullmatch(commit_sha)):
        raise ValueError("Remote Gallery verification requires an exact lowercase commit SHA.")
    for relative_path, expected in records.items():
        packaged_path = package_root / relative_path
        if packaged_path.is_file():
            if _file_record(packaged_path) != expected:
                raise ValueError(f"Packaged Gallery hash mismatch: {relative_path}")
            continue
        remote_url = remote_assets.get(relative_path)
        if not isinstance(remote_url, str) or _remote_path(remote_url, commit_sha or "") != relative_path:
            raise ValueError(f"Missing verified packaged Gallery artifact: {relative_path}")
    unexpected_remote = set(remote_assets) - set(records)
    if unexpected_remote:
        raise ValueError(
            f"Remote Gallery manifest contains unverified targets: {sorted(unexpected_remote)}"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Write or verify Gallery artifact hashes.")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--package-root", type=Path)
    parser.add_argument("--commit-sha")
    args = parser.parse_args(argv)
    if args.write:
        if args.package_root or args.commit_sha:
            parser.error("--write cannot be combined with package verification options")
        write_gallery_artifact_manifest()
    else:
        verify_gallery_artifacts(
            package_root=args.package_root,
            commit_sha=args.commit_sha,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
