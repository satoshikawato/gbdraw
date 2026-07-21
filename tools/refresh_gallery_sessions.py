#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from xml.etree import ElementTree as ET


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_io import (  # noqa: E402
    CURRENT_SESSION_VERSION,
    PROTEIN_IDENTITY_MANIFEST_SCHEMA,
    classify_raw_losat_cache_entry,
    load_session,
    session_mode,
    validate_session,
    write_session_json,
)
from gbdraw.render.formats import INTERACTIVE_SVG_FORMAT  # noqa: E402


GALLERY_ROOT = REPO_ROOT / "gbdraw" / "web" / "gallery"
SESSION_ROOT = GALLERY_ROOT / "sessions"
SESSION_PROMOTER = REPO_ROOT / "tools" / "promote_gallery_session.mjs"

GALLERY_SESSION_FILES = (
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_basic_circular.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "tobacco-chloroplast.gbdraw-session.json",
    "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "lambda_basic_linear.gbdraw-session.json",
)


def _validate_gallery_session_inventory() -> None:
    """Keep the physical, refresh-tool, and Gallery example inventories aligned."""

    from tools.prepare_interactive_gallery_assets import EXAMPLES

    configured = set(GALLERY_SESSION_FILES)
    if len(configured) != len(GALLERY_SESSION_FILES):
        raise ValueError("GALLERY_SESSION_FILES contains duplicate session names")
    physical = {
        path.name
        for pattern in ("*.gbdraw-session.json", "*.gbdraw-session.json.gz")
        for path in SESSION_ROOT.glob(pattern)
    }
    examples = {example.session_path.name for example in EXAMPLES}
    if configured != physical or configured != examples:
        raise ValueError(
            "Gallery session inventory mismatch: "
            f"missing-physical={sorted(configured - physical)}, "
            f"unconfigured-physical={sorted(physical - configured)}, "
            f"missing-examples={sorted(configured - examples)}, "
            f"unconfigured-examples={sorted(examples - configured)}"
        )


def _gallery_mutation_targets(
    session_names: tuple[str, ...],
    *,
    include_assets: bool,
) -> tuple[Path, ...]:
    targets = {_session_path(name) for name in session_names}
    if include_assets:
        from tools.prepare_interactive_gallery_assets import (
            EXAMPLES,
            EXAMPLE_ROOT,
            GALLERY_ROOT as ASSET_GALLERY_ROOT,
            SOURCE_ROOT,
            THUMBNAIL_ROOT,
        )

        for example in EXAMPLES:
            targets.update(
                {
                    example.session_path,
                    example.source_svg_path,
                    example.gallery_svg_path,
                    example.thumbnail_path,
                }
            )
        targets.add(ASSET_GALLERY_ROOT / "examples.json")
        for root, pattern in (
            (EXAMPLE_ROOT, "*.svg"),
            (SOURCE_ROOT, "*.svg"),
            (THUMBNAIL_ROOT, "*.webp"),
        ):
            if root.exists():
                targets.update(root.glob(pattern))
    return tuple(sorted(targets, key=lambda path: str(path)))


@contextmanager
def _gallery_file_transaction(paths: tuple[Path, ...]):
    """Restore every Gallery output if a refresh phase raises an exception."""

    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-backup-") as tmpdir:
        backup_root = Path(tmpdir)
        snapshots: list[tuple[Path, Path | None]] = []
        for index, path in enumerate(paths):
            if path.is_file():
                backup_path = backup_root / f"{index:04d}.bak"
                shutil.copy2(path, backup_path)
                snapshots.append((path, backup_path))
            else:
                snapshots.append((path, None))
        try:
            yield
        except BaseException:
            for path, backup_path in snapshots:
                if backup_path is None:
                    if path.is_file():
                        path.unlink()
                    continue
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(backup_path, path)
            raise


def _session_path(name_or_id: str) -> Path:
    name = name_or_id
    if name.endswith((".gbdraw-session.json", ".gbdraw-session.json.gz")):
        return SESSION_ROOT / name
    compressed_path = SESSION_ROOT / f"{name}.gbdraw-session.json.gz"
    if compressed_path.exists():
        return compressed_path
    return SESSION_ROOT / f"{name}.gbdraw-session.json"


def _session_cli_invocation(session: Mapping[str, Any]) -> Mapping[str, Any] | None:
    cli = session.get("cliInvocation")
    if isinstance(cli, Mapping) and isinstance(cli.get("args"), list):
        return cli
    config = session.get("config")
    if isinstance(config, Mapping):
        cli = config.get("cliInvocation")
        if isinstance(cli, Mapping) and isinstance(cli.get("args"), list):
            return cli
    return None


def _with_interactive_svg_format(args: list[Any]) -> list[str]:
    updated: list[str] = []
    index = 0
    found_format = False
    while index < len(args):
        token = str(args[index])
        if token in {"-f", "--format"}:
            updated.append(token)
            updated.append("interactive_svg")
            index += 2 if index + 1 < len(args) else 1
            found_format = True
            continue
        if token.startswith("--format="):
            updated.append("--format=interactive_svg")
            index += 1
            found_format = True
            continue
        updated.append(token)
        index += 1
    if not found_format:
        updated.extend(["-f", "interactive_svg"])
    return updated


def _preserve_gallery_cli_invocation(
    source_session: Mapping[str, Any],
    refreshed_session: dict[str, Any],
    *,
    mode: str,
) -> bool:
    source_cli = _session_cli_invocation(source_session)
    if source_cli is None:
        return False

    preserved_cli = copy.deepcopy(dict(source_cli))
    preserved_cli["schema"] = 1
    preserved_cli["mode"] = mode
    preserved_cli["args"] = _with_interactive_svg_format(list(source_cli["args"]))
    preserved_cli["renderFormats"] = [INTERACTIVE_SVG_FORMAT]
    preserved_cli.setdefault("fileBindings", [])
    preserved_cli.setdefault("generatedBy", "gbdraw")
    refreshed_session["cliInvocation"] = preserved_cli
    return True


def _promote_gallery_session(
    source_path: Path,
    output_path: Path,
    *,
    env: Mapping[str, str],
) -> None:
    node = shutil.which("node")
    if node is None:
        raise RuntimeError(
            "Gallery session promotion requires Node.js because the browser and "
            "refresh command share the canonical session projection code."
        )
    if not SESSION_PROMOTER.is_file():
        raise FileNotFoundError(
            f"Missing gallery session promoter: {SESSION_PROMOTER.relative_to(REPO_ROOT)}"
        )
    subprocess.run(
        [
            node,
            "--experimental-default-type=module",
            str(SESSION_PROMOTER),
            str(source_path),
            str(output_path),
        ],
        cwd=REPO_ROOT,
        env=dict(env),
        check=True,
    )
    load_session(output_path)


def _merge_refreshed_gallery_artifacts(
    promoted_session: Mapping[str, Any],
    refreshed_session: Mapping[str, Any],
) -> dict[str, Any]:
    """Keep the promoted request authoritative while accepting fresh render artifacts."""

    merged = copy.deepcopy(dict(promoted_session))
    for key in (
        "version",
        "createdAt",
        "results",
        "features",
        "editorState",
        "orthogroupState",
        "losatCache",
        "losatDerivedCache",
        "proteinIdentityManifest",
        "legacyArtifacts",
    ):
        if key in refreshed_session:
            merged[key] = copy.deepcopy(refreshed_session[key])
        elif key in {
            "editorState",
            "losatCache",
            "losatDerivedCache",
            "proteinIdentityManifest",
            "legacyArtifacts",
        }:
            merged.pop(key, None)

    refreshed_resources = refreshed_session.get("resources")
    promoted_resources = promoted_session.get("resources")
    merged["resources"] = {
        **(
            copy.deepcopy(dict(refreshed_resources))
            if isinstance(refreshed_resources, Mapping)
            else {}
        ),
        **(
            copy.deepcopy(dict(promoted_resources))
            if isinstance(promoted_resources, Mapping)
            else {}
        ),
    }
    return merged


def _validate_staged_gallery_session(
    session_path: Path,
    session: dict[str, Any],
) -> None:
    from tools.prepare_interactive_gallery_assets import (
        EXAMPLES,
        _validate_session_interactive_orthogroups,
    )

    validate_session(session)
    if session.get("version") != CURRENT_SESSION_VERSION:
        raise ValueError(
            f"{session_path.name} has session version {session.get('version')}; "
            f"expected {CURRENT_SESSION_VERSION}"
        )
    request = session.get("renderRequest")
    if not isinstance(request, Mapping) or request.get("schema") != 3:
        raise ValueError(f"{session_path.name} has no canonical schema-3 render request")
    resources = session.get("resources")
    if not isinstance(resources, Mapping):
        raise ValueError(f"{session_path.name} has no canonical resources")

    manifest = session.get("proteinIdentityManifest")
    if (
        not isinstance(manifest, Mapping)
        or manifest.get("schema") != PROTEIN_IDENTITY_MANIFEST_SCHEMA
    ):
        raise ValueError(
            f"{session_path.name} has no protein identity manifest schema 1"
        )
    raw_cache = session.get("losatCache")
    raw_entries = raw_cache.get("entries", []) if isinstance(raw_cache, Mapping) else []
    if not isinstance(raw_entries, list) or any(
        classify_raw_losat_cache_entry(entry)
        not in {"protein-current", "nucleotide-current"}
        for entry in raw_entries
    ):
        raise ValueError(f"{session_path.name} contains a non-current raw LOSAT artifact")
    derived_cache = session.get("losatDerivedCache")
    derived_entries = (
        derived_cache.get("entries", [])
        if isinstance(derived_cache, Mapping)
        else []
    )
    if not isinstance(derived_entries, list) or any(
        not isinstance(entry, Mapping) or entry.get("schema") != 2
        for entry in derived_entries
    ):
        raise ValueError(
            f"{session_path.name} contains a non-current derived LOSATP artifact"
        )
    legacy_artifacts = session.get("legacyArtifacts")
    candidates = (
        legacy_artifacts.get("proteinRawCandidates")
        if isinstance(legacy_artifacts, Mapping)
        else None
    )
    if isinstance(candidates, Mapping) and candidates.get("entries"):
        raise ValueError(
            f"{session_path.name} retained legacy protein raw candidates after refresh"
        )
    protein_artifacts = {
        key: session.get(key)
        for key in (
            "results",
            "features",
            "editorState",
            "orthogroupState",
            "losatCache",
            "losatDerivedCache",
            "proteinIdentityManifest",
        )
    }
    if "p_r_" in json.dumps(protein_artifacts, ensure_ascii=False):
        raise ValueError(
            f"{session_path.name} contains unresolved legacy protein identifiers"
        )

    def referenced_resource_ids(value: object):
        if isinstance(value, Mapping):
            resource_id = value.get("resourceId")
            if isinstance(resource_id, str) and resource_id:
                yield resource_id
            for nested in value.values():
                yield from referenced_resource_ids(nested)
        elif isinstance(value, list):
            for nested in value:
                yield from referenced_resource_ids(nested)

    for resource_id in set(referenced_resource_ids(request)):
        resource = resources.get(resource_id)
        if not isinstance(resource, Mapping) or not str(resource.get("data") or ""):
            raise ValueError(
                f"{session_path.name} references missing or empty resource {resource_id}"
            )

    results = session.get("results")
    svg_content = next(
        (
            result.get("content")
            for result in results
            if isinstance(result, Mapping)
            and isinstance(result.get("content"), str)
            and "<svg" in result["content"]
        ),
        None,
    ) if isinstance(results, list) else None
    if not isinstance(svg_content, str):
        raise ValueError(f"{session_path.name} has no generated SVG result")
    ET.fromstring(svg_content)

    example = next(
        (item for item in EXAMPLES if item.session_path.name == session_path.name),
        None,
    )
    if example is None:
        return
    _validate_session_interactive_orthogroups(example, session)


def _refresh_one_session(
    session_path: Path,
    *,
    destination_path: Path | None = None,
) -> None:
    session = load_session(session_path)
    mode = session_mode(session)
    if mode not in {"circular", "linear"}:
        raise RuntimeError(f"Could not determine gallery session mode: {session_path}")
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not env.get("PYTHONPATH")
        else f"{REPO_ROOT}{os.pathsep}{env['PYTHONPATH']}"
    )
    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-session-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        source_session = tmpdir_path / f"source-{session_path.name}"
        write_session_json(source_session, session)
        render_session = tmpdir_path / "promoted.gbdraw-session.json"
        _promote_gallery_session(source_session, render_session, env=env)
        refreshed_session = tmpdir_path / session_path.name
        subprocess.run(
            [
                sys.executable,
                "-m",
                "gbdraw.cli",
                mode,
                "--session",
                str(render_session),
                "-f",
                "interactive_svg",
                "-o",
                "out",
                "--session_output",
                str(refreshed_session),
            ],
            cwd=tmpdir_path,
            env=env,
            check=True,
        )
        promoted_payload = load_session(render_session)
        refreshed_payload = _merge_refreshed_gallery_artifacts(
            promoted_payload,
            load_session(refreshed_session),
        )
        _preserve_gallery_cli_invocation(
            session,
            refreshed_payload,
            mode=mode,
        )
        write_session_json(refreshed_session, refreshed_payload)
        load_session(refreshed_session)
        shutil.move(str(refreshed_session), destination_path or session_path)


def refresh_gallery_sessions(
    session_names: tuple[str, ...] = GALLERY_SESSION_FILES,
) -> None:
    if session_names == GALLERY_SESSION_FILES:
        _validate_gallery_session_inventory()
    session_paths = [_session_path(session_name) for session_name in session_names]
    for session_path in session_paths:
        if not session_path.exists():
            raise FileNotFoundError(
                f"Missing gallery session: {session_path.relative_to(REPO_ROOT)}"
            )

    # Render every requested session successfully before replacing any tracked session.
    # A failed promotion therefore cannot leave a half-refreshed gallery behind.
    with tempfile.TemporaryDirectory(
        prefix=".gbdraw-gallery-refresh-",
        dir=SESSION_ROOT,
    ) as staging_dir:
        staging_root = Path(staging_dir)
        staged_paths: list[tuple[Path, Path]] = []
        for session_path in session_paths:
            print(f"Refreshing gallery session: {session_path.relative_to(REPO_ROOT)}")
            staged_path = staging_root / session_path.name
            _refresh_one_session(session_path, destination_path=staged_path)
            staged_session = load_session(staged_path)
            _validate_staged_gallery_session(session_path, staged_session)
            staged_paths.append((session_path, staged_path))

        for session_path, staged_path in staged_paths:
            staged_path.replace(session_path)


def prepare_gallery_assets() -> None:
    from tools.prepare_interactive_gallery_assets import (
        prepare_gallery_assets as prepare_assets,
    )

    prepare_assets(refresh_sources=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate web gallery session JSON files with the current CLI/session schema. "
            "By default, gallery SVG sources, thumbnails, and examples.json are refreshed afterwards."
        )
    )
    parser.add_argument(
        "--session",
        action="append",
        default=None,
        help=(
            "Refresh one gallery session by file name or gallery id. "
            "May be repeated; defaults to all gallery sessions."
        ),
    )
    parser.add_argument(
        "--no-assets",
        action="store_true",
        help="Only refresh session JSON files; skip gallery SVG/thumbnail/examples.json preparation.",
    )
    parser.add_argument(
        "--skip-session-refresh",
        action="store_true",
        help="Only prepare gallery SVG/thumbnail/examples.json assets from the existing session files.",
    )
    args = parser.parse_args(argv)

    session_names = tuple(args.session or GALLERY_SESSION_FILES)
    targets = _gallery_mutation_targets(
        () if args.skip_session_refresh else session_names,
        include_assets=not args.no_assets,
    )
    with _gallery_file_transaction(targets):
        if not args.skip_session_refresh:
            refresh_gallery_sessions(session_names)
        if not args.no_assets:
            prepare_gallery_assets()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
