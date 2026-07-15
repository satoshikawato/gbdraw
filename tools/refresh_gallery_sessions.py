#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw.session_io import load_session, session_mode, write_session_json  # noqa: E402
from gbdraw.render.formats import INTERACTIVE_SVG_FORMAT  # noqa: E402


GALLERY_ROOT = REPO_ROOT / "gbdraw" / "web" / "gallery"
SESSION_ROOT = GALLERY_ROOT / "sessions"

GALLERY_SESSION_FILES = (
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_basic_circular.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "Vnig_TUMSAT-TG-2018.gbdraw-session.json",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json",
    "majanivirus_orthogroup.gbdraw-session.json",
    "lambda_basic_linear.gbdraw-session.json",
)

def _session_path(name_or_id: str) -> Path:
    name = name_or_id
    if not name.endswith(".gbdraw-session.json"):
        name = f"{name}.gbdraw-session.json"
    return SESSION_ROOT / name


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


def _refresh_one_session(session_path: Path) -> None:
    session = load_session(session_path)
    mode = session_mode(session)
    if mode not in {"circular", "linear"}:
        raise RuntimeError(f"Could not determine gallery session mode: {session_path}")
    if _session_cli_invocation(session) is None:
        print(
            "Skipping gallery session without cliInvocation; "
            f"leaving existing output unchanged: {session_path.relative_to(REPO_ROOT)}"
        )
        return

    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(REPO_ROOT)
        if not env.get("PYTHONPATH")
        else f"{REPO_ROOT}{os.pathsep}{env['PYTHONPATH']}"
    )
    with tempfile.TemporaryDirectory(prefix="gbdraw-gallery-session-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        refreshed_session = tmpdir_path / session_path.name
        subprocess.run(
            [
                sys.executable,
                "-m",
                "gbdraw.cli",
                mode,
                "--session",
                str(session_path),
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
        refreshed_payload = load_session(refreshed_session)
        _preserve_gallery_cli_invocation(
            session,
            refreshed_payload,
            mode=mode,
        )
        write_session_json(refreshed_session, refreshed_payload)
        load_session(refreshed_session)
        shutil.move(str(refreshed_session), session_path)


def refresh_gallery_sessions(session_names: tuple[str, ...] = GALLERY_SESSION_FILES) -> None:
    for session_name in session_names:
        session_path = _session_path(session_name)
        if not session_path.exists():
            raise FileNotFoundError(f"Missing gallery session: {session_path.relative_to(REPO_ROOT)}")
        print(f"Refreshing gallery session: {session_path.relative_to(REPO_ROOT)}")
        _refresh_one_session(session_path)


def prepare_gallery_assets() -> None:
    from tools.prepare_interactive_gallery_assets import prepare_gallery_assets as prepare_assets

    prepare_assets()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate web gallery .gbdraw-session.json files with the current CLI/session schema. "
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

    if not args.skip_session_refresh:
        refresh_gallery_sessions(tuple(args.session or GALLERY_SESSION_FILES))
    if not args.no_assets:
        prepare_gallery_assets()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
