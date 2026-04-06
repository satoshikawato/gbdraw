#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import hashlib
import io
import re
import shutil
import tempfile
import zipfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DIST_ROOT = REPO_ROOT / "dist"
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
CONFIG_PATH = WEB_ROOT / "js" / "config.js"
TOP_LEVEL_BROWSER_WHEEL_RE = re.compile(r"^gbdraw/web/[^/]+\.whl$")
CONFIG_WHEEL_NAME_RE = re.compile(r'export const GBDRAW_WHEEL_NAME = ".*?";')
CONFIG_CACHE_BUST_RE = re.compile(r'export const GBDRAW_WHEEL_CACHE_BUST = ".*?";')
BROWSER_WHEEL_EMBEDDED_CACHE_BUST = "embedded-browser-wheel"


def _project_version() -> str:
    version_ns: dict[str, str] = {}
    exec((REPO_ROOT / "gbdraw" / "version.py").read_text(encoding="utf-8"), version_ns)
    version = version_ns.get("__version__")
    if not version:
        raise RuntimeError("Could not determine gbdraw version from gbdraw/version.py")
    return version


def _canonical_browser_wheel_name() -> str:
    return f"gbdraw-{_project_version()}-py3-none-any.whl"


def _default_outer_wheel_path() -> Path:
    canonical = DIST_ROOT / _canonical_browser_wheel_name()
    if canonical.exists():
        return canonical
    candidates = sorted(DIST_ROOT.glob("gbdraw-*.whl"), key=lambda path: path.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(
            "No built wheel found under dist/. Run `python -m build --wheel` before syncing the browser wheel."
        )
    return candidates[0]


def _top_level_repo_browser_wheels() -> list[Path]:
    return sorted(path for path in WEB_ROOT.glob("*.whl") if path.is_file())


def _read_zip_entries(zip_path: Path) -> list[tuple[str, bytes]]:
    with zipfile.ZipFile(zip_path) as zf:
        return [(info.filename, zf.read(info.filename)) for info in zf.infolist() if not info.is_dir()]


def _urlsafe_sha256(data: bytes) -> str:
    digest = hashlib.sha256(data).digest()
    encoded = base64.urlsafe_b64encode(digest).rstrip(b"=")
    return "sha256=" + encoded.decode("ascii")


def _rewrite_wheel(entries: list[tuple[str, bytes]]) -> bytes:
    filtered = [(name, data) for name, data in entries if not name.endswith(".dist-info/RECORD")]
    record_name = next((name for name, _ in entries if name.endswith(".dist-info/RECORD")), None)
    if record_name is None:
        raise RuntimeError("Wheel does not contain a .dist-info/RECORD entry")

    rows: list[tuple[str, str, str]] = []
    for name, data in filtered:
        rows.append((name, _urlsafe_sha256(data), str(len(data))))
    rows.append((record_name, "", ""))

    record_buffer = io.StringIO()
    writer = csv.writer(record_buffer, lineterminator="\n")
    writer.writerows(rows)
    filtered.append((record_name, record_buffer.getvalue().encode("utf-8")))

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in filtered:
            info = zipfile.ZipInfo(filename=name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            zf.writestr(info, data)
    return buffer.getvalue()


def _build_browser_wheel_bytes(outer_wheel: Path) -> bytes:
    entries = _read_zip_entries(outer_wheel)
    pruned_entries = [
        (
            name,
            _update_config_text(
                data.decode("utf-8"),
                wheel_name=_canonical_browser_wheel_name(),
                cache_bust=BROWSER_WHEEL_EMBEDDED_CACHE_BUST,
            ).encode("utf-8")
            if name == "gbdraw/web/js/config.js"
            else data,
        )
        for name, data in entries
        if not TOP_LEVEL_BROWSER_WHEEL_RE.match(name)
    ]
    return _rewrite_wheel(pruned_entries)


def _cache_bust_value(wheel_bytes: bytes) -> str:
    return hashlib.sha256(wheel_bytes).hexdigest()[:12]


def _update_config_text(config_text: str, *, wheel_name: str, cache_bust: str) -> str:
    if CONFIG_WHEEL_NAME_RE.search(config_text) is None:
        raise RuntimeError("Could not update GBDRAW_WHEEL_NAME in gbdraw/web/js/config.js")
    config_text = CONFIG_WHEEL_NAME_RE.sub(
        f'export const GBDRAW_WHEEL_NAME = "{wheel_name}";',
        config_text,
        count=1,
    )
    if CONFIG_CACHE_BUST_RE.search(config_text) is None:
        raise RuntimeError("Could not update GBDRAW_WHEEL_CACHE_BUST in gbdraw/web/js/config.js")
    return CONFIG_CACHE_BUST_RE.sub(
        f'export const GBDRAW_WHEEL_CACHE_BUST = "{cache_bust}";',
        config_text,
        count=1,
    )


def _write_bytes_atomic(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as tmp:
        tmp.write(data)
        tmp_path = Path(tmp.name)
    tmp_path.replace(path)


def _remove_obsolete_repo_browser_wheels(canonical_name: str) -> list[Path]:
    removed: list[Path] = []
    for path in _top_level_repo_browser_wheels():
        if path.name == canonical_name:
            continue
        path.unlink()
        removed.append(path)
    return removed


def _patch_outer_wheel(
    outer_wheel: Path,
    *,
    browser_wheel_name: str,
    browser_wheel_bytes: bytes,
    config_bytes: bytes,
) -> None:
    entries = _read_zip_entries(outer_wheel)
    patched_entries: list[tuple[str, bytes]] = []
    nested_entry_name = f"gbdraw/web/{browser_wheel_name}"

    for name, data in entries:
        if TOP_LEVEL_BROWSER_WHEEL_RE.match(name) and name != nested_entry_name:
            continue
        if name == nested_entry_name:
            patched_entries.append((name, browser_wheel_bytes))
            continue
        if name == "gbdraw/web/js/config.js":
            patched_entries.append((name, config_bytes))
            continue
        patched_entries.append((name, data))

    if not any(name == nested_entry_name for name, _ in patched_entries):
        patched_entries.append((nested_entry_name, browser_wheel_bytes))

    _write_bytes_atomic(outer_wheel, _rewrite_wheel(patched_entries))


def sync_browser_wheel(outer_wheel: Path) -> dict[str, str]:
    if not outer_wheel.exists():
        raise FileNotFoundError(outer_wheel)
    outer_wheel = outer_wheel.resolve()

    browser_wheel_name = _canonical_browser_wheel_name()
    browser_wheel_bytes = _build_browser_wheel_bytes(outer_wheel)
    browser_wheel_path = WEB_ROOT / browser_wheel_name
    _write_bytes_atomic(browser_wheel_path, browser_wheel_bytes)
    removed_paths = _remove_obsolete_repo_browser_wheels(browser_wheel_name)

    cache_bust = _cache_bust_value(browser_wheel_bytes)
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    updated_config = _update_config_text(
        config_text,
        wheel_name=browser_wheel_name,
        cache_bust=cache_bust,
    )
    CONFIG_PATH.write_text(updated_config, encoding="utf-8")

    _patch_outer_wheel(
        outer_wheel,
        browser_wheel_name=browser_wheel_name,
        browser_wheel_bytes=browser_wheel_bytes,
        config_bytes=updated_config.encode("utf-8"),
    )

    return {
        "outer_wheel": str(outer_wheel.relative_to(REPO_ROOT)),
        "browser_wheel": str(browser_wheel_path.relative_to(REPO_ROOT)),
        "cache_bust": cache_bust,
        "removed": ", ".join(str(path.relative_to(REPO_ROOT)) for path in removed_paths) or "(none)",
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate the browser wheel referenced by gbdraw/web/js/config.js from a built outer wheel, "
            "remove obsolete top-level browser wheels, update the cache-bust token, and patch the outer wheel "
            "to embed the same browser wheel for gbdraw gui."
        )
    )
    parser.add_argument(
        "outer_wheel",
        nargs="?",
        type=Path,
        help="Optional path to the built outer wheel. Defaults to dist/gbdraw-<version>-py3-none-any.whl.",
    )
    args = parser.parse_args()

    try:
        outer_wheel = args.outer_wheel or _default_outer_wheel_path()
        result = sync_browser_wheel(outer_wheel)
    except (FileNotFoundError, RuntimeError) as exc:
        print(str(exc))
        return 1

    print(f"Synchronized browser wheel: {result['browser_wheel']}")
    print(f"Patched outer wheel: {result['outer_wheel']}")
    print(f"Updated cache bust: {result['cache_bust']}")
    print(f"Removed obsolete wheels: {result['removed']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
