from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"
WEB_ROOT = Path(__file__).resolve().parent / "web"
CONFIG_PATH = WEB_ROOT / "js" / "config.js"

BROWSER_WHEEL_BUILD_ENV = "GBDRAW_BUILDING_BROWSER_WHEEL"

_PROJECT_VERSION_RE = re.compile(r'^version\s*=\s*"([^"]+)"\s*$', re.MULTILINE)
_WHEEL_NAME_RE = re.compile(r'^(export const GBDRAW_WHEEL_NAME\s*=\s*")[^"]+(";\s*)$', re.MULTILINE)
_SYNCED = False

_BASE_PACKAGE_DATA = [
    "data/color_palettes.toml",
    "data/config.toml",
    "data/*.ttf",
    "web/index.html",
    "web/open-source-notices.html",
    "web/js/*.js",
    "web/js/app/*.js",
    "web/js/app/*/*.js",
    "web/js/workers/*.js",
    "web/js/services/*.js",
    "web/js/utils/*.js",
    "web/presets/*.txt",
    "web/presets/*.tsv",
    "web/vendor/fonts/*.ttf",
    "web/vendor/*/*.css",
    "web/vendor/*/*.js",
    "web/vendor/*/*.json",
    "web/vendor/*/*.svg",
    "web/vendor/*/*.ttf",
    "web/vendor/*/*.whl",
    "web/vendor/*/*.woff",
    "web/vendor/*/*.woff2",
    "web/vendor/*/*/*.css",
    "web/vendor/*/*/*.js",
    "web/vendor/*/*/*.json",
    "web/vendor/*/*/*.mjs",
    "web/vendor/*/*/*.svg",
    "web/vendor/*/*/*.ttf",
    "web/vendor/*/*/*.whl",
    "web/vendor/*/*/*.woff",
    "web/vendor/*/*/*.woff2",
    "web/vendor/*/*/*.zip",
    "web/vendor/*/*/*.wasm",
    "web/vendor/*/*/*/*.js",
    "web/vendor/*/*/*/*.json",
    "web/vendor/*/*/*/*.mjs",
    "web/vendor/*/*/*/*.whl",
    "web/vendor/*/*/*/*.zip",
    "web/vendor/*/*/*/*.wasm",
    "web/wasm/*/*.md",
    "web/wasm/*/*.wasm",
]


def is_browser_wheel_build() -> bool:
    return os.environ.get(BROWSER_WHEEL_BUILD_ENV) == "1"


def get_package_data_patterns(*, include_browser_wheel: bool) -> list[str]:
    patterns = list(_BASE_PACKAGE_DATA)
    if include_browser_wheel:
        patterns.insert(3, "web/*.whl")
    return patterns


def get_excluded_package_data_patterns(*, include_browser_wheel: bool) -> list[str]:
    if include_browser_wheel:
        return []
    return ["web/*.whl"]


def read_project_version() -> str:
    pyproject_text = PYPROJECT_PATH.read_text(encoding="utf-8")
    match = _PROJECT_VERSION_RE.search(pyproject_text)
    if match is None:
        raise RuntimeError(f"Could not determine project version from {PYPROJECT_PATH}")
    return match.group(1)


def expected_browser_wheel_name(version: str | None = None) -> str:
    if version is None:
        version = read_project_version()
    return f"gbdraw-{version}-py3-none-any.whl"


def remove_browser_wheels_from_build_dir(build_lib: str | Path) -> None:
    build_web_root = Path(build_lib) / "gbdraw" / "web"
    if not build_web_root.exists():
        return
    for wheel_path in build_web_root.glob("gbdraw-*.whl"):
        wheel_path.unlink()


def _update_config_wheel_name(wheel_name: str) -> None:
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    updated_text, replacements = _WHEEL_NAME_RE.subn(
        lambda match: f'{match.group(1)}{wheel_name}{match.group(2)}',
        config_text,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError(f"Could not update GBDRAW_WHEEL_NAME in {CONFIG_PATH}")
    if updated_text != config_text:
        CONFIG_PATH.write_text(updated_text, encoding="utf-8")


def sync_browser_wheel() -> Path | None:
    global _SYNCED
    if _SYNCED or is_browser_wheel_build():
        return None

    version = read_project_version()
    wheel_name = expected_browser_wheel_name(version)

    with tempfile.TemporaryDirectory(prefix="gbdraw-browser-wheel-") as tmpdir:
        tmpdir_path = Path(tmpdir)
        staging_root = tmpdir_path / "source-tree"
        wheelhouse = tmpdir_path / "wheelhouse"
        shutil.copytree(
            REPO_ROOT,
            staging_root,
            ignore=shutil.ignore_patterns(
                ".git",
                ".agents",
                ".codex",
                ".pytest_cache",
                ".ruff_cache",
                "__pycache__",
                "*.egg-info",
                "build",
                "dist",
            ),
        )
        for staged_wheel in (staging_root / "gbdraw" / "web").glob("gbdraw-*.whl"):
            staged_wheel.unlink()

        env = os.environ.copy()
        env[BROWSER_WHEEL_BUILD_ENV] = "1"
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pip",
                "wheel",
                ".",
                "--no-deps",
                "--no-build-isolation",
                "--wheel-dir",
                str(wheelhouse),
            ],
            cwd=staging_root,
            env=env,
            check=True,
        )

        built_wheel_path = wheelhouse / wheel_name
        if not built_wheel_path.exists():
            available = ", ".join(path.name for path in sorted(wheelhouse.glob("gbdraw-*.whl")))
            raise FileNotFoundError(
                f"Expected browser wheel {wheel_name} was not produced. Available wheels: {available or 'none'}"
            )

        for existing_wheel in WEB_ROOT.glob("gbdraw-*.whl"):
            existing_wheel.unlink()

        target_path = WEB_ROOT / wheel_name
        shutil.copy2(built_wheel_path, target_path)

    _update_config_wheel_name(wheel_name)
    _SYNCED = True
    return WEB_ROOT / wheel_name
