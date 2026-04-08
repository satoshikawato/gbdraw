from __future__ import annotations

import os
import re
from datetime import datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"
WEB_ROOT = Path(__file__).resolve().parent / "web"
CONFIG_PATH = WEB_ROOT / "js" / "config.js"

BROWSER_WHEEL_BUILD_ENV = "GBDRAW_BUILDING_BROWSER_WHEEL"

_PROJECT_VERSION_RE = re.compile(r'^version\s*=\s*"([^"]+)"\s*$', re.MULTILINE)
_WHEEL_NAME_RE = re.compile(r'^(export const GBDRAW_WHEEL_NAME\s*=\s*")[^"]+(";\s*)$', re.MULTILINE)
_WHEEL_CACHE_BUST_RE = re.compile(r'^(export const GBDRAW_WHEEL_CACHE_BUST\s*=\s*")[^"]+(";\s*)$', re.MULTILINE)

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


def read_browser_wheel_name_from_config() -> str:
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    match = _WHEEL_NAME_RE.search(config_text)
    if match is None:
        raise RuntimeError(f"Could not determine GBDRAW_WHEEL_NAME from {CONFIG_PATH}")
    _, value = match.group(0).split("=", 1)
    return value.strip().strip(";").strip().strip('"').strip("'")


def generate_cache_bust_token() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def update_browser_wheel_config(*, wheel_name: str, cache_bust: str | None = None) -> None:
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    updated_text, replacements = _WHEEL_NAME_RE.subn(
        lambda match: f'{match.group(1)}{wheel_name}{match.group(2)}',
        config_text,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError(f"Could not update GBDRAW_WHEEL_NAME in {CONFIG_PATH}")
    if cache_bust is not None:
        updated_text, cache_bust_replacements = _WHEEL_CACHE_BUST_RE.subn(
            lambda match: f'{match.group(1)}{cache_bust}{match.group(2)}',
            updated_text,
            count=1,
        )
        if cache_bust_replacements != 1:
            raise RuntimeError(f"Could not update GBDRAW_WHEEL_CACHE_BUST in {CONFIG_PATH}")
    if updated_text != config_text:
        CONFIG_PATH.write_text(updated_text, encoding="utf-8")


def validate_browser_wheel_prepared() -> Path:
    expected_name = expected_browser_wheel_name()
    expected_path = WEB_ROOT / expected_name
    existing_wheels = sorted(WEB_ROOT.glob("gbdraw-*.whl"))

    if not expected_path.exists():
        raise FileNotFoundError(
            f"Missing browser wheel {expected_path.relative_to(REPO_ROOT)}. "
            "Run `python tools/prepare_browser_wheel.py`."
        )

    unexpected_wheels = [path.name for path in existing_wheels if path.name != expected_name]
    if unexpected_wheels:
        raise RuntimeError(
            "Unexpected browser wheel files in gbdraw/web: "
            f"{', '.join(unexpected_wheels)}. Run `python tools/prepare_browser_wheel.py` to refresh them."
        )

    configured_name = read_browser_wheel_name_from_config()
    if configured_name != expected_name:
        raise RuntimeError(
            f"gbdraw/web/js/config.js points to {configured_name}, expected {expected_name}. "
            "Run `python tools/prepare_browser_wheel.py`."
        )

    return expected_path
