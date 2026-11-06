"""The vendored Web assets: their files, source URLs, versions and notices."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
VENDOR_ROOT = WEB_ROOT / "vendor"
PYODIDE_RUNTIME_DIR = VENDOR_ROOT / "pyodide" / "v0.29.0" / "full"

# The Vue build that index.html loads: the global production build, which
# still contains the template compiler that the in-DOM templates need.
VUE_RUNTIME_FILE = Path("vendor") / "vue" / "vue.global.prod.js"


@dataclass(frozen=True)
class VendoredDownload:
    """A file that `tools/verify_gui_offline.py vendor-assets` downloads, and the SHA-256 its bytes must have."""

    url: str
    sha256: str


# Each SHA-256 is that of the download whose files are the vendored ones; an npm
# tarball's bytes also match the integrity value the npm registry publishes.
ASSET_DOWNLOADS = {
    "vue": VendoredDownload(
        f"https://unpkg.com/vue@3.5.25/dist/{VUE_RUNTIME_FILE.name}",
        "78133fb37d7da6ab9fb0ede6c8ed446f38773595d1faf52a3c1b77d539b06b61",
    ),
    "pyodide_core": VendoredDownload(
        "https://github.com/pyodide/pyodide/releases/download/0.29.0/pyodide-core-0.29.0.tar.bz2",
        "fe85edd0d5c9828598da14c98d2f79cea45aac5a9e6448ca65dcc2d58d4ed033",
    ),
    "micropip": VendoredDownload(
        "https://cdn.jsdelivr.net/pyodide/v0.29.0/full/micropip-0.11.0-py3-none-any.whl",
        "d98e6100df5e4145c1b11d559a8620e92b98d994cc1b52381a8815cebf3e7128",
    ),
    "tzdata": VendoredDownload(
        "https://cdn.jsdelivr.net/pyodide/v0.29.0/full/tzdata-2025.2-py2.py3-none-any.whl",
        "efe840e3b9cc04707b6fbca5c145910e50916f801a61a8db6730ba4c92abba86",
    ),
    "browser_wasi_shim": VendoredDownload(
        "https://registry.npmjs.org/@bjorn3/browser_wasi_shim/-/browser_wasi_shim-0.4.2.tgz",
        "9c0281520d0e99f027ec7c1c79b4036c0f8168ed9bf98aba19db4737a1333782",
    ),
    "jspdf": VendoredDownload(
        "https://registry.npmjs.org/jspdf/-/jspdf-3.0.3.tgz",
        "e3615356fb5d2408e05215a848dec3924f67fb83f0a44956bdbe55c5398f36f3",
    ),
    "svg2pdf": VendoredDownload(
        "https://registry.npmjs.org/svg2pdf.js/-/svg2pdf.js-2.6.0.tgz",
        "c242ebf7a9b8cbc64a7cc4cd6058026485f11ba95e55c3a5b57fc23a3acf0edd",
    ),
    "dompurify": VendoredDownload(
        "https://registry.npmjs.org/dompurify/-/dompurify-3.2.7.tgz",
        "e6e0410255ad080c4efc43d797348f372b0c2d94218321c8ed0b83ba57f42e43",
    ),
    "phosphor": VendoredDownload(
        "https://registry.npmjs.org/@phosphor-icons/web/-/web-2.1.2.tgz",
        "40e3096099ca818c047979cece6a3764944d5200c67b64092f7bcd8bcd1d2b08",
    ),
    "inter": VendoredDownload(
        "https://registry.npmjs.org/@fontsource/inter/-/inter-5.2.8.tgz",
        "7004eee208371e46125e999c67b094ba47f3606b1983331c6f8e77b39a2f825c",
    ),
    "noto_sans_jp": VendoredDownload(
        "https://registry.npmjs.org/@fontsource/noto-sans-jp/-/noto-sans-jp-5.2.9.tgz",
        "e4d4ba8825aebe6293e19400a07d45b8f111ff31532c8dfbdd77b986c0123652",
    ),
}

# The files vendor-assets takes from a package folder; the package holds others.
BROWSER_WASI_SHIM_FILES = (
    "debug.js",
    "fd.js",
    "fs_mem.js",
    "fs_opfs.js",
    "index.js",
    "strace.js",
    "wasi.js",
    "wasi_defs.js",
)
PHOSPHOR_REGULAR_FILES = ("style.css", "Phosphor.svg", "Phosphor.ttf", "Phosphor.woff", "Phosphor.woff2")

_NPM_DOWNLOAD_URLS = (
    re.compile(r"https://registry\.npmjs\.org/(?P<name>(?:@[^/]+/)?[^/]+)/-/[^/]+-(?P<version>[0-9][^/]*)\.tgz"),
    re.compile(r"https://unpkg\.com/(?P<name>(?:@[^/]+/)?[^@/]+)@(?P<version>[0-9][^/]*)/.+"),
)
_NPM_HOSTS = ("https://registry.npmjs.org/", "https://unpkg.com/")


def npm_package(url: str) -> tuple[str, str] | None:
    """The npm package name and version a download URL names; None for a download from elsewhere."""

    for pattern in _NPM_DOWNLOAD_URLS:
        match = pattern.fullmatch(url)
        if match:
            return match["name"], match["version"]
    if url.startswith(_NPM_HOSTS):
        raise ValueError(f"{url} names no npm package version")
    return None


# No tool installs this manifest. GitHub's dependency graph reads it, so
# Dependabot reports an advisory for a vendored npm package.
NPM_ADVISORY_MANIFEST = Path("tools") / "vendored-web-assets" / "package.json"


def npm_advisory_manifest_text() -> str:
    """The text of NPM_ADVISORY_MANIFEST: each vendored npm package at the version vendored."""

    packages = sorted(
        package for package in (npm_package(download.url) for download in ASSET_DOWNLOADS.values()) if package
    )
    manifest = {
        "name": "gbdraw-vendored-web-assets",
        "private": True,
        "description": (
            "The npm packages vendored under gbdraw/web/vendor/, for advisories only; nothing installs them. "
            "Written by `python tools/verify_gui_offline.py vendor-assets` from gbdraw/_web_assets.py."
        ),
        "dependencies": dict(packages),
    }
    return json.dumps(manifest, indent=2) + "\n"

# The files of PYODIDE_RUNTIME_DIR that Pyodide's loader downloads at start,
# in that order (config.js PYODIDE_RUNTIME_FILES, the page's prefetch).
PYODIDE_LOADER_FILES = (
    "pyodide.asm.js",
    "pyodide.asm.wasm",
    "python_stdlib.zip",
    "pyodide-lock.json",
)

PYODIDE_CORE_FILES = {
    "package.json",
    "pyodide.js",
    "pyodide.mjs",
    *PYODIDE_LOADER_FILES,
}

PYODIDE_RUNTIME_PACKAGE_WHEELS = {
    "micropip": "micropip-0.11.0-py3-none-any.whl",
    "tzdata": "tzdata-2025.2-py2.py3-none-any.whl",
}

PYODIDE_LOCAL_WHEELS = (
    Path("vendor") / "pyodide-wheels" / "six-1.17.0-py2.py3-none-any.whl",
    Path("vendor") / "pyodide-wheels" / "python_dateutil-2.9.0.post0-py2.py3-none-any.whl",
    Path("vendor") / "pyodide-wheels" / "pytz-2025.2-py2.py3-none-any.whl",
    Path("vendor") / "pyodide-wheels" / "numpy-2.2.5-cp313-cp313-pyodide_2025_0_wasm32.whl",
    Path("vendor") / "pyodide-wheels" / "biopython-1.85-cp313-cp313-pyodide_2025_0_wasm32.whl",
    Path("vendor") / "pyodide-wheels" / "fonttools-4.56.0-py3-none-any.whl",
    Path("vendor") / "pyodide-wheels" / "svgwrite-1.4.3-py3-none-any.whl",
    Path("vendor") / "pyodide-wheels" / "bcbio_gff-0.7.1-py3-none-any.whl",
    Path("vendor") / "pyodide-wheels" / "pandas-2.3.2-cp313-cp313-pyodide_2025_0_wasm32.whl",
)

PYODIDE_RUNTIME_WHEELS = tuple(
    Path("vendor") / "pyodide" / "v0.29.0" / "full" / filename
    for filename in PYODIDE_RUNTIME_PACKAGE_WHEELS.values()
)

PYODIDE_BUNDLED_WHEELS = PYODIDE_LOCAL_WHEELS + PYODIDE_RUNTIME_WHEELS

UI_FONT_ASSETS = {
    "inter": (
        "inter-latin-300-normal.woff2",
        "inter-latin-400-normal.woff2",
        "inter-latin-500-normal.woff2",
        "inter-latin-600-normal.woff2",
        "inter-latin-700-normal.woff2",
        "inter-latin-ext-300-normal.woff2",
        "inter-latin-ext-400-normal.woff2",
        "inter-latin-ext-500-normal.woff2",
        "inter-latin-ext-600-normal.woff2",
        "inter-latin-ext-700-normal.woff2",
        "inter-cyrillic-300-normal.woff2",
        "inter-cyrillic-400-normal.woff2",
        "inter-cyrillic-500-normal.woff2",
        "inter-cyrillic-600-normal.woff2",
        "inter-cyrillic-700-normal.woff2",
        "inter-vietnamese-300-normal.woff2",
        "inter-vietnamese-400-normal.woff2",
        "inter-vietnamese-500-normal.woff2",
        "inter-vietnamese-600-normal.woff2",
        "inter-vietnamese-700-normal.woff2",
    ),
    "noto-sans-jp": (
        "noto-sans-jp-japanese-300-normal.woff2",
        "noto-sans-jp-japanese-400-normal.woff2",
        "noto-sans-jp-japanese-500-normal.woff2",
        "noto-sans-jp-japanese-700-normal.woff2",
    ),
}

REQUIRED_UI_FONT_FILES = tuple(
    Path("vendor") / "fonts" / family / filename
    for family, filenames in UI_FONT_ASSETS.items()
    for filename in filenames
)


@dataclass(frozen=True)
class WebAssetNotice:
    display_name: str
    version: str
    license_expression: str
    source_url: str
    bundled_path: str
    notice: str
    license_text_anchors: tuple[str, ...] = ()


WEB_ASSET_NOTICES = (
    WebAssetNotice(
        display_name="Vue.js",
        version="3.5.25",
        license_expression="MIT",
        source_url=ASSET_DOWNLOADS["vue"].url,
        bundled_path=VUE_RUNTIME_FILE.as_posix(),
        notice="Copyright (c) 2018-present Yuxi (Evan) You and Vue contributors.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="Tailwind CSS (generated stylesheet)",
        version="3.4.17",
        license_expression="MIT",
        source_url="https://tailwindcss.com",
        bundled_path="css/tailwind.css",
        notice=(
            "Copyright (c) Tailwind Labs, Inc. The stylesheet is generated with the "
            "Tailwind CSS compiler and contains its base styles (Preflight)."
        ),
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="Pyodide core runtime",
        version="0.29.0",
        license_expression="MPL-2.0",
        source_url=ASSET_DOWNLOADS["pyodide_core"].url,
        bundled_path="vendor/pyodide/v0.29.0/full/",
        notice="The local Pyodide package.json declares MPL-2.0.",
        license_text_anchors=("license-mpl-2.0",),
    ),
    WebAssetNotice(
        display_name="browser_wasi_shim",
        version="0.4.2",
        license_expression="MIT OR Apache-2.0",
        source_url=ASSET_DOWNLOADS["browser_wasi_shim"].url,
        bundled_path="vendor/browser_wasi_shim/dist/",
        notice="The package metadata identifies the author as bjorn3.",
        license_text_anchors=("license-mit", "license-apache-2.0"),
    ),
    WebAssetNotice(
        display_name="jsPDF",
        version="3.0.3",
        license_expression="MIT",
        source_url=ASSET_DOWNLOADS["jspdf"].url,
        bundled_path="vendor/jspdf/jspdf.umd.min.js",
        notice="Copyright (c) 2010-2025 James Hall and contributors.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="svg2pdf.js",
        version="2.6.0",
        license_expression="MIT",
        source_url=ASSET_DOWNLOADS["svg2pdf"].url,
        bundled_path="vendor/svg2pdf.js/svg2pdf.umd.min.js",
        notice="Authored by the yFiles for HTML Support Team at yWorks.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="DOMPurify",
        version="3.2.7",
        license_expression="MPL-2.0 OR Apache-2.0",
        source_url=ASSET_DOWNLOADS["dompurify"].url,
        bundled_path="vendor/dompurify/purify.min.js",
        notice="The distributed header identifies Cure53 and other contributors.",
        license_text_anchors=("license-mpl-2.0", "license-apache-2.0"),
    ),
    WebAssetNotice(
        display_name="Phosphor Icons Web",
        version="2.1.2",
        license_expression="MIT",
        source_url=ASSET_DOWNLOADS["phosphor"].url,
        bundled_path="vendor/phosphor-icons/regular/",
        notice="Package author: rektdeckard.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="Inter font package",
        version="Fontsource package 5.2.8",
        license_expression="OFL-1.1",
        source_url=ASSET_DOWNLOADS["inter"].url,
        bundled_path="vendor/fonts/inter/",
        notice="Font package metadata identifies Google Inc. as the author.",
        license_text_anchors=("license-ofl-1.1",),
    ),
    WebAssetNotice(
        display_name="Noto Sans JP font package",
        version="Fontsource package 5.2.9",
        license_expression="OFL-1.1",
        source_url=ASSET_DOWNLOADS["noto_sans_jp"].url,
        bundled_path="vendor/fonts/noto-sans-jp/",
        notice="Font package metadata identifies Google Inc. as the author.",
        license_text_anchors=("license-ofl-1.1",),
    ),
)
