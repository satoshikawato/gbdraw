from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
VENDOR_ROOT = WEB_ROOT / "vendor"
PYODIDE_RUNTIME_DIR = VENDOR_ROOT / "pyodide" / "v0.29.0" / "full"

# The Vue build that index.html loads: the global production build, which
# still contains the template compiler that the in-DOM templates need.
VUE_RUNTIME_FILE = Path("vendor") / "vue" / "vue.global.prod.js"

ASSET_URLS = {
    "vue": f"https://unpkg.com/vue@3.5.25/dist/{VUE_RUNTIME_FILE.name}",
    "pyodide_core": "https://github.com/pyodide/pyodide/releases/download/0.29.0/pyodide-core-0.29.0.tar.bz2",
    "micropip": "https://cdn.jsdelivr.net/pyodide/v0.29.0/full/micropip-0.11.0-py3-none-any.whl",
    "tzdata": "https://cdn.jsdelivr.net/pyodide/v0.29.0/full/tzdata-2025.2-py2.py3-none-any.whl",
    "browser_wasi_shim": "https://registry.npmjs.org/@bjorn3/browser_wasi_shim/-/browser_wasi_shim-0.4.2.tgz",
    "jspdf": "https://registry.npmjs.org/jspdf/-/jspdf-3.0.3.tgz",
    "svg2pdf": "https://registry.npmjs.org/svg2pdf.js/-/svg2pdf.js-2.6.0.tgz",
    "dompurify": "https://registry.npmjs.org/dompurify/-/dompurify-3.2.7.tgz",
    "phosphor": "https://registry.npmjs.org/@phosphor-icons/web/-/web-2.1.2.tgz",
    "inter": "https://registry.npmjs.org/@fontsource/inter/-/inter-5.2.8.tgz",
    "noto_sans_jp": "https://registry.npmjs.org/@fontsource/noto-sans-jp/-/noto-sans-jp-5.2.9.tgz",
}

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
        source_url=ASSET_URLS["vue"],
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
        source_url=ASSET_URLS["pyodide_core"],
        bundled_path="vendor/pyodide/v0.29.0/full/",
        notice="The local Pyodide package.json declares MPL-2.0.",
        license_text_anchors=("license-mpl-2.0",),
    ),
    WebAssetNotice(
        display_name="browser_wasi_shim",
        version="0.4.2",
        license_expression="MIT OR Apache-2.0",
        source_url=ASSET_URLS["browser_wasi_shim"],
        bundled_path="vendor/browser_wasi_shim/dist/",
        notice="The package metadata identifies the author as bjorn3.",
        license_text_anchors=("license-mit", "license-apache-2.0"),
    ),
    WebAssetNotice(
        display_name="jsPDF",
        version="3.0.3",
        license_expression="MIT",
        source_url=ASSET_URLS["jspdf"],
        bundled_path="vendor/jspdf/jspdf.umd.min.js",
        notice="Copyright (c) 2010-2025 James Hall and contributors.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="svg2pdf.js",
        version="2.6.0",
        license_expression="MIT",
        source_url=ASSET_URLS["svg2pdf"],
        bundled_path="vendor/svg2pdf.js/svg2pdf.umd.min.js",
        notice="Authored by the yFiles for HTML Support Team at yWorks.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="DOMPurify",
        version="3.2.7",
        license_expression="MPL-2.0 OR Apache-2.0",
        source_url=ASSET_URLS["dompurify"],
        bundled_path="vendor/dompurify/purify.min.js",
        notice="The distributed header identifies Cure53 and other contributors.",
        license_text_anchors=("license-mpl-2.0", "license-apache-2.0"),
    ),
    WebAssetNotice(
        display_name="Phosphor Icons Web",
        version="2.1.2",
        license_expression="MIT",
        source_url=ASSET_URLS["phosphor"],
        bundled_path="vendor/phosphor-icons/regular/",
        notice="Package author: rektdeckard.",
        license_text_anchors=("license-mit",),
    ),
    WebAssetNotice(
        display_name="Inter font package",
        version="Fontsource package 5.2.8",
        license_expression="OFL-1.1",
        source_url=ASSET_URLS["inter"],
        bundled_path="vendor/fonts/inter/",
        notice="Font package metadata identifies Google Inc. as the author.",
        license_text_anchors=("license-ofl-1.1",),
    ),
    WebAssetNotice(
        display_name="Noto Sans JP font package",
        version="Fontsource package 5.2.9",
        license_expression="OFL-1.1",
        source_url=ASSET_URLS["noto_sans_jp"],
        bundled_path="vendor/fonts/noto-sans-jp/",
        notice="Font package metadata identifies Google Inc. as the author.",
        license_text_anchors=("license-ofl-1.1",),
    ),
)
