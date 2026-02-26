from __future__ import annotations

import re
from importlib import resources
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INDEX_HTML = PROJECT_ROOT / "gbdraw" / "web" / "index.html"
LOSAT_JS = PROJECT_ROOT / "gbdraw" / "web" / "js" / "services" / "losat.js"
PYODIDE_JS = PROJECT_ROOT / "gbdraw" / "web" / "js" / "app" / "pyodide.js"
WEB_CONFIG_JS = PROJECT_ROOT / "gbdraw" / "web" / "js" / "config.js"
TAILWIND_SOURCE_CSS = PROJECT_ROOT / "gbdraw" / "web" / "css" / "tailwind.source.css"
TAILWIND_OUTPUT_CSS = PROJECT_ROOT / "gbdraw" / "web" / "vendor" / "tailwind" / "app.min.css"
VENDOR_DIR = PROJECT_ROOT / "gbdraw" / "web" / "vendor"
PYODIDE_WHEELS_DIR = VENDOR_DIR / "pyodide-wheels"
PYODIDE_CORE_DIR = VENDOR_DIR / "pyodide" / "v0.29.0" / "full"
FONTS_DIR = VENDOR_DIR / "fonts"

RUNTIME_VENDOR_FILES = [
    VENDOR_DIR / "vue" / "vue.global.js",
    VENDOR_DIR / "phosphor" / "index.js",
    VENDOR_DIR / "jspdf" / "jspdf.umd.min.js",
    VENDOR_DIR / "svg2pdf" / "svg2pdf.umd.min.js",
    VENDOR_DIR / "dompurify" / "purify.min.js",
    VENDOR_DIR / "browser_wasi_shim" / "index.js",
    VENDOR_DIR / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
]

PYODIDE_CORE_REQUIRED_FILES = [
    "ffi.d.ts",
    "package.json",
    "pyodide-lock.json",
    "pyodide.asm.js",
    "pyodide.asm.wasm",
    "pyodide.d.ts",
    "pyodide.js",
    "pyodide.mjs",
    "python",
    "python_cli_entry.mjs",
    "python_stdlib.zip",
    "micropip-0.11.0-py3-none-any.whl",
]

BROWSER_WASI_SHIM_REQUIRED_FILES = [
    "index.js",
    "fd.js",
    "fs_mem.js",
    "fs_opfs.js",
    "wasi.js",
    "wasi_defs.js",
    "strace.js",
    "debug.js",
]

VENDORED_FONT_FILES = [
    "Inter-Variable.ttf",
    "NotoSansJP-Variable.ttf",
    "OFL-Inter.txt",
    "OFL-NotoSansJP.txt",
]


def _parse_required_wheels(config_content: str) -> list[str]:
    match = re.search(
        r"PYODIDE_REQUIRED_WHEELS\s*=\s*\[(?P<body>.*?)\]",
        config_content,
        flags=re.S,
    )
    assert match is not None
    return re.findall(r'"([^"]+\.whl)"', match.group("body"))


def test_index_html_has_no_cdn_runtime_dependencies() -> None:
    content = INDEX_HTML.read_text(encoding="utf-8")
    forbidden_markers = [
        "https://cdn",
        "https://unpkg",
        "https://cdnjs",
        "https://fonts.googleapis",
    ]
    for marker in forbidden_markers:
        assert marker not in content


def test_index_html_uses_static_tailwind_css() -> None:
    content = INDEX_HTML.read_text(encoding="utf-8")
    assert "./vendor/tailwind/app.min.css" in content
    assert "tailwindcss.cdn.js" not in content
    assert "@apply" not in content


def test_tailwind_source_and_generated_css_contract() -> None:
    source = TAILWIND_SOURCE_CSS.read_text(encoding="utf-8")
    generated = TAILWIND_OUTPUT_CSS.read_text(encoding="utf-8")

    assert "@apply" in source
    assert "@apply" not in generated
    assert ".fixed{position:fixed}" in generated
    assert ".flex{display:flex}" in generated


def test_vendored_fonts_exist_and_tailwind_references_them() -> None:
    source = TAILWIND_SOURCE_CSS.read_text(encoding="utf-8")
    generated = TAILWIND_OUTPUT_CSS.read_text(encoding="utf-8")

    for filename in VENDORED_FONT_FILES:
        assert (FONTS_DIR / filename).is_file()

    assert "font-family: \"Inter\", \"Noto Sans JP\", sans-serif;" in source
    assert "../fonts/Inter-Variable.ttf" in generated
    assert "../fonts/NotoSansJP-Variable.ttf" in generated


def test_vendor_runtime_assets_are_real_files_not_placeholders() -> None:
    for asset in RUNTIME_VENDOR_FILES:
        content = asset.read_text(encoding="utf-8")
        assert "Placeholder vendor asset" not in content

    phosphor_loader = (VENDOR_DIR / "phosphor" / "index.js").read_text(encoding="utf-8")
    assert "https://cdn.jsdelivr.net" not in phosphor_loader
    assert "./vendor/phosphor/regular/style.css" in phosphor_loader


def test_losat_does_not_import_wasi_shim_from_cdn() -> None:
    content = LOSAT_JS.read_text(encoding="utf-8")
    assert "https://unpkg.com/@bjorn3/browser_wasi_shim" not in content
    assert "../../vendor/browser_wasi_shim/index.js" in content


def test_pyodide_uses_local_index_url_and_explicit_wheel_list() -> None:
    config_content = WEB_CONFIG_JS.read_text(encoding="utf-8")
    pyodide_content = PYODIDE_JS.read_text(encoding="utf-8")
    assert "PYODIDE_INDEX_URL" in config_content
    assert "PYODIDE_WHEELS_BASE_URL" in config_content
    assert "PYODIDE_REQUIRED_WHEELS" in config_content
    assert "loadPyodide({ indexURL:" in pyodide_content
    assert "ensureDependencyWheelsAvailable" in pyodide_content
    assert "micropip.install(dependencyWheelUrls)" in pyodide_content


def test_required_pyodide_wheels_are_present_and_valid_archives() -> None:
    config_content = WEB_CONFIG_JS.read_text(encoding="utf-8")
    required_wheels = _parse_required_wheels(config_content)
    assert required_wheels

    for wheel_name in required_wheels:
        wheel_path = PYODIDE_WHEELS_DIR / wheel_name
        assert wheel_path.is_file()
        header = wheel_path.read_bytes()[:2]
        assert header == b"PK"


def test_pyodide_core_required_files_exist() -> None:
    for filename in PYODIDE_CORE_REQUIRED_FILES:
        assert (PYODIDE_CORE_DIR / filename).is_file()


def test_browser_wasi_shim_required_modules_exist() -> None:
    shim_dir = VENDOR_DIR / "browser_wasi_shim"
    for filename in BROWSER_WASI_SHIM_REQUIRED_FILES:
        assert (shim_dir / filename).is_file()


def test_web_vendor_manifest_is_packaged_resource() -> None:
    vendor_manifest = resources.files("gbdraw").joinpath("web", "vendor", "manifest.json")
    assert vendor_manifest.is_file()


def test_losat_wasm_is_packaged_resource() -> None:
    losat_wasm = resources.files("gbdraw").joinpath("web", "wasm", "losat", "losat.wasm")
    assert losat_wasm.is_file()


def test_pyodide_core_assets_are_packaged_resources() -> None:
    core_root = resources.files("gbdraw").joinpath("web", "vendor", "pyodide", "v0.29.0", "full")
    for filename in PYODIDE_CORE_REQUIRED_FILES:
        assert core_root.joinpath(filename).is_file()


def test_required_pyodide_wheels_are_packaged_resources() -> None:
    config_content = WEB_CONFIG_JS.read_text(encoding="utf-8")
    required_wheels = _parse_required_wheels(config_content)
    wheels_root = resources.files("gbdraw").joinpath("web", "vendor", "pyodide-wheels")
    for wheel_name in required_wheels:
        assert wheels_root.joinpath(wheel_name).is_file()


def test_vendored_fonts_are_packaged_resources() -> None:
    fonts_root = resources.files("gbdraw").joinpath("web", "vendor", "fonts")
    for filename in VENDORED_FONT_FILES:
        assert fonts_root.joinpath(filename).is_file()
