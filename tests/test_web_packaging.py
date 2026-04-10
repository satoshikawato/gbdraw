from __future__ import annotations

import importlib.util
import socket
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
README_PATH = REPO_ROOT / "README.md"
ABOUT_PATH = REPO_ROOT / "docs" / "ABOUT.md"
CITATION_PATH = REPO_ROOT / "CITATION.cff"

PREPRINT_TITLE = "gbdraw: a genome diagram generator for microbes and organelles"
PREPRINT_DOI = "10.64898/2026.04.07.716863"


def _load_verify_module():
    module_path = REPO_ROOT / "tools" / "verify_gui_offline.py"
    spec = importlib.util.spec_from_file_location("verify_gui_offline", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load verification module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _can_bind_loopback() -> bool:
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.bind(("127.0.0.1", 0))
    except OSError:
        return False
    return True


def _run_prepare_browser_wheel(*args: str) -> None:
    subprocess.run(
        [sys.executable, "tools/prepare_browser_wheel.py", *args],
        cwd=REPO_ROOT,
        check=True,
    )


def ensure_prepared_browser_wheel():
    verify_module = _load_verify_module()
    if importlib.util.find_spec("build") is None:
        pytest.skip("python -m build is not available in this environment")

    try:
        browser_wheel_path = verify_module.BUILD_SUPPORT.validate_browser_wheel_prepared()
    except (FileNotFoundError, RuntimeError):
        _run_prepare_browser_wheel()
        browser_wheel_path = verify_module.BUILD_SUPPORT.validate_browser_wheel_prepared()
    return verify_module, browser_wheel_path


def test_web_offline_assets_can_be_prepared_for_packaging() -> None:
    verify_module, expected_wheel_path = ensure_prepared_browser_wheel()
    expected_wheel_name = "gbdraw-0.9.2-py3-none-any.whl"
    assert verify_module._parse_wheel_name() == expected_wheel_name
    assert expected_wheel_path.name == expected_wheel_name
    verify_module.assert_browser_wheel_is_not_recursive(expected_wheel_path)
    verify_module._assert_packaged_assets()


def test_index_links_to_open_source_notices() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "./open-source-notices.html" in index_html
    assert "Open Source Notices" in index_html


def test_index_includes_preprint_citation() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "How to cite" in index_html
    assert PREPRINT_TITLE in index_html
    assert PREPRINT_DOI in index_html


def test_project_docs_and_citation_metadata_include_preprint_doi() -> None:
    assert PREPRINT_DOI in README_PATH.read_text(encoding="utf-8")
    assert PREPRINT_DOI in ABOUT_PATH.read_text(encoding="utf-8")
    citation_cff = CITATION_PATH.read_text(encoding="utf-8")
    assert PREPRINT_DOI in citation_cff
    assert "preferred-citation:" in citation_cff


def test_web_run_analysis_wires_scale_and_tick_font_size_options() -> None:
    source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    assert "tick_label_font_size" in source
    assert '"tick_label_font_size": "--tick_label_font_size" in _source' in source
    assert "args.push('--tick_label_font_size', adv.tick_label_font_size);" in source
    assert "if (form.scale_style === 'ruler')" in source
    assert "args.push('--ruler_label_font_size', adv.scale_font_size);" in source
    assert "args.push('--scale_font_size', adv.scale_font_size);" in source


@pytest.mark.slow
def test_build_py_copies_offline_gui_assets(tmp_path: Path) -> None:
    verify_module, _ = ensure_prepared_browser_wheel()
    build_root = tmp_path / "build_lib"
    subprocess.run(
        [sys.executable, "setup.py", "build_py", "--build-lib", str(build_root)],
        cwd=REPO_ROOT,
        check=True,
    )

    required = [
        build_root / "gbdraw" / "web" / "index.html",
        build_root / "gbdraw" / "web" / "open-source-notices.html",
        build_root / "gbdraw" / "web" / verify_module._parse_wheel_name(),
        build_root / "gbdraw" / "web" / "vendor" / "vue" / "vue.global.js",
        build_root / "gbdraw" / "web" / "vendor" / "tailwindcss" / "tailwindcss-play.js",
        build_root / "gbdraw" / "web" / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
        build_root / "gbdraw" / "web" / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.wasm",
        build_root / "gbdraw" / "web" / "vendor" / "browser_wasi_shim" / "dist" / "index.js",
        build_root / "gbdraw" / "web" / "vendor" / "phosphor-icons" / "regular" / "style.css",
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat.wasm",
        *(build_root / "gbdraw" / "web" / path for path in verify_module.REQUIRED_UI_FONT_FILES),
        *(build_root / "gbdraw" / "web" / path for path in verify_module._parse_local_wheel_paths()),
    ]
    missing = [str(path.relative_to(build_root)) for path in required if not path.exists()]
    assert not missing, "build_py did not copy required offline GUI assets:\n" + "\n".join(missing)
    copied_wheels = sorted(path.name for path in (build_root / "gbdraw" / "web").glob("gbdraw-*.whl"))
    assert copied_wheels == [verify_module._parse_wheel_name()]


@pytest.mark.slow
def test_built_wheel_contains_offline_gui_assets(tmp_path: Path) -> None:
    if importlib.util.find_spec("build") is None:
        pytest.skip("python -m build is not available in this environment")
    if importlib.util.find_spec("wheel") is None:
        pytest.skip("wheel is not available in this environment")

    verify_module, _ = ensure_prepared_browser_wheel()
    dist_dir = tmp_path / "dist"
    subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--no-isolation", "--outdir", str(dist_dir)],
        cwd=REPO_ROOT,
        check=True,
    )

    wheel_path = next(dist_dir.glob("gbdraw-*.whl"))
    assert wheel_path.name == "gbdraw-0.9.2-py3-none-any.whl"
    subprocess.run(
        [sys.executable, "tools/verify_gui_offline.py", "inspect-wheel", str(wheel_path)],
        cwd=REPO_ROOT,
        check=True,
    )
    verify_module.assert_embedded_browser_wheel_is_not_recursive(wheel_path)

    with zipfile.ZipFile(wheel_path) as outer_wheel:
        browser_wheel_member = f"gbdraw/web/{verify_module._parse_wheel_name()}"
        browser_wheels = sorted(
            name
            for name in outer_wheel.namelist()
            if name.startswith("gbdraw/web/gbdraw-") and name.endswith(".whl")
        )
        assert browser_wheels == [browser_wheel_member]


@pytest.mark.slow
def test_offline_gui_smoke_test_covers_palette_preview_behavior() -> None:
    if importlib.util.find_spec("playwright") is None:
        pytest.skip("playwright is not available in this environment")
    if not _can_bind_loopback():
        pytest.skip("loopback sockets are not permitted in this environment")

    ensure_prepared_browser_wheel()
    subprocess.run(
        [sys.executable, "tools/verify_gui_offline.py", "smoke-test"],
        cwd=REPO_ROOT,
        check=True,
    )
