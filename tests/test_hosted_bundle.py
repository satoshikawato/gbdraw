"""The hosted bundle's cache contract (tools/prepare_cloudflare_pages.py).

The bundle names the gbdraw wheel by its content hash, so a cached copy of an
older wheel is never served under the new URL.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

import pytest

from gbdraw import _build_support as build_support

REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"


def test_source_config_carries_no_wheel_cache_bust() -> None:
    # Only the hosted bundle carries a token; the package (`gbdraw gui`) and the
    # CI wheel jobs load the wheel by its plain name.
    config_text = (WEB_ROOT / "js" / "config.js").read_text(encoding="utf-8")
    assert "export const GBDRAW_WHEEL_CACHE_BUST = null;" in config_text


def test_wheel_cache_bust_is_the_wheel_content_hash(tmp_path: Path) -> None:
    wheel = tmp_path / "gbdraw-0.0.0-py3-none-any.whl"
    wheel.write_bytes(b"wheel bytes")
    token = build_support.browser_wheel_cache_bust_token(wheel)
    assert token == hashlib.sha256(b"wheel bytes").hexdigest()[:16]

    config_path = tmp_path / "config.js"
    shutil.copyfile(WEB_ROOT / "js" / "config.js", config_path)
    build_support.write_wheel_cache_bust(config_path, token)
    config_text = config_path.read_text(encoding="utf-8")
    assert f'export const GBDRAW_WHEEL_CACHE_BUST = "{token}";' in config_text
    assert config_text.count("GBDRAW_WHEEL_CACHE_BUST =") == 1


def _load_cloudflare_module():
    from importlib.util import module_from_spec, spec_from_file_location

    spec = spec_from_file_location(
        "prepare_cloudflare_pages_for_cache_tests",
        REPO_ROOT / "tools" / "prepare_cloudflare_pages.py",
    )
    assert spec is not None and spec.loader is not None
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bundle_headers_cache_only_the_immutable_vendor_paths() -> None:
    cloudflare = _load_cloudflare_module()
    headers = cloudflare.render_bundle_headers()
    immutable = "Cache-Control: public, max-age=31536000, immutable"
    blocks = {
        block.splitlines()[0]: block for block in headers.strip().split("\n\n")
    }
    for path in (
        "/vendor/pyodide/*",
        "/vendor/pyodide-wheels/*",
        "/vendor/fonts/*",
        "/vendor/phosphor-icons/*",
    ):
        assert immutable in blocks[path]
    # Mutable URLs keep the platform default (revalidate every load).
    assert headers.count(immutable) == 4
    assert "Cross-Origin-Embedder-Policy: require-corp" in blocks["/*"]
    assert "frame-ancestors 'self'" in blocks["/gallery/examples/*"]


def test_immutable_vendor_files_match_their_manifest() -> None:
    cloudflare = _load_cloudflare_module()
    cloudflare.verify_immutable_assets(WEB_ROOT)


def test_immutable_vendor_file_changed_in_place_fails(tmp_path: Path) -> None:
    cloudflare = _load_cloudflare_module()
    root = tmp_path / "web"
    shutil.copytree(WEB_ROOT / "vendor" / "fonts", root / "vendor" / "fonts")
    shutil.copytree(WEB_ROOT / "vendor" / "phosphor-icons", root / "vendor" / "phosphor-icons")
    shutil.copytree(WEB_ROOT / "vendor" / "pyodide", root / "vendor" / "pyodide")
    shutil.copytree(WEB_ROOT / "vendor" / "pyodide-wheels", root / "vendor" / "pyodide-wheels")
    cloudflare.verify_immutable_assets(root)

    changed = root / "vendor" / "pyodide-wheels" / "six-1.17.0-py2.py3-none-any.whl"
    changed.write_bytes(changed.read_bytes() + b"\0")
    added = root / "vendor" / "fonts" / "inter" / "inter-latin-800-normal.woff2"
    added.write_bytes(b"font")
    with pytest.raises(RuntimeError) as failure:
        cloudflare.verify_immutable_assets(root)
    message = str(failure.value)
    assert "vendor/pyodide-wheels/six-1.17.0-py2.py3-none-any.whl" in message
    assert "vendor/fonts/inter/inter-latin-800-normal.woff2" in message
