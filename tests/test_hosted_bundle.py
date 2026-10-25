"""The hosted bundle's cache contract (tools/prepare_cloudflare_pages.py).

The bundle names the gbdraw wheel by its content hash, so a cached copy of an
older wheel is never served under the new URL.
"""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path

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
