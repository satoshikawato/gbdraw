#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "dist" / "cloudflare-pages"
DEFAULT_ANALYTICS_TOKEN = "e4dc2e66d09549868f5a5ac7d7a6e633"
ANALYTICS_TOKEN_ENV = "CLOUDFLARE_WEB_ANALYTICS_TOKEN"
SCRIPT_MARKER = "<!-- CLOUDFLARE_WEB_ANALYTICS_SCRIPT -->"
NOTICE_MARKER = "<!-- CLOUDFLARE_WEB_ANALYTICS_NOTICE -->"


def _load_prepare_browser_wheel_module():
    module_path = REPO_ROOT / "tools" / "prepare_browser_wheel.py"
    spec = spec_from_file_location("prepare_browser_wheel", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load browser wheel preparation module from {module_path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _replace_once(source: str, old: str, new: str) -> str:
    if old not in source:
        raise RuntimeError(f"Expected text not found while preparing Cloudflare bundle: {old}")
    return source.replace(old, new, 1)


def _render_analytics_script(token: str) -> str:
    return (
        "    <!-- Cloudflare Web Analytics --><script defer "
        "src=\"https://static.cloudflareinsights.com/beacon.min.js\" "
        f"data-cf-beacon='{{\"token\": \"{token}\"}}'></script><!-- End Cloudflare Web Analytics -->"
    )


def _render_analytics_notice() -> str:
    return """                            <div class="bg-amber-50 text-amber-900 p-2.5 rounded-lg border border-amber-200 leading-relaxed">
                                <strong class="flex items-center gap-1 mb-1"><i class="ph ph-chart-line"></i> Hosted Site Analytics</strong>
                                The hosted <strong>gbdraw.app</strong> deployment uses Cloudflare Web Analytics for aggregate page-visit and performance metrics. Uploaded genome files are still processed locally in your browser and are not sent to Cloudflare by gbdraw itself.
                            </div>"""


def build_cloudflare_pages_bundle(
    *,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    analytics_token: str = DEFAULT_ANALYTICS_TOKEN,
) -> Path:
    if output_root.exists():
        shutil.rmtree(output_root)
    shutil.copytree(WEB_ROOT, output_root)

    index_path = output_root / "index.html"
    index_html = index_path.read_text(encoding="utf-8")
    index_html = _replace_once(index_html, SCRIPT_MARKER, _render_analytics_script(analytics_token))
    index_html = _replace_once(index_html, NOTICE_MARKER, _render_analytics_notice())
    index_html = _replace_once(
        index_html,
        "script-src 'self' 'unsafe-inline' 'unsafe-eval';",
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://static.cloudflareinsights.com;",
    )
    index_html = _replace_once(
        index_html,
        "connect-src 'self';",
        "connect-src 'self' https://cloudflareinsights.com;",
    )
    index_path.write_text(index_html, encoding="utf-8")
    return output_root


def prepare_cloudflare_pages(
    *,
    refresh_cache_bust: bool = False,
    analytics_token: str | None = None,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
) -> Path:
    prepare_browser_wheel_module = _load_prepare_browser_wheel_module()
    prepare_browser_wheel_module.prepare_browser_wheel(refresh_cache_bust=refresh_cache_bust)
    token = analytics_token or os.environ.get(ANALYTICS_TOKEN_ENV) or DEFAULT_ANALYTICS_TOKEN
    return build_cloudflare_pages_bundle(output_root=output_root, analytics_token=token)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare the Cloudflare Pages bundle in dist/cloudflare-pages with the browser wheel, "
            "Cloudflare Web Analytics snippet, and hosted-site privacy notice."
        )
    )
    parser.add_argument(
        "--refresh-cache-bust",
        action="store_true",
        help="Refresh GBDRAW_WHEEL_CACHE_BUST while preparing the browser wheel.",
    )
    parser.add_argument(
        "--analytics-token",
        help=(
            "Cloudflare Web Analytics token. Defaults to the value of "
            f"{ANALYTICS_TOKEN_ENV} or the repository default."
        ),
    )
    args = parser.parse_args(argv)
    output_root = prepare_cloudflare_pages(
        refresh_cache_bust=args.refresh_cache_bust,
        analytics_token=args.analytics_token,
    )
    print(f"Prepared Cloudflare Pages bundle: {output_root.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
