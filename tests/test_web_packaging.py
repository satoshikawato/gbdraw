from __future__ import annotations

import base64
import contextlib
import functools
import gzip
import html
import importlib.util
import json
import re
import shutil
import socket
import subprocess
import sys
import threading
import zipfile
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace

import pytest
from PIL import Image

from gbdraw.session_io import CURRENT_SESSION_VERSION, load_session


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
GALLERY_ROOT = WEB_ROOT / "gallery"
README_PATH = REPO_ROOT / "README.md"
ABOUT_PATH = REPO_ROOT / "docs" / "ABOUT.md"
CITATION_PATH = REPO_ROOT / "CITATION.cff"

PREPRINT_TITLE = "gbdraw: a genome diagram generator for microbes and organelles"
PREPRINT_DOI = "10.64898/2026.04.07.716863"
BROWSER_WHEEL_FORBIDDEN_PREFIXES = (
    "gbdraw/web/assets/",
    "gbdraw/web/gallery/",
    "gbdraw/web/js/",
    "gbdraw/web/presets/",
    "gbdraw/web/vendor/",
    "gbdraw/web/wasm/",
)
BROWSER_WHEEL_FORBIDDEN_FILES = {
    "gbdraw/web/index.html",
    "gbdraw/web/open-source-notices.html",
}
BROWSER_WHEEL_REQUIRED_RUNTIME_DATA = {
    "gbdraw/data/color_palettes.toml",
    "gbdraw/data/config.toml",
}
GALLERY_SESSION_FILES = [
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_basic_circular.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "tobacco-chloroplast.gbdraw-session.json",
    "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "lambda_basic_linear.gbdraw-session.json",
]
GALLERY_MULTI_RECORD_LINEAR_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
}
GALLERY_EDITOR_STATE_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "HmmtDNA_ATskew.gbdraw-session.json",
    "WSSV_genome_comparison.gbdraw-session.json",
}
GALLERY_LOSAT_CACHE_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
    "WSSV_genome_comparison.gbdraw-session.json",
    "hepatoplasmataceae_collinear.gbdraw-session.json.gz",
    "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
    "majanivirus_orthogroup.gbdraw-session.json.gz",
    "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
}
GALLERY_LOSAT_DERIVED_CACHE_SESSION_FILES = {
    "BGC0000708-BGC0000713.gbdraw-session.json",
}
def _load_verify_module():
    module_path = REPO_ROOT / "tools" / "verify_gui_offline.py"
    spec = importlib.util.spec_from_file_location("verify_gui_offline", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load verification module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _read_session_prefix(path: Path, length: int = 256) -> str:
    with path.open("rb") as session_file:
        is_gzip = session_file.read(2) == b"\x1f\x8b"
    if is_gzip:
        with gzip.open(path, mode="rt", encoding="utf-8") as session_handle:
            return session_handle.read(length)
    with path.open(encoding="utf-8") as session_handle:
        return session_handle.read(length)


def _load_prepare_browser_wheel_module():
    module_path = REPO_ROOT / "tools" / "prepare_browser_wheel.py"
    spec = importlib.util.spec_from_file_location("prepare_browser_wheel", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load browser wheel preparation module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_prepare_cloudflare_pages_module():
    module_path = REPO_ROOT / "tools" / "prepare_cloudflare_pages.py"
    spec = importlib.util.spec_from_file_location("prepare_cloudflare_pages", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load Cloudflare packaging module from {module_path}")
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


class _WebTestRequestHandler(SimpleHTTPRequestHandler):
    extensions_map = {
        **SimpleHTTPRequestHandler.extensions_map,
        ".css": "text/css; charset=utf-8",
        ".data": "application/octet-stream",
        ".html": "text/html; charset=utf-8",
        ".js": "text/javascript; charset=utf-8",
        ".json": "application/json; charset=utf-8",
        ".mjs": "text/javascript; charset=utf-8",
        ".svg": "image/svg+xml",
        ".wasm": "application/wasm",
        ".whl": "application/octet-stream",
    }

    def log_message(self, format: str, *args: object) -> None:  # noqa: A002
        return


@contextlib.contextmanager
def _serve_repo_root():
    handler = functools.partial(_WebTestRequestHandler, directory=str(REPO_ROOT))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        host, port = server.server_address
        yield f"http://{host}:{port}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


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


def _gui_search_field_ids(source: str) -> list[str]:
    block = source.split("FEATURE_SEARCH_FIELD_DEFINITIONS = Object.freeze([", 1)[1].split("]);", 1)[0]
    return re.findall(r"\{\s*id:\s*'([^']+)'", block)


def _standalone_search_field_ids(source: str) -> list[str]:
    block = source.split("var searchFieldOptions = [", 1)[1].split("];", 1)[0]
    return re.findall(r"\['([^']+)',\s*'[^']+'\]", block)


def _standalone_interactivity_source() -> str:
    service_root = WEB_ROOT / "js" / "services"
    return "\n".join(
        [
            (service_root / "standalone-interactivity.js").read_text(encoding="utf-8"),
            (service_root / "standalone-interactivity-assets.js").read_text(encoding="utf-8"),
        ]
    )


def _gallery_svg_metadata(svg_source: str) -> dict[str, object]:
    metadata_match = re.search(
        r'<metadata(?P<attributes>[^>]*id="gbdraw-interactive-feature-metadata"[^>]*)>'
        r'(?P<payload>.*?)</metadata>',
        svg_source,
        re.S,
    )
    assert metadata_match, "missing interactive feature metadata"
    payload = html.unescape(metadata_match.group("payload"))
    if 'data-encoding="gzip-base64"' in metadata_match.group("attributes"):
        payload = gzip.decompress(base64.b64decode(payload)).decode("utf-8")
    return json.loads(payload)


def _assert_white_gallery_thumbnail(path: Path) -> None:
    image = Image.open(path).convert("RGB")
    assert image.size == (640, 360)

    width, height = image.size
    corners = [
        image.getpixel((0, 0)),
        image.getpixel((width - 1, 0)),
        image.getpixel((0, height - 1)),
        image.getpixel((width - 1, height - 1)),
    ]
    assert all(min(pixel) >= 245 for pixel in corners)

    border_pixels = []
    for x in range(width):
        border_pixels.append(image.getpixel((x, 0)))
        border_pixels.append(image.getpixel((x, height - 1)))
    for y in range(height):
        border_pixels.append(image.getpixel((0, y)))
        border_pixels.append(image.getpixel((width - 1, y)))
    average_luminance = sum(sum(pixel) / 3 for pixel in border_pixels) / len(border_pixels)
    assert average_luminance >= 245


def test_web_offline_assets_can_be_prepared_for_packaging() -> None:
    verify_module, expected_wheel_path = ensure_prepared_browser_wheel()
    expected_wheel_name = "gbdraw-0.14.0b0-py3-none-any.whl"
    assert verify_module._parse_wheel_name() == expected_wheel_name
    assert expected_wheel_path.name == expected_wheel_name
    verify_module.assert_browser_wheel_is_not_recursive(expected_wheel_path)
    verify_module.BUILD_SUPPORT.inspect_browser_wheel_payload(expected_wheel_path)
    verify_module._assert_packaged_assets()


def test_browser_wheel_excludes_hosted_web_assets() -> None:
    verify_module, browser_wheel_path = ensure_prepared_browser_wheel()
    with zipfile.ZipFile(browser_wheel_path) as zf:
        names = set(zf.namelist())

    forbidden = sorted(
        name
        for name in names
        if name in BROWSER_WHEEL_FORBIDDEN_FILES
        or any(name.startswith(prefix) for prefix in BROWSER_WHEEL_FORBIDDEN_PREFIXES)
        or (name.startswith("gbdraw/web/gbdraw-") and name.endswith(".whl"))
    )
    assert forbidden == []
    assert BROWSER_WHEEL_REQUIRED_RUNTIME_DATA <= names
    assert browser_wheel_path.stat().st_size <= verify_module.BUILD_SUPPORT.BROWSER_WHEEL_MAX_BYTES


def test_local_web_package_data_excludes_gallery_assets() -> None:
    build_support = _load_verify_module().BUILD_SUPPORT
    package_data_patterns = build_support.get_package_data_patterns(include_browser_wheel=True)

    assert all("web/gallery" not in pattern for pattern in package_data_patterns)
    assert "gbdraw/web/gallery" not in (REPO_ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "gbdraw/web/gallery/" in build_support._BROWSER_WHEEL_FORBIDDEN_PREFIXES


def test_index_links_to_open_source_notices() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "./open-source-notices.html" in index_html
    assert "Open Source Notices" in index_html


def test_index_links_to_hosted_interactive_gallery() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "Interactive Gallery" in index_html
    assert "./gallery/" not in index_html

    match = re.search(r'<a href="https://gbdraw\.app/gallery/"(?P<attrs>[^>]*)>Interactive Gallery</a>', index_html)
    assert match is not None
    assert 'target="_blank"' in match.group("attrs")
    assert 'rel="noopener noreferrer"' in match.group("attrs")


def test_linear_record_selector_source_contract() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    helper_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    app_setup_js = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    watcher_js = (WEB_ROOT / "js" / "app" / "watchers.js").read_text(encoding="utf-8")

    assert '<select\n                                                v-model="seq.region_record_id"' in index_html
    assert '<input type="text" v-model="seq.region_record_id"' not in index_html
    assert "Record (optional)" in index_html
    assert "linearRecordOptions(seq)" in index_html
    assert "linearRecordSelectorDisabled(seq)" in index_html
    assert "def list_sequence_records(path, format):" in helper_js
    assert "def list_gff_fasta_records(gff_path, fasta_path):" in helper_js
    assert "load_gff_fasta(" in helper_js
    assert 'format_map = {"genbank": "genbank", "fasta": "fasta"}' in helper_js
    assert "def list_genbank_records(" not in helper_js
    assert (WEB_ROOT / "js" / "app" / "record-discovery.js").is_file()
    assert (WEB_ROOT / "js" / "app" / "linear-record-selector.js").is_file()
    assert "createLinearRecordSelector" in app_setup_js
    assert "refreshLinearRecordSelectors" in watcher_js


def test_web_run_info_tab_source_contract() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")

    assert "Preview" in index_html
    assert "Run info" in index_html
    assert "copyRunCommand" in index_html
    assert "lastRunInfo" in index_html
    assert "const resultPanelTab = ref('preview');" in state_js
    assert "const lastRunInfo = ref(null);" in state_js
    assert "buildRunInfo({" in run_analysis_js
    assert "resultPanelTab.value = 'preview';" in run_analysis_js
    assert "if (activePaletteName !== 'default') args.push('-p', activePaletteName);" in run_analysis_js
    assert "const dContent = buildDefaultColorOverrideTsv({" in run_analysis_js
    assert "if (dContent.trim() !== '') {" in run_analysis_js
    assert "cliInvocation: exportableCliInvocation" in config_js


def test_web_run_info_helper_builds_display_commands() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/run-info.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_cli_arg_helpers_omit_default_values() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/cli-args.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_canonical_session_request_codec() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/session-request.test.mjs"], check=True, cwd=REPO_ROOT)
    subprocess.run(
        [
            node,
            "tests/web/session-request.test.mjs",
            "--project-session",
            str(
                GALLERY_ROOT
                / "sessions"
                / "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz"
            ),
        ],
        check=True,
        cwd=REPO_ROOT,
    )
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(
        encoding="utf-8"
    )
    assert "files: serializedFiles" not in config_source


def test_web_losat_settings_preserve_requested_thread_count() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/losat-settings.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_file_import_helpers_preserve_specific_rule_order() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/file-imports.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_feature_color_caption_scope_updates_specific_rule() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/feature-color-actions.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_feature_visibility_helpers() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/feature-visibility.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_feature_selector_helpers() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/feature-selector.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_specific_rule_qualifier_accepts_suggestions_and_custom_values() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(
        encoding="utf-8"
    )

    assert '<datalist id="specific-rule-qualifier-options">' in index_html
    assert index_html.count('list="specific-rule-qualifier-options"') == 2
    assert 'v-for="value in specificRuleQualifierSuggestions"' in index_html
    assert 'v-model="newSpecRule.qual"' in index_html
    assert '<select v-model="newSpecRule.qual"' not in index_html
    assert "collectSpecificColorQualifierSuggestions" in state_source
    assert "specificRuleQualifierSuggestions" in app_setup_source


def test_web_preview_runtime_helpers() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/preview-runtime.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_session_feature_metadata_helpers() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run([node, "tests/web/session-feature-metadata.test.mjs"], check=True, cwd=REPO_ROOT)


def test_web_feature_visibility_table_uses_matching_exclusion_mode() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    helper_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    worker_js = (WEB_ROOT / "js" / "workers" / "diagram-generation-worker.js").read_text(encoding="utf-8")
    svg_actions_js = (WEB_ROOT / "js" / "app" / "feature-editor" / "svg-actions.js").read_text(encoding="utf-8")

    assert '<option value="exclude_matching">Exclude from matching' in index_html
    assert '<option value="suppress">Suppress' not in index_html
    assert "featureVisibilityFeatureSuggestions" in index_html
    assert "Feature Visibility Scope" in index_html
    assert "feature-visibility-feature-type-options" not in index_html
    assert "'transcript'" in state_js
    assert "{svg_id: 'on' | 'off' | 'exclude_matching'}" in state_js
    assert "/web_feature_visibility_table.tsv" in run_analysis_js
    assert "args.push('--feature_visibility_table', featureVisibilityTablePath);" in run_analysis_js
    assert "serializeFeatureVisibilityRules(featureVisibilityRules?.value || [])" in run_analysis_js
    assert "featureVisibility: featureVisibilityCacheKey" in run_analysis_js
    assert "featureVisibilityTsv: featureVisibilityCacheKey" in run_analysis_js
    assert "featureVisibilityTablePath || null" in worker_js
    assert "feature_visibility_table_path=None" in helper_js
    assert "extract_features_from_genbank(gb_path, region_spec=None, record_selector=None, reverse_flag=None, selected_features=None, feature_visibility_table_path=None)" in helper_js
    assert "feature_visibility_rules=feature_visibility_rules" in helper_js
    assert "if (mode === 'off')" in svg_actions_js
    assert "mode === 'off' || mode === 'suppress'" not in svg_actions_js


def test_web_losatp_orthogroup_membership_uses_anchor_core_model() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    normalization_js = (WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")

    assert '<option value="distribution_split">Distribution split</option>' not in index_html
    assert '<option value="family_merge">Family merge</option>' not in index_html
    assert "orthogroupMembershipMode: 'anchor_core_v1'" in state_js
    assert "outparalog_split: 'anchor_core_v1'" in normalization_js
    assert "normalizeOrthogroupMembershipMode" in config_js
    assert "normalizeOrthogroupMembershipMode" in run_analysis_js


def test_web_linear_definition_line_styles_contract() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    app_setup_js = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    cli_args_js = (WEB_ROOT / "js" / "app" / "cli-args.js").read_text(encoding="utf-8")

    assert "Definition Line Styles" in index_html
    assert "Subtitle / title (optional)" in index_html
    assert "v-model=\"seq.record_subtitle\"" in index_html
    assert "Name / Species" in app_setup_js
    assert "Subtitle" in app_setup_js
    assert "Length / Coord." in app_setup_js
    assert ">Normal</button>" in index_html
    assert "data-definition-line-kind" in state_js
    assert "record_subtitle: String(source.record_subtitle ?? '')" in state_js
    assert "linear_definition_line_styles: createDefaultLinearDefinitionLineStyles()" in state_js
    assert "record_subtitle: seq.record_subtitle ?? ''" in config_js
    assert "normalizeDefinitionLineStyleState" in config_js
    assert "buildDefinitionLineStyleAssignments" in run_analysis_js
    assert "args.push('--record_subtitle', subtitle);" in run_analysis_js
    assert "definition_line_style" in run_analysis_js
    assert "args.push('--definition_line_style', assignment);" in run_analysis_js
    assert "'subtitle'" in cli_args_js
    assert "color=${style.fill}" in cli_args_js


def test_web_collinear_blocks_use_rbh_evidence_scope_ui() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    normalizer_js = (WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    helper_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")

    assert "Edge mode" not in index_html
    assert 'v-model="losat.blastp.collinearAnchorMode"' not in index_html
    assert "Top1" not in index_html
    assert "All hits" not in index_html
    assert "Evidence scope" in index_html
    assert "ribbons are still emitted only between adjacent display rows" in index_html
    assert "Merge conflicts" in index_html
    assert "export const normalizeCollinearAnchorMode = (_value) => 'rbh';" in normalizer_js
    assert "delete cloned.blastp.collinearAnchorMode;" in config_js
    assert "collinearBlockMergeGap" not in state_js
    assert "collinearSingletonMergeGap" not in state_js
    assert "delete state.losat.blastp.collinearBlockMergeGap;" in config_js
    assert "delete state.losat.blastp.collinearSingletonMergeGap;" in config_js
    assert "losat.blastp.collinearAnchorMode = normalizeCollinearAnchorMode" in run_analysis_js
    assert "losat.blastp.collinearMaxConflictsInMergeGap" in run_analysis_js
    assert "max_conflicts=_collinear_int(collinear_max_conflicts_in_merge_gap, 1)" in helper_js
    assert "collinear_block_merge_gap=50" not in helper_js
    assert "collinear_singleton_merge_gap=25" not in helper_js
    assert "normalized_collinear_anchor_mode = \"rbh\"" in helper_js
    assert "normalized_collinear_anchor_mode," in helper_js
    assert "'data-group-kind'" in state_js
    assert "'data-group-scope'" in state_js
    assert "'data-collinear-group-scope'" in state_js


def test_web_orthogroup_payload_serializes_record_local_scope() -> None:
    helpers_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    member = SimpleNamespace(
        orthogroup_id="og_2",
        protein_id="a0",
        source_protein_id=None,
        record_index=0,
        record_id="record_a",
        feature_index=0,
        label="a0",
        feature_svg_id="feature-a0",
        start=0,
        end=9,
        strand=1,
        representative=True,
        role="local_paralog",
        confidence="high",
        assignment_reason="record-local reciprocal paralog cluster",
        supporting_edges=("a0->a1:record_local_paralog",),
        best_core_support=10.0,
        second_best_core_support=0.0,
        gene=None,
        product=None,
        note=None,
        locus_tag=None,
        gene_id=None,
        old_locus_tag=None,
    )
    orthogroups = SimpleNamespace(
        orthogroups={"og_2": [member]},
        names_by_orthogroup_id={},
        descriptions_by_orthogroup_id={},
        name_candidates_by_orthogroup_id={},
        confidence_by_orthogroup_id={},
        rbh_orthogroups={"og_2": ("a0",)},
        ortholog_edges_by_orthogroup_id={},
        ortholog_paths_by_orthogroup_id={},
        related_edges_by_orthogroup_id={},
        scope_by_orthogroup_id={"og_2": "record_local"},
        source_record_index_by_orthogroup_id={"og_2": 0},
    )

    payload = namespace["_serialize_orthogroups_payload"](orthogroups)

    assert payload[0]["scope"] == "record_local"
    assert payload[0]["source_record_index"] == 0
    assert "sourceRecordIndex" not in payload[0]
    assert payload[0]["member_count"] == 1
    assert payload[0]["record_coverage_count"] == 1
    assert payload[0]["members"][0]["role"] == "local_paralog"
    assert payload[0]["members"][0]["strand"] == "+"


def test_web_losatp_orthogroup_members_use_absolute_region_coordinates() -> None:
    helpers_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    groups = [
        {
            "members": [
                {
                    "proteinId": "p1",
                    "start": 0,
                    "end": 457,
                    "strand": "-",
                }
            ]
        }
    ]
    metadata = namespace["_web_orthogroup_member_display_metadata"](
        [
            {
                "protein_map": {
                    "p1": {
                        "protein_id": "p1",
                        "start": 0,
                        "end": 457,
                        "strand": -1,
                        "coord_base": 10000,
                        "coord_step": 1,
                        "coord_length": 1000,
                    }
                },
                "view_transform": {"length": 1000, "reverse": False},
            }
        ]
    )

    namespace["_apply_web_orthogroup_member_display_metadata"](groups, metadata)

    assert groups[0]["members"][0]["start"] == 9999
    assert groups[0]["members"][0]["end"] == 10456
    assert groups[0]["members"][0]["strand"] == "-"


def test_web_collinear_payload_splits_local_groups_from_global_orthogroups() -> None:
    helpers_js = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")

    assert 'group_scope = "global_collinear" if search_scope == "all" else "adjacent_local"' in helpers_js
    assert 'converted["group_kind"] = "orthogroup" if group_scope == "global_collinear" else "collinear_gene_group"' in helpers_js
    assert 'converted["collinear_group_scope"] = group_scope' in helpers_js
    assert 'orthogroup_groups = serialized_groups if group_scope == "global_collinear" else []' in helpers_js
    assert 'collinear_groups = serialized_groups if group_scope == "adjacent_local" else []' in helpers_js
    assert '"orthogroups": orthogroup_groups' in helpers_js
    assert '"collinearGroups": collinear_groups' in helpers_js
    assert '"collinearGroupScope": group_scope' in helpers_js


def test_web_record_local_orthogroup_scope_survives_state_and_ui_layers() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    orthogroups_js = (WEB_ROOT / "js" / "app" / "orthogroups.js").read_text(encoding="utf-8")
    normalizer_js = (WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8")

    assert "state.orthogroups.value = groups;" in config_js
    assert "orthogroupScope" in config_js
    assert "orthogroupScope" in run_analysis_js
    assert "Species-specific orthogroup" in normalizer_js
    assert "Collinearity-backed global evidence" in normalizer_js
    assert "Local collinear group" in normalizer_js
    assert "groupMetadataScopeLabel(orthogroupScope(groupOrId))" in orthogroups_js
    assert "orthogroupScopeLabel(group)" in index_html
    assert "orthogroupScopeLabel(selectedOrthogroup)" in index_html


def test_web_run_analysis_orthogroup_top_label_mode_is_wired() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")

    assert '<option value="orthogroup_top"' in index_html
    assert "Top Orthogroup Record" in index_html
    assert "losat.blastp.mode === 'orthogroup' || losat.blastp.mode === 'collinear'" in index_html
    assert "if (form.show_labels_linear === 'orthogroup_top') args.push('orthogroup_top');" in run_analysis_js
    assert "show_labels_linear: 'none'" in state_js
    assert "form: state.form" in config_js


def test_web_losatp_derived_payload_cache_is_persisted_separately() -> None:
    state_js = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")

    assert "const losatDerivedCache = ref(new Map());" in state_js
    assert "losatDerivedCache:" in config_js
    assert "serializeLosatDerivedCache()" in config_js
    assert "applyLosatDerivedCache(data.losatDerivedCache?.entries)" in config_js
    assert "kind: 'derived-losatp-payload'" in config_js
    assert "buildLosatDerivedPayloadCachePayload({" in run_analysis_js
    assert "getLosatDerivedCacheEntry(derivedCacheMap, derivedCacheKey)" in run_analysis_js
    assert "setLosatDerivedCacheEntry(derivedCacheMap, derivedCacheKey" in run_analysis_js


def test_web_losatp_orthogroup_and_collinear_blastp_omit_hsp_cap() -> None:
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")

    assert (
        "if (!useOrthogroupBlastp && !useCollinearBlastp) {\n"
        "              pushArg(args, '--max-hsps-per-subject', 1);\n"
        "            }"
    ) in run_analysis_js
    assert "pushArg(args, '--max-target-seqs', getBlastpCandidateLimit());" in run_analysis_js


def test_public_web_html_entrypoints_are_not_gitignored() -> None:
    required_html_paths = [
        "gbdraw/web/index.html",
        "gbdraw/web/open-source-notices.html",
        "gbdraw/web/gallery/index.html",
    ]
    gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")

    for path in required_html_paths:
        assert f"!{path}" in gitignore

    if shutil.which("git") and (REPO_ROOT / ".git").exists():
        for path in required_html_paths:
            result = subprocess.run(
                ["git", "check-ignore", "-q", path],
                cwd=REPO_ROOT,
                check=False,
            )
            assert result.returncode == 1, f"{path} must be commit-visible for hosted builds"


def test_gallery_csp_allows_same_origin_tutorial_media() -> None:
    gallery_index = (WEB_ROOT / "gallery" / "index.html").read_text(encoding="utf-8")
    assert "media-src 'self';" in gallery_index
    assert "img-src 'self' data:;" in gallery_index


def test_open_source_notices_are_generated() -> None:
    subprocess.run(
        [sys.executable, "tools/generate_open_source_notices.py", "--check"],
        cwd=REPO_ROOT,
        check=True,
    )


def test_open_source_notices_include_current_project_version() -> None:
    from gbdraw._build_support import read_project_version

    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    assert read_project_version() in notices_html


def test_open_source_notices_omit_internal_generation_details() -> None:
    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    assert "Distribution Summary" not in notices_html
    assert "Inventory Sources" not in notices_html
    assert "Project version:" not in notices_html
    assert "tools/generate_open_source_notices.py" not in notices_html


def test_open_source_notices_include_local_pyodide_wheels() -> None:
    from gbdraw._web_assets import PYODIDE_LOCAL_WHEELS

    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    for wheel_path in PYODIDE_LOCAL_WHEELS:
        assert wheel_path.as_posix() in notices_html


def test_web_config_pyodide_local_wheels_match_shared_manifest() -> None:
    from gbdraw._web_assets import PYODIDE_LOCAL_WHEELS

    config_js = (WEB_ROOT / "js" / "config.js").read_text(encoding="utf-8")
    match = re.search(r"export const PYODIDE_LOCAL_WHEELS\s*=\s*\[(.*?)\];", config_js, re.DOTALL)
    assert match is not None
    configured_wheels = tuple(
        Path(raw.lstrip("./"))
        for raw in re.findall(r"""["']([^"']+\.whl)["']""", match.group(1))
    )
    assert configured_wheels == PYODIDE_LOCAL_WHEELS


def test_index_uses_title_logo_separately_from_icon_assets() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert './assets/gbdraw-logo-title.png' in index_html
    assert '<link rel="icon" href="./assets/gbdraw-logo.svg" type="image/svg+xml">' in index_html


def test_circular_gff3_input_renders_single_gff3_uploader() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert index_html.count('label="GFF3 File" v-model="files.c_gff"') == 1
    assert '<file-uploader label="GenBank/DDBJ File" v-model="files.c_gb"' in index_html
    assert '<file-uploader label="GFF3 File" v-model="files.c_gff" accept=".gff,.gff3,.txt,text/plain,text/*"></file-uploader>' in index_html
    assert '<file-uploader label="FASTA File" v-model="files.c_fasta" accept=".fa,.fas,.fasta,.fna,.ffn,.faa,.txt,text/plain,text/*"></file-uploader>' in index_html


def test_meta_csp_omits_frame_ancestors_header_only_directive() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    notices_html = (WEB_ROOT / "open-source-notices.html").read_text(encoding="utf-8")
    gallery_html = (GALLERY_ROOT / "index.html").read_text(encoding="utf-8")
    assert "frame-ancestors" not in index_html
    assert "frame-ancestors" not in notices_html
    assert "frame-ancestors" not in gallery_html


def test_interactive_gallery_shell_is_static_and_sandboxed() -> None:
    gallery_html = (GALLERY_ROOT / "index.html").read_text(encoding="utf-8")
    gallery_js = (GALLERY_ROOT / "gallery.js").read_text(encoding="utf-8")
    gallery_css = (GALLERY_ROOT / "gallery.css").read_text(encoding="utf-8")

    assert "default-src 'self';" in gallery_html
    assert "script-src 'self';" in gallery_html
    assert "style-src 'self';" in gallery_html
    assert "frame-src 'self';" in gallery_html
    assert '<script type="module" src="./gallery.js"></script>' in gallery_html
    assert 'sandbox="allow-scripts allow-downloads"' in gallery_html
    assert "allow-same-origin" not in gallery_html
    assert 'id="demo-frame"' in gallery_html
    assert 'id="session-link"' in gallery_html
    assert 'id="tag-filter-list"' in gallery_html
    assert 'id="clear-tag-filters"' in gallery_html
    assert "fetch('./examples.json'" in gallery_js
    assert "activeTagFilters" in gallery_js
    assert "filteredExamples" in gallery_js
    assert "frame.src = sample.svg" in gallery_js
    assert "sample.session" in gallery_js
    assert "Vnig_TUMSAT-TG-2018" in gallery_js
    assert "lambda-phage-linear" not in gallery_js
    assert "lambda-phage-linear.svg" not in gallery_html
    assert "hepatoplasmataceae-comparison.svg" not in gallery_html
    assert re.search(r"\.viewer-panel\s*>\s*\*\s*\{[^}]*min-width:\s*0;", gallery_css, re.S)
    assert re.search(
        r"\.frame-wrap\s*\{[^}]*width:\s*100%;[^}]*max-width:\s*100%;[^}]*min-width:\s*0;",
        gallery_css,
        re.S,
    )
    assert re.search(
        r"\.demo-frame\s*\{[^}]*width:\s*100%;[^}]*max-width:\s*100%;[^}]*min-width:\s*0;",
        gallery_css,
        re.S,
    )
    combined = "\n".join([gallery_html, gallery_js, gallery_css]).lower()
    assert "pyodide" not in combined
    assert "vue" not in combined
    assert "tailwind" not in combined


def test_interactive_gallery_examples_are_wired() -> None:
    expected_ids = [
        "HmmtDNA_basic_circular",
        "lambda_basic_linear",
        "HmmtDNA_ATskew",
        "tobacco-chloroplast",
        "Vnig_TUMSAT-TG-2018",
        "hepatoplasmataceae_collinear",
        "vibrio-harveyi-group-collinear",
        "hepatoplasmataceae_orthogroup",
        "BGC0000708-BGC0000713",
        "majanivirus_orthogroup",
        "WSSV_genome_comparison",
    ]
    examples = json.loads((GALLERY_ROOT / "examples.json").read_text(encoding="utf-8"))
    expected_tags = {
        "HmmtDNA_basic_circular": ["Circular", "Interactive SVG"],
        "lambda_basic_linear": ["Linear", "Interactive SVG"],
        "HmmtDNA_ATskew": ["Circular", "Interactive SVG"],
        "tobacco-chloroplast": ["Circular", "Interactive SVG"],
        "Vnig_TUMSAT-TG-2018": ["Circular", "Multi-record", "Interactive SVG"],
        "hepatoplasmataceae_collinear": [
            "Linear",
            "Collinear groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "vibrio-harveyi-group-collinear": [
            "Linear",
            "Multi-record",
            "Collinear groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "hepatoplasmataceae_orthogroup": [
            "Linear",
            "Similarity groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "BGC0000708-BGC0000713": [
            "Linear",
            "Similarity groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "majanivirus_orthogroup": [
            "Linear",
            "Similarity groups",
            "LOSAT",
            "Interactive SVG",
        ],
        "WSSV_genome_comparison": ["Circular", "LOSAT", "Interactive SVG"],
    }

    assert [entry["id"] for entry in examples] == expected_ids
    assert [entry["displayOrder"] for entry in examples] == sorted(
        entry["displayOrder"] for entry in examples
    )
    assert [entry["title"] for entry in examples] == [
        "Human mitochondrial genome: first circular figure",
        "Lambda phage: first linear figure",
        "Human mitochondrial genome (AT skew)",
        "<i>Nicotiana tabacum</i> chloroplast genome regions",
        "<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018",
        "Hepatoplasmataceae collinear protein-match blocks",
        "<i>Vibrio</i> Harveyi group multi-record collinearity",
        "Hepatoplasmataceae CDS protein-similarity links",
        "Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.",
        "Majanivirus CDS protein-similarity links",
        "White spot syndrome virus nucleotide-similarity rings",
    ]
    for entry in examples:
        assert entry["title"]
        assert entry["description"]
        assert "difficulty" not in entry
        assert entry["workflow"]
        assert entry["inputSummary"]
        assert "estimatedTime" not in entry
        assert isinstance(entry["displayOrder"], int)
        assert entry["commandKind"] in {"runnable", "provenance"}
        assert entry["commandNote"]
        assert not entry.get("interactiveStep")
        assert entry["tags"] == expected_tags[entry["id"]]
        assert entry["command"].startswith("gbdraw ")
        assert "interactive-svg" not in entry["command"]
        assert entry["fileSizeLabel"]
        assert entry["sourceNote"]
        assert entry["featureSources"]
        assert entry["svg"].startswith("./examples/")
        assert entry["session"].startswith("./sessions/")
        assert entry["thumbnail"].startswith("./thumbnails/")
        assert entry["tutorial"].startswith("./tutorials/")
        assert entry["tutorialStatus"] == "ready"
        assert entry["sourceSession"].startswith("gbdraw/web/gallery/sessions/")
        assert entry["sourceOutput"].startswith("gbdraw/web/gallery/examples/")
        assert entry["sourceFigure"].startswith("gbdraw/web/gallery/sources/")

        svg_path = GALLERY_ROOT / entry["svg"].removeprefix("./")
        session_path = GALLERY_ROOT / entry["session"].removeprefix("./")
        source_figure_path = REPO_ROOT / entry["sourceFigure"]
        thumbnail_path = GALLERY_ROOT / entry["thumbnail"].removeprefix("./")
        assert svg_path.exists()
        assert session_path.exists()
        assert source_figure_path.exists()
        assert thumbnail_path.exists()

        svg_source = svg_path.read_text(encoding="utf-8")
        session_prefix = _read_session_prefix(session_path)
        thumbnail_header = thumbnail_path.read_bytes()[:16]

        assert svg_path.stat().st_size > 1024
        assert session_path.stat().st_size > 1024
        assert session_path.stat().st_size < 100_000_000
        assert source_figure_path.stat().st_size > 1024
        assert thumbnail_path.stat().st_size > 1024
        assert '"format":"gbdraw-session"' in session_prefix
        version_match = re.search(r'"version":(\d+)', session_prefix)
        assert version_match is not None
        assert int(version_match.group(1)) == CURRENT_SESSION_VERSION
        assert "gbdraw-gallery-interactive-script" not in svg_source
        assert "data-gbdraw-gallery" not in svg_source
        assert "window.parent" not in svg_source
        assert "parent.postMessage" not in svg_source
        assert "window.top" not in svg_source
        assert "window.opener" not in svg_source

        if entry.get("svgType") == "static":
            assert 'data-gbdraw-interactive-svg="true"' not in svg_source
            assert "gbdraw-interactive-feature-script" not in svg_source
        else:
            assert 'data-gbdraw-interactive-svg="true"' in svg_source
            assert "gbdraw-interactive-feature-metadata" in svg_source
            assert "gbdraw-interactive-feature-script" in svg_source
            assert 'data-popup-mode="rich"' in svg_source
            assert "data-gbdraw-original-viewbox" in svg_source
            payload = _gallery_svg_metadata(svg_source)
            assert payload["schema"] in {
                "gbdraw-interactive-feature-popup-v1",
                "gbdraw-interactive-feature-popup-v2",
            }
            assert payload["popup_mode"] == "rich"
            features = payload["features"]
            assert features
            assert any(feature.get("qualifiers") for feature in features)
            assert any(feature.get("location_parts") for feature in features)
            assert all(feature.get("nucleotide_sequence") for feature in features)

        assert thumbnail_header.startswith(b"RIFF")
        assert b"WEBP" in thumbnail_header
        _assert_white_gallery_thumbnail(thumbnail_path)

    provenance = [entry for entry in examples if entry["commandKind"] == "provenance"]
    assert [entry["id"] for entry in provenance] == ["WSSV_genome_comparison"]
    assert "not directly runnable" in provenance[0]["commandNote"]
    collinear = next(entry for entry in examples if entry["id"] == "hepatoplasmataceae_collinear")
    assert collinear["command"].count("--losatp_threads") == 1


def test_runnable_gallery_support_downloads_exist() -> None:
    expected_local_downloads = {
        "HmmtDNA_ATskew": {"./files/HmmtDNA_qualifier_priority.tsv"},
        "tobacco-chloroplast": {
            "./files/chloroplast_specific_table.tsv",
            "./files/qualifier_priority.tsv",
            "./files/nicotiana-tabacum-regions.tsv",
        },
        "BGC0000708-BGC0000713": {
            "./files/BGC0000708-BGC0000713_default_colors.tsv",
            "./files/BGC0000708-BGC0000713_specific_colors.tsv",
            "./files/BGC0000708-BGC0000713_qualifier_priority.tsv",
        },
        "majanivirus_orthogroup": {
            "./files/modified_default_colors.tsv",
            "./files/majani_custom_color_table.tsv",
        },
    }

    for example_id, hrefs in expected_local_downloads.items():
        tutorial = json.loads(
            (GALLERY_ROOT / "tutorials" / f"{example_id}.json").read_text(encoding="utf-8")
        )
        tutorial_hrefs = {download.get("href") for download in tutorial.get("downloads", [])}
        assert hrefs <= tutorial_hrefs
        for href in hrefs:
            assert (GALLERY_ROOT / href.removeprefix("./")).is_file()


def test_index_cloaks_vue_template_until_mount() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert '<div id="app" v-cloak' in index_html
    assert "[v-cloak] { display: none !important; }" in index_html
    assert 'id="app-boot-splash"' in index_html
    assert "Initializing gbdraw..." in index_html
    assert "#app:not([v-cloak]) + #app-boot-splash" in index_html


def test_download_helpers_stop_synthetic_link_clicks_from_closing_popups() -> None:
    source_paths = [
        WEB_ROOT / "js" / "app" / "app-setup.js",
        WEB_ROOT / "js" / "app" / "orthogroups.js",
        WEB_ROOT / "js" / "app" / "run-analysis.js",
        WEB_ROOT / "js" / "app" / "feature-editor" / "label-actions.js",
        WEB_ROOT / "js" / "services" / "export.js",
        WEB_ROOT / "js" / "services" / "standalone-interactivity.js",
    ]
    for source_path in source_paths:
        source = source_path.read_text(encoding="utf-8")
        for match in re.finditer(r"link\.click\(\);", source):
            preceding_source = source[max(0, match.start() - 260):match.start()]
            assert "link.addEventListener('click'" in preceding_source, source_path
            assert "event.stopPropagation();" in preceding_source, source_path


def test_index_includes_preprint_citation() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "How to cite" in index_html
    assert PREPRINT_TITLE in index_html
    assert PREPRINT_DOI in index_html


def test_feature_popup_metadata_ui_is_wired_without_new_dependencies() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    svg_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "svg-actions.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    helper_source = (WEB_ROOT / "js" / "app" / "python-helpers.js").read_text(encoding="utf-8")
    feature_metadata_source = (REPO_ROOT / "gbdraw" / "web_support" / "feature_metadata.py").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")

    assert "rich_feature_popup: true" in state_source
    assert 'v-model="adv.rich_feature_popup"' in index_html
    assert "Rich Feature Popup" in index_html
    assert "feature-popup--simple" in index_html
    assert "feature-hover-summary" in index_html
    assert "!adv.rich_feature_popup || clickedFeature.activeTab === 'edit'" in index_html
    assert "adv?.rich_feature_popup === false ? 440 : 720" in svg_actions_source
    assert "const scheduleHoverSummary = (feat, featureElement, eventLike) => {" in svg_actions_source
    assert "svg.addEventListener('mousemove', handleMouseMove);" in svg_actions_source
    assert "function hideHoverSummary()" in svg_actions_source
    assert "state.adv.rich_feature_popup = data?.adv?.rich_feature_popup !== false;" in config_source
    assert "clickedFeature.activeTab" in index_html
    assert "Details" in index_html
    assert "Qualifiers" in index_html
    assert "Sequence" in index_html
    assert "ph ph-copy" in index_html
    assert "copyText(row.value)" in index_html
    assert "clickedFeature.nucleotideFasta" in index_html
    assert "clickedFeature.aminoAcidFasta" in index_html
    assert "qualifierRows" in svg_actions_source
    assert "locationParts" in svg_actions_source
    assert "nucleotideSequence" in svg_actions_source
    assert "aminoAcidSequence" in svg_actions_source
    assert "buildFeatureSequenceFastas" in svg_actions_source
    assert "import { getFeatureCaption, normalizeStringArray, resolveDisplayProteinId } from '../feature-utils.js';" in svg_actions_source
    assert "const proteinId = resolveDisplayProteinId(feat, member);" in svg_actions_source
    assert "label: 'Protein ID', value: proteinId" in svg_actions_source
    assert "document.elementsFromPoint(eventLike.clientX, eventLike.clientY)" in svg_actions_source
    assert "label: 'Source protein ID'" not in svg_actions_source
    assert "label: 'SVG ID'" not in svg_actions_source
    assert "label: 'Record index'" not in svg_actions_source
    assert "label: 'Strand'" not in svg_actions_source
    assert "navigator.clipboard?.writeText" in app_setup_source
    assert "from gbdraw.web_support.feature_metadata import (" in helper_source
    assert "extract_features_from_genbank_json" in helper_source
    assert "extract_features_from_gff_fasta_json" in helper_source
    assert "return extract_features_from_genbank_json(" in helper_source
    assert "return extract_features_from_gff_fasta_json(" in helper_source
    assert "location_parts" in feature_metadata_source
    assert "nucleotide_sequence" in feature_metadata_source
    assert "amino_acid_sequence" in feature_metadata_source
    assert '"organism": organism' in feature_metadata_source
    assert '"source_protein_id": _first_qualifier_value(feat.qualifiers, "protein_id")' in feature_metadata_source
    assert '"gene_id": _first_qualifier_value(feat.qualifiers, "gene_id")' in feature_metadata_source
    assert '"old_locus_tag": _first_qualifier_value(feat.qualifiers, "old_locus_tag")' in feature_metadata_source
    assert "sanitizeExtractedFeaturesForSession(state.extractedFeatures.value)" in config_source


def test_web_session_feature_metadata_recovery_source_contract() -> None:
    extraction_path = WEB_ROOT / "js" / "app" / "feature-metadata-extraction.js"
    recovery_path = WEB_ROOT / "js" / "app" / "session-feature-metadata.js"
    run_analysis_js = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    config_js = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    extraction_js = extraction_path.read_text(encoding="utf-8")
    recovery_js = recovery_path.read_text(encoding="utf-8")

    assert extraction_path.exists()
    assert "export const extractFeatureMetadataForPreview" in extraction_js
    assert "export const makeLinearRenderedFeatureId" in extraction_js
    assert "export const buildLinearRegionExtractionContext" in extraction_js
    assert "import { extractFeatureMetadataForPreview } from './feature-metadata-extraction.js';" in run_analysis_js
    assert "extractFeatureMetadataForPreview({" in run_analysis_js
    assert "const makeLinearRenderedFeatureId" not in run_analysis_js

    assert recovery_path.exists()
    assert "export const classifyFeatureMetadataState" in recovery_js
    assert "export const collectRenderedFeatureIdentitiesFromSvg" in recovery_js
    assert "export const collectRenderedFeatureIdsFromSvg" in recovery_js
    assert "export const alignRecoveredFeatureIdsToRenderedSvg" in recovery_js
    assert "export const buildSessionFeatureRecoveryPlan" in recovery_js
    assert "export const buildFeatureOverrideMigration" in recovery_js
    assert "export const migrateFeatureOverrideState" in recovery_js
    assert "buildSessionFeatureRecoveryPlan" in config_js
    assert "classifyFeatureMetadataState" in config_js

    import_session_source = config_js.split("export const importSession", 1)[1]
    recovery_call = import_session_source.index("await recoverSessionFeatureMetadataIfNeeded")
    assert import_session_source.index("applyOrthogroupStateData(") < recovery_call
    assert import_session_source.index("applyEditorStateData(data.editorState);") < recovery_call

    export_session_source = config_js.split("export const exportSession", 1)[1]
    assert export_session_source.index("await guardSessionFeatureMetadataForExport();") < export_session_source.index("const sessionData = {")

    assert "READY_MATCH_RATIO" not in recovery_js
    assert "exactMatchingCount" in recovery_js
    assert "aliasMatchingCount" in recovery_js
    assert "missingExactCount" in recovery_js
    assert "featureStableCandidate" in recovery_js
    assert "stableRecordKey" in recovery_js
    assert "data-gbdraw-stable-feature-id" in recovery_js
    assert "data-gbdraw-record-index" in recovery_js
    assert "data-gbdraw-record-id" in recovery_js
    assert "mode === 'circular' && cInputType === 'gb'" in recovery_js
    assert "mode === 'linear' && lInputType === 'gb'" in recovery_js
    assert "featureVisibilityTsv" in recovery_js
    assert "nextFeatureState.featureColorOverrides = rewriteOverrideMap(" in recovery_js
    assert "featureIdMaps" in recovery_js
    assert "nextFeatureState.featureVisibilityOverrides = rewriteOverrideMap(" in recovery_js
    assert "svgIdMaps" in recovery_js
    assert "nextEditorState.featureStrokes.overrides = rewriteOverrideMap(" in recovery_js


def test_gallery_sessions_ship_resumable_state_without_duplicate_files() -> None:
    for session_name in GALLERY_SESSION_FILES:
        session_path = GALLERY_ROOT / "sessions" / session_name
        session = load_session(session_path)
        features = session.get("features", {}).get("extractedFeatures", [])
        svg_text = "\n".join(result.get("content", "") for result in session.get("results", []))
        rendered_feature_ids = {
            re.sub(r"__part\d+$", "", match)
            for match in re.findall(r"data-gbdraw-feature-id=[\"']([^\"']+)[\"']", svg_text)
        }
        feature_ids = {
            candidate
            for feature in features
            for candidate in (
                str(feature.get("svg_id") or ""),
                re.sub(r"_record_\d+$", "", str(feature.get("svg_id") or "")),
            )
            if candidate
        }
        pairwise_ids = set(re.findall(r"data-gbdraw-pairwise-match-id=[\"']([^\"']+)[\"']", svg_text))
        collinearity_ids = set(re.findall(r"data-collinearity-block-id=[\"']([^\"']+)[\"']", svg_text))

        assert session.get("version") == CURRENT_SESSION_VERSION, session_name
        assert "files" not in session, session_name
        assert features, session_name
        assert session.get("results"), session_name
        assert "orthogroupState" in session, session_name
        assert rendered_feature_ids, session_name
        assert rendered_feature_ids <= feature_ids, session_name
        if session_name in GALLERY_EDITOR_STATE_SESSION_FILES:
            assert session.get("editorState"), session_name
        if session_name in GALLERY_LOSAT_CACHE_SESSION_FILES:
            assert session.get("losatCache", {}).get("entries"), session_name
        if session_name in GALLERY_LOSAT_DERIVED_CACHE_SESSION_FILES:
            assert session.get("losatDerivedCache", {}).get("entries"), session_name
        if session_name in GALLERY_MULTI_RECORD_LINEAR_SESSION_FILES:
            assert pairwise_ids or collinearity_ids, session_name


def test_tobacco_gallery_session_keeps_chloroplast_region_annotations() -> None:
    session = json.loads(
        (GALLERY_ROOT / "sessions" / "tobacco-chloroplast.gbdraw-session.json").read_text(
            encoding="utf-8"
        )
    )
    annotation_sets = (
        session.get("renderRequest", {})
        .get("diagramOptions", {})
        .get("annotations", {})
        .get("sets", [])
    )
    assert len(annotation_sets) == 1
    assert annotation_sets[0]["id"] == "plastome_regions"
    annotations = annotation_sets[0]["annotations"]
    assert [annotation["label"] for annotation in annotations] == ["LSC", "IRb", "SSC", "IRa"]
    actual_spans = [
        (annotation["target"]["start"], annotation["target"]["end"])
        for annotation in annotations
    ]
    assert actual_spans == [
        (1, 86686),
        (86687, 112029),
        (112030, 130600),
        (130601, 155943),
    ]


def test_vnig_gallery_session_multirecord_positions_are_restoreable() -> None:
    session = load_session(
        GALLERY_ROOT
        / "sessions"
        / "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz"
    )
    expected_positions = ["#1@1", "#2@1", "#3@2", "#4@2", "#5@2", "#6@2"]
    config_positions = session.get("config", {}).get("adv", {}).get("multi_record_positions")
    if isinstance(config_positions, list) and config_positions:
        actual_positions = [
            f"{entry.get('selector')}@{entry.get('row')}"
            for entry in config_positions
        ]
    else:
        args = session.get("cliInvocation", {}).get("args", [])
        actual_positions = [
            str(args[index + 1])
            for index, arg in enumerate(args[:-1])
            if arg == "--multi_record_position"
        ]

    assert actual_positions == expected_positions

    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    assert "hydrateMissingMultiRecordPositionsFromCliInvocation(data.config, data.cliInvocation)" in config_source


def test_feature_sequence_fasta_formatter_uses_ncbi_style_headers(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "feature-sequence-fasta.js"
    feature_utils_path = tmp_path / "feature-utils.mjs"
    feature_utils_path.write_text((WEB_ROOT / "js" / "app" / "feature-utils.js").read_text(encoding="utf-8"), encoding="utf-8")
    module_path = tmp_path / "feature-sequence-fasta.mjs"
    module_path.write_text(
        source_path.read_text(encoding="utf-8").replace("./feature-utils.js", "./feature-utils.mjs"),
        encoding="utf-8",
    )
    check_path = tmp_path / "check-feature-sequence-fasta.mjs"
    check_path.write_text(
        f"""
        import {{ buildFeatureSequenceFastas }} from {module_path.as_uri()!r};

        const assert = (condition, message) => {{
          if (!condition) throw new Error(message);
        }};

        const feature = {{
          record_id: 'NC_000001.1',
          type: 'CDS',
          start: 0,
          end: 9,
          strand: '+',
          organism: 'Example organism',
          qualifiers: {{
            product: ['example protein'],
            protein_id: ['WP_000001.1'],
            locus_tag: ['ABC_0001']
          }},
          nucleotide_sequence: 'ATGAAATAA',
          amino_acid_sequence: 'MK'
        }};
        const fasta = buildFeatureSequenceFastas(feature);
        assert(
          fasta.nucleotideFasta === '>NC_000001.1:1-9 example protein [Example organism]\\nATGAAATAA',
          fasta.nucleotideFasta
        );
        assert(
          fasta.aminoAcidFasta === '>WP_000001.1 example protein [Example organism]\\nMK',
          fasta.aminoAcidFasta
        );

        const fallback = buildFeatureSequenceFastas({{
          record_id: 'seq1',
          start: 0,
          end: 6,
          strand: '-',
          qualifiers: {{
            locus_tag: ['LOC_1'],
            product: ['fallback protein']
          }},
          nucleotide_sequence: 'ATGAAA',
          amino_acid_sequence: 'M'.repeat(61)
        }});
        assert(
          fallback.nucleotideFasta === '>seq1:c6-1 fallback protein\\nATGAAA',
          fallback.nucleotideFasta
        );
        assert(fallback.aminoAcidFasta.startsWith('>LOC_1 fallback protein\\n'), fallback.aminoAcidFasta);
        assert(fallback.aminoAcidFasta.endsWith('\\nM'), fallback.aminoAcidFasta);
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)


def test_interactive_svg_export_decouples_interactivity_from_rich_popup_payload() -> None:
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")
    standalone_source = _standalone_interactivity_source()
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "import { enrichSvgWithStandaloneInteractivity, stripEditorOnlyCursorStyles } from './standalone-interactivity.js';" in export_source
    assert "gbdraw-interactive-feature-popup-v1" not in export_source
    assert "gbdraw-feature-search-controls" not in export_source
    assert "const STANDALONE_INTERACTIVE_SCRIPT" not in export_source
    assert "enrichSvgWithStandaloneFeaturePopup" not in export_source

    standalone_needles = [
        "gbdraw-interactive-feature-popup-v2",
        "gbdraw-interactive-feature-metadata",
        "gbdraw-interactive-feature-script",
        "gbdraw-feature-search-controls",
        "gbdraw-interactive-feature--match",
        "gbdraw-interactive-feature--active-match",
        "gbdraw-interactive-feature--dimmed",
        "gbdraw-interactive-feature-match-glow",
        "filter: url(#gbdraw-interactive-feature-match-glow);",
        "stroke-opacity: 0.6;",
        "stroke-opacity: 1;",
        "function normalizeSearchText(value)",
        "function compileSearchMatcher(query, useRegex)",
        "function buildPreparedSearchIndex()",
        "function preparedFeatureSearchMatches(document, matcher, field, qualifierKey)",
        "matchDetails: {}",
        "data-search-match-detail",
        "gfs-button--clear",
        "gfs-match-detail",
        "['orthogroup', 'Orthogroup']",
        "['nucleotide', 'Nucleotide']",
        "['amino-acid', 'Amino acid']",
        "var NUCLEOTIDE_IUPAC = {",
        "var AMINO_ACID_IUPAC = {",
        "function buildIupacQueryPattern(query, alphabet)",
        "function supportsStandaloneControls()",
        "function setSearchState(nextState)",
        "var pendingSearchState = {",
        "function setPendingSearchState(nextState)",
        "queryInput.addEventListener('input', function () {\n      setPendingSearchState({ query: queryInput.value });",
        "fieldSelect.addEventListener('change', function () {\n      setPendingSearchState({ field: fieldSelect.value });",
        "searchButton.addEventListener('click', function () {\n      setSearchState({",
        "query: pendingSearchState.query",
        "openButton.addEventListener('click', function () {\n      openActiveMatchPopup();",
        "function applySearchResults()",
        "function setActiveMatch(index, options)",
        "function clearSearch()",
        "Search match",
        "Orthogroup members",
        "function scheduleInitialViewportRefresh()",
        "var initialView = copyViewRect(getViewRect());",
        "rectsNearlyEqual(getViewRect(), initialView)",
        "scheduleInitialViewportRefresh();",
        "var targetRect = fitRectToAspect({\n      x: bounds.x + bounds.width / 2 - targetWidth / 2,",
        "}, homeViewRect.width / homeViewRect.height);\n    setSvgViewRect(targetRect);",
        "visibleView.x + visibleView.width - (controlWidth * unit) - margin",
        "visibleView.y + margin + (searchControlsOffsetCss.y * unit)",
        "gfi-og-members-table",
        "Coordinates (+/-)",
        "Product / note",
        "displayProteinId(null, member)",
        "function displayProteinId(feature, member, fallback)",
        "function firstDisplayText()",
        "display_label",
        "search_labels",
        "orthogroup_id",
        "protein_id",
        "const buildStandaloneOrthogroupPayloads = (features, context) => {",
        "const orthogroups = buildStandaloneOrthogroupPayloads(features, context);",
        "data-gbdraw-interactive-feature",
        "data-gbdraw-original-viewbox",
        "data-gbdraw-original-width",
        "data-gbdraw-original-height",
        "export const enrichSvgWithStandaloneInteractivity = (svg, options = {}) => {\n  if (!svg) return false;",
        "const features = buildStandaloneFeaturePayloads(svg, {\n    ...context,\n    popupMode: normalizedPopupMode\n  });",
        "svg.setAttribute('width', '100vw');",
        "svg.setAttribute('height', '100vh');",
        "svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');",
        "svg.style.setProperty('width', '100vw');",
        "svg.style.setProperty('height', '100vh');",
        "function parseOriginalViewRectFromSvg()",
        "function getViewportAspect()",
        "function fitRectToAspect(rect, targetAspect)",
        "var homeViewRect = fitRectToAspect(originalViewRect, getViewportAspect());",
        "homeViewRect.width / maxZoom",
        "homeViewRect.width * nextScale",
        "{ action: 'reset', label: 'Original', title: 'Return to original view', width: 62 }",
        "function ensureStickyLegendBackground(legend, bbox)",
        "gbdraw-sticky-legend-background",
        "setClassToken(svg, 'gbdraw-interactive-pan-enabled', true);",
        "function refitViewportToWindow()",
        "function scheduleViewportRefit()",
        "popup_mode: normalizedPopupMode",
        "orthogroups",
        "var popupMode = payload.popup_mode === 'simple' ? 'simple' : 'rich';",
        "var orthogroups = Array.isArray(payload.orthogroups) ? payload.orthogroups : [];",
        "var richSearchFields = {",
        "if (popupMode === 'simple') {\n      searchFieldOptions = searchFieldOptions.filter",
        "function renderSimplePopup(feature)",
        "if (normalizedPopupMode === 'rich') {\n      Object.assign(payload, {\n        qualifiers,",
        "nucleotide_sequence",
        "amino_acid_sequence",
        "function featureFasta(feature, sequenceKind)",
        "function featureAminoAcidSequence(feature)",
        "getVisibleViewRect()",
        "var visibleView = getVisibleViewRect();",
        "window.addEventListener('scroll', updateViewportControlsPosition, { passive: true });",
        "window.addEventListener('resize', scheduleViewportRefit);",
        "window.visualViewport.addEventListener('scroll', updateViewportControlsPosition, { passive: true });",
        "popupCssWidth",
        "getPopupCssMetrics",
        "var effectiveScaleX = safeScaleX * metrics.zoomScale;",
        "var marginCss = metrics.margin;",
        "var dragZoomScale = getBrowserZoomScale(getViewportClientRect());",
        "var updateActivePopupViewportMetrics = null;",
        "function refreshActivePopupForViewport()",
        "updateActivePopupViewportMetrics();",
        "updateActivePopupViewportMetrics = function () {",
        "gbdraw-interactive-feature-glow",
        "gbdraw-interactive-feature-match-glow",
        "gbdraw-interactive-feature--hover",
        "gbdraw-interactive-orthogroup-link--hover",
        "function setOrthogroupHover(orthogroupId, highlight)",
        "member_rows: orthogroupMemberRows",
        "function renderMatchMemberTable(section, rows)",
        "function memberFeatureSvgId(memberOrRow)",
        "function featureForMember(memberOrRow)",
        "function featureFasta(feature, sequenceKind)",
        "function memberFasta(memberOrRow, sequenceKind)",
        "function memberSequenceFilename(memberOrRow, sequenceKind, orthogroupId)",
        "function renderMemberSequenceActions(memberOrRow, orthogroupId)",
        "renderMemberSequenceActions(row, orthogroupId)",
        "gfi-seq-actions",
        "<th>Seq</th>",
        "function groupFasta(memberRows, sequenceKind)",
        "function groupSequenceFilename(orthogroupId, displayName, sequenceKind)",
        "function downloadText(filename, text, mimeType)",
        "var safeFilename = filename || 'download.txt';",
        "var link = document.createElementNS(XHTML_NS, 'a');",
        "link.setAttribute('download', safeFilename);",
        "data-download-index",
        "var downloadValues = [];",
        "downloadText(payload.filename, payload.text, payload.type);",
        "featureSvgId: standaloneMemberFeatureSvgId(member)",
        "activePopupDrag",
        "activeSearchControlsDrag",
        "gbdraw-feature-hover-popup",
        "function scheduleHoverPopup(feature, svgId, event)",
        "function renderHoverPopupHtml(feature, svgId)",
        "svg.addEventListener('mousemove'",
        "const collectRenderedFeatureEntries = (svg) => {",
        "const buildFallbackStandaloneFeaturePayload = (svgId, entry, captionsByColor) => {",
        "function startSearchControlsDrag(event, root)",
        "document.addEventListener('mouseup', onEnd, true);",
        "window.addEventListener('blur', onEnd);",
        'data-drag-handle="true"',
        "function startPopupDrag(event)",
        "setFeatureHighlight",
        "svg.addEventListener('mouseover'",
        "root.style.transform = 'scale('",
        "overscroll-behavior: contain;",
        "root.addEventListener('wheel', function (rootEvent) {\n      rootEvent.stopPropagation();\n    }, { passive: true });",
    ]
    for needle in standalone_needles:
        assert needle in standalone_source

    standalone_absent = [
        "['sequence', 'Sequence']",
        "queryInput.addEventListener('input', function () {\n      setSearchState({ query: queryInput.value });",
        "fieldSelect.addEventListener('change', function () {\n      setSearchState({ field: fieldSelect.value });",
        "setActiveMatch(searchState.activeIndex < 0 ? 0 : searchState.activeIndex, { center: true",
        "var yOffset = 42 * unit",
        "['Source protein ID'",
        "addStandaloneMatchRow(rows, 'Source protein ID'",
        "addStandaloneMatchRow(rows, 'Feature SVG ID'",
        "addStandaloneMatchRow(summaryRows, 'Match style'",
        "enrichSvgWithStandaloneFeaturePopup",
        "if (!svg || state.adv.rich_feature_popup === false) return false;",
        "{ action: 'pan', label: 'Pan'",
        "{ action: 'legend', label: 'Legend'",
        "--gfi-text-scale",
        "--gfi-font-size",
        "getPopupTextScale",
        "function setPopupTextScale",
        "['SVG ID'",
        "['Record index'",
        "['Strand'",
        "memberNtFasta",
        "memberAaFasta",
        "ntFasta",
        "aaFasta",
    ]
    for needle in standalone_absent:
        assert needle not in standalone_source

    assert "popupMode: state.adv.rich_feature_popup === false ? 'simple' : 'rich'" in export_source
    assert "features: state.extractedFeatures.value" in export_source
    assert "editableLabels: state.editableLabels.value" in export_source
    assert "orthogroups: state.orthogroups.value" in export_source
    assert "if (!svg || state.adv.rich_feature_popup === false) return false;" not in export_source

    zoom_block = standalone_source.split("  function zoomViewBy", 1)[1].split("  function resetViewport", 1)[0]
    assert "closePopup();" not in zoom_block
    popup_resize_block = standalone_source.split("    function startPopupResize", 1)[1].split("    function startPopupDrag", 1)[0]
    popup_drag_block = standalone_source.split("    function startPopupDrag", 1)[1].split("    function redraw", 1)[0]
    assert "document.addEventListener('mousemove', onMove, true);" in popup_resize_block
    assert "document.addEventListener('mouseup', onEnd, true);" in popup_resize_block
    assert "window.addEventListener('mouseup', onEnd, true);" in popup_resize_block
    assert "window.addEventListener('blur', onEnd);" in popup_resize_block
    assert "typeof moveEvent.buttons === 'number' && (moveEvent.buttons & 1) !== 1" in popup_resize_block
    assert "document.addEventListener('mousemove', onMove, true);" in popup_drag_block
    assert "document.addEventListener('mouseup', onEnd, true);" in popup_drag_block
    assert "window.addEventListener('mouseup', onEnd, true);" in popup_drag_block
    assert "window.addEventListener('blur', onEnd);" in popup_drag_block
    assert "typeof moveEvent.buttons === 'number' && (moveEvent.buttons & 1) !== 1" in popup_drag_block

    assert "export const downloadSVG = () => {\n  const svgString = getCurrentSvgString();" in export_source
    assert "export const downloadInteractiveSVG = () => {\n  const svgString = getCurrentSvgString({ interactive: true });" in export_source
    assert "export const downloadPNG = () => {\n  const svgString = getCurrentSvgString();" in export_source
    assert "export const downloadPDF = async () => {\n  const svgString = getCurrentSvgString();" in export_source
    assert "export const downloadSVG = () => {\n  const svgString = getCurrentSvgString({ interactive: true });" not in export_source
    assert "export const downloadPNG = () => {\n  const svgString = getCurrentSvgString({ interactive: true });" not in export_source
    assert "export const downloadPDF = async () => {\n  const svgString = getCurrentSvgString({ interactive: true });" not in export_source
    assert "downloadInteractiveSVG" in app_setup_source
    assert '@click="downloadInteractiveSVG"' in index_html
    assert "Interactive SVG" in index_html
    assert "Browser-oriented SVG with embedded controls. Rich Feature Popup controls how much feature detail is embedded." in index_html


def test_gui_preview_feature_search_is_wired_and_kept_export_transient() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")
    standalone_source = _standalone_interactivity_source()
    search_core_source = (WEB_ROOT / "js" / "app" / "feature-search" / "search-core.js").read_text(encoding="utf-8")
    preview_svg_source = (WEB_ROOT / "js" / "app" / "feature-search" / "preview-svg.js").read_text(encoding="utf-8")

    assert "import { createPreviewFeatureSearch } from './feature-search/preview-actions.js';" in app_setup_source
    assert "const previewFeatureSearch = createPreviewFeatureSearch({" in app_setup_source
    assert "openFeatureEditorForFeature: featureActions.openFeatureEditorForFeature" in app_setup_source
    assert "previewFeatureSearchInput = ref('')" in state_source
    assert "previewFeatureSearchField = ref('all')" in state_source
    assert 'placeholder="Search features"' in index_html
    assert 'v-model="previewFeatureSearchField"' in index_html
    assert 'v-for="field in previewFeatureSearchFieldOptions"' in index_html
    assert '@keydown.enter.prevent="applyPreviewFeatureSearch"' in index_html
    assert '@click.stop="applyPreviewFeatureSearch"' in index_html
    assert "applyPreviewFeatureSearch: previewFeatureSearch.applySearch" in app_setup_source
    assert "goToNextPreviewFeatureSearchMatch" in index_html
    assert "openPreviewFeatureSearchActiveMatch" in index_html
    assert "gbdraw-preview-feature-search-match" in index_html
    assert "gbdraw-preview-feature-search-active-match" in index_html
    assert "gbdraw-preview-feature-search-dimmed" in index_html

    assert _gui_search_field_ids(search_core_source) == _standalone_search_field_ids(standalone_source)
    assert "RICH_FEATURE_SEARCH_FIELD_IDS = Object.freeze([" in search_core_source
    assert "'qualifier-key'" in search_core_source
    assert "'qualifier-value'" in search_core_source
    assert "'nucleotide'" in search_core_source
    assert "'amino-acid'" in search_core_source
    assert "buildIupacQueryPattern" in search_core_source
    assert "featureSearchMatches" in search_core_source
    assert "formatSearchMatchDetail" in search_core_source

    assert "stripPreviewFeatureSearchClasses" in preview_svg_source
    assert "resolvePreviewSvg" in preview_svg_source
    assert "centerPreviewFeature" in preview_svg_source
    assert "import { stripPreviewFeatureSearchClasses } from '../app/feature-search/preview-svg.js';" in export_source
    assert "stripPreviewFeatureSearchClasses(clone);" in export_source
    assert "gbdraw-feature-search-controls" in standalone_source
    assert "export const downloadSVG = () => {\n  const svgString = getCurrentSvgString();" in export_source
    assert "export const downloadPNG = () => {\n  const svgString = getCurrentSvgString();" in export_source
    assert "export const downloadPDF = async () => {\n  const svgString = getCurrentSvgString();" in export_source
    assert "export const downloadInteractiveSVG = () => {\n  const svgString = getCurrentSvgString({ interactive: true });" in export_source


def test_feature_search_core_matches_labels_qualifiers_and_sequence_aliases(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "feature-search" / "search-core.js"
    module_path = tmp_path / "search-core.mjs"
    feature_utils_path = tmp_path / "feature-utils.mjs"
    feature_utils_path.write_text((WEB_ROOT / "js" / "app" / "feature-utils.js").read_text(encoding="utf-8"), encoding="utf-8")
    module_path.write_text(
        source_path.read_text(encoding="utf-8").replace("../feature-utils.js", "./feature-utils.mjs"),
        encoding="utf-8",
    )
    check_path = tmp_path / "check-feature-search.mjs"
    check_path.write_text(
        f"""
        import {{
          buildFeatureSearchIndex,
          featureSearchItems,
          runFeatureSearch
        }} from {module_path.as_uri()!r};

        const assert = (condition, message) => {{
          if (!condition) throw new Error(message);
        }};

        const feature = {{
          svg_id: 'fabc12345',
          displayLabel: 'Edited beta subunit',
          gene: 'rpoB',
          locus_tag: 'b3987',
          product: 'DNA-directed RNA polymerase subunit beta',
          recordId: 'NC_000913.3',
          type: 'CDS',
          start: 10,
          end: 120,
          strand: 'positive',
          location_parts: [{{ display: 'join(11..40, 80..120)' }}],
          qualifiers: {{
            gene: ['rpoB'],
            product: ['DNA-directed RNA polymerase subunit beta'],
            note: ['core enzyme']
          }},
          nucleotideSequence: 'ATGGCN',
          aminoAcidSequence: 'MXX'
        }};
        const renderedFeatureIds = new Set(['fabc12345']);
        const searchIndex = buildFeatureSearchIndex({{
          features: [feature],
          popupMode: 'rich'
        }});
        assert(searchIndex.featureOrder.join(',') === 'fabc12345', 'Prepared feature order missing');
        assert(searchIndex.byId.has('fabc12345'), 'Prepared document missing');

        const productSearch = runFeatureSearch({{
          features: [feature],
          renderedFeatureIds,
          query: 'polymerase',
          field: 'all',
          popupMode: 'rich',
          searchIndex
        }});
        assert(productSearch.matches.join(',') === 'fabc12345', `Product search failed: ${{JSON.stringify(productSearch)}}`);
        assert(
          productSearch.matchDetails.fabc12345.some((detail) => detail.value.includes('polymerase')),
          `Product match details missing product: ${{JSON.stringify(productSearch.matchDetails)}}`
        );

        const recordItems = featureSearchItems(feature, 'record-id', '', {{ popupMode: 'rich' }}).map((item) => item.value);
        assert(recordItems.includes('NC_000913.3'), `recordId alias missing: ${{JSON.stringify(recordItems)}}`);

        const locationItems = featureSearchItems(feature, 'location', '', {{ popupMode: 'rich' }}).map((item) => item.value);
        assert(locationItems.includes('join(11..40, 80..120)'), `location_parts fallback missing: ${{JSON.stringify(locationItems)}}`);

        const nucleotideSearch = runFeatureSearch({{
          features: [feature],
          renderedFeatureIds,
          query: 'ATGGNN',
          field: 'nucleotide',
          popupMode: 'rich'
        }});
        assert(nucleotideSearch.matches.join(',') === 'fabc12345', `IUPAC nucleotide search failed: ${{JSON.stringify(nucleotideSearch)}}`);

        const simpleSearch = runFeatureSearch({{
          features: [feature],
          renderedFeatureIds,
          query: 'core enzyme',
          field: 'all',
          popupMode: 'simple'
        }});
        assert(simpleSearch.matches.length === 0, `Simple popup mode should not search rich qualifier payloads: ${{JSON.stringify(simpleSearch)}}`);

        const translationOnly = {{
          svg_id: 'ftranslation',
          type: 'CDS',
          qualifiers: {{ translation: ['MPEPTIDE'] }}
        }};
        const translationSearch = runFeatureSearch({{
          features: [translationOnly],
          renderedFeatureIds: new Set(['ftranslation']),
          query: 'MPEPTIDE',
          field: 'amino-acid',
          popupMode: 'rich'
        }});
        assert(translationSearch.matches.join(',') === 'ftranslation', `Translation fallback failed: ${{JSON.stringify(translationSearch)}}`);
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)

    standalone_source = _standalone_interactivity_source()
    assert "feature && feature.displayLabel" in standalone_source
    assert "feature && feature.product" in standalone_source
    assert "feature && feature.searchLabels" in standalone_source
    assert "buildFeatureLocation(feature)" in standalone_source
    assert "feature && (feature.nucleotide_sequence || feature.nucleotideSequence)" in standalone_source
    assert "feature && (feature.amino_acid_sequence || feature.aminoAcidSequence)" in standalone_source


def test_orthogroup_match_popup_payload_uses_orthogroup_summary(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    feature_utils_path = tmp_path / "feature-utils.mjs"
    feature_utils_path.write_text((WEB_ROOT / "js" / "app" / "feature-utils.js").read_text(encoding="utf-8"), encoding="utf-8")
    sequence_fasta_path = tmp_path / "feature-sequence-fasta.mjs"
    sequence_fasta_path.write_text(
        (WEB_ROOT / "js" / "app" / "feature-sequence-fasta.js").read_text(encoding="utf-8").replace("./feature-utils.js", "./feature-utils.mjs"),
        encoding="utf-8",
    )
    normalization_path = tmp_path / "losat-normalization.mjs"
    normalization_path.write_text((WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8"), encoding="utf-8")
    match_sequences_path = tmp_path / "match-sequences.mjs"
    match_sequences_path.write_text(
        (WEB_ROOT / "js" / "app" / "match-sequences.js")
        .read_text(encoding="utf-8")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs"),
        encoding="utf-8",
    )
    source_path = WEB_ROOT / "js" / "app" / "pairwise-match-popup.js"
    module_path = tmp_path / "pairwise-match-popup.mjs"
    module_path.write_text(
        source_path.read_text(encoding="utf-8")
        .replace("./feature-utils.js", "./feature-utils.mjs")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs")
        .replace("./losat-normalization.js", "./losat-normalization.mjs")
        .replace("./match-sequences.js", "./match-sequences.mjs"),
        encoding="utf-8",
    )
    check_path = tmp_path / "check-pairwise-popup.mjs"
    check_path.write_text(
        f"""
        import {{ buildPairwiseMatchHoverRows, buildPairwiseMatchPayload }} from {module_path.as_uri()!r};

        const assert = (condition, message) => {{
          if (!condition) throw new Error(message);
        }};

        const attrs = new Map(Object.entries({{
          'data-gbdraw-pairwise-match-id': 'edge_1',
          'data-match-kind': 'orthogroup',
          'data-orthogroup-id': 'og_1',
          'data-query-record-id': 'record_a',
          'data-subject-record-id': 'record_b',
          'data-qstart': '10',
          'data-qend': '40',
          'data-sstart': '90',
          'data-send': '130',
          'data-pairwise-match-style': 'ribbon',
          'data-query-feature-svg-id': 'fq',
          'data-subject-feature-svg-id': 'fs',
          'data-query-protein-id': 'p_internal_query',
          'data-subject-protein-id': 'p_internal_subject',
          'data-identity': '99.0'
        }}));
        const element = {{
          style: {{}},
          getAttribute: (name) => attrs.get(name) || ''
        }};
        const featureLookup = new Map([
          ['fq', {{
            svg_id: 'fq',
            record_id: 'record_a',
            orthogroupId: 'og_1',
            orthogroupMemberCount: 2,
            orthogroupRecordCoverage: 2,
            proteinId: 'p_internal_query',
            sourceProteinId: 'WP_000001.1',
            qualifiers: {{ protein_id: ['WP_000001.1'] }},
            product: 'query product'
          }}],
          ['fs', {{
            svg_id: 'fs',
            record_id: 'record_b',
            orthogroupId: 'og_1',
            orthogroupMemberCount: 2,
            orthogroupRecordCoverage: 2,
            proteinId: 'p_internal_subject',
            qualifiers: {{ protein_id: ['WP_000002.1'] }},
            product: 'subject product'
          }}]
        ]);
        const payload = buildPairwiseMatchPayload(element, {{
          featureLookup,
          orthogroups: [{{
            id: 'legacy_og_1',
            name: 'rpoB',
            member_count: 2,
            record_coverage_count: 2,
            members: [
              {{
                recordId: 'record_a',
                featureSvgId: 'fq',
                start: 9,
                end: 40,
                strand: '+',
                proteinId: 'p_internal_query',
                sourceProteinId: 'WP_000001.1',
                product: 'query product'
              }},
              {{
                recordId: 'record_b',
                featureSvgId: 'fs',
                start: 89,
                end: 130,
                strand: '-',
                proteinId: 'p_internal_subject',
                sourceProteinId: 'WP_000002.1',
                product: 'subject product'
              }}
            ]
          }}]
        }});

        assert(payload.title === 'og_1:rpoB', `Unexpected title: ${{payload.title}}`);
        assert(payload.subtitle === '', `Orthogroup popup should not duplicate subtitle: ${{payload.subtitle}}`);
        assert(payload.sections.map((section) => section.title).join(',') === 'Summary', JSON.stringify(payload.sections));
        const labels = payload.sections.flatMap((section) => section.rows.map((row) => row.label));
        assert(!labels.includes('Match style'), `Match style leaked: ${{JSON.stringify(labels)}}`);
        assert(!labels.includes('Feature SVG ID'), `Feature SVG ID leaked: ${{JSON.stringify(labels)}}`);
        assert(!labels.includes('Source protein ID'), `Source protein ID leaked: ${{JSON.stringify(labels)}}`);
        assert(!labels.includes('Query record'), `Query record should be omitted: ${{JSON.stringify(labels)}}`);
        assert(!labels.includes('Subject record'), `Subject record should be omitted: ${{JSON.stringify(labels)}}`);
        assert(!labels.includes('Query interval'), `Query interval should be omitted: ${{JSON.stringify(labels)}}`);
        assert(!labels.includes('Subject interval'), `Subject interval should be omitted: ${{JSON.stringify(labels)}}`);
        assert(labels.includes('Orthogroup ID'), `Orthogroup ID missing: ${{JSON.stringify(labels)}}`);
        assert(labels.includes('Display name'), `Display name missing: ${{JSON.stringify(labels)}}`);
        const summary = payload.sections[0];
        assert(summary.memberRows.length === 2, JSON.stringify(summary));
        assert(summary.memberRows.map((row) => row.proteinId).join(',') === 'WP_000001.1,WP_000002.1', JSON.stringify(summary.memberRows));
        assert(summary.memberCopyText.includes('Record\\tCoordinates (+/-)\\tProtein ID\\tRole\\tConfidence\\tAssignment reason\\tProduct / note'), summary.memberCopyText);
        const hoverLabels = buildPairwiseMatchHoverRows(payload).map((row) => row.label);
        assert(hoverLabels.includes('Orthogroup'), `Hover orthogroup row missing: ${{JSON.stringify(hoverLabels)}}`);
        assert(!hoverLabels.includes('Query'), `Hover query row should be omitted: ${{JSON.stringify(hoverLabels)}}`);
        assert(!hoverLabels.includes('Subject'), `Hover subject row should be omitted: ${{JSON.stringify(hoverLabels)}}`);

        const fallbackPayload = buildPairwiseMatchPayload(element, {{
          featureLookup,
          orthogroups: []
        }});
        assert(fallbackPayload.title === 'og_1:query product', `Feature fallback title failed: ${{fallbackPayload.title}}`);
        const fallbackLabels = fallbackPayload.sections.flatMap((section) => section.rows.map((row) => row.label));
        assert(fallbackLabels.includes('Members'), `Feature fallback members missing: ${{JSON.stringify(fallbackLabels)}}`);
        assert(fallbackLabels.includes('Record coverage'), `Feature fallback coverage missing: ${{JSON.stringify(fallbackLabels)}}`);
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)


def test_collinear_adjacent_popup_labels_local_collinear_groups(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    feature_utils_path = tmp_path / "feature-utils.mjs"
    feature_utils_path.write_text((WEB_ROOT / "js" / "app" / "feature-utils.js").read_text(encoding="utf-8"), encoding="utf-8")
    sequence_fasta_path = tmp_path / "feature-sequence-fasta.mjs"
    sequence_fasta_path.write_text(
        (WEB_ROOT / "js" / "app" / "feature-sequence-fasta.js").read_text(encoding="utf-8").replace("./feature-utils.js", "./feature-utils.mjs"),
        encoding="utf-8",
    )
    normalization_path = tmp_path / "losat-normalization.mjs"
    normalization_path.write_text((WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8"), encoding="utf-8")
    match_sequences_path = tmp_path / "match-sequences.mjs"
    match_sequences_path.write_text(
        (WEB_ROOT / "js" / "app" / "match-sequences.js")
        .read_text(encoding="utf-8")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs"),
        encoding="utf-8",
    )
    source_path = WEB_ROOT / "js" / "app" / "pairwise-match-popup.js"
    module_path = tmp_path / "pairwise-match-popup.mjs"
    module_path.write_text(
        source_path.read_text(encoding="utf-8")
        .replace("./feature-utils.js", "./feature-utils.mjs")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs")
        .replace("./losat-normalization.js", "./losat-normalization.mjs")
        .replace("./match-sequences.js", "./match-sequences.mjs"),
        encoding="utf-8",
    )
    check_path = tmp_path / "check-collinear-popup.mjs"
    check_path.write_text(
        f"""
        import {{ buildPairwiseMatchHoverRows, buildPairwiseMatchPayload }} from {module_path.as_uri()!r};

        const assert = (condition, message) => {{
          if (!condition) throw new Error(message);
        }};

        const attrs = new Map(Object.entries({{
          'data-gbdraw-pairwise-match-id': 'block_path_1',
          'data-match-kind': 'collinear',
          'data-collinearity-block-id': 'block_1',
          'data-collinearity-block-kind': 'syntenic',
          'data-collinear-group-scope': 'adjacent_local',
          'data-group-kind': 'collinear_gene_group',
          'data-orthogroup-id': 'og_local_1;og_local_2',
          'data-query-record-id': 'record_a',
          'data-subject-record-id': 'record_b',
          'data-qstart': '10',
          'data-qend': '40',
          'data-sstart': '90',
          'data-send': '130',
          'data-identity': '88.5'
        }}));
        const element = {{
          style: {{}},
          getAttribute: (name) => attrs.get(name) || ''
        }};
        const payload = buildPairwiseMatchPayload(element, {{ featureLookup: new Map(), orthogroups: [] }});
        const sectionTitles = payload.sections.map((section) => section.title);
        assert(sectionTitles.includes('Local collinear groups'), JSON.stringify(sectionTitles));
        assert(!sectionTitles.includes('Orthogroups covered'), JSON.stringify(sectionTitles));
        const groupSection = payload.sections.find((section) => section.title === 'Local collinear groups');
        const groupLabels = groupSection.rows.map((row) => row.label);
        assert(groupLabels.includes('Number of local collinear groups'), JSON.stringify(groupLabels));
        assert(payload.blockOrthogroupCount === 2, `Unexpected group count: ${{payload.blockOrthogroupCount}}`);
        const detailLabels = payload.blockOrthogroups[0].detailRows.map((row) => row.label);
        assert(detailLabels.includes('Collinear group ID'), JSON.stringify(detailLabels));
        assert(!detailLabels.includes('Orthogroup ID'), JSON.stringify(detailLabels));
        const hoverLabels = buildPairwiseMatchHoverRows(payload).map((row) => row.label);
        assert(hoverLabels.includes('Collinear groups'), JSON.stringify(hoverLabels));
        assert(!hoverLabels.includes('Orthogroups'), JSON.stringify(hoverLabels));
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)


def test_plain_svg_export_strips_editor_only_cursor_affordances() -> None:
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")
    standalone_source = _standalone_interactivity_source()

    assert "export const stripEditorOnlyCursorStyles = (svg) => {" in standalone_source
    assert "svg.querySelectorAll('[style]').forEach((element) => {" in standalone_source
    assert "if (!style || !/\\bcursor\\s*:/i.test(style)) return;" in standalone_source
    assert "element.style.removeProperty('cursor');" in standalone_source
    assert "if (!element.getAttribute('style')?.trim()) {" in standalone_source
    assert "import { enrichSvgWithStandaloneInteractivity, stripEditorOnlyCursorStyles } from './standalone-interactivity.js';" in export_source
    assert "  } else {\n    stripEditorOnlyCursorStyles(clone);\n  }\n  return new XMLSerializer().serializeToString(clone);" in export_source
    assert "export const downloadInteractiveSVG = () => {\n  const svgString = getCurrentSvgString({ interactive: true });" in export_source


def test_layout_reposition_mode_gates_preview_dragging() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    watcher_source = (WEB_ROOT / "js" / "app" / "watchers.js").read_text(encoding="utf-8")
    diagram_drag_source = (WEB_ROOT / "js" / "app" / "legend-layout" / "diagram-drag.js").read_text(
        encoding="utf-8"
    )
    legend_drag_source = (WEB_ROOT / "js" / "app" / "legend" / "drag-actions.js").read_text(
        encoding="utf-8"
    )
    ui_source = (WEB_ROOT / "js" / "app" / "ui.js").read_text(encoding="utf-8")

    assert "const layoutRepositionMode = ref(false);" in state_source
    assert "layoutRepositionMode," in state_source
    assert '@click="layoutRepositionMode = !layoutRepositionMode"' in index_html
    assert 'v-if="svgContent && layoutRepositionMode"' in index_html
    assert "if (!isLayoutRepositionModeEnabled()) return;" in diagram_drag_source
    assert "const refreshDiagramDragAffordances = () => {" in diagram_drag_source
    assert "const refreshLegendDragAffordances = () => {" in legend_drag_source
    assert "() => layoutRepositionMode.value" in watcher_source
    assert "refreshLegendDragAffordances();" in watcher_source
    assert "refreshDiagramDragAffordances();" in watcher_source
    assert "isSvgEditingTarget(target)" in ui_source
    assert "isLayoutRepositionModeEnabled() && target.closest('svg')" in ui_source


def test_local_index_keeps_google_analytics_as_deploy_only() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "googletagmanager.com/gtag/js" not in index_html
    assert "static.cloudflareinsights.com" not in index_html
    assert "cloudflareinsights.com" not in index_html
    assert "GOOGLE_ANALYTICS_SCRIPT" in index_html
    assert "GOOGLE_ANALYTICS_NOTICE" in index_html
    assert "GBDRAW_HOSTED_BUILD_LABEL" in index_html


def test_web_losat_threaded_browser_wiring() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "config.js").read_text(encoding="utf-8")
    losat_source = (WEB_ROOT / "js" / "services" / "losat.js").read_text(encoding="utf-8")
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")

    assert 'LOSAT_THREADED_WASM_URL = "./wasm/losat/losat-threaded.wasm"' in config_source
    assert "losat-threaded-worker.js" in losat_source
    assert "wasi_thread_start" in (WEB_ROOT / "js" / "workers" / "losat-wasi-thread-worker.js").read_text(encoding="utf-8")
    assert 'v-model="losat.parallelWorkers"' in index_html
    assert 'v-model="losat.threadsPerJob"' in index_html
    assert "executionMode: 'auto'" in (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    assert "buildLosatCachePayload" in run_source
    assert "losatRuntimeCompatibility" not in run_source
    assert "losatThreadsPerJob" not in run_source
    assert "onRuntimeStatus" in run_source


def test_web_losat_thread_count_options_are_contiguous() -> None:
    source = (WEB_ROOT / "js" / "app" / "losat-settings.js").read_text(encoding="utf-8")
    service_source = (WEB_ROOT / "js" / "services" / "losat.js").read_text(encoding="utf-8")
    worker_source = (WEB_ROOT / "js" / "workers" / "losat-threaded-worker.js").read_text(encoding="utf-8")
    assert "const createPositiveIntegerOptions = (maxValue) =>" in source
    assert "const appendRequestedIntegerOption = (options, requestedValue) =>" in source
    assert "return createPositiveIntegerOptions(losatHardwareThreads.value);" in source
    assert "createPositiveIntegerOptions(maxThreads)," in source
    assert "const perJobSlots = losatEffectiveThreadsPerJob.value;" in source
    assert "const workersPerThreadedJob = Math.max(1, threadsPerJob);" in service_source
    assert "getChildWorkerCount(effectiveThreads)" in worker_source
    assert "Array.from({ length: childWorkerCount }" in worker_source


def test_prepare_browser_wheel_refreshes_open_source_notices(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prepare_module = _load_prepare_browser_wheel_module()
    repo_root = tmp_path / "repo"
    web_root = repo_root / "gbdraw" / "web"
    web_root.mkdir(parents=True)
    expected_name = "gbdraw-0.14.0b0-py3-none-any.whl"
    calls: list[object] = []

    def fake_run(args: list[str], *, cwd: Path, env: dict[str, str], check: bool) -> None:
        calls.append("build")
        assert cwd == repo_root
        assert check is True
        assert env["GBDRAW_BUILDING_BROWSER_WHEEL"] == "1"
        outdir = Path(args[args.index("--outdir") + 1])
        outdir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(outdir / expected_name, "w") as zf:
            zf.writestr("gbdraw/__init__.py", "")

    def fake_validate_browser_wheel_prepared() -> Path:
        calls.append("validate")
        return web_root / expected_name

    build_support = SimpleNamespace(
        BROWSER_WHEEL_BUILD_ENV="GBDRAW_BUILDING_BROWSER_WHEEL",
        expected_browser_wheel_name=lambda: expected_name,
        refresh_open_source_notices=lambda: calls.append("notices"),
        generate_cache_bust_token=lambda: "cache-token",
        update_browser_wheel_config=lambda **kwargs: calls.append(("config", kwargs)),
        validate_browser_wheel_prepared=fake_validate_browser_wheel_prepared,
    )
    monkeypatch.setattr(prepare_module, "REPO_ROOT", repo_root)
    monkeypatch.setattr(prepare_module, "WEB_ROOT", web_root)
    monkeypatch.setattr(prepare_module, "_load_build_support_module", lambda: build_support)
    monkeypatch.setattr(prepare_module.subprocess, "run", fake_run)

    assert prepare_module.prepare_browser_wheel(refresh_cache_bust=True) == 0
    assert (web_root / expected_name).exists()
    assert calls == [
        "notices",
        "build",
        ("config", {"wheel_name": expected_name, "cache_bust": "cache-token"}),
        "validate",
    ]


def test_cloudflare_bundle_includes_google_analytics_and_hosted_notice(tmp_path: Path) -> None:
    from gbdraw._build_support import read_project_version

    ensure_prepared_browser_wheel()
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    output_root = tmp_path / "cloudflare-pages"
    commit_sha = "abcdef1234567890abcdef1234567890abcdef12"
    remote_base = (
        "https://raw.githubusercontent.com/satoshikawato/gbdraw/"
        f"{commit_sha}/gbdraw/web/"
    )
    bundle_path = cloudflare_module.build_cloudflare_pages_bundle(
        output_root=output_root,
        commit_sha=commit_sha,
    )

    index_html = (bundle_path / "index.html").read_text(encoding="utf-8")
    assert "https://www.googletagmanager.com/gtag/js?id=G-GG6JMKM02Y" in index_html
    assert "gtag('config', 'G-GG6JMKM02Y');" in index_html
    assert "static.cloudflareinsights.com" not in index_html
    assert "cloudflareinsights.com" not in index_html
    assert "Hosted Site Analytics" in index_html
    assert "uses Google Analytics 4 for aggregate page-usage metrics" in index_html
    assert "Uploaded genome files and generated diagrams are still processed locally in your browser" in index_html
    assert "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://*.googletagmanager.com;" in index_html
    assert (
        "img-src 'self' data: blob: https://*.google-analytics.com "
        "https://*.googletagmanager.com;"
    ) in index_html
    assert (
        "connect-src 'self' https://*.google-analytics.com "
        "https://*.analytics.google.com https://*.googletagmanager.com;"
    ) in index_html
    assert "GOOGLE_ANALYTICS_SCRIPT" not in index_html
    assert "GOOGLE_ANALYTICS_NOTICE" not in index_html
    assert "GBDRAW_HOSTED_BUILD_LABEL" not in index_html
    assert f"Version: v{read_project_version()}+abcdef1" in index_html
    assert f'title="Commit {commit_sha}"' in index_html
    headers = (bundle_path / "_headers").read_text(encoding="utf-8")
    assert "Cross-Origin-Opener-Policy: same-origin" in headers
    assert "Cross-Origin-Embedder-Policy: require-corp" in headers
    assert "Cross-Origin-Resource-Policy: same-origin" in headers
    assert "Content-Security-Policy: frame-ancestors 'none'" in headers
    assert "/gallery/examples/*" in headers
    assert "! Content-Security-Policy" in headers
    assert "frame-ancestors 'self'" in headers
    assert "gallery/media/**/*" in cloudflare_module.GALLERY_REMOTE_ASSET_PATTERNS
    assert "gallery/sessions/*.gbdraw-session.json.gz" in cloudflare_module.GALLERY_REMOTE_ASSET_PATTERNS
    remote_assets = json.loads((bundle_path / "gallery" / "remote-assets.json").read_text(encoding="utf-8"))
    assert "gallery/examples/Vnig_TUMSAT-TG-2018.svg" not in remote_assets
    assert "gallery/sessions/Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz" not in remote_assets
    assert (
        remote_assets[
            "gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
        ]
        == f"{remote_base}gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
    )
    assert (
        remote_assets["gallery/examples/vibrio-harveyi-group-collinear.svg"]
        == f"{remote_base}gallery/examples/vibrio-harveyi-group-collinear.svg"
    )
    assert all("/main/" not in url for url in remote_assets.values())
    assert (bundle_path / "gallery" / "examples" / "Vnig_TUMSAT-TG-2018.svg").exists()
    assert not (
        bundle_path
        / "gallery"
        / "examples"
        / "vibrio-harveyi-group-collinear.svg"
    ).exists()
    assert (
        bundle_path
        / "gallery"
        / "sessions"
        / "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz"
    ).exists()
    assert not (
        bundle_path
        / "gallery"
        / "sessions"
        / "vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
    ).exists()
    assert (bundle_path / "gallery" / "examples" / "majanivirus_orthogroup.svg").exists()


def test_cloudflare_gallery_remote_base_rejects_mutable_refs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    monkeypatch.delenv(cloudflare_module.GALLERY_REMOTE_BASE_ENV, raising=False)
    monkeypatch.setenv(cloudflare_module.GALLERY_REMOTE_REF_ENV, "main")

    with pytest.raises(RuntimeError, match="full 40-character commit SHA"):
        cloudflare_module._default_gallery_remote_base()


def test_cloudflare_gallery_remote_base_uses_resolved_commit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    commit_sha = "1234567890abcdef1234567890abcdef12345678"
    monkeypatch.delenv(cloudflare_module.GALLERY_REMOTE_BASE_ENV, raising=False)
    monkeypatch.delenv(cloudflare_module.GALLERY_REMOTE_REF_ENV, raising=False)
    monkeypatch.setattr(
        cloudflare_module,
        "_load_stamp_web_build_module",
        lambda: SimpleNamespace(resolve_commit_sha=lambda: commit_sha),
    )

    assert cloudflare_module._default_gallery_remote_base() == (
        "https://raw.githubusercontent.com/satoshikawato/gbdraw/"
        f"{commit_sha}/gbdraw/web/"
    )


def test_cloudflare_prepare_refreshes_gallery_only_when_requested(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    calls: list[tuple[str, object]] = []
    output_root = tmp_path / "cloudflare-pages"

    monkeypatch.setattr(
        cloudflare_module,
        "_load_prepare_browser_wheel_module",
        lambda: SimpleNamespace(
            prepare_browser_wheel=lambda refresh_cache_bust=False: calls.append(
                ("wheel", refresh_cache_bust)
            )
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_build_support_module",
        lambda: SimpleNamespace(
            refresh_open_source_notices=lambda: calls.append(("notices", None))
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_refresh_gallery_sessions_module",
        lambda: pytest.fail("Gallery refresh should be opt-in for Cloudflare Pages."),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "build_cloudflare_pages_bundle",
        lambda *,
        output_root=cloudflare_module.DEFAULT_OUTPUT_ROOT,
        google_analytics_measurement_id=cloudflare_module.DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID,
        gallery_remote_base=None: calls.append(("bundle", output_root)) or output_root,
    )

    assert cloudflare_module.prepare_cloudflare_pages(output_root=output_root) == output_root
    assert calls == [("wheel", False), ("notices", None), ("bundle", output_root)]

    calls.clear()
    gallery_module = SimpleNamespace(
        refresh_gallery_sessions=lambda: calls.append(("gallery", "sessions")),
        prepare_gallery_assets=lambda: calls.append(("gallery", "assets")),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_refresh_gallery_sessions_module",
        lambda: gallery_module,
    )

    assert (
        cloudflare_module.prepare_cloudflare_pages(
            refresh_cache_bust=True,
            refresh_gallery_sessions=True,
            output_root=output_root,
        )
        == output_root
    )
    assert calls == [
        ("gallery", "sessions"),
        ("gallery", "assets"),
        ("wheel", True),
        ("notices", None),
        ("bundle", output_root),
    ]


def test_cloudflare_prepare_refreshes_open_source_notices_before_copy(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    calls: list[str] = []
    output_root = tmp_path / "cloudflare-pages"

    monkeypatch.setattr(
        cloudflare_module,
        "_load_prepare_browser_wheel_module",
        lambda: SimpleNamespace(
            prepare_browser_wheel=lambda refresh_cache_bust=False: calls.append("wheel")
        ),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "_load_build_support_module",
        lambda: SimpleNamespace(refresh_open_source_notices=lambda: calls.append("notices")),
    )
    monkeypatch.setattr(
        cloudflare_module,
        "build_cloudflare_pages_bundle",
        lambda *,
        output_root=cloudflare_module.DEFAULT_OUTPUT_ROOT,
        google_analytics_measurement_id=cloudflare_module.DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID,
        gallery_remote_base=None: calls.append("copy") or output_root,
    )

    assert cloudflare_module.prepare_cloudflare_pages(output_root=output_root) == output_root
    assert calls == [
        "wheel",
        "notices",
        "copy",
    ]


def test_wrangler_uses_cloudflare_bundle_directory() -> None:
    wrangler_toml = (REPO_ROOT / "wrangler.toml").read_text(encoding="utf-8")
    assert 'main = "./gbdraw/web/cloudflare-worker.js"' in wrangler_toml
    assert 'directory = "./dist/cloudflare-pages"' in wrangler_toml
    assert 'binding = "ASSETS"' in wrangler_toml
    assert 'not_found_handling = "single-page-application"' in wrangler_toml
    assert '"/gallery/*"' in wrangler_toml


def test_cloudflare_worker_proxies_remote_gallery_assets() -> None:
    source = (WEB_ROOT / "cloudflare-worker.js").read_text(encoding="utf-8")
    assert "/gallery/remote-assets.json" in source
    assert "/gallery/examples/" in source
    assert "/gallery/sessions/" in source
    assert "/gallery/media/" in source
    assert "video/mp4" in source
    assert "video/webm" in source
    assert "video/ogg" in source
    assert "image/webp" in source
    assert "env.ASSETS.fetch" in source
    assert "Content-Security-Policy" in source
    assert "Cross-Origin-Embedder-Policy" in source
    assert "Cross-Origin-Opener-Policy" in source
    assert "Cross-Origin-Resource-Policy" in source
    assert "isGalleryViewRoute" in source
    assert "'/gallery/'" in source
    assert "PALETTE_EXPLORER_ASSET_PATH = '/gallery/palettes/'" in source


def test_cloudflare_worker_routes_gallery_pages_and_static_assets() -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    subprocess.run(
        [node, "tests/web/cloudflare-worker.test.mjs"],
        check=True,
        cwd=REPO_ROOT,
    )


def test_project_docs_and_citation_metadata_include_preprint_doi() -> None:
    readme = README_PATH.read_text(encoding="utf-8")
    assert PREPRINT_DOI in readme
    assert "./gbdraw/web/assets/gbdraw-logo-title.png" in readme
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


def test_web_linear_run_ignores_hidden_circular_species_strain_args() -> None:
    source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    assert (
        "if (mode.value === 'circular') {\n"
        "        if (form.species) args.push('--species', form.species);\n"
        "        if (form.strain) args.push('--strain', form.strain);\n"
        "      }"
    ) in source
    assert "if (form.species) args.push('--species', form.species);\n      if (form.strain)" not in source


def test_web_feature_lookup_uses_stable_data_attribute_with_dom_id_fallback() -> None:
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    feature_dom_source = (WEB_ROOT / "js" / "app" / "feature-dom.js").read_text(encoding="utf-8")
    svg_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "svg-actions.js").read_text(encoding="utf-8")
    label_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "label-actions.js").read_text(encoding="utf-8")
    color_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "color-actions.js").read_text(encoding="utf-8")
    stroke_actions_source = (WEB_ROOT / "js" / "app" / "legend" / "stroke-actions.js").read_text(encoding="utf-8")
    svg_styles_source = (WEB_ROOT / "js" / "app" / "svg-styles.js").read_text(encoding="utf-8")
    orthogroups_source = (WEB_ROOT / "js" / "app" / "orthogroups.js").read_text(encoding="utf-8")
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")
    standalone_source = _standalone_interactivity_source()

    assert "'data-gbdraw-feature-id'" in state_source
    assert "FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id'" in feature_dom_source
    assert "FEATURE_PART_ATTRIBUTE = 'data-gbdraw-feature-part'" in feature_dom_source
    assert "normalizeFeatureIdentity" in svg_actions_source
    assert "FEATURE_SELECTOR = [" in feature_dom_source
    assert "`path[${FEATURE_ID_ATTRIBUTE}]`" in feature_dom_source
    assert "element?.getAttribute?.(FEATURE_ID_ATTRIBUTE)" in feature_dom_source
    assert "isFeatureFillTarget" in feature_dom_source
    assert "svg.querySelectorAll(FEATURE_SELECTOR)" in svg_actions_source
    assert "getFeatureIdentity(element)" in svg_actions_source
    assert "getFeatureHoverKey(getFeatureIdentity(relatedFeature))" in svg_actions_source
    assert "getFeatureIdentity(el)" in label_actions_source
    assert "getFeatureFillElements(svg, feat.svg_id)" in color_actions_source
    assert "getFeatureElements(svg, svgId)" in color_actions_source
    assert "getFeatureElements(svg, svgId)" in stroke_actions_source
    assert "getFeatureIdentity(path)" in svg_styles_source
    assert "getFeatureElements(svg, featureId)" in orthogroups_source
    assert "import { enrichSvgWithStandaloneInteractivity, stripEditorOnlyCursorStyles } from './standalone-interactivity.js';" in export_source
    assert "FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id'" in standalone_source
    assert "normalizeFeatureElementId" in standalone_source
    assert "function getElementFeatureId(element)" in standalone_source
    assert "var svgId = getElementFeatureId(featureElement);" in standalone_source
    assert "const id = getElementFeatureId(element);" in standalone_source


def test_web_linear_custom_track_slots_are_wired() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    module_source = (WEB_ROOT / "js" / "app" / "linear-track-slots.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")

    assert "linear_track_slots_enabled: false" in state_source
    assert "linear_track_slots_schema_version: 1" in state_source
    assert "createDefaultLinearTrackSlots" in state_source
    assert "<span class=\"truncate\">Custom Track Slots</span>" in index_html
    assert '@click="setLinearTrackSlotsEnabled(!adv.linear_track_slots_enabled)"' in index_html
    assert 'aria-controls="linear-custom-track-slots-panel"' in index_html
    assert "Linear Custom Track Slots" not in index_html
    assert "Linear Custom Track Slots" not in config_source
    assert 'v-model="adv.linear_track_slots_enabled"' not in index_html
    assert 'v-model.number="entry.slot.params.track_index"' in index_html
    assert "--linear_track_slot" in run_source
    assert "buildLinearTrackSlotSpec" in run_source
    assert "linearSlotNeedsDepth" in run_source
    assert "validateImportedLinearTrackSlots" in config_source
    assert "LINEAR_TRACK_SLOT_SCHEMA_VERSION = 1" in config_source
    assert "createLinearTrackSlotEditor" in module_source
    assert "linearTrackStackEntries" in app_setup_source
    assert "linearTrackSlotUsesPresetGeometry(entry.slot)" in index_html
    assert "linearTrackSlotUsesPresetGeometry: linearTrackSlotEditor.linearTrackSlotUsesPresetGeometry" in app_setup_source
    assert "args.push('--ruler_label_font_size', adv.scale_font_size);" in run_source
    assert "args.push('--scale_font_size', adv.scale_font_size);" in run_source


def test_web_wires_gc_content_percent_options() -> None:
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "gc_content_mode: 'deviation'" in state_source
    assert "gc_content_min_percent: 0" in state_source
    assert "gc_content_max_percent: 100" in state_source
    assert "gc_content_tick_interval: 20" in state_source
    assert "state.adv.gc_content_mode" in config_source
    assert "state.adv.gc_content_min_percent" in config_source
    assert "state.adv.gc_content_max_percent" in config_source
    assert "appendGcContentPercentArgs" in run_source
    assert "if (adv.gc_content_mode !== 'percent') return;" in run_source
    assert "args.push('--gc_content_mode', 'percent');" in run_source
    assert "args.push('--gc_content_min_percent', String(minPercent));" in run_source
    assert "args.push('--gc_content_max_percent', String(maxPercent));" in run_source
    assert "args.push('--gc_content_large_tick_interval', adv.gc_content_tick_interval);" in run_source
    assert "args.push('--hide_gc_content_axis');" in run_source
    assert "args.push('--hide_gc_content_ticks');" in run_source
    assert 'v-model="adv.gc_content_mode"' in index_html
    assert "GC Content Mode" in index_html
    assert "adv.gc_content_mode === 'percent'" in index_html
    assert "GC Min %" in index_html


def test_web_wires_addable_depth_tracks() -> None:
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    depth_tracks_source = (WEB_ROOT / "js" / "app" / "depth-tracks.js").read_text(encoding="utf-8")
    depth_state_source = (WEB_ROOT / "js" / "app" / "depth-track-state.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "depth_tracks: []" in state_source
    assert "normalizeDepthTracks(state.adv.depth_tracks, state.adv)" in config_source
    assert "getDepthTrackFileBaseName" in depth_tracks_source
    assert "getDepthTrackLabelFromFile" in depth_tracks_source
    assert "isDepthTrackAutoLabel" in depth_tracks_source
    assert "depthTrackRows" in app_setup_source
    assert "circularDepthTrackRows" in app_setup_source
    assert "linearDepthTrackRows" in app_setup_source
    assert "ensureDepthTrackConfigShape" in app_setup_source
    assert "ensureDepthTrackConfigShape" in depth_state_source
    rows_block = app_setup_source.split("const rowsForDepthTrackCount", 1)[1].split(
        "const linearDepthTrackUiCount",
        1,
    )[0]
    assert "ensureDepthTrackEditableConfigCount(normalizedCount);" in rows_block
    assert "ensureDepthTrackConfigCount(" not in rows_block
    assert "hasCircularDepthFiles" in app_setup_source
    assert "hasLinearDepthFiles" in app_setup_source
    assert "depthTrackCountLabel" in app_setup_source
    assert "getDepthTrackLegendLabelForSlot" in app_setup_source
    assert "setDepthTrackLegendLabelForSlot" in app_setup_source
    assert "syncDepthTrackSlotLabel" in app_setup_source
    assert "addCircularDepthTrack" in app_setup_source
    assert "addLinearDepthTrack" in app_setup_source
    assert "setCircularDepthFile" in app_setup_source
    assert "setLinearDepthFile" in app_setup_source
    assert "updateDepthTrackLabelFromFile(idx, file, previousFile);" in app_setup_source
    assert "Depth TSV tracks" in index_html
    assert 'v-if="!hasCircularDepthFiles"' in index_html
    assert 'v-if="!hasLinearDepthFiles(seq)"' in index_html
    assert "Add TSV" in index_html
    assert "Per-track settings" in index_html
    assert ':value="getDepthTrackLabel(track.index)"' in index_html
    assert '@input="setDepthTrackLabel(track.index, $event.target.value)"' in index_html
    assert ':value="getDepthTrackLegendLabelForSlot(entry.slot)"' in index_html
    assert '@input="setDepthTrackLegendLabelForSlot(entry.slot, $event.target.value)"' in index_html
    assert "v-model.number=\"track.config.height\"" in index_html
    assert "v-model.number=\"track.config.large_tick_interval\"" in index_html
    assert "v-model.number=\"track.config.small_tick_interval\"" in index_html
    assert "v-model.number=\"track.config.tick_font_size\"" in index_html
    assert ':value="optionalNumberInputValue(adv.plot_title_font_size)"' in index_html
    assert "@input=\"setOptionalNumberInputValue(adv, 'plot_title_font_size', $event.target.value)\"" in index_html
    assert ':value="optionalNumberInputValue(adv.def_font_size)"' in index_html
    assert "@input=\"setOptionalNumberInputValue(adv, 'def_font_size', $event.target.value)\"" in index_html
    definition_size_block = app_setup_source.split("const setDefinitionLineStyleSize", 1)[1].split(
        "const getDefinitionLineStyleWeight",
        1,
    )[0]
    assert "setOptionalNumberInputValue" in definition_size_block
    assert "Number(value)" not in definition_size_block
    assert "depthTrackConfigAt(index, file)" in run_source
    assert "syncDepthSlotLegendLabelsFromTrackConfigs" in run_source
    assert "syncDepthSlotLabels({" in run_source
    assert "slot.params.legend_label = label;" in depth_state_source
    assert "args.push('--depth_track_label', ...labels);" in run_source
    assert "args.push('--depth_track_color', ...colors);" in run_source
    assert "args.push('--depth_track_height', ...heights);" in run_source
    assert "args.push('--depth_track_large_tick_interval', ...largeTicks);" in run_source
    assert "args.push('--depth_track_small_tick_interval', ...smallTicks);" in run_source
    assert "args.push('--depth_track_tick_font_size', ...tickFontSizes);" in run_source
    assert "depthPaths.push('');" in run_source


def test_web_run_analysis_wires_circular_track_slot_options() -> None:
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    slot_source = (WEB_ROOT / "js" / "app" / "circular-track-slots.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "circular_track_slots_enabled" in state_source
    assert "circular_track_slots_axis_index" in state_source
    assert "center_reserved_radius: null" in state_source
    assert "createDefaultCircularTrackSlots()" in state_source
    assert "state.adv.center_reserved_radius = normalizeNonNegativeNumberOrNull(state.adv.center_reserved_radius);" in config_source
    assert "inferLegacyAxisIndexFromFeature(normalizedSlots, state.form.track_type)" in config_source
    assert '"circular_track_slot": "--circular_track_slot" in _source' in run_source
    assert '"circular_track_axis_index": "--circular_track_axis_index" in _source' in run_source
    assert '"center_reserved_radius": "--center_reserved_radius" in _source' in run_source
    assert "args.push('--track_type', form.track_type);" in run_source
    assert "args.push('--center_reserved_radius', String(normalizedCenterReservedRadius));" in run_source
    assert "applyCircularSuppressControlsToSlots" in run_source
    assert "if (form.suppress_gc) args.push('--suppress_gc');" in run_source
    assert "if (form.suppress_skew) args.push('--suppress_skew');" in run_source
    assert "args.push('--circular_track_axis_index', String(adv.circular_track_slots_axis_index));" in run_source
    assert "buildCircularTrackSlotSpec(slot, adv.nt, form.track_type, {" in run_source
    assert "applyCircularTrackOrderPlacements(" in run_source
    assert "if (useCircularTrackSlots)" in run_source
    assert "hasEnabledCircularTrackRenderer(circularTrackSlots, 'depth')" in run_source
    assert "Custom Track Slots" in index_html
    assert '@click="setCircularTrackSlotsEnabled(!adv.circular_track_slots_enabled)"' in index_html
    assert 'aria-controls="circular-custom-track-slots-panel"' in index_html
    assert "Track Preset" in index_html
    assert 'v-if="!adv.circular_track_slots_enabled"' in index_html
    assert "Reset to Tuckin" in index_html
    assert "Reset to Middle" in index_html
    assert "Reset to Spreadout" in index_html
    assert "Center Reserved Radius" in index_html
    assert "Apply Tuckin" not in index_html
    assert "arrows reorder within the current side" in index_html
    assert "Radial track stack" in index_html
    assert "circularTrackStackEntries()" in index_html
    assert "Use Move outside or Move inside to cross the Axis" in index_html
    assert "Move outside Axis" in index_html
    assert "Move inside Axis" in index_html
    assert "entry.onAxis" in index_html
    assert "Add track" in index_html
    assert "Outer tracks" not in slot_source
    assert "On-axis tracks" not in slot_source
    assert "Inner tracks" not in slot_source
    assert "circularTrackSlots" in slot_source
    assert "axisIndexForSlots" in slot_source
    assert "effectiveSlotPlacement" in slot_source
    assert "wouldCircularTrackSlotMoveCrossAxis" in slot_source
    assert "syncSlotPlacementFromSide" in slot_source
    assert "moveCircularTrackSlotOutside" in slot_source
    assert "moveCircularTrackSlotInside" in slot_source
    assert "moveCircularTrackSlotToAxis" in slot_source
    assert "Feature Layout" not in index_html
    assert "params.axis" not in slot_source
    assert "axis=true" not in slot_source
    assert "side = null" in slot_source
    assert "isLegacyDefaultWebSlotShape" in slot_source
    assert "ensureCircularTrackDepthSlot" in slot_source
    assert "setCircularGcSuppressed" in slot_source
    assert "setCircularSkewSuppressed" in slot_source
    assert "override the custom track settings" in slot_source
    assert "Replace the current custom circular track slots with this preset" not in slot_source
    assert "setCircularTrackSlotsEnabled" in slot_source
    assert "setCircularTrackSlotEnabled: circularTrackSlotEditor.setCircularTrackSlotEnabled" in app_setup_source
    assert "circularTrackSlotHiddenBySuppress: circularTrackSlotEditor.circularTrackSlotHiddenBySuppress" in app_setup_source
    assert "const templateSlots = createDefaultCircularTrackSlots" in slot_source
    assert "const suppressed = applyCircularSuppressControlsToSlots(normalized, state.form);" in slot_source
    assert "state.adv.circular_track_slots.splice(" in slot_source
    assert '@change="setCircularGcSuppressed($event.target.checked, $event)"' in index_html
    assert '@change="setCircularSkewSuppressed($event.target.checked, $event)"' in index_html
    assert "circularTrackSlotSuppressMessage(entry.slot)" in index_html
    assert "depthFileSlotsFromValue(files.c_depth).length" in app_setup_source
    assert "circularTrackSlotEditor.normalizeCircularTrackSlots();" in app_setup_source
    assert "circularTrackSlotEditor.ensureCircularTrackDepthSlot();" in app_setup_source
    assert "resetCircularTrackSlotsToPreset: circularTrackSlotEditor.resetCircularTrackSlotsToPreset" in app_setup_source
    assert "setCircularTrackSlotsEnabled: circularTrackSlotEditor.setCircularTrackSlotsEnabled" in app_setup_source
    assert "moveCircularTrackSlotOutside: circularTrackSlotEditor.moveCircularTrackSlotOutside" in app_setup_source
    assert "moveCircularTrackSlotInside: circularTrackSlotEditor.moveCircularTrackSlotInside" in app_setup_source
    assert "canMoveCircularTrackSlotOutside: circularTrackSlotEditor.canMoveCircularTrackSlotOutside" in app_setup_source
    assert "canMoveCircularTrackSlotInside: circularTrackSlotEditor.canMoveCircularTrackSlotInside" in app_setup_source


def test_web_wires_circular_conservation_options() -> None:
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    slot_source = (WEB_ROOT / "js" / "app" / "circular-track-slots.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    components_source = (WEB_ROOT / "js" / "components.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "const circularConservation = reactive" in state_source
    assert "c_conservation_blasts: []" in state_source
    assert "c_conservation_fastas: []" in state_source
    assert "losat_program: 'blastn'" in state_source
    assert "subject_gencode: 1" in state_source
    assert "losat_gencode: normalizePositiveInteger(entry.losat_gencode, 1)" in config_source
    assert "series: []" in state_source
    assert "'data-source-index'" in state_source
    assert "'data-reference-record-id'" in state_source
    assert "'data-track-color'" in state_source
    assert "circularConservation" in app_setup_source
    assert "circularConservationSeriesRows" in app_setup_source
    assert "syncCircularConservationSeries" in app_setup_source
    assert "syncCircularConservationEnabled" in app_setup_source
    assert "addCircularConservationComparisonFile" in app_setup_source
    assert "removeCircularConservationSource" in app_setup_source
    assert "circularConservation: state.circularConservation" in config_source
    assert "normalizeCircularConservationSeries" in config_source
    assert "normalizeCircularConservationLosatProgram" in config_source
    assert "c_conservation_blasts: await serializeFileArray" in config_source
    assert "normalizeCircularConservationReference" in config_source
    assert "Pairwise Comparisons" in index_html
    assert 'v-model="circularConservation.losat_program"' in index_html
    assert 'v-model.number="circularConservation.series[row.index].losat_gencode"' in index_html
    assert 'v-model.number="circularConservation.query_gencode"' not in index_html
    assert "Default subject gencode" not in index_html
    assert "comparisonEntry?.losat_gencode" in run_source
    assert "TLOSATX" in index_html
    assert "Conservation Rings" not in index_html
    assert 'v-model="circularConservation.enabled"' not in index_html
    assert "BLAST outfmt 6/7 files" in index_html
    assert "Comparison FASTA files" in index_html
    assert "@click=\"openCircularConservationComparisonFilePicker\"" in index_html
    assert "@click=\"removeCircularConservationSource(row.index)\"" in index_html
    assert "type=\"color\" v-model=\"circularConservation.series[row.index].color\"" in index_html
    assert ":multiple=\"true\"" in index_html
    assert "props: ['label', 'accept', 'modelValue', 'small', 'multiple']" in components_source
    assert '"conservation_blast": "--conservation_blast" in _source' in run_source
    assert '"records_table": "--records_table" in _source' in run_source
    assert '"conservation_table": "--conservation_table" in _source' in run_source
    assert '"circular_track_table": "--circular_track_table" in _source' in run_source
    assert '"conservation_colors": "--conservation_colors" in _source' in run_source
    assert "runCircularLosatConservation" in run_source
    assert "buildConservationSeries" in run_source
    assert "program: circularLosatProgram" in run_source
    assert "circularLosatProgram === 'tblastx'" in run_source
    assert "args.push('--conservation_colors', ...colors);" in run_source
    assert "args.push('--conservation_blast', ...conservationBlastPaths);" in run_source
    assert "args.push('--conservation_reference', conservationReference);" in run_source
    assert "flow: 'circular-conservation'" in run_source
    assert "conservationReference = 'subject';" in run_source
    assert "'sequence_conservation'" in slot_source


def test_web_collinear_orientation_identity_mode_is_wired() -> None:
    index_source = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    color_source = (WEB_ROOT / "js" / "app" / "color-utils.js").read_text(encoding="utf-8")
    svg_styles_source = (WEB_ROOT / "js" / "app" / "svg-styles.js").read_text(encoding="utf-8")

    assert '<option value="orientation_identity">Orientation + identity</option>' in index_source
    assert "['average_identity', 'orientation', 'orientation_identity']" in run_source
    assert "['average_identity', 'orientation', 'orientation_identity']" in config_source
    assert "normalizedMode === 'orientation_identity'" in color_source
    assert "collinear_block_plus_min" in color_source
    assert "collinear_block_minus_min" in color_source
    assert "identityFactor: Number.isFinite(metadataFactor) ? metadataFactor : null" in svg_styles_source


def test_web_collinear_orientation_identity_recoloring_uses_identity_factor(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "color-utils.js"
    module_path = tmp_path / "color-utils.mjs"
    module_path.write_text(source_path.read_text(encoding="utf-8"), encoding="utf-8")
    check_path = tmp_path / "check-color-utils.mjs"
    check_path.write_text(
        f"""
        import {{
          buildDefaultColorOverrideTsv,
          buildPaletteColorOverrideRows,
          interpolateColor,
          normalizePaletteColors,
          resolveCollinearMatchColor,
          resolvePairwiseLegendGradientColorKeys
        }} from {module_path.as_uri()!r};

        const colors = {{
          pairwise_match_min: '#ffffff',
          pairwise_match_max: '#000000',
          collinear_block_plus_min: '#eeeeee',
          collinear_block_plus: '#808080',
          collinear_block_minus_min: '#ffeeee',
          collinear_block_minus: '#ff0000'
        }};
        const graded = resolveCollinearMatchColor({{
          blockId: 'block_1',
          colorMode: 'orientation_identity',
          orientation: 'minus',
          identityFactor: 0.5,
          colors
        }});
        const expected = interpolateColor('#ffeeee', '#ff0000', 0.5);
        if (graded !== expected) {{
          throw new Error(`Expected orientation identity ramp ${{expected}}, got ${{graded}}`);
        }}
        const fixed = resolveCollinearMatchColor({{
          blockId: 'block_1',
          colorMode: 'orientation',
          orientation: 'minus',
          identityFactor: 0.5,
          colors
        }});
        if (fixed !== '#ff0000') {{
          throw new Error(`Orientation mode should keep a fixed endpoint color, got ${{fixed}}`);
        }}
        const averageIdentity = resolveCollinearMatchColor({{
          blockId: 'block_1',
          colorMode: 'average_identity',
          orientation: 'minus',
          identityFactor: 0.5,
          colors
        }});
        if (averageIdentity !== null) {{
          throw new Error(`Average identity mode should fall back to pairwise recoloring, got ${{averageIdentity}}`);
        }}
        const plusLegendKeys = resolvePairwiseLegendGradientColorKeys('Collinear');
        if (
          plusLegendKeys.minKey !== 'collinear_block_plus_min' ||
          plusLegendKeys.maxKey !== 'collinear_block_plus'
        ) {{
          throw new Error(`Collinear legend should use plus gradient keys, got ${{JSON.stringify(plusLegendKeys)}}`);
        }}
        const minusLegendKeys = resolvePairwiseLegendGradientColorKeys('Inverted');
        if (
          minusLegendKeys.minKey !== 'collinear_block_minus_min' ||
          minusLegendKeys.maxKey !== 'collinear_block_minus'
        ) {{
          throw new Error(`Inverted legend should use minus gradient keys, got ${{JSON.stringify(minusLegendKeys)}}`);
        }}
        const aliasPalette = normalizePaletteColors({{ collinear_block_plus_max: '#123456' }});
        if (aliasPalette.collinear_block_plus !== '#123456') {{
          throw new Error(`collinear_block_plus_max should alias collinear_block_plus, got ${{aliasPalette.collinear_block_plus}}`);
        }}
        const paletteColors = normalizePaletteColors({{
          CDS: '#54bcf8',
          rRNA: '#71ee7d',
          pairwise_match_min: '#FFE7E7',
          pairwise_match_max: '#FF7272'
        }});
        const noOverrides = buildPaletteColorOverrideRows({{
          colors: {{
            CDS: '#54bcf8',
            rRNA: '#71EE7D',
            pairwise_match_min: '#ffe7e7',
            pairwise_match_max: '#ff7272'
          }},
          paletteColors
        }});
        if (noOverrides.length !== 0) {{
          throw new Error(`Palette-equivalent colors should not produce -d overrides, got ${{JSON.stringify(noOverrides)}}`);
        }}
        const overrideTsv = buildDefaultColorOverrideTsv({{
          colors: {{ CDS: '#000000', rRNA: '#71ee7d', custom_feature: '#abcdef' }},
          paletteColors
        }});
        if (overrideTsv !== 'CDS\\t#000000\\ncustom_feature\\t#abcdef') {{
          throw new Error(`Expected only changed/default-missing color overrides, got ${{JSON.stringify(overrideTsv)}}`);
        }}
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)


def test_linear_track_slot_axis_sync_actions_and_specs(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "linear-track-slots.js"
    module_path = tmp_path / "linear-track-slots.mjs"
    (tmp_path / "package.json").write_text('{"type":"module"}', encoding="utf-8")
    for dependency in ["depth-track-state.js", "color-utils.js", "track-slot-colors.js", "track-slot-display.js"]:
        dep_path = WEB_ROOT / "js" / "app" / dependency
        (tmp_path / dependency).write_text(dep_path.read_text(encoding="utf-8"), encoding="utf-8")
    module_path.write_text(source_path.read_text(encoding="utf-8"), encoding="utf-8")
    check_path = tmp_path / "check-linear-track-slots.mjs"
    check_path.write_text(
        f"""
        import {{
          applyLinearTrackOrderPlacements,
          buildLinearTrackSlotSpec,
          createDefaultLinearTrackSlots,
          createLinearTrackSlotEditor
        }} from {module_path.as_uri()!r};

        const assert = (condition, message) => {{
          if (!condition) throw new Error(message);
        }};

        const defaultState = {{
          adv: {{
            nt: 'GC',
            linear_track_slots_enabled: true,
            linear_track_slots_axis_index: null,
            linear_track_slots: createDefaultLinearTrackSlots({{
              showDepth: true,
              depthTrackCount: 1,
              showGc: true,
              showSkew: true,
              trackLayout: 'middle'
            }})
          }},
          form: {{
            linear_track_layout: 'middle',
            show_depth: true,
            show_gc: true,
            show_skew: true
          }},
          linearSeqs: [{{ depth: [{{ name: 'depth.tsv' }}] }}]
        }};
        const defaultEditor = createLinearTrackSlotEditor({{ state: defaultState }});
        defaultEditor.normalizeLinearTrackSlots();
        const defaultEntries = defaultEditor.linearTrackStackEntries();
        assert(defaultState.adv.linear_track_slots_axis_index === 0, 'Middle default should place Axis at the feature row');
        assert(defaultEntries[0]?.kind === 'slot' && defaultEntries[0]?.onAxis === true, 'Middle default feature should embed the Axis row');
        assert(defaultState.adv.linear_track_slots[0].side === 'overlay', 'Middle default feature should use side=overlay');

        const state = {{
          adv: {{
            nt: 'GC',
            linear_track_slots_enabled: true,
            linear_track_slots_axis_index: 1,
            linear_track_slots: [
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'above', params: {{ nt: 'GC' }} }},
              {{ id: 'features', renderer: 'features', side: 'below', params: {{}} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'below', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            linear_track_layout: 'middle',
            show_depth: false,
            show_gc: true,
            show_skew: true
          }},
          linearSeqs: [],
          currentColors: {{ value: {{ skew_high: '#123456', skew_low: '#654321' }} }},
          paletteDefinitions: {{
            value: {{
              default: {{ skew_high: '#6dded3', skew_low: '#ad72e3' }},
              custom: {{ skew_high: '#abcdef', skew_low: '#fedcba' }}
            }}
          }},
          selectedPalette: {{ value: 'custom' }}
        }};
        const editor = createLinearTrackSlotEditor({{ state }});
        const stateSkewSlot = state.adv.linear_track_slots.find((slot) => slot.id === 'gc_skew');
        assert(editor.linearTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') === '#123456', 'Linear inherited positive skew color should use currentColors.skew_high');
        assert(editor.linearTrackSlotSkewColorValue(stateSkewSlot, 'negative_color') === '#654321', 'Linear inherited negative skew color should use currentColors.skew_low');
        state.currentColors.value.skew_high = '#abc';
        state.currentColors.value.skew_low = 'tomato';
        assert(editor.linearTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') === '#aabbcc', 'Linear inherited positive swatch should update with currentColors');
        assert(editor.linearTrackSlotSkewColorValue(stateSkewSlot, 'negative_color') === '#ff6347', 'Linear inherited negative swatch should normalize named currentColors');
        assert(!('positive_color' in stateSkewSlot.params) && !('negative_color' in stateSkewSlot.params), 'Inherited linear swatches should not add slot params');
        editor.setLinearTrackSlotSkewColor(stateSkewSlot, 'positive_color', '#010203');
        editor.setLinearTrackSlotSkewColor(stateSkewSlot, 'negative_color', '#040506');
        assert(editor.linearTrackSlotHasSkewColorOverride(stateSkewSlot, 'positive_color'), 'Linear positive override should be detected');
        assert(editor.linearTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') === '#010203', 'Linear explicit positive color should override inherited color');
        const explicitLinearSpec = editor.linearTrackSlotCliSpec(stateSkewSlot);
        assert(explicitLinearSpec.includes('positive_color=#010203') && explicitLinearSpec.includes('negative_color=#040506'), `Linear explicit skew colors were not serialized: ${{explicitLinearSpec}}`);
        editor.clearLinearTrackSlotSkewColor(stateSkewSlot, 'positive_color');
        editor.clearLinearTrackSlotSkewColor(stateSkewSlot, 'negative_color');
        assert(!editor.linearTrackSlotHasSkewColorOverride(stateSkewSlot, 'positive_color'), 'Linear clear should remove positive override');
        assert(editor.linearTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') === '#aabbcc', 'Linear clear should return positive swatch to inherited color');
        const clearedLinearSpec = editor.linearTrackSlotCliSpec(stateSkewSlot);
        assert(!clearedLinearSpec.includes('positive_color=') && !clearedLinearSpec.includes('negative_color='), `Linear cleared skew colors should not serialize params: ${{clearedLinearSpec}}`);
        assert(!editor.canMoveLinearTrackSlot(0, 1), 'Normal down arrow should not cross Axis');
        editor.moveLinearTrackSlot(0, 1);
        assert(state.adv.linear_track_slots.map((slot) => slot.id).join(',') === 'gc_content,features,gc_skew', 'Arrow move crossed Axis unexpectedly');

        editor.moveLinearTrackSlotToAxis(1);
        assert(state.adv.linear_track_slots_axis_index === 1, 'Feature on-Axis move should set axis index to the feature row');
        assert(state.adv.linear_track_slots[1].side === 'overlay', 'Feature on-Axis move should use side=overlay');
        assert(editor.linearTrackStackEntries()[1]?.onAxis === true, 'On-Axis feature should be embedded in stack entries');

        editor.moveLinearTrackSlotAbove(2);
        assert(state.adv.linear_track_slots.map((slot) => slot.id).join(',') === 'gc_content,gc_skew,features', 'Move above Axis inserted in the wrong position');
        assert(state.adv.linear_track_slots_axis_index === 2, 'Move above Axis should shift the boundary after the moved slot');
        assert(state.adv.linear_track_slots.map((slot) => slot.side).join(',') === 'above,above,overlay', 'Move above Axis did not sync sides');

        editor.moveLinearTrackSlotBelow(0);
        assert(state.adv.linear_track_slots.map((slot) => slot.id).join(',') === 'gc_skew,features,gc_content', 'Move below Axis inserted in the wrong position');
        assert(state.adv.linear_track_slots_axis_index === 1, 'Move below Axis should keep the on-Axis feature boundary');
        assert(state.adv.linear_track_slots.map((slot) => slot.side).join(',') === 'above,overlay,below', 'Move below Axis did not sync sides');

        const feature = state.adv.linear_track_slots.find((slot) => slot.id === 'features');
        editor.updateLinearTrackSlotPlacement(feature, 'below');
        let movedFeature = state.adv.linear_track_slots.find((slot) => slot.id === 'features');
        assert(movedFeature?.side === 'below' && state.adv.linear_track_slots_axis_index === 1, 'Feature below-axis selection should move the row, not just mutate side');
        editor.updateLinearTrackSlotPlacement(movedFeature, 'above');
        movedFeature = state.adv.linear_track_slots.find((slot) => slot.id === 'features');
        assert(movedFeature?.side === 'above' && state.adv.linear_track_slots_axis_index === 2, 'Feature above-axis selection should move above the boundary');
        editor.updateLinearTrackSlotPlacement(movedFeature, 'overlay');
        movedFeature = state.adv.linear_track_slots.find((slot) => slot.id === 'features');
        assert(movedFeature?.side === 'overlay' && state.adv.linear_track_slots_axis_index === 1, 'Feature on-axis selection should restore overlay placement');

        const duplicateAxis = applyLinearTrackOrderPlacements(
          [
            {{ id: 'features_a', renderer: 'features', side: 'overlay', params: {{}} }},
            {{ id: 'features_b', renderer: 'features', side: 'overlay', params: {{}} }},
            {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'below', params: {{ nt: 'GC' }} }}
          ],
          1
        );
        const onAxisIds = duplicateAxis.filter((slot) => slot.side === 'overlay').map((slot) => slot.id).join(',');
        assert(onAxisIds === 'features_b', 'Only the feature at the Axis index should remain on-axis');
        assert(duplicateAxis[0].side === 'above', 'Demoted duplicate on-axis feature should match row placement');

        const gcSpec = buildLinearTrackSlotSpec(
          {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'above', params: {{ nt: 'GC' }} }},
          {{ includeSide: false }}
        );
        assert(!gcSpec.includes('side='), 'Axis-derived above/below side should be omitted from CLI spec');
        const defaultSkewSpec = buildLinearTrackSlotSpec(
          {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'below', params: {{ nt: 'GC' }} }},
          {{ includeSide: false }}
        );
        assert(defaultSkewSpec === 'gc_skew:dinucleotide_skew@nt=GC', `Default skew slot should serialize exactly as before: ${{defaultSkewSpec}}`);
        const coloredSkewSpec = buildLinearTrackSlotSpec(
          {{ id: 'at_skew', renderer: 'dinucleotide_skew', side: 'below', params: {{ nt: 'AT', high_color: 'tomato', low_color: '#2a9d8f' }} }},
          {{ includeSide: false }}
        );
        assert(coloredSkewSpec.includes('positive_color=#FF6347'), `Skew alias high_color was not canonicalized: ${{coloredSkewSpec}}`);
        assert(coloredSkewSpec.includes('negative_color=#2a9d8f'), `Skew alias low_color was not canonicalized: ${{coloredSkewSpec}}`);
        const featureSpec = buildLinearTrackSlotSpec(
          {{ id: 'features', renderer: 'features', side: 'overlay', params: {{}} }},
          {{ includeSide: false }}
        );
        assert(featureSpec.includes('side=overlay'), 'On-Axis feature CLI spec should keep side=overlay');

        const reconcileState = {{
          adv: {{
            nt: 'AT',
            linear_track_slots_enabled: true,
            linear_track_slots_axis_index: 0,
            linear_track_slots: [
              {{ id: 'features', renderer: 'features', side: 'overlay', params: {{}} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'below', params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'below', params: {{ nt: 'GC', positive_color: '#e76f51' }} }},
              {{ id: 'custom_gc', renderer: 'dinucleotide_content', side: 'below', height: '22px', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            linear_track_layout: 'middle',
            show_depth: false,
            show_gc: false,
            show_skew: false
          }},
          linearSeqs: [{{ depth: [{{ name: 'one.tsv' }}, {{ name: 'two.tsv' }}] }}]
        }};
        const reconcileEditor = createLinearTrackSlotEditor({{ state: reconcileState }});
        reconcileEditor.syncLinearNumericSlotsFromSimpleControls();
        assert(!reconcileState.adv.linear_track_slots.some((slot) => slot.id === 'gc_content'), 'Show GC off should remove only the default GC slot');
        assert(reconcileState.adv.linear_track_slots.some((slot) => slot.id === 'custom_gc'), 'Show GC off should preserve customized duplicate GC slots');
        assert(reconcileState.adv.linear_track_slots.some((slot) => slot.id === 'gc_skew' && slot.params?.positive_color === '#e76f51'), 'Show Skew off should preserve color-overridden gc_skew as custom');
        reconcileState.form.show_gc = true;
        reconcileState.form.show_skew = true;
        reconcileEditor.syncLinearNumericSlotsFromSimpleControls();
        const restoredGc = reconcileState.adv.linear_track_slots.find((slot) => slot.id === 'gc_content');
        assert(restoredGc?.params?.nt === 'AT', 'Show GC on should restore default GC using the simple nt control');
        assert(reconcileState.adv.linear_track_slots.filter((slot) => slot.id === 'gc_skew').length === 1, 'Show Skew on should not add a duplicate default over color-overridden gc_skew');
        reconcileState.form.show_depth = true;
        reconcileEditor.syncLinearNumericSlotsFromSimpleControls();
        const depthIndexes = reconcileState.adv.linear_track_slots
          .filter((slot) => slot.renderer === 'depth')
          .map((slot) => Number(slot.params?.track_index))
          .join(',');
        assert(depthIndexes === '0,1', 'Show Depth on should materialize one slot per uploaded depth series');
        reconcileState.form.linear_track_layout = 'below';
        reconcileEditor.applyLinearTrackLayoutPreset('below');
        const reconciledFeature = reconcileState.adv.linear_track_slots.find((slot) => slot.id === 'features');
        assert(reconciledFeature?.side === 'below' && reconcileState.adv.linear_track_slots_axis_index === 0, 'Track Layout Below should move feature below Axis');
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)


def test_circular_track_slot_axis_crossing_actions_keep_neighbor_sides(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "circular-track-slots.js"
    module_path = tmp_path / "circular-track-slots.mjs"
    (tmp_path / "package.json").write_text('{"type":"module"}', encoding="utf-8")
    for dependency in ["conservation-series.js", "color-utils.js", "depth-track-state.js", "track-slot-display.js", "track-slot-colors.js"]:
        dep_path = WEB_ROOT / "js" / "app" / dependency
        (tmp_path / dependency).write_text(dep_path.read_text(encoding="utf-8"), encoding="utf-8")
    module_path.write_text(source_path.read_text(encoding="utf-8"), encoding="utf-8")
    check_path = tmp_path / "check-circular-track-slots.mjs"
    check_path.write_text(
        f"""
        import {{
          applyCircularSuppressControlsToSlots,
          applyCircularTrackOrderPlacements,
          buildCircularTrackSlotSpec,
          createDefaultCircularTrackSlots,
          createCircularTrackSlotEditor,
          estimateCircularConservationLayoutWarning
        }} from {module_path.as_uri()!r};
        import {{ normalizeConservationSeriesColor }} from './conservation-series.js';
        import {{ formatCircularWidthValue }} from './track-slot-display.js';

        const defaultSlots = createDefaultCircularTrackSlots({{ preset: 'tuckin' }});
        const defaultTick = defaultSlots.find((slot) => slot.id === 'ticks');
        if (defaultTick?.params?.tick_label_layout !== 'label_in_tick_out') {{
          throw new Error(`Default Tick layout should point labels inward when Tick is inside Feature: ${{JSON.stringify(defaultTick)}}`);
        }}
        const gapSpec = buildCircularTrackSlotSpec(
          {{
            id: 'gc_content',
            renderer: 'dinucleotide_content',
            side: 'inside',
            inner_gap_px: '4',
            outer_gap_px: '6',
            params: {{ nt: 'GC' }}
          }},
          'GC',
          'tuckin',
          {{ includeSide: false }}
        );
        if (!gapSpec.includes('inner_gap_px=4') || !gapSpec.includes('outer_gap_px=6') || gapSpec.includes('spacing=')) {{
          throw new Error(`Circular gap spec did not serialize physical px gaps: ${{gapSpec}}`);
        }}
        const pixelWidthSpec = buildCircularTrackSlotSpec(
          {{
            id: 'gc_content',
            renderer: 'dinucleotide_content',
            side: 'inside',
            width: '50px',
            params: {{ nt: 'GC' }}
          }},
          'GC',
          'tuckin',
          {{ includeSide: false }}
        );
        if (!pixelWidthSpec.includes('w=50px')) {{
          throw new Error(`Circular px width was not preserved: ${{pixelWidthSpec}}`);
        }}
        const factorWidthSpec = buildCircularTrackSlotSpec(
          {{
            id: 'gc_content',
            renderer: 'dinucleotide_content',
            side: 'inside',
            width: '0.15',
            params: {{ nt: 'GC' }}
          }},
          'GC',
          'tuckin',
          {{ includeSide: false }}
        );
        if (!factorWidthSpec.includes('w=0.15') || factorWidthSpec.includes('w=0.15px')) {{
          throw new Error(`Circular unitless width factor was not preserved: ${{factorWidthSpec}}`);
        }}
        if (formatCircularWidthValue('50px') !== '50 px' || formatCircularWidthValue('0.15') !== '0.15 R') {{
          throw new Error('Circular width display did not distinguish px and radius factors');
        }}
        const legacySpacingSpec = buildCircularTrackSlotSpec(
          {{
            id: 'gc_skew',
            renderer: 'dinucleotide_skew',
            side: 'inside',
            spacing: '5',
            params: {{ nt: 'GC' }}
          }},
          'GC',
          'tuckin',
          {{ includeSide: false }}
        );
        if (!legacySpacingSpec.includes('inner_gap_px=5') || !legacySpacingSpec.includes('outer_gap_px=5') || legacySpacingSpec.includes('spacing=')) {{
          throw new Error(`Legacy circular spacing was not converted to physical gaps: ${{legacySpacingSpec}}`);
        }}
        const defaultSkewSpec = buildCircularTrackSlotSpec(
          {{
            id: 'gc_skew',
            renderer: 'dinucleotide_skew',
            side: 'inside',
            params: {{ nt: 'GC' }}
          }},
          'GC',
          'tuckin',
          {{ includeSide: false }}
        );
        if (defaultSkewSpec !== 'gc_skew:dinucleotide_skew') {{
          throw new Error(`Default circular skew slot should serialize exactly as before: ${{defaultSkewSpec}}`);
        }}
        const coloredSkewSpec = buildCircularTrackSlotSpec(
          {{
            id: 'at_skew',
            renderer: 'dinucleotide_skew',
            side: 'inside',
            params: {{ nt: 'AT', high_color: 'tomato', low_color: '#2a9d8f' }}
          }},
          'GC',
          'tuckin',
          {{ includeSide: false }}
        );
        if (!coloredSkewSpec.includes('positive_color=#FF6347') || !coloredSkewSpec.includes('negative_color=#2a9d8f')) {{
          throw new Error(`Circular skew alias colors were not canonicalized: ${{coloredSkewSpec}}`);
        }}

        const state = {{
          adv: {{
            nt: 'GC',
            circular_track_slots: [
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'outside', params: {{ nt: 'GC' }} }},
              {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }},
          currentColors: {{ value: {{ skew_high: '#0f1e2d', skew_low: '#2d1e0f' }} }},
          paletteDefinitions: {{
            value: {{
              default: {{ skew_high: '#6dded3', skew_low: '#ad72e3' }},
              custom: {{ skew_high: '#102030', skew_low: '#405060' }}
            }}
          }},
          selectedPalette: {{ value: 'custom' }}
        }};

        const defaultState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_enabled: false,
            circular_track_slots_axis_index: null,
            circular_track_slots: createDefaultCircularTrackSlots({{ preset: 'tuckin' }})
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};
        const defaultEditor = createCircularTrackSlotEditor({{ state: defaultState }});
        defaultEditor.setCircularTrackSlotsEnabled(true);
        const resetTick = defaultState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        if (resetTick?.side !== 'inside' || resetTick?.params?.tick_label_layout !== 'label_in_tick_out') {{
          throw new Error(`Enabling custom slots did not use the inward default Tick layout: ${{JSON.stringify(defaultState.adv.circular_track_slots)}}`);
        }}
        defaultEditor.resetCircularTrackSlotsToPreset('middle');
        const middleFeature = defaultState.adv.circular_track_slots.find((slot) => slot.id === 'features');
        const middleFeatureSpec = buildCircularTrackSlotSpec(
          middleFeature,
          defaultState.adv.nt,
          defaultState.form.track_type,
          {{ includeSide: false, forceSplitLane: true }}
        );
        if (defaultState.form.track_type !== 'middle' || middleFeature?.side !== 'overlay' || middleFeature?.params?.lane_direction !== 'split') {{
          throw new Error(`Reset to Middle did not put Feature on the Axis: ${{JSON.stringify(defaultState.adv.circular_track_slots)}}`);
        }}
        if (middleFeatureSpec !== 'features:features@lane_direction=split') {{
          throw new Error(`Middle feature CLI spec must keep lane_direction=split with circular axis index: ${{middleFeatureSpec}}`);
        }}

        const suppressState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_enabled: true,
            circular_track_slots_axis_index: 0,
            circular_track_slots: [
              {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', enabled: true, params: {{ nt: 'GC' }} }},
              {{ id: 'manual_gc', renderer: 'dinucleotide_content', side: 'inside', enabled: false, params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', enabled: true, params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};
        const suppressEditor = createCircularTrackSlotEditor({{ state: suppressState }});
        let confirmCalls = 0;
        globalThis.confirm = (message) => {{
          confirmCalls += 1;
          if (!String(message).includes('override the custom track settings')) {{
            throw new Error(`Suppress confirmation did not describe custom override: ${{message}}`);
          }}
          return false;
        }};
        const cancelEvent = {{ target: {{ checked: true }} }};
        suppressEditor.setCircularGcSuppressed(true, cancelEvent);
        if (confirmCalls !== 1 || suppressState.form.suppress_gc !== false || cancelEvent.target.checked !== false) {{
          throw new Error(`Cancelled suppress did not restore checkbox/form state: calls=${{confirmCalls}} form=${{suppressState.form.suppress_gc}} event=${{cancelEvent.target.checked}}`);
        }}
        const uncancelledGc = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'gc_content');
        if (uncancelledGc?.enabled !== true) {{
          throw new Error(`Cancelled suppress should not disable custom GC slot: ${{JSON.stringify(suppressState.adv.circular_track_slots)}}`);
        }}

        globalThis.confirm = () => true;
        const proceedEvent = {{ target: {{ checked: true }} }};
        suppressEditor.setCircularGcSuppressed(true, proceedEvent);
        const suppressedGc = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'gc_content');
        const manualGc = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'manual_gc');
        const activeSkew = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'gc_skew');
        if (suppressState.form.suppress_gc !== true || proceedEvent.target.checked !== true) {{
          throw new Error('Confirmed suppress did not set form/checkbox state');
        }}
        if (suppressedGc?.enabled !== false || suppressedGc?.params?._suppressed_by_global !== 'gc_content') {{
          throw new Error(`Confirmed suppress did not mark enabled GC slot as globally hidden: ${{JSON.stringify(suppressedGc)}}`);
        }}
        if (manualGc?.enabled !== false || manualGc?.params?._suppressed_by_global) {{
          throw new Error(`Manually disabled GC slot should stay disabled without global marker: ${{JSON.stringify(manualGc)}}`);
        }}
        if (activeSkew?.enabled !== true || suppressEditor.circularTrackSlotEffectiveEnabled(suppressedGc)) {{
          throw new Error(`Suppress GC should not affect skew and hidden GC should not be effectively enabled: ${{JSON.stringify(suppressState.adv.circular_track_slots)}}`);
        }}
        if (!suppressEditor.circularTrackSlotSuppressMessage(suppressedGc).includes('Hide GC Content')) {{
          throw new Error(`Suppress message did not name the controlling checkbox: ${{suppressEditor.circularTrackSlotSuppressMessage(suppressedGc)}}`);
        }}
        suppressEditor.normalizeCircularTrackSlots();
        const normalizedManualGc = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'manual_gc');
        if (normalizedManualGc?.params?._suppressed_by_global) {{
          throw new Error(`Re-normalizing while hidden should not mark manually disabled GC slots: ${{JSON.stringify(normalizedManualGc)}}`);
        }}
        const forcedSkewSlots = applyCircularSuppressControlsToSlots(
          [{{ id: 'gc_skew', renderer: 'dinucleotide_skew', enabled: true, params: {{ nt: 'GC' }} }}],
          {{ suppress_skew: true }}
        );
        if (forcedSkewSlots[0]?.enabled !== false || forcedSkewSlots[0]?.params?._suppressed_by_global !== 'gc_skew') {{
          throw new Error(`Run-time suppress guard did not disable skew slot: ${{JSON.stringify(forcedSkewSlots)}}`);
        }}
        suppressEditor.setCircularGcSuppressed(false, {{ target: {{ checked: false }} }});
        const restoredGc = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'gc_content');
        const stillManualGc = suppressState.adv.circular_track_slots.find((slot) => slot.id === 'manual_gc');
        if (suppressState.form.suppress_gc !== false || restoredGc?.enabled !== true || restoredGc?.params?._suppressed_by_global) {{
          throw new Error(`Unhiding GC did not restore globally hidden slot: ${{JSON.stringify(restoredGc)}}`);
        }}
        if (stillManualGc?.enabled !== false) {{
          throw new Error(`Unhiding GC should not restore manually disabled duplicate: ${{JSON.stringify(stillManualGc)}}`);
        }}

        const multiDepthState = {{
          adv: {{
            nt: 'GC',
            depth_tracks: [{{ label: 'Mars' }}, {{ label: 'DRR271272' }}],
            circular_track_slots_enabled: false,
            circular_track_slots_axis_index: null,
            circular_track_slots: createDefaultCircularTrackSlots({{ preset: 'tuckin' }})
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: true,
            suppress_gc: false,
            suppress_skew: false
          }},
          files: {{
            c_depth: [
              {{ name: 'Mars.depth.tsv', size: 10, lastModified: 100 }},
              {{ name: 'DRR271272.depth.tsv', size: 20, lastModified: 200 }}
            ]
          }}
        }};
        const multiDepthEditor = createCircularTrackSlotEditor({{ state: multiDepthState }});
        multiDepthEditor.setCircularTrackSlotsEnabled(true);
        const multiDepthSlots = multiDepthState.adv.circular_track_slots.filter((slot) => slot.renderer === 'depth');
        const multiDepthIndexes = multiDepthSlots.map((slot) => Number(slot.params?.track_index ?? 0)).join(',');
        if (multiDepthSlots.length !== 2 || multiDepthIndexes !== '0,1') {{
          throw new Error(`Multiple circular depth tracks were not materialized as slots: ${{JSON.stringify(multiDepthState.adv.circular_track_slots)}}`);
        }}
        const secondDepthSpec = multiDepthEditor.circularTrackSlotCliSpec(multiDepthSlots[1]);
        if (!secondDepthSpec.includes('track_index=1')) {{
          throw new Error(`Second depth slot CLI spec did not pin track_index=1: ${{secondDepthSpec}}`);
        }}
        multiDepthState.files.c_depth.push({{ name: 'third.depth.tsv', size: 30, lastModified: 300 }});
        multiDepthEditor.ensureCircularTrackDepthSlot();
        const expandedDepthSlots = multiDepthState.adv.circular_track_slots.filter((slot) => slot.renderer === 'depth');
        const expandedDepthIndexes = expandedDepthSlots.map((slot) => Number(slot.params?.track_index ?? 0)).join(',');
        if (expandedDepthSlots.length !== 3 || expandedDepthIndexes !== '0,1,2') {{
          throw new Error(`Adding a circular depth file did not add a custom depth slot: ${{JSON.stringify(multiDepthState.adv.circular_track_slots)}}`);
        }}

        const conservationState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_enabled: false,
            circular_track_slots_axis_index: null,
            circular_track_slots: createDefaultCircularTrackSlots({{ preset: 'tuckin' }})
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }},
          files: {{
            c_conservation_blasts: [],
            c_conservation_fastas: [
              {{ name: 'alpha.fa', size: 10, lastModified: 100 }},
              {{ name: 'beta.fa', size: 20, lastModified: 200 }}
            ]
          }},
          circularConservation: {{
            enabled: true,
            source: 'losat',
            series: [
              {{ sourceKey: 'beta.fa|20|200|0', fileName: 'beta.fa', label: 'Beta', color: '#f28e2b' }},
              {{ sourceKey: 'alpha.fa|10|100|0', fileName: 'alpha.fa', label: 'Alpha', color: '#4e79a7' }}
            ]
          }}
        }};
        const conservationEditor = createCircularTrackSlotEditor({{ state: conservationState }});
        conservationEditor.setCircularTrackSlotsEnabled(true);
        const conservationSlots = conservationState.adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation');
        if (conservationSlots.length !== 2 || conservationSlots[0].params.label !== 'Beta' || conservationSlots[1].params.label !== 'Alpha') {{
          throw new Error(`Conservation slots were not materialized in series order: ${{JSON.stringify(conservationState.adv.circular_track_slots)}}`);
        }}
        const betaSpec = conservationEditor.circularTrackSlotCliSpec(conservationSlots[0]);
        if (!betaSpec.includes('track_index=1') || !betaSpec.includes('source_index=0')) {{
          throw new Error(`Conservation slot CLI spec did not pin the ordered source: ${{betaSpec}}`);
        }}
        if (betaSpec.includes('w=')) {{
          throw new Error(`Default conservation slot should keep auto width, got: ${{betaSpec}}`);
        }}
        if (conservationEditor.circularTrackSlotDisplayLabel(conservationSlots[0]) !== 'Beta') {{
          throw new Error(`Conservation slot display label was not series label: ${{conservationEditor.circularTrackSlotDisplayLabel(conservationSlots[0])}}`);
        }}
        const warningState = {{
          mode: {{ value: 'circular' }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }},
          files: {{
            c_conservation_blasts: [],
            c_conservation_fastas: Array.from({{ length: 6 }}, (_, idx) => ({{
              name: `comparison_${{idx + 1}}.fa`,
              size: 10 + idx,
              lastModified: 100 + idx
            }}))
          }},
          circularConservation: {{
            enabled: true,
            source: 'losat',
            series: []
          }}
        }};
        const conservationWarning = estimateCircularConservationLayoutWarning(warningState);
        if (!conservationWarning.includes('auto-compress')) {{
          throw new Error(`Expected conservation layout warning, got: ${{conservationWarning}}`);
        }}
        const softenedDefaultColor = normalizeConservationSeriesColor(null, 0);
        if (softenedDefaultColor !== '#6e91b7') {{
          throw new Error(`Expected softened default conservation color #6e91b7, got: ${{softenedDefaultColor}}`);
        }}
        const explicitColor = normalizeConservationSeriesColor('#4e79a7', 0);
        if (explicitColor !== '#4e79a7') {{
          throw new Error(`Explicit conservation colors should be preserved, got: ${{explicitColor}}`);
        }}

        const outerTickSlots = applyCircularTrackOrderPlacements(
          [
            {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'label_in_tick_out' }} }},
            {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
            {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }}
          ],
          'GC',
          'tuckin',
          0
        );
        const outerTick = outerTickSlots.find((slot) => slot.id === 'ticks');
        if (outerTick?.params?.tick_label_layout !== 'label_out_tick_in') {{
          throw new Error(`Tick outside Feature did not flip labels outward: ${{JSON.stringify(outerTickSlots)}}`);
        }}

        const manualTickSlots = applyCircularTrackOrderPlacements(
          [
            {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
            {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
            {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }}
          ],
          'GC',
          'tuckin',
          0
        );
        const manualTick = manualTickSlots.find((slot) => slot.id === 'ticks');
        if (manualTick?.params?.tick_label_layout !== 'tick_only') {{
          throw new Error(`Manual Tick layout should not be auto-oriented: ${{JSON.stringify(manualTickSlots)}}`);
        }}

        const editor = createCircularTrackSlotEditor({{ state }});
        const stateSkewSlot = state.adv.circular_track_slots.find((slot) => slot.id === 'gc_skew');
        if (editor.circularTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') !== '#0f1e2d') {{
          throw new Error('Circular inherited positive skew color should use currentColors.skew_high');
        }}
        if (editor.circularTrackSlotSkewColorValue(stateSkewSlot, 'negative_color') !== '#2d1e0f') {{
          throw new Error('Circular inherited negative skew color should use currentColors.skew_low');
        }}
        state.currentColors.value = {{}};
        if (editor.circularTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') !== '#102030') {{
          throw new Error('Circular inherited positive skew color should fall back to the selected palette');
        }}
        state.selectedPalette.value = 'missing';
        if (editor.circularTrackSlotSkewColorValue(stateSkewSlot, 'negative_color') !== '#ad72e3') {{
          throw new Error('Circular inherited negative skew color should fall back to the default palette');
        }}
        editor.setCircularTrackSlotSkewColor(stateSkewSlot, 'positive_color', '#010203');
        if (!editor.circularTrackSlotHasSkewColorOverride(stateSkewSlot, 'positive_color')) {{
          throw new Error('Circular positive skew override should be detected');
        }}
        if (editor.circularTrackSlotSkewColorValue(stateSkewSlot, 'positive_color') !== '#010203') {{
          throw new Error('Circular explicit positive color should override inherited color');
        }}
        const explicitCircularSpec = editor.circularTrackSlotCliSpec(stateSkewSlot);
        if (!explicitCircularSpec.includes('positive_color=#010203')) {{
          throw new Error(`Circular explicit skew color was not serialized: ${{explicitCircularSpec}}`);
        }}
        editor.clearCircularTrackSlotSkewColor(stateSkewSlot, 'positive_color');
        if (editor.circularTrackSlotHasSkewColorOverride(stateSkewSlot, 'positive_color')) {{
          throw new Error('Circular clear should remove positive override');
        }}
        const clearedCircularSpec = editor.circularTrackSlotCliSpec(stateSkewSlot);
        if (clearedCircularSpec.includes('positive_color=') || clearedCircularSpec.includes('negative_color=')) {{
          throw new Error(`Circular cleared skew colors should not serialize params: ${{clearedCircularSpec}}`);
        }}
        if (editor.canMoveCircularTrackSlot(1, -1)) {{
          throw new Error('Feature up arrow should not cross Axis');
        }}

        editor.moveCircularTrackSlot(1, 0);
        if (state.adv.circular_track_slots[0].id !== 'gc_content') {{
          throw new Error('Arrow move crossed Axis unexpectedly');
        }}

        editor.moveCircularTrackSlotOutside(1);
        const ids = state.adv.circular_track_slots.map((slot) => slot.id).join(',');
        const sides = Object.fromEntries(state.adv.circular_track_slots.map((slot) => [slot.id, slot.side]));
        if (ids !== 'gc_content,features,ticks,gc_skew') {{
          throw new Error(`Unexpected slot order after Move outside: ${{ids}}`);
        }}
        if (sides.gc_content !== 'outside' || sides.features !== 'outside' || sides.ticks !== 'inside' || sides.gc_skew !== 'inside') {{
          throw new Error(`Unexpected sides after Move outside: ${{JSON.stringify(sides)}}`);
        }}
        if (state.adv.circular_track_slots[1].params.lane_direction !== 'outside') {{
          throw new Error('Feature lane_direction was not updated to outside');
        }}

        const tickState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots: [
              {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};

        const tickEditor = createCircularTrackSlotEditor({{ state: tickState }});
        tickEditor.moveCircularTrackSlotOutside(1);
        if (tickState.adv.circular_track_slots[0].id !== 'ticks' || tickState.adv.circular_track_slots[0].side !== 'outside') {{
          throw new Error(`Tick did not move outside: ${{JSON.stringify(tickState.adv.circular_track_slots)}}`);
        }}
        tickEditor.moveCircularTrackSlotInside(0);
        const tickIds = tickState.adv.circular_track_slots.map((slot) => slot.id).join(',');
        const tickSlot = tickState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        if (tickIds !== 'ticks,features,gc_content,gc_skew' || tickSlot?.side !== 'inside') {{
          throw new Error(`Tick did not return inside: ${{tickIds}} ${{JSON.stringify(tickSlot)}}`);
        }}
        if (tickSlot.params.tick_label_layout !== 'tick_only') {{
          throw new Error(`Tick layout did not remain tick_only: ${{JSON.stringify(tickSlot.params)}}`);
        }}

        const outsideFeatureState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots: [
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'outside', params: {{ nt: 'GC' }} }},
              {{ id: 'features', renderer: 'features', side: 'outside', params: {{ lane_direction: 'outside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};

        const outsideFeatureEditor = createCircularTrackSlotEditor({{ state: outsideFeatureState }});
        outsideFeatureEditor.moveCircularTrackSlotOutside(2);
        const outsideFeatureIds = outsideFeatureState.adv.circular_track_slots.map((slot) => slot.id).join(',');
        const outsideFeatureTick = outsideFeatureState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        if (outsideFeatureIds !== 'gc_skew,ticks,features,gc_content' || outsideFeatureTick?.side !== 'outside') {{
          throw new Error(`Tick did not move outside with outside Feature: ${{outsideFeatureIds}} ${{JSON.stringify(outsideFeatureTick)}}`);
        }}

        const laneSelectState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_axis_index: 0,
            circular_track_slots: [
              {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};

        const laneSelectEditor = createCircularTrackSlotEditor({{ state: laneSelectState }});
        laneSelectEditor.updateCircularTrackFeatureLane(laneSelectState.adv.circular_track_slots[0], 'outside');
        const movedFeature = laneSelectState.adv.circular_track_slots.find((slot) => slot.id === 'features');
        if (laneSelectState.adv.circular_track_slots_axis_index !== 1 || movedFeature?.side !== 'outside' || movedFeature?.params?.lane_direction !== 'outside') {{
          throw new Error(`Feature outside selection did not move outside Axis: ${{JSON.stringify(laneSelectState.adv)}}`);
        }}

        laneSelectEditor.updateCircularTrackFeatureLane(movedFeature, 'split');
        const stackEntries = laneSelectEditor.circularTrackStackEntries();
        if (laneSelectState.adv.circular_track_slots_axis_index !== 0 || stackEntries[0]?.kind !== 'slot' || stackEntries[0]?.onAxis !== true) {{
          throw new Error(`Feature on-axis selection was not embedded in Axis entry: ${{JSON.stringify(stackEntries)}}`);
        }}

        const tickAxisState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_axis_index: 0,
            circular_track_slots: [
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
              {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};

        const tickAxisEditor = createCircularTrackSlotEditor({{ state: tickAxisState }});
        tickAxisEditor.moveCircularTrackSlotToAxis(0);
        const axisTick = tickAxisState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        const tickStackEntries = tickAxisEditor.circularTrackStackEntries();
        if (axisTick?.side !== 'overlay' || axisTick?.params?.tick_label_layout !== 'tick_only') {{
          throw new Error(`Tick on-axis layout was not preserved: ${{JSON.stringify(axisTick)}}`);
        }}
        if (tickStackEntries[0]?.kind !== 'slot' || tickStackEntries[0]?.onAxis !== true) {{
          throw new Error(`Tick on-axis selection was not embedded in Axis entry: ${{JSON.stringify(tickStackEntries)}}`);
        }}
        if (!tickAxisEditor.canMoveCircularTrackSlotInside(0)) {{
          throw new Error('Tick on-axis slot should be movable inside Axis');
        }}
        axisTick.params.tick_label_layout = 'label_in_tick_out';
        const tickAxisSpec = tickAxisEditor.circularTrackSlotCliSpec(axisTick);
        if (!tickAxisSpec.includes('side=overlay') || !tickAxisSpec.includes('tick_label_layout=label_in_tick_out')) {{
          throw new Error(`Tick on-axis CLI spec lost inverted layout: ${{tickAxisSpec}}`);
        }}
        tickAxisEditor.moveCircularTrackSlotInside(0);
        const demotedAxisTick = tickAxisState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        const demotedTickStackEntries = tickAxisEditor.circularTrackStackEntries();
        if (
          tickAxisState.adv.circular_track_slots_axis_index !== 0 ||
          demotedAxisTick?.side !== 'inside' ||
          demotedTickStackEntries[0]?.kind !== 'axis' ||
          demotedTickStackEntries[1]?.slot?.id !== 'ticks'
        ) {{
          throw new Error(`Tick on-axis slot did not move inside Axis: ${{JSON.stringify(tickAxisState.adv)}} entries=${{JSON.stringify(demotedTickStackEntries)}}`);
        }}

        const axisSwapState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_axis_index: 0,
            circular_track_slots: [
              {{ id: 'features', renderer: 'features', side: 'overlay', params: {{ lane_direction: 'split' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ tick_label_layout: 'tick_only' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};

        const axisSwapEditor = createCircularTrackSlotEditor({{ state: axisSwapState }});
        axisSwapEditor.moveCircularTrackSlotToAxis(1);
        const axisSwapIds = axisSwapState.adv.circular_track_slots.map((slot) => slot.id).join(',');
        const axisSwapOnAxis = axisSwapState.adv.circular_track_slots.filter((slot) => slot.side === 'overlay').map((slot) => slot.id).join(',');
        const swappedFeature = axisSwapState.adv.circular_track_slots.find((slot) => slot.id === 'features');
        const swappedTick = axisSwapState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        if (axisSwapIds !== 'ticks,features,gc_content,gc_skew' || axisSwapOnAxis !== 'ticks') {{
          throw new Error(`Axis swap did not replace the existing on-axis slot: ${{axisSwapIds}} onAxis=${{axisSwapOnAxis}}`);
        }}
        if (axisSwapState.adv.circular_track_slots_axis_index !== 0 || swappedFeature?.side !== 'inside' || swappedFeature?.params?.lane_direction !== 'inside') {{
          throw new Error(`Previous on-axis feature was not demoted inside: ${{JSON.stringify(axisSwapState.adv)}}`);
        }}
        if (swappedTick?.params?.tick_label_layout !== 'tick_only') {{
          throw new Error(`New on-axis tick layout was not preserved after swap: ${{JSON.stringify(swappedTick)}}`);
        }}

        axisSwapEditor.updateCircularTrackFeatureLane(swappedFeature, 'split');
        const laneSwapIds = axisSwapState.adv.circular_track_slots.map((slot) => slot.id).join(',');
        const laneSwapOnAxis = axisSwapState.adv.circular_track_slots.filter((slot) => slot.side === 'overlay').map((slot) => slot.id).join(',');
        const laneSwappedTick = axisSwapState.adv.circular_track_slots.find((slot) => slot.id === 'ticks');
        if (laneSwapIds !== 'features,ticks,gc_content,gc_skew' || laneSwapOnAxis !== 'features') {{
          throw new Error(`Feature lane on-axis selection did not swap with existing tick: ${{laneSwapIds}} onAxis=${{laneSwapOnAxis}}`);
        }}
        if (laneSwappedTick?.side !== 'inside' || laneSwappedTick?.params?.tick_label_layout !== 'tick_only') {{
          throw new Error(`Previous on-axis tick was not demoted inside: ${{JSON.stringify(laneSwappedTick)}}`);
        }}

        const featureInsideAcrossAxisState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots_axis_index: 1,
            circular_track_slots: [
              {{ id: 'features', renderer: 'features', side: 'outside', params: {{ lane_direction: 'outside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'overlay', params: {{ tick_label_layout: 'label_out_tick_in' }} }},
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }},
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'inside', params: {{ nt: 'GC' }} }}
            ]
          }},
          form: {{
            track_type: 'tuckin',
            show_depth: false,
            suppress_gc: false,
            suppress_skew: false
          }}
        }};

        const featureInsideAcrossAxisEditor = createCircularTrackSlotEditor({{ state: featureInsideAcrossAxisState }});
        featureInsideAcrossAxisEditor.moveCircularTrackSlotInside(0);
        const featureInsideAcrossAxisIds = featureInsideAcrossAxisState.adv.circular_track_slots.map((slot) => slot.id).join(',');
        const featureInsideAcrossAxisFeature = featureInsideAcrossAxisState.adv.circular_track_slots.find((slot) => slot.id === 'features');
        const featureInsideAcrossAxisEntries = featureInsideAcrossAxisEditor.circularTrackStackEntries();
        if (
          featureInsideAcrossAxisIds !== 'ticks,features,gc_content,gc_skew' ||
          featureInsideAcrossAxisState.adv.circular_track_slots_axis_index !== 0 ||
          featureInsideAcrossAxisFeature?.side !== 'inside' ||
          featureInsideAcrossAxisFeature?.params?.lane_direction !== 'inside' ||
          featureInsideAcrossAxisEntries[0]?.slot?.id !== 'ticks' ||
          featureInsideAcrossAxisEntries[0]?.onAxis !== true ||
          featureInsideAcrossAxisEntries[1]?.slot?.id !== 'features'
        ) {{
          throw new Error(`Feature did not move inside across an on-axis tick: ${{JSON.stringify(featureInsideAcrossAxisState.adv)}} entries=${{JSON.stringify(featureInsideAcrossAxisEntries)}}`);
        }}
        const featureInsideAcrossAxisNormalized = applyCircularTrackOrderPlacements(
          featureInsideAcrossAxisState.adv.circular_track_slots,
          'GC',
          'tuckin',
          featureInsideAcrossAxisState.adv.circular_track_slots_axis_index
        );
        const normalizedFeatureInsideAcrossAxisFeature = featureInsideAcrossAxisNormalized.find((slot) => slot.id === 'features');
        if (normalizedFeatureInsideAcrossAxisFeature?.side !== 'inside' || normalizedFeatureInsideAcrossAxisFeature?.params?.lane_direction !== 'inside') {{
          throw new Error(`Feature inside placement did not survive final normalization: ${{JSON.stringify(featureInsideAcrossAxisNormalized)}}`);
        }}

        const duplicateAxisSlots = applyCircularTrackOrderPlacements(
          [
            {{ id: 'features', renderer: 'features', side: 'overlay', params: {{ lane_direction: 'split' }} }},
            {{ id: 'ticks', renderer: 'ticks', side: 'overlay', params: {{ tick_label_layout: 'label_out_tick_in' }} }},
            {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', params: {{ nt: 'GC' }} }}
          ],
          'GC',
          'tuckin',
          1
        );
        const normalizedOnAxis = duplicateAxisSlots.filter((slot) => slot.side === 'overlay').map((slot) => slot.id).join(',');
        const normalizedFeature = duplicateAxisSlots.find((slot) => slot.id === 'features');
        if (normalizedOnAxis !== 'ticks' || normalizedFeature?.side !== 'outside' || normalizedFeature?.params?.lane_direction !== 'outside') {{
          throw new Error(`Duplicate on-axis normalization did not keep only the preferred slot: ${{JSON.stringify(duplicateAxisSlots)}}`);
        }}
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)


def test_web_config_rejects_obsolete_circular_track_slot_import_shapes() -> None:
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")

    assert "circular_track_slots_schema_version: 4" in state_source
    assert "const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 4;" in config_source
    assert "const LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;" in config_source
    assert "adv.circular_track_slots_schema_version !== CIRCULAR_TRACK_SLOT_SCHEMA_VERSION" in config_source
    assert "validateImportedCircularTrackSlots(data);" in config_source
    assert "isLegacyConfigPayload(data)" in config_source
    assert "validateImportedCircularTrackSlots(data.config);" in config_source
    assert "Legacy configuration loaded. Save as a session to use the current format." in config_source
    assert "Failed to load session: ${message}" in config_source

    for obsolete_key in [
        "gapAfter",
        "gap_after",
        "innerRadius",
        "inner_radius",
        "outerRadius",
        "outer_radius",
        "placement",
    ]:
        assert f"'{obsolete_key}'" in config_source

    for obsolete_param_key in ["side", "radius", "width"]:
        assert f"'{obsolete_param_key}'" in config_source


def test_web_session_uses_structured_depth_file_codec(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "services" / "depth-file-codec.js"
    module_path = tmp_path / "depth-file-codec.mjs"
    module_path.write_text(source_path.read_text(encoding="utf-8"), encoding="utf-8")
    check_path = tmp_path / "check-depth-file-codec.mjs"
    check_path.write_text(
        f"""
        import {{ DEPTH_FILE_ENCODING, decodeDepthText, encodeDepthText }} from {module_path.as_uri()!r};

        const depthText = [
          'reference_name\\tposition\\tdepth',
          'refA\\t1\\t10',
          'refA\\t2\\t20.5',
          'refA\\t4\\t0',
          'refB\\t10\\t1e2',
          ''
        ].join('\\n');
        const payload = encodeDepthText(depthText);
        if (!payload) throw new Error('Depth codec did not encode a valid TSV.');
        if (payload.records.length !== 2) throw new Error(`Expected two reference records: ${{JSON.stringify(payload)}}`);
        if (payload.records[0].id !== 'refA') throw new Error(`Reference ID was not stored once per record: ${{JSON.stringify(payload.records[0])}}`);
        const refARuns = payload.records[0].runs;
        if (refARuns.length !== 2) throw new Error(`Gapped coordinates should create two runs: ${{JSON.stringify(refARuns)}}`);
        if (JSON.stringify(refARuns[0].slice(0, 3)) !== JSON.stringify([1, 1, 2])) {{
          throw new Error(`Unexpected start/step/count for first run: ${{JSON.stringify(refARuns[0])}}`);
        }}
        if (decodeDepthText(payload) !== depthText) throw new Error('Depth codec did not round-trip text.');
        if (encodeDepthText('refA\\t1\\t10\\textra\\n') !== null) throw new Error('Extra TSV columns should fall back to base64 storage.');

        const largeDepthText = Array.from({{ length: 5000 }}, (_, idx) =>
          `refA\\t${{idx + 1}}\\t${{(idx * 17) % 80}}`
        ).join('\\n') + '\\n';
        const largePayload = encodeDepthText(largeDepthText);
        const structuredJson = JSON.stringify({{ encoding: DEPTH_FILE_ENCODING, data: largePayload }});
        const base64Json = JSON.stringify({{ data: Buffer.from(largeDepthText).toString('base64') }});
        if (structuredJson.length >= base64Json.length * 0.5) {{
          throw new Error(`Structured depth codec was not compact enough: ${{structuredJson.length}} vs ${{base64Json.length}}`);
        }}
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)

    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    assert "depth: await serializeDepthFile(seq.depth)" in config_source
    assert "c_depth: await serializeDepthFile(state.files.c_depth)" in config_source
    assert "await downloadCompressedSession(sessionData, sessionFilename);" in config_source


def test_web_config_persists_manual_qualifier_priority_rules() -> None:
    source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    assert "qualifierPriorityRules: cloneQualifierPriorityRules(state.manualPriorityRules)" in source
    assert "replaceQualifierPriorityRules(data.qualifierPriorityRules)" in source
    assert "replaceQualifierPriorityRules(data.priorityRules)" in source


def test_web_session_json_is_single_visible_settings_workflow() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")

    assert "Save Config" not in index_html
    assert "Load Config" not in index_html
    assert 'ref="configInput"' not in index_html
    assert "exportConfig" not in app_setup_source
    assert "importConfig" not in app_setup_source
    assert "editorState: buildEditorStateData()" in config_source
    assert "applyEditorStateData(data.editorState)" in config_source
    assert "featureStrokeOverrides" in config_source


def test_pyodide_palette_init_treats_comparison_only_colors_as_uninitialized() -> None:
    source = (WEB_ROOT / "js" / "app" / "pyodide.js").read_text(encoding="utf-8")

    assert "const currentHasPaletteColors = hasPaletteColorEntries(currentColors.value);" in source
    assert "if (!currentHasPaletteColors) {\n            currentColors.value = resolvedColors;" in source
    assert "currentHasPaletteColors ? currentColors.value : resolvedColors" in source
    assert "const currentHasColors = hasColorEntries(currentColors.value);" not in source


def test_conda_build_prepares_browser_wheel_before_install() -> None:
    build_sh = (REPO_ROOT / "recipe" / "build.sh").read_text(encoding="utf-8")
    meta_yaml = (REPO_ROOT / "recipe" / "meta.yaml").read_text(encoding="utf-8")

    prepare_index = build_sh.index("$PYTHON tools/prepare_browser_wheel.py")
    install_index = build_sh.index("$PYTHON -m pip install . --no-deps --ignore-installed -vv")
    assert prepare_index < install_index
    assert "python-build" in meta_yaml
    assert re.search(r"^\s+- wheel\s*$", meta_yaml, re.MULTILINE)


def test_hosted_web_build_refreshes_gallery_sessions_before_copy() -> None:
    deploy_yml = (REPO_ROOT / ".github" / "workflows" / "deploy_web.yml").read_text(encoding="utf-8")
    cloudflare_source = (REPO_ROOT / "tools" / "prepare_cloudflare_pages.py").read_text(encoding="utf-8")

    assert 'python -m pip install -e ".[dev]"' in deploy_yml
    refresh_index = deploy_yml.index("python tools/refresh_gallery_sessions.py")
    copy_index = deploy_yml.index("cp -r gbdraw/web/* public/")
    stamp_index = deploy_yml.index("python tools/stamp_web_build.py public")
    assert refresh_index < copy_index
    assert copy_index < stamp_index
    assert "refresh_gallery_sessions: bool = False" in cloudflare_source
    assert '"--refresh-gallery"' in cloudflare_source
    assert "refresh_gallery_sessions=args.refresh_gallery" in cloudflare_source
    assert "refresh_gallery_sessions_module.refresh_gallery_sessions()" in cloudflare_source
    assert "refresh_gallery_sessions_module.prepare_gallery_assets()" in cloudflare_source


def test_conda_recipe_does_not_copy_entire_web_tree() -> None:
    build_sh = (REPO_ROOT / "recipe" / "build.sh").read_text(encoding="utf-8")

    assert "cp -r gbdraw/web/*" not in build_sh
    assert "gbdraw/web/gallery" not in build_sh
    assert "gbdraw/web/index.html" in build_sh
    assert "gbdraw/web/open-source-notices.html" in build_sh
    assert "gbdraw/web/gbdraw-*.whl" in build_sh
    assert "for web_asset_dir in assets js presets vendor wasm" in build_sh


def test_setup_commands_refresh_open_source_notices() -> None:
    setup_source = (REPO_ROOT / "setup.py").read_text(encoding="utf-8")

    assert "class build_py(_build_py):" in setup_source
    assert "class sdist(_sdist):" in setup_source
    assert setup_source.count("_build_support.refresh_open_source_notices()") == 2
    assert 'cmdclass={"build_py": build_py, "sdist": sdist}' in setup_source


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
        build_root / "gbdraw" / "web" / "assets" / "favicon.ico",
        build_root / "gbdraw" / "web" / "assets" / "gbdraw-logo.svg",
        build_root / "gbdraw" / "web" / "assets" / "gbdraw-logo-title.svg",
        build_root / "gbdraw" / "web" / "assets" / "gbdraw-logo-title.png",
        build_root / "gbdraw" / "web" / verify_module._parse_wheel_name(),
        build_root / "gbdraw" / "web" / "vendor" / "vue" / "vue.global.js",
        build_root / "gbdraw" / "web" / "vendor" / "tailwindcss" / "tailwindcss-play.js",
        build_root / "gbdraw" / "web" / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
        build_root / "gbdraw" / "web" / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.wasm",
        build_root / "gbdraw" / "web" / "vendor" / "browser_wasi_shim" / "dist" / "index.js",
        build_root / "gbdraw" / "web" / "vendor" / "phosphor-icons" / "regular" / "style.css",
        build_root / "gbdraw" / "web" / "js" / "workers" / "losat-threaded-worker.js",
        build_root / "gbdraw" / "web" / "js" / "workers" / "losat-wasi-thread-worker.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "record-discovery.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "record-options.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "linear-record-selector.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "annotations" / "record-catalog.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "annotations" / "record-selector.js",
        build_root / "gbdraw" / "web" / "js" / "app" / "annotations" / "validation.js",
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat.wasm",
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat-threaded.wasm",
        *(build_root / "gbdraw" / "web" / path for path in verify_module.REQUIRED_UI_FONT_FILES),
        *(build_root / "gbdraw" / "web" / path for path in verify_module._parse_local_wheel_paths()),
    ]
    missing = [str(path.relative_to(build_root)) for path in required if not path.exists()]
    assert not missing, "build_py did not copy required offline GUI assets:\n" + "\n".join(missing)
    assert not (build_root / "gbdraw" / "web" / "gallery").exists()
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
    assert wheel_path.name == "gbdraw-0.14.0b0-py3-none-any.whl"
    subprocess.run(
        [sys.executable, "tools/verify_gui_offline.py", "inspect-wheel", str(wheel_path)],
        cwd=REPO_ROOT,
        check=True,
    )
    verify_module.assert_embedded_browser_wheel_is_not_recursive(wheel_path)

    with zipfile.ZipFile(wheel_path) as outer_wheel:
        outer_names = outer_wheel.namelist()
        browser_wheel_member = f"gbdraw/web/{verify_module._parse_wheel_name()}"
        browser_wheels = sorted(
            name
            for name in outer_names
            if name.startswith("gbdraw/web/gbdraw-") and name.endswith(".whl")
        )
        gallery_members = sorted(name for name in outer_names if name.startswith("gbdraw/web/gallery/"))
        assert browser_wheels == [browser_wheel_member]
        assert gallery_members == []
        assert "gbdraw/web/js/app/record-discovery.js" in outer_names
        assert "gbdraw/web/js/app/record-options.js" in outer_names
        assert "gbdraw/web/js/app/linear-record-selector.js" in outer_names
        assert "gbdraw/web/js/app/annotations/record-catalog.js" in outer_names
        assert "gbdraw/web/js/app/annotations/record-selector.js" in outer_names
        assert "gbdraw/web/js/app/annotations/validation.js" in outer_names


@pytest.mark.slow
def test_linear_record_selector_browser_flow(tmp_path: Path) -> None:
    playwright_sync_api = pytest.importorskip(
        "playwright.sync_api",
        reason="playwright is not available in this environment",
    )
    if not _can_bind_loopback():
        pytest.skip("loopback sockets are not permitted in this environment")

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    def write_genbank(path: Path, specs: list[tuple[str, int]]) -> None:
        records = []
        for index, (record_id, length) in enumerate(specs):
            record = SeqRecord(
                Seq(("ATGC" * ((length + 3) // 4))[:length]),
                id=record_id,
                name=record_id,
                description=f"Record {index + 1}",
            )
            record.annotations["molecule_type"] = "DNA"
            record.features = [
                SeqFeature(
                    FeatureLocation(0, length),
                    type="source",
                    qualifiers={"organism": [f"Organism {index + 1}"]},
                ),
                SeqFeature(FeatureLocation(0, min(80, length)), type="CDS"),
            ]
            records.append(record)
        SeqIO.write(records, path, "genbank")

    unique_gbk = tmp_path / "unique.gbk"
    duplicate_gbk = tmp_path / "duplicate.gbk"
    gff_path = tmp_path / "records.gff3"
    fasta_path = tmp_path / "records.fasta"
    session_path = tmp_path / "record-selector-session.json"
    write_genbank(unique_gbk, [("RecA", 1234), ("RecB", 567)])
    write_genbank(duplicate_gbk, [("RecA", 1234), ("RecA", 1100)])
    gff_path.write_text("##gff-version 3\nFastaA\t.\tgene\t1\t4\t.\t+\t.\tID=gene1\n", encoding="utf-8")
    fasta_path.write_text(">FastaA\nATGC\n>FastaB\nATGCGT\n", encoding="utf-8")

    ensure_prepared_browser_wheel()
    with _serve_repo_root() as base_url, playwright_sync_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.goto(f"{base_url}/gbdraw/web/index.html", wait_until="domcontentloaded")
        page.wait_for_function("() => window.__GBDRAW_APP__")
        page.evaluate(
            """() => {
                window.__gbdrawDialogMessages = [];
                window.alert = (message) => window.__gbdrawDialogMessages.push(String(message || ''));
            }"""
        )
        page.locator("button.app-mode-button").filter(has_text="Linear").click()
        page.wait_for_function("() => window.__GBDRAW_APP__?.pyodideReady", timeout=180_000)

        page.locator('input[type="file"][accept^=".gb,"]').first.set_input_files(str(unique_gbk))
        selector = page.locator("select[data-record-selector-uid]").first
        page.wait_for_function(
            """() => {
                const select = document.querySelector('select[data-record-selector-uid]');
                return select && !select.disabled && select.options.length === 3;
            }""",
            timeout=60_000,
        )
        assert selector.locator("option").all_text_contents() == [
            "Automatic (no explicit selector)",
            "RecA (1,234 bp)",
            "RecB (567 bp)",
        ]
        selector.select_option("RecB")
        assert page.evaluate("() => window.__GBDRAW_APP__.linearSeqs[0].region_record_id") == "RecB"

        result = page.evaluate("async () => await window.__GBDRAW_APP__.runAnalysis()")
        assert result["status"] == "ok"
        svg_content = page.evaluate("() => String(window.__GBDRAW_APP__.results[0]?.content || '')")
        assert "567 bp" in svg_content
        assert "1,234 bp" not in svg_content

        page.evaluate(
            """() => {
                const seq = window.__GBDRAW_APP__.linearSeqs[0];
                seq.region_start = 10;
                seq.region_end = 20;
            }"""
        )
        region_result = page.evaluate("async () => await window.__GBDRAW_APP__.runAnalysis()")
        assert region_result["status"] == "ok"
        region_svg = page.evaluate("() => String(window.__GBDRAW_APP__.results[0]?.content || '')")
        assert "10-20" in region_svg

        page.locator('input[type="file"][accept^=".gb,"]').first.set_input_files(str(duplicate_gbk))
        page.wait_for_function(
            """() => {
                const select = document.querySelector('select[data-record-selector-uid]');
                return select && !select.disabled &&
                    Array.from(select.options).some((option) => option.value === '#2');
            }""",
            timeout=60_000,
        )
        duplicate_labels = selector.locator("option").all_text_contents()
        assert "RecA (1,234 bp) [#1]" in duplicate_labels
        assert "RecA (1,100 bp) [#2]" in duplicate_labels
        selector.select_option("#2")

        page.evaluate("() => { window.__GBDRAW_APP__.sessionTitle = 'Record selector test'; }")
        with page.expect_download() as download_info:
            page.get_by_role("button", name="Save Session", exact=True).click()
        download_info.value.save_as(session_path)
        assert download_info.value.suggested_filename.endswith(".gbdraw-session.json.gz")
        assert session_path.read_bytes().startswith(b"\x1f\x8b")

        page.locator('input[type="file"][accept^=".gb,"]').first.set_input_files(str(unique_gbk))
        page.locator('input[type="file"][accept^=".json,"]').set_input_files(str(session_path))
        page.wait_for_function(
            """() => {
                const app = window.__GBDRAW_APP__;
                const select = document.querySelector('select[data-record-selector-uid]');
                return app?.linearSeqs?.[0]?.region_record_id === '#2' &&
                    select && !select.disabled && select.value === '#2';
            }""",
            timeout=60_000,
        )

        page.locator('input[type="radio"][value="gff"]').check()
        page.locator('input[type="file"][accept^=".gff,"]').first.set_input_files(str(gff_path))
        page.locator('input[type="file"][accept^=".fa,"]').first.set_input_files(str(fasta_path))
        page.wait_for_function(
            """() => {
                const select = document.querySelector('select[data-record-selector-uid]');
                return select && !select.disabled &&
                    Array.from(select.options).some((option) => option.value === 'FastaB');
            }""",
            timeout=60_000,
        )
        assert selector.locator("option").all_text_contents() == [
            "Automatic (no explicit selector)",
            "#2 (not found in current file)",
            "FastaA (4 bp)",
            "FastaB (6 bp)",
        ]
        browser.close()


@pytest.mark.slow
def test_linear_gff_feature_click_and_selection_smoke(tmp_path: Path) -> None:
    playwright_sync_api = pytest.importorskip(
        "playwright.sync_api",
        reason="playwright is not available in this environment",
    )
    if not _can_bind_loopback():
        pytest.skip("loopback sockets are not permitted in this environment")

    gff_path = tmp_path / "selectable.gff3"
    fasta_path = tmp_path / "selectable.fna"
    gff_path.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "##sequence-region GffRecord 1 90",
                "GffRecord\ttest\tgene\t1\t30\t.\t+\t.\tID=gene1;Name=example",
                "GffRecord\ttest\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=gene1;gene=example;product=selectable%20protein",
                "",
            ]
        ),
        encoding="utf-8",
    )
    fasta_path.write_text(">GffRecord\n" + ("ATG" * 30) + "\n", encoding="utf-8")

    ensure_prepared_browser_wheel()
    with _serve_repo_root() as base_url, playwright_sync_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.goto(f"{base_url}/gbdraw/web/index.html", wait_until="domcontentloaded")
        page.wait_for_function("() => window.__GBDRAW_APP__?.pyodideReady", timeout=180_000)
        page.locator("button.app-mode-button").filter(has_text="Linear").click()
        page.locator('input[type="radio"][value="gff"]').check()
        page.locator('input[type="file"][accept^=".gff,"]').first.set_input_files(str(gff_path))
        page.locator('input[type="file"][accept^=".fa,"]').first.set_input_files(str(fasta_path))

        result = page.evaluate("async () => await window.__GBDRAW_APP__.runAnalysis()")
        assert result["status"] == "ok"
        page.wait_for_function(
            "() => Array.isArray(window.__GBDRAW_APP__?.extractedFeatures) && window.__GBDRAW_APP__.extractedFeatures.length > 0",
            timeout=60_000,
        )
        target = page.evaluate(
            """async () => {
                const app = window.__GBDRAW_APP__;
                const featureIds = new Set(app.extractedFeatures.map((feature) => String(feature.svg_id || '')));
                const element = Array.from(document.querySelectorAll('[data-gbdraw-feature-id]'))
                    .find((candidate) => featureIds.has(String(candidate.getAttribute('data-gbdraw-feature-id') || '')));
                if (!element) return null;
                element.scrollIntoView({ block: 'center', inline: 'center' });
                await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
                const rect = element.getBoundingClientRect();
                return {
                    id: String(element.getAttribute('data-gbdraw-feature-id') || ''),
                    x: rect.left + rect.width / 2,
                    y: rect.top + rect.height / 2
                };
            }"""
        )
        assert target is not None

        page.evaluate(
            """({ id, x, y, ctrlKey }) => {
                const element = Array.from(document.querySelectorAll('[data-gbdraw-feature-id]'))
                    .find((candidate) => candidate.getAttribute('data-gbdraw-feature-id') === id);
                element.dispatchEvent(new MouseEvent('click', {
                    bubbles: true,
                    cancelable: true,
                    view: window,
                    clientX: x,
                    clientY: y,
                    ctrlKey
                }));
            }""",
            {**target, "ctrlKey": False},
        )
        page.wait_for_function(
            "(id) => window.__GBDRAW_APP__?.clickedFeature?.svg_id === id",
            arg=target["id"],
            timeout=20_000,
        )

        page.evaluate("() => { window.__GBDRAW_APP__.clickedFeature = null; }")
        page.evaluate(
            """({ id, x, y }) => {
                const element = Array.from(document.querySelectorAll('[data-gbdraw-feature-id]'))
                    .find((candidate) => candidate.getAttribute('data-gbdraw-feature-id') === id);
                element.dispatchEvent(new MouseEvent('click', {
                    bubbles: true,
                    cancelable: true,
                    view: window,
                    clientX: x,
                    clientY: y,
                    ctrlKey: true
                }));
            }""",
            target,
        )
        page.wait_for_function(
            "(id) => window.__GBDRAW_APP__?.selectedFeatureCount === 1 && window.__GBDRAW_APP__.selectedFeatureIds?.has(id)",
            arg=target["id"],
            timeout=20_000,
        )
        browser.close()


@pytest.mark.slow
def test_gallery_session_restore_smoke() -> None:
    playwright_sync_api = pytest.importorskip(
        "playwright.sync_api",
        reason="playwright is not available in this environment",
    )
    if not _can_bind_loopback():
        pytest.skip("loopback sockets are not permitted in this environment")

    ensure_prepared_browser_wheel()
    bad_console_terms = (
        "Definition update error",
        "float() argument",
        "JsNull",
        "JsUndefined",
    )

    with _serve_repo_root() as base_url, playwright_sync_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        console_errors: list[str] = []

        page.on(
            "console",
            lambda message: console_errors.append(message.text) if message.type == "error" else None,
        )
        page.on("pageerror", lambda error: console_errors.append(str(error)))

        for session_name in GALLERY_SESSION_FILES:
            page.goto(f"{base_url}/gbdraw/web/index.html", wait_until="domcontentloaded")
            page.wait_for_function("() => window.__GBDRAW_APP__")
            page.evaluate(
                """() => {
                    window.__gbdrawDialogMessages = [];
                    window.alert = (message) => {
                        window.__gbdrawDialogMessages.push(String(message || ''));
                    };
                }"""
            )
            page.locator('input[accept^=".json,"]').first.set_input_files(
                str(GALLERY_ROOT / "sessions" / session_name)
            )
            page.wait_for_function(
                "() => Array.isArray(window.__gbdrawDialogMessages) && window.__gbdrawDialogMessages.length > 0",
                timeout=180_000,
            )
            dialog_message = page.evaluate("() => window.__gbdrawDialogMessages.at(-1)")
            page.wait_for_function(
                """() => {
                    const app = window.__GBDRAW_APP__;
                    return Array.isArray(app?.results) &&
                        app.results.length > 0 &&
                        !app.featureExtractionPending;
                }""",
                timeout=180_000,
            )

            assert dialog_message == "Session loaded successfully!"
            summary = page.evaluate(
                """() => {
                    const app = window.__GBDRAW_APP__;
                    const index = Number.isInteger(app.selectedResultIndex) ? app.selectedResultIndex : 0;
                    const selected = app.results?.[index] || app.results?.[0] || {};
                    return {
                        resultCount: Array.isArray(app.results) ? app.results.length : 0,
                        hasSvg: String(selected.content || '').includes('<svg'),
                        status: app.featureEditorStatus?.status || '',
                        featureExtractionError: app.featureExtractionError || null,
                        extractedCount: Array.isArray(app.extractedFeatures) ? app.extractedFeatures.length : 0
                    };
                }"""
            )
            assert summary["resultCount"] > 0, session_name
            assert summary["hasSvg"], session_name
            assert summary["status"] == "summary-ready", session_name
            assert summary["featureExtractionError"] in (None, ""), session_name
            assert summary["extractedCount"] > 0, session_name

            if session_name == "Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz":
                multi_record_positions = page.evaluate(
                    """() => {
                        const app = window.__GBDRAW_APP__;
                        return (Array.isArray(app?.adv?.multi_record_positions) ? app.adv.multi_record_positions : [])
                            .map((entry) => `${entry.selector}@${entry.row}`);
                    }"""
                )
                assert multi_record_positions == [
                    "#1@1",
                    "#2@1",
                    "#3@2",
                    "#4@2",
                    "#5@2",
                    "#6@2",
                ]

            if session_name in GALLERY_MULTI_RECORD_LINEAR_SESSION_FILES:
                target = page.evaluate(
                    """async () => {
                        const app = window.__GBDRAW_APP__;
                        const svg = document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement ||
                            document.querySelector('svg');
                        const featuresBySvgId = new Map(
                            (Array.isArray(app.extractedFeatures) ? app.extractedFeatures : [])
                                .map((feature) => [String(feature?.svg_id || '').trim(), feature])
                                .filter(([id]) => id)
                        );
                        const candidates = Array.from(svg?.querySelectorAll('[data-gbdraw-feature-id]') || [])
                            .filter((candidate) => {
                                const id = String(candidate.getAttribute('data-gbdraw-feature-id') || candidate.id || '')
                                    .replace(/__part\\d+$/, '');
                                return featuresBySvgId.has(id);
                            });
                        const element = candidates.find((candidate) => {
                            const rect = candidate.getBoundingClientRect();
                            return rect.width > 0 && rect.height > 0;
                        }) || candidates[0];
                        if (!element) return null;
                        element.scrollIntoView({ block: 'center', inline: 'center' });
                        await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
                        const rect = element.getBoundingClientRect();
                        return {
                            id: String(element.getAttribute('data-gbdraw-feature-id') || element.id || '')
                                .replace(/__part\\d+$/, ''),
                            x: rect.left + rect.width / 2,
                            y: rect.top + rect.height / 2
                        };
                    }"""
                )
                assert target is not None, session_name
                page.evaluate(
                    """(featureId) => {
                        const svg = document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement ||
                            document.querySelector('svg');
                        const element = Array.from(svg?.querySelectorAll('[data-gbdraw-feature-id]') || [])
                            .find((candidate) => {
                                const id = String(candidate.getAttribute('data-gbdraw-feature-id') || candidate.id || '')
                                    .replace(/__part\\d+$/, '');
                                return id === featureId;
                            });
                        if (!element) return;
                        const rect = element.getBoundingClientRect();
                        element.dispatchEvent(new MouseEvent('click', {
                            bubbles: true,
                            cancelable: true,
                            view: window,
                            clientX: rect.left + rect.width / 2,
                            clientY: rect.top + rect.height / 2
                        }));
                    }""",
                    target["id"],
                )
                try:
                    page.wait_for_function(
                        """(featureId) => {
                            const clicked = window.__GBDRAW_APP__?.clickedFeature;
                            return clicked?.svg_id === featureId || Boolean(document.querySelector('.feature-popup'));
                        }""",
                        arg=target["id"],
                        timeout=20_000,
                    )
                except Exception as error:
                    debug = page.evaluate(
                        """(featureId) => {
                            const app = window.__GBDRAW_APP__;
                            return {
                                featureId,
                                clickedFeature: app?.clickedFeature?.svg_id || '',
                                popupVisible: Boolean(document.querySelector('.feature-popup')),
                                status: app?.featureEditorStatus?.status || '',
                                extractedCount: Array.isArray(app?.extractedFeatures) ? app.extractedFeatures.length : 0,
                                hasFeature: (Array.isArray(app?.extractedFeatures) ? app.extractedFeatures : [])
                                    .some((feature) => String(feature?.svg_id || '').trim() === featureId)
                            };
                        }""",
                        target["id"],
                    )
                    raise AssertionError(f"{session_name} feature click did not resolve: {debug}") from error
                clicked = page.evaluate(
                    "() => window.__GBDRAW_APP__?.clickedFeature ? { svg_id: window.__GBDRAW_APP__.clickedFeature.svg_id } : null"
                )
                assert clicked is not None, session_name
                assert clicked["svg_id"] == target["id"], session_name

            matching_console_errors = [
                message
                for message in console_errors
                if any(term in message for term in bad_console_terms)
            ]
            assert matching_console_errors == [], session_name

        browser.close()


@pytest.mark.slow
def test_circular_radial_label_browser_flow_and_reflow() -> None:
    playwright_sync_api = pytest.importorskip(
        "playwright.sync_api",
        reason="playwright is not available in this environment",
    )
    if not _can_bind_loopback():
        pytest.skip("loopback sockets are not permitted in this environment")

    ensure_prepared_browser_wheel()
    input_path = REPO_ROOT / "tests" / "test_inputs" / "MjeNMV.gbk"
    with _serve_repo_root() as base_url, playwright_sync_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.goto(f"{base_url}/gbdraw/web/index.html", wait_until="domcontentloaded")
        page.wait_for_function("() => window.__GBDRAW_APP__?.pyodideReady", timeout=180_000)
        assert page.evaluate("() => window.__GBDRAW_APP__?.adv?.circular_label_placement") == "horizontal"
        page.evaluate(
            """async () => {
                window.__GBDRAW_APP__.form.labels_mode = 'out';
                await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
            }"""
        )

        controls = page.evaluate(
            """() => {
                const app = window.__GBDRAW_APP__;
                const label = Array.from(document.querySelectorAll('label'))
                    .find((candidate) => candidate.textContent.includes('Circular Label Placement'));
                const select = label?.parentElement?.querySelector('select');
                return {
                    state: app?.adv?.circular_label_placement,
                    value: select?.value,
                    options: Array.from(select?.options || []).map((option) => option.value),
                    geometryVisible: Array.from(document.querySelectorAll('h4'))
                        .some((heading) => heading.textContent.includes('LABEL GEOMETRY'))
                };
            }"""
        )
        assert controls == {
            "state": "horizontal",
            "value": "horizontal",
            "options": ["horizontal", "radial"],
            "geometryVisible": True,
        }

        page.locator('input[type="file"][accept^=".gb,"]').first.set_input_files(str(input_path))
        page.evaluate(
            """() => {
                const app = window.__GBDRAW_APP__;
                app.form.labels_mode = 'both';
                app.form.legend = 'none';
                app.adv.label_rendering = 'external_only';
                app.adv.circular_label_placement = 'radial';
            }"""
        )
        page.wait_for_function(
            """() => !Array.from(document.querySelectorAll('h4'))
                .some((heading) => heading.textContent.includes('LABEL GEOMETRY'))"""
        )

        result = page.evaluate("async () => await window.__GBDRAW_APP__.runAnalysis()")
        assert result["status"] == "ok"
        page.wait_for_function(
            """() => window.__GBDRAW_APP__?.featureEditorStatus?.status === 'summary-ready' &&
                Boolean(document.querySelector('#label_text text'))""",
            timeout=180_000,
        )
        page.evaluate("() => window.__GBDRAW_APP__.syncLabelEditor()")
        page.wait_for_function(
            "() => Array.isArray(window.__GBDRAW_APP__?.editableLabels) && window.__GBDRAW_APP__.editableLabels.length > 0",
            timeout=180_000,
        )

        def radial_stats() -> dict:
            return page.evaluate(
                """() => {
                    const content = String(window.__GBDRAW_APP__.results?.[0]?.content || '');
                    const document = new DOMParser().parseFromString(content, 'image/svg+xml');
                    return {
                        textCount: document.querySelectorAll('#label_text text[transform^="rotate("]').length,
                        lineCount: document.querySelectorAll('#label_leaders line').length,
                        content
                    };
                }"""
            )

        before = radial_stats()
        assert before["textCount"] > 0
        assert before["lineCount"] == 2 * before["textCount"]

        edited_text = "radial browser reflow label with a deliberately long value"
        edit_started = page.evaluate(
            """async (nextText) => {
                const app = window.__GBDRAW_APP__;
                const entry = app.editableLabels.find((candidate) => candidate.featureId);
                if (!entry) return false;
                await app.requestLabelTextChangeByFeatureId(entry.featureId, nextText);
                await app.handleLabelTextScopeChoice('single');
                return true;
            }""",
            edited_text,
        )
        assert edit_started
        page.wait_for_function(
            """(nextText) => {
                const app = window.__GBDRAW_APP__;
                const content = String(app?.results?.[0]?.content || '');
                return !app?.processing && !app?.labelReflowProcessing && content.includes(nextText);
            }""",
            arg=edited_text,
            timeout=180_000,
        )
        after = radial_stats()
        assert edited_text in after["content"]
        assert after["textCount"] == before["textCount"]
        assert after["lineCount"] == 2 * after["textCount"]
        browser.close()


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
