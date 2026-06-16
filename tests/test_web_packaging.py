from __future__ import annotations

import importlib.util
import shutil
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
    expected_wheel_name = "gbdraw-0.12.0-py3-none-any.whl"
    assert verify_module._parse_wheel_name() == expected_wheel_name
    assert expected_wheel_path.name == expected_wheel_name
    verify_module.assert_browser_wheel_is_not_recursive(expected_wheel_path)
    verify_module._assert_packaged_assets()


def test_index_links_to_open_source_notices() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "./open-source-notices.html" in index_html
    assert "Open Source Notices" in index_html


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
    assert "frame-ancestors" not in index_html
    assert "frame-ancestors" not in notices_html


def test_index_cloaks_vue_template_until_mount() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert '<div id="app" v-cloak' in index_html
    assert "[v-cloak] { display: none !important; }" in index_html
    assert 'id="app-boot-splash"' in index_html
    assert "Initializing gbdraw..." in index_html
    assert "#app:not([v-cloak]) + #app-boot-splash" in index_html


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
    assert "qualifierRows" in svg_actions_source
    assert "locationParts" in svg_actions_source
    assert "nucleotideSequence" in svg_actions_source
    assert "aminoAcidSequence" in svg_actions_source
    assert "const displayProteinId = (feat)" in svg_actions_source
    assert "label: 'Protein ID', value: proteinId" in svg_actions_source
    assert "label: 'Source protein ID'" not in svg_actions_source
    assert "label: 'SVG ID'" not in svg_actions_source
    assert "label: 'Record index'" not in svg_actions_source
    assert "label: 'Strand'" not in svg_actions_source
    assert "navigator.clipboard?.writeText" in app_setup_source
    assert "location_parts" in helper_source
    assert "nucleotide_sequence" in helper_source
    assert "amino_acid_sequence" in helper_source
    assert "sanitizeExtractedFeaturesForSession(state.extractedFeatures.value)" in config_source


def test_interactive_svg_export_decouples_interactivity_from_rich_popup_payload() -> None:
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "gbdraw-interactive-feature-popup-v1" in export_source
    assert "gbdraw-interactive-feature-metadata" in export_source
    assert "gbdraw-interactive-feature-script" in export_source
    assert "gbdraw-feature-search-controls" in export_source
    assert "gbdraw-interactive-feature--match" in export_source
    assert "gbdraw-interactive-feature--active-match" in export_source
    assert "gbdraw-interactive-feature--dimmed" in export_source
    assert "gbdraw-interactive-feature-match-glow" in export_source
    assert "filter: url(#gbdraw-interactive-feature-match-glow);" in export_source
    assert "stroke-opacity: 0.6;" in export_source
    assert "stroke-opacity: 1;" in export_source
    assert "function normalizeSearchText(value)" in export_source
    assert "function featureSearchValues(feature, field, qualifierKey)" in export_source
    assert "function featureSearchMatches(feature, matcher, field, qualifierKey)" in export_source
    assert "function compileSearchMatcher(query, useRegex)" in export_source
    assert "function featureMatchesSearch(feature, matcher, field, qualifierKey)" in export_source
    assert "matchDetails: {}" in export_source
    assert "data-search-match-detail" in export_source
    assert "gfs-button--clear" in export_source
    assert "gfs-match-detail" in export_source
    assert "['orthogroup', 'Orthogroup']" in export_source
    assert "['nucleotide', 'Nucleotide']" in export_source
    assert "['amino-acid', 'Amino acid']" in export_source
    assert "['sequence', 'Sequence']" not in export_source
    assert "var NUCLEOTIDE_IUPAC = {" in export_source
    assert "var AMINO_ACID_IUPAC = {" in export_source
    assert "function buildIupacQueryPattern(query, alphabet)" in export_source
    assert "function supportsStandaloneControls()" in export_source
    assert "function setSearchState(nextState)" in export_source
    assert "function applySearchResults()" in export_source
    assert "function setActiveMatch(index, options)" in export_source
    assert "function clearSearch()" in export_source
    assert "Search match" in export_source
    assert "Orthogroup members" in export_source
    assert "function scheduleInitialViewportRefresh()" in export_source
    assert "var initialView = copyViewRect(getViewRect());" in export_source
    assert "rectsNearlyEqual(getViewRect(), initialView)" in export_source
    assert "scheduleInitialViewportRefresh();" in export_source
    assert "visibleView.x + visibleView.width - (controlWidth * unit) - margin" in export_source
    assert "visibleView.y + margin + (searchControlsOffsetCss.y * unit)" in export_source
    assert "var yOffset = 42 * unit" not in export_source
    assert "gfi-og-members-table" in export_source
    assert "Coordinates (+/-)" in export_source
    assert "Product / note" in export_source
    assert "displayProteinId(null, member)" in export_source
    assert "function displayProteinId(feature, member)" in export_source
    assert "['Source protein ID'" not in export_source
    assert "display_label" in export_source
    assert "search_labels" in export_source
    assert "orthogroup_id" in export_source
    assert "protein_id" in export_source
    assert "const buildStandaloneOrthogroupPayloads = (features) => {" in export_source
    assert "const orthogroups = buildStandaloneOrthogroupPayloads(features);" in export_source
    assert "data-gbdraw-interactive-feature" in export_source
    assert "data-gbdraw-original-viewbox" in export_source
    assert "data-gbdraw-original-width" in export_source
    assert "data-gbdraw-original-height" in export_source
    assert "enrichSvgWithStandaloneFeaturePopup" not in export_source
    assert "const enrichSvgWithStandaloneInteractivity = (svg, { popupMode = 'rich' } = {}) => {\n  if (!svg) return false;" in export_source
    assert "if (!svg || state.adv.rich_feature_popup === false) return false;" not in export_source
    assert "popupMode: state.adv.rich_feature_popup === false ? 'simple' : 'rich'" in export_source
    assert "buildStandaloneFeaturePayloads(svg, { popupMode: normalizedPopupMode })" in export_source
    assert "svg.setAttribute('width', '100vw');" in export_source
    assert "svg.setAttribute('height', '100vh');" in export_source
    assert "svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');" in export_source
    assert "svg.style.setProperty('width', '100vw');" in export_source
    assert "svg.style.setProperty('height', '100vh');" in export_source
    assert "function parseOriginalViewRectFromSvg()" in export_source
    assert "function getViewportAspect()" in export_source
    assert "function fitRectToAspect(rect, targetAspect)" in export_source
    assert "var homeViewRect = fitRectToAspect(originalViewRect, getViewportAspect());" in export_source
    assert "homeViewRect.width / maxZoom" in export_source
    assert "homeViewRect.width * nextScale" in export_source
    assert "{ action: 'reset', label: 'Original', title: 'Return to original view', width: 62 }" in export_source
    assert "{ action: 'pan', label: 'Pan'" not in export_source
    assert "{ action: 'legend', label: 'Legend'" not in export_source
    assert "function ensureStickyLegendBackground(legend, bbox)" in export_source
    assert "gbdraw-sticky-legend-background" in export_source
    assert "setClassToken(svg, 'gbdraw-interactive-pan-enabled', true);" in export_source
    assert "function refitViewportToWindow()" in export_source
    assert "function scheduleViewportRefit()" in export_source
    assert "popup_mode: normalizedPopupMode" in export_source
    assert "orthogroups" in export_source
    assert "var popupMode = payload.popup_mode === 'simple' ? 'simple' : 'rich';" in export_source
    assert "var orthogroups = Array.isArray(payload.orthogroups) ? payload.orthogroups : [];" in export_source
    assert "var richSearchFields = {" in export_source
    assert "if (popupMode === 'simple') {\n      searchFieldOptions = searchFieldOptions.filter" in export_source
    assert "function renderSimplePopup(feature)" in export_source
    assert "if (normalizedPopupMode === 'rich') {\n      Object.assign(payload, {\n        qualifiers: normalizeQualifierMap(feature?.qualifiers)," in export_source
    assert "nucleotide_sequence" in export_source
    assert "amino_acid_sequence" in export_source
    assert "getVisibleViewRect()" in export_source
    assert "var visibleView = getVisibleViewRect();" in export_source
    assert "window.addEventListener('scroll', updateViewportControlsPosition, { passive: true });" in export_source
    assert "window.addEventListener('resize', scheduleViewportRefit);" in export_source
    assert "window.visualViewport.addEventListener('scroll', updateViewportControlsPosition, { passive: true });" in export_source
    assert "popupCssWidth" in export_source
    assert "getPopupCssMetrics" in export_source
    assert "var effectiveScaleX = safeScaleX * metrics.zoomScale;" in export_source
    assert "var marginCss = metrics.margin;" in export_source
    assert "var dragZoomScale = getBrowserZoomScale(getViewportClientRect());" in export_source
    assert "--gfi-text-scale" in export_source
    assert "getPopupTextScale" in export_source
    assert "root.style.setProperty('--gfi-text-scale'" in export_source
    assert "gbdraw-interactive-feature-glow" in export_source
    assert "gbdraw-interactive-feature-match-glow" in export_source
    assert "gbdraw-interactive-feature--hover" in export_source
    assert "gbdraw-interactive-orthogroup-link--hover" in export_source
    assert "function setOrthogroupHover(orthogroupId, highlight)" in export_source
    assert "activePopupDrag" in export_source
    assert "activeSearchControlsDrag" in export_source
    assert "gbdraw-feature-hover-popup" in export_source
    assert "function scheduleHoverPopup(feature, svgId, event)" in export_source
    assert "function renderHoverPopupHtml(feature, svgId)" in export_source
    assert "svg.addEventListener('mousemove'" in export_source
    assert "const collectRenderedFeatureEntries = (svg) => {" in export_source
    assert "const buildFallbackStandaloneFeaturePayload = (svgId, entry, captionsByColor) => {" in export_source
    assert "function startSearchControlsDrag(event, root)" in export_source
    assert "document.addEventListener('mouseup', onEnd, true);" in export_source
    assert "window.addEventListener('blur', onEnd);" in export_source
    assert 'data-drag-handle="true"' in export_source
    assert "function startPopupDrag(event)" in export_source
    assert "setFeatureHighlight" in export_source
    assert "svg.addEventListener('mouseover'" in export_source
    assert "['SVG ID'" not in export_source
    assert "['Record index'" not in export_source
    assert "['Strand'" not in export_source
    assert "root.style.transform = 'scale('" in export_source
    assert "overscroll-behavior: contain;" in export_source
    assert "root.addEventListener('wheel', function (rootEvent) {\n      rootEvent.stopPropagation();\n    }, { passive: true });" in export_source
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


def test_plain_svg_export_strips_editor_only_cursor_affordances() -> None:
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")

    assert "const stripEditorOnlyCursorStyles = (svg) => {" in export_source
    assert "svg.querySelectorAll('[style]').forEach((element) => {" in export_source
    assert "if (!style || !/\\bcursor\\s*:/i.test(style)) return;" in export_source
    assert "element.style.removeProperty('cursor');" in export_source
    assert "if (!element.getAttribute('style')?.trim()) {" in export_source
    assert "  } else {\n    stripEditorOnlyCursorStyles(clone);\n  }\n  return new XMLSerializer().serializeToString(clone);" in export_source
    assert "export const downloadInteractiveSVG = () => {\n  const svgString = getCurrentSvgString({ interactive: true });" in export_source


def test_local_index_keeps_cloudflare_analytics_as_deploy_only() -> None:
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")
    assert "static.cloudflareinsights.com" not in index_html
    assert "CLOUDFLARE_WEB_ANALYTICS_SCRIPT" in index_html
    assert "CLOUDFLARE_WEB_ANALYTICS_NOTICE" in index_html


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
    assert "losatRuntimeCompatibility" in run_source
    assert "onRuntimeStatus" in run_source


def test_web_losat_thread_count_options_are_contiguous() -> None:
    source = (WEB_ROOT / "js" / "app" / "losat-settings.js").read_text(encoding="utf-8")
    service_source = (WEB_ROOT / "js" / "services" / "losat.js").read_text(encoding="utf-8")
    worker_source = (WEB_ROOT / "js" / "workers" / "losat-threaded-worker.js").read_text(encoding="utf-8")
    assert "const createPositiveIntegerOptions = (maxValue) =>" in source
    assert "return createPositiveIntegerOptions(losatHardwareThreads.value);" in source
    assert "return createPositiveIntegerOptions(maxThreads);" in source
    assert "const perJobSlots = losatEffectiveThreadsPerJob.value;" in source
    assert "const workersPerThreadedJob = Math.max(1, threadsPerJob);" in service_source
    assert "getChildWorkerCount(effectiveThreads)" in worker_source
    assert "Array.from({ length: childWorkerCount }" in worker_source


def test_cloudflare_bundle_includes_analytics_and_hosted_notice(tmp_path: Path) -> None:
    ensure_prepared_browser_wheel()
    cloudflare_module = _load_prepare_cloudflare_pages_module()
    output_root = tmp_path / "cloudflare-pages"
    bundle_path = cloudflare_module.build_cloudflare_pages_bundle(output_root=output_root)

    index_html = (bundle_path / "index.html").read_text(encoding="utf-8")
    assert "https://static.cloudflareinsights.com/beacon.min.js" in index_html
    assert 'data-cf-beacon=\'{"token": "e4dc2e66d09549868f5a5ac7d7a6e633"}\'' in index_html
    assert "Hosted Site Analytics" in index_html
    assert "Uploaded genome files are still processed locally in your browser" in index_html
    assert "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://static.cloudflareinsights.com;" in index_html
    assert "connect-src 'self' https://cloudflareinsights.com;" in index_html
    assert "CLOUDFLARE_WEB_ANALYTICS_SCRIPT" not in index_html
    assert "CLOUDFLARE_WEB_ANALYTICS_NOTICE" not in index_html
    headers = (bundle_path / "_headers").read_text(encoding="utf-8")
    assert "Cross-Origin-Opener-Policy: same-origin" in headers
    assert "Cross-Origin-Embedder-Policy: require-corp" in headers
    assert "Cross-Origin-Resource-Policy: same-origin" in headers
    assert "Content-Security-Policy: frame-ancestors 'none'" in headers


def test_wrangler_uses_cloudflare_bundle_directory() -> None:
    wrangler_toml = (REPO_ROOT / "wrangler.toml").read_text(encoding="utf-8")
    assert 'directory = "./dist/cloudflare-pages"' in wrangler_toml
    assert 'not_found_handling = "single-page-application"' in wrangler_toml


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
    svg_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "svg-actions.js").read_text(encoding="utf-8")
    label_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "label-actions.js").read_text(encoding="utf-8")
    color_actions_source = (WEB_ROOT / "js" / "app" / "feature-editor" / "color-actions.js").read_text(encoding="utf-8")
    stroke_actions_source = (WEB_ROOT / "js" / "app" / "legend" / "stroke-actions.js").read_text(encoding="utf-8")
    svg_styles_source = (WEB_ROOT / "js" / "app" / "svg-styles.js").read_text(encoding="utf-8")
    orthogroups_source = (WEB_ROOT / "js" / "app" / "orthogroups.js").read_text(encoding="utf-8")
    export_source = (WEB_ROOT / "js" / "services" / "export.js").read_text(encoding="utf-8")

    assert "FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id'" in svg_actions_source
    assert "FEATURE_SELECTOR = [" in svg_actions_source
    assert "`path[${FEATURE_ID_ATTRIBUTE}]`" in svg_actions_source
    assert "element?.getAttribute?.(FEATURE_ID_ATTRIBUTE)" in svg_actions_source
    assert "svg.querySelectorAll(FEATURE_SELECTOR)" in svg_actions_source
    assert "getFeatureIdentity(element)" in svg_actions_source
    assert "getFeatureHoverKey(getFeatureIdentity(relatedFeature))" in svg_actions_source
    assert "getFeatureIdentity(el)" in label_actions_source
    assert "getFeatureElements(svg, feat.svg_id)" in color_actions_source
    assert "getFeatureElements(svg, svgId)" in color_actions_source
    assert "getFeatureElements(svg, svgId)" in stroke_actions_source
    assert "getFeatureIdentity(path)" in svg_styles_source
    assert "getFeatureElements(svg, featureId)" in orthogroups_source
    assert "FEATURE_ID_ATTRIBUTE = 'data-gbdraw-feature-id'" in export_source
    assert "function getElementFeatureId(element)" in export_source
    assert "var svgId = getElementFeatureId(featureElement);" in export_source
    assert "const id = getElementFeatureId(element);" in export_source


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
    assert '>Custom Track Slots <help-tip text="Edit the explicit linear track stack.' in index_html
    assert "Linear Custom Track Slots" not in index_html
    assert "Linear Custom Track Slots" not in config_source
    assert 'v-model="adv.linear_track_slots_enabled"' in index_html
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
    for dependency in ["depth-track-state.js"]:
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
          linearSeqs: []
        }};
        const editor = createLinearTrackSlotEditor({{ state }});
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
        reconcileState.form.show_gc = true;
        reconcileEditor.syncLinearNumericSlotsFromSimpleControls();
        const restoredGc = reconcileState.adv.linear_track_slots.find((slot) => slot.id === 'gc_content');
        assert(restoredGc?.params?.nt === 'AT', 'Show GC on should restore default GC using the simple nt control');
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
    for dependency in ["conservation-series.js", "color-utils.js", "depth-track-state.js"]:
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
          }}
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
    assert "validateImportedCircularTrackSlots(data.config);" in config_source
    assert "Failed to load config: ${message}" in config_source
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
    assert "downloadJson(sessionData, sessionFilename, { pretty: false });" in config_source


def test_web_config_persists_manual_qualifier_priority_rules() -> None:
    source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    assert "qualifierPriorityRules: cloneQualifierPriorityRules(state.manualPriorityRules)" in source
    assert "replaceQualifierPriorityRules(data.qualifierPriorityRules)" in source
    assert "replaceQualifierPriorityRules(data.priorityRules)" in source


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
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat.wasm",
        build_root / "gbdraw" / "web" / "wasm" / "losat" / "losat-threaded.wasm",
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
    assert wheel_path.name == "gbdraw-0.12.0-py3-none-any.whl"
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
