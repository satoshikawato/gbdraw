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
    expected_wheel_name = "gbdraw-0.12.0b0-py3-none-any.whl"
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
    assert "args.push('--ruler_label_font_size', adv.scale_font_size);" in source
    assert "args.push('--scale_font_size', adv.scale_font_size);" in source


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
    assert "args.push('--circular_track_axis_index', String(adv.circular_track_slots_axis_index));" in run_source
    assert "buildCircularTrackSlotSpec(slot, adv.nt, form.track_type, {" in run_source
    assert "applyCircularTrackOrderPlacements(" in run_source
    assert "if (!useCircularTrackSlots)" in run_source
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
    assert "Replace the current custom circular track slots with this preset" not in slot_source
    assert "setCircularTrackSlotsEnabled" in slot_source
    assert "const templateSlots = createDefaultCircularTrackSlots" in slot_source
    assert "state.adv.circular_track_slots.splice(0, state.adv.circular_track_slots.length, ...normalized);" in slot_source
    assert "() => [adv.circular_track_slots_enabled, form.show_depth]" in app_setup_source
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
    assert "c_conservation_blasts: await serializeFileArray" in config_source
    assert "normalizeCircularConservationReference" in config_source
    assert "Pairwise Comparisons" in index_html
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


def test_circular_track_slot_axis_crossing_actions_keep_neighbor_sides(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "circular-track-slots.js"
    module_path = tmp_path / "circular-track-slots.mjs"
    (tmp_path / "package.json").write_text('{"type":"module"}', encoding="utf-8")
    for dependency in ["conservation-series.js", "color-utils.js"]:
        dep_path = WEB_ROOT / "js" / "app" / dependency
        (tmp_path / dependency).write_text(dep_path.read_text(encoding="utf-8"), encoding="utf-8")
    module_path.write_text(source_path.read_text(encoding="utf-8"), encoding="utf-8")
    check_path = tmp_path / "check-circular-track-slots.mjs"
    check_path.write_text(
        f"""
        import {{
          applyCircularTrackOrderPlacements,
          createDefaultCircularTrackSlots,
          createCircularTrackSlotEditor,
          estimateCircularConservationLayoutWarning
        }} from {module_path.as_uri()!r};

        const defaultSlots = createDefaultCircularTrackSlots({{ preset: 'tuckin' }});
        const defaultTick = defaultSlots.find((slot) => slot.id === 'ticks');
        if (defaultTick?.params?.tick_label_layout !== 'label_in_tick_out') {{
          throw new Error(`Default Tick layout should point labels inward when Tick is inside Feature: ${{JSON.stringify(defaultTick)}}`);
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

    assert "circular_track_slots_schema_version: 3" in state_source
    assert "const CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 3;" in config_source
    assert "const LEGACY_CIRCULAR_TRACK_SLOT_SCHEMA_VERSION = 2;" in config_source
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
    assert wheel_path.name == "gbdraw-0.12.0b0-py3-none-any.whl"
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
