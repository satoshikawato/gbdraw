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
    expected_wheel_name = "gbdraw-0.11.0-py3-none-any.whl"
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
    assert '<img src="./assets/gbdraw-logo.svg" alt="" class="animate-spin w-16 h-16 mb-6">' in index_html


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


def test_web_run_analysis_wires_circular_track_slot_options() -> None:
    run_source = (WEB_ROOT / "js" / "app" / "run-analysis.js").read_text(encoding="utf-8")
    state_source = (WEB_ROOT / "js" / "state.js").read_text(encoding="utf-8")
    config_source = (WEB_ROOT / "js" / "services" / "config.js").read_text(encoding="utf-8")
    slot_source = (WEB_ROOT / "js" / "app" / "circular-track-slots.js").read_text(encoding="utf-8")
    app_setup_source = (WEB_ROOT / "js" / "app" / "app-setup.js").read_text(encoding="utf-8")
    index_html = (WEB_ROOT / "index.html").read_text(encoding="utf-8")

    assert "circular_track_slots_enabled" in state_source
    assert "circular_track_slots_axis_index" in state_source
    assert "createDefaultCircularTrackSlots()" in state_source
    assert "inferLegacyAxisIndexFromFeature(normalizedSlots, state.form.track_type)" in config_source
    assert '"circular_track_slot": "--circular_track_slot" in _source' in run_source
    assert '"circular_track_axis_index": "--circular_track_axis_index" in _source' in run_source
    assert "args.push('--track_type', form.track_type);" in run_source
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


def test_circular_track_slot_axis_crossing_actions_keep_neighbor_sides(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    source_path = WEB_ROOT / "js" / "app" / "circular-track-slots.js"
    module_path = tmp_path / "circular-track-slots.mjs"
    module_path.write_text(source_path.read_text(encoding="utf-8"), encoding="utf-8")
    check_path = tmp_path / "check-circular-track-slots.mjs"
    check_path.write_text(
        f"""
        import {{ createCircularTrackSlotEditor }} from {module_path.as_uri()!r};

        const state = {{
          adv: {{
            nt: 'GC',
            circular_track_slots: [
              {{ id: 'gc_content', renderer: 'dinucleotide_content', side: 'outside', params: {{ nt: 'GC' }} }},
              {{ id: 'features', renderer: 'features', side: 'inside', params: {{ lane_direction: 'inside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ label_side: 'inside', tick_side: 'inside' }} }},
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
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ label_side: 'inside', tick_side: 'inside' }} }},
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
        if (tickSlot.params.label_side !== 'inside' || tickSlot.params.tick_side !== 'inside') {{
          throw new Error(`Tick params did not resync inside: ${{JSON.stringify(tickSlot.params)}}`);
        }}

        const outsideFeatureState = {{
          adv: {{
            nt: 'GC',
            circular_track_slots: [
              {{ id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'outside', params: {{ nt: 'GC' }} }},
              {{ id: 'features', renderer: 'features', side: 'outside', params: {{ lane_direction: 'outside' }} }},
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ label_side: 'inside', tick_side: 'inside' }} }},
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
              {{ id: 'ticks', renderer: 'ticks', side: 'inside', params: {{ label_side: 'inside', tick_side: 'inside' }} }},
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
    assert wheel_path.name == "gbdraw-0.11.0-py3-none-any.whl"
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
