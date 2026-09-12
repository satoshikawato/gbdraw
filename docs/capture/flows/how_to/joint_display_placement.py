"""Capture rotation and multipart placement using the existing plastome recipe."""
from __future__ import annotations

import gzip
import json
from pathlib import Path

from playwright.sync_api import expect

from flows.tutorials.gui_annotated_chloroplast import (
    GUI_ANNOTATION_GENBANK_PATH, GUI_ANNOTATION_GENBANK_SHA256,
    GUI_ANNOTATION_GENBANK_SIZE, GUI_ANNOTATION_TABLE_PATH,
    _configure_gallery_presentation, _configure_gallery_slots,
    _inspect_tracks_svg, _assert_safe_svg, _fit_circular_preview,
    _validate_complete_record, _validate_annotation_fixture, CaptureResult,
)
from flows.web_capture import (
    assert_fixture_identity, assert_output_paths, generate_and_inspect,
    open_browser_capture, wait_for_app_shell,
)

NAMES = ('01-record-start.png', '02-feature-placement.png')


def capture_joint_display_placement(browser_type, base_url, output_paths, download_dir):
    """Use real controls; keep the annotated plastome and saved replay evidence."""
    report = _validate_complete_record(
        GUI_ANNOTATION_GENBANK_PATH, expected_size=GUI_ANNOTATION_GENBANK_SIZE,
        expected_sha256=GUI_ANNOTATION_GENBANK_SHA256,
        expected_id='NC_001879.2', expected_length=155943,
    )
    provenance = json.loads((Path(__file__).parents[2] / 'joint-source-verification.json').read_text())
    assert_fixture_identity(GUI_ANNOTATION_GENBANK_PATH,
                            expected_size=provenance['mirror_size'],
                            expected_sha256=provenance['mirror_sha256'])
    _validate_annotation_fixture()
    assert_output_paths(output_paths, NAMES, 'H-GUI-16')
    download_dir.mkdir(parents=True, exist_ok=True)
    capture = open_browser_capture(browser_type, base_url, device_scale_factor=2)
    page = capture.page
    page.on('console', lambda message: print('Capture console:', message.text) if message.type == 'error' else None)
    page.on('pageerror', lambda error: print('Capture page error:', str(error)))
    try:
        page.goto(base_url, wait_until='domcontentloaded')
        wait_for_app_shell(page)
        page.get_by_label('GenBank/DDBJ File', exact=True).set_input_files(GUI_ANNOTATION_GENBANK_PATH)
        _configure_gallery_presentation(page)
        page.get_by_label('Separate Strands', exact=True).uncheck()
        page.get_by_label('Resolve Overlaps', exact=True).check()
        annotations = page.get_by_label('Region Annotations', exact=True)
        annotations.click()
        page.get_by_label('Import TSV', exact=True).set_input_files(GUI_ANNOTATION_TABLE_PATH)
        slots = _configure_gallery_slots(page)
        generate_and_inspect(page, _inspect_tracks_svg, _assert_safe_svg)
        start = page.get_by_role('spinbutton', name='Display start NC_001879.2 #1', exact=True)
        start.fill('5500')
        start.press('Tab')
        expect(start).to_have_value('5500')
        # Existing row carries the biological record identity; no accessible row name.
        row = start.locator('xpath=ancestor::div[@data-record-rotation]')
        row.screenshot(path=str(output_paths[NAMES[0]]))
        tolerance = page.get_by_role('spinbutton', name='Feature overlap tolerance (bp)', exact=True)
        tolerance.fill('1')
        tolerance.press('Tab')
        _fit_circular_preview(page, target_zoom='70%', pan_left_ratio=0.0)
        feature = page.evaluate('''() => window.__GBDRAW_APP__.extractedFeatures.find(f =>
            JSON.stringify(f.qualifiers || f).includes('NP_054479.1'))''')
        assert feature, 'rps16 source feature missing'
        # The existing drawer toggle and search input lack accessible names.
        page.locator('.drawer-toggle').click()
        page.get_by_placeholder('Search by feature or annotation...').fill('ribosomal protein S16')
        edit = page.locator('.right-drawer').get_by_role('button', name='Edit', exact=True)
        expect(edit).to_have_count(1)
        edit.click()
        placement = page.get_by_role('combobox', name='Feature placement', exact=True)
        expect(placement).to_be_visible()
        placement.select_option('outward')
        expect(placement).to_have_value('outward')
        dialog = page.get_by_role('dialog', name='Feature details: ribosomal protein S16', exact=True)
        expect(dialog).to_be_visible()
        box = dialog.bounding_box()
        placement_box = placement.bounding_box()
        assert box and placement_box
        # Keep the source identity, control and Generate instruction in one operation crop.
        page.screenshot(path=str(output_paths[NAMES[1]]), clip={
            **box, 'height': placement_box['y'] + placement_box['height'] + 26 - box['y']})
        page.get_by_role('button', name='Close feature popup', exact=True).click()
        page.locator('.drawer-toggle').click()
        final = generate_and_inspect(page, _inspect_tracks_svg, _assert_safe_svg)
        values = page.evaluate('''() => ({request: window.__GBDRAW_APP__.lastSuccessfulRunRequest,
            svg: window.__GBDRAW_APP__.results[0].content})''')
        (download_dir / 'joint.svg').write_text(values['svg'])
        # The next action remains Save/Load; formal joint browser tests exercise fresh Load and replay.
        page.on('dialog', lambda dialog: (print('Capture dialog:', dialog.type, dialog.message), dialog.accept('Rotated placed plastome')))
        with page.expect_download() as pending:
            page.get_by_role('button', name='Save Session', exact=True).click()
        session = pending.value
        session_path = download_dir / session.suggested_filename
        session.save_as(session_path)
        saved = json.loads(gzip.decompress(session_path.read_bytes()))
        assert saved['version'] == 42 and saved['renderRequest']['schema'] == 7
        assert saved['renderRequest']['records'][0]['display']['startCoordinate'] == 5500
        assert saved['renderRequest']['diagramOptions']['configOverrides']['canvas.feature_overlap_tolerance_bp'] == 1
        assert len(saved['renderRequest']['diagramOptions']['featurePlacements']) == 1
        values['request'] = saved['renderRequest']
        (download_dir / 'joint-controls.json').write_text(json.dumps({'source': report, 'feature': feature, 'values': values}, indent=2))
        capture.assert_clean()
        return CaptureResult(
            screenshot_bytes={name: output_paths[name].stat().st_size for name in NAMES},
            final_svg_semantics=final, download={'filename': session.suggested_filename, 'bytes': (download_dir / session.suggested_filename).stat().st_size},
            fixture_report=report, track_slots=slots,
        )
    finally:
        capture.close()
