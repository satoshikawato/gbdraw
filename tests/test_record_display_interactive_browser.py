"""Real Preview binder and standalone actions on typed-plan display fragments."""
from __future__ import annotations

import functools
import json
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest
from playwright.sync_api import sync_playwright

from gbdraw.linear_comparison import LinearComparison
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.render.interactive_svg import enrich_svg
from tests.test_record_display_comparisons import hit_frame, record, request, svg


pytestmark = pytest.mark.browser
ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def source_server():
    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, *_args):
            pass
    server = ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Handler, directory=str(ROOT)))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{server.server_port}"
    server.shutdown()
    server.server_close()
    thread.join()


def interactive_fixture(mode, reverse, start=41):
    items = [record(), record(120, "subject")]
    options = ({"linear_comparisons": [LinearComparison(0, 1, hit_frame())]} if mode == "linear" else
               {"conservation_dataframes": [hit_frame()], "conservation_reference": "query"})
    plan, root = svg(request(mode, items if mode == "linear" else items[:1],
                             [start, 51] if mode == "linear" else [start],
                             [reverse, False] if mode == "linear" else [reverse], **options))
    context = build_interactive_svg_context(
        plan.records, mode=mode, linear_rendered_feature_ids=mode == "linear",
        comparison_sequence_records=[items[1:]] if mode == "circular" else (),
    )
    standalone = enrich_svg(ET.tostring(root, encoding="unicode"), context)
    catalog = json.loads(next(n.text for n in ET.fromstring(standalone).iter()
                              if n.get("id") == "gbdraw-interactive-feature-metadata"))
    return ET.tostring(root, encoding="unicode"), standalone, catalog


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("reverse", [False, True])
def test_display_fragments_preview_and_standalone_source_actions(source_server, tmp_path, mode, reverse):
    source, standalone, catalog = interactive_fixture(mode, reverse)
    (tmp_path / "source.svg").write_text(source)
    standalone_path = tmp_path / "standalone.svg"
    standalone_path.write_text(standalone)
    (tmp_path / "catalog.json").write_text(json.dumps(catalog, indent=2))
    expected_count = len([n for n in ET.fromstring(source).iter() if n.get("data-gbdraw-match-id")])
    # Independent source oracle: raw HSP inclusive endpoints, existing materialized RC sequence.
    query_sequence = str(record().seq.reverse_complement() if reverse else record().seq)[20:60]
    subject_sequence = str(record(120, "subject").seq)[30:90]
    requests = []
    errors = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        context = browser.new_context(viewport={"width": 1440, "height": 1000}, accept_downloads=True)
        def route(request_route):
            url = request_route.request.url
            requests.append(url)
            if url.startswith(source_server) or url.startswith("file:"):
                request_route.continue_()
            else:
                request_route.abort()
        context.route("**/*", route)
        page = context.new_page()
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(source_server)
        page.add_script_tag(url=f"{source_server}/gbdraw/web/vendor/dompurify/purify.min.js")
        page.evaluate("""async ({source, catalog, mode}) => {
          const {createFeatureSvgActions} = await import('/gbdraw/web/js/app/feature-editor/svg-actions.js');
          const {admitFeatureCatalog} = await import('/gbdraw/web/js/services/feature-catalog.js');
          const {createSequenceSourceRegistry} = await import('/gbdraw/web/js/app/match-sequences.js');
          const {copyTextToClipboard} = await import('/gbdraw/web/js/utils/clipboard.js');
          const {downloadTextFile} = await import('/gbdraw/web/js/services/text-download.js');
          const ref = value => ({value});
          const results = [{name: 'interactive.svg', content: source}];
          const admission = admitFeatureCatalog(catalog, results, {mode});
          document.body.innerHTML = '<div id="preview"></div>';
          const host = document.querySelector('#preview');
          host.innerHTML = DOMPurify.sanitize(source, {USE_PROFILES: {svg: true, svgFilters: true}});
          const state = {
            results: ref(results), selectedResultIndex: ref(0), orthogroups: ref([]), collinearGroups: ref([]),
            orthogroupNameOverrides: {}, orthogroupDescriptionOverrides: {},
            extractedFeatures: ref(admission.featureState.extractedFeatures),
            biologicalFeatures: ref(admission.featureState.biologicalFeatures),
            featuresBySvgId: ref(new Map(admission.featureState.extractedFeatures.map(feature => [feature.svg_id, feature]))),
            featureColorOverrides: {}, featureVisibilityOverrides: {}, svgContainer: ref(host),
            clickedFeature: ref(null), clickedFeaturePos: {}, clickedPairwiseMatch: ref(null), clickedPairwiseMatchPos: {},
            matchSequenceRegistry: createSequenceSourceRegistry(catalog.items[0].sequenceSources),
            selectedAnnotation: ref(null), featurePopupSize: {}, featureSelectionDrag: {active: false},
            skipCaptureBaseConfig: ref(false), adv: {rich_feature_popup: true},
          };
          window.fixtureState = state;
          window.fixtureActions = createFeatureSvgActions({state, getFeatureColor: () => '#123456', getEffectiveLegendCaption: () => ''});
          window.fixtureActions.attachSvgFeatureHandlers();
          window.fixtureClipboard = [];
          Object.defineProperty(navigator, 'clipboard', {value: {writeText: async text => window.fixtureClipboard.push(text)}, configurable: true});
          window.fixtureCopy = copyTextToClipboard;
          window.fixtureDownload = downloadTextFile;
          window.fixtureMount = content => {
            results[0].content = content;
            host.innerHTML = DOMPurify.sanitize(content, {USE_PROFILES: {svg: true, svgFilters: true}});
            window.fixtureActions.attachSvgFeatureHandlers();
          };
        }""", {"source": source, "catalog": catalog, "mode": mode})
        paths = page.locator('[data-gbdraw-match-id]')
        assert paths.count() == expected_count
        copied = []
        for index in range(expected_count):
            paths.nth(index).dispatch_event("mouseover")
            assert page.locator('[data-gbdraw-match-id][data-gbdraw-hover-opacity]').count() == expected_count
            paths.nth(index).dispatch_event("click")
            payload = page.evaluate("fixtureState.clickedPairwiseMatch.value")
            assert payload is not None
            bundle = payload["sequenceBundle"]
            assert len(bundle["entries"]) == 2
            assert all(entry["available"] for entry in bundle["entries"])
            fasta = bundle["combinedFasta"]
            assert fasta.count(">") == 2
            bodies = ["".join(part.splitlines()[1:]) for part in fasta.split(">")[1:]]
            assert bodies == [query_sequence, subject_sequence]
            await_expression = "fixtureCopy(fixtureState.clickedPairwiseMatch.value.sequenceBundle.combinedFasta)"
            page.evaluate(await_expression)
            assert page.evaluate("fixtureClipboard.at(-1)") == fasta
            with page.expect_download() as download_info:
                page.evaluate("fixtureDownload('match.fna', fixtureState.clickedPairwiseMatch.value.sequenceBundle.combinedFasta)")
            downloaded = download_info.value
            assert Path(downloaded.path()).read_text() == fasta
            copied.append(fasta)
        assert len(set(copied)) == 1
        page.screenshot(path=str(tmp_path / "preview.png"), full_page=False)

        feature_paths = page.locator('path[data-gbdraw-feature-id][fill="#54bcf8"]')
        feature_id = feature_paths.first.get_attribute('data-gbdraw-feature-id')
        feature_paths = page.locator(f'path[data-gbdraw-feature-id="{feature_id}"]:not([fill="none"])')
        assert feature_paths.count() == 2
        feature_payloads = []
        for index in range(feature_paths.count()):
            feature_paths.nth(index).dispatch_event('click')
            payload = page.evaluate('fixtureState.clickedFeature.value')
            fasta = payload['nucleotideFasta']
            assert ''.join(fasta.splitlines()[1:]) == str(record().features[0].extract(record().seq))
            assert fasta.count('>') == 1
            page.evaluate('fixtureCopy(fixtureState.clickedFeature.value.nucleotideFasta)')
            assert page.evaluate('fixtureClipboard.at(-1)') == fasta
            with page.expect_download() as download_info:
                page.evaluate("fixtureDownload('feature.fna', fixtureState.clickedFeature.value.nucleotideFasta)")
            assert Path(download_info.value.path()).read_text() == fasta
            feature_payloads.append(payload)
        assert feature_payloads[0] == feature_payloads[1]

        b_source, _, _ = interactive_fixture(mode, reverse, start=11)
        for content in (source, b_source, source):
            page.evaluate('fixtureMount', content)
            current_paths = page.locator('[data-gbdraw-match-id]')
            for index in range(current_paths.count()):
                current_paths.nth(index).dispatch_event('mouseover')
                assert page.locator('[data-gbdraw-match-id][data-gbdraw-hover-opacity]').count() == current_paths.count()
                current_paths.nth(index).dispatch_event('click')
                assert page.evaluate('fixtureState.clickedPairwiseMatch.value.sequenceBundle.combinedFasta') == copied[0]
        assert page.locator('[data-gbdraw-match-id]').count() == expected_count

        # A separate document runs the exact packaged standalone runtime.
        page.goto(standalone_path.as_uri())
        page.evaluate("""() => {
          window.fixtureClipboard = [];
          Object.defineProperty(navigator, 'clipboard', {value: {writeText: async text => window.fixtureClipboard.push(text)}, configurable: true});
        }""")
        standalone_paths = page.locator('[data-gbdraw-interactive-match="true"]')
        assert standalone_paths.count() == expected_count
        for index in range(expected_count):
            standalone_paths.nth(index).dispatch_event("click")
            assert page.locator('.gbdraw-interactive-pairwise-match--selected').count() == expected_count
            popup = page.locator('#gbdraw-feature-popup')
            assert popup.is_visible()
            assert '21' in popup.text_content() and '60' in popup.text_content()
            button = popup.locator('.gfi-block-actions').filter(has_text="Both spans").get_by_role("button", name="Copy", exact=True)
            button.click()
            page.wait_for_function("fixtureClipboard.length > 0")
            fasta = page.evaluate("fixtureClipboard.at(-1)")
            with page.expect_download() as download_info:
                popup.locator('.gfi-block-actions').filter(has_text="Both spans").get_by_role("button", name="FASTA", exact=True).click()
            assert Path(download_info.value.path()).read_text() == fasta
            assert fasta.count(">") == 2
            assert ["".join(part.splitlines()[1:]) for part in fasta.split(">")[1:]] == [query_sequence, subject_sequence]
        feature_paths = page.locator(f'path[data-gbdraw-feature-id="{feature_id}"]:not([fill="none"])')
        assert feature_paths.count() == 2
        for index in range(feature_paths.count()):
            feature_paths.nth(index).dispatch_event('click')
            popup = page.locator('#gbdraw-feature-popup')
            popup.locator('[data-tab="sequence"]').click()
            block = popup.locator('.gfi-block').filter(has=page.locator('.gfi-block-title', has_text='Nucleotide'))
            block.locator('.gfi-copy').filter(has_text='Copy').click()
            fasta = page.evaluate('fixtureClipboard.at(-1)')
            assert ''.join(fasta.splitlines()[1:]) == str(record().features[0].extract(record().seq))
            assert fasta.count('>') == 1
            with page.expect_download() as download_info:
                block.locator('[data-download-index]').click()
            assert Path(download_info.value.path()).read_text() == fasta
        page.screenshot(path=str(tmp_path / "standalone.png"), full_page=False)
        page.set_viewport_size({"width": 390, "height": 844})
        standalone_paths.first.dispatch_event('click')
        assert page.locator('#gbdraw-feature-popup').is_visible()
        page.screenshot(path=str(tmp_path / "standalone-narrow.png"), full_page=False)
        assert not errors
        (tmp_path / "browser-evidence.json").write_text(json.dumps({
            "mode": mode, "reverse": reverse, "browser": browser.version,
            "logical_matches": 1, "fragments": expected_count, "preview_copy_payloads": copied, "preview_feature_payloads": feature_payloads,
            "requests": requests, "errors": errors,
        }, indent=2))
        browser.close()


@pytest.mark.parametrize("conflict", ["endpoint", "fragment_index", "feature_claim"])
def test_standalone_rejects_conflicting_fragment_or_source_identity(tmp_path, conflict):
    _, standalone, _ = interactive_fixture("linear", False)
    root = ET.fromstring(standalone)
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    if conflict == "endpoint":
        paths[1].set("data-qstart", "22")
    elif conflict == "fragment_index":
        paths[1].set("data-gbdraw-match-fragment", "0")
    else:
        metadata = next(n for n in root.iter() if n.get("id") == "gbdraw-interactive-feature-metadata")
        catalog = json.loads(metadata.text)
        catalog["items"][0]["comparisonMatches"][0]["queryBiologicalFeatureId"] = "missing-feature"
        metadata.text = json.dumps(catalog)
    path = tmp_path / "conflicting.svg"
    path.write_text(ET.tostring(root, encoding="unicode"))
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page()
        page.goto(path.as_uri())
        page.locator('[data-gbdraw-match-id]').first.dispatch_event('click')
        if conflict == "feature_claim":
            popup = page.locator('#gbdraw-feature-popup')
            assert 'Match feature endpoint identity is invalid.' in popup.text_content()
            assert popup.locator('.gfi-block-actions').filter(has_text="Both spans").count() == 0
        else:
            assert page.locator('.gbdraw-interactive-pairwise-match--selected').count() == 0
            assert not page.locator('#gbdraw-feature-popup').is_visible()
        browser.close()


def test_circular_grid_match_actions_resolve_duplicate_record_instances(tmp_path):
    from dataclasses import replace
    from gbdraw.api.options import CircularMultiRecordOptions
    items = [record(), record()]
    req = request("circular", items, [41, 41], [False, True],
                  conservation_dataframes=[hit_frame()], conservation_reference="query")
    req = replace(req, grouping="grid", layout=CircularMultiRecordOptions())
    plan, root = svg(req)
    context = build_interactive_svg_context(plan.records, mode="circular", record_transforms=plan.transforms, comparison_sequence_records=[[record(120, "subject")]])
    enriched = enrich_svg(ET.tostring(root, encoding="unicode"), context)
    path = tmp_path / "grid.svg"
    path.write_text(enriched)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        page = browser.new_page()
        page.goto(path.as_uri())
        page.evaluate("""() => {
          window.copied = [];
          Object.defineProperty(navigator, 'clipboard', {value: {writeText: async text => copied.push(text)}, configurable: true});
        }""")
        paths = page.locator('[data-gbdraw-interactive-match="true"]')
        assert paths.count() == 4
        for index in range(paths.count()):
            element = paths.nth(index)
            record_index = int(element.get_attribute('data-query-record-index'))
            element.dispatch_event('click')
            assert page.locator('.gbdraw-interactive-pairwise-match--selected').count() == 2
            page.locator('.gfi-block-actions').filter(has_text="Both spans").get_by_role('button', name='Copy', exact=True).click()
            fasta = page.evaluate('copied.at(-1)')
            expected = str(items[0].seq.reverse_complement() if record_index else items[0].seq)[20:60]
            assert ''.join(fasta.split('>')[1].splitlines()[1:]) == expected
            assert fasta.count('>') == 2
        page.screenshot(path=str(tmp_path / 'grid.png'))
        browser.close()
