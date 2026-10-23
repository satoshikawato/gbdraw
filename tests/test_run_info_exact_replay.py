"""Issue #469: execute the helper downloaded after the original History sequence."""
from __future__ import annotations

import base64
import functools
import gzip
import hashlib
import http.server
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import threading
import xml.etree.ElementTree as ET
import zipfile

import pytest

from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA
from tests.utils.svg_compare import compare_svgs

pytestmark = pytest.mark.browser
ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/sessions/HmmtDNA_basic_circular.issue-469.json.gz"
SOURCE_SHA = "f2e922c26561d37a8d3b410ba972526c709cad5635f598bbdcfc234f02642c92"
OUTPUT_SHA = "e0503332f9a48545417cc75da80e57fd50af7726aa1fb46b6985995a0eb81b90"
# The replay before the feature label style (dev d84352d8a, hash
# 76e2e0a1250b90d63fcd55ef3a96119c8ab3fafc5c890594033fa27aa94379d0), with only
# the overlayObstacles array of the composition removed. The leader width went
# from 1.5 to 1.0, which shrinks each obstacle box by half that change per
# side, so the boxes cannot be peeled and are left out of the comparison.
PRE_FEATURE_LABEL_STYLE_COMPARABLE_SHA = (
    "c42139da6f913db165da9bfaa5de8f42e694b81ab8d7b8296570f647932c425a"
)


def _without_feature_label_style(svg: bytes) -> bytes:
    """Peel what the Circular feature label style adds or changes."""
    svg = re.sub(rb'<use aria-hidden="true"[^>]*xlink:href="#label_text_[0-9a-f]+" />', b"", svg)
    svg = re.sub(rb' id="label_text_[0-9a-f]+"', b"", svg)
    svg = re.sub(rb'(<line\b[^>]*?) stroke-width="1\.0"', rb'\1 stroke-width="1.5"', svg)
    return re.sub(
        rb'&quot;overlayObstacles&quot;:\[.*?\],&quot;overlayPolicy',
        b"&quot;overlayPolicy",
        svg,
        count=1,
    )


def test_downloaded_exact_replay_after_original_history(tmp_path: Path) -> None:
    # Keep browser imports inside the browser test for non-browser collection.
    from playwright.sync_api import sync_playwright

    raw = gzip.decompress(FIXTURE.read_bytes())
    assert hashlib.sha256(raw).hexdigest() == (
        "69786dd18f7a431441085fad3e7d20ed9b54ed2f0b5aeac9685dbb3a704e6da1"
    )
    original = json.loads(raw)
    original["title"] = "lazy-normal"
    original["resources"]["unused-lazy-contract"] = {
        "kind": "web-file", "name": "unused-lazy-contract.txt", "type": "text/plain",
        "size": 7, "lastModified": 0, "encoding": "base64", "data": "dW51c2VkCg==",
    }
    source = base64.b64decode(original["resources"]["record-1-genbank"]["data"])
    assert hashlib.sha256(source).hexdigest() == SOURCE_SHA
    fixture = tmp_path / "input.json"
    fixture.write_text(json.dumps(original), encoding="utf-8")
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(ROOT))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    origin = f"http://127.0.0.1:{server.server_port}"
    snapshots = {}
    external_requests = []
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            page = browser.new_page(viewport={"width": 1900, "height": 1100})
            page.set_default_timeout(180_000)
            page.on("dialog", lambda dialog: dialog.dismiss())

            def route_request(route):
                if route.request.url.startswith(origin + "/"):
                    route.continue_()
                else:
                    external_requests.append(route.request.url)
                    route.abort()

            page.route("**/*", route_request)
            page.goto(origin + "/gbdraw/web/index.html")
            page.wait_for_function("Boolean(window.__GBDRAW_APP__)")
            page.locator(
                'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
            ).set_input_files(fixture)
            page.wait_for_function("""() => {
                const app = window.__GBDRAW_APP__;
                return app.sessionTitle === 'lazy-normal' && app.results.length > 0
                    && app.sessionImportPending === false;
            }""")

            def snapshot(name):
                value = page.evaluate("""async () => {
                    const a = window.__GBDRAW_APP__, h = window.__GBDRAW_HISTORY__;
                    const { getCommittedCanonicalSession } = await import('./js/services/config.js');
                    return { info: a.lastRunInfo, result: a.results[a.selectedResultIndex],
                        selectedResultIndex: a.selectedResultIndex, undo: h.getUndoCount(),
                        redo: h.getRedoCount(), committed: getCommittedCanonicalSession() };
                }""")
                snapshots[name] = value
                (tmp_path / f"{name}.svg").write_text(value["result"]["content"], encoding="utf-8")
                return value

            snapshot("loaded")
            assert page.evaluate("window.__GBDRAW_APP__.runAnalysis()") == {"status": "ok"}
            generated = snapshot("first")
            assert page.evaluate("window.__GBDRAW_HISTORY__.undo()") is True
            snapshot("undo-first")
            assert page.evaluate("window.__GBDRAW_HISTORY__.redo()") is True
            snapshot("redo-first")
            # Exact operation/arguments from the retained Issue diagnostic and the
            # current lazy-materialization contract; this does not edit the legend.
            width = page.evaluate("""async () => {
                const m = await import('./js/services/diagram-generation.js');
                const helper = await m.runDiagramHelperOperation(
                    m.DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES,
                    { features: [], rules: [], kind: 'color' });
                return Array.isArray(helper.result.matches);
            }""")
            assert width is True
            assert page.evaluate("window.__GBDRAW_APP__.runAnalysis()") == {"status": "ok"}
            second = snapshot("second")
            assert second["committed"] == generated["committed"]
            edit = page.evaluate("""async () => {
                const a = window.__GBDRAW_APP__;
                const feature = a.extractedFeatures[0];
                const color = String(a.getFeatureColorValue(feature) || '').toLowerCase()
                    === '#123456' ? '#654321' : '#123456';
                return { color, applied: await a.setFeatureColorValue(feature, color) };
            }""")
            assert edit == {"color": "#123456", "applied": True}
            snapshot("edited")
            for name, action in [("undo-edit", "undo"), ("undo-generate", "undo"),
                                 ("redo-generate", "redo"), ("redo-edit", "redo")]:
                assert page.evaluate(f"window.__GBDRAW_HISTORY__.{action}()") is True
                snapshot(name)
            selected = snapshots["redo-edit"]
            assert selected["committed"] == generated["committed"]
            assert (selected["undo"], selected["redo"], selected["selectedResultIndex"]) == (2, 0, 0)
            page.get_by_role("button", name=re.compile("Run info", re.I)).click()
            info = selected["info"]
            for kind in ("sourceRecipe", "exactReplay"):
                assert info[kind]["command"] in page.locator("pre").all_text_contents()
            with page.expect_download() as download:
                page.get_by_role("button", name=re.compile("Download reproducibility files")).click()
            archive = tmp_path / "helpers.zip"
            download.value.save_as(archive)
            page.screenshot(path=str(tmp_path / "run-info.png"))
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join()
    assert external_requests == []
    (tmp_path / "history.json").write_text(json.dumps(snapshots), encoding="utf-8")
    with zipfile.ZipFile(archive) as bundle:
        assert bundle.namelist() == ["out.gbdraw-session.json"]
        session = json.loads(bundle.read("out.gbdraw-session.json"))
    assert session["version"] == CURRENT_SESSION_VERSION
    assert session["renderRequest"]["schema"] == CANONICAL_REQUEST_SCHEMA
    assert session["renderRequest"] == generated["committed"]["renderRequest"]
    assert session["resources"] == generated["committed"]["resources"]
    assert session["results"] == [generated["result"]]
    assert list(session["resources"]) == ["record-1-genbank"]
    assert hashlib.sha256(base64.b64decode(
        session["resources"]["record-1-genbank"]["data"]
    )).hexdigest() == SOURCE_SHA

    env = {key: value for key, value in os.environ.items() if key not in ("PYTHONPATH", "PYTHONHOME")}
    for kind in ("exactReplay", "sourceRecipe"):
        work = tmp_path / kind
        work.mkdir()
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(work)
        if kind == "sourceRecipe":
            # The UI explicitly requires the original upload for Source recipe.
            (work / "HmmtDNA.gbk").write_bytes(source)
        command = shlex.split(info[kind]["command"])
        result = subprocess.run(command, cwd=work, env=env, capture_output=True, text=True)
        (work / "cli.log").write_text(result.stdout + result.stderr, encoding="utf-8")
        assert result.returncode == 0, result.stderr
        output = work / "out.svg"
        assert output.stat().st_size > 0
        assert hashlib.sha256(output.read_bytes()).hexdigest() == OUTPUT_SHA
        assert hashlib.sha256(_without_feature_label_style(output.read_bytes())).hexdigest() == (
            PRE_FEATURE_LABEL_STYLE_COMPARABLE_SHA
        )
        # DOMPurify removes only this non-rendering SVG profile declaration.
        # Retain every geometry, style, label, and composition attribute.
        tree = ET.parse(output)
        assert tree.getroot().attrib.pop("baseProfile") == "full"
        comparable = work / "comparable.svg"
        tree.write(comparable, encoding="unicode")
        comparison = compare_svgs(tmp_path / "first.svg", comparable)
        assert comparison.equal, comparison.differences
