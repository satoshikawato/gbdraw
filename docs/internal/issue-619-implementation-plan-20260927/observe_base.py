import argparse
import functools
import hashlib
import json
import threading
import subprocess
from importlib.metadata import version
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from playwright.sync_api import sync_playwright

parser = argparse.ArgumentParser(
    description="Observe frozen issue #619 baseline without modifying runtime."
)
parser.add_argument("--source-root", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
ROOT = args.source_root.resolve()
BASE_SHA = "f5f86634459e0dcd46c1a452e9219fbba635d429"
repo_root = Path(__file__).resolve().parents[3]
runtime_sources = [
    "gbdraw/web/index.html",
    "gbdraw/web/js/app/circular-track-slots.js",
    "gbdraw/web/js/app/track-slot-display.js",
    "gbdraw/web/js/app/app-setup.js",
    "gbdraw/web/js/services/session-request.js",
    "gbdraw/web/js/services/config.js",
]
source_hashes = {}
for source in runtime_sources:
    actual = (ROOT / source).read_bytes()
    expected = subprocess.check_output(
        ["git", "-C", str(repo_root), "show", f"{BASE_SHA}:{source}"]
    )
    assert actual == expected, f"Source differs from frozen base: {source}"
    source_hashes[source] = hashlib.sha256(actual).hexdigest()


class QuietHandler(SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


server = ThreadingHTTPServer(
    ("127.0.0.1", 0), functools.partial(QuietHandler, directory=str(ROOT))
)
threading.Thread(target=server.serve_forever, daemon=True).start()
url = f"http://127.0.0.1:{server.server_port}/gbdraw/web/index.html"
observations = []
try:
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        for width in (1440, 390):
            page = browser.new_page(viewport={"width": width, "height": 900})
            page.set_default_timeout(180_000)
            page_errors = []
            console_errors = []
            page.on(
                "console",
                lambda e: console_errors.append(e.text) if e.type == "error" else None,
            )
            page.on("pageerror", lambda e: page_errors.append(str(e)))
            page.route(
                "**/*",
                lambda route: (
                    route.continue_()
                    if route.request.url.startswith(
                        ("http://127.0.0.1:", "data:", "blob:")
                    )
                    else route.abort()
                ),
            )
            page.add_init_script(
                """window.__issue619Workers=[]; const NativeWorker=window.Worker; window.Worker=new Proxy(NativeWorker,{construct(target,args){window.__issue619Workers.push(String(args[0]));return Reflect.construct(target,args);}});"""
            )
            page.goto(url)
            page.wait_for_function(
                "window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions).length>0"
            )
            page.locator(
                'input[type=file][accept=".json,.gz,application/json,application/gzip"]'
            ).set_input_files(
                ROOT
                / "gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json"
            )
            try:
                page.wait_for_function(
                    "window.__GBDRAW_APP__.adv.circular_track_slots.some(s=>s.id==='plastome_regions') && !window.__GBDRAW_APP__.sessionImportPending"
                )
            except Exception:
                print(
                    page.evaluate(
                        "JSON.stringify({pending:window.__GBDRAW_APP__.sessionImportPending,error:window.__GBDRAW_APP__.errorLog,slots:window.__GBDRAW_APP__.adv.circular_track_slots.map(s=>s.id),workers:window.__issue619Workers})"
                    )
                )
                print(console_errors, page_errors)
                raise
            page.locator("summary[aria-label=Layout]").click()
            page.get_by_title("Open Custom Track Slots", exact=True).click()
            page.locator("input[title=Width]").first.wait_for(state="visible")
            before = page.evaluate(
                """() => ({slots:JSON.parse(JSON.stringify(window.__GBDRAW_APP__.adv.circular_track_slots)),results:JSON.parse(JSON.stringify(window.__GBDRAW_APP__.results)),undo:window.__GBDRAW_HISTORY__.getUndoCount(),workers:window.__issue619Workers.slice()})"""
            )
            fields = page.evaluate(
                """() => [...document.querySelectorAll('input[title="Width"],input[title="Radius"]')].map(el=>({title:el.title,text:el.value,suffix:el.parentElement.querySelector('.track-slot-unit-suffix')?.textContent.trim()||'',auto:el.parentElement.querySelector('.track-slot-auto-placeholder')?.textContent.trim()||'',visible:!!el.getClientRects().length}))"""
            )
            page.locator("input[title=Width]").nth(1).focus()
            page.locator("input[title=Radius]").nth(1).focus()
            page.locator("input[title=Width]").first.focus()
            after = page.evaluate(
                """() => ({slots:JSON.parse(JSON.stringify(window.__GBDRAW_APP__.adv.circular_track_slots)),results:JSON.parse(JSON.stringify(window.__GBDRAW_APP__.results)),undo:window.__GBDRAW_HISTORY__.getUndoCount(),workers:window.__issue619Workers})"""
            )
            assert before["slots"] == after["slots"]
            assert before["results"] == after["results"]
            assert before["undo"] == after["undo"]
            assert before["workers"] == after["workers"]
            assert sum(f["text"] == "[object Object]" for f in fields) == 4
            assert not page_errors, page_errors
            observations.append(
                {
                    "viewportWidth": width,
                    "fields": fields,
                    "slotsPreservedOnFocusBlur": True,
                    "resultPreservedOnFocusBlur": True,
                    "undoBefore": before["undo"],
                    "undoAfter": after["undo"],
                    "workerUrls": after["workers"],
                    "pageErrors": page_errors,
                    "consoleErrors": console_errors,
                }
            )
            page.close()
        browser.close()
finally:
    server.shutdown()
    server.server_close()
result = {
    "baseSha": BASE_SHA,
    "runtimeSourceMatchesBase": True,
    "runtimeSourceSha256": source_hashes,
    "playwrightVersion": version("playwright"),
    "browserVersion": browser.version,
    "wheelSha256": hashlib.sha256(
        (ROOT / "gbdraw/web/gbdraw-0.14.0-py3-none-any.whl").read_bytes()
    ).hexdigest(),
    "fixtureSha256": hashlib.sha256(
        (
            ROOT / "gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json"
        ).read_bytes()
    ).hexdigest(),
    "observations": observations,
    "limits": [
        "No Generate/Save/Load round trip or full-app Undo/Redo after a fix was performed.",
        "Only the frozen baseline was observed; numeric/unit controls were not implemented.",
    ],
}
args.output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
