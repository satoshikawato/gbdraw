"""S07 real Collinear Generate/export/save/reload, using the current local wheel."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse

from playwright.sync_api import sync_playwright
from tests import run_losat_cache_browser_acceptance as acceptance

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    fixture = ROOT / "gbdraw/web/gallery/sessions/hepatoplasmataceae_collinear.gbdraw-session.json.gz"
    wheel, = (ROOT / "gbdraw/web").glob("gbdraw-*-py3-none-any.whl")
    report = {"fixtureSha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
              "wheelSha256": hashlib.sha256(wheel.read_bytes()).hexdigest(), "runs": []}
    artifacts = ROOT / ".venv/s07-browser-artifacts"
    artifacts.mkdir(parents=True, exist_ok=True)
    with acceptance._serve_repo() as base, sync_playwright() as p:
        browser = p.chromium.launch()
        try:
            for width, height in ((1280, 720), (390, 844)):
                checks = acceptance.AcceptanceChecks()
                run = {"viewport": [width, height], "browser": browser.version,
                       "externalRequests": [], "pageErrors": [], "localPaths": set(), "generations": []}
                context = browser.new_context(accept_downloads=True, viewport={"width": width, "height": height})
                def route(r):
                    url = urlparse(r.request.url)
                    if url.hostname not in {"127.0.0.1", "localhost"}:
                        run["externalRequests"].append(r.request.url)
                        r.abort()
                    else:
                        run["localPaths"].add(url.path)
                        r.continue_()
                context.route("**/*", route)
                page = context.new_page()
                page.on("pageerror", lambda error: run["pageErrors"].append(str(error)))
                page.set_default_timeout(120_000)
                page.add_init_script("""window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
                  window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
                    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__++;
                    throw new Error('S07 must reuse the verified Gallery raw evidence');
                  };""")
                page.goto(base + "/gbdraw/web/index.html", wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                acceptance._import_session(page, fixture, checks)
                page.evaluate("""async () => {
                  const app = window.__GBDRAW_APP__;
                  app.losat.blastp.mode = 'collinear';
                  await Vue.nextTick();
                  app.losat.blastp.collinearInferOrthogroups = true;
                }""")
                def generate():
                    checks.require(page.evaluate("() => window.__GBDRAW_APP__.losat.blastp.mode") == "collinear", "Wrong active analysis mode")
                    page.wait_for_function("() => Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length > 0")
                    observed = acceptance._generate(page)
                    checks.require(observed["result"] == {"status": "ok"}, str(observed))
                    checks.require(observed["executorCalls"] == 0, "Unexpected new raw search")
                    acceptance._settle_app_render(page)
                    exported = acceptance._download_svg(page, checks)
                    acceptance._assert_svg_geometry_parity(page, exported, checks)
                    checks.require(page.evaluate("""async () => {
                      const { state } = await import('/gbdraw/web/js/state.js');
                      return Array.from(state.svgContainer.value.querySelectorAll('[data-collinearity-block-id]'))
                        .some(element => element.getAttribute('data-collinearity-block-id'));
                    }"""), "No nonempty Collinear block IDs rendered")
                    geometry = acceptance._inspect_layout(page)["geometrySignature"]
                    observed["geometrySha256"] = hashlib.sha256(geometry.encode()).hexdigest()
                    observed["svgSha256"] = hashlib.sha256(exported.read_bytes()).hexdigest()
                    run["generations"].append(observed)
                    return geometry, exported.read_bytes()
                first, svg = generate()
                (artifacts / f"{width}.svg").write_bytes(svg)
                checks.require(generate()[0] == first, "Repeated Generate changed geometry")
                saved = acceptance._save_session(page, checks)
                saved_data = saved.read_bytes()
                saved_path = artifacts / f"{width}.gbdraw-session.json"
                saved_path.write_bytes(saved_data)
                run["sessionSha256"] = hashlib.sha256(saved_data).hexdigest()
                page.reload(wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                acceptance._import_session(page, saved_path, checks)
                checks.require(generate()[0] == first, "Save/reload/regeneration changed geometry")
                page.screenshot(path=str(artifacts / f"{width}.png"), full_page=True)
                checks.require(not run["externalRequests"] and not run["pageErrors"], str(run))
                run["localPaths"] = sorted(run["localPaths"])
                run["assertions"] = checks.count
                report["runs"].append(run)
                context.close()
        finally:
            browser.close()
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"runs": len(report["runs"]), "assertions": sum(r["assertions"] for r in report["runs"])}))


if __name__ == "__main__":
    main()
