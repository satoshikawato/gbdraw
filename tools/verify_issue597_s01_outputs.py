"""Check existing Save normalization and Python/CLI replay for S01 measurement files."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

from measure_issue597_s01 import ROOT, document, summary


def compare_svg(expected, actual_path):
    from tests.utils.svg_compare import compare_svgs, parse_svg

    a = parse_svg(expected)
    b = parse_svg(actual_path.read_text())
    if a.attrib.get("baseProfile") == "full" and "baseProfile" not in b.attrib:
        a.attrib.pop("baseProfile")
    result = compare_svgs(ET.tostring(a, encoding="unicode"), str(actual_path))
    return {
        "equal": result.equal,
        "message": result.message,
        "differences": result.differences[:10],
    }


def main():
    sys.path.insert(0, str(ROOT))
    import gbdraw
    from gbdraw.api import (
        load_session_document,
        materialize_session,
        session_to_request,
    )

    assert Path(gbdraw.__file__).resolve().parent == ROOT / "gbdraw"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--saved", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--fresh-load", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source, saved = document(args.source), document(args.saved)
    a, b = summary(args.source), summary(args.saved)
    exact_keys = [
        "renderRequest",
        "resources",
        "proteinIdentityManifest",
        "losatDerivedCache",
    ]
    checks = {k: a["hashes"].get(k) == b["hashes"].get(k) for k in exact_keys}
    checks["sourceEditorStateFields"] = all(
        saved["editorState"].get(key) == value
        for key, value in source["editorState"].items()
    )
    checks["featureCatalog"] = (
        source["editorState"]["featureCatalog"]
        == saved["editorState"]["featureCatalog"]
    )
    entries_a, entries_b = (
        source["losatCache"]["entries"],
        saved["losatCache"]["entries"],
    )
    checks["rawCacheSourceFields"] = len(entries_a) == len(entries_b) and all(
        all(after.get(key) == value for key, value in before.items())
        for before, after in zip(entries_a, entries_b)
    )
    added_cache_fields = sorted(
        {
            key
            for before, after in zip(entries_a, entries_b)
            for key in after.keys() - before.keys()
        }
    )
    allowed_cache_fields = {
        "queryIndex",
        "subjectIndex",
        "queryUid",
        "subjectUid",
        "ordinal",
        "edgeKey",
    }
    checks["cacheOnlyExistingRestoreMetadataAdded"] = (
        set(added_cache_fields) <= allowed_cache_fields
    )
    svg_path = args.output / "saved.svg"
    svg_path.write_text(saved["results"][0]["content"])
    svg = compare_svg(source["results"][0]["content"], svg_path)
    checks["savedSvgSemanticEqual"] = svg["equal"]
    result = {
        "source": a,
        "saved": b,
        "checks": checks,
        "svg": svg,
        "addedCacheFields": added_cache_fields,
        "addedEditorFields": {
            k: v
            for k, v in saved["editorState"].items()
            if k not in source["editorState"]
        },
        "webBindingsSchema": saved.get("webFiles", {})
        .get("bindings", {})
        .get("schema"),
    }
    if args.replay:
        session = load_session_document(args.saved)
        with materialize_session(
            session, output_directory=args.output / "python-materialized"
        ) as materialized:
            request = session_to_request(materialized)
            result["pythonReader"] = {
                "version": session.version,
                "mode": session.mode,
                "records": len(request.records),
            }
        prefix = args.output.resolve() / "cli-replay"
        command = [
            sys.executable,
            "-m",
            "gbdraw.cli",
            "linear",
            "--session",
            str(args.saved.resolve()),
            "-o",
            str(prefix),
            "-f",
            "svg",
            "--overwrite",
        ]
        with (args.output / "cli-replay.log").open("w") as log:
            replay = subprocess.run(
                command,
                cwd=ROOT,
                env={
                    **os.environ,
                    "XDG_CACHE_HOME": str(args.output.resolve() / ".cache"),
                },
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        result["cliReplayExit"] = replay.returncode
        result["cliSvgBytes"] = (
            prefix.with_suffix(".svg").stat().st_size
            if replay.returncode == 0
            else None
        )
        checks["pythonReaderRecords"] = (
            result["pythonReader"]["records"] == a["records"]
        )
        checks["cliReplay"] = replay.returncode == 0 and result["cliSvgBytes"] > 0
        # Browser serialization/sanitization can omit presentation-only SVG metadata.
        # CLI replay is checked as a successful scientific conversion, separately from Save SVG equivalence.
    if args.fresh_load:
        from functools import partial
        from http.server import ThreadingHTTPServer
        import threading

        from playwright.sync_api import sync_playwright
        from measure_issue597_s01 import Handler, configure_page, IMPORT

        server = ThreadingHTTPServer(
            ("127.0.0.1", 0), partial(Handler, directory=str(ROOT))
        )
        threading.Thread(target=server.serve_forever, daemon=True).start()
        try:
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch()
                context, page, errors, dialogs, external = configure_page(
                    browser, f"http://127.0.0.1:{server.server_port}"
                )
                page.locator(IMPORT).set_input_files(str(args.saved.resolve()))
                page.wait_for_function(
                    'window.__s01.events.some(e=>e.name==="interactiveReady") && !window.__GBDRAW_APP__.sessionImportPending'
                )
                ready = page.evaluate("""async () => {
                    const {state}=await import('./js/state.js');
                    const request=(await import('./js/services/config.js')).getCommittedCanonicalRenderRequest();
                    return {request, results:state.results.value.length, cacheEntries:state.losatCache.value.size,
                      features:state.featureCatalog.value.items.reduce((sum,item)=>sum+item.biologicalFeatures.length,0),
                      mounted:Boolean(document.querySelector('.shadow-xl.origin-top > svg')), workers:window.__s01.workers};
                }""")
                checks["freshLoadCanonicalRequest"] = (
                    ready.pop("request") == saved["renderRequest"]
                )
                checks["freshLoadReady"] = ready["mounted"] and ready["results"] == len(
                    saved["results"]
                )
                checks["freshLoadCatalogAndCache"] = (
                    ready["features"] == b["features"]
                    and ready["cacheEntries"] == b["cacheEntries"]
                )
                checks["freshLoadNoWorkerOrExternalError"] = (
                    not ready["workers"] and not errors and not external
                )
                result["freshLoad"] = {
                    **ready,
                    "pageErrors": errors,
                    "externalRequests": external,
                    "dialogs": dialogs,
                }
                context.close()
                browser.close()
        finally:
            server.shutdown()
            server.server_close()
    (args.output / "semantic-checks.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(checks))
    assert all(checks.values()), checks


if __name__ == "__main__":
    main()
