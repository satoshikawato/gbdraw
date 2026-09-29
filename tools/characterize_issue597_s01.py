"""Record missing S01 discovery baselines; expectations describe the current runtime."""

import argparse
from functools import partial
from http.server import ThreadingHTTPServer
import json
import gzip
import tempfile
from pathlib import Path
import platform
import subprocess
import threading

from playwright.sync_api import sync_playwright
from measure_issue597_s01 import (
    ROOT,
    Handler,
    IMPORT,
    configure_page,
    sha,
    free_port,
    transport,
)

SNAPSHOT = """async () => {
  const {state}=await import('./js/state.js');
  const {getCommittedCanonicalRenderRequest}=await import('./js/services/config.js');
  return {status:state.circularRecordDiscovery.status,error:state.circularRecordDiscovery.error,
    records:state.circularRecordList.value,workers:window.__s01.workers,
    inputType:state.cInputType.value,primaryName:state.circularRecordDiscovery.primaryFile?.name,
    resultCount:state.results.value.length,committedRecords:getCommittedCanonicalRenderRequest()?.records?.length ?? null,
    controls:document.querySelectorAll('[aria-label^="Display start"]').length,
    panelOpen:document.querySelector('[data-circular-record-presentation]')?.open};
}"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    server = ThreadingHTTPServer(
        ("127.0.0.1", 0), partial(Handler, directory=str(ROOT))
    )
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_port}"
    evidence = {
        "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "platform": platform.platform(),
        "observations": {},
        "fixtures": {},
    }
    inputs = [
        "tests/fixtures/sessions/cli-web-mito.gb",
        "tests/test_inputs/MjeNMV.gbk",
        "tests/test_inputs/NC_013668.gff3",
        "tests/test_inputs/NC_013668.fasta",
        "gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json",
        "tests/fixtures/sessions/settings-only.v42.json.gz",
    ]
    evidence["fixtures"] = {p: sha(ROOT / p) for p in inputs}

    def observe(name, page):
        value = page.evaluate(SNAPSHOT)
        evidence["observations"][name] = value
        args.output.write_text(json.dumps(evidence, indent=2) + "\n")
        print(name, value["status"], len(value["records"]), flush=True)
        return value

    try:
        with sync_playwright() as p:
            debug_port = free_port()
            browser = p.chromium.launch(args=[f"--remote-debugging-port={debug_port}"])
            evidence["browser"] = browser.version
            context, page, errors, dialogs, external = configure_page(browser, url)
            upload = page.get_by_label("GenBank/DDBJ File", exact=True)
            source = (ROOT / inputs[0]).read_bytes()
            upload.set_input_files(
                {
                    "name": "duplicate.gb",
                    "mimeType": "text/plain",
                    "buffer": source + b"\n" + source,
                }
            )
            page.wait_for_function(
                "window.__GBDRAW_APP__.circularRecordList.length===2"
            )
            observe("duplicateNative", page)
            assert [
                r["selector"]
                for r in evidence["observations"]["duplicateNative"]["records"]
            ] == ["#1", "#2"]
            assert evidence["observations"]["duplicateNative"]["workers"] == []
            upload.set_input_files(str(ROOT / inputs[1]))
            page.wait_for_function(
                'window.__GBDRAW_APP__.circularRecordList[0]?.record_id==="LC738868.1"'
            )
            observe("ddbjAccessionNative", page)
            assert evidence["observations"]["ddbjAccessionNative"]["workers"] == []
            # Delay only this native File's first text read. Exercise the actual watcher/version guard.
            page.evaluate(
                """async text => {
              const {state}=await import('./js/state.js');
              const slow=new File([text],'slow-old.gb');
              const pending=new Promise(resolve=>window.__s01ReleaseOld=()=>resolve(new TextEncoder().encode(text).buffer));
              slow.arrayBuffer=()=>pending;
              state.files.c_gb=slow;
            }""",
                source.decode(),
            )
            page.wait_for_function("Boolean(window.__s01ReleaseOld)")
            page.evaluate(
                """async text => {
              const {state}=await import('./js/state.js');
              state.files.c_gb=new File([text],'fast-new.gb');
            }""",
                (ROOT / inputs[1]).read_text(),
            )
            page.wait_for_function(
                'window.__GBDRAW_APP__.circularRecordList[0]?.record_id==="LC738868.1"'
            )
            page.evaluate("window.__s01ReleaseOld()")
            page.wait_for_timeout(200)
            observe("rapidReplaceStaleSettlement", page)
            assert (
                evidence["observations"]["rapidReplaceStaleSettlement"]["records"][0][
                    "record_id"
                ]
                == "LC738868.1"
            )
            # The delayed read is a watcher-only probe. Isolate ordinary UI/History checks.
            context.close()
            context, page, more_errors, more_dialogs, more_external = configure_page(
                browser, url
            )
            first_errors, first_dialogs, first_external = errors, dialogs, external
            errors, dialogs, external = more_errors, more_dialogs, more_external
            upload = page.get_by_label("GenBank/DDBJ File", exact=True)
            upload.set_input_files(str(ROOT / inputs[1]))
            page.wait_for_function(
                "window.__GBDRAW_APP__.circularRecordList.length===1"
            )
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            ).locator("button").click()
            page.wait_for_function(
                "window.__GBDRAW_APP__.circularRecordList.length===0"
            )
            observe("remove", page)
            page.locator("input[type=radio][value=gff]").first.check()
            page.get_by_label("GFF3 File", exact=True).set_input_files(
                str(ROOT / inputs[2])
            )
            page.wait_for_timeout(200)
            observe("incompleteGffPair", page)
            assert evidence["observations"]["incompleteGffPair"]["records"] == []
            page.get_by_label("FASTA File", exact=True).set_input_files(
                str(ROOT / inputs[3])
            )
            page.wait_for_function(
                "window.__GBDRAW_APP__.circularRecordList.length===1"
            )
            observe("completeGffPair", page)
            assert evidence["observations"]["completeGffPair"]["workers"] == []
            page.locator(IMPORT).set_input_files(str(ROOT / inputs[4]))
            page.wait_for_function(
                "window.__GBDRAW_APP__.sessionImportPending===false && window.__GBDRAW_APP__.results.length===1"
            )
            observe("savedPreview", page)
            assert evidence["observations"]["savedPreview"]["workers"] == []
            page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
                str(ROOT / inputs[1])
            )
            page.wait_for_function(
                'window.__GBDRAW_APP__.circularRecordList[0]?.record_id==="LC738868.1"'
            )
            observe("draftDiffersFromArtifact", page)
            assert (
                evidence["observations"]["draftDiffersFromArtifact"]["resultCount"] == 1
            )
            # Force the existing parser/helper boundary by rejecting only this source's text read.
            helper = page.evaluate("""async () => {
              const {discoverSequenceRecords}=await import('./js/app/record-discovery.js');
              const source=window.__GBDRAW_APP__.files.c_gb;
              return await discoverSequenceRecords({file:source,format:'genbank',
                readText:async()=>{throw new Error('S01 forced fast-path miss');}});
            }""")
            evidence["observations"]["existingHelperSettlement"] = {
                "records": helper,
                "workers": page.evaluate("window.__s01.workers"),
            }
            assert helper[0]["recordId"] == "LC738868.1"
            page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
                {
                    "name": "invalid.gb",
                    "mimeType": "text/plain",
                    "buffer": b"not a biological record",
                }
            )
            page.wait_for_function(
                "async () => (await import('./js/state.js')).state.circularRecordDiscovery.status==='error'"
            )
            observe("invalidSource", page)
            assert evidence["observations"]["invalidSource"]["resultCount"] == 1
            context.close()
            context, page, settings_errors, settings_dialogs, settings_external = (
                configure_page(browser, url)
            )
            page.locator(IMPORT).set_input_files(str(ROOT / inputs[5]))
            page.wait_for_function(
                'window.__s01.events.some(e=>e.name==="interactiveReady") && !window.__GBDRAW_APP__.sessionImportPending'
            )
            observe("settingsOnly", page)
            assert evidence["observations"]["settingsOnly"]["resultCount"] == 0
            assert evidence["observations"]["settingsOnly"]["workers"] == []
            evidence.update(
                {
                    "pageErrors": first_errors + errors + settings_errors,
                    "externalRequests": first_external + external + settings_external,
                    "dialogs": first_dialogs + dialogs + settings_dialogs,
                }
            )
            assert not evidence["pageErrors"] and not evidence["externalRequests"]
            context.close()
            # Auxiliary codec stress, not a replacement for the real biological fixture.
            with tempfile.TemporaryDirectory(
                prefix="issue597-s01-unicode-"
            ) as temporary:
                out = Path(temporary)
                text = json.dumps(
                    {
                        "resources": {
                            "huge": "x" * 131071 + "😀\ud800\udfff" + "β" * 200000
                        },
                        "array": [None, {}, [], [1, False, "𝄞"]],
                        "__proto__": {"value": "retained"},
                    },
                    ensure_ascii=True,
                ).encode()
                evidence["unicodeTransport"] = {}
                for name, payload in [
                    ("json", text),
                    ("gzip", gzip.compress(text, mtime=0)),
                ]:
                    fixture = out / (name + ".data")
                    fixture.write_bytes(payload)
                    for mode in ["whole", "bounded"]:
                        observed = transport(
                            browser, url, debug_port, fixture, out, name, mode, 0
                        )
                        evidence["unicodeTransport"][name + "-" + mode] = {
                            "wholeDocumentEqual": observed["wholeDocumentEqual"],
                            "maximumMessageUpperBytes": observed["transport"][
                                "maximumMessageUpperBytes"
                            ],
                            "transferredBytes": observed["transport"][
                                "transferredBytes"
                            ],
                        }
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
        args.output.write_text(json.dumps(evidence, indent=2) + "\n")
    print("S01 discovery characterization assertions passed")


if __name__ == "__main__":
    main()
