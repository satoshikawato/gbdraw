"""Run S01 probes on an isolated checkout/server/browser; no application writes."""

import argparse
from functools import partial
import gzip
import hashlib
from importlib.metadata import version
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import platform
import socket
import subprocess
import threading
import time
import urllib.request

import psutil
from playwright.sync_api import sync_playwright
import websocket

ROOT = Path(__file__).resolve().parents[1]
OBSERVER = ROOT / "tests/web/fixtures/issue597-session-measurement.js"
VIBRIO = (
    ROOT
    / "gbdraw/web/gallery/sessions/vibrio-harveyi-group-collinear.gbdraw-session.json.gz"
)
IMPORT = 'input[type=file][accept*="application/json"][accept*="application/gzip"]'
TIMEOUT = 1_800_000


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def document(path):
    data = path.read_bytes()
    return json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)


def summary(path):
    data = document(path)

    def normalise_numbers(value):
        # JSON has one numeric type. CLI 1.0 and browser 1 are the same value;
        # booleans retain their distinct JSON type. Byte SHA remains separate.
        if isinstance(value, float) and value.is_integer():
            return int(value)
        if isinstance(value, dict):
            return {key: normalise_numbers(item) for key, item in value.items()}
        if isinstance(value, list):
            return [normalise_numbers(item) for item in value]
        return value

    def digest(v):
        return hashlib.sha256(
            json.dumps(
                normalise_numbers(v),
                sort_keys=True,
                ensure_ascii=True,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()

    return {
        "sha256": sha(path),
        "bytes": path.stat().st_size,
        "version": data.get("version"),
        "requestSchema": data.get("renderRequest", {}).get("schema"),
        "records": len(data.get("renderRequest", {}).get("records", [])),
        "features": sum(
            len(i.get("biologicalFeatures", []))
            for i in (data.get("editorState", {}).get("featureCatalog") or {}).get(
                "items", []
            )
        ),
        "cacheEntries": len(data.get("losatCache", {}).get("entries", [])),
        "hashes": {k: digest(v) for k, v in data.items()},
        "wholeHash": digest(data),
    }


class Handler(SimpleHTTPRequestHandler):
    def log_message(self, *_args):
        pass

    def do_GET(self):
        if self.path == "/__issue597_probe.html":
            body = b'<!doctype html><meta charset="utf-8"><input id="file" type="file"><title>S01 disposable transport</title>'
            self.send_response(200)
            self.send_header("Content-Type", "text/html")
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()


def free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


class MemorySampler:
    """Sample dedicated Chromium descendants and V8 isolates, without touching other browsers."""

    def __init__(self, debug_port):
        self.port = debug_port
        self.stop = threading.Event()
        self.samples = []
        self.errors = []
        self.connections = {}
        self.sequence = 0
        self.pid = None
        self.rss_samples = []
        self.rss_errors = []
        self.rss_thread = threading.Thread(target=self.sample_rss, daemon=True)
        self.thread = threading.Thread(target=self.run, daemon=True)

    def start(self):
        self.rss_thread.start()
        self.thread.start()

    def sample_rss(self):
        while not self.stop.is_set():
            try:
                if self.pid is None:
                    for proc in psutil.process_iter(["pid", "cmdline"]):
                        if f"--remote-debugging-port={self.port}" in (
                            proc.info["cmdline"] or []
                        ):
                            self.pid = proc.pid
                            break
                if self.pid:
                    processes = [
                        psutil.Process(self.pid),
                        *psutil.Process(self.pid).children(recursive=True),
                    ]
                    rss = sum(p.memory_info().rss for p in processes if p.is_running())
                    self.rss_samples.append({"time": time.time(), "rssBytes": rss})
            except psutil.Error as error:
                self.rss_errors.append(type(error).__name__ + ": " + str(error))
            self.stop.wait(0.1)

    def run(self):
        while not self.stop.is_set():
            sample = {"time": time.time(), "rssBytes": 0, "isolates": {}}
            try:
                targets = json.load(
                    urllib.request.urlopen(
                        f"http://127.0.0.1:{self.port}/json/list", timeout=2
                    )
                )
                for target in targets:
                    if target["type"] not in ("page", "worker"):
                        continue
                    key = target["id"]
                    if key not in self.connections:
                        self.connections[key] = websocket.create_connection(
                            target["webSocketDebuggerUrl"],
                            timeout=2,
                            suppress_origin=True,
                        )
                    ws = self.connections[key]
                    self.sequence += 1
                    request_id = self.sequence
                    ws.send(
                        json.dumps({"id": request_id, "method": "Runtime.getHeapUsage"})
                    )
                    while True:
                        reply = json.loads(ws.recv())
                        if reply.get("id") == request_id:
                            break
                    sample["isolates"][target["type"] + ":" + key] = reply.get(
                        "result", {}
                    )
            except Exception as error:
                self.errors.append(type(error).__name__ + ": " + str(error))
            self.samples.append(sample)
            self.stop.wait(0.1)

    def finish(self):
        self.stop.set()
        self.thread.join(5)
        self.rss_thread.join(5)
        for ws in self.connections.values():
            ws.close()
        peaks = {}
        for s in self.samples:
            for key, value in s["isolates"].items():
                kind = key.split(":")[0]
                peaks[kind] = max(peaks.get(kind, 0), value.get("usedSize", 0))
        return {
            "sampleCount": len(self.samples),
            "rssPeakBytes": max(
                (s["rssBytes"] for s in self.rss_samples), default=None
            ),
            "rssErrors": sorted(set(self.rss_errors)),
            "rssStartBytes": next(
                (s["rssBytes"] for s in self.rss_samples if s["rssBytes"]), None
            ),
            "v8UsedHeapPeaksBytes": peaks,
            "errors": sorted(set(self.errors)),
            "samples": self.samples,
            "rssSamples": self.rss_samples,
        }


def compact(raw):
    gaps = sorted(raw.pop("gaps"))
    tasks = raw.pop("longTasks")
    heaps = raw.pop("samples")
    receiver = raw.pop("receiverTasks")
    raw["heartbeat"] = {
        "intervalMs": 100,
        "count": len(gaps),
        "p95Ms": gaps[min(len(gaps) - 1, int(len(gaps) * 0.95))] if gaps else None,
        "maxMs": max(gaps, default=None),
    }
    raw["longTasks"] = {
        "count": len(tasks),
        "totalMs": sum(t["duration"] for t in tasks),
        "maxMs": max((t["duration"] for t in tasks), default=0),
        "entries": tasks,
    }
    raw["mainSampledHeapPeakBytes"] = max(
        [
            raw.get("heapBefore") or 0,
            raw.get("heapAfter") or 0,
            *[s["heap"] or 0 for s in heaps],
        ]
    )
    raw["receiverTaskMaxMs"] = max(receiver, default=0)
    events = raw["events"]
    starts = {}
    stages = {}
    for e in events:
        name = e["name"]
        if name.endswith("-start"):
            starts[name[:-6]] = e["timestamp"]
        if name.endswith("-end") and name[:-4] in starts:
            stages[name[:-4]] = e["timestamp"] - starts[name[:-4]]
    times = {e["name"]: e["timestamp"] for e in events}
    if "current-session-preflight-end" in times and "svg-admission-start" in times:
        stages["restoreBeforeSvgAdmission"] = (
            times["svg-admission-start"] - times["current-session-preflight-end"]
        )
    raw["stagesMs"] = stages
    return raw


def configure_page(browser, url, app=True):
    context = browser.new_context(
        viewport={"width": 1440, "height": 1000}, accept_downloads=True
    )
    page = context.new_page()
    page.set_default_timeout(TIMEOUT)
    errors, dialogs, external = [], [], []
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.on(
        "dialog",
        lambda d: (dialogs.append({"type": d.type, "message": d.message}), d.accept()),
    )

    def network(route):
        if route.request.url.startswith(url) or route.request.url.startswith(
            ("blob:", "data:")
        ):
            route.continue_()
        else:
            external.append(route.request.url)
            route.abort()

    page.route("**/*", network)
    page.add_init_script(OBSERVER.read_text())
    page.goto(url + ("/gbdraw/web/index.html" if app else "/__issue597_probe.html"))
    if app:
        page.wait_for_function(
            "window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length"
        )
    return context, page, errors, dialogs, external


def codec_observers(page):
    page.evaluate("""() => {
      window.__s01Codec = {decodeCpuMs: 0, encodeCpuMs: 0, encodedBytes: 0, decodedBytes: 0,streams:[]};
      for (const [prototype, name, field] of [[TextDecoder.prototype,'decode','decodeCpuMs'],
                                            [TextEncoder.prototype,'encodeInto','encodeCpuMs'],
                                            [TextEncoder.prototype,'encode','encodeCpuMs']]) {
        const original = prototype[name];
        prototype[name] = function(...args) {
          const start = performance.now();
          const result = original.apply(this, args);
          window.__s01Codec[field] += performance.now()-start;
          if (name==='decode') window.__s01Codec.decodedBytes += args[0]?.byteLength || 0;
          else window.__s01Codec.encodedBytes += name==='encode' ? result.byteLength : result.written;
          return result;
        };
      }
      const original = DOMPurify.sanitize;
      DOMPurify.sanitize = function(...args) {
        const start = performance.now(); const result = original.apply(this,args);
        window.__s01.events.push({name:'sanitize-call',timestamp:start,durationMs:performance.now()-start});
        return result;
      };
    }""")


def pipeline(browser, url, port, fixture, out, label):
    context, page, errors, dialogs, external = configure_page(browser, url)
    codec_observers(page)
    sampler = MemorySampler(port)
    sampler.start()
    page.evaluate("window.__s01.reset()")
    page.locator(IMPORT).set_input_files(str(fixture))
    page.wait_for_function(
        'window.__s01.events.some(e=>e.name==="interactiveReady") && window.__GBDRAW_APP__.sessionImportPending===false'
    )
    page.wait_for_timeout(200)
    load = compact(page.evaluate("window.__s01.finish()"))
    load["codec"] = page.evaluate("window.__s01Codec")
    load["memory"] = sampler.finish()
    ready = page.evaluate("""async () => {
      const {state}=await import('./js/state.js');
      return {resultCount: state.results.value.length, records: state.linearSeqs.length,
        catalog: state.featureCatalog?.value?.items?.length,
        requestRecords: (await import('./js/services/config.js')).getCommittedCanonicalRenderRequest()?.records?.length,
        cacheEntries: state.losatCache.value.size,
        mounted: Boolean(document.querySelector('.shadow-xl.origin-top > svg'))};
    }""")
    load["ready"] = ready
    if not ready["mounted"]:
        raise AssertionError(f"Load failed: {dialogs}")
    page.evaluate("""() => {
      window.__GBDRAW_APP__.sessionTitle='S01 measurement';
      window.__s01Codec={decodeCpuMs:0,encodeCpuMs:0,encodedBytes:0,decodedBytes:0,streams:[]};
      window.__s01.reset();
    }""")
    sampler = MemorySampler(port)
    sampler.start()
    saved = out / f"{label}-saved.gbdraw-session.json.gz"
    with page.expect_download(timeout=TIMEOUT) as download:
        outcome = page.evaluate("""async () => {
          const outcome=await window.__GBDRAW_APP__.saveSessionWithTitle();
          return {status:outcome.status, compressedBytes:outcome.blob?.size};
        }""")
    download.value.save_as(saved)
    page.wait_for_timeout(200)
    save = compact(page.evaluate("window.__s01.finish()"))
    save["codec"] = page.evaluate("window.__s01Codec")
    save["memory"] = sampler.finish()
    save["outcome"] = outcome
    context.close()
    source_summary, saved_summary = summary(fixture), summary(saved)
    keys = [
        "renderRequest",
        "resources",
        "webFiles",
        "editorState",
        "losatCache",
        "losatDerivedCache",
        "proteinIdentityManifest",
    ]
    equal = {
        k: source_summary["hashes"].get(k) == saved_summary["hashes"].get(k)
        for k in keys
    }
    return {
        "load": load,
        "save": save,
        "source": source_summary,
        "saved": saved_summary,
        "sectionEquality": equal,
        "pageErrors": errors,
        "dialogs": dialogs,
        "externalRequests": external,
    }


def transport(browser, url, port, fixture, out, label, mode, repetition):
    context, page, errors, dialogs, external = configure_page(browser, url, app=False)
    page.locator("#file").set_input_files(str(fixture))
    sampler = MemorySampler(port)
    sampler.start()
    raw = page.evaluate(
        """async mode => await window.__s01Transport(document.querySelector('#file').files[0], mode)""",
        mode,
    )
    raw["memory"] = sampler.finish()
    result = compact(raw)
    # Hash/equivalence serialization is outside responsiveness interval.
    with page.expect_download(timeout=TIMEOUT) as download:
        page.evaluate("""async () => {
          const {compressSessionData}=await import('/gbdraw/web/js/services/session-file.js');
          const blob=await compressSessionData(window.__s01Candidate);
          const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='probe.json.gz';a.click();
          setTimeout(()=>URL.revokeObjectURL(a.href),0);
        }""")
    reply = out / f"{label}-{mode}-{repetition}.json.gz"
    download.value.save_as(reply)
    result["wholeDocumentEqual"] = (
        summary(fixture)["wholeHash"] == summary(reply)["wholeHash"]
    )
    result["pageErrors"], result["externalRequests"] = errors, external
    context.close()
    assert result["wholeDocumentEqual"], "Transport payload changed"
    assert not errors and not external, "Transport page/network failure"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fixture", type=Path)
    parser.add_argument("--mode", choices=["pipeline", "transport"], required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer(
        ("127.0.0.1", 0), partial(Handler, directory=str(ROOT))
    )
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_port}"
    port = free_port()
    evidence = {
        "head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpuCount": os.cpu_count(),
        "cpuModel": next(
            (
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            platform.machine(),
        ),
        "versions": {
            name: version(name)
            for name in ("playwright", "psutil", "websocket-client", "biopython")
        },
        "workerProbeSha256": sha(
            ROOT / "tests/web/fixtures/issue597-session-transport.worker.js"
        ),
        "ramBytes": psutil.virtual_memory().total,
        "url": url,
        "debugPort": port,
        "launchArgs": [
            "--enable-precise-memory-info",
            f"--remote-debugging-port={port}",
        ],
        "observerSha256": sha(OBSERVER),
        "fixtures": {},
    }
    fixtures = {"vibrio": VIBRIO}
    if args.fixture:
        fixtures["real"] = args.fixture.resolve()
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(args=evidence["launchArgs"])
            evidence["browser"] = browser.version
            for label, fixture in fixtures.items():
                result = []
                for i in range(args.repetitions):
                    print(f"{label} {args.mode} {i + 1}/{args.repetitions}", flush=True)
                    if args.mode == "pipeline":
                        result.append(
                            pipeline(browser, url, port, fixture, out, label + f"-{i}")
                        )
                    else:
                        result.append(
                            {
                                mode: transport(
                                    browser, url, port, fixture, out, label, mode, i
                                )
                                for mode in ["whole", "bounded"]
                            }
                        )
                    evidence["fixtures"][label] = result
                    (out / f"{args.mode}-metrics.json").write_text(
                        json.dumps(evidence, indent=2) + "\n"
                    )
            browser.close()
    except Exception as error:
        evidence["failure"] = {"type": type(error).__name__, "message": str(error)}
        (out / f"{args.mode}-metrics.json").write_text(
            json.dumps(evidence, indent=2) + "\n"
        )
        raise
    finally:
        server.shutdown()
        server.server_close()


if __name__ == "__main__":
    main()
