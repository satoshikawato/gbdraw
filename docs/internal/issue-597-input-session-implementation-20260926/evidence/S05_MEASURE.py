"""Measure the production S05 transport and pipeline with S01 metric definitions."""

import argparse
from functools import partial
import hashlib
from http.server import ThreadingHTTPServer
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import threading
import time
import urllib.request

import psutil
from playwright.sync_api import sync_playwright
import websocket

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools"))
import measure_issue597_s01 as s01  # noqa: E402

# These test-served wrappers count codec work, never log uploaded content or pause execution.
CODEC_PROBE = r"""
const probe={decodeCpuMs:0,decodedBytes:0,streams:[]};
const originalDecode=TextDecoder.prototype.decode;
TextDecoder.prototype.decode=function(...args){
  const start=performance.now();const result=originalDecode.apply(this,args);
  probe.decodeCpuMs+=performance.now()-start;probe.decodedBytes+=args[0]?.byteLength||0;
  return result;
};
const observe=(stream,kind)=>{
  const metric={kind,start:performance.now(),bytes:0,chunks:0};probe.streams.push(metric);
  return stream.pipeThrough(new TransformStream({
    transform(chunk,controller){metric.bytes+=chunk.byteLength;metric.chunks++;controller.enqueue(chunk);},
    flush(){metric.wallMs=performance.now()-metric.start;}
  }));
};
const originalStream=Blob.prototype.stream;
Blob.prototype.stream=function(...args){return observe(originalStream.apply(this,args),'file-read');};
const NativeDecompressionStream=DecompressionStream;
self.DecompressionStream=class extends NativeDecompressionStream{
  constructor(...args){super(...args);return {writable:this.writable,readable:observe(this.readable,'decompress')};}
};
const originalPost=self.postMessage.bind(self);
self.postMessage=reply=>{if(reply.timings)reply.timings.codecProbe=probe;return originalPost(reply);};
"""
TRACE = r"""() => {
  const NativeWorker=Worker;
  window.__s05Workers=[];
  window.Worker=class extends NativeWorker {
    constructor(...args){super(...args);this.entry={url:String(args[0]),stack:new Error().stack,
      messages:[],listeners:0,terminations:0};window.__s05Workers.push(this.entry);}
    postMessage(message,...args){this.entry.messages.push({type:message?.type||null,
      operation:message?.operation||null,operationId:message?.operationId||null,
      fileBytes:message?.file?.size??null,stack:new Error().stack});return super.postMessage(message,...args);}
    addEventListener(...args){this.entry.listeners++;return super.addEventListener(...args);}
    removeEventListener(...args){this.entry.listeners--;return super.removeEventListener(...args);}
    terminate(){this.entry.terminations++;return super.terminate();}
  };
  const finish=window.__s01.finish.bind(window.__s01);
  window.__s01.finish=()=>({...finish(),workerTrace:window.__s05Workers});
}"""

original_configure = s01.configure_page


def configure(browser, url, app=True):
    context, page, errors, dialogs, external = original_configure(browser, url, app)
    page.route("**/workers/session-import-worker.js", lambda route: route.fulfill(
        content_type="text/javascript",
        body=CODEC_PROBE + (ROOT / "gbdraw/web/js/workers/session-import-worker.js").read_text(),
    ))
    page.evaluate(TRACE)
    return context, page, errors, dialogs, external


class Sampler(s01.MemorySampler):
    def __init__(self, port):
        super().__init__(port)
        self.targets = {}

    def run(self):
        while not self.stop.is_set():
            sample = {"time": time.time(), "isolates": {}}
            try:
                targets = json.load(urllib.request.urlopen(
                    f"http://127.0.0.1:{self.port}/json/list", timeout=2))
                for target in targets:
                    if target["type"] not in ("page", "worker"):
                        continue
                    key = target["id"]
                    self.targets[key] = {"type": target["type"], "url": target["url"]}
                    if key not in self.connections:
                        self.connections[key] = websocket.create_connection(
                            target["webSocketDebuggerUrl"], timeout=2, suppress_origin=True)
                    ws = self.connections[key]
                    self.sequence += 1
                    ws.send(json.dumps({"id": self.sequence, "method": "Runtime.getHeapUsage"}))
                    while True:
                        reply = json.loads(ws.recv())
                        if reply.get("id") == self.sequence:
                            break
                    sample["isolates"][target["type"] + ":" + key] = reply.get("result", {})
            except Exception as error:
                self.errors.append(type(error).__name__ + ": " + str(error))
            self.samples.append(sample)
            self.stop.wait(0.1)

    def finish(self):
        result = super().finish()
        result["targets"] = self.targets
        return result


def transport(browser, url, port, fixture, output, label):
    context, page, errors, dialogs, external = configure(browser, url, app=False)
    page.locator("#file").set_input_files(str(fixture))
    sampler = Sampler(port)
    sampler.start()
    raw = page.evaluate("""async () => {
      const {importSessionFile}=await import('/gbdraw/web/js/services/session-import-client.js');
      window.__s01.reset();
      const reply=await importSessionFile(document.querySelector('#file').files[0]);
      window.__s05Candidate=reply.data;
      const completionMs=performance.now()-window.__s01.start;
      // Same 200 ms terminal observation as S01 pipeline; excluded from completionMs.
      await new Promise(resolve=>setTimeout(resolve,200));
      return {...window.__s01.finish(),completionMs,workerCodec:reply.timings,characters:reply.characters};
    }""")
    raw["memory"] = sampler.finish()
    result = s01.compact(raw)
    result["receiverTaskMaxMs"] = None  # native reply deserialization is not directly observed
    assert "codecProbe" in result["workerCodec"]
    result["fileBytes"] = fixture.stat().st_size
    result["transferredBufferBytes"] = 0  # the production postMessage has no transfer list
    result["nativeStructuredCloneWireBytes"] = None
    result["nativeStructuredCloneCopyBytes"] = None
    # Complete graph comparison/serialization follows the measured responsiveness interval.
    with page.expect_download(timeout=s01.TIMEOUT) as download:
        page.evaluate("""async () => {
          const {compressSessionData}=await import('/gbdraw/web/js/services/session-file.js');
          const blob=await compressSessionData(window.__s05Candidate);
          const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='reply.json.gz';a.click();
          setTimeout(()=>URL.revokeObjectURL(a.href),0);
        }""")
    saved = output / f"{label}-reply.json.gz"
    download.value.save_as(saved)
    context.close()
    result["wholeDocumentEqual"] = s01.summary(fixture)["wholeHash"] == s01.summary(saved)["wholeHash"]
    result["pageErrors"], result["externalRequests"] = errors, external
    assert result["wholeDocumentEqual"] and not errors and not external
    return result


def limits(browser, url, fixture, output, label):
    import gzip
    import shutil
    plain = output / 'real-expanded-rejection.json'
    with gzip.open(fixture, 'rb') as source, plain.open('wb') as target:
        shutil.copyfileobj(source, target)
    assert plain.stat().st_size > 200 * 1024 * 1024
    bomb = output / 'expanded-limit-rejection.json.gz'
    with gzip.open(bomb, 'wb') as target:
        block = b' ' * (1024 * 1024)
        for _ in range(513):
            target.write(block)
    context, page, errors, dialogs, external = configure(browser, url, app=False)
    rows = []
    for path, expected in [(plain, 'Session file is too large.'),
                           (bomb, 'Expanded session file is too large.')]:
        page.locator('#file').set_input_files(str(path))
        reply = page.evaluate("""async () => {
          const {importSessionFile}=await import('/gbdraw/web/js/services/session-import-client.js');
          try { await importSessionFile(document.querySelector('#file').files[0]); return {status:'ok'}; }
          catch(error){return {status:'error',message:error.message,stage:error.stage,code:error.code};}
        }""")
        assert reply['status'] == 'error' and reply['message'] == expected
        rows.append({'path': str(path), 'fileBytes': path.stat().st_size, 'reply': reply})
    trace = page.evaluate('window.__s05Workers')
    context.close()
    assert all(w['terminations'] == 1 and w['listeners'] == 0 for w in trace)
    return {'rejections': rows, 'workerTrace': trace, 'pageErrors': errors, 'externalRequests': external}


def journey(browser, url, fixture, output, label):
    context, page, errors, dialogs, external = configure(browser, url)
    page.locator(s01.IMPORT).set_input_files(str(fixture))
    page.wait_for_function('window.__s01.events.some(e=>e.name==="interactiveReady") && !window.__GBDRAW_APP__.sessionImportPending')
    before = page.evaluate("""async () => {
      const service=await import('./js/services/config.js');
      return {request:service.getCommittedCanonicalRenderRequest(),workers:window.__s05Workers,
        preview:Boolean(document.querySelector('.shadow-xl.origin-top > svg'))};
    }""")
    outcome = page.evaluate("""async () => {
      const app=window.__GBDRAW_APP__;
      const result=await app.runAnalysis();
      return {status:result?.status,error:String(app.errorLog||''),results:app.results.length,
        preview:Boolean(document.querySelector('.shadow-xl.origin-top > svg')),workers:window.__s05Workers};
    }""")
    saved = output / f"{label}-generated.json.gz"
    page.evaluate("window.__GBDRAW_APP__.sessionTitle='S05 generated journey'")
    with page.expect_download(timeout=s01.TIMEOUT) as download:
        status = page.evaluate("async () => (await window.__GBDRAW_APP__.saveSessionWithTitle()).status")
    download.value.save_as(saved)
    context.close()
    result = {"beforeGenerate": before, "generate": outcome, "saveStatus": status,
              "saved": str(saved), "pageErrors": errors, "externalRequests": external}
    result["pythonPreviewWorkerCount"] = sum("diagram-generation-worker.js" in w["url"] for w in before["workers"])
    result["previewWorkerFreeAcceptance"] = result["pythonPreviewWorkerCount"] == 0
    result["generateAcceptance"] = outcome["status"] == "ok" and outcome["preview"] and not outcome["error"]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=("transport", "pipeline", "journey", "limits"), required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(s01.Handler, directory=str(ROOT)))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_port}"
    port = s01.free_port()
    s01.configure_page = configure
    s01.MemorySampler = Sampler
    result = {"mode": args.mode, "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "platform": platform.platform(), "python": platform.python_version(),
              "cpuCount": os.cpu_count(), "ramBytes": psutil.virtual_memory().total,
              "cpuModel": next((line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name")), platform.machine()),
              "versions": {name: version(name) for name in ("playwright", "psutil", "websocket-client", "biopython")},
              "fixtures": {}, "codecProbeSha256": hashlib.sha256(CODEC_PROBE.encode()).hexdigest(),
              "traceSha256": hashlib.sha256(TRACE.encode()).hexdigest(),
              "observerSha256": s01.sha(s01.OBSERVER),
              "recipeSha256": s01.sha(Path(__file__)), "s01RecipeSha256": s01.sha(Path(s01.__file__)),
              "url": url, "debugPort": port,
              "launchArgs": ["--enable-precise-memory-info", f"--remote-debugging-port={port}"],
              "startedUtc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    fixtures = {"vibrio": s01.VIBRIO, "real": args.fixture.resolve()}
    if args.mode == "limits":
        fixtures = {"real": args.fixture.resolve()}
    if args.mode == "transport":
        fixtures["small-plain"] = ROOT / "gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json"
    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(args=["--enable-precise-memory-info", f"--remote-debugging-port={port}"])
            result["browser"] = browser.version
            for label, fixture in fixtures.items():
                rows = []
                for i in range(args.repetitions):
                    print(f"{label} {args.mode} {i + 1}/{args.repetitions}", flush=True)
                    name = f"{label}-{i}"
                    if args.mode == "transport":
                        row = transport(browser, url, port, fixture, out, name)
                    elif args.mode == "pipeline":
                        row = s01.pipeline(browser, url, port, fixture, out, name)
                        row["load"]["receiverTaskMaxMs"] = None
                        row["load"]["transferredBufferBytes"] = 0
                        row["load"]["nativeStructuredCloneWireBytes"] = None
                        row["load"]["nativeStructuredCloneCopyBytes"] = None
                    elif args.mode == "limits":
                        row = limits(browser, url, fixture, out, name)
                    else:
                        row = journey(browser, url, fixture, out, name)
                    rows.append(row)
                    result["fixtures"][label] = {"sha256": s01.sha(fixture), "results": rows}
                    (out / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
            browser.close()
    except Exception as error:
        result["failure"] = {"type": type(error).__name__, "message": str(error)}
        (out / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
        raise
    finally:
        server.shutdown()
        server.server_close()


if __name__ == "__main__":
    main()
