// Observers for the unchanged application and a disposable transport page.
window.__s01 = { events: [], metrics: {}, workers: [], gaps: [], longTasks: [], samples: [],
  receiverTasks: [], reset() {
    this.events = []; this.metrics = {}; this.gaps = []; this.longTasks = [];
    this.samples = []; this.receiverTasks = []; this.start = performance.now();
    this.last = this.start; this.heapBefore=performance.memory?.usedJSHeapSize ?? null; this.active = true;
  }, finish() { this.active = false; return {
    wallMs: performance.now() - this.start, heapBefore: this.heapBefore, heapAfter:performance.memory?.usedJSHeapSize ?? null, events: this.events, metrics: this.metrics,
    gaps: this.gaps, longTasks: this.longTasks, samples: this.samples,
    workers: this.workers, receiverTasks: this.receiverTasks
  }; }
};
const probe = window.__s01;
setInterval(() => {
  if (!probe.active) return;
  const now = performance.now();
  probe.gaps.push(now - probe.last); probe.last = now;
  probe.samples.push({ timestamp: now, heap: performance.memory?.usedJSHeapSize ?? null });
}, 100);
new PerformanceObserver(list => {
  for (const entry of list.getEntries()) {
    if (entry.startTime >= probe.start) probe.longTasks.push({ start: entry.startTime, duration: entry.duration });
  }
}).observe({ type: 'longtask', buffered: true });
window.__GBDRAW_TEST_HOOKS__ = {
  onSessionLifecycleEvent: e => probe.events.push(e),
  onStructuralMetric: e => probe.metrics[e.name] = (probe.metrics[e.name] || 0) + e.value
};
const NativeWorker = window.Worker;
window.Worker = class extends NativeWorker {
  constructor(...args) { super(...args); probe.workers.push(String(args[0])); }
};
window.__s01Transport = async (file, mode) => {
  let candidate;
  let stringChunks;
  let output;
  const worker = new Worker('/tests/web/fixtures/issue597-session-transport.worker.js', {type: 'module'});
  const locate = path => path.reduce((value, key) => value[key], candidate);
  const set = (path, value) => {
    if (!path.length) candidate = value;
    else Object.defineProperty(locate(path.slice(0, -1)), path.at(-1),
      { value, configurable: true, enumerable: true, writable: true });
  };
  probe.reset();
  try {
    output = await new Promise((resolve, reject) => {
      worker.onerror = e => reject(new Error(e.message));
      worker.onmessage = ({data}) => {
        if (data.id !== 1) return;
        if (data.kind === 'error') { reject(new Error(data.message)); return; }
        if (data.kind === 'parsed') { probe.events.push({name: 'worker-parsed', timestamp: performance.now()}); return; }
        if (data.kind === 'done') { resolve(data); return; }
        const start = performance.now();
        try {
          if (data.kind === 'value') set(data.path, data.value);
          if (data.kind === 'batch') {
            const array = locate(data.path);
            data.value.forEach((value, i) => array[data.index + i] = value);
          }
          if (data.kind === 'string-start') { stringChunks = []; }
          if (data.kind === 'string-chunk') {
            for (let index=0; index<data.units.length; index+=8192)
              stringChunks.push(String.fromCharCode(...data.units.subarray(index,index+8192)));
          }
          if (data.kind === 'string-end') { set(data.path, stringChunks.join('')); stringChunks = null; }
          probe.receiverTasks.push(performance.now() - start);
          worker.postMessage({ack: true});
        } catch (error) { reject(error); }
      };
      worker.postMessage({id: 1, file, mode});
    });
    // Include the timer/longtask observation for the final message.
    await new Promise(resolve => setTimeout(resolve, 120));
    return { ...probe.finish(), transport: output };
  } finally {
    worker.terminate();
    // The untrusted plain candidate is not routed into the application.
    window.__s01Candidate = candidate;
  }
};

// Stream spans overlap decode/compression; they are not additive exclusive CPU stages.
window.__s01Codec = {decodeCpuMs:0,encodeCpuMs:0,encodedBytes:0,decodedBytes:0,streams:[]};
const observeReadable = (stream, kind, start) => {
  const metric={kind,start,bytes:0,chunks:0};
  window.__s01Codec.streams.push(metric);
  return stream.pipeThrough(new TransformStream({
    transform(chunk,controller) {metric.bytes+=chunk.byteLength;metric.chunks++;controller.enqueue(chunk);},
    flush() {metric.end=performance.now();metric.wallMs=metric.end-start;}
  }));
};
const nativeBlobStream=Blob.prototype.stream;
Blob.prototype.stream=function(...args) {
  const stream=nativeBlobStream.apply(this,args);
  return probe.active ? observeReadable(stream,'file-read',performance.now()) : stream;
};
const NativeDecompress=window.DecompressionStream;
window.DecompressionStream=class extends NativeDecompress {
  constructor(...args) {
    super(...args);
    if (probe.active) return {writable:this.writable,
      readable:observeReadable(this.readable,'decompress',performance.now())};
  }
};
