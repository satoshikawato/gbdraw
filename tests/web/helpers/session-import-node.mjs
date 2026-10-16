import { Worker as NodeWorker } from 'node:worker_threads';

// Browser event API around a real Node isolate running the production handler.
// The production client starts a Worker per import and terminates it when the
// import ends. Starting a Node isolate costs about 50 ms, so the adapter keeps
// one started isolate per test file: `terminate()` after the import finished
// parks it for the next import; `terminate()` while an import still runs (a
// cancel) stops it, as in the browser.
/** @type {{ url: string, worker: NodeWorker, ready: Promise<void> } | null} */
let parked = null;

const startIsolate = (url) => {
  const worker = new NodeWorker(`
      const {parentPort}=require('node:worker_threads');
      global.self={addEventListener:(_,fn)=>parentPort.on('message',data=>fn({data})),
        postMessage:data=>parentPort.postMessage(data)};
      import(${JSON.stringify(url)}).then(()=>parentPort.postMessage({ready:true}));
    `, { eval: true });
  const ready = new Promise((resolve) => {
    const started = (data) => {
      if (!data?.ready) return;
      worker.off('message', started);
      resolve();
    };
    worker.on('message', started);
  });
  // A parked isolate has no import listening; one that fails or exits is not reused.
  const unpark = () => { if (parked?.worker === worker) parked = null; };
  worker.on('error', unpark);
  worker.once('exit', unpark);
  return { url, worker, ready };
};

export class NodeSessionImportWorker {
  constructor(url) {
    if (!String(url).endsWith('/session-import-worker.js')) {
      throw new Error('The Node import test adapter does not supply a diagram Worker.');
    }
    const isolate = parked?.url === String(url) ? parked : startIsolate(String(url));
    if (isolate === parked) parked = null;
    isolate.worker.ref();
    this.worker = isolate.worker;
    this.isolate = isolate;
    this.ready = isolate.ready;
    this.listeners = new Map();
    // The operation sent and not yet answered with `ok` or `error`.
    this.running = null;
    this.forward = {
      message: (data) => {
        if (data?.ready) return;
        if (data?.operationId === this.running && (data.status === 'ok' || data.status === 'error')) this.running = null;
        this.emit('message', { data });
      },
      error: (error) => this.emit('error', { message: error.message }),
      messageerror: () => this.emit('messageerror', {})
    };
    for (const [name, fn] of Object.entries(this.forward)) this.worker.on(name, fn);
  }
  emit(name, event) { for (const fn of this.listeners.get(name) || []) fn(event); }
  addEventListener(name, fn) {
    if (!this.listeners.has(name)) this.listeners.set(name, new Set());
    this.listeners.get(name).add(fn);
  }
  removeEventListener(name, fn) { this.listeners.get(name)?.delete(fn); }
  postMessage(message) {
    if (!message.ack) this.running = message.operationId;
    this.ready.then(() => this.worker.postMessage(message));
  }
  terminate() {
    this.listeners.clear();
    for (const [name, fn] of Object.entries(this.forward)) this.worker.off(name, fn);
    if (this.running === null && this.worker.threadId !== -1 && !parked) {
      this.worker.unref();
      parked = this.isolate;
    } else {
      this.worker.terminate();
    }
  }
}

export const installSessionImportWorker = () => { globalThis.Worker = NodeSessionImportWorker; };

// A Session older than 40 reads its features again through the diagram
// Worker. Node supplies none, so here its Load lists that recovery; in the
// browser it loads without a row (OV-751). This Load port accepts only that
// row: any other row cancels the Load, so a test still sees an unexpected
// repair.
export const acceptNodeFeatureRecoveryRow = async ({ items }) => items.every(
  (line) => line.startsWith('Features: the feature information of the saved diagram could not be restored')
);
