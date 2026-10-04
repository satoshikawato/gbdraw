import { Worker as NodeWorker } from 'node:worker_threads';

// Browser event API around a real Node isolate running the production handler.
export class NodeSessionImportWorker {
  constructor(url) {
    if (!String(url).endsWith('/session-import-worker.js')) {
      throw new Error('The Node import test adapter does not supply a diagram Worker.');
    }
    this.listeners = new Map();
    this.worker = new NodeWorker(`
      const {parentPort}=require('node:worker_threads');
      global.self={addEventListener:(_,fn)=>parentPort.on('message',data=>fn({data})),
        postMessage:data=>parentPort.postMessage(data)};
      import(${JSON.stringify(String(url))}).then(()=>parentPort.postMessage({ready:true}));
    `, { eval: true });
    this.ready = new Promise(resolve => {
      this.worker.on('message', data => {
        if (data.ready) resolve();
        else this.emit('message', { data });
      });
    });
    this.worker.on('error', error => this.emit('error', { message: error.message }));
    this.worker.on('messageerror', () => this.emit('messageerror', {}));
  }
  emit(name, event) { for (const fn of this.listeners.get(name) || []) fn(event); }
  addEventListener(name, fn) {
    if (!this.listeners.has(name)) this.listeners.set(name, new Set());
    this.listeners.get(name).add(fn);
  }
  removeEventListener(name, fn) { this.listeners.get(name)?.delete(fn); }
  postMessage(message) { this.ready.then(() => this.worker.postMessage(message)); }
  terminate() { this.listeners.clear(); this.worker.terminate(); }
}

export const installSessionImportWorker = () => { globalThis.Worker = NodeSessionImportWorker; };
