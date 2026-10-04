import assert from 'node:assert/strict';
import test from 'node:test';
import { sendBoundedJson } from '../../gbdraw/web/js/services/bounded-json-transport.js';
import { EXPECTED_WEB_RUNTIME_CAPABILITIES } from '../../gbdraw/web/js/services/runtime-capabilities.js';

globalThis.location = { href: 'https://example.test/gbdraw/web/' };
class ReplyWorker {
  static instances = [];
  static malformed = false;
  constructor() { this.listeners = new Map(); this.acks = 0; ReplyWorker.instances.push(this); }
  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, new Set());
    this.listeners.get(type).add(listener);
  }
  removeEventListener(type, listener) { this.listeners.get(type)?.delete(listener); }
  emit(data) { this.listeners.get('message')?.forEach((listener) => listener({ data })); }
  terminate() { this.terminated = true; }
  postMessage(message) {
    if (message.type === 'init') {
      queueMicrotask(() => this.emit({ type: 'init', id: message.id, ok: true,
        capabilities: EXPECTED_WEB_RUNTIME_CAPABILITIES }));
    } else if (message.type === 'auxiliary-ack') {
      this.acks++;
      this.ack?.();
    } else if (message.type === 'helper') {
      const envelope = { type: 'helper', requestId: message.requestId };
      queueMicrotask(async () => {
        if (ReplyWorker.malformed) {
          this.emit({ ...envelope, status: 'part', kind: 'batch', path: [], index: 7, value: [] });
          return;
        }
        this.emit({ ...envelope, requestId: message.requestId + 1, status: 'part', kind: 'value', path: [], value: 'stale' });
        await sendBoundedJson(message.payload.expected, (part, transfer = []) => new Promise((resolve) => {
          this.ack = resolve;
          this.emit(structuredClone({ ...envelope, status: 'part', ...part }, { transfer }));
        }));
        this.emit({ ...envelope, ok: true });
        this.emit({ ...envelope, ok: true });
      });
    }
  }
}
globalThis.Worker = ReplyWorker;
const { runDiagramHelperOperation, disposeDiagramGenerationWorker, DIAGRAM_HELPER_OPERATIONS } = await import('../../gbdraw/web/js/services/diagram-generation.js');

test('the real helper client assembles large replies privately and releases listeners once', async () => {
  const expected = { text: 'α😀\ud800'.repeat(80000), rows: Array.from({ length: 2000 }, (_, index) => ({ index })) };
  try {
    const response = await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EXTRACT_FIRST_FASTA, { expected });
    assert.deepEqual(response.result, expected);
    const worker = ReplyWorker.instances.at(-1);
    assert.ok(worker.acks > 3);
    for (const listeners of worker.listeners.values()) assert.equal(listeners.size, 0);
  } finally { disposeDiagramGenerationWorker(); }
});

test('a malformed helper part rejects without admitting a partial result or acknowledging it', async () => {
  ReplyWorker.malformed = true;
  try {
    await assert.rejects(runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.EXTRACT_FIRST_FASTA), { code: 'RESULT_INVALID', stage: 'result-admission' });
    const worker = ReplyWorker.instances.at(-1);
    assert.equal(worker.acks, 0);
    assert.equal(worker.terminated, true, 'release the sender waiting for a part acknowledgement');
    for (const listeners of worker.listeners.values()) assert.equal(listeners.size, 0);
  } finally { ReplyWorker.malformed = false; disposeDiagramGenerationWorker(); }
});
