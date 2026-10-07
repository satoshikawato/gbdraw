// @ts-check
import { createBoundedJsonReceiver } from './bounded-json-transport.js';
import { recordSessionLifecycleEvent } from './runtime-test-hooks.js';
import { diagnosticError } from './error-normalization.js';

let nextOperationId = 1;

// A Worker that cannot start, crashes or replies unreadably leaves no import path.
const unavailable = () => diagnosticError('SESSION_IMPORT_UNAVAILABLE', {}, { stage: 'transport' });

// This owns transport lifetime only. The caller still admits an untrusted candidate.
/**
 * @param {File} file
 * @param {{ signal?: AbortSignal }} [options]
 */
export const importSessionFile = (file, { signal } = {}) => new Promise((resolve, reject) => {
  const operationId = nextOperationId++;
  /** @type {Worker | null} */
  let worker = null;
  let settled = false;
  const assembly = createBoundedJsonReceiver();
  const settle = (error, reply) => {
    if (settled) return;
    settled = true;
    signal?.removeEventListener('abort', cancel);
    if (worker) {
      worker.removeEventListener('message', receive);
      worker.removeEventListener('error', crash);
      worker.removeEventListener('messageerror', unreadable);
      worker.terminate();
      recordSessionLifecycleEvent('session-import-worker-terminated', { operationId });
      worker = null;
    }
    recordSessionLifecycleEvent('session-import-worker-settled', {
      operationId, status: error ? 'error' : 'ok'
    });
    if (error) reject(error);
    else resolve(reply);
  };
  const cancel = () => settle(Object.assign(
    new Error('Session loading was canceled.'), { name: 'AbortError', code: 'SESSION_IMPORT_CANCELED', stage: 'transport' }
  ));
  const crash = (event) => {
    event.preventDefault?.();
    settle(unavailable());
  };
  const unreadable = () => settle(unavailable());
  const receive = ({ data: reply }) => {
    if (reply?.operationId !== operationId || settled) return;
    if (reply.status === 'part') {
      try {
        assembly.receivePart(reply);
        // `settled` is false here, and only `settle` clears `worker`, so the started Worker is still held.
        /** @type {Worker} */ (worker).postMessage({ operationId, ack: true });
      } catch { unreadable(); }
    } else if (reply.status === 'error') {
      const { code = 'UNKNOWN', stage = 'transport', context = {} } = reply.error || {};
      settle(diagnosticError(code, context, { stage }));
    } else if (reply.status === 'ok' && Number.isSafeInteger(reply.characters) && reply.characters >= 0) {
      try {
        settle(null, { data: assembly.getValue(), characters: reply.characters, timings: reply.timings });
      } catch { unreadable(); }
    } else {
      unreadable();
    }
  };
  if (signal?.aborted) { cancel(); return; }
  if (typeof Worker !== 'function') {
    settle(unavailable());
    return;
  }
  try {
    worker = new Worker(new URL('../workers/session-import-worker.js', import.meta.url), { type: 'module' });
    worker.addEventListener('message', receive);
    worker.addEventListener('error', crash);
    worker.addEventListener('messageerror', unreadable);
    signal?.addEventListener('abort', cancel, { once: true });
    recordSessionLifecycleEvent('session-import-worker-start', { operationId, fileBytes: file.size });
    worker.postMessage({ operationId, file });
  } catch {
    settle(unavailable());
  }
});
