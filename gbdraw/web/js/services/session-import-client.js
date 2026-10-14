// @ts-check
import { createBoundedJsonReceiver } from './bounded-json-transport.js';
import { recordSessionLifecycleEvent } from './runtime-test-hooks.js';
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { BoundedJsonPart, BoundedJsonValue } from './bounded-json-transport.js' */

/**
 * What the Worker measured, in milliseconds.
 * @typedef {{ readMs: number, parseMs: number, replyMs: number }} SessionImportTimings
 */
/**
 * The reassembled, still unvalidated JSON of a Session file.
 * @typedef {{ data: BoundedJsonValue, characters: number, timings: SessionImportTimings }} SessionImportCandidate
 */
/**
 * A message of the Session import Worker: a transport part, the final result, or a diagnostic.
 * @typedef {(BoundedJsonPart & { operationId: number, status: 'part' })
 *   | { operationId: number, status: 'ok', characters: number, timings: SessionImportTimings }
 *   | { operationId: number, status: 'error', error?: { code?: string, stage?: string, context?: Record<string, unknown> } }
 * } SessionImportReply
 */

let nextOperationId = 1;

// A Worker that cannot start, crashes or replies unreadably leaves no import path.
const unavailable = () => diagnosticError('SESSION_IMPORT_UNAVAILABLE', {}, { stage: 'transport' });

// This owns transport lifetime only. The caller still admits an untrusted candidate.
/**
 * @param {File} file
 * @param {{ signal?: AbortSignal }} [options]
 * @returns {Promise<SessionImportCandidate>}
 */
export const importSessionFile = (file, { signal } = {}) => new Promise((resolve, reject) => {
  const operationId = nextOperationId++;
  /** @type {Worker | null} */
  let worker = null;
  let settled = false;
  const assembly = createBoundedJsonReceiver();
  // Ends the operation once; false when it had already ended.
  const release = () => {
    if (settled) return false;
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
    return true;
  };
  /** @param {Error} error */
  const fail = (error) => {
    if (!release()) return;
    recordSessionLifecycleEvent('session-import-worker-settled', { operationId, status: 'error' });
    reject(error);
  };
  /** @param {SessionImportCandidate} reply */
  const succeed = (reply) => {
    if (!release()) return;
    recordSessionLifecycleEvent('session-import-worker-settled', { operationId, status: 'ok' });
    resolve(reply);
  };
  const cancel = () => fail(Object.assign(
    new Error('Session loading was canceled.'), { name: 'AbortError', code: 'SESSION_IMPORT_CANCELED', stage: 'transport' }
  ));
  /** @param {ErrorEvent} event */
  const crash = (event) => {
    event.preventDefault?.();
    fail(unavailable());
  };
  const unreadable = () => fail(unavailable());
  /** @param {MessageEvent<SessionImportReply>} event */
  const receive = ({ data: reply }) => {
    if (reply?.operationId !== operationId || settled) return;
    if (reply.status === 'part') {
      try {
        assembly.receivePart(reply);
        // `settled` is false here, and only `release` clears `worker`, so the started Worker is still held.
        /** @type {Worker} */ (worker).postMessage({ operationId, ack: true });
      } catch { unreadable(); }
    } else if (reply.status === 'error') {
      const { code = 'UNKNOWN', stage = 'transport', context = {} } = reply.error || {};
      fail(diagnosticError(code, context, { stage }));
    } else if (reply.status === 'ok' && Number.isSafeInteger(reply.characters) && reply.characters >= 0) {
      try {
        succeed({ data: assembly.getValue(), characters: reply.characters, timings: reply.timings });
      } catch { unreadable(); }
    } else {
      unreadable();
    }
  };
  if (signal?.aborted) { cancel(); return; }
  if (typeof Worker !== 'function') {
    fail(unavailable());
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
    fail(unavailable());
  }
});
