import { recordSessionLifecycleEvent } from './runtime-test-hooks.js';

let nextOperationId = 1;

const importError = (message, code, stage, name = 'Error') => Object.assign(
  new Error(message), { name, code, stage }
);

// This owns transport lifetime only. The caller still admits an untrusted candidate.
export const importSessionFile = (file, { signal } = {}) => new Promise((resolve, reject) => {
  const operationId = nextOperationId++;
  let worker = null;
  let settled = false;
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
  const cancel = () => settle(importError(
    'Session loading was canceled.', 'SESSION_IMPORT_CANCELED', 'transport', 'AbortError'
  ));
  const crash = (event) => {
    event.preventDefault?.();
    settle(importError('Session import Worker failed.', 'SESSION_IMPORT_CRASH', 'transport'));
  };
  const unreadable = () => settle(importError(
    'Session import Worker reply could not be read.', 'SESSION_IMPORT_REPLY_FAILED', 'transport'
  ));
  const receive = ({ data: reply }) => {
    if (reply?.operationId !== operationId || settled) return;
    if (reply.status === 'error') {
      settle(importError(
        reply.error?.message || 'Session import failed.',
        reply.error?.code || 'SESSION_IMPORT_FAILED',
        reply.error?.stage || 'transport',
        reply.error?.name || 'Error'
      ));
    } else if (reply.status === 'ok' && Object.hasOwn(reply, 'data')
      && Number.isSafeInteger(reply.characters) && reply.characters >= 0) {
      settle(null, { data: reply.data, characters: reply.characters, timings: reply.timings });
    } else {
      unreadable();
    }
  };
  if (signal?.aborted) { cancel(); return; }
  if (typeof Worker !== 'function') {
    settle(importError(
      'This browser does not support Session import Workers.',
      'SESSION_IMPORT_UNAVAILABLE', 'transport'
    ));
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
    settle(importError(
      'Session import Worker could not be started.', 'SESSION_IMPORT_START_FAILED', 'transport'
    ));
  }
});
