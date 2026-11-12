import { MAX_EXPANDED_SESSION_BYTES, readSessionText } from '../services/session-file.js';
import { assertSafeObjectKeys, assertUniqueJsonKeys } from '../services/safe-object-keys.js';
import { isSessionContainerFile, readSessionContainer } from '../services/session-container.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';

import { sendBoundedJson } from '../services/bounded-json-transport.js';

let acknowledge = null;
const sendPart = (operationId, part, transfers = []) => new Promise((resolve) => {
  acknowledge = { operationId, resolve };
  self.postMessage({ operationId, status: 'part', ...part }, transfers);
});

self.addEventListener('message', async ({ data: { operationId, file, ack } }) => {
  if (ack) {
    if (acknowledge?.operationId === operationId) {
      const current = acknowledge;
      acknowledge = null;
      current.resolve();
    }
    return;
  }
  let stage = 'read';
  try {
    const started = performance.now();
    // Content decides (receipt session.file-container): a container's
    // directory and document are read here, and its entries are opened by the
    // page from the same file; any other file is a JSON Session, plain or gzip.
    const container = await isSessionContainerFile(file)
      ? await readSessionContainer(file, {
        // The browser loads at most 512 MiB of settings and results; the inputs have no cap.
        admitDocument: (sessionBytes) => {
          if (sessionBytes > MAX_EXPANDED_SESSION_BYTES) {
            throw diagnosticError('SESSION_SIZE_LIMIT', { sessionBytes, sessionByteLimit: MAX_EXPANDED_SESSION_BYTES },
              { stage: 'read', operation: 'session-load' });
          }
        }
      })
      : null;
    const text = container ? '' : await readSessionText(file);
    const readCompleted = performance.now();
    stage = 'parse';
    const data = container ? container.document : JSON.parse(text);
    // Python refuses a Session with one key twice as invalid JSON (OV-2694); the container reader checked its document.
    if (!container) assertUniqueJsonKeys(text);
    const parseCompleted = performance.now();
    // Session Load's one check of unsafe own keys (OV-952), off the main thread; its diagnostic names its stage.
    assertSafeObjectKeys(data);
    stage = 'reply';
    const replyStarted = performance.now();
    await sendBoundedJson(data, (part, transfers) => sendPart(operationId, part, transfers), { consume: true });
    self.postMessage({
      operationId, status: 'ok', characters: container ? container.documentEntry.size : text.length,
      // The entries the page streams into the input copies, by resource id (plain data).
      container: container ? { resources: container.resources, inputBytes: container.inputBytes } : null,
      timings: { readMs: readCompleted - started, parseMs: parseCompleted - readCompleted,
        replyMs: performance.now() - replyStarted }
    });
  } catch (error) {
    // The final diagnostic code and bounded context; never the cause message,
    // which can quote uploaded text.
    const known = normalizeUserFacingError(error, { stage });
    self.postMessage({
      operationId, status: 'error',
      error: known.code !== 'UNKNOWN' ? { code: known.code, stage: known.stage, context: known.context }
        : stage === 'read' ? { code: 'INPUT_UNREADABLE', stage, context: { field: 'schema' } }
          : stage === 'parse' ? { code: 'INPUT_INVALID', stage, context: { field: 'schema', reason: 'JSON_FORMAT' } }
            : { code: 'SESSION_IMPORT_UNAVAILABLE', stage: 'transport', context: {} }
    });
  }
});
