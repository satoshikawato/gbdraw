import { readSessionText } from '../services/session-file.js';
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
    const text = await readSessionText(file);
    const readCompleted = performance.now();
    stage = 'parse';
    const data = JSON.parse(text);
    const parseCompleted = performance.now();
    stage = 'reply';
    const replyStarted = performance.now();
    await sendBoundedJson(data, (part, transfers) => sendPart(operationId, part, transfers), [], { consume: true });
    self.postMessage({
      operationId, status: 'ok', characters: text.length,
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
