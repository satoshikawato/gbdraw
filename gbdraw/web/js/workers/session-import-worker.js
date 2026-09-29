import { readSessionText } from '../services/session-file.js';

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
    self.postMessage({
      operationId, status: 'error',
      error: {
        name: error?.name || 'Error',
        code: `SESSION_IMPORT_${stage.toUpperCase()}_FAILED`, stage,
        // JSON SyntaxError messages can quote uploaded text.
        message: stage === 'parse' ? 'Invalid JSON structure.'
          : stage === 'reply' ? 'Session import Worker reply failed.'
            : error?.message || 'Session file could not be read.'
      }
    });
  }
});
