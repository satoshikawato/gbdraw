import { readSessionText } from '../services/session-file.js';

self.addEventListener('message', async ({ data: { operationId, file } }) => {
  let stage = 'read';
  try {
    const started = performance.now();
    const text = await readSessionText(file);
    const readCompleted = performance.now();
    stage = 'parse';
    const data = JSON.parse(text);
    const parseCompleted = performance.now();
    stage = 'reply';
    self.postMessage({
      operationId, status: 'ok', data, characters: text.length,
      timings: { readMs: readCompleted - started, parseMs: parseCompleted - readCompleted }
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
