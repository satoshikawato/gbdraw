import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import test from 'node:test';

// G-G(3) (Web GUI audit 2026-09-30): the shared Generate and Save helpers fail
// on an UNKNOWN diagnostic, an uncaught page error, or NaN/undefined in Run
// Info, and accept an opt-out only when it names an audit ID.
const require = createRequire(import.meta.url);
const { assertOperationHealth } = require('./helpers/app-lifecycle.cjs');

const fakePage = (app) => {
  const listeners = new Map();
  return {
    app,
    on: (type, listener) => listeners.set(type, listener),
    emit: (type, value) => listeners.get(type)?.(value),
    evaluate: async (callback) => {
      globalThis.window = { __GBDRAW_APP__: app };
      try {
        return callback();
      } finally {
        delete globalThis.window;
      }
    }
  };
};

test('a clean operation passes and each health failure is reported', async () => {
  const page = fakePage({ errorLog: null, lastRunInfo: { args: ['-w', '1000'], windowSize: 1000 } });
  assert.deepEqual(await assertOperationHealth(page), { errorCode: null, pageErrors: [], invalidRunInfo: [] });

  page.app.errorLog = { code: 'UNKNOWN', operation: 'generate', stage: 'request-validation', summary: 'failed' };
  await assert.rejects(assertOperationHealth(page), /must not report an UNKNOWN diagnostic/);
  page.app.errorLog = null;

  page.app.lastRunInfo = { args: ['-w', 'NaN'], stepSize: Number.NaN, nested: [{ label: 'undefined' }] };
  await assert.rejects(assertOperationHealth(page), (error) => {
    assert.match(error.message, /must not record NaN or undefined in Run Info/);
    assert.match(error.message, /lastRunInfo\.args\[1\]/);
    assert.match(error.message, /lastRunInfo\.stepSize=NaN/);
    assert.match(error.message, /lastRunInfo\.nested\[0\]\.label/);
    return true;
  });
  page.app.lastRunInfo = null;

  page.emit('pageerror', new Error('boom'));
  await assert.rejects(assertOperationHealth(page), /must not raise an uncaught page error/);
});

test('an opt-out names the audit ID of a known defect', async () => {
  const page = fakePage({
    errorLog: { code: 'UNKNOWN', operation: 'save-session', stage: 'save', summary: 'failed' },
    lastRunInfo: { windowSize: Number.NaN }
  });
  const health = await assertOperationHealth(page, { allowUnknown: 'IN-08', allowInvalidRunInfo: 'TR-04' });
  assert.equal(health.errorCode, 'UNKNOWN');
  assert.deepEqual(health.invalidRunInfo, ['lastRunInfo.windowSize=NaN']);
  for (const value of [true, 'unknown', 'IN08']) {
    await assert.rejects(assertOperationHealth(page, { allowUnknown: value }), /must name the audit ID/);
  }
});

test('an alert left by an earlier operation is not attributed to Save', async () => {
  const errorLog = { code: 'UNKNOWN', operation: 'generate', stage: 'render', summary: 'earlier' };
  const page = fakePage({ errorLog, lastRunInfo: null });
  const errorSignatureBefore = JSON.stringify([errorLog.code, errorLog.operation, errorLog.stage, errorLog.summary]);
  assert.equal((await assertOperationHealth(page, { errorSignatureBefore })).errorCode, null);
  page.app.errorLog = { ...errorLog, operation: 'save-session', summary: 'new' };
  await assert.rejects(assertOperationHealth(page, { errorSignatureBefore }), /UNKNOWN diagnostic/);
});
