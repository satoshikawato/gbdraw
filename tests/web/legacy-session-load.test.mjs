// Load of a Session 27-44 (DESIGN 4.3): the stored figure shows at once,
// Python moves the file to a Session 49 in one Worker call, and Load reads
// Python's document as it reads a Session 49 file, with the stored Results
// a Session 49 cannot hold shown in memory. A Worker failure installs nothing.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { evaluatePython, pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { body: { appendChild: () => {} }, createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null }) };
globalThis.alert = () => {};
const { installFakeSvgDom } = await import('./fake-svg-dom.mjs');
installFakeSvgDom();
installSessionImportWorker();

const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { setLegacySessionUpgrader } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}.gbdraw-session.json.gz`, import.meta.url))));
const pythonUpgrade = async (document) => JSON.parse(await evaluatePython(`
from gbdraw.web_support.legacy_session import legacy_session_to_drawings_json
result = legacy_session_to_drawings_json(payload)
`, JSON.stringify(document)));
const workerCalls = [];
setLegacySessionUpgrader(pythonLegacySessionUpgrader(workerCalls));
const load = (session) => importSession({
  target: { files: [new Blob([JSON.stringify(session)], { type: 'application/json' })], value: 'selected' }
}, {
  confirmLoadRepairs: async () => true,
  // The composition root's port; these Results need no change.
  transformLegacyResultSvg: Object.assign(() => false, { appliesToContent: () => false })
});
const drawingsShown = () => state.drawingList().map(({ id, mode }) => ({ id, mode }));

test('a Session 44 Load makes one upgradeLegacySession call and installs the drawings of Python\'s Session 49', async () => {
  const session = fixture('HmmtDNA_basic_circular.v44-schema8');
  const expected = await pythonUpgrade(session);
  workerCalls.length = 0;
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(workerCalls.length, 1);
  assert.deepEqual(drawingsShown(), expected.document.drawings.map(({ id, mode }) => ({ id, mode })));
  assert.equal(state.activeDrawingId.value, expected.document.ui.activeDrawingId);
  assert.equal(state.results.value.length, 1);
});

test('a Session 39 Load shows the stored Results Python returns beside its document, in the drawing that opens', async () => {
  const session = fixture('BGC0000708-BGC0000713.v39');
  const expected = await pythonUpgrade(session);
  workerCalls.length = 0;
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.equal(workerCalls.length, 1);
  assert.equal(state.activeDrawingId.value, expected.document.ui.activeDrawingId);
  assert.deepEqual(Object.keys(expected.storedResults), [expected.document.ui.activeDrawingId]);
  assert.equal(state.results.value.length, expected.storedResults[expected.document.ui.activeDrawingId].length);
});

test('a Worker failure of a Session 27-44 Load installs nothing, withdraws the stored figure and reports one diagnostic', async () => {
  const before = await load(fixture('HmmtDNA_basic_circular.v44-schema8'));
  assert.equal(before.status, 'ok');
  const drawings = drawingsShown();
  const results = state.results.value;
  setLegacySessionUpgrader(async () => { throw { code: 'HELPER_PROTOCOL', stage: 'helper', operation: 'upgradeLegacySession' }; });
  try {
    const loaded = await load(fixture('BGC0000708-BGC0000713.v39'));
    assert.equal(loaded.status, 'error');
    assert.equal(loaded.error.code, 'HELPER_PROTOCOL');
    assert.equal(loaded.error.operation, 'upgradeLegacySession');
    assert.equal(state.errorLog.value, loaded.error);
    assert.deepEqual(drawingsShown(), drawings);
    assert.equal(state.results.value, results);
  } finally {
    setLegacySessionUpgrader(pythonLegacySessionUpgrader(workerCalls));
  }
});
