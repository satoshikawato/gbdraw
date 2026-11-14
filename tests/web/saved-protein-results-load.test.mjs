// Load of a Session 46-50 whose Linear request holds a LOSATP search (OV-3545): the browser's render has no LOSAT,
// so Load gives the committed request the results of the LOSAT rows the file saved, in one Worker call, as a
// Generate commits them; a redraw before the first Generate then needs no search. Other Sessions make no call.
import assert from 'node:assert/strict';
import { test } from 'node:test';
import { readFileSync } from 'node:fs';
import { gunzipSync } from 'node:zlib';
import { installSessionImportWorker } from './helpers/session-import-node.mjs';
import { callerSession } from './helpers/session-caller.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';
import { pythonLegacySessionUpgrader, validateWebConfigOverrides } from './helpers/python-evaluator.mjs';
import { sessionLoadEvent } from './helpers/session-documents.mjs';

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
const { setLegacySessionUpgrader, readSavedProteinResults } = await import('../../gbdraw/web/js/services/legacy-session-upgrade.js');
sessionConfig.setUnmanagedConfigOverrideValidator(async (payload) => ({ result: await validateWebConfigOverrides(payload) }));
const { state } = await import('../../gbdraw/web/js/state.js');
const { importSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));

const fixture = (name) => JSON.parse(gunzipSync(readFileSync(new URL(`../fixtures/sessions/${name}.gbdraw-session.json.gz`, import.meta.url))));
// The composition root's port of Load's Python steps, with the operation of each call.
const workerCalls = [];
setLegacySessionUpgrader(pythonLegacySessionUpgrader([], workerCalls));
const load = async (session) => importSession(await sessionLoadEvent(session), {
  confirmLoadRepairs: async () => true,
  transformLegacyResultSvg: Object.assign(() => false, { appliesToContent: () => false })
});
const comparisonKinds = (request) => request.comparisons.map(({ kind, mode }) => [kind, mode ?? null]);

test('a CLI Session 47 Load commits its LOSATP search with the results of its saved rows', async () => {
  const session = fixture('BGC0000711-BGC0000713-losatp-cli.v47');
  assert.deepEqual(comparisonKinds(session.renderRequest), [['generatedProteinComparison', 'orthogroup']]);
  workerCalls.length = 0;
  const loaded = await load(session);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(workerCalls, ['storeSavedProteinResults']);
  const committed = sessionConfig.getCommittedCanonicalSession();
  assert.deepEqual(comparisonKinds(committed.renderRequest), [
    ['precomputedProteinComparison', null], ['orthogroupResult', null], ['generatedProteinComparison', 'none']
  ]);
  // The tables the request names are held with the Session's inputs, which keep their IDs.
  const named = committed.renderRequest.comparisons.flatMap(({ resourceId }) => (resourceId ? [resourceId] : []));
  assert.equal(named.length, 2);
  for (const resourceId of [...named, ...Object.keys(session.resources)]) assert.ok(committed.resources[resourceId], resourceId);
});

test('a Session 46-50 without a LOSATP search, or whose search has no saved rows, makes no Worker call', async () => {
  // A Linear request without protein comparisons, and a Pairwise search saved before any run (OV-3541's case).
  for (const name of ['cli-linear-blacklist-file.v47', 'python-api-linear-pairs.v50']) {
    workerCalls.length = 0;
    const loaded = await load(fixture(name));
    assert.equal(loaded.status, 'ok', `${name}: ${JSON.stringify(loaded.error)}`);
    assert.deepEqual(workerCalls, [], name);
  }
});

test('a reply of another shape, one that adds a resource the file holds, or one naming a drawing Load did not send fails',
  async () => {
    const session = fixture('BGC0000711-BGC0000713-losatp-cli.v47');
    const failure = { code: 'RESULT_INVALID', stage: 'helper', operation: 'storeSavedProteinResults' };
    const reply = (result) => ({ upgrade: async () => ({ result }) });
    await assert.rejects(readSavedProteinResults(session, reply({
      requests: { linear: {} }, resources: { [Object.keys(session.resources)[0]]: {} }
    })), failure);
    await assert.rejects(readSavedProteinResults(session, reply(null)), failure);
    await assert.rejects(readSavedProteinResults(session, { upgrade: null }), { ...failure, code: 'HELPER_PROTOCOL' });
    setLegacySessionUpgrader(async () => ({ result: { requests: { circular: {} }, resources: {} } }));
    try {
      const loaded = await load(session);
      assert.equal(loaded.status, 'error');
      assert.deepEqual({ code: loaded.error.code, operation: loaded.error.operation }, { code: 'RESULT_INVALID', operation: 'storeSavedProteinResults' });
    } finally {
      setLegacySessionUpgrader(pythonLegacySessionUpgrader([], workerCalls));
    }
  });
