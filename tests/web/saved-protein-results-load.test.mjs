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
// `dialogs` collects each Load repair dialog the Load shows.
const load = async (session, dialogs = []) => importSession(await sessionLoadEvent(session), {
  confirmLoadRepairs: async (dialog) => { dialogs.push(dialog); return true; },
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
  // Sessions with saved LOSAT rows but no LOSATP search: committed protein results (mode none), a collinear
  // result, a nucleotide BLAST table, and a Circular drawing.
  for (const name of [
    'BGC0000708-BGC0000713.v47', 'majanivirus_collinear-rc.v47', 'two-mode-linear-shown.v48', 'web-losat-ring-copies.v50'
  ]) {
    const session = fixture(name);
    assert.ok(session.losatCache.entries.length > 0, name);
    workerCalls.length = 0;
    await load(session);
    assert.deepEqual(workerCalls, [], name);
  }
});

test('a saved-rows step Python refuses keeps the search and lists a Load row; a protocol failure stops the Load', async () => {
  // load_session_document is stricter than the Web reader: without the protein identity manifest Python refuses
  // the file, which the Web loads as before, with its search.
  const session = fixture('BGC0000711-BGC0000713-losatp-cli.v47');
  delete session.proteinIdentityManifest;
  workerCalls.length = 0;
  const dialogs = [];
  const loaded = await load(session, dialogs);
  assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
  assert.deepEqual(workerCalls, ['storeSavedProteinResults']);
  assert.deepEqual(comparisonKinds(sessionConfig.getCommittedCanonicalSession().renderRequest), [['generatedProteinComparison', 'orthogroup']]);
  assert.deepEqual(dialogs.flatMap(({ items }) => items), [
    'Saved diagram: the LOSAT results saved for its protein comparison (losatCache) could not be read. Click Generate to run the comparison again before you edit the diagram.'
  ]);
  setLegacySessionUpgrader(async () => { throw { code: 'HELPER_PROTOCOL', stage: 'helper', operation: 'storeSavedProteinResults' }; });
  try {
    const failed = await load(fixture('BGC0000711-BGC0000713-losatp-cli.v47'));
    assert.equal(failed.status, 'error');
    assert.equal(failed.error.code, 'HELPER_PROTOCOL');
  } finally {
    setLegacySessionUpgrader(pythonLegacySessionUpgrader([], workerCalls));
  }
});

test('the Load dialog lists what Load could not restore of the stored request the step replaced', async () => {
  // Python's reply for the file as written; the file Load reads holds an alignment the Web reader cannot read.
  const session = fixture('BGC0000711-BGC0000713-losatp-cli.v47');
  const reply = await pythonLegacySessionUpgrader()('storeSavedProteinResults', { document: session });
  session.renderRequest.layout.similarityAlignment = 'unreadable';
  setLegacySessionUpgrader(async () => reply);
  try {
    const dialogs = [];
    const loaded = await load(session, dialogs);
    assert.equal(loaded.status, 'ok', JSON.stringify(loaded.error));
    assert.deepEqual(comparisonKinds(sessionConfig.getCommittedCanonicalSession().renderRequest)[0], ['precomputedProteinComparison', null]);
    const items = dialogs.flatMap((dialog) => dialog.items);
    assert.deepEqual(items, ['Settings: 1 entry of renderRequest.layout could not be read and was dropped.']);
  } finally {
    setLegacySessionUpgrader(pythonLegacySessionUpgrader([], workerCalls));
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
