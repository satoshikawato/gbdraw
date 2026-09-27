import assert from 'node:assert/strict';
import test from 'node:test';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const root = await mkdtemp(join(tmpdir(), 'gbdraw-session-operations-'));
await cp(join(process.cwd(), 'gbdraw/web/js'), join(root, 'js'), { recursive: true });
await writeFile(join(root, 'package.json'), '{"type":"module"}');
globalThis.window = { Vue: {
  ref: value => ({ value }), reactive: value => value,
  computed: getter => ({ get value() { return getter(); } }), nextTick: async () => {}
} };
const load = path => import(pathToFileURL(join(root, 'js', path)));
const { state, sessionOperationAvailability } = await load('state.js');
const { createHistoryManager } = await load('services/history.js');
const { createAnnotationEditor } = await load('app/annotations.js');
const { createLinearRecordSelector } = await load('app/linear-record-selector.js');
const { createRulePreparation, ruleMatchesFeature } = await load('app/rule-matching.js');
const gate = () => {
  let release;
  const promise = new Promise(resolve => { release = resolve; });
  return { promise, release };
};
const saving = { status: 'busy', reason: 'Saving session. Retry after saving finishes.' };

test('availability derives Session exclusivity and preparation from existing owners', () => {
  assert.equal(sessionOperationAvailability(), null);
  state.processing.value = true;
  assert.match(sessionOperationAvailability('save').reason, /Generating/);
  assert.match(sessionOperationAvailability('load').reason, /Generating/);
  assert.equal(sessionOperationAvailability(), null);
  state.processing.value = false;
  state.labelReflowProcessing.value = true;
  assert.match(sessionOperationAvailability('save').reason, /Updating/);
  state.labelReflowProcessing.value = false;
  state.sessionSavePending.value = true;
  assert.deepEqual(sessionOperationAvailability(), saving);
  state.sessionImportPending.value = true;
  assert.match(sessionOperationAvailability('save').reason, /Loading/);
  state.sessionImportPending.value = false;
  state.sessionSavePending.value = false;
});

const historyFor = (options = {}) => createHistoryManager({
  buildIntent: async () => ({ value: 1 }), applyIntent: async () => {},
  buildCheckpoint: async () => ({ value: 1 }), applyCheckpoint: async () => {},
  mutationAvailability: sessionOperationAvailability, ...options
});

test('History admits no late command or action after an awaited capture', async () => {
  const capture = gate();
  const history = historyFor({ buildIntent: () => capture.promise });
  let mutations = 0;
  const action = history.runUndoable('late edit', () => { mutations += 1; });
  state.sessionSavePending.value = true;
  capture.release({ value: 1 });
  assert.deepEqual(await action, saving);
  assert.equal(mutations, 0);
  assert.equal(history.getUndoCount(), 0);
  state.sessionSavePending.value = false;
  const command = gate();
  const pending = history.runUndoableCommand('late command', () => command.promise);
  state.sessionSavePending.value = true;
  command.release({ changes: [] });
  assert.deepEqual(await pending, saving);
  state.sessionSavePending.value = false;
});

test('Save and Load report an active asynchronous History mutation until settlement', async () => {
  const history = historyFor();
  state.sessionPreparationBusyReason = () => history.mutationPending()
    ? 'Applying an edit. Retry after the edit finishes.' : '';
  const work = gate();
  let started = false;
  const pending = history.runUndoable('asynchronous edit', async () => {
    started = true;
    await work.promise;
  });
  while (!started) await Promise.resolve();
  assert.match(sessionOperationAvailability('save').reason, /Applying/);
  assert.match(sessionOperationAvailability('load').reason, /Applying/);
  work.release();
  await pending;
  assert.equal(sessionOperationAvailability('save'), null);
  state.sessionPreparationBusyReason = null;
});

test('focused input intent does not occupy Session until asynchronous mutation starts', async () => {
  const history = historyFor();
  const tx = await history.begin('focused input', { source: 'input-adapter' });
  assert.equal(history.mutationPending(), false);
  await history.commit(tx);
});

test('canceled Load baseline does not clear the previous History', async () => {
  let value = 1;
  let read = null;
  const history = historyFor({ buildIntent: () => read ? read.promise : { value } });
  await history.runUndoable('existing edit', () => { value = 2; });
  assert.equal(history.getUndoCount(), 1);
  read = gate();
  let current = true;
  const pending = history.initializeIntentBaseline('Loaded session', { isCurrent: () => current });
  current = false;
  read.release({ value: 3 });
  assert.equal(await pending, false);
  assert.equal(history.getUndoCount(), 1);
});

test('late annotation reads return busy without replacing the existing rows', async () => {
  const read = gate();
  const sets = [];
  const owner = createAnnotationEditor({
    state: { annotationSets: sets, adv: {}, sessionOperationAvailability },
    getRecordCatalog: () => ({ status: 'ready', records: [] })
  });
  const file = { name: 'annotations.tsv', text: () => read.promise };
  const input = { files: [file], value: 'annotations.tsv' };
  const pending = owner.importAnnotationTableFile({ target: input });
  state.sessionSavePending.value = true;
  read.release('set_id\tid\tstart\tend\nreview\twindow\t1\t10\n');
  assert.deepEqual(await pending, saving);
  assert.deepEqual(sets, []);
  assert.equal(input.value, '');
  state.sessionSavePending.value = false;
});

test('late source discovery cannot expand the live source set during Session', async () => {
  const read = gate();
  let expanded = 0;
  const file = {};
  const source = { uid: 'source', gb: file, region_record_id: '' };
  const selector = createLinearRecordSelector({
    state: { mode: { value: 'linear' }, lInputType: { value: 'gb' },
      linearSeqs: [source], sessionOperationAvailability },
    reactive: value => value, recordReader: () => read.promise,
    onRecordsDiscovered: () => { expanded += 1; }
  });
  const pending = selector.refresh();
  state.sessionSavePending.value = true;
  read.release([{ selector: '#1', recordId: 'first', recordLength: 100 }]);
  assert.deepEqual(await pending, saving);
  assert.equal(expanded, 0);
  assert.equal(selector.statusFor(source), 'deferred');
  state.sessionSavePending.value = false;
  await selector.refresh();
  assert.equal(selector.statusFor(source), 'ready');
  assert.equal(expanded, 1);
});

test('late rule helper replies cannot publish cache entries while Session is pending', async () => {
  const reply = gate();
  const feature = { type: 'CDS', record: 'record', id: 'gene', qualifiers: {} };
  const rule = { feat: '*', qual: 'product', val: 'gene' };
  const owner = createRulePreparation({
    state: { extractedFeatures: { value: [feature] }, manualSpecificRules: [rule], sessionOperationAvailability },
    evaluate: () => reply.promise
  });
  const pending = owner.prepare([rule]);
  state.sessionSavePending.value = true;
  reply.release({ matches: [[0]], priorities: [0] });
  assert.equal(await pending, false);
  assert.equal(ruleMatchesFeature(feature, rule), null);
  state.sessionSavePending.value = false;
});
