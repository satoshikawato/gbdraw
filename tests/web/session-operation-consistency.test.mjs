import assert from 'node:assert/strict';
import test from 'node:test';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { withDrawings } from './helpers/drawing-state.mjs';

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
const { createRulePreparation } = await load('app/rule-matching.js');
const { ruleMatchesFeature } = await load('services/rule-matchers.js');
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
  // PD-OI-051 (D-37): draft edits stay available while Generate runs. Undo and
  // Redo use History availability instead (D-28, next test).
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

const trackedRefs = () => {
  let reads = null;
  const makeRef = (initial) => {
    let current = initial;
    const box = {
      writes: 0,
      get value() { reads?.add(box); return current; },
      set value(next) { current = next; box.writes += 1; }
    };
    return box;
  };
  const dependencies = (read) => {
    reads = new Set();
    try { read(); return [...reads]; } finally { reads = null; }
  };
  return { makeRef, dependencies };
};
const writesOf = refs => refs.reduce((total, box) => total + box.writes, 0);
const generating = { status: 'busy', reason: 'Generating diagram. Retry after generation finishes.' };
const applying = { status: 'busy', reason: 'Applying an edit. Retry after the edit finishes.' };

test('D-28: Undo and Redo are busy while Generate replaces the artifact; draft edits stay allowed', async (t) => {
  let value = 1;
  let fingerprint = 0;
  const refs = trackedRefs();
  const history = historyFor({
    makeRef: refs.makeRef,
    buildIntent: async () => ({ value }),
    applyIntent: async (intent) => { value = intent.value; },
    captureGeneratedArtifactHandle: () => ({
      retainedBytes: 1, identity: { fingerprint: String(fingerprint), compactSignature: '' }
    }),
    restoreGeneratedArtifactHandle: async () => {}
  });
  await history.initializeIntentBaseline();
  await history.runUndoable('first edit', () => { value = 2; });
  await history.runUndoable('second edit', () => { value = 3; });
  await history.undo();
  assert.deepEqual([history.getUndoCount(), history.getRedoCount(), value], [1, 1, 2]);
  const dependencies = refs.dependencies(() => [history.canUndo(), history.canRedo()]);
  const writesBefore = writesOf(dependencies);

  const render = gate();
  let started = false;
  state.processing.value = true;
  t.after(() => { state.processing.value = false; });
  const generate = history.runUndoableArtifactReplacement('Generate diagram', async () => {
    started = true;
    const result = await render.promise;
    fingerprint += 1;
    return result;
  }, { shouldCommit: result => result?.status === 'ok' });
  while (!started) await Promise.resolve();
  assert.ok(writesOf(dependencies) > writesBefore, 'canUndo/canRedo must be reactive to the open replacement');
  assert.equal(history.canUndo(), false);
  assert.equal(history.canRedo(), false);
  assert.deepEqual(await history.undo(), generating);
  assert.deepEqual(await history.redo(), generating);
  assert.deepEqual([history.getUndoCount(), history.getRedoCount(), value], [1, 1, 2]);
  assert.equal(sessionOperationAvailability(), null);
  await history.runUndoable('draft edit during Generate', () => { value = 4; });
  assert.equal(value, 4);

  const writesWhileOpen = writesOf(dependencies);
  render.release({ status: 'ok' });
  await generate;
  state.processing.value = false;
  assert.ok(writesOf(dependencies) > writesWhileOpen, 'closing the replacement must refresh availability');
  assert.equal(history.getUndoCount(), 3);
  assert.equal(history.canUndo(), true);
  assert.equal(await history.undo(), true);
  assert.equal(await history.undo(), true);
  assert.equal(value, 2);
});

test('D-28: Undo and Redo are busy while a History checkpoint is open', async () => {
  let value = 1;
  const history = historyFor({
    buildIntent: async () => ({ value }),
    applyIntent: async (intent) => { value = intent.value; },
    buildCheckpoint: () => ({ value }),
    applyCheckpoint: async (checkpoint) => { value = checkpoint.value; }
  });
  await history.captureBaseline();
  await history.runUndoable('edit', () => { value = 2; });
  const work = gate();
  let started = false;
  const checkpoint = history.runUndoableCheckpoint('Reset settings', async () => {
    started = true;
    await work.promise;
    value = 5;
  });
  while (!started) await Promise.resolve();
  assert.equal(history.canUndo(), false);
  assert.deepEqual(await history.undo(), applying);
  assert.equal(value, 2);
  work.release();
  await checkpoint;
  assert.equal(history.canUndo(), true);
  assert.equal(await history.undo(), true);
  assert.equal(value, 2);
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

// R13: the composition root passes the busy source of an edit still applying
// as an argument; availability keeps no registered callback on `state`.
test('Save and Load report an active asynchronous History mutation until settlement', async () => {
  const history = historyFor();
  const preparationBusyReason = () => history.mutationPending()
    ? 'Applying an edit. Retry after the edit finishes.' : '';
  const work = gate();
  let started = false;
  const pending = history.runUndoable('asynchronous edit', async () => {
    started = true;
    await work.promise;
  });
  while (!started) await Promise.resolve();
  assert.match(sessionOperationAvailability('save', preparationBusyReason).reason, /Applying/);
  assert.match(sessionOperationAvailability('load', preparationBusyReason).reason, /Applying/);
  // Draft edits and callers without the source do not read it.
  assert.equal(sessionOperationAvailability('mutation', preparationBusyReason), null);
  assert.equal(sessionOperationAvailability('save'), null);
  // A Generate in progress is worded before the edit still applying.
  state.processing.value = true;
  assert.match(sessionOperationAvailability('save', preparationBusyReason).reason, /Generating/);
  state.processing.value = false;
  work.release();
  await pending;
  assert.equal(sessionOperationAvailability('save', preparationBusyReason), null);
  assert.equal('sessionPreparationBusyReason' in state, false);
});

test('Save and Load consult the availability the composition root passes (R13)', async () => {
  const { exportSession, importSession } = await load('services/config.js');
  const calls = [];
  const busyFrom = (operation, call) => (requested) => {
    calls.push(requested);
    return requested === operation && calls.length >= call ? applying : null;
  };
  // Save checks before it starts and again when its turn comes.
  assert.deepEqual(await exportSession(null, { availability: busyFrom('save', 2) }), applying);
  assert.deepEqual(calls, ['save', 'save']);
  assert.equal(state.sessionSavePending.value, false);
  calls.length = 0;
  const input = { files: [{ name: 'session.gbdraw' }], value: 'session.gbdraw' };
  assert.deepEqual(await importSession({ target: input }, { availability: busyFrom('load', 1) }), applying);
  assert.deepEqual(calls, ['load']);
  assert.equal(input.value, '');
  assert.equal(state.sessionImportPending.value, false);
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
    state: withDrawings({ extractedFeatures: { value: [feature] }, manualSpecificRules: [rule], sessionOperationAvailability }),
    evaluate: () => reply.promise
  });
  const pending = owner.prepare([rule]);
  state.sessionSavePending.value = true;
  reply.release({ matches: [[0]], priorities: [0] });
  assert.equal(await pending, false);
  assert.equal(ruleMatchesFeature(feature, rule), null);
  state.sessionSavePending.value = false;
});
