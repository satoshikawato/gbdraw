import assert from 'node:assert/strict';
import test from 'node:test';
import { createRecordDisplayControls, RECORD_TARGET_NOT_DISCOVERED } from '../../gbdraw/web/js/app/record-display-options.js';
import {
  bindRecordDisplayDrafts, buildRecordDisplayRows, parseRecordDisplayStart, reconcileRecordDisplayDrafts, effectiveRecordReverseComplement, migrateLegacyRecordDisplayDrafts, recordDisplaySurface, requestedRecordTransform, readRecordDisplayDrafts, validateAnchorIntent, validateRecordDisplayDrafts
} from '../../gbdraw/web/js/services/record-display-model.js';
import { circularDiscoveryForInput, parseSequenceRecordText } from '../../gbdraw/web/js/services/record-discovery.js';
import { isCurrentFeature } from '../../gbdraw/web/js/services/feature-identity.js';

import { createFeaturePlacementActions } from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { adoptCurrentSessionResources, createCombinedSessionResourceFileView, createSessionResourceFileView } from '../../gbdraw/web/js/services/session-resource-backing.js';
import { readFileBytes } from '../../gbdraw/web/js/services/file-content-cache.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { withDrawings, withModeDrawings } from './helpers/drawing-state.mjs';

const source = {};
const records = [
  { selector: '#1', recordId: 'same', recordLength: 100, detectedTopology: 'circular' },
  { selector: '#2', recordId: 'same', recordLength: 100, detectedTopology: 'linear' }
];
const rows = buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-1', source, records });

test('discovery observes topology without inferring it from the diagram or record name', () => {
  const parsed = parseSequenceRecordText('LOCUS       same 100 bp DNA circular\n//\nLOCUS       same 100 bp DNA linear\n//\nLOCUS       circular 100 bp DNA\n//\n', 'genbank');
  assert.deepEqual(parsed.map(r => r.detectedTopology), ['circular', 'linear', 'unknown']);
  assert.equal(parseSequenceRecordText('>circular\nACGT\n', 'fasta')[0].detectedTopology, 'unknown');
});

test('ALL rows and exact-one use instance selectors, with no duplicate-ID fallback', () => {
  assert.equal(rows.length, 2);
  assert.notEqual(rows[0].key, rows[1].key);
  assert.deepEqual(buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-1', source, records, selector: '#2' }), [rows[1]]);
  for (const selector of ['same', 'missing', '#3']) {
    assert.throws(() => buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-1', records, selector }), /selector/);
  }
});

test('source replacement purges drafts; inactive selectors survive rediscovery and reorder', () => {
  const drafts = rows.map(r => ({ scope: r.scope, sourceUid: r.sourceUid, selector: r.selector, topologyOverride: null, startCoordinate: 25 }));
  assert.deepEqual(reconcileRecordDisplayDrafts(drafts, [...rows].reverse()), drafts);
  assert.deepEqual(reconcileRecordDisplayDrafts(drafts, rows, ['card-1']), []);
  assert.deepEqual(reconcileRecordDisplayDrafts(drafts, rows.map(r => ({ ...r, sourceUid: 'replacement' }))), []);
});

// OV-278: Load projects a request's record-ID selection into a draft with the
// placeholder selector '#1'; it belongs to the row with that record ID.
test('a draft that names its record by ID moves to that record once the rows are known', () => {
  const shared = buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card-2', source, records: [
    { selector: '#1', recordId: 'AAA.1', recordLength: 60, detectedTopology: 'circular' },
    { selector: '#2', recordId: 'BBB.1', recordLength: 60, detectedTopology: 'circular' }
  ] });
  const draft = (selector, recordId, startCoordinate) => ({ sourceUid: 'card-2', selector, recordId,
    topologyOverride: null, startCoordinate, reverseComplementOverride: null, anchorIntent: null });
  assert.deepEqual(bindRecordDisplayDrafts([draft('#1', 'BBB.1', 11)], shared), [draft('#2', 'BBB.1', 11)]);
  // Unchanged (the same array): a draft on its record, records not yet listed,
  // a duplicate ID, a record the file lacks, and a row that already has a draft.
  for (const [drafts, known] of [
    [[draft('#1', 'AAA.1', 7)], shared],
    [[draft('#1', 'BBB.1', 11)], []],
    [[draft('#2', 'same', 25)], rows.map((row) => ({ ...row, sourceUid: 'card-2' }))],
    [[draft('#1', 'CCC.1', 11)], shared],
    [[draft('#1', 'BBB.1', 11), draft('#2', 'BBB.1', 9)], shared]
  ]) assert.equal(bindRecordDisplayDrafts(drafts, known), drafts);
});

test('blank remains null, explicit 1 survives, malformed input is never clamped', () => {
  for (const value of [null, '', '  ']) assert.equal(parseRecordDisplayStart(value), null);
  for (const value of [1, 100, '25']) assert.equal(parseRecordDisplayStart(value), Number(value));
  for (const value of [true, false, undefined, 0, -1, 1.5, 'oops', Infinity, [], {}]) {
    assert.throws(() => parseRecordDisplayStart(value), /integer/);
  }
});

const anchorIntent = {
  schema: 1,
  recordKey: 'record-key',
  biologicalFeatureId: 'feature-key',
  placement: 'anchor',
  anchor: 'five-prime',
  offsetBp: 0,
  orientForward: true
};

test('record display drafts use one exact transform and provenance contract', () => {
  // A drawing's draft: its mode is the drawing's (PR-1).
  const draft = {
    sourceUid: 'card-1', selector: '#1', recordId: 'same',
    topologyOverride: null, startCoordinate: 25, reverseComplementOverride: true,
    anchorIntent
  };
  const drafts = [draft];
  assert.equal(validateRecordDisplayDrafts(drafts), drafts);
  assert.equal(validateAnchorIntent(anchorIntent), anchorIntent);
  // A Session 41-44 draft names its mode until the Session 46 split.
  const migrated = migrateLegacyRecordDisplayDrafts([{
    scope: 'linear', sourceUid: 'card-1', selector: '#1', recordId: 'same',
    topologyOverride: null, startCoordinate: 25
  }]);
  assert.deepEqual(migrated[0], { scope: 'linear', ...draft, reverseComplementOverride: null, anchorIntent: null });
  assert.equal(validateRecordDisplayDrafts(migrated, { scoped: true }), migrated);
  assert.throws(() => validateRecordDisplayDrafts(migrated));
  for (const invalid of [
    { ...draft, extra: true },
    { ...draft, scope: 'linear' },
    { ...draft, reverseComplementOverride: 'true' },
    { ...draft, anchorIntent: { ...anchorIntent, schema: 2 } }
  ]) assert.throws(() => validateRecordDisplayDrafts([invalid]), JSON.stringify(invalid));
  // OV-575, OV-639: a row's selector is text, as gbdraw/session_io.py:853-857
  // requires, in every Session version: Save refuses another row, and Load
  // reads the others (the type is checked before the `#N` pattern).
  for (const selector of [['#1'], 1]) {
    assert.throws(() => validateRecordDisplayDrafts([{ ...draft, selector }]), /^Error: INPUT_INVALID\/FIELDS$/, JSON.stringify(selector));
  }
  assert.throws(() => validateRecordDisplayDrafts({}), /^Error: INPUT_INVALID\/FIELDS$/);
  const read = readRecordDisplayDrafts([{ ...draft, selector: ['#1'] }, draft, draft, 'x']);
  assert.deepEqual(read, { value: [draft], lost: true, dropped: 3 });
  assert.deepEqual(readRecordDisplayDrafts({ 0: draft }), { value: [], lost: true });
  assert.deepEqual(readRecordDisplayDrafts(null), { value: [], lost: false });
  assert.deepEqual(readRecordDisplayDrafts(migrated, { scoped: true }), { value: migrated, lost: false });
});

test('complete orientation override has one owner while cropped reverse stays region-owned', () => {
  const row = { ...rows[0], reverse: false, cropped: false };
  assert.equal(effectiveRecordReverseComplement(row, { reverseComplementOverride: true }), true);
  assert.equal(effectiveRecordReverseComplement(
    { ...row, reverse: true, cropped: true },
    { reverseComplementOverride: false }
  ), true);
  assert.deepEqual(requestedRecordTransform(row, {
    topologyOverride: null, startCoordinate: 25, reverseComplementOverride: true
  }), {
    display: { isCircular: null, startCoordinate: 25 },
    reverseComplement: true
  });
});

test('detected topology, nullable reset, crop lock, and RC default stay separate from raw draft', () => {
  const draft = { topologyOverride: false, startCoordinate: 25 };
  assert.equal(recordDisplaySurface(rows[0], draft).startEnabled, false);
  assert.equal(draft.startCoordinate, 25);
  assert.equal(recordDisplaySurface(rows[0], { ...draft, topologyOverride: null }).startEnabled, true);
  assert.equal(recordDisplaySurface(rows[1]).startEnabled, false);
  assert.equal(recordDisplaySurface(rows[1], { topologyOverride: true }).startEnabled, true);
  assert.equal(recordDisplaySurface(rows[0], {}, { cropped: true }).startEnabled, false);
  assert.equal(recordDisplaySurface(rows[0], {}, { reverse: true }).currentStart, 100);
});

const descriptor = (name, text) => ({ kind: 'genbank', name, type: 'text/plain',
  lastModified: 0, size: Buffer.byteLength(text), encoding: 'base64', data: btoa(text) });
const compositeControls = ({ linear = false, discovered = true, createHistory = null } = {}) => {
  const resources = Object.fromEntries(['first', 'second'].map((id) => [id,
    descriptor(`${id}.gb`, `LOCUS       ${id} 100 bp DNA circular\n//\n`)]));
  const table = adoptCurrentSessionResources(resources);
  const components = Object.keys(resources).map((resourceId) => ({ resourceId }));
  const makeFile = () => createCombinedSessionResourceFileView(table, components);
  const file = makeFile();
  let committed = { resources, renderRequest: { mode: linear ? 'linear' : 'circular', records: components.map(({ resourceId }, i) => ({
    recordKey: `record-${i + 1}`, cardinality: 'exactly_one', source: { kind: 'genbank', resourceId },
    selector: null })) } };
  const state = { sessionOperationAvailability: () => null, mode: { value: linear ? 'linear' : 'circular' }, cInputType: { value: 'gb' },
    lInputType: { value: 'gb' }, files: { c_gb: file }, linearSeqs: [],
    form: createDefaultForm(), adv: createDefaultAdv('circular'),
    recordDisplayDrafts: [], featurePlacementOverrides: {},
    circularRecordList: { value: discovered ? components.map(({ resourceId }, i) => ({
      record_id: resourceId, record_length: 100, selector: `#${i + 1}`, detectedTopology: 'circular' })) : [] },
    // A loaded Session has read no record bytes, so Circular discovery is not current.
    circularRecordDiscovery: discovered
      ? { status: 'ready', error: '', inputType: 'gb', primaryFile: file, pairedFile: null }
      : { status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null },
    featureCatalog: { value: { items: [{ recordKeys: ['record-1', 'record-2'] }] } } };
  state.form.track_type = 'middle';
  if (linear) state.linearSeqs = components.map(({resourceId}, index) => ({
    uid:`record-${index+1}`,gb:createSessionResourceFileView(table,resourceId),region_record_id:'',region_reverse:false
  }));
  const getCommittedRequest = () => committed.renderRequest;
  const history = createHistory ? createHistory(state) : { runUndoable: (_label, fn) => fn() };
  const linearDiscovery = { status: discovered ? 'ready' : 'loading', error: '' };
  const controls = createRecordDisplayControls({ state, computed: (fn) => ({ get value() { return fn(); } }),
    watch: () => {}, linearRecordsFor: seq => discovered ? [{selector:'#1',recordId:components[Number(seq.uid.slice(-1))-1].resourceId,recordLength:100,detectedTopology:'circular'}] : [],
    linearRecordStatusFor: () => linearDiscovery.status, linearRecordErrorFor: () => linearDiscovery.error,
    runUndoable: history.runUndoable, getCommittedRequest, getCommittedSession: () => committed });
  // The composition root's port: the pure check over the record display's binding.
  const actions = createFeaturePlacementActions({ state: withDrawings(state), runUndoable: history.runUndoable, getCommittedRequest,
    isCurrentFeature: (feature) => isCurrentFeature(feature, controls.sourceBinding()) });
  const feature = { scope: state.mode.value, record_key: 'record-2', biological_feature_id: 'logical-feature' };
  return { state, actions, controls, feature, file, makeFile, linearDiscovery, history,
    enabled: () => actions.choices([feature]).filter((choice) => choice.enabled).map((choice) => choice.value),
    // The committed request with every record naming one saved resource id.
    recommitSources: (resourceId) => {
      committed = { ...committed, renderRequest: { ...committed.renderRequest,
        records: committed.renderRequest.records.map((record) => ({ ...record, source: { kind: 'genbank', resourceId } })) } };
    },
    commitCombined: async () => {
      const bytes = await readFileBytes(file);
      const combined = { ...descriptor(file.name, ''), size: bytes.byteLength,
        data: Buffer.from(bytes).toString('base64') };
      committed = { resources: { combined }, renderRequest: { ...committed.renderRequest,
        records: committed.renderRequest.records.map((record, i) => ({ ...record,
          source: { kind: 'genbank', resourceId: 'combined' }, selector: { kind: 'recordIndex', index: i } })) } };
    } };
};

test('same composite retains semantic capability when the committed resource becomes combined', async () => {
  const model = compositeControls();
  const all = ['auto', 'main', 'outward', 'inward'];
  assert.deepEqual(model.enabled(), all);
  model.actions.setPlacement([model.feature], 'main');
  assert.equal(model.actions.valueFor(model.feature), 'main');
  await model.commitCombined();
  assert.deepEqual(model.enabled(), all);
  assert.equal(model.actions.valueFor(model.feature), 'main');
});

test('manual start and orientation edits clear feature provenance', () => {
  const model = compositeControls();
  const row = model.controls.allRows.value[0];
  model.controls.setResolvedTransform(row, {
    startCoordinate: 25, reverseComplement: true, anchorIntent
  });
  assert.deepEqual(model.state.recordDisplayDrafts[0].anchorIntent, anchorIntent);
  model.controls.setStart(row, 10);
  assert.equal(model.state.recordDisplayDrafts[0].anchorIntent, null);
  model.controls.setResolvedTransform(row, {
    startCoordinate: 25, reverseComplement: true, anchorIntent
  });
  model.controls.setReverseComplement(row, false);
  assert.equal(model.state.recordDisplayDrafts[0].anchorIntent, null);
});

test('explicit popup target and its draft checkpoint stay record-bound', () => {
  const model = compositeControls();
  const { row, target } = model.controls.targetForFeature(model.feature);
  assert.equal(target.recordKey, 'record-2');
  assert.equal(target.canonicalRecordKey, 'record-2');
  assert.equal(target.source.resourceId, 'second');
  assert.equal(target.recordLength, 100);
  assert.equal(target.effectiveCircular, true);
  assert.equal(target.committedReverseComplement, false);
  assert.equal(target.members.length, 1);
  const before = model.controls.captureTargetDraft(row);
  model.controls.commitResolvedTransform(row, {
    startCoordinate: 25,
    reverseComplement: true,
    anchorIntent
  });
  assert.equal(model.state.recordDisplayDrafts.length, 1);
  model.controls.restoreTargetDraft(before);
  assert.deepEqual(model.state.recordDisplayDrafts, []);
  assert.throws(
    () => model.controls.targetForFeature({ ...model.feature, record_key: 'missing' }),
    /stale or ambiguous/
  );
});

test('a popup target whose records are not read yet is not stale (Session Load, 776a2f93)', () => {
  for (const linear of [false, true]) {
    const model = compositeControls({ linear, discovered: false });
    assert.throws(() => model.controls.targetForFeature(model.feature), (error) => {
      assert.equal(error.kind, RECORD_TARGET_NOT_DISCOVERED);
      assert.doesNotMatch(error.message, /stale|ambiguous/);
      return true;
    });
  }
  const failed = compositeControls({ linear: true, discovered: false });
  failed.linearDiscovery.status = 'error';
  failed.linearDiscovery.error = 'Records could not be loaded: malformed LOCUS line.';
  assert.throws(() => failed.controls.targetForFeature(failed.feature), (error) => {
    assert.notEqual(error.kind, RECORD_TARGET_NOT_DISCOVERED);
    assert.equal(error.message, 'Records could not be loaded: malformed LOCUS line.');
    return true;
  });
  const discovered = compositeControls({ linear: true });
  assert.throws(
    () => discovered.controls.targetForFeature({ ...discovered.feature, record_key: 'missing' }),
    (error) => error.kind !== RECORD_TARGET_NOT_DISCOVERED && /stale or ambiguous/.test(error.message)
  );
});

test('a Linear rotation writes the File card orientation and leaves no row override (CO-03, N-09)', () => {
  const model = compositeControls({ linear: true });
  const row = model.controls.rows.value.find((entry) => entry.sourceUid === 'record-2');
  const sequence = model.state.linearSeqs[1];
  const before = model.controls.captureTargetDraft(row);
  model.controls.commitResolvedTransform(row, { startCoordinate: 25, reverseComplement: true, anchorIntent });
  assert.equal(sequence.region_reverse, true, 'the card checkbox shows the applied orientation');
  assert.equal(model.state.recordDisplayDrafts[0].reverseComplementOverride, null);
  assert.equal(model.state.recordDisplayDrafts[0].startCoordinate, 25);
  sequence.region_reverse = false;
  assert.equal(model.controls.rows.value.find((entry) => entry.sourceUid === 'record-2').reverse, false,
    'no override can mask a later card checkbox edit');
  sequence.region_reverse = true;
  model.controls.restoreTargetDraft(before);
  assert.equal(sequence.region_reverse, false, 'a failed apply restores the card orientation');
  assert.deepEqual(model.state.recordDisplayDrafts, []);
  model.controls.setReverseComplement(row, true);
  assert.deepEqual([sequence.region_reverse, model.state.recordDisplayDrafts[0].reverseComplementOverride], [true, null]);
});

// The real History manager over the intent the app captures for these owners:
// the record display drafts and the Linear File card orientation.
const draftHistory = (state) => {
  const capture = () => ({
    recordDisplayDrafts: structuredClone(state.recordDisplayDrafts),
    linearReverse: state.linearSeqs.map((seq) => Boolean(seq.region_reverse))
  });
  const restore = (intent) => {
    state.recordDisplayDrafts.splice(0, state.recordDisplayDrafts.length,
      ...structuredClone(intent.recordDisplayDrafts));
    intent.linearReverse.forEach((value, index) => { state.linearSeqs[index].region_reverse = value; });
  };
  return createHistoryManager({
    buildIntent: capture, applyIntent: restore, buildCheckpoint: capture, applyCheckpoint: restore
  });
};

test('Apply on Generate is one undoable draft step, and re-staging replaces the draft (PD-OI-085)', async () => {
  for (const linear of [false, true]) {
    const model = compositeControls({ linear, createHistory: draftHistory });
    const pending = () => model.controls.targetForFeature(model.feature).target.pendingTransform;
    const reverse = () => (linear ? model.state.linearSeqs[1].region_reverse : null);
    const { row } = model.controls.targetForFeature(model.feature);
    assert.equal(pending(), null);

    await model.controls.setResolvedTransform(row, {
      startCoordinate: 25, reverseComplement: true, anchorIntent
    });
    assert.equal(model.history.getUndoCount(), 1, 'one History step');
    assert.equal(model.history.undoLabel(), 'Rotate record to feature on Generate');
    assert.deepEqual(pending(), { startCoordinate: 25, reverseComplement: true });
    assert.equal(reverse(), linear ? true : null, 'a single-record Linear card owns the orientation');

    const forwardIntent = { ...anchorIntent, orientForward: false };
    await model.controls.setResolvedTransform(row, {
      startCoordinate: 40, reverseComplement: false, anchorIntent: forwardIntent
    });
    assert.equal(model.state.recordDisplayDrafts.length, 1, 'the later value replaces the draft');
    assert.deepEqual(model.state.recordDisplayDrafts[0].anchorIntent, forwardIntent);
    assert.deepEqual(pending(), { startCoordinate: 40, reverseComplement: false });
    assert.equal(model.history.getUndoCount(), 2);

    await model.history.undo();
    assert.deepEqual(pending(), { startCoordinate: 25, reverseComplement: true });
    assert.deepEqual(model.state.recordDisplayDrafts[0].anchorIntent, anchorIntent);
    assert.equal(reverse(), linear ? true : null);
    await model.history.undo();
    assert.deepEqual(model.state.recordDisplayDrafts, []);
    assert.equal(pending(), null);
    assert.equal(reverse(), linear ? false : null);
    await model.history.redo();
    assert.deepEqual(pending(), { startCoordinate: 25, reverseComplement: true });
    assert.equal(reverse(), linear ? true : null);
  }
});

test('same-name, same-content source replacement does not retain old feature capability', () => {
  const model = compositeControls();
  assert.equal(model.enabled().length, 4);
  model.state.files.c_gb = model.makeFile();
  assert.equal(model.state.files.c_gb.name, model.file.name);
  assert.notEqual(model.state.files.c_gb, model.file);
  assert.deepEqual(model.enabled(), []);
  assert.throws(() => model.actions.setPlacement([model.feature], 'main'), /A selected feature is not in the displayed Result\./);
});

test('isCurrentFeature is a pure check of the committed record, its bound source, and the saved resource', () => {
  const genbank = {};
  const fasta = {};
  const gff = {};
  const request = { mode: 'linear', records: [
    { recordKey: 'seq-a', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'a' }, selector: null },
    { recordKey: 'seq-b', cardinality: 'all', source: { kind: 'gffFasta', gffResourceId: 'b-gff', fastaResourceId: 'b-fasta' },
      selector: null },
    { recordKey: 'seq-c:2', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'c' },
      selector: { kind: 'recordIndex', index: 1 } }
  ] };
  const sources = [
    { scope: 'linear', sourceUid: 'seq-a', source: genbank, paired: null },
    { scope: 'linear', sourceUid: 'seq-b', source: gff, paired: fasta },
    { scope: 'linear', sourceUid: 'seq-c', source: genbank, paired: null }
  ];
  const saved = new Map([['a', genbank], ['b-gff', gff], ['b-fasta', fasta], ['c', genbank]]);
  const binding = { request, sources, boundSources: sources.map((entry) => ({ ...entry })),
    matchesSavedSource: (file, resourceId) => saved.get(resourceId) === file };
  const current = (recordKey, changes = {}) => isCurrentFeature({ record_key: recordKey }, { ...binding, ...changes });
  assert.equal(current('seq-a'), true);
  assert.equal(current('seq-b:3'), true, 'one record of a source drawn whole');
  assert.equal(current('seq-c:2'), true, 'a record named by its index');
  assert.equal(current('seq-b:0'), false);
  assert.equal(current('seq-missing'), false);
  assert.equal(current('seq-a', { request: null }), false);
  assert.equal(current('seq-a', { request: { ...request, mode: 'circular' } }), false);
  // A source replaced after the request was committed is not the bound one.
  assert.equal(current('seq-a', { sources: [{ ...sources[0], source: {} }, ...sources.slice(1)] }), false);
  assert.equal(current('seq-a', { boundSources: binding.boundSources.slice(1) }), false);
  // Both the primary input and its paired FASTA match the committed Session.
  assert.equal(current('seq-a', { matchesSavedSource: () => false }), false);
  assert.equal(current('seq-b:1', { matchesSavedSource: (file, resourceId) => resourceId !== 'b-fasta' }), false);
});

// A feature row of a Session before 40 is stored unchecked: it can lack `record_key` or hold a value that is
// not text. Such a key names no committed record, also when the request draws a source whole.
test('a feature whose record key is not text is not current, with a source drawn whole (OV-618)', () => {
  const file = {};
  const request = { mode: 'linear', records: [
    { recordKey: 'card', cardinality: 'all', source: { kind: 'genbank', resourceId: 'a' }, selector: null }
  ] };
  const sources = [{ scope: 'linear', sourceUid: 'card', source: file, paired: null }];
  const binding = { request, sources, boundSources: sources.map((entry) => ({ ...entry })), matchesSavedSource: () => true };
  assert.equal(isCurrentFeature({ record_key: 'card:1' }, binding), true);
  for (const feature of [{}, { record_key: null }, { recordKey: 'card:1' }, { record_key: 5 }, { record_key: ['card:1'] },
    { record_key: new String('card:1') }]) {
    assert.equal(isCurrentFeature(feature, binding), false, JSON.stringify(feature));
  }
  // The placement control reads a feature's identity from either spelling; the check reads `record_key`.
  const state = { sessionOperationAvailability: () => null, mode: { value: 'linear' }, form: createDefaultForm(), adv: createDefaultAdv('linear'),
    featurePlacementOverrides: {}, featureCatalog: { value: { items: [{ recordKeys: ['card:1'] }] } } };
  const actions = createFeaturePlacementActions({ state: withDrawings(state), runUndoable: (_label, fn) => fn(),
    getCommittedRequest: () => request, isCurrentFeature: (feature) => isCurrentFeature(feature, binding) });
  const enabled = (feature) => actions.choices([feature]).filter((choice) => choice.enabled).map((choice) => choice.value);
  assert.ok(enabled({ record_key: 'card:1', biological_feature_id: 'gene' }).includes('auto'));
  assert.deepEqual(enabled({ recordKey: 'card:1', biologicalFeatureId: 'gene' }), []);
});

test('the record display binding is refreshed for the committed request before the check reads it', async () => {
  const model = compositeControls();
  const before = model.controls.sourceBinding();
  assert.equal(before.boundSources[0].source, model.file);
  assert.equal(isCurrentFeature(model.feature, before), true);
  model.state.files.c_gb = model.makeFile();
  assert.equal(isCurrentFeature(model.feature, model.controls.sourceBinding()), false);
  // A new committed request binds the inputs it was drawn from.
  await model.commitCombined();
  const after = model.controls.sourceBinding();
  assert.notEqual(after.request, before.request);
  assert.equal(after.boundSources[0].source, model.state.files.c_gb);
  assert.equal(isCurrentFeature(model.feature, after), true);
});

// OV-557: a resource id is Session data. A record that names an id the Session
// does not hold has no saved source, also when the id names an inherited member.
test('a record that names a missing resource id matches no saved source', () => {
  const model = compositeControls();
  const { matchesSavedSource } = model.controls.sourceBinding();
  // The combined view of two resources is no saved resource by itself.
  for (const missing of ['missing', 'constructor', 'toString', '__proto__', undefined]) {
    assert.equal(matchesSavedSource(model.file, missing), false, String(missing));
  }
  assert.equal(isCurrentFeature(model.feature, model.controls.sourceBinding()), true);
  model.recommitSources('constructor');
  assert.equal(isCurrentFeature(model.feature, model.controls.sourceBinding()), false,
    'the feature of a record whose saved source is missing is not current');
  assert.deepEqual(model.enabled(), []);
});

test('unsupported draft lanes remain unavailable independently of current placement value', () => {
  const model = compositeControls();
  model.actions.setPlacement([model.feature], 'outward');
  model.state.form.track_type = 'tuckin';
  assert.equal(model.actions.valueFor(model.feature), 'outward');
  assert.deepEqual(model.enabled(), ['auto', 'main']);
  assert.throws(() => model.actions.setPlacement([model.feature], 'inward'), /The current track settings have no feature track for this placement\./);
});


test('alignment intent changes only its target direction and preserves pending display settings', () => {
  const model=compositeControls({linear:true});
  const rows=model.controls.allRows.value.filter(row=>row.scope==='linear');
  model.controls.setStart(rows[0],25);model.controls.setStart(rows[1],40);
  const pending=structuredClone(model.state.recordDisplayDrafts);
  const directions=[{recordKey:'record-1',reverseComplement:true}];
  const checkpoint=model.controls.captureAlignmentOrientationIntent(directions);
  model.controls.commitAlignmentOrientations(directions);
  assert.equal(model.state.linearSeqs[0].region_reverse,true);
  assert.equal(model.state.linearSeqs[1].region_reverse,false);
  assert.equal(model.state.recordDisplayDrafts.length,pending.length);
  assert.equal(model.state.recordDisplayDrafts[0].startCoordinate,25);
  assert.deepEqual(model.state.recordDisplayDrafts[1],pending[1]);
  assert.equal(effectiveRecordReverseComplement({...rows[0],reverse:true},model.state.recordDisplayDrafts[0]),true);
  // A later normal Reverse checkbox remains effective, with no forced override.
  model.state.linearSeqs[0].region_reverse=false;
  assert.equal(effectiveRecordReverseComplement({...rows[0],reverse:false},model.state.recordDisplayDrafts[0]),false);
  model.controls.restoreAlignmentOrientationIntent(checkpoint);
  assert.deepEqual(model.state.recordDisplayDrafts,pending);
});

// OV-579: each mode has its own drawing, and a draft names no mode. The
// alignment checkpoint is of the Linear drawing's drafts.
test('an alignment checkpoint restores the Linear drafts and leaves the Circular drawing as it is', () => {
  const model = compositeControls({ linear: true });
  const circular = { form: createDefaultForm(), adv: createDefaultAdv('circular'), featurePlacementOverrides: {},
    // Not yet on its record's row (OV-278): only an edit of the Circular drawing moves it.
    recordDisplayDrafts: [{ sourceUid: 'circular', selector: '#1', recordId: 'second', topologyOverride: null,
      startCoordinate: 7, reverseComplementOverride: null, anchorIntent: null }] };
  withModeDrawings(model.state, { circular, linear: model.state });
  const circularBefore = structuredClone(circular.recordDisplayDrafts);
  const rows = model.controls.allRows.value.filter((row) => row.scope === 'linear');
  model.controls.setStart(rows[0], 25);
  // A row override, as a card that drew several records wrote it.
  model.state.recordDisplayDrafts[0].reverseComplementOverride = true;
  const pending = structuredClone(model.state.recordDisplayDrafts);
  const directions = [{ recordKey: 'record-1', reverseComplement: true }];
  const checkpoint = model.controls.captureAlignmentOrientationIntent(directions);
  model.controls.commitAlignmentOrientations(directions);
  assert.equal(model.state.linearSeqs[0].region_reverse, true);
  assert.equal(model.state.recordDisplayDrafts[0].reverseComplementOverride, null);
  model.controls.restoreAlignmentOrientationIntent(checkpoint);
  assert.equal(model.state.linearSeqs[0].region_reverse, false);
  assert.deepEqual({ linear: model.state.recordDisplayDrafts, circular: circular.recordDisplayDrafts },
    { linear: pending, circular: circularBefore });
});


test('fresh loaded alignment binds canonical sources before discovery and restores pending orientation', () => {
  const model = compositeControls({linear:true, discovered:false});
  model.state.linearSeqs[0].region_reverse = true;
  const target = {recordKey:'record-1', reverseComplement:true};
  assert.equal(model.controls.allRows.value.filter(row => row.scope === 'linear').length, 0);
  const checkpoint = model.controls.captureAlignmentOrientationIntent([target]);
  model.controls.commitAlignmentOrientations([target]);
  assert.equal(model.state.linearSeqs[0].region_reverse, true);
  assert.equal(model.state.linearSeqs[1].region_reverse, false);
  // Artifact rollback restores its canonical direction first; the same owner
  // then restores the pre-operation pending input, without needing discovery.
  model.state.linearSeqs[0].region_reverse = false;
  model.controls.restoreAlignmentOrientationIntent(checkpoint);
  assert.equal(model.state.linearSeqs[0].region_reverse, true);
  model.state.linearSeqs[0].gb = new Blob(['replacement']);
  assert.throws(() => model.controls.captureAlignmentOrientationIntent([target]), /source binding changed/);
});

test('Circular discovery metadata belongs only to the current source and complete input pair', () => {
  const source = {};
  const state = { sessionOperationAvailability: () => null, cInputType: { value: 'gb' }, files: { c_gb: source },
    circularRecordList: { value: [{ selector: '#1', record_id: 'committed' }] },
    circularRecordDiscovery: { status: 'ready', inputType: 'gb', primaryFile: source, pairedFile: null, error: '' } };
  assert.equal(circularDiscoveryForInput(state).status, 'ready');
  assert.equal(circularDiscoveryForInput(state).records.length, 1);
  state.files.c_gb = {};
  assert.equal(circularDiscoveryForInput(state).status, 'deferred');
  assert.deepEqual(circularDiscoveryForInput(state).records, []);
  state.files.c_gb = source;
  state.circularRecordDiscovery.status = 'deferred';
  assert.equal(circularDiscoveryForInput(state).status, 'deferred');
  assert.deepEqual(circularDiscoveryForInput(state).records, []);
  state.cInputType.value = 'gff';
  state.files.c_gff = {};
  assert.equal(circularDiscoveryForInput(state).status, 'idle');
  state.files.c_fasta = {};
  assert.equal(circularDiscoveryForInput(state).status, 'deferred');
  assert.deepEqual(circularDiscoveryForInput(state).records, []);
});

// OV-278: Load of a Session without a Web draft writes a draft selected by
// record ID on '#1'. The row of its record reads it, and an edit there leaves
// one draft for that record.
test('a loaded draft selected by record ID is read and edited on its own record row', () => {
  const model = compositeControls();
  const loaded = (recordId, startCoordinate) => ({ sourceUid: 'circular', selector: '#1', recordId, startCoordinate,
    topologyOverride: null, reverseComplementOverride: null, anchorIntent: null });
  model.state.recordDisplayDrafts.push(loaded('first', 7), loaded('second', 11));
  const [first, second] = model.controls.allRows.value;
  assert.equal(model.controls.draftFor(first).startCoordinate, 7);
  assert.equal(model.controls.draftFor(second).startCoordinate, 11);
  model.controls.setStart(second, 12);
  assert.deepEqual(model.state.recordDisplayDrafts.map(({ selector, recordId, startCoordinate }) => [selector, recordId, startCoordinate]),
    [['#1', 'first', 7], ['#2', 'second', 12]]);
});

// OV-579, review U1: a draft that names its record by ID and is not yet on its row (as a Load projects a
// request that selects by record ID, OV-278) is settled before the checkpoint looks for it.
test('an alignment checkpoint holds a draft that its own capture moves onto its row (OV-579)', () => {
  const resources = { both: descriptor('both.gb', 'LOCUS       first 100 bp DNA circular\n//\nLOCUS       second 100 bp DNA circular\n//\n') };
  const table = adoptCurrentSessionResources(resources);
  const committed = { resources, renderRequest: { mode: 'linear', records: [{ recordKey: 'record-1', cardinality: 'exactly_one',
    source: { kind: 'genbank', resourceId: 'both' }, selector: { kind: 'recordId', value: 'second' } }] } };
  const state = { sessionOperationAvailability: () => null, mode: { value: 'linear' }, cInputType: { value: 'gb' }, lInputType: { value: 'gb' }, files: { c_gb: null },
    linearSeqs: [{ uid: 'record-1', gb: createSessionResourceFileView(table, 'both'), region_record_id: 'second', region_reverse: false }],
    form: createDefaultForm(), adv: createDefaultAdv('linear'),
    recordDisplayDrafts: [{ sourceUid: 'record-1', selector: '#1', recordId: 'second', topologyOverride: null,
      startCoordinate: 25, reverseComplementOverride: true, anchorIntent: null }],
    featurePlacementOverrides: {}, circularRecordList: { value: [] },
    circularRecordDiscovery: { status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null },
    featureCatalog: { value: { items: [] } } };
  const controls = createRecordDisplayControls({ state: withDrawings(state), computed: (fn) => ({ get value() { return fn(); } }),
    watch: () => {}, linearRecordsFor: () => [
      { selector: '#1', recordId: 'first', recordLength: 100, detectedTopology: 'circular' },
      { selector: '#2', recordId: 'second', recordLength: 100, detectedTopology: 'circular' }],
    linearRecordStatusFor: () => 'ready', linearRecordErrorFor: () => '',
    runUndoable: (_label, fn) => fn(), getCommittedRequest: () => committed.renderRequest, getCommittedSession: () => committed });
  const directions = [{ recordKey: 'record-1', reverseComplement: true }];
  const checkpoint = controls.captureAlignmentOrientationIntent(directions);
  assert.equal(state.recordDisplayDrafts[0].selector, '#2', 'the capture settled the draft onto its row');
  controls.commitAlignmentOrientations(directions);
  assert.equal(state.recordDisplayDrafts[0].reverseComplementOverride, null);
  controls.restoreAlignmentOrientationIntent(checkpoint);
  assert.equal(state.recordDisplayDrafts[0].reverseComplementOverride, true, 'the row override comes back');
  assert.equal(state.linearSeqs[0].region_reverse, false);
});

// A Circular record is cropped only where Generate crops it: one record offered
// alone, with a region Generate can draw. Emptied region fields name no region,
// a region Generate refuses crops nothing, and separate diagrams or a
// Multi-Record Canvas draw whole records. Display start is unavailable exactly
// for a cropped record (OV-1093).
test('a Circular row is cropped, and loses Display start, where Generate crops its record', () => {
  const model = compositeControls();
  const circularRows = () => model.controls.rows.value.filter((row) => row.sourceUid === 'circular');
  const cropped = () => circularRows().map((row) => row.cropped);
  const startEnabled = () => circularRows().map((row) => model.controls.surfaceFor(row).startEnabled);
  const { form } = model.state.activeDrawing();
  form.multi_record_canvas = false;
  form.circular_region_start = 1;
  form.circular_region_end = 50;
  assert.deepEqual(cropped(), [false, false], 'two records offered: separate diagrams draw whole records');
  assert.deepEqual(startEnabled(), [true, true], 'Display start counts on the whole records');
  form.circular_record_selector = 'first';
  assert.deepEqual(cropped(), [true], 'the record chosen alone is cropped');
  assert.deepEqual(startEnabled(), [false], 'Display start is unavailable for the cropped record');
  form.circular_region_end = '';
  assert.deepEqual(cropped(), [false], 'Start without End crops nothing: Generate refuses the region');
  form.circular_region_end = 200;
  assert.deepEqual(cropped(), [true], 'an End beyond the 100 bp record is clamped: Generate crops to the record end');
  assert.deepEqual(startEnabled(), [false]);
  form.circular_region_start = 150;
  assert.deepEqual(cropped(), [false], 'a Start beyond the 100 bp record crops nothing: Generate refuses the region');
  assert.deepEqual(startEnabled(), [true]);
  form.circular_region_start = 1;
  form.circular_region_end = 50;
  form.multi_record_canvas = true;
  assert.deepEqual(cropped(), [false], 'a Multi-Record Canvas draws whole records');
  form.multi_record_canvas = false;
  form.circular_region_start = '';
  form.circular_region_end = '';
  assert.deepEqual(cropped(), [false], 'emptied fields name no region');
});

// A2 W5 (R10): the drafts of a replaced source go in the transition that
// replaced it (`changeRecordSources`), not in a watcher of the sources.
const sourceControls = () => {
  /** @type {Function[]} */
  const callbacks = [];
  const draft = (sourceUid) => ({ sourceUid, selector: '#1', recordId: sourceUid, topologyOverride: 'linear',
    startCoordinate: null, reverseComplementOverride: null, anchorIntent: null });
  const circular = { form: createDefaultForm(), adv: createDefaultAdv('circular'), recordDisplayDrafts: [draft('circular')] };
  const linear = { form: createDefaultForm(), adv: createDefaultAdv('linear'), recordDisplayDrafts: [draft('record-1'), draft('record-2')] };
  const state = withModeDrawings({ sessionOperationAvailability: () => null, mode: { value: 'linear' }, cInputType: { value: 'gb' }, lInputType: { value: 'gb' },
    files: { c_gb: new Blob(['c']) }, circularRecordList: { value: [] }, featureCatalog: { value: null },
    linearSeqs: ['record-1', 'record-2'].map((uid) => ({ uid, gb: new Blob([uid]), region_record_id: '', region_reverse: false })) },
  { circular, linear });
  const controls = createRecordDisplayControls({ state, computed: (fn) => ({ get value() { return fn(); } }),
    watch: (_source, callback) => { callbacks.push(callback); }, linearRecordsFor: () => [],
    linearRecordStatusFor: () => 'ready', linearRecordErrorFor: () => '', runUndoable: (_label, fn) => fn(),
    getCommittedRequest: () => null, getCommittedSession: () => null });
  const uids = (drawing) => drawing.recordDisplayDrafts.map(({ sourceUid }) => sourceUid);
  return { state, controls, callbacks, circular, linear, uids };
};

test('no watcher of the record sources drops drafts (A2 W5)', () => {
  const { state, callbacks, circular, linear, uids } = sourceControls();
  const tuple = (uid, source) => ({ scope: 'linear', sourceUid: uid, source, paired: null });
  const previous = [tuple('record-1', state.linearSeqs[0].gb), tuple('record-2', state.linearSeqs[1].gb)];
  const current = [tuple('record-1', new Blob(['x'])), previous[1]];
  callbacks.forEach((callback) => callback(current, previous));
  assert.deepEqual(uids(linear), ['record-1', 'record-2']);
  assert.deepEqual(uids(circular), ['circular']);
});

test('a replaced Linear source drops its own drafts only, in its transition (A2 W5)', () => {
  const { state, controls, circular, linear, uids } = sourceControls();
  controls.changeRecordSources(() => { state.linearSeqs[0].gb = new Blob(['x']); });
  assert.deepEqual(uids(linear), ['record-2']);
  assert.deepEqual(uids(circular), ['circular']);
  controls.changeRecordSources(() => { state.cInputType.value = 'gff'; });
  assert.deepEqual(uids(circular), [], 'the input type switch replaced the Circular source');
  assert.deepEqual(uids(linear), ['record-2']);
});
