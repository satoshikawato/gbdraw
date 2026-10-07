import assert from 'node:assert/strict';
import test from 'node:test';
import {
  buildRecordDisplayRows, createRecordDisplayControls, parseRecordDisplayStart, reconcileRecordDisplayDrafts,
  RECORD_TARGET_NOT_DISCOVERED,
  effectiveRecordReverseComplement, migrateLegacyRecordDisplayDrafts, recordDisplaySurface,
  requestedRecordTransform, selectedFeatureDisplayStart, validateAnchorIntent,
  validateRecordDisplayDrafts
} from '../../gbdraw/web/js/app/record-display-options.js';
import { circularDiscoveryForInput, parseSequenceRecordText } from '../../gbdraw/web/js/app/record-discovery.js';
import { isCurrentFeature } from '../../gbdraw/web/js/services/feature-identity.js';

import { createFeaturePlacementActions } from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { adoptCurrentSessionResources, createCombinedSessionResourceFileView, createSessionResourceFileView } from '../../gbdraw/web/js/services/session-resource-backing.js';
import { readFileBytes } from '../../gbdraw/web/js/services/file-content-cache.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { withDrawings } from './helpers/drawing-state.mjs';

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
  const draft = {
    scope: 'linear', sourceUid: 'card-1', selector: '#1', recordId: 'same',
    topologyOverride: null, startCoordinate: 25, reverseComplementOverride: true,
    anchorIntent
  };
  const drafts = [draft];
  assert.equal(validateRecordDisplayDrafts(drafts), drafts);
  assert.equal(validateAnchorIntent(anchorIntent), anchorIntent);
  assert.deepEqual(migrateLegacyRecordDisplayDrafts([{
    scope: 'linear', sourceUid: 'card-1', selector: '#1', recordId: 'same',
    topologyOverride: null, startCoordinate: 25
  }])[0], { ...draft, reverseComplementOverride: null, anchorIntent: null });
  for (const invalid of [
    { ...draft, extra: true },
    { ...draft, reverseComplementOverride: 'true' },
    { ...draft, anchorIntent: { ...anchorIntent, schema: 2 } }
  ]) assert.throws(() => validateRecordDisplayDrafts([invalid]));
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

for (const strand of ['+', '-']) {
  for (const parts of [[[0, 1]], [[99, 100]], [[10, 14]], [[90, 100], [0, 5]], [[10, 13], [40, 43]]]) {
    test(`covered-base shortcut oracle ${strand} ${JSON.stringify(parts)}`, () => {
      const ordered = strand === '-' ? [...parts].reverse() : parts;
      const bases = ordered.flatMap(([start, end]) => {
        const sequence = Array.from({ length: end - start }, (_, i) => start + i + 1);
        return strand === '-' ? sequence.reverse() : sequence;
      });
      const feature = { record_key: 'instance', biological_feature_id: 'feature', strand,
        anchorProfile: { precision: 'exact', operator: parts.length === 1 ? 'single' : 'join',
          partOrder: 'biological', strand },
        location_parts: ordered.map(([start, end]) => ({ start, end, strand })) };
      const args = { row: rows[0], committedRow: { ...rows[0], recordKey: 'instance' }, feature };
      assert.equal(selectedFeatureDisplayStart({ ...args, shortcut: 'five-prime' }), bases[0]);
      assert.equal(selectedFeatureDisplayStart({ ...args, shortcut: 'midpoint' }), bases[Math.floor((bases.length - 1) / 2)]);
    });
  }
}

test('shortcuts reject stale, unbound, cross-record, empty, and unknown/mixed-strand selections', () => {
  const feature = { record_key: 'instance', biological_feature_id: 'feature', strand: '+',
    anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' },
    location_parts: [{ start: 5, end: 10, strand: '+' }] };
  const args = { row: rows[0], committedRow: { ...rows[0], recordKey: 'instance' }, feature, shortcut: 'midpoint' };
  for (const changed of [
    { committedRow: null }, { committedRow: { ...args.committedRow, source: {} } },
    { committedRow: { ...args.committedRow, selector: '#2' } },
    { feature: null },
    { feature: { ...feature, record_key: 'other' } },
    { feature: { ...feature, anchorProfile: { ...feature.anchorProfile, strand: 'mixed' } } },
    { feature: { ...feature, location_parts: [] } },
    { feature: { ...feature, location_parts: [{ start: 5, end: 10, strand: '-' }] } }
  ]) assert.throws(() => selectedFeatureDisplayStart({ ...args, ...changed }));
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
  const state = { mode: { value: linear ? 'linear' : 'circular' }, cInputType: { value: 'gb' },
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
  assert.throws(() => model.actions.setPlacement([model.feature], 'main'), /Unavailable/);
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

test('unsupported draft lanes remain unavailable independently of current placement value', () => {
  const model = compositeControls();
  model.actions.setPlacement([model.feature], 'outward');
  model.state.form.track_type = 'tuckin';
  assert.equal(model.actions.valueFor(model.feature), 'outward');
  assert.deepEqual(model.enabled(), ['auto', 'main']);
  assert.throws(() => model.actions.setPlacement([model.feature], 'inward'), /Unavailable/);
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
  const state = { cInputType: { value: 'gb' }, files: { c_gb: source },
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
