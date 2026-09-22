import assert from 'node:assert/strict';
import test from 'node:test';
import {
  buildRecordDisplayRows, createRecordDisplayControls, parseRecordDisplayStart, reconcileRecordDisplayDrafts,
  effectiveRecordReverseComplement, migrateLegacyRecordDisplayDrafts, recordDisplaySurface,
  requestedRecordTransform, selectedFeatureDisplayStart, validateAnchorIntent,
  validateRecordDisplayDrafts
} from '../../gbdraw/web/js/app/record-display-options.js';
import { parseSequenceRecordText } from '../../gbdraw/web/js/app/record-discovery.js';

import { createFeaturePlacementActions } from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { adoptCurrentSessionResources, createCombinedSessionResourceFileView } from '../../gbdraw/web/js/services/session-resource-backing.js';
import { readFileBytes } from '../../gbdraw/web/js/services/file-content-cache.js';

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
const compositeControls = () => {
  const resources = Object.fromEntries(['first', 'second'].map((id) => [id,
    descriptor(`${id}.gb`, `LOCUS       ${id} 100 bp DNA circular\n//\n`)]));
  const table = adoptCurrentSessionResources(resources);
  const components = Object.keys(resources).map((resourceId) => ({ resourceId }));
  const makeFile = () => createCombinedSessionResourceFileView(table, components);
  const file = makeFile();
  let committed = { resources, renderRequest: { mode: 'circular', records: components.map(({ resourceId }, i) => ({
    recordKey: `record-${i + 1}`, cardinality: 'exactly_one', source: { kind: 'genbank', resourceId },
    selector: null })) } };
  const state = { mode: { value: 'circular' }, cInputType: { value: 'gb' },
    lInputType: { value: 'gb' }, files: { c_gb: file }, linearSeqs: [],
    form: createDefaultForm(), adv: createDefaultAdv('circular'),
    recordDisplayDrafts: [], featurePlacementOverrides: {},
    circularRecordList: { value: components.map(({ resourceId }, i) => ({
      record_id: resourceId, record_length: 100, selector: `#${i + 1}`, detectedTopology: 'circular' })) },
    featureCatalog: { value: { items: [{ recordKeys: ['record-1', 'record-2'] }] } } };
  state.form.track_type = 'middle';
  const getCommittedRequest = () => committed.renderRequest;
  const history = { runUndoable: (_label, fn) => fn() };
  const controls = createRecordDisplayControls({ state, computed: (fn) => ({ get value() { return fn(); } }),
    watch: () => {}, linearRecordSelector: { recordsFor: () => [] }, history,
    getCommittedRequest, getCommittedSession: () => committed });
  const actions = createFeaturePlacementActions({ state, history, getCommittedRequest,
    isCurrentFeature: controls.isCurrentFeature });
  const feature = { record_key: 'record-2', biological_feature_id: 'logical-feature' };
  return { state, actions, controls, feature, file, makeFile,
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

test('same-name, same-content source replacement does not retain old feature capability', () => {
  const model = compositeControls();
  assert.equal(model.enabled().length, 4);
  model.state.files.c_gb = model.makeFile();
  assert.equal(model.state.files.c_gb.name, model.file.name);
  assert.notEqual(model.state.files.c_gb, model.file);
  assert.deepEqual(model.enabled(), []);
  assert.throws(() => model.actions.setPlacement([model.feature], 'main'), /Unavailable/);
});

test('unsupported draft lanes remain unavailable independently of current placement value', () => {
  const model = compositeControls();
  model.actions.setPlacement([model.feature], 'outward');
  model.state.form.track_type = 'tuckin';
  assert.equal(model.actions.valueFor(model.feature), 'outward');
  assert.deepEqual(model.enabled(), ['auto', 'main']);
  assert.throws(() => model.actions.setPlacement([model.feature], 'inward'), /Unavailable/);
});
