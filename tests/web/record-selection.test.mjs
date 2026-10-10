// The record selection owner (app/record-selection.js, D-02..D-10): toggles,
// bulk actions on the rows shown, the question before the last drawn record
// of a file goes OFF (D-06), Delete settings (D-08), and the list an upload
// opens (D-04). History, the dialog choice, and the root transitions are
// fakes here; history-config-restore.test.mjs runs Delete settings through the
// real History.
import assert from 'node:assert/strict';
import test from 'node:test';
import { createRecordSelection } from '../../gbdraw/web/js/app/record-selection.js';
import {
  RECORD_SETTINGS_DEFAULTS,
  linearCardHasSettings,
  recordsWithEdits,
  removeRecordEdits
} from '../../gbdraw/web/js/services/record-draw-selection.js';

const computed = (getter) => ({ get value() { return getter(); } });
const reactive = (value) => value;

const makeOwner = ({ sources, recordsOff = { linear: [], circular: [] }, requests = [] } = {}) => {
  const log = [];
  const offLists = { linear: [...recordsOff.linear], circular: [...recordsOff.circular] };
  const choice = { pending: false };
  const owner = createRecordSelection({
    reactive,
    computed,
    recordsOff: (mode) => offLists[mode],
    source: (mode, key) => sources[mode]?.find((entry) => entry.key === key) || null,
    runUndoable: async (label, change) => { log.push(`step:${label}`); return change(); },
    withDialogChoice: (label, handler, cancel) => (value) => {
      if (value === 'cancel') return cancel();
      log.push(`step:${label()}`);
      return handler(value);
    },
    closeAfterDialogChoice: (close) => close(),
    afterRecordSetChange: (mode) => log.push(`invalidate:${mode}`),
    removeSource: (mode, key) => { log.push(`remove:${mode}:${key}`); return true; },
    deleteRecordSettings: (mode, keys) => { log.push(`delete:${mode}:${keys.join(',')}`); return true; },
    autoOpenRequests: requests
  });
  return { owner, offLists, log, choice };
};

const linearFile = (key, count, name = `${key}.gb`) => ({
  key, name,
  records: Array.from({ length: count }, (_, index) => ({
    key: `${key}-${index + 1}`, recordId: `contig_${index + 1}`, length: (index + 1) * 100
  }))
});
const checkbox = (checked) => ({ target: { checked } });

test('a record checkbox turns one record OFF and back ON, one step each', async () => {
  const { owner, offLists, log } = makeOwner({ sources: { linear: [linearFile('f', 3)] } });
  const event = checkbox(false);
  await owner.toggleRecord('linear', 'f', 'f-2', event);
  assert.deepEqual(offLists.linear, ['f-2']);
  assert.equal(event.target.checked, false);
  assert.equal(owner.isDrawn('linear', 'f-2'), false);
  assert.equal(owner.drawnCount('linear', ['f-1', 'f-2', 'f-3']), 2);
  await owner.toggleRecord('linear', 'f', 'f-2', checkbox(true));
  assert.deepEqual(offLists.linear, []);
  assert.deepEqual(log, ['step:Leave out records', 'invalidate:linear', 'step:Draw records', 'invalidate:linear']);
});

test('the last drawn record of a file asks to remove the file; Cancel changes nothing (D-06)', async () => {
  const { owner, offLists, log } = makeOwner({
    sources: { linear: [linearFile('f', 2, 'draft.gb')] }, recordsOff: { linear: ['f-1'], circular: [] }
  });
  const event = checkbox(false);
  const outcome = await owner.toggleRecord('linear', 'f', 'f-2', event);
  assert.deepEqual(outcome, { status: 'question' });
  // Nothing written, and the checkbox shows the record ON again.
  assert.deepEqual(offLists.linear, ['f-1']);
  assert.equal(event.target.checked, true);
  assert.deepEqual({ ...owner.removeFileDialog }, { open: true, mode: 'linear', sourceKey: 'f', name: 'draft.gb' });
  owner.resolveRemoveFile('cancel');
  assert.equal(owner.removeFileDialog.open, false);
  assert.deepEqual(offLists.linear, ['f-1']);
  assert.deepEqual(log, []);
});

test('Remove File removes the whole file in one step and closes its list (D-06)', async () => {
  const { owner, log } = makeOwner({
    sources: { circular: [linearFile('circular', 2, 'genome.gb')] }, recordsOff: { linear: [], circular: ['circular-1'] }
  });
  owner.openRecordList('circular', 'circular');
  await owner.toggleRecord('circular', 'circular', 'circular-2', checkbox(false));
  assert.equal(owner.removeFileDialog.open, true);
  owner.resolveRemoveFile('remove');
  assert.deepEqual(log, ['step:Remove File', 'remove:circular:circular']);
  assert.equal(owner.removeFileDialog.open, false);
  assert.equal(owner.recordListView.value.open, false);
});

test('Select all and Select none act on the rows shown; a bulk change that empties the file applies nothing', async () => {
  const { owner, offLists, log } = makeOwner({ sources: { linear: [linearFile('f', 12)] } });
  owner.openRecordList('linear', 'f');
  owner.setRecordListQuery('contig_1');
  // contig_1, contig_10, contig_11, contig_12
  assert.equal(owner.recordListView.value.rows.length, 4);
  await owner.selectShown(false);
  assert.deepEqual([...offLists.linear].sort(), ['f-1', 'f-10', 'f-11', 'f-12']);
  assert.equal(owner.recordListView.value.drawn, 8);
  owner.setRecordListQuery('');
  await owner.selectShown(false);
  // Every record of the file would be OFF: ask, write nothing.
  assert.equal(owner.removeFileDialog.open, true);
  assert.equal(offLists.linear.length, 4);
  owner.resolveRemoveFile('cancel');
  await owner.selectShown(true);
  assert.deepEqual(offLists.linear, []);
  assert.deepEqual(log.filter((entry) => entry.startsWith('step:')), [
    'step:Leave out records', 'step:Draw records'
  ]);
});

test('the list sorts for display only, numbers as numbers (D-10)', () => {
  const file = { key: 'f', name: 'f.gb', records: [
    { key: 'a', recordId: 'contig_10', length: 5 },
    { key: 'b', recordId: 'contig_2', length: 50 },
    { key: 'c', recordId: 'contig_1', length: 20 }
  ] };
  const { owner } = makeOwner({ sources: { linear: [file] } });
  owner.openRecordList('linear', 'f');
  const order = () => owner.recordListView.value.rows.map((row) => row.recordId);
  assert.deepEqual(order(), ['contig_10', 'contig_2', 'contig_1']);
  owner.setRecordListSort('id-asc');
  assert.deepEqual(order(), ['contig_1', 'contig_2', 'contig_10']);
  owner.setRecordListSort('id-desc');
  assert.deepEqual(order(), ['contig_10', 'contig_2', 'contig_1']);
  owner.setRecordListSort('length-desc');
  assert.deepEqual(order(), ['contig_2', 'contig_1', 'contig_10']);
  owner.setRecordListSort('length-asc');
  assert.deepEqual(order(), ['contig_10', 'contig_1', 'contig_2']);
  owner.setRecordListSort('not-a-sort');
  assert.equal(owner.recordList.sort, 'length-asc');
  owner.setRecordListSort('file');
  assert.deepEqual(order(), ['contig_10', 'contig_2', 'contig_1']);
  // The file keeps its order.
  assert.deepEqual(file.records.map((record) => record.key), ['a', 'b', 'c']);
});

test('Delete settings reaches only OFF records, as one step (D-08)', async () => {
  const file = linearFile('f', 3);
  file.records[0].hasSettings = true;
  file.records[1].hasSettings = true;
  const { owner, log } = makeOwner({ sources: { linear: [file] }, recordsOff: { linear: ['f-1', 'f-3'], circular: [] } });
  assert.equal(await owner.deleteSettings('linear', ['f-2']), false);
  await owner.deleteSettings('linear', ['f-1', 'f-2']);
  owner.openRecordList('linear', 'f');
  // The OFF rows shown with settings: f-1 only (f-3 has none, f-2 is ON).
  assert.deepEqual(owner.recordListView.value.offWithSettings, ['f-1']);
  await owner.deleteShownOffSettings();
  assert.deepEqual(log, ['step:Delete record settings', 'delete:linear:f-1', 'step:Delete record settings', 'delete:linear:f-1']);
});

test('an upload with many records opens its list; several open one at a time (D-04)', () => {
  const requests = [];
  const { owner } = makeOwner({ sources: { linear: [linearFile('a', 25), linearFile('b', 30)] }, requests });
  assert.equal(owner.recordListView.value.open, false);
  requests.push({ mode: 'linear', sourceKey: 'gone' }, { mode: 'linear', sourceKey: 'a' }, { mode: 'linear', sourceKey: 'b' });
  assert.equal(owner.recordListView.value.sourceKey, 'a');
  assert.equal(owner.recordListView.value.drawn, 25);
  owner.closeRecordList();
  assert.equal(owner.recordListView.value.sourceKey, 'b');
  // Choose records… on another file shows that one first; the request waits.
  owner.openRecordList('linear', 'a');
  assert.equal(owner.recordListView.value.sourceKey, 'a');
  owner.closeRecordList();
  assert.equal(owner.recordListView.value.sourceKey, 'b');
  owner.closeRecordList();
  assert.equal(owner.recordListView.value.open, false);
  assert.deepEqual(requests, []);
});

test('Delete settings removes the record edits, annotations, pairs, and rotation of its records only', () => {
  const featureKey = (recordKey, id) => JSON.stringify([recordKey, id]);
  const data = {
    featureOverrides: {
      [featureKey('uid-a', 'f1')]: { recordKey: 'uid-a', biologicalFeatureId: 'f1' },
      [featureKey('uid-a:2', 'f2')]: { recordKey: 'uid-a:2', biologicalFeatureId: 'f2' },
      [featureKey('uid-b', 'f1')]: { recordKey: 'uid-b', biologicalFeatureId: 'f1' }
    },
    featurePlacementOverrides: {
      [featureKey('uid-a', 'f3')]: { recordKey: 'uid-a', biologicalFeatureId: 'f3', placement: 'lane' }
    },
    featureStrokeOverrides: { 'uid-a\0f1': { color: '#000' }, 'uid-b\0f1': { color: '#111' } },
    annotationSets: [{ annotations: [
      { id: 'bound-a', target: { kind: 'coordinateSpan' }, metadata: { binding: 'cat-a' } },
      { id: 'feature-a', target: { kind: 'featureIdentity', recordKey: 'uid-a' }, metadata: {} },
      { id: 'bound-b', target: { kind: 'coordinateSpan' }, metadata: { binding: 'cat-b' } },
      { id: 'unbound', target: { kind: 'coordinateSpan' }, metadata: {} }
    ] }],
    recordDisplayDrafts: [{ sourceUid: 'uid-a', selector: '#1' }, { sourceUid: 'uid-b', selector: '#1' }],
    comparisonEdges: [
      { queryUid: 'uid-a', subjectUid: 'uid-b' }, { queryUid: 'uid-b', subjectUid: 'uid-c' }, { queryUid: 'uid-c', subjectUid: 'uid-a' }
    ],
    annotationBindingField: 'binding'
  };
  const owners = [{
    key: 'uid-a', requestKeys: ['uid-a'], ownsExpansions: true, bindingKeys: ['cat-a'],
    displaySource: 'uid-a', displaySelector: null
  }];
  assert.deepEqual([...recordsWithEdits(owners, data)], ['uid-a']);
  assert.equal(removeRecordEdits(owners, data), 9);
  assert.deepEqual(Object.values(data.featureOverrides).map((row) => row.recordKey), ['uid-b']);
  assert.deepEqual(data.featurePlacementOverrides, {});
  assert.deepEqual(Object.keys(data.featureStrokeOverrides), ['uid-b\0f1']);
  assert.deepEqual(data.annotationSets[0].annotations.map((item) => item.id), ['bound-b', 'unbound']);
  assert.deepEqual(data.recordDisplayDrafts, [{ sourceUid: 'uid-b', selector: '#1' }]);
  assert.deepEqual(data.comparisonEdges, [{ queryUid: 'uid-b', subjectUid: 'uid-c' }]);
  assert.equal(recordsWithEdits(owners, data).size, 0);

  // A Circular record: its `record-N` and single-record keys and its own rotation row.
  const circular = {
    featureOverrides: {
      [featureKey('record-2', 'f1')]: { recordKey: 'record-2' },
      [featureKey('circular-chr2-_2', 'f1')]: { recordKey: 'circular-chr2-_2' },
      [featureKey('record-1', 'f1')]: { recordKey: 'record-1' },
      // `record-2` is not an expandable key, so `record-2:1` is another record.
      [featureKey('record-2:1', 'f1')]: { recordKey: 'record-2:1' }
    },
    recordDisplayDrafts: [{ sourceUid: 'circular', selector: '#1' }, { sourceUid: 'circular', selector: '#2' }],
    annotationBindingField: 'binding'
  };
  removeRecordEdits([{
    key: '#2', requestKeys: ['record-2', 'circular-chr2-_2'], bindingKeys: [],
    displaySource: 'circular', displaySelector: '#2'
  }], circular);
  assert.deepEqual(Object.values(circular.featureOverrides).map((row) => row.recordKey), ['record-1', 'record-2:1']);
  assert.deepEqual(circular.recordDisplayDrafts, [{ sourceUid: 'circular', selector: '#1' }]);
});

test('a card holds settings when a field Delete settings resets differs from its default', () => {
  assert.equal(linearCardHasSettings({ uid: 'a', definition: '', region_start: null, losat_gencode: 1 }), false);
  assert.equal(linearCardHasSettings({ ...RECORD_SETTINGS_DEFAULTS, record_subtitle: 'x' }), true);
  assert.equal(linearCardHasSettings({ ...RECORD_SETTINGS_DEFAULTS, losat_gencode: 11 }), true);
  assert.equal(linearCardHasSettings({ ...RECORD_SETTINGS_DEFAULTS, region_reverse: true }), true);
});
