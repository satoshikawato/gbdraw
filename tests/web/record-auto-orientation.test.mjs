// "Orient records from comparisons": the writer reports what it did, and the
// button's status never says it changed a record that it did not write.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createRecordDisplayControls } from '../../gbdraw/web/js/app/record-display-options.js';
import { createRecordAutoOrientation } from '../../gbdraw/web/js/app/record-auto-orientation.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { adoptCurrentSessionResources, createSessionResourceFileView } from '../../gbdraw/web/js/services/session-resource-backing.js';
import { withModeDrawings } from './helpers/drawing-state.mjs';

const genbank = (name, text) => ({
  kind: 'genbank', name, type: 'text/plain', lastModified: 0, size: Buffer.byteLength(text), encoding: 'base64', data: btoa(text)
});
const resources = {
  a: genbank('a.gb', 'LOCUS       first 100 bp DNA linear\n//\n'),
  b: genbank('b.gb', 'LOCUS       second 100 bp DNA linear\n//\n')
};
const committed = {
  resources,
  renderRequest: {
    mode: 'linear',
    comparisons: [{ kind: 'nucleotideBlast' }],
    records: [
      { recordKey: 'ca', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'a' } },
      { recordKey: 'cb', cardinality: 'exactly_one', source: { kind: 'genbank', resourceId: 'b' } }
    ]
  }
};
const DECISIONS = [{ recordKey: 'ca', reverseComplement: false }, { recordKey: 'cb', reverseComplement: true }];

const controlsFor = ({ mode = 'linear', runUndoable, sessionBusy = null }) => {
  const table = adoptCurrentSessionResources(resources);
  const state = {
    mode: { value: 'linear' }, cInputType: { value: 'gb' }, lInputType: { value: 'gb' }, files: { c_gb: null },
    linearSeqs: ['a', 'b'].map((id) => ({
      uid: `c${id}`, gb: createSessionResourceFileView(table, id), region_record_id: '', region_reverse: false
    })),
    form: createDefaultForm(), adv: createDefaultAdv('linear'), recordDisplayDrafts: [], featurePlacementOverrides: {},
    circularRecordList: { value: [] },
    circularRecordDiscovery: { status: 'idle', error: '', inputType: '', primaryFile: null, pairedFile: null },
    featureCatalog: { value: { items: [] } },
    sessionOperationAvailability: () => sessionBusy
  };
  const drawing = (id) => ({
    id, mode: id, form: createDefaultForm(), adv: createDefaultAdv(id), recordDisplayDrafts: [], featurePlacementOverrides: {}
  });
  const wrapped = withModeDrawings(state, { circular: drawing('circular'), linear: drawing('linear') });
  state.mode.value = mode;
  const controls = createRecordDisplayControls({
    state: wrapped, computed: (fn) => ({ get value() { return fn(); } }), watch: () => {},
    linearRecordsFor: (seq) => [{
      selector: '#1', recordId: seq.uid === 'ca' ? 'first' : 'second', recordLength: 100, detectedTopology: 'linear'
    }],
    linearRecordStatusFor: () => 'ready', linearRecordErrorFor: () => '',
    runUndoable, getCommittedRequest: () => committed.renderRequest, getCommittedSession: () => committed
  });
  return { controls, reversed: () => state.linearSeqs.map((seq) => seq.region_reverse) };
};
const applies = async (_label, fn) => fn();

test('the writer writes the changed records in one step and names them', async () => {
  const labels = [];
  const { controls, reversed } = controlsFor({ runUndoable: (label, fn) => { labels.push(label); return applies(label, fn); } });
  assert.deepEqual(await controls.writeDecidedOrientations(DECISIONS, 'Orient'), { status: 'written', recordIds: ['second'] });
  assert.deepEqual(reversed(), [false, true]);
  assert.deepEqual(labels, ['Orient']);
});

test('a step History refuses is busy with its reason, and nothing is written', async () => {
  const reason = 'Saving session. Retry after saving finishes.';
  const { controls, reversed } = controlsFor({ runUndoable: async () => ({ status: 'busy', reason }) });
  assert.deepEqual(await controls.writeDecidedOrientations(DECISIONS, 'Orient'), { status: 'busy', reason });
  assert.deepEqual(reversed(), [false, false]);
});

test('a pending Session operation is busy before History is asked', async () => {
  const reason = 'Loading session. Retry after loading finishes.';
  let asked = false;
  const { controls, reversed } = controlsFor({
    sessionBusy: { status: 'busy', reason }, runUndoable: (label, fn) => { asked = true; return applies(label, fn); }
  });
  assert.deepEqual(await controls.writeDecidedOrientations(DECISIONS, 'Orient'), { status: 'busy', reason });
  assert.equal(asked, false);
  assert.deepEqual(reversed(), [false, false]);
});

test('a shown drawing that is not Linear has its own outcome', async () => {
  const { controls, reversed } = controlsFor({ mode: 'circular', runUndoable: applies });
  assert.deepEqual(await controls.writeDecidedOrientations(DECISIONS, 'Orient'), { status: 'notLinear' });
  assert.deepEqual(reversed(), [false, false]);
});

test('the status line follows the writer: a refused write never reads as changed', async () => {
  const ref = (value) => ({ value });
  const statusAfter = async (written) => {
    const orientation = createRecordAutoOrientation({
      computed: (fn) => ({ get value() { return fn(); } }), ref,
      getCommittedSession: () => committed, readCommittedRevision: () => 1, busy: () => null,
      decideOrientation: async () => ({ result: { records: DECISIONS } }),
      writeDecidedOrientations: async () => written, errorMessage: () => 'failed'
    });
    await orientation.orient();
    return orientation.status.value;
  };
  assert.match(await statusAfter({ status: 'written', recordIds: ['second'] }), /changed second/);
  assert.equal(await statusAfter({ status: 'busy', reason: 'Saving session. Retry after saving finishes.' }),
    'Saving session. Retry after saving finishes.');
  assert.match(await statusAfter({ status: 'notLinear' }), /shown drawing changed/);
  assert.match(await statusAfter({ status: 'unchanged' }), /current orientation/);
});
