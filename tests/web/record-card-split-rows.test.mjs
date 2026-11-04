import assert from 'node:assert/strict';
import test from 'node:test';
import {
  buildCanonicalRequestState, buildCanonicalRenderRequest, projectCommittedRecordTransform
} from '../../gbdraw/web/js/services/session-request.js';
import { createDefaultForm, createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { buildRecordDisplayRows } from '../../gbdraw/web/js/services/record-display-model.js';
import { resolveLinearComparisonPlan } from '../../gbdraw/web/js/services/linear-comparisons.js';

// OV-2500, Owner decision R17-41 (Q-RC-1 = B): a split of a File card that
// draws several records leaves every member's row empty while Arrange in rows
// is off, so the planner gives each record its own row; while it is on, the
// members keep the card's row.

const text = 'LOCUS       same 100 bp DNA circular\n//\nLOCUS       same 100 bp DNA circular\n//\n';
const file = { name: 'same.gb', size: text.length, type: 'text/plain', lastModified: 0, data: btoa(text) };
const other = { ...file, name: 'other.gb' };

/** The canonical Session of a Linear drawing whose `card` draws two records, with these starts. */
const buildFor = (arranged, starts) => {
  const { state, drawing } = buildCanonicalRequestState({
    session: { renderRequest: { records: [] } },
    projection: { mode: 'linear', inputType: 'gb', files: {}, config: {} },
    config: { form: createDefaultForm(), adv: createDefaultAdv('linear'), colors: {} }
  });
  Object.assign(state, drawing);
  state.linearRecordLayoutEnabled.value = arranged;
  state.linearRecordRows.splice(0, state.linearRecordRows.length, { uid: 'first', row: 1 }, { uid: 'card', row: 2 });
  const rows = buildRecordDisplayRows({ scope: 'linear', sourceUid: 'card', source: {},
    records: [1, 2].map(() => ({ recordId: 'same', recordLength: 100, detectedTopology: 'circular' })) });
  state.recordDisplayDrafts = starts.map((startCoordinate, index) => ({ sourceUid: 'card',
    selector: rows[index].selector, recordId: rows[index].recordId, topologyOverride: null, startCoordinate,
    reverseComplementOverride: null, anchorIntent: null }));
  const filesData = { linearSeqs: [{ uid: 'first', gb: other, region_record_id: 'same' }, { uid: 'card', gb: file }] };
  return buildCanonicalRenderRequest({ state, drawing: state, filesData, recordDisplayRows: rows,
    comparisonPlanSnapshot: resolveLinearComparisonPlan({ plan: state.linearComparisonPlan,
      sequences: filesData.linearSeqs, layout: [], losatProgram: 'blastn', blastpMode: 'orthogroup' }) });
};
// Different starts make the request builder split the card.
const splitRequest = (arranged) => buildFor(arranged, [1, 71]).renderRequest;
const rowsOf = (request) => request.records.map((record) => [record.recordKey, record.presentation.gridRow]);

test('the request builder\'s split leaves rows to the planner while Arrange in rows is off', () => {
  assert.deepEqual(rowsOf(splitRequest(false)), [['first', null], ['card:1', null], ['card:2', null]]);
});

test('the request builder\'s split keeps the card\'s row while Arrange in rows is on', () => {
  assert.deepEqual(rowsOf(splitRequest(true)), [['first', 1], ['card:1', 2], ['card:2', 2]]);
});

// The same starts keep the card whole; the Rotate popup then splits it.
const rotate = (arranged) => {
  const committed = buildFor(arranged, [null, null]);
  const card = committed.renderRequest.records[1];
  assert.equal(card.cardinality, 'all');
  const members = [1, 2].map((index) => ({
    canonicalRecordKey: 'card', recordKey: `card:${index}`, selector: `#${index}`, recordId: 'same',
    recordLength: 100, committedDisplay: structuredClone(card.display), committedReverseComplement: false
  }));
  return projectCommittedRecordTransform({
    committed,
    target: { ...members[1], scope: 'linear', source: structuredClone(card.source),
      effectiveCircular: true, cropped: false, members },
    transform: { recordLength: 100, startCoordinate: 40, reverseComplement: false }
  }).canonical.renderRequest;
};

test('the Rotate popup\'s split leaves rows to the planner while Arrange in rows is off', () => {
  assert.deepEqual(rowsOf(rotate(false)), [['first', null], ['card:1', null], ['card:2', null]]);
});

test('the Rotate popup\'s split keeps the card\'s row while Arrange in rows is on', () => {
  assert.deepEqual(rowsOf(rotate(true)), [['first', 1], ['card:1', 2], ['card:2', 2]]);
});
