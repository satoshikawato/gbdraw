// Per-record draw ON/OFF (record-selection D-01..D-10): the state-free
// answers of services/record-draw-selection.js.
import assert from 'node:assert/strict';
import test from 'node:test';
import {
  circularRecordDrawKey,
  circularRecordsToDraw,
  compareRecordIdsNaturally,
  drawnLinearSequences,
  isRecordDrawKey,
  nextRecordsOff,
  offChangeLeavesSourceEmpty,
  omittedLinearUids,
  pruneRecordsOff,
  recordListRows
} from '../../gbdraw/web/js/services/record-draw-selection.js';

const cards = ['a', 'b', 'c', 'd'].map((uid) => ({ uid }));

test('the drawn Linear cards keep file order and leave out OFF cards', () => {
  assert.deepEqual(drawnLinearSequences(cards, []).map((card) => card.uid), ['a', 'b', 'c', 'd']);
  assert.deepEqual(drawnLinearSequences(cards, ['c', 'a']).map((card) => card.uid), ['b', 'd']);
  // The same card objects, not copies.
  assert.equal(drawnLinearSequences(cards, ['a'])[0], cards[1]);
  assert.deepEqual([...omittedLinearUids(cards, ['c', 'gone'])], ['c']);
});

test('the Circular records to draw keep their source index', () => {
  const entries = [0, 1, 2].map((sourceIndex) => ({ selector: `#${sourceIndex + 1}`, sourceIndex }));
  const drawn = circularRecordsToDraw(entries, ['#2']);
  assert.deepEqual(drawn.map((entry) => entry.sourceIndex), [0, 2]);
  assert.equal(circularRecordDrawKey({ sourceIndex: 4 }), '#5');
  assert.equal(isRecordDrawKey('circular', '#3'), true);
  assert.equal(isRecordDrawKey('circular', '#0'), false);
  assert.equal(isRecordDrawKey('circular', 'contig_1'), false);
  assert.equal(isRecordDrawKey('linear', 'linear-seq-1'), true);
  assert.equal(isRecordDrawKey('linear', ''), false);
});

test('turning records ON or OFF and pruning keep one entry per key', () => {
  assert.deepEqual(nextRecordsOff(['a'], ['b', 'a'], false), ['a', 'b']);
  assert.deepEqual(nextRecordsOff(['a', 'b'], ['a'], true), ['b']);
  assert.deepEqual(pruneRecordsOff(['a', 'gone', 'b'], ['a', 'b', 'c']), ['a', 'b']);
});

test('the last ON record of a source cannot go OFF without asking (D-06)', () => {
  const sourceKeys = ['a', 'b', 'c'];
  assert.equal(offChangeLeavesSourceEmpty({ sourceKeys, recordsOff: ['a', 'b'], keys: ['c'], drawn: false }), true);
  assert.equal(offChangeLeavesSourceEmpty({ sourceKeys, recordsOff: ['a'], keys: ['c'], drawn: false }), false);
  // A bulk Select none on the whole source.
  assert.equal(offChangeLeavesSourceEmpty({ sourceKeys, recordsOff: [], keys: sourceKeys, drawn: false }), true);
  assert.equal(offChangeLeavesSourceEmpty({ sourceKeys, recordsOff: ['a', 'b', 'c'], keys: ['a'], drawn: true }), false);
  // The Circular source: its records' selectors.
  assert.equal(offChangeLeavesSourceEmpty({ sourceKeys: ['#1', '#2'], recordsOff: ['#1'], keys: ['#2'], drawn: false }), true);
});

const records = [
  { key: 'k1', recordId: 'contig_10', length: 500 },
  { key: 'k2', recordId: 'contig_2', length: 1500 },
  { key: 'k3', recordId: 'Plasmid_A', length: null },
  { key: 'k4', recordId: 'contig_1', length: 900, hasSettings: true }
];
const ids = (rows) => rows.map((row) => row.recordId);

test('the record list filters by record ID and sorts for display only (D-10)', () => {
  assert.deepEqual(ids(recordListRows({ records, recordsOff: [] })), ['contig_10', 'contig_2', 'Plasmid_A', 'contig_1']);
  assert.deepEqual(ids(recordListRows({ records, recordsOff: [], sort: 'id-asc' })), ['contig_1', 'contig_2', 'contig_10', 'Plasmid_A']);
  assert.deepEqual(ids(recordListRows({ records, recordsOff: [], sort: 'id-desc' })), ['Plasmid_A', 'contig_10', 'contig_2', 'contig_1']);
  assert.deepEqual(ids(recordListRows({ records, recordsOff: [], sort: 'length-desc' })), ['contig_2', 'contig_1', 'contig_10', 'Plasmid_A']);
  assert.deepEqual(ids(recordListRows({ records, recordsOff: [], sort: 'length-asc' })), ['contig_10', 'contig_1', 'contig_2', 'Plasmid_A']);
  const filtered = recordListRows({ records, recordsOff: ['k2'], query: 'CONTIG_1', sort: 'file' });
  assert.deepEqual(filtered.map(({ recordId, drawn, hasSettings, position }) => ({ recordId, drawn, hasSettings, position })), [
    { recordId: 'contig_10', drawn: true, hasSettings: false, position: 0 },
    { recordId: 'contig_1', drawn: true, hasSettings: true, position: 3 }
  ]);
  assert.equal(recordListRows({ records, recordsOff: ['k2'] })[1].drawn, false);
  // The input order is untouched.
  assert.deepEqual(records.map((record) => record.key), ['k1', 'k2', 'k3', 'k4']);
});

test('record IDs compare numbers as numbers', () => {
  assert.ok(compareRecordIdsNaturally('contig_2', 'contig_10') < 0);
  assert.ok(compareRecordIdsNaturally('seq9b', 'seq10a') < 0);
  assert.ok(compareRecordIdsNaturally('a', 'B') < 0);
  assert.equal(compareRecordIdsNaturally('x1', 'x1'), 0);
});
