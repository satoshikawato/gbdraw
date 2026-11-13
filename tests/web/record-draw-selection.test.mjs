// Per-record draw ON/OFF (record-selection D-01..D-10): the state-free
// answers of services/record-draw-selection.js.
import assert from 'node:assert/strict';
import test from 'node:test';
import { readFileSync } from 'node:fs';
import {
  circularRecordDrawKey,
  circularRecordsToDraw,
  compareRecordIdsNaturally,
  drawnLinearSequences,
  isRecordDrawKey,
  nextRecordsOff,
  offChangeLeavesSourceEmpty,
  offRecordRequestRecords,
  omittedLinearUids,
  pruneRecordsOff,
  recordIndexNumber,
  recordListRows,
  recordOrderRows,
  recordSelectorNumber
} from '../../gbdraw/web/js/services/record-draw-selection.js';

const cards = ['a', 'b', 'c', 'd'].map((uid) => ({ uid }));

test('the drawn Linear cards keep file order and leave out OFF cards', () => {
  assert.deepEqual(drawnLinearSequences(cards, []).map((card) => card.uid), ['a', 'b', 'c', 'd']);
  assert.deepEqual(drawnLinearSequences(cards, ['c', 'a']).map((card) => card.uid), ['b', 'd']);
  // The same card objects, not copies.
  assert.equal(drawnLinearSequences(cards, ['a'])[0], cards[1]);
  assert.deepEqual([...omittedLinearUids(cards, ['c', 'gone'])], ['c']);
});

test('the Circular records to draw are named by their identity', () => {
  const entries = [0, 1, 2].map((sourceIndex) => ({ identity: `card-a:${sourceIndex + 1}`, sourceIndex }));
  const drawn = circularRecordsToDraw(entries, ['card-a:2']);
  assert.deepEqual(drawn.map((entry) => entry.sourceIndex), [0, 2]);
  assert.equal(circularRecordDrawKey({ identity: 'card-b' }), 'card-b');
  // One shape of key in both modes: a card uid or a record identity.
  assert.equal(isRecordDrawKey('circular', 'card-a:3'), true);
  assert.equal(isRecordDrawKey('circular', ' card-a'), false);
  assert.equal(isRecordDrawKey('linear', 'linear-seq-1'), true);
  assert.equal(isRecordDrawKey('linear', ''), false);
});

test('turning records ON or OFF and pruning keep one entry per key', () => {
  assert.deepEqual(nextRecordsOff(['a'], ['b', 'a'], false), ['a', 'b']);
  assert.deepEqual(nextRecordsOff(['a', 'b'], ['a'], true), ['b']);
  assert.deepEqual(pruneRecordsOff(['a', 'gone', 'b'], ['a', 'b', 'c']), ['a', 'b']);
});

test('the request keys of OFF records follow the request builder', () => {
  // A Linear OFF card keeps its rows under its uid and its `<uid>:<k>` identities.
  // A Circular OFF record keeps them under its identity: one rule for both modes.
  assert.deepEqual(offRecordRequestRecords(['uid-2', 'uid-2']), [{ recordKey: 'uid-2', cardinality: 'all' }]);
  assert.deepEqual(offRecordRequestRecords(['card-a:2', 'card-b']), [
    { recordKey: 'card-a:2', cardinality: 'all' }, { recordKey: 'card-b', cardinality: 'all' }
  ]);
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

// OV-2960: a `#N` record selector reads as Python's renderer reads it (`record_selector_number` and
// `record_index_number` in gbdraw/io/record_select.py, whose pattern the generated contracts carry;
// tests/test_web_contracts_generated.py checks the same texts in Python).
test('a #N record selector names the record the renderer draws', () => {
  /** @type {[string, number | null, number | null][]} text, record, N with 0 */
  const cases = [
    ['#1', 1, 1], ['#01', 1, 1], [' # 2 ', 2, 2], ['#0', null, 0], ['#00', null, 0],
    ['#x', null, null], ['#', null, null], ['#1e1', null, null], ['#+1', 1, 1], ['# +02', 2, 2], ['#+ 1', null, null],
    ['#++1', null, null], ['#-1', null, null], ['#+0', null, 0], ['+#1', null, null], ['#0x2', null, null],
    ['#1.0', null, null], ['BGC0000708', null, null], ['#\u00853', 3, 3], ['#\u001c3', 3, 3], ['#\ufeff3', null, null]
  ];
  for (const [text, record, number] of cases) {
    assert.equal(recordSelectorNumber(text), record, JSON.stringify(text));
    assert.equal(recordIndexNumber(text), number, JSON.stringify(text));
  }
  assert.equal(recordSelectorNumber(1), null);
});

// OV-3313: the row of a Record Order token (`<selector>@<row>`, the text after its last `@`) reads as Python's
// renderer reads it (`parse_record_row_position` in gbdraw/layout/record_placement.py, whose pattern the generated
// contracts carry; tests/fixtures/record_order_row_vectors.json holds the tokens, which
// tests/test_web_contracts_generated.py checks in Python). `1.0`, `0x1` and `1e0`,
// which `Number()` read as 1, are no rows; a row is read from the token as stored (U+FEFF, which `trim()` strips and
// Python's `strip()` does not, stays and refuses it).
test('a Record Order token\'s row is the row the renderer reads', () => {
  /** @type {{ tokens: { token: string, row: number | null }[] }} */
  const vectors = JSON.parse(readFileSync(new URL('../fixtures/record_order_row_vectors.json', import.meta.url), 'utf8'));
  assert.ok(vectors.tokens.length > 20);
  for (const { token, row } of vectors.tokens) assert.deepEqual(recordOrderRows([token]), [row], JSON.stringify(token));
  assert.deepEqual(recordOrderRows([1, '#1']), [null, null]);
});

// OV-3313 review 1 L1: the renderer reads rows of up to 4300 digits exactly and only their order (`resolve_record_row_positions`
// ranks them; a longer row is refused). The rows of one list read as safe integers in the same order: each row's value
// when every row is a safe integer, else each row's rank, so `9007199254740993` stays after `9007199254740992`.
test('the rows of a Record Order list keep the order the renderer draws, beyond 2^53 and up to 4300 digits', () => {
  /** @type {{ rowLists: { tokens: string[], order: number[] | null }[] }} */
  const { rowLists } = JSON.parse(readFileSync(new URL('../fixtures/record_order_row_vectors.json', import.meta.url), 'utf8'));
  assert.ok(rowLists.length >= 8);
  for (const { tokens, order } of rowLists) {
    const rows = recordOrderRows(tokens);
    const at = tokens.map((token) => token.slice(0, 24)).join(' ');
    if (order === null) {
      assert.ok(rows.includes(null), at);
      continue;
    }
    assert.ok(rows.every((row) => Number.isSafeInteger(row) && row >= 1), at);
    const drawn = tokens.map((_token, index) => index).sort((left, right) => rows[left] - rows[right] || left - right);
    assert.deepEqual(drawn, order, at);
  }
  assert.deepEqual(recordOrderRows(['#1@9007199254740991', '#2@3']), [9007199254740991, 3]);
  assert.deepEqual(recordOrderRows(['#1@9007199254740993', '#2@3', '#3@3', '#4@x']), [2, 1, 1, null]);
});
