// The name the Circular record panel shows for a Record select value or a Record
// Order row (app/circular-record-panel.js, OV-2949, OV-2959, D-98, D-100): a read
// record and an unread record identity that the shown Result names
// (`shownRecordIds`, tests/web/shown-drawing-readers.test.mjs) are named from one
// entry (its card's file name when the drawing has several cards, its record ID,
// and its `#k` when the ID repeats in its file, is unsafe, or, unread, is not
// proven unique); a Record Order row of a one-card drawing is `#k (<record ID>)`;
// a record the Result does not name is named by its card's file name (with `#k`
// for the k-th record of a file), never by the identity itself; any other value
// shows as it is.
// tests/web/circular-record-presentation.playwright.spec.js checks a loaded
// Session in the browser.
import assert from 'node:assert/strict';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

globalThis.window = { Vue: loadPageVue() };
const { createCircularRecordPresentation } = await import('../../gbdraw/web/js/app/circular-record-panel.js');

const file = (name) => ({ fileId: `id-${name}`, name });
const ONE = 'circular-source-1111';
const TWO = 'circular-source-2222';

/**
 * The panel over one Circular drawing whose records are not read (a loaded Session) or are read; its cards
 * hold the files `files` names by card uid.
 * `named` gives the record ID the Result names a record by; a record its committed request selected by that ID is
 * named by `{ recordId, selectedById: true }`.
 * @param {{ selector: string, named?: Record<string, string | { recordId: string, selectedById: boolean }>, read?: object[] | null, files?: Record<string, string> }} input
 */
const panelOf = ({ selector, named = {}, read = null, files = { [ONE]: 'one.gb', [TWO]: 'two.gb' } }) => {
  const sources = Object.entries(files).map(([uid, name]) => ({ uid, gb: file(name), gff: null, fasta: null }));
  const drawing = {
    id: 'circular', mode: 'circular', inputType: { value: 'gb' }, circularSources: sources,
    form: { circular_record_selector: selector, multi_record_canvas: false }, adv: { circular_grouping_intent: '' },
    recordsOff: [], circularRecordList: { value: read ?? [] },
    circularRecordDiscovery: read
      ? { status: 'ready', error: '', inputType: 'gb', sources: sources.map((card) => ({ uid: card.uid, primaryFile: card.gb, pairedFile: null, error: '' })) }
      : { status: 'idle', error: '', inputType: '', sources: [] }
  };
  const state = { activeDrawing: () => drawing, semanticFileWatchersSuppressed: { value: false } };
  const panel = createCircularRecordPresentation({
    state, mode: { value: 'circular' }, semanticMutationAvailable: { value: true }, sessionOperationAvailability: () => null,
    shownDrawingOfMode: () => drawing, refreshCircularRecordOrder: () => undefined,
    shownRecordIds: (shown) => new Map(shown.id === drawing.id ? Object.entries(named).map(([key, name]) => [
      key, typeof name === 'string' ? { recordId: name, selectedById: false } : name
    ]) : [])
  });
  return { panel, synthetic: () => panel.circularRecordPresentationOptions.value.find((option) => option.synthetic)?.label };
};

const byId = (recordId) => ({ recordId, selectedById: true });
const readRecord = (uid, k, recordId, identity = `${uid}:${k}`) => ({
  uid, k, identity, selector: `#${k}`, record_id: recordId, record_length: 100, detectedTopology: 'circular'
});

test('an unread record identity shows the record ID its Result names it by, else its file name', () => {
  // The k-th record of a file: its `#k` unless the request selected it by its record ID (review 1 of #1778).
  assert.equal(panelOf({ selector: `${TWO}:2`, named: { [`${TWO}:2`]: byId('BGC0000709') } }).synthetic(),
    'two.gb: BGC0000709 (not inspected)');
  assert.equal(panelOf({ selector: `${TWO}:2`, named: { [`${TWO}:2`]: 'BGC0000709' } }).synthetic(),
    'two.gb: BGC0000709 [#2] (not inspected)');
  // One card: the record ID the Result drew it with (OV-2959); the file only while no Result names it.
  const one = { [ONE]: 'presentation.gb' };
  assert.equal(panelOf({ selector: ONE, named: { [ONE]: 'NC_012920.1' }, files: one }).synthetic(), 'NC_012920.1 (not inspected)');
  assert.equal(panelOf({ selector: ONE, files: one }).synthetic(), 'presentation.gb (not inspected)');
  // No Result: the k-th record of the file.
  assert.equal(panelOf({ selector: `${TWO}:2` }).synthetic(), 'two.gb #2 (not inspected)');
  // A value that is not an identity (a kept `#N`, a record ID) shows as it is.
  assert.equal(panelOf({ selector: '#3' }).synthetic(), '#3 (not inspected)');
});
test('a read identity that no read record has is not found by its file name, in the option and the error', () => {
  const read = [{ uid: ONE, k: 1, identity: ONE, selector: '#1', record_id: 'NC_1', record_length: 100, detectedTopology: 'circular' }];
  const { panel, synthetic } = panelOf({ selector: `${TWO}:3`, read });
  assert.equal(synthetic(), 'two.gb #3 (not found)');
  assert.equal(panel.circularRecordPresentationError.value,
    'No record in the current input matches the record selector \'two.gb #3\'.');
});

test('a Record Order row of an unread record is named like the Record select', () => {
  const { panel } = panelOf({ selector: '' });
  assert.equal(panel.getCircularRecordOrderLabel(`${TWO}:2`), 'two.gb #2');
  assert.equal(panelOf({ selector: '', named: { [`${TWO}:2`]: byId('BGC0000709') } }).panel.getCircularRecordOrderLabel(`${TWO}:2`),
    'two.gb: BGC0000709');
  assert.equal(panel.getCircularRecordOrderLabel('#4'), '#4');
});

// D-98: the rows a loaded Session's Record Order shows before its cards are read tell the records apart as the
// read rows do. Each case names what 6ac46afe9 (before #1770) showed.
test('unread records whose IDs repeat are told apart by file or by #k', () => {
  // Two one-record cards with one ID (6ac46afe9: `a.gbk`, `b.gbk`).
  const pair = panelOf({
    selector: '', files: { [ONE]: 'a.gbk', [TWO]: 'b.gbk' }, named: { [ONE]: 'BGC0000708', [TWO]: 'BGC0000708' }
  }).panel;
  assert.deepEqual([ONE, TWO].map(pair.getCircularRecordOrderLabel), ['a.gbk: BGC0000708', 'b.gbk: BGC0000708']);
  // One file whose records repeat an ID (6ac46afe9: `dup.gb #1`, `dup.gb #2`): with one card a Record Order row is
  // `#k (<record ID>)`, as a read row is; the Record select names it with the `#k` its options show.
  const dupNamed = { [`${ONE}:1`]: 'X', [`${ONE}:2`]: 'X', [`${ONE}:3`]: byId('Y') };
  const dup = panelOf({ selector: `${ONE}:2`, files: { [ONE]: 'dup.gb' }, named: dupNamed });
  assert.deepEqual([1, 2, 3].map((k) => dup.panel.getCircularRecordOrderLabel(`${ONE}:${k}`)), ['#1 (X)', '#2 (X)', '#3 (Y)']);
  assert.equal(dup.synthetic(), 'X [#2] (not inspected)');
  // Only one repeat drawn (single-record mode, or the other repeat OFF): the Result names one `X`, which the Web's
  // request selected by `#2`, so the name keeps its `#2` (6ac46afe9: `dup.gb #2 (not inspected)`).
  const oneRepeat = { [`${ONE}:2`]: 'X' };
  assert.equal(panelOf({ selector: `${ONE}:2`, files: { [ONE]: 'dup.gb' }, named: oneRepeat }).synthetic(), 'X [#2] (not inspected)');
  const twoCards = panelOf({ selector: '', files: { [ONE]: 'dup.gb', [TWO]: 'b.gb' }, named: oneRepeat }).panel;
  assert.equal(twoCards.getCircularRecordOrderLabel(`${ONE}:2`), 'dup.gb: X [#2]');
  // An unsafe record ID takes its `#k` too.
  const unsafe = panelOf({ selector: `${ONE}:1`, files: { [ONE]: 'u.gb' }, named: { [`${ONE}:1`]: '#2', [`${ONE}:2`]: byId('Z') } });
  assert.equal(unsafe.synthetic(), '#2 [#1] (not inspected)');
  assert.equal(unsafe.panel.getCircularRecordOrderLabel(`${ONE}:2`), '#2 (Z)');
});

// Review 1 of #1778: a Record Order row names its record from the same entry before and after the records are read.
test('a Record Order row names a record alike before and after its card is read', () => {
  const files = { [ONE]: 'dup.gb', [TWO]: 'b.gb' };
  const named = { [`${ONE}:1`]: 'X', [`${ONE}:2`]: 'X', [TWO]: 'Z' };
  const rows = [`${ONE}:1`, `${ONE}:2`, TWO];
  const unread = panelOf({ selector: '', files, named }).panel;
  const read = panelOf({ selector: '', files, named, read: [readRecord(ONE, 1, 'X'), readRecord(ONE, 2, 'X'), readRecord(TWO, 1, 'Z', TWO)] }).panel;
  assert.deepEqual(rows.map(read.getCircularRecordOrderLabel), ['dup.gb: X [#1]', 'dup.gb: X [#2]', 'b.gb: Z']);
  assert.deepEqual(rows.map(unread.getCircularRecordOrderLabel), rows.map(read.getCircularRecordOrderLabel));
  // One card: `#k (<record ID>)` read and unread.
  const one = { [ONE]: 'dup.gb' };
  const oneRead = panelOf({ selector: '', files: one, read: [readRecord(ONE, 1, 'X'), readRecord(ONE, 2, 'X'), readRecord(ONE, 3, 'Y')] }).panel;
  const oneUnread = panelOf({ selector: '', files: one, named: { [`${ONE}:1`]: 'X', [`${ONE}:2`]: 'X', [`${ONE}:3`]: byId('Y') } }).panel;
  const oneRows = [1, 2, 3].map((k) => `${ONE}:${k}`);
  assert.deepEqual(oneRows.map(oneRead.getCircularRecordOrderLabel), ['#1 (X)', '#2 (X)', '#3 (Y)']);
  assert.deepEqual(oneRows.map(oneUnread.getCircularRecordOrderLabel), oneRows.map(oneRead.getCircularRecordOrderLabel));
});

test('the three one-record cards of a Web Session keep their file names (Sakai)', () => {
  const files = { a: 'NC_002695.gb', b: 'NC_002128.gb', c: 'NC_002127.gb' };
  const named = { a: 'NC_002695.2', b: 'NC_002128.1', c: 'NC_002127.1' };
  // 6ac46afe9: `NC_002695.gb`, `NC_002128.gb`, `NC_002127.gb`.
  const { panel } = panelOf({ selector: '', files, named });
  assert.deepEqual(['a', 'b', 'c'].map(panel.getCircularRecordOrderLabel),
    ['NC_002695.gb: NC_002695.2', 'NC_002128.gb: NC_002128.1', 'NC_002127.gb: NC_002127.1']);
});
