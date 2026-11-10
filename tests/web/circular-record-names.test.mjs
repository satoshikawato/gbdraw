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
 * named by `{ recordId, idSuffices: true }`, as is a record of a file its request drew whole.
 * @param {{ selector: string, named?: Record<string, string | { recordId: string, idSuffices: boolean }>, read?: object[] | null, files?: Record<string, string>,
 *   recordsOff?: string[], canvas?: boolean }} input
 */
const panelOf = ({ selector, named = {}, read = null, files = { [ONE]: 'one.gb', [TWO]: 'two.gb' }, recordsOff = [], canvas = false }) => {
  const sources = Object.entries(files).map(([uid, name]) => ({ uid, gb: file(name), gff: null, fasta: null }));
  const drawing = {
    id: 'circular', mode: 'circular', inputType: { value: 'gb' }, circularSources: sources,
    form: { circular_record_selector: selector, multi_record_canvas: canvas }, adv: { circular_grouping_intent: '' },
    recordsOff, circularRecordList: { value: read ?? [] },
    circularRecordDiscovery: read
      ? { status: 'ready', error: '', inputType: 'gb', sources: sources.map((card) => ({ uid: card.uid, primaryFile: card.gb, pairedFile: null, error: '' })) }
      : { status: 'idle', error: '', inputType: '', sources: [] }
  };
  const state = { activeDrawing: () => drawing, semanticFileWatchersSuppressed: { value: false } };
  const panel = createCircularRecordPresentation({
    state, mode: { value: 'circular' }, semanticMutationAvailable: { value: true }, sessionOperationAvailability: () => null,
    shownDrawingOfMode: () => drawing, refreshCircularRecordOrder: () => undefined,
    shownRecordIds: (shown) => new Map(shown.id === drawing.id ? Object.entries(named).map(([key, name]) => [
      key, typeof name === 'string' ? { recordId: name, idSuffices: false } : name
    ]) : [])
  });
  return { panel, synthetic: () => panel.circularRecordPresentationOptions.value.find((option) => option.synthetic)?.label };
};

const byId = (recordId) => ({ recordId, idSuffices: true });
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
  // A Session the CLI or the Python API wrote draws each card's whole file, so every ID of the file is named
  // (`idSuffices`): a unique ID shows no `#k` before the card is read, and a repeated one keeps it (review 2 of #1778).
  const cliRows = [`${ONE}:1`, `${ONE}:2`, TWO];
  for (const [first, second] of [['A', 'B'], ['X', 'X']]) {
    const cliNamed = { [`${ONE}:1`]: byId(first), [`${ONE}:2`]: byId(second), [TWO]: 'Z' };
    const cliRead = panelOf({ selector: '', files: { [ONE]: 'm.gb', [TWO]: 'b.gb' }, named: cliNamed,
      read: [readRecord(ONE, 1, first), readRecord(ONE, 2, second), readRecord(TWO, 1, 'Z', TWO)] }).panel;
    const cliUnread = panelOf({ selector: '', files: { [ONE]: 'm.gb', [TWO]: 'b.gb' }, named: cliNamed }).panel;
    assert.deepEqual(cliRows.map(cliUnread.getCircularRecordOrderLabel), cliRows.map(cliRead.getCircularRecordOrderLabel));
  }
  assert.deepEqual(cliRows.map(panelOf({ selector: '', files: { [ONE]: 'm.gb', [TWO]: 'b.gb' },
    named: { [`${ONE}:1`]: byId('A'), [`${ONE}:2`]: byId('B'), [TWO]: 'Z' } }).panel.getCircularRecordOrderLabel),
  ['m.gb: A', 'm.gb: B', 'b.gb: Z']);
});

test('the three one-record cards of a Web Session keep their file names (Sakai)', () => {
  const files = { a: 'NC_002695.gb', b: 'NC_002128.gb', c: 'NC_002127.gb' };
  const named = { a: 'NC_002695.2', b: 'NC_002128.1', c: 'NC_002127.1' };
  // 6ac46afe9: `NC_002695.gb`, `NC_002128.gb`, `NC_002127.gb`.
  const { panel } = panelOf({ selector: '', files, named });
  assert.deepEqual(['a', 'b', 'c'].map(panel.getCircularRecordOrderLabel),
    ['NC_002695.gb: NC_002695.2', 'NC_002128.gb: NC_002128.1', 'NC_002127.gb: NC_002127.1']);
});

// OV-3340: the status line counts the records the next Generate draws (`circularDrawingRecordSet`): the record the
// Record select names alone, else the ON records. The automatic option describes the all-records choice, so it keeps
// the ON count.
test('the drawn count is the record set the drawing draws; the automatic option keeps the ON count (OV-3340)', () => {
  const read = [readRecord(TWO, 1, 'TESTA'), readRecord(TWO, 2, 'TESTB')];
  const files = { [TWO]: 'two.gb' };
  const automatic = (panel) => panel.circularRecordPresentationOptions.value[0].label;
  const chosen = panelOf({ selector: `${TWO}:2`, read, files }).panel;
  assert.equal(chosen.circularDrawnRecordCount.value, 1, 'TESTB alone');
  assert.equal(automatic(chosen), 'All records (separate diagrams)');
  const all = panelOf({ selector: '', read, files });
  assert.equal(all.panel.circularDrawnRecordCount.value, 2);
  const off = panelOf({ selector: '', read, files, recordsOff: [`${TWO}:1`] }).panel;
  assert.equal(off.circularDrawnRecordCount.value, 1);
  assert.equal(automatic(off), '1 of 2 records (separate diagrams)');
  const chosenOff = panelOf({ selector: `${TWO}:2`, read, files, recordsOff: [`${TWO}:1`] }).panel;
  assert.equal(chosenOff.circularDrawnRecordCount.value, 1);
  assert.equal(automatic(chosenOff), '1 of 2 records (separate diagrams)', 'the all-records choice draws the ON records');
  const canvas = panelOf({ selector: `${TWO}:2`, read, files, canvas: true }).panel;
  assert.equal(canvas.circularDrawnRecordCount.value, 2, 'a Multi-Record Canvas draws every ON record');
});
