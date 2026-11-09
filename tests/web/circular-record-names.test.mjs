// The name the Circular record panel shows for a Record select value whose
// record is not read (app/circular-record-panel.js, OV-2949, OV-2959, D-98): a
// record identity that the shown Result names (`shownRecordIds`,
// tests/web/shown-drawing-readers.test.mjs) is named as a read record is (its
// card's file name when the drawing has several cards, its record ID, and its
// `#k` when the ID repeats in its file or is unsafe); one the Result does not
// name by its card's file name (with `#k` for the k-th record of a file), never
// by the identity itself; any other value shows as it is.
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
 * @param {{ selector: string, named?: Record<string, string>, read?: object[] | null, files?: Record<string, string> }} input
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
    shownRecordIds: (shown) => new Map(shown.id === drawing.id ? Object.entries(named) : [])
  });
  return { panel, synthetic: () => panel.circularRecordPresentationOptions.value.find((option) => option.synthetic)?.label };
};

test('an unread record identity shows the record ID its Result names it by, else its file name', () => {
  assert.equal(panelOf({ selector: `${TWO}:2`, named: { [`${TWO}:2`]: 'BGC0000709' } }).synthetic(),
    'two.gb: BGC0000709 (not inspected)');
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
  assert.equal(panelOf({ selector: '', named: { [`${TWO}:2`]: 'BGC0000709' } }).panel.getCircularRecordOrderLabel(`${TWO}:2`),
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
  // One file whose records repeat an ID (6ac46afe9: `dup.gb #1`, `dup.gb #2`): the `#k` the Record select shows.
  const dup = panelOf({ selector: '', files: { [ONE]: 'dup.gb' }, named: { [`${ONE}:1`]: 'X', [`${ONE}:2`]: 'X', [`${ONE}:3`]: 'Y' } }).panel;
  assert.deepEqual([1, 2, 3].map((k) => dup.getCircularRecordOrderLabel(`${ONE}:${k}`)), ['X [#1]', 'X [#2]', 'Y']);
  // An unsafe record ID takes its `#k` too.
  const unsafe = panelOf({ selector: '', files: { [ONE]: 'u.gb' }, named: { [`${ONE}:1`]: '#2', [`${ONE}:2`]: 'Z' } }).panel;
  assert.equal(unsafe.getCircularRecordOrderLabel(`${ONE}:1`), '#2 [#1]');
});

test('the three one-record cards of a Web Session keep their file names (Sakai)', () => {
  const files = { a: 'NC_002695.gb', b: 'NC_002128.gb', c: 'NC_002127.gb' };
  const named = { a: 'NC_002695.2', b: 'NC_002128.1', c: 'NC_002127.1' };
  // 6ac46afe9: `NC_002695.gb`, `NC_002128.gb`, `NC_002127.gb`.
  const { panel } = panelOf({ selector: '', files, named });
  assert.deepEqual(['a', 'b', 'c'].map(panel.getCircularRecordOrderLabel),
    ['NC_002695.gb: NC_002695.2', 'NC_002128.gb: NC_002128.1', 'NC_002127.gb: NC_002127.1']);
});
