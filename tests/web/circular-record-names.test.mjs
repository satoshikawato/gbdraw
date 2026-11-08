// The name the Circular record panel shows for a Record select value whose
// record is not read (app/circular-record-panel.js, OV-2949): a record
// identity names its record by the record ID the shown Result's request names
// it by, else by its card's file name (with `#k` for the k-th record of a
// file), never by the identity itself; any other value shows as it is.
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
 * The panel over one Circular drawing whose records are not read (a loaded Session) or are read.
 * @param {{ selector: string, request?: object | null, read?: object[] | null }} input
 */
const panelOf = ({ selector, request = null, read = null }) => {
  const sources = [{ uid: ONE, gb: file('one.gb'), gff: null, fasta: null }, { uid: TWO, gb: file('two.gb'), gff: null, fasta: null }];
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
    shownCommittedRequest: (shown) => (shown.id === drawing.id ? request : null)
  });
  return { panel, synthetic: () => panel.circularRecordPresentationOptions.value.find((option) => option.synthetic)?.label };
};

test('an unread record identity shows the record ID its Result request names, else its file name', () => {
  const request = { mode: 'circular', records: [{ recordKey: `${TWO}:2`, selector: { kind: 'recordId', value: 'BGC0000709' } }] };
  assert.equal(panelOf({ selector: `${TWO}:2`, request }).synthetic(), 'BGC0000709 (not inspected)');
  // A file's only record has no selector in the request: the file names it.
  assert.equal(panelOf({ selector: ONE, request: { mode: 'circular', records: [{ recordKey: ONE, selector: null }] } }).synthetic(),
    'one.gb (not inspected)');
  // No Result: the k-th record of the file.
  assert.equal(panelOf({ selector: `${TWO}:2` }).synthetic(), 'two.gb #2 (not inspected)');
  // A request of another mode names no Circular record.
  assert.equal(panelOf({ selector: `${TWO}:2`, request: { mode: 'linear', records: request.records } }).synthetic(),
    'two.gb #2 (not inspected)');
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
  assert.equal(panel.getCircularRecordOrderLabel('#4'), '#4');
});
