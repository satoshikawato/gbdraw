// PD-OI-044 scenario revision 5: the crop a CLI or Python Session gave a record of a Circular grid or batch is held
// per record, by its card and its name there; it draws while the record is drawn in a grid or batch (the
// app transitions are in drawing-list.test.mjs).
import assert from 'node:assert/strict';
import { test } from 'node:test';

const {
  circularDrawingCroppedRecords, circularRecordEntries, circularRecordRange
} = await import('../../gbdraw/web/js/services/record-options.js');
const { readCircularRecordRanges } = await import('../../gbdraw/web/js/services/config-field-readers.js');

const window = { start: 11, end: 90, reverseComplement: true, column: null, acrossOrigin: false };
const list = [
  { uid: 'one', k: 1, identity: 'one', record_id: 'solo', record_length: 100 },
  { uid: 'two', k: 1, identity: 'two:1', record_id: 'A', record_length: 100 },
  { uid: 'two', k: 2, identity: 'two:2', record_id: 'B', record_length: 100 }
];
const entries = circularRecordEntries(list);
const drawing = (ranges, form = {}) => ({
  form: { circular_record_selector: '', multi_record_canvas: false, circular_region_start: null,
    circular_region_end: null, circular_windows: [], circular_record_ranges: ranges, ...form },
  adv: { circular_grouping_intent: 'batch' }, recordsOff: []
});

test('a stored range names its card, a record and windows with both ends; anything else is left out', () => {
  assert.deepEqual(readCircularRecordRanges([
    { source: 'two', record: 'B', windows: [window] },
    { source: '', record: 'A', windows: [window] },
    { source: 'two', record: 'A', windows: [] },
    { source: 'two', record: 'A', windows: [{ start: 3, end: null }] },
    'junk'
  ]), [{ source: 'two', record: 'B', windows: [window] }]);
  assert.deepEqual(readCircularRecordRanges({ source: 'two' }), []);
});

test('a range finds its record by card and by the name the request gave it there', () => {
  const ranges = readCircularRecordRanges([
    { source: 'one', record: '', windows: [window] }, { source: 'two', record: '#2', windows: [window] }
  ]);
  assert.deepEqual(entries.map((entry) => circularRecordRange(ranges, entries, entry)?.source ?? null), ['one', null, 'two']);
  const byId = readCircularRecordRanges([{ source: 'two', record: 'A', windows: [window] }]);
  assert.deepEqual(entries.map((entry) => circularRecordRange(byId, entries, entry) !== null), [false, true, false]);
  // Another card's record of that name is not it.
  const otherCard = readCircularRecordRanges([{ source: 'one', record: 'A', windows: [window] }]);
  assert.deepEqual(entries.map((entry) => circularRecordRange(otherCard, entries, entry)), [null, null, null]);
});

test('a grid or batch draws its ranged records cropped; a record offered alone draws the fields instead', () => {
  const ranges = readCircularRecordRanges([{ source: 'two', record: 'B', windows: [window] }]);
  assert.deepEqual([...circularDrawingCroppedRecords(list, drawing(ranges))], ['two:2']);
  const alone = drawing(ranges, { circular_record_selector: 'two:2' });
  alone.adv.circular_grouping_intent = 'single';
  assert.deepEqual([...circularDrawingCroppedRecords(list, alone)], []);
});
