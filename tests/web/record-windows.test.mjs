// Windows of one record in one Linear row (PR-W, D-01): the Web's slot, request writer, stored-request reader,
// axis entries, search spans, record selectors, and the card's read-only window list.
import assert from 'node:assert/strict';
import test from 'node:test';

import { CANONICAL_REQUEST_SCHEMA } from '../../gbdraw/web/js/contracts/schema.generated.js';
import { createLinearSeq } from '../../gbdraw/web/js/services/drawing/inputs.js';
import {
  promoteCanonicalRenderRequestToCurrent,
  requestRecordPlan
} from '../../gbdraw/web/js/services/session-request.js';
import { linearRecordWindowRows } from '../../gbdraw/web/js/services/record-options.js';
import { readRegionAxisRecords } from '../../gbdraw/web/js/services/region-axis.js';
import { searchSpans } from '../../gbdraw/web/js/services/sequence-search.js';
import { resultRecordBindings } from '../../gbdraw/web/js/app/annotations/record-selector.js';

const WINDOWS = [
  { start: 100, end: 200, reverseComplement: false, column: 3 },
  { start: 900, end: 950, reverseComplement: true, column: 1 }
];

test('a slot keeps the windows of a record drawn in several and reads one window as its own crop', () => {
  const several = createLinearSeq({ uid: 'a', windows: WINDOWS, region_start: 5, region_end: 9, region_reverse: true });
  // The slot marks each window's "Across the origin" box; a stored window reads it from its span (D-45).
  assert.deepEqual(several.windows, WINDOWS.map((window) => ({ ...window, acrossOrigin: false })));
  assert.deepEqual([several.region_start, several.region_end, several.region_reverse], [null, null, false]);
  const one = createLinearSeq({ uid: 'b', windows: [{ start: 7, end: 70, reverseComplement: true, column: null }] });
  assert.deepEqual([one.region_start, one.region_end, one.region_reverse, one.windows], [7, 70, true, []]);
  const crop = createLinearSeq({ uid: 'c', region_start: 3, region_end: 30 });
  assert.deepEqual([crop.region_start, crop.region_end, crop.windows], [3, 30, []]);
});

/** A Linear drawing as the record plan reads it. @param {boolean} arranged */
const linearDrawing = (arranged) => /** @type {any} */ ({
  mode: 'linear',
  inputType: { value: 'genbank' },
  linearRecordRows: arranged ? [{ uid: 'a', row: 1 }] : [],
  linearRecordLayoutEnabled: { value: arranged }
});
const GB = { name: 'a.gb', size: 3, type: 'text/plain', data: 'YWJj', encoding: 'base64' };

test('the request writes a card crop as one window and the windows of a record drawn in several (schema 10)', () => {
  const crop = requestRecordPlan({
    drawing: linearDrawing(false),
    filesData: /** @type {any} */ ({ linearSeqs: [{ uid: 'a', gb: GB, region_record_id: 'NC_1', region_start: 50, region_end: 10 }] })
  });
  assert.ok(crop);
  assert.deepEqual(crop[0].selector, { kind: 'recordId', value: 'NC_1' });
  assert.deepEqual(crop[0].windows, [{ start: 10, end: 50, reverseComplement: true, column: null }]);
  assert.equal(crop[0].presentation.reverseComplement, false);
  assert.equal(Object.hasOwn(crop[0], 'region'), false);

  for (const arranged of [false, true]) {
    const several = requestRecordPlan({
      drawing: linearDrawing(arranged),
      filesData: /** @type {any} */ ({ linearSeqs: [{ uid: 'a', gb: GB, region_record_id: 'NC_1', windows: WINDOWS }] })
    });
    assert.ok(several);
    assert.deepEqual(several[0].windows, WINDOWS.map((window) => ({ ...window, column: arranged ? window.column : null })));
    assert.equal(several[0].presentation.gridColumn, null);
    assert.equal(several[0].cardinality, 'exactly_one');
  }
});

test('a stored schema-9 region reads as one window, its selector as the record selector', () => {
  const promoted = promoteCanonicalRenderRequestToCurrent({
    schema: 9,
    mode: 'linear',
    records: [
      { recordKey: 'a', selector: null, region: { selector: { kind: 'recordId', value: 'NC_1' }, start: 5, end: 9, reverseComplement: true } },
      { recordKey: 'b', selector: { kind: 'recordIndex', index: 1 }, region: null }
    ],
    diagramOptions: { featureOverrides: [{ keep: true }] }
  });
  assert.equal(promoted.schema, CANONICAL_REQUEST_SCHEMA);
  assert.deepEqual(promoted.records.map(({ selector, windows }) => ({ selector, windows })), [
    { selector: { kind: 'recordId', value: 'NC_1' }, windows: [{ start: 5, end: 9, reverseComplement: true, column: null }] },
    { selector: { kind: 'recordIndex', index: 1 }, windows: [] }
  ]);
  assert.equal(promoted.records.some((record) => Object.hasOwn(record, 'region')), false);
  // Schema 9 kept its per-feature edits.
  assert.deepEqual(promoted.diagramOptions.featureOverrides, [{ keep: true }]);
});

test("the card's window list follows the windows' columns", () => {
  assert.deepEqual(linearRecordWindowRows(WINDOWS, null), [
    { number: 1, span: '900–950', length: '51 bp', reverseComplement: true },
    { number: 2, span: '100–200', length: '101 bp', reverseComplement: false }
  ]);
});

test('a window across the origin of a circular record counts the bases it draws, once the record length is known', () => {
  const across = [{ start: 16000, end: 100, reverseComplement: false, column: null }];
  assert.deepEqual(linearRecordWindowRows(across, 16569), [
    { number: 1, span: '16,000–100', length: '670 bp', reverseComplement: false }
  ]);
  assert.equal(linearRecordWindowRows(across, null)[0].length, '');
});

const AXIS = {
  kind: 'linear',
  rows: [{ xStartPx: 0, xEndPx: 10, yPx: 0, bandTopPx: -1, bandBottomPx: 1, displayStart: 0, displayEnd: 10 }],
  source: { anchorBase: 1, recordBase: 1, step: 1, sourceLength: 1000, drawnLength: 10, isCircular: false }
};

test('each window of a record is its own axis record, found by the record key and the window index', () => {
  const geometry = {
    records: [
      { resultIndex: 0, recordIndex: 0, recordKey: 'a', windowIndex: 0, recordId: 'NC_1', axis: AXIS },
      { resultIndex: 0, recordIndex: 1, recordKey: 'a', windowIndex: 1, recordId: 'NC_1', axis: AXIS },
      { resultIndex: 0, recordIndex: 2, recordKey: 'b', recordId: 'NC_2', axis: AXIS },
      { resultIndex: 0, recordIndex: 3, recordKey: 'c', windowIndex: -1, recordId: 'NC_3', axis: AXIS },
      { resultIndex: 0, recordIndex: 4, recordKey: 'd', recordId: 'NC_2', axis: AXIS }
    ]
  };
  const records = readRegionAxisRecords(geometry, 0);
  assert.deepEqual(records.map(({ recordKey, windowIndex }) => [recordKey, windowIndex]), [['a', 0], ['a', 1], ['b', null], ['d', null]]);
  // A search reads each record once (D-11); a `#n` selector counts records, not windows (I-6).
  assert.deepEqual(searchSpans(records, [{ reach: 30 }]).map(({ recordKey }) => recordKey), ['a', 'b', 'd']);
  // Two records named NC_2 are told apart by number: b is the second record, d the fourth (c, whose axis entry is
  // malformed, is still a record the Result draws).
  const bindings = resultRecordBindings(records, []);
  assert.deepEqual([bindings.selectorOf('b'), bindings.selectorOf('d')],
    [{ kind: 'recordIndex', index: 1 }, { kind: 'recordIndex', index: 3 }]);
});
