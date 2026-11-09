// The readers of the shown drawing (app/shown-drawing-readers.js): Wrap Rows as the shown Result drew it compares
// the records the next Generate draws with Python's input summary by the windows their cards write, and the shown
// Result names its records by record ID.
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    computed: (definition) => ({ get value() { return definition(); } })
  }
};

const { createShownDrawingReaders } = await import('../../gbdraw/web/js/app/shown-drawing-readers.js');
const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { linearWindowsPayload } = await import('../../gbdraw/web/js/services/session-request.js');

/** A shown Linear drawing of one card with a wrapped Result, its input summary naming `summaryRegion`. */
const readersFor = (seq, summaryRegion) => {
  const drawing = {
    id: 'd1', mode: 'linear', form: { linear_wrap_length: 500 },
    drawnLinearSeqs: { value: [seq] }, linearSeqs: [seq], linearRecordLayoutEnabled: { value: false }
  };
  const request = { mode: 'linear', diagramOptions: { configOverrides: { 'canvas.linear.wrap_length': 500 } }, records: [] };
  const state = {
    activeDrawing: () => drawing,
    generatedDrawingId: { value: 'd1' },
    results: { value: [{ name: 'r' }] },
    selectedResultIndex: { value: 0 },
    inputSummaries: { value: [{
      resultIndex: 0, resultName: 'r', lines: [], warnings: [],
      records: [{ sourceLength: 16569, region: summaryRegion }],
      linearWrap: { rowLength: 500, rowCount: 2, requested: 500, recordRows: [2], drawnAxisLength: 500 }
    }] }
  };
  return createShownDrawingReaders({
    state, shownDrawingOfMode: () => drawing,
    linearRecordsFor: () => [{ recordLength: 16569 }], linearRecordStatusFor: () => 'ready',
    refreshLinearRecords: () => ({ busy: true }), committedLinearRecordCards: () => null,
    getCommittedCanonicalRenderRequest: () => request,
    windowsOf: (shown) => (card) => linearWindowsPayload(card, { arranged: Boolean(shown.linearRecordLayoutEnabled.value) })
  });
};

test('a card with a window across the origin is the record the wrapped Result drew', () => {
  const seq = createLinearSeq({ windows: [{ start: 16000, end: 100, reverseComplement: false, column: null }] });
  const { shownLinearWrapReading } = readersFor(seq, { start: 16000, end: 100 });
  assert.equal(shownLinearWrapReading.value.kind, 'rows');
});

test('a reversed crop typed End before Start still matches the ordered span Python names', () => {
  const seq = createLinearSeq({ region_start: 300, region_end: 50 });
  const { shownLinearWrapReading } = readersFor(seq, { start: 50, end: 300 });
  assert.equal(shownLinearWrapReading.value.kind, 'rows');
});

// OV-2959: the record ID the shown Result names each record by, by record key: the one its run metadata says the
// record was drawn with (a Session 51 Result), else the record ID its committed request selects (an older Result,
// whose run metadata keys no record); nothing while the drawing shows no Result of its own.
test('the shown Result names its records by the record IDs it drew them with, else by its request', () => {
  const axis = {
    kind: 'circular', radiusPx: 390, bandInnerPx: 370, bandOuterPx: 410, startAngleDeg: -90, degreesPerBp: 0.02,
    source: { anchorBase: 1, recordBase: 1, step: 1, sourceLength: 16569, drawnLength: 16569, isCircular: true }
  };
  const drawn = (recordIndex, recordKey, recordId) => ({ resultIndex: 0, recordIndex, recordKey, recordId, axis });
  const readersOf = ({ geometry, records, shown = 'c1' }) => {
    const drawing = { id: 'c1', mode: 'circular' };
    const state = {
      activeDrawing: () => drawing, generatedDrawingId: { value: shown }, results: { value: [{ name: 'r' }] },
      selectedResultIndex: { value: 0 }, inputSummaries: { value: [] }, trackSlotResolvedGeometry: { value: geometry }
    };
    const { shownRecordIds } = createShownDrawingReaders({
      state, shownDrawingOfMode: () => drawing, linearRecordsFor: () => [], linearRecordStatusFor: () => 'ready',
      refreshLinearRecords: () => ({ busy: true }), committedLinearRecordCards: () => null,
      getCommittedCanonicalRenderRequest: () => ({ mode: 'circular', records }), windowsOf: () => () => []
    });
    return Object.fromEntries(shownRecordIds(drawing));
  };
  const records = [
    { recordKey: 'card-a', selector: null },
    { recordKey: 'card-b:2', selector: { kind: 'recordId', value: ' BGC0000709 ' } },
    { recordKey: 'card-b:3', selector: { kind: 'recordIndex', index: 2 } }
  ];
  const geometry = { schema: 1, mode: 'circular', records: [drawn(0, 'card-a', 'NC_012920.1'), drawn(1, 'card-b:3', 'BGC3')] };
  assert.deepEqual(readersOf({ geometry, records }),
    { 'card-a': 'NC_012920.1', 'card-b:2': 'BGC0000709', 'card-b:3': 'BGC3' });
  // A Result saved before Session 51: its run metadata keys no record, so only the request's record IDs name them.
  const older = { schema: 1, mode: 'circular', records: [{ recordIndex: 0, recordId: 'NC_012920.1' }] };
  assert.deepEqual(readersOf({ geometry: older, records }), { 'card-b:2': 'BGC0000709' });
  // Another drawing's Result names none of this drawing's records.
  assert.deepEqual(readersOf({ geometry, records, shown: 'c2' }), {});
});
