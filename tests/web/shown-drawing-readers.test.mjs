// The readers of the shown drawing (app/shown-drawing-readers.js): Wrap Rows as the shown Result drew it compares
// the records the next Generate draws with Python's input summary by the windows their cards write.
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
