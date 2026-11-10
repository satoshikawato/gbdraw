// The record a Linear card draws has one owner, `linearCardRecord` (app/shown-drawing-readers.js, OV-3227): the
// record list (app/record-settings.js), the File name field's default and Generate's raw LOSAT TSV names
// (app/run-analysis.js, tests/web/run-analysis-simple-path.test.mjs) and the card's windows
// (app/linear-card-window-readers.js) read it there, so the record list never names a record Generate does not draw.
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}, watch: () => {}
  }
};

const { createShownDrawingReaders } = await import('../../gbdraw/web/js/app/shown-drawing-readers.js');
const { createRecordSettingsSources } = await import('../../gbdraw/web/js/app/record-settings.js');
const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const { createDrawingState } = await import('../../gbdraw/web/js/state.js');
const { linearCardWindowPlanRecords, linearDrawingWindowsOf } = await import('../../gbdraw/web/js/services/session-request.js');

test('the record list names each Linear card\'s record as the card draws it, from the one owner', () => {
  const read = (recordId) => ({ recordId, selector: '', recordLength: 8 });
  // Each card's records as the record selector read them; `unread` has no list yet (a loaded Session's card).
  const files = new Map([['one', [read('ONLY')]], ['gone', [read('ONLY')]], ['multi', [read('A'), read('B')]], ['picked', [read('A'), read('B')]]]);
  const drawing = createDrawingState('linear-record-owner', 'linear');
  const cards = [
    createLinearSeq({ uid: 'one' }),
    // A selector that names no record of its one-record file: Generate refuses it ("did not match any record ID"), so
    // the card draws no record and the list shows the record the card names, not the file's.
    createLinearSeq({ uid: 'gone', region_record_id: 'GONE' }),
    // No selector over a file of several records: the card draws every record and names none.
    createLinearSeq({ uid: 'multi' }),
    createLinearSeq({ uid: 'picked', region_record_id: 'B' }),
    createLinearSeq({ uid: 'unread', region_record_id: 'LATER' })
  ];
  const { linearCardRecord } = createShownDrawingReaders({
    state: { activeDrawing: () => drawing, inputSummaries: { value: [] }, selectedResultIndex: { value: 0 }, results: { value: [] } },
    shownDrawingOfMode: () => drawing,
    linearRecordsFor: (_drawing, seq) => files.get(seq.uid) ?? [],
    linearRecordStatusFor: (_drawing, seq) => (files.has(seq.uid) ? 'ready' : 'deferred'),
    refreshLinearRecords: () => ({ busy: true }), committedLinearRecordCards: () => null,
    getCommittedCanonicalRenderRequest: () => null,
    windowsOf: linearDrawingWindowsOf, cardWindowPlanRecords: linearCardWindowPlanRecords,
    planCardWindows: async () => ({ result: null }), recordDisplayRows: () => []
  });
  const { recordSourceOf } = createRecordSettingsSources({
    shownDrawingOfMode: () => drawing,
    linearSourceGroups: { value: [{ uid: 'file', records: cards.map((sequence) => ({ sequence })) }] },
    linearSourceFileName: () => 'cards.gb',
    linearCardRecord,
    linearAnnotationSource: () => ({ sourceKey: 'file', records: [], selector: '' }),
    getCircularRecordDiscoveryState: () => ({ status: 'idle', records: [] }),
    circularAnnotationSourceKey: () => '',
    takeBackRuleRows: () => {}
  });
  const shown = recordSourceOf('linear', 'file')?.records.map(({ key, recordId, length }) => [key, recordId, length]);
  assert.deepEqual(shown, [
    ['one', 'ONLY', 8], ['gone', 'GONE', null], ['multi', '', null], ['picked', 'B', 8], ['unread', 'LATER', null]
  ]);
  // The list says what the owner says of every card that draws one read record.
  for (const card of cards) {
    const record = linearCardRecord(drawing, card);
    if (record.count === 1) assert.equal(shown.find(([key]) => key === card.uid)?.[1], record.id, card.uid);
  }
});
