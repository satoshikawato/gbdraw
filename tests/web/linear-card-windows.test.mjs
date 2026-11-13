// The windows of a Linear record card (PD-OI-124): the Region block and the Windows list edit one owner
// (services/linear-card-windows.js); an End typed before the Start without "Across the origin" is stored in order with
// Reverse complement ticked (OV-2884); Generate fills an empty window with the whole record and refuses a window
// marked across the origin whose Start is not after its End; the Record Layout has one cell per window; Add window and
// Remove go through the record rows' one transition, which drops the stored comparisons whose unit indexes they shift
// (F7).
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    // The stub runs an immediate watcher once, at creation: the readers ask Python for the plan of the state they see.
    watch: (source, callback, options) => { if (options?.immediate) callback(source()); }
  }
};

const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  addCardWindow, cardWindowList, editCardWindow, orientLinearCard, windowSpanOf,
  removeCardWindow, UNKNOWN_CARD_RECORD, cardRecordFacts, catalogCardRecord
} = await import('../../gbdraw/web/js/services/linear-card-windows.js');
const {
  linearCardRequestedCircular, linearCardWindowPlanRecords, linearDrawingUnitPlaces, linearDrawingWindowsOf, linearWindowsPayload
} = await import('../../gbdraw/web/js/services/session-request.js');
const {
  layoutCellKey, linearLayoutCells, linearLayoutPlaces, moveLinearLayoutCell, placeNewLayoutWindow, placeRecordAtRowEnd, setLinearRecordRow
} = await import('../../gbdraw/web/js/services/linear-record-layout.js');
const { createLinearRecordEdits } = await import('../../gbdraw/web/js/app/linear-source-edits.js');
const { linearComparisonEdgeKey } = await import('../../gbdraw/web/js/services/linear-comparisons.js');
const { windowBases } = await import('../../gbdraw/web/js/services/record-options.js');
const { RECORD_SETTINGS_DEFAULTS, linearCardHasSettings } = await import('../../gbdraw/web/js/services/record-draw-selection.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const { pythonCardWindowPlan, pythonLinearUnitPlan } = await import('./helpers/linear-unit-plan.mjs');
const { createLinearCardWindowReaders } = await import('../../gbdraw/web/js/app/linear-card-window-readers.js');
const { linearUnitPlanPayload } = await import('../../gbdraw/web/js/services/linear-units.js');

const span = (start, end, reverseComplement = false, column = null) => ({ start, end, reverseComplement, column, acrossOrigin: false });
const crop = (card) => [card.region_start, card.region_end, card.region_reverse, card.windows.length];

test('an End typed before the Start is stored in order with Reverse complement ticked (OV-2884)', () => {
  const card = createLinearSeq({ uid: 'a', region_start: 500, region_end: null });
  assert.equal(editCardWindow(card, 0, { end: 100 }), true);
  assert.deepEqual(crop(card), [100, 500, true, 0]);
  // The figure stays the reversed span 100..500; the request writes it as before.
  assert.deepEqual(linearWindowsPayload(card, { arranged: false, ordinal: 1, record: UNKNOWN_CARD_RECORD }), [{ start: 100, end: 500, reverseComplement: true, column: null }]);
  // In a window of the list the same holds.
  const listed = createLinearSeq({ uid: 'b', windows: [span(1, 50), span(200, 300)] });
  editCardWindow(listed, 1, { start: 900 });
  assert.deepEqual(listed.windows[1], span(300, 900, true));
});

test('with "Across the origin" ticked a Start after the End stays a window across the origin', () => {
  const card = createLinearSeq({ uid: 'a' });
  editCardWindow(card, 0, { acrossOrigin: true });
  editCardWindow(card, 0, { start: 16000 });
  editCardWindow(card, 0, { end: 100 });
  assert.deepEqual(crop(card), [null, null, false, 1]);
  assert.deepEqual(card.windows, [{ start: 16000, end: 100, reverseComplement: false, column: null, acrossOrigin: true }]);
  assert.equal(windowBases(16000, 100, 16569), 670);
  assert.equal(windowBases(16000, 100, null), null);
  // A window across the origin that New drawing from region or Find in other files writes arrives with the box ticked.
  const made = createLinearSeq({ uid: 'm', windows: [{ start: 16000, end: 100, reverseComplement: false, column: null }] });
  assert.equal(made.windows[0].acrossOrigin, true);
  assert.deepEqual(linearWindowsPayload(card, { arranged: false, ordinal: 1, record: UNKNOWN_CARD_RECORD }), [{ start: 16000, end: 100, reverseComplement: false, column: null }]);
  // Clearing the box keeps the fields as they are, now the reversed span.
  editCardWindow(card, 0, { acrossOrigin: false });
  assert.deepEqual(crop(card), [100, 16000, true, 0]);
});

test('one rule reads Start, End and "Across the origin" for the card, the request and Find in other files', async () => {
  const { referenceWindowOf } = await import('../../gbdraw/web/js/services/find-in-other-files.js');
  const cases = [[8100, 1501, false], [8100, 1501, true], [100, 1501, true], [1501, 8100, false]];
  for (const [start, end, acrossOrigin] of cases) {
    const spanned = windowSpanOf(start, end, acrossOrigin);
    const row = referenceWindowOf({ start: String(start), end: String(end), acrossOrigin });
    assert.deepEqual('window' in row ? row.window : null, spanned);
    const card = createLinearSeq({ uid: 'a', windows: [{ start, end, reverseComplement: false, column: null, acrossOrigin }, span(1, 10)] });
    if (spanned === null) {
      assert.throws(() => linearWindowsPayload(card, { arranged: false, ordinal: 1, record: UNKNOWN_CARD_RECORD }));
    } else {
      assert.deepEqual(linearWindowsPayload(card, { arranged: false, ordinal: 1, record: UNKNOWN_CARD_RECORD })[0],
        { ...spanned, column: null });
    }
  }
});

test('Generate refuses a window across the origin whose Start is not after its End, naming the card and the record', () => {
  const card = createLinearSeq({ uid: 'a', windows: [{ ...span(100, 900), acrossOrigin: true }] });
  const record = { length: 16569, id: 'NC_012920.1', count: 1 };
  assert.throws(() => linearWindowsPayload(card, { arranged: false, ordinal: 2, record }), (error) => {
    const message = normalizeUserFacingError(error, { stage: 'request-validation' }).summary;
    assert.match(message, /Sequence 2\./);
    // It names the record (S2).
    assert.match(message, /Record: NC_012920\.1\./);
    assert.match(message, /Across the origin is ticked: give a Start greater than the End, or clear the box\./);
    return true;
  });
});

test('Add window appends an empty window; Remove leaves the last window; one window goes back to the crop fields', () => {
  const card = createLinearSeq({ uid: 'a', region_start: 10, region_end: 90 });
  addCardWindow(card, null);
  assert.deepEqual(card.windows, [span(10, 90), { start: null, end: null, reverseComplement: false, column: null, acrossOrigin: false }]);
  assert.deepEqual(crop(card).slice(0, 3), [null, null, false]);
  assert.equal(removeCardWindow(card, 0), true);
  // The empty window is now the card's Region block: the whole record.
  assert.deepEqual(crop(card), [null, null, false, 0]);
  assert.equal(removeCardWindow(card, 0), false);
  assert.deepEqual(cardWindowList(card), [{ start: null, end: null, reverseComplement: false, column: null, acrossOrigin: false }]);
});

test('an empty window draws the whole record for every reader of the card\'s windows (B1)', () => {
  const card = createLinearSeq({ uid: 'a', windows: [span(1, 50), span(null, null)] });
  const catalog = { records: [{ sourceIndex: 0, recordLength: 3000, recordId: 'NC_1', sourceRecordCount: 1 },
    { sourceIndex: 1, recordLength: 'x', sourceRecordCount: 3 }] };
  // `fileCount`: the records of the card's file, by which the request names the record (`cliRecordSelector`).
  assert.deepEqual(catalogCardRecord(catalog, 0), { length: 3000, id: 'NC_1', count: 1, fileCount: 1 });
  assert.deepEqual(catalogCardRecord(catalog, 1), { length: null, id: '', count: 1, fileCount: 3 });
  assert.deepEqual(catalogCardRecord(null, 0), UNKNOWN_CARD_RECORD);
  // The request, the drawn-unit plan, the LOSAT estimate and the extraction all read the windows here.
  const record = catalogCardRecord(catalog, 0);
  assert.deepEqual(linearWindowsPayload(card, { arranged: false, ordinal: 2, record }).map(({ start, end }) => [start, end]),
    [[1, 50], [1, 3000]]);
  const drawing = { linearRecordLayoutEnabled: { value: false }, linearSeqs: [createLinearSeq({ uid: 'off' }), card] };
  assert.deepEqual(linearDrawingWindowsOf(drawing, () => record)(card).map(({ start, end }) => [start, end]), [[1, 50], [1, 3000]]);
  // The card is untouched: it still shows an empty window.
  assert.deepEqual(card.windows[1], span(null, null));
  // Before the record is read, the refusal names the card, counted among all cards as the card list shows them.
  assert.throws(() => linearDrawingWindowsOf(drawing, () => UNKNOWN_CARD_RECORD)(card), (error) => {
    assert.equal(error.context.inputOrdinal, 2);
    assert.match(normalizeUserFacingError(error, { stage: 'request-validation' }).summary, /Sequence 2\./);
    return true;
  });
});

test('the extraction reads the topology with an empty window as the request writes it (B1)', () => {
  const card = createLinearSeq({ uid: 'a', windows: [span(1, 50), span(null, null)] });
  const drawing = { linearRecordLayoutEnabled: { value: false }, linearSeqs: [card], recordDisplayDrafts: [] };
  const windows = linearDrawingWindowsOf(drawing, () => ({ length: 3000, id: 'NC_1', count: 1 }))(card);
  assert.equal(linearCardRequestedCircular({ drawing, card, index: 0, recordDisplayRows: [], windows }), null);
});

/**
 * The card views and Record Layout cells of a Linear drawing of `cards`, each record `length` bp and `topology`, as
 * the readers show them once Python (the real `plan_card_windows_json`) has answered.
 */
const shownWindows = async (cards, { length, topology = 'linear', entries = null, recordsOf = null }) => {
  const record = { recordLength: length, detectedTopology: topology, recordId: 'NC_1' };
  const drawnRecords = (_drawing, card) => recordsOf?.(card) ?? [record];
  const drawing = {
    mode: 'linear', linearSeqs: cards, drawnLinearSeqs: { value: cards }, recordDisplayDrafts: [],
    linearRecordLayoutEnabled: { value: entries !== null }, linearRecordRows: entries || []
  };
  let asked = Promise.resolve();
  const readers = createLinearCardWindowReaders({
    state: { activeDrawing: () => drawing }, shownDrawingOfMode: () => drawing,
    linearSeqRecord: (_drawing, card) => (drawnRecords(_drawing, card).length === 1 ? drawnRecords(_drawing, card)[0] : null),
    linearCardRecord: (_drawing, card) => cardRecordFacts(drawnRecords(_drawing, card), null),
    planCardWindows: (payload) => { asked = pythonCardWindowPlan(payload); return asked; },
    cardWindowPlanRecords: linearCardWindowPlanRecords, recordDisplayRows: () => []
  });
  await asked;
  await new Promise((resolve) => setImmediate(resolve));
  return {
    cards: cards.map((card) => readers.linearCardWindowView(card).windows.map(({ number, gapText }) => [number, gapText])),
    cells: readers.linearLayoutCellRows().map(({ uid, window, number }) => [uid, window, number])
  };
};

test('the card shows the drawing\'s window numbers and omitted lengths, from Python (OV-2992)', async () => {
  const between = (bp, first, second) => `${bp} bp of the record lie between windows ${first} and ${second} and are left out.`;
  // Two windows of a linear record with 200 bp between them.
  assert.deepEqual((await shownWindows([createLinearSeq({ uid: 'a', windows: [span(1, 100), span(301, 400)] })],
    { length: 10000 })).cards, [[[1, between(200, 1, 2)], [2, '']]]);
  // Overlapping windows of a circular record: the drawing writes no length, nor does the card.
  assert.deepEqual((await shownWindows([createLinearSeq({ uid: 'a', windows: [span(1, 600), span(500, 900)] })],
    { length: 1000, topology: 'circular' })).cards, [[[1, ''], [2, '']]]);
  // A third window covers the interval between windows 1 and 2.
  assert.deepEqual((await shownWindows([createLinearSeq({ uid: 'a', windows: [span(1, 100), span(301, 400), span(50, 350)] })],
    { length: 10000 })).cards, [[[1, ''], [2, ''], [3, '']]]);
  // Up on the third window puts it second in the row: the numbers follow the row and no window pair is next to its
  // record neighbour, so the drawing writes no length.
  const cards = [createLinearSeq({ uid: 'a', windows: [span(1, 100), span(201, 300), span(401, 500)] })];
  const entries = moveLinearLayoutCell(cards, [{ uid: 'a', row: 1 }], { uid: 'a', window: 2 }, -1);
  assert.ok(entries);
  const moved = await shownWindows(cards, { length: 10000, entries });
  assert.deepEqual(moved.cards, [[[1, ''], [3, ''], [2, '']]]);
  assert.deepEqual(moved.cells, [['a', 0, 1], ['a', 2, 2], ['a', 1, 3]]);
});

test('a card drawn whole beside a windowed card does not hide its numbers, and each record it draws is counted', async () => {
  // m reads a file of two records, no record chosen: it draws both whole, with no length the card knows. a is drawn
  // in two windows. On the base the card waited for m's length and never showed a number.
  const m = createLinearSeq({ uid: 'm' });
  const a = createLinearSeq({ uid: 'a', windows: [span(1, 100), span(301, 400)] });
  const recordsOf = (card) => (card.uid === 'm' ? [{ recordLength: 500, recordId: 'm1' }, { recordLength: 700, recordId: 'm2' }] : null);
  // Rows of their own: a's windows are 1 and 2 of its row.
  const between = (first, second) => `200 bp of the record lie between windows ${first} and ${second} and are left out.`;
  assert.deepEqual((await shownWindows([m, a], { length: 10000, recordsOf })).cards[1], [[1, between(1, 2)], [2, '']]);
  // One row: m's two records come first, so the drawing numbers a's windows 3 and 4.
  const entries = [{ uid: 'm', row: 1 }, { uid: 'a', row: 1 }];
  assert.deepEqual((await shownWindows([m, a], { length: 10000, recordsOf, entries })).cards[1], [[3, between(3, 4)], [4, '']]);
});

test('until Python can plan the windows the card shows no number and no omitted length', async () => {
  // An End without a Start is a window Generate refuses: no plan is asked for.
  const card = createLinearSeq({ uid: 'a', windows: [span(1, 100), span(null, 400)] });
  assert.deepEqual((await shownWindows([card], { length: 10000 })).cards, [[[null, ''], [null, '']]]);
});

// The places the request writes, read by the real Python codec: tests/web/linear-card-windows-request.test.mjs.

/** The cards and Record Layout of an arranged drawing: a in row 1 drawn in two windows, b in row 1, c in row 2. */
const arrangedLayout = () => ({
  cards: [
    createLinearSeq({ uid: 'a', windows: [span(1, 100), span(301, 400)] }),
    createLinearSeq({ uid: 'b' }),
    createLinearSeq({ uid: 'c' })
  ],
  entries: [{ uid: 'a', row: 1 }, { uid: 'b', row: 1 }, { uid: 'c', row: 2 }]
});
const cellsOf = (cells) => cells.map(({ uid, window, row, position }) => [uid, window, row, position]);

test('the Record Layout has one cell per window; Up and Down move a cell within its row', () => {
  const { cards, entries } = arrangedLayout();
  assert.deepEqual(cellsOf(linearLayoutCells(cards, entries)), [['a', 0, 1, 1], ['a', 1, 1, 2], ['b', null, 1, 3], ['c', null, 2, 1]]);
  // No cell has a place yet: the request writes none.
  assert.equal(linearLayoutPlaces(linearLayoutCells(cards, entries)), null);
  const moved = moveLinearLayoutCell(cards, entries, { uid: 'b', window: null }, -1);
  assert.ok(moved);
  assert.deepEqual(cellsOf(linearLayoutCells(cards, moved)), [['a', 0, 1, 1], ['b', null, 1, 2], ['a', 1, 1, 3], ['c', null, 2, 1]]);
  // b's whole-card cell passed a window of a: the cards swap in the list, so their "#n" swap as on dev (PD-OI-018).
  assert.deepEqual(cards.map(({ uid }) => uid), ['b', 'a', 'c']);
  assert.deepEqual(cards[1].windows.map(({ column }) => column), [1, 3]);
  // Row 1's cells now have their places and c in row 2 has none: the request writes a place for every cell, as
  // Python reads a place for every cell or none.
  const places = linearLayoutPlaces(linearLayoutCells(cards, moved));
  assert.ok(places);
  assert.deepEqual([['a', 0], ['b', null], ['a', 1], ['c', null]].map(([uid, window]) => places.get(layoutCellKey(uid, window))), [1, 2, 3, 1]);
  assert.equal(moveLinearLayoutCell(cards, moved, { uid: 'a', window: 0 }, -1), null);
});

test('Up and Down depend on the cell: a window moves only its place, a whole card swaps with its neighbour\'s card', () => {
  const { cards, entries } = arrangedLayout();
  // A window cell passing the whole card b: places only, the list stays.
  const moved = moveLinearLayoutCell(cards, entries, { uid: 'a', window: 1 }, 1);
  assert.ok(moved);
  assert.deepEqual(cards.map(({ uid }) => uid), ['a', 'b', 'c']);
  assert.deepEqual(cellsOf(linearLayoutCells(cards, moved)), [['a', 0, 1, 1], ['b', null, 1, 2], ['a', 1, 1, 3], ['c', null, 2, 1]]);
  // The whole card b passing window 0 of a: the cards swap, the places follow the cells.
  const back = moveLinearLayoutCell(cards, moved, { uid: 'b', window: null }, -1);
  assert.ok(back);
  assert.deepEqual(cards.map(({ uid }) => uid), ['b', 'a', 'c']);
  assert.deepEqual(cellsOf(linearLayoutCells(cards, back)), [['b', null, 1, 1], ['a', 0, 1, 2], ['a', 1, 1, 3], ['c', null, 2, 1]]);
});

test('a Row change puts the record\'s windows at the end of the new row', () => {
  // a (two windows) in row 1, b in row 2, no places: Row 2 for a puts a's windows after b.
  const cards = [createLinearSeq({ uid: 'a', windows: [span(1, 100), span(301, 400)] }), createLinearSeq({ uid: 'b' })];
  const entries = [{ uid: 'a', row: 1 }, { uid: 'b', row: 2 }];
  setLinearRecordRow(entries, 'a', 2);
  const placed = placeRecordAtRowEnd(cards, entries, 'a');
  assert.ok(placed);
  assert.deepEqual(cellsOf(linearLayoutCells(cards, placed)), [['b', null, 2, 1], ['a', 0, 2, 2], ['a', 1, 2, 3]]);
  // A card drawn whole keeps the rule it had: its place goes with its row.
  assert.equal(placeRecordAtRowEnd(cards, placed, 'b'), null);
});

test('a cell without a place among cells with one gets its position in the request', () => {
  const { cards, entries } = arrangedLayout();
  cards[0].windows[1].column = 1;
  const places = linearLayoutPlaces(linearLayoutCells(cards, entries));
  assert.ok(places);
  assert.deepEqual([layoutCellKey('a', 1), layoutCellKey('a', 0), layoutCellKey('b', null)].map((key) => places.get(key)), [1, 2, 3]);
});

test('Add window places the new window right after the card\'s last cell', () => {
  const { cards, entries } = arrangedLayout();
  const placed = placeNewLayoutWindow(cards, entries, 'a');
  assert.ok(placed);
  assert.deepEqual([placed.column, placed.lastColumn], [3, 2]);
  addCardWindow(cards[0], placed.column);
  assert.deepEqual(cellsOf(linearLayoutCells(cards, placed.entries)),
    [['a', 0, 1, 1], ['a', 1, 1, 2], ['a', 2, 1, 3], ['b', null, 1, 4], ['c', null, 2, 1]]);
});

/** A Linear record edits owner over one drawing whose rows `cards` are, with the transition `applied` records. */
const recordEditsOf = (cards, places = new Map()) => {
  const drawing = {
    linearSeqs: cards, linearRecordRows: [], linearRecordLayoutEnabled: { value: false },
    linearComparisonPlan: { mode: 'none', defaultSource: 'losat', edges: [] }
  };
  const applied = [];
  const steps = [];
  const cleared = [];
  const edits = createLinearRecordEdits({
    sessionOperationAvailability: () => null,
    shownDrawingOfMode: () => drawing,
    runUndoable: (label, action) => { steps.push(label); action(); return true; },
    clearAlignmentForMutation: (reason) => { cleared.push(String(reason)); },
    linearRecordsFor: () => [],
    linearRecordStatusFor: () => ({}),
    shownUnitPlaces: () => places,
    rebindAnnotationRecords: () => {},
    applyLinearSeqMutation: (target, items, options) => {
      applied.push({ items: items.map((item) => ({ uid: item.uid, windows: item.windows.length })), options });
      target.linearSeqs = items;
    }
  });
  return { drawing, applied, steps, edits, cleared };
};

test('Reverse complement turns a window, so it keeps the alignment; Start, End and Across the origin crop (OV-3026)', () => {
  const cards = [createLinearSeq({ uid: 'a', region_start: 1, region_end: 100, region_reverse: true }),
    createLinearSeq({ uid: 'b', windows: [span(1, 50), span(200, 300)] })];
  const { edits, cleared } = recordEditsOf(cards);
  // The Region block's Reverse and a listed window's Reverse are the orientation the alignment records.
  assert.equal(edits.setLinearCardWindow(cards[0], 0, 'reverseComplement', false), true);
  assert.equal(edits.setLinearCardWindow(cards[1], 1, 'reverseComplement', true), true);
  assert.deepEqual([cards[0].region_reverse, cards[1].windows[1].reverseComplement], [false, true]);
  assert.deepEqual(cleared, []);
  // An edit that changes nothing clears nothing.
  assert.equal(edits.setLinearCardWindow(cards[0], 0, 'start', 1), false);
  assert.deepEqual(cleared, []);
  // The span is the crop: each of its fields clears the alignment once.
  edits.setLinearCardWindow(cards[0], 0, 'end', 90);
  edits.setLinearCardWindow(cards[1], 0, 'start', 5);
  edits.setLinearCardWindow(cards[1], 1, 'acrossOrigin', true);
  assert.deepEqual(cleared, ['record crop changed.', 'record crop changed.', 'record crop changed.']);
  // An End typed before the Start also ticks Reverse complement (OV-2884): still a crop.
  edits.setLinearCardWindow(cards[1], 0, 'end', 2);
  assert.deepEqual(cards[1].windows[0], span(2, 5, true));
  assert.equal(cleared.length, 4);
});

// OV-2742: Remove window moves the names typed on window pairs with their windows, inside its own History step.
test('Remove window moves the window-pair names of the plan with their windows', () => {
  const cards = [createLinearSeq({ uid: 'a', windows: [span(1, 100), span(101, 200), span(201, 300)] }), createLinearSeq({ uid: 'b' })];
  // The drawn Result's units, as Python keys them (`window_unit_key`).
  const places = new Map([['a~w1', { recordKey: 'a', window: 1 }], ['a~w2', { recordKey: 'a', window: 2 }],
    ['a~w3', { recordKey: 'a', window: 3 }], ['b', { recordKey: 'b', window: null }]]);
  const { drawing, steps, edits } = recordEditsOf(cards, places);
  drawing.linearComparisonPlan.edges.push({ id: 'e1', queryUid: 'a', subjectUid: 'b', windowLosatFilenames: {
    [linearComparisonEdgeKey('a~w2', 'b')]: 'second.tsv', [linearComparisonEdgeKey('a~w1', 'b')]: 'first.tsv' } });
  edits.removeLinearCardWindow(cards[0], 0);
  assert.deepEqual(steps.map(String), ['Remove window']);
  assert.deepEqual(drawing.linearComparisonPlan.edges[0].windowLosatFilenames, { [linearComparisonEdgeKey('a~w1', 'b')]: 'second.tsv' });
  // Down to one window, the record is drawn whole: its window pairs and their names go.
  edits.removeLinearCardWindow(drawing.linearSeqs[0], 0);
  edits.removeLinearCardWindow(drawing.linearSeqs[0], 0);
  assert.equal(Object.hasOwn(drawing.linearComparisonPlan.edges[0], 'windowLosatFilenames'), false);
});

// OV-2742: with no drawn Result the Web cannot place a window-pair name (Python forms the unit keys), so Remove window
// drops the edited record's names rather than leave one on another window's span; another record pair keeps its own.
test('Remove window with no drawn Result drops the edited record\'s window-pair names', () => {
  const cards = [createLinearSeq({ uid: 'a', windows: [span(1, 100), span(101, 200), span(201, 300)] }),
    createLinearSeq({ uid: 'b' }), createLinearSeq({ uid: 'c', windows: [span(1, 50), span(51, 100)] })];
  const { drawing, steps, edits } = recordEditsOf(cards);
  const othersNames = { [linearComparisonEdgeKey('c~w1', 'b')]: 'other.tsv' };
  drawing.linearComparisonPlan.edges.push(
    { id: 'e1', queryUid: 'a', subjectUid: 'b', windowLosatFilenames: { [linearComparisonEdgeKey('a~w3', 'b')]: 'third.tsv' } },
    { id: 'e2', queryUid: 'c', subjectUid: 'b', windowLosatFilenames: { ...othersNames } });
  edits.removeLinearCardWindow(cards[0], 0);
  assert.deepEqual(steps.map(String), ['Remove window']);
  assert.equal(Object.hasOwn(drawing.linearComparisonPlan.edges[0], 'windowLosatFilenames'), false);
  assert.deepEqual(drawing.linearComparisonPlan.edges[1].windowLosatFilenames, othersNames);
});

test('Add window and Remove are one History step each through the record rows\' transition (F7)', () => {
  const cards = [createLinearSeq({ uid: 'a', region_start: 1, region_end: 100 }), createLinearSeq({ uid: 'b' })];
  const { drawing, applied, steps, edits } = recordEditsOf(cards);
  edits.addLinearCardWindow(cards[0]);
  edits.removeLinearCardWindow(drawing.linearSeqs[0], 0);
  assert.deepEqual(steps.map(String), ['Add window', 'Remove window']);
  assert.deepEqual(applied.map(({ items }) => items), [[{ uid: 'a', windows: 2 }, { uid: 'b', windows: 0 }], [{ uid: 'a', windows: 0 }, { uid: 'b', windows: 0 }]]);
  // The transition keeps a LOSAT cache-info row's file; that it empties the stored comparisons (their unit indexes
  // shift) is checked through the app's own transition in linear-card-windows-state.test.mjs.
  assert.equal(applied.every(({ options }) => options.preserveLosatCacheInfo === true), true);
  // Editing a Start keeps the unit plan, so it stays a field edit.
  edits.setLinearCardWindow(drawing.linearSeqs[0], 0, 'start', 5);
  assert.equal(applied.length, 2);
  assert.equal(drawing.linearSeqs[0].region_start, 5);
});

test('the card owns every write of its windows: orientation, Reset Settings and Delete settings values', () => {
  const card = createLinearSeq({ uid: 'a', windows: [span(1, 100), span(301, 400)] });
  assert.equal(orientLinearCard(card, true), true);
  assert.deepEqual(card.windows.map(({ reverseComplement }) => reverseComplement), [true, true]);
  assert.equal(orientLinearCard(card, true), false);
  // Delete settings resets the window fields to the whole record, and a card without a window has no settings
  // (OV-2993: `[] !== []` made every card report settings).
  const cleared = createLinearSeq({ ...card, ...RECORD_SETTINGS_DEFAULTS });
  assert.deepEqual(crop(cleared), [null, null, false, 0]);
  assert.equal(linearCardHasSettings(createLinearSeq({ uid: 'b' })), false);
  assert.equal(linearCardHasSettings(card), true);
});

test('Up and Down on a window cell change only the places, so the unit indexes stored comparisons hold stay (F7)', () => {
  const { cards, entries } = arrangedLayout();
  const unitsOf = () => pythonLinearUnitPlan(linearUnitPlanPayload(cards, { edges: [] }, linearDrawingUnitPlaces(
    { linearRecordLayoutEnabled: { value: true }, linearRecordRows: entries, linearSeqs: cards }, () => UNKNOWN_CARD_RECORD)));
  const before = unitsOf();
  moveLinearLayoutCell(cards, entries, { uid: 'a', window: 1 }, -1);
  assert.notDeepEqual(cards[0].windows.map(({ column }) => column), [null, null]);
  assert.deepEqual(unitsOf(), before);
});
