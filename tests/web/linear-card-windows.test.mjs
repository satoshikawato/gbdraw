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
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {}, watch: () => {}
  }
};

const { createLinearSeq } = await import('../../gbdraw/web/js/services/drawing/inputs.js');
const {
  addCardWindow, cardWindowGaps, cardWindowList, catalogRecordLength, editCardWindow, originWindowLength,
  removeCardWindow, resolveWholeRecordWindows
} = await import('../../gbdraw/web/js/services/linear-card-windows.js');
const { linearWindowsPayload } = await import('../../gbdraw/web/js/services/session-request.js');
const {
  layoutCellKey, linearLayoutCells, linearLayoutPlaces, moveLinearLayoutCell, placeNewLayoutWindow
} = await import('../../gbdraw/web/js/services/linear-record-layout.js');
const { createLinearRecordEdits } = await import('../../gbdraw/web/js/app/linear-source-edits.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const { pythonLinearUnitPlan } = await import('./helpers/linear-unit-plan.mjs');
const { linearUnitPlanPayload } = await import('../../gbdraw/web/js/services/linear-units.js');

const span = (start, end, reverseComplement = false, column = null) => ({ start, end, reverseComplement, column, acrossOrigin: false });
const crop = (card) => [card.region_start, card.region_end, card.region_reverse, card.windows.length];

test('an End typed before the Start is stored in order with Reverse complement ticked (OV-2884)', () => {
  const card = createLinearSeq({ uid: 'a', region_start: 500, region_end: null });
  assert.equal(editCardWindow(card, 0, { end: 100 }), true);
  assert.deepEqual(crop(card), [100, 500, true, 0]);
  // The figure stays the reversed span 100..500; the request writes it as before.
  assert.deepEqual(linearWindowsPayload(card, { arranged: false }), [{ start: 100, end: 500, reverseComplement: true, column: null }]);
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
  assert.equal(originWindowLength(card.windows[0], 16569), 670);
  assert.equal(originWindowLength(card.windows[0], null), null);
  // A window across the origin that New drawing from region or Find in other files writes arrives with the box ticked.
  const made = createLinearSeq({ uid: 'm', windows: [{ start: 16000, end: 100, reverseComplement: false, column: null }] });
  assert.equal(made.windows[0].acrossOrigin, true);
  assert.deepEqual(linearWindowsPayload(card, { arranged: false }), [{ start: 16000, end: 100, reverseComplement: false, column: null }]);
  // Clearing the box keeps the fields as they are, now the reversed span.
  editCardWindow(card, 0, { acrossOrigin: false });
  assert.deepEqual(crop(card), [100, 16000, true, 0]);
});

test('Generate refuses a window across the origin whose Start is not after its End, naming the card', () => {
  const card = createLinearSeq({ uid: 'a', windows: [{ ...span(100, 900), acrossOrigin: true }] });
  assert.throws(() => linearWindowsPayload(card, { arranged: false, ordinal: 2 }), (error) => {
    const message = normalizeUserFacingError(error, { stage: 'request-validation' }).summary;
    assert.match(message, /Sequence 2\./);
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

test('an empty window draws the whole record at Generate, from the record catalog length', () => {
  const card = createLinearSeq({ uid: 'a', windows: [span(1, 50), span(null, null)] });
  assert.throws(() => linearWindowsPayload(card, { arranged: false }), /region/i);
  const catalog = { records: [{ sourceIndex: 0, recordLength: 3000 }, { sourceIndex: 1, recordLength: 'x' }] };
  assert.equal(catalogRecordLength(catalog, 0), 3000);
  assert.equal(catalogRecordLength(catalog, 1), null);
  assert.equal(catalogRecordLength(null, 0), null);
  const resolved = resolveWholeRecordWindows(card, catalogRecordLength(catalog, 0));
  assert.deepEqual(linearWindowsPayload(resolved, { arranged: false }).map(({ start, end }) => [start, end]), [[1, 50], [1, 3000]]);
  assert.equal(resolveWholeRecordWindows(card, null), card);
  assert.deepEqual(card.windows[1], span(null, null));
});

test('the bases between two windows next to each other in the same direction are counted', () => {
  const linear = { recordLength: 10000, circular: false };
  assert.deepEqual(cardWindowGaps([span(1, 100), span(301, 400)], linear), [{ after: 0, bp: 200 }]);
  // Reversed windows read in record order from the later one back.
  assert.deepEqual(cardWindowGaps([span(301, 400, true), span(1, 100, true)], linear), [{ after: 0, bp: 200 }]);
  // Different directions, overlapping windows, or another window between them: no count.
  assert.deepEqual(cardWindowGaps([span(1, 100), span(301, 400, true)], linear), []);
  assert.deepEqual(cardWindowGaps([span(1, 100), span(50, 400)], linear), []);
  assert.deepEqual(cardWindowGaps([span(1, 100), span(501, 600), span(201, 300)], linear), []);
  assert.deepEqual(cardWindowGaps([span(1, 100), span(301, 400)], { recordLength: null, circular: false }), []);
  // A circular record's interval may run through the origin.
  assert.deepEqual(cardWindowGaps([span(9901, 10000), span(51, 100)], { recordLength: 10000, circular: true }), [{ after: 0, bp: 50 }]);
});

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
  assert.deepEqual(cards[0].windows.map(({ column }) => column), [1, 3]);
  // Row 1's cells now have their places and c in row 2 has none: the request writes a place for every cell, as
  // Python reads a place for every cell or none.
  const places = linearLayoutPlaces(linearLayoutCells(cards, moved));
  assert.ok(places);
  assert.deepEqual([['a', 0], ['b', null], ['a', 1], ['c', null]].map(([uid, window]) => places.get(layoutCellKey(uid, window))), [1, 2, 3, 1]);
  assert.equal(moveLinearLayoutCell(cards, moved, { uid: 'a', window: 0 }, -1), null);
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
const recordEditsOf = (cards) => {
  const drawing = { linearSeqs: cards, linearRecordRows: [], linearRecordLayoutEnabled: { value: false } };
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

test('Add window and Remove are one History step each through the record rows\' transition (F7)', () => {
  const cards = [createLinearSeq({ uid: 'a', region_start: 1, region_end: 100 }), createLinearSeq({ uid: 'b' })];
  const { drawing, applied, steps, edits } = recordEditsOf(cards);
  edits.addLinearCardWindow(cards[0]);
  edits.removeLinearCardWindow(drawing.linearSeqs[0], 0);
  assert.deepEqual(steps.map(String), ['Add window', 'Remove window']);
  assert.deepEqual(applied.map(({ items }) => items), [[{ uid: 'a', windows: 2 }, { uid: 'b', windows: 0 }], [{ uid: 'a', windows: 0 }, { uid: 'b', windows: 0 }]]);
  // The transition drops the stored comparisons (their unit indexes shift); a LOSAT cache-info row keeps its file.
  assert.equal(applied.every(({ options }) => options.preserveLosatCacheInfo === true), true);
  // Editing a Start keeps the unit plan, so it stays a field edit.
  edits.setLinearCardWindow(drawing.linearSeqs[0], 0, 'start', 5);
  assert.equal(applied.length, 2);
  assert.equal(drawing.linearSeqs[0].region_start, 5);
});

test('Up and Down change only the places, so the unit indexes stored comparisons hold stay (F7)', () => {
  const { cards, entries } = arrangedLayout();
  const unitsOf = () => pythonLinearUnitPlan(linearUnitPlanPayload(cards, { edges: [] },
    (card) => linearWindowsPayload(card, { arranged: true })));
  const before = unitsOf();
  moveLinearLayoutCell(cards, entries, { uid: 'a', window: 1 }, -1);
  assert.notDeepEqual(cards[0].windows.map(({ column }) => column), [null, null]);
  assert.deepEqual(unitsOf(), before);
});
