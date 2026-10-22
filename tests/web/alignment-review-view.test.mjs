import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
  ALIGNMENT_REFERENCE_KEY, alignmentOverlayRecords, buildAlignmentReviewView, layoutAlignmentStrips, recordsNeedingChoice
} from '../../gbdraw/web/js/services/alignment-review-view.js';

const span = (left, right) => ({ left, right });
const close = (actual, expected, message) => assert.ok(
  Math.abs(actual - expected) < 1e-9, `${message}: ${actual} !== ${expected}`
);
const centre = ({ x, width }) => x + width / 2;

// Three records in figure units. The reference feature's centre is at 120.
const geometry = () => ({
  referenceCenter: 120,
  records: [
    { recordKey: 'ref', extent: span(100, 300), features: [span(110, 130), span(200, 240)],
      candidates: [{ key: 'ref-gene', span: span(110, 130) }] },
    { recordKey: 'two', extent: span(0, 200), features: [span(20, 60), span(150, 170)],
      candidates: [{ key: 'left', span: span(20, 60) }, { key: 'right', span: span(150, 170) }] },
    { recordKey: 'still', extent: span(60, 160), features: [span(70, 90)], candidates: [] }
  ]
});
const states = (overrides = {}) => [
  { recordKey: 'ref', lineKey: 'ref-gene', mirrorShift: null },
  { recordKey: 'two', lineKey: 'right', mirrorShift: null, ...overrides.two },
  { recordKey: 'still', lineKey: '', mirrorShift: null, ...overrides.still }
];
const rowOf = (layout, recordKey) => layout.rows.find((row) => row.recordKey === recordKey);
const boxOf = (row, key) => row.candidates.find((candidate) => candidate.key === key);

test('the reference feature and each chosen candidate have their centres on the line', () => {
  const layout = layoutAlignmentStrips({ geometry: geometry(), records: states(), width: 350 });
  close(centre(boxOf(rowOf(layout, 'ref'), 'ref-gene')), layout.lineX, 'reference');
  close(centre(boxOf(rowOf(layout, 'two'), 'right')), layout.lineX, 'chosen candidate');
  // The other candidate keeps its distance from the chosen one, at the shared scale.
  const two = rowOf(layout, 'two');
  close(centre(boxOf(two, 'right')) - centre(boxOf(two, 'left')), 120 * layout.scale, 'candidate distance');
});

test('a record without a candidate on the line keeps its present offset from the reference feature', () => {
  const layout = layoutAlignmentStrips({ geometry: geometry(), records: states(), width: 350 });
  // The reference line is 60 units right of this record's left edge today.
  close(layout.lineX - rowOf(layout, 'still').baseline.x, 60 * layout.scale, 'no candidate');
  const needsChoice = layoutAlignmentStrips({
    geometry: geometry(), records: states({ two: { lineKey: '' } }), width: 350
  });
  close(needsChoice.lineX - rowOf(needsChoice, 'two').baseline.x, 120 * needsChoice.scale, 'needs a choice');
});

test('one scale fits every placed extent into the strip width', () => {
  const layout = layoutAlignmentStrips({ geometry: geometry(), records: states(), width: 350 });
  // Left of the line: ref 20, two 160, still 60. Right of it: ref 180, two 40, still 40.
  close(layout.scale, 350 / (160 + 180), 'scale');
  close(layout.lineX, 160 * layout.scale, 'line');
  layout.rows.forEach((row) => {
    assert.ok(row.baseline.x >= -1e-9, row.recordKey);
    assert.ok(row.baseline.x + row.baseline.width <= 350 + 1e-9, row.recordKey);
  });
  close(rowOf(layout, 'two').baseline.x, 0, 'the widest left side starts at 0');
  close(rowOf(layout, 'ref').baseline.x + rowOf(layout, 'ref').baseline.width, 350, 'the widest right side ends at the width');
});

test('a choice of another candidate moves its box to the line and the record with it', () => {
  const layout = layoutAlignmentStrips({
    geometry: geometry(), records: states({ two: { lineKey: 'left' } }), width: 350
  });
  const two = rowOf(layout, 'two');
  close(centre(boxOf(two, 'left')), layout.lineX, 'new choice');
  assert.ok(centre(boxOf(two, 'right')) > layout.lineX, 'the other candidate is right of the line');
  close(layout.lineX - two.baseline.x, 40 * layout.scale, 'the record starts 40 units left of the line');
});

test('a reversed record is mirrored about the axis the projection gives, not about its drawn features', () => {
  // Python reverses the record about its axis centre. The chosen centre is at 160 and moves by
  // -120 (afterCenterX - beforeCenterX), to 40, so the axis is at 100: the middle of the record's
  // extent 0..200, while its features span 20..170 (middle 95).
  const layout = layoutAlignmentStrips({
    geometry: geometry(), records: states({ two: { mirrorShift: -120 } }), width: 350
  });
  const two = rowOf(layout, 'two');
  close(centre(boxOf(two, 'right')), layout.lineX, 'the chosen candidate stays on the line');
  // Mirrored, the other candidate (centre 40 -> 160) is 120 units right of the chosen one.
  close(centre(boxOf(two, 'left')) - layout.lineX, 120 * layout.scale, 'the other candidate changes side');
  // The extent 0..200 maps onto itself, so the record starts 40 units left of the line.
  close(layout.lineX - two.baseline.x, 40 * layout.scale, 'extent after the mirror');
  // A feature 20..60 becomes 140..180.
  const first = two.features[0];
  close(first.x - two.baseline.x, 140 * layout.scale, 'feature after the mirror');
});

test('a reversed reference record keeps its feature on the line', () => {
  const layout = layoutAlignmentStrips({
    geometry: geometry(),
    records: [{ recordKey: 'ref', lineKey: 'ref-gene', mirrorShift: 160 }, ...states().slice(1)],
    width: 350
  });
  const ref = rowOf(layout, 'ref');
  close(centre(boxOf(ref, 'ref-gene')), layout.lineX, 'reference');
  // Axis at 120 + 80 = 200, the middle of 100..300; the feature 200..240 becomes 160..200.
  close(ref.features[1].x - ref.baseline.x, 60 * layout.scale, 'feature after the mirror');
});

test('small boxes stay visible and centred, and a record without geometry has no row', () => {
  const layout = layoutAlignmentStrips({
    geometry: { referenceCenter: 5, records: [
      { recordKey: 'ref', extent: span(0, 1000), features: [span(4, 6)], candidates: [{ key: 'g', span: span(4, 6) }] }
    ] },
    records: [{ recordKey: 'ref', lineKey: 'g', mirrorShift: null }, { recordKey: 'undrawn', lineKey: '', mirrorShift: null }],
    width: 100
  });
  assert.deepEqual(layout.rows.map(({ recordKey }) => recordKey), ['ref']);
  const box = boxOf(layout.rows[0], 'g');
  assert.equal(box.width, 4);
  close(centre(box), layout.lineX, 'widened box');
  assert.ok(layout.rows[0].features[0].width >= 0.6);
});

test('no geometry, no width or an empty extent gives no strips', () => {
  assert.equal(layoutAlignmentStrips({ geometry: null, records: states(), width: 350 }), null);
  assert.equal(layoutAlignmentStrips({ geometry: geometry(), records: states(), width: 0 }), null);
  assert.equal(layoutAlignmentStrips({
    geometry: { referenceCenter: 3, records: [{ recordKey: 'ref', extent: span(3, 3), features: [], candidates: [] }] },
    records: [{ recordKey: 'ref', lineKey: '', mirrorShift: null }], width: 350
  }), null);
});

// The panel's cards: three records below the reference, in figure order.
const reviewRows = () => [
  { recordKey: 'two', candidates: [{ key: 'right' }, { key: 'undrawn' }, { key: 'left' }], choice: null },
  { recordKey: 'still', candidates: [], choice: { kind: 'skip', candidateKey: null } },
  { recordKey: 'one', candidates: [{ key: 'only' }], choice: { kind: 'select', candidateKey: 'only' } }
];
const directions = () => [
  { recordKey: 'two', status: 'ambiguous', reversed: false },
  { recordKey: 'ref', status: 'reference', reversed: false },
  { recordKey: 'still', status: 'skipped', reversed: false },
  { recordKey: 'one', status: 'aligned', reversed: false }
];
const reviewGeometry = () => {
  const drawn = geometry();
  drawn.records[0].candidates[0].key = ALIGNMENT_REFERENCE_KEY;
  drawn.records.push({ recordKey: 'one', extent: span(40, 90), features: [span(50, 60)], candidates: [{ key: 'only', span: span(50, 60) }] });
  return drawn;
};

const sectionsOf = (view) => view.groups.slice(1);
const referenceOf = (view) => view.groups[0].cards[0];

test('the reference card comes first and each record is in the section it had when the review opened', () => {
  const rows = reviewRows();
  const openedWithoutChoice = recordsNeedingChoice(rows);
  assert.deepEqual(openedWithoutChoice, ['two']);
  const view = buildAlignmentReviewView({ rows, directions: directions(), geometry: reviewGeometry(), openedWithoutChoice, width: 350 });
  assert.deepEqual(view.groups.map(({ id, title, cards }) => [id, title, cards.map(({ recordKey }) => recordKey)]), [
    ['reference', '', ['ref']],
    ['choose', 'Choose a gene', ['two']],
    ['single', 'One candidate', ['one']],
    ['none', 'No usable candidate', ['still']]
  ]);
  assert.equal(referenceOf(view).kind, 'reference');
  // The summary counts what the headings count.
  assert.equal(view.summary, '4 records, 1 record awaiting a choice, 1 record with one candidate, 1 record without a usable candidate.');
  const [two] = sectionsOf(view)[0].cards;
  assert.deepEqual({ kind: two.kind, needsChoice: two.needsChoice, skipped: two.skipped, chosenKey: two.chosenKey },
    { kind: 'several', needsChoice: true, skipped: false, chosenKey: '' });

  // A choice removes the amber frame; the card keeps its section and its layout.
  rows[0] = { ...rows[0], choice: { kind: 'select', candidateKey: 'left' } };
  const chosen = buildAlignmentReviewView({ rows, directions: directions(), geometry: reviewGeometry(), openedWithoutChoice, width: 350 });
  assert.deepEqual(chosen.groups.map(({ id }) => id), ['reference', 'choose', 'single', 'none']);
  const [after] = sectionsOf(chosen)[0].cards;
  assert.deepEqual({ needsChoice: after.needsChoice, chosenKey: after.chosenKey }, { needsChoice: false, chosenKey: 'left' });
  close(centre(boxOf(after.strip, 'left')), chosen.lineX, 'the new choice is on the line');
  // A record that opened with a choice and several candidates is under "Several candidates".
  const preselected = buildAlignmentReviewView({ rows, directions: directions(), geometry: reviewGeometry(), openedWithoutChoice: [] });
  assert.deepEqual(sectionsOf(preselected).map(({ id, title, cards }) => [id, title, cards.length]),
    [['several', 'Several candidates', 1], ['single', 'One candidate', 1], ['none', 'No usable candidate', 1]]);
  assert.equal(preselected.stripWidth, 350);
  assert.equal(preselected.summary, '4 records, 1 record with several candidates, 1 record with one candidate, 1 record without a usable candidate.');
});

test('candidates are listed left to right as the figure draws them, undrawn ones last', () => {
  const rows = reviewRows();
  const view = buildAlignmentReviewView({
    rows, directions: directions(), geometry: reviewGeometry(), openedWithoutChoice: recordsNeedingChoice(rows), width: 350
  });
  const [two] = sectionsOf(view)[0].cards;
  assert.deepEqual(two.candidates.map(({ key }) => key), ['left', 'right', 'undrawn']);
  assert.equal(two.candidates[0], rows[0].candidates[2], 'the draft candidate itself');
});

test('the cards carry the strips of one layout, and a reversed record its mirror shift', () => {
  const rows = reviewRows();
  rows[0] = { ...rows[0], choice: { kind: 'select', candidateKey: 'right' } };
  const reversed = directions();
  reversed[0] = { ...reversed[0], status: 'aligned', reversed: true, beforeCenterX: 160, afterCenterX: 40 };
  const view = buildAlignmentReviewView({ rows, directions: reversed, geometry: reviewGeometry(), openedWithoutChoice: [], width: 350 });
  close(centre(boxOf(referenceOf(view).strip, ALIGNMENT_REFERENCE_KEY)), view.lineX, 'reference');
  const [two] = sectionsOf(view)[0].cards;
  close(centre(boxOf(two.strip, 'right')), view.lineX, 'chosen');
  assert.ok(centre(boxOf(two.strip, 'left')) > view.lineX, 'mirrored: the other candidate is right of the line');

  const withoutFigure = buildAlignmentReviewView({ rows, directions: reversed, geometry: null, openedWithoutChoice: [], width: 350 });
  assert.equal(withoutFigure.lineX, null);
  assert.equal(referenceOf(withoutFigure).strip, null);
  assert.deepEqual(sectionsOf(withoutFigure)[0].cards[0].candidates.map(({ key }) => key), ['right', 'undrawn', 'left']);
});

test('the figure overlay marks each row by its choice and names a ring without a number', () => {
  const rows = [
    { recordKey: 'two', recordLabel: 'Record two', choice: null, candidates: [
      { key: 'left', anchor: { recordKey: 'two', biologicalFeatureId: 'f1' }, name: 'geneA', label: 'geneA (P1)', coordinates: '1..300 bp' }] },
    { recordKey: 'one', recordLabel: 'Record one', choice: { kind: 'select', candidateKey: 'only' }, candidates: [
      { key: 'only', anchor: { recordKey: 'one', biologicalFeatureId: 'f2' }, name: 'geneB', label: 'geneB', coordinates: '5..90 bp' }] },
    { recordKey: 'still', recordLabel: 'Record still', choice: { kind: 'skip', candidateKey: null }, candidates: [] }
  ];
  assert.deepEqual(alignmentOverlayRecords(rows), [
    { recordKey: 'two', chosenKey: '', needsChoice: true, skipped: false, candidates: [
      { key: 'left', anchor: rows[0].candidates[0].anchor, name: 'geneA', accessibleName: 'Select geneA (P1), 1..300 bp for Record two' }] },
    { recordKey: 'one', chosenKey: 'only', needsChoice: false, skipped: false, candidates: [
      { key: 'only', anchor: rows[1].candidates[0].anchor, name: 'geneB', accessibleName: 'Select geneB, 5..90 bp for Record one' }] },
    { recordKey: 'still', chosenKey: '', needsChoice: false, skipped: true, candidates: [] }
  ]);
});
