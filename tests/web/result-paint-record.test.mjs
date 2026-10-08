import assert from 'node:assert/strict';
import test from 'node:test';

import { createResultPaintRecord } from '../../gbdraw/web/js/app/result-paint-record.js';

// An editor paint state as app-setup reads it (`currentEditorProjectionState`).
/** @param {string} step @returns {import('../../gbdraw/web/js/app/result-paint-record.js').EditorPaintState} */
const paint = (step) => ({
  colors: [{ CDS: step }, `rules ${step}`],
  visibility: `visibility ${step}`,
  strokes: `strokes ${step}`,
  labels: 'labels',
  legendOrder: ''
});
/** @param {ReturnType<typeof paint>} shown @param {ReturnType<typeof paint>} current */
const changed = (shown, current) => ({
  colors: shown.colors[0] !== current.colors[0] || shown.colors[1] !== current.colors[1],
  visibility: shown.visibility !== current.visibility,
  strokes: shown.strokes !== current.strokes
});
const NOTHING = { colors: false, visibility: false, strokes: false };
const EVERY_PAINT_DOMAIN = { colors: true, visibility: true, strokes: true };

// U2BFIX review #1: Save keeps a batch Result's bytes as stored, so a Result
// not displayed since the last edits is saved without them. A Load (or a
// History restore of an artifact) records every Result but the mounted one as
// showing an unknown paint state, so its first display shows every paint
// domain; a Generate draws every Result with the edits.
test('a Result restored from saved bytes shows every paint domain on its first display', () => {
  const edited = paint('edited');
  const loaded = createResultPaintRecord();
  loaded.commit(['a', 'b'], 'a', edited, { restored: true });
  const shownB = loaded.display('b', edited);
  assert.deepEqual(changed(shownB, edited), EVERY_PAINT_DOMAIN, 'the Result saved without the edits');
  assert.equal(shownB.labels, edited.labels);
  assert.equal(shownB.legendOrder, edited.legendOrder);
  loaded.shown('b', edited, shownB);
  assert.deepEqual(changed(loaded.display('a', edited), edited), NOTHING, 'the Result mounted by the Load showed them');
  assert.deepEqual(changed(loaded.display('b', edited), edited), NOTHING, 'shown once, the Result holds them');

  const generated = createResultPaintRecord();
  generated.commit(['a', 'b'], 'a', edited);
  assert.deepEqual(changed(generated.display('b', edited), edited), NOTHING, 'Generate drew the edits');
});

// U2BFIX review #2: a display whose palette and rules projection declined
// (rules changed meanwhile, or a Session operation) showed neither the fills
// nor the visibility; the Result is recorded behind on both, also after it
// leaves, so its next display shows them.
test('a display whose palette and rules projection declined leaves the Result behind on fills and visibility', () => {
  const before = paint('before');
  const after = paint('after');
  const record = createResultPaintRecord();
  record.commit(['a', 'b'], 'a', before);
  const shownB = record.display('b', after);
  assert.deepEqual(changed(shownB, after), EVERY_PAINT_DOMAIN);
  record.shown('b', after, shownB, { declined: true });
  record.display('a', after);
  assert.deepEqual(changed(record.display('b', after), after), { colors: true, visibility: true, strokes: false },
    'the declined domains are shown on the next display');
  record.shown('b', after, record.display('b', after));
  record.display('a', after);
  assert.deepEqual(changed(record.display('b', after), after), NOTHING);
  assert.equal(record.depart(after), 'b');
  assert.deepEqual(changed(record.display('b', after), after), NOTHING, 'a Result that left after a full display');
});
