// The drawing ID rule (B2-D4): `allocateDrawingId` in
// services/drawing-view.js and Python's
// `gbdraw.session_drawings.allocate_drawing_id`
// (tests/test_session_drawings_vectors.py) read the same `ids` vectors of
// tests/fixtures/sessions/drawings-vectors.json.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { allocateDrawingId } from '../../gbdraw/web/js/services/drawing-view.js';

/** @typedef {{ mode: 'circular' | 'linear', taken: string[], id: string }} DrawingIdVector */

/** @type {{ schema: number, ids: DrawingIdVector[] }} */
const VECTORS = JSON.parse(
  readFileSync(new URL('../fixtures/sessions/drawings-vectors.json', import.meta.url), 'utf8'),
);

test('the vectors declare schema 1', () => {
  assert.equal(VECTORS.schema, 1);
});

for (const { mode, taken, id } of VECTORS.ids) {
  test(`allocateDrawingId(${mode}, [${taken.join(', ')}]) is ${id}`, () => {
    assert.equal(allocateDrawingId(mode, taken), id);
  });
}
