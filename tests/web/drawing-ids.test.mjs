// The drawing ID rule (B2-D4): `allocateDrawingId` in
// services/drawing-view.js and Python's
// `gbdraw.session_drawings.allocate_drawing_id`
// (tests/test_session_drawings_vectors.py) read the same `ids` vectors of
// tests/fixtures/sessions/drawings-vectors.json. Its `names` vectors pin a
// new drawing's name (`defaultDrawingName`; Python's
// `allocate_session_drawing`) and the name rule (`drawingNamesProblem`;
// Python's `_identity_problem`).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { allocateDrawingId, defaultDrawingName, drawingNamesProblem } from '../../gbdraw/web/js/services/drawing-view.js';

/** @typedef {{ mode: 'circular' | 'linear', taken: string[], id: string }} DrawingIdVector */

/** @typedef {DrawingIdVector & { name: string }} DrawingNameVector */
/** @typedef {{ names: string[], valid: boolean }} DrawingIdentityVector */

/** @type {{ schema: number, ids: DrawingIdVector[], names: { defaults: DrawingNameVector[], identities: DrawingIdentityVector[] } }} */
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

for (const { mode, taken, id, name } of VECTORS.names.defaults) {
  test(`a new ${mode} drawing beside [${taken.join(', ')}] is ${id}, named ${name}`, () => {
    const drawingId = allocateDrawingId(mode, taken);
    assert.equal(drawingId, id);
    assert.equal(defaultDrawingName(mode, drawingId), name);
  });
}

for (const { names, valid } of VECTORS.names.identities) {
  test(`the names [${names.join(' | ')}] are ${valid ? '' : 'not '}valid`, () => {
    assert.equal(drawingNamesProblem(names) === null, valid);
  });
}
