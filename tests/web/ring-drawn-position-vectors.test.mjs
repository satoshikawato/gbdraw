// The Web's side of the shared ring drawn-position vectors (tests/fixtures/ring_drawn_position_vectors.json, written by
// tools/generate_ring_drawn_position_vectors.py from Python's `_drawn_ring_colors`). `ringDrawnPositions` gives each
// ring Python draws its place in Python's drawn order, outermost first, and the rings Python does not draw the
// positions after them, in series order (the Web's own rule: such a ring draws nothing).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { ringDrawnPositions } from '../../gbdraw/web/js/services/conservation-series.js';

const vectors = JSON.parse(readFileSync(new URL('../fixtures/ring_drawn_position_vectors.json', import.meta.url), 'utf8'));

test('the Web gives each ring the drawn position Python draws it at', () => {
  assert.ok(vectors.cases.length > 0);
  for (const { id, count, slots, expected } of vectors.cases) {
    const positions = ringDrawnPositions(count, slots);
    expected.drawnPlaces.forEach((place, position) => assert.equal(positions[place], position, `${id}: ring ${place + 1}`));
    const undrawn = Array.from({ length: count }, (_, place) => place).filter((place) => !expected.drawnPlaces.includes(place));
    assert.deepEqual(undrawn.map((place) => positions[place]),
      undrawn.map((_, index) => expected.drawnPlaces.length + index), `${id}: the rings Python does not draw follow`);
  }
});
