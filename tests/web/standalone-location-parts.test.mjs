// The Web export writer and the Python CLI writer both normalize the location
// parts of a feature into the interactive SVG (gbdraw/render/interactive_svg.py,
// services/standalone-interactivity.js). They read the vectors in
// ../fixtures/location-parts-vectors.json
// (tests/test_interactive_svg_cli_format.py is the Python side).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';
import { normalizeLocationParts } from '../../gbdraw/web/js/services/standalone-interactivity.js';

const vectors = JSON.parse(readFileSync(new URL('../fixtures/location-parts-vectors.json', import.meta.url), 'utf8'));

test('location parts normalize to the vectors the Python writer shares', () => {
  for (const { name, input, expected } of vectors.cases) {
    assert.deepEqual(normalizeLocationParts(structuredClone(input)), expected, name);
  }
});
