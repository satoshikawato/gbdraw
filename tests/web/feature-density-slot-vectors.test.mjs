// The rows Python adds for feature density tracks equal the Web's copy
// (vector:tests/fixtures/feature_density_slot_vectors.json, tools/twins.json).
// Each case runs `insertFeatureDensitySlots` over the listed rows, or the Web's
// default stack builder, with the same track ids and compares the rows and the
// Axis index. tools/generate_feature_density_slot_vectors.py writes the cases;
// tests/test_feature_density_slot_vectors.py checks them against Python.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';

import { featureDensitySlotId, insertFeatureDensitySlots } from '../../gbdraw/web/js/services/feature-density-tracks.js';
import { createDefaultCircularTrackSlots } from '../../gbdraw/web/js/services/circular-track-slot-model.js';
import { createDefaultLinearTrackSlots } from '../../gbdraw/web/js/services/linear-track-slot-model.js';

const { cases: CASES } = JSON.parse(readFileSync(new URL('../fixtures/feature_density_slot_vectors.json', import.meta.url), 'utf8'));

const asRow = (slot) => {
  const row = { id: slot.id, renderer: slot.renderer };
  if (slot.renderer === 'feature_density') {
    row.track = slot.params.track;
    row.side = slot.side;
  }
  return row;
};

const makeSlot = (mode) => (track) => ({
  id: featureDensitySlotId(track),
  renderer: 'feature_density',
  side: mode === 'linear' ? 'below' : null,
  params: { track }
});

const givenSlot = ({ id, renderer, track, side }) => ({ id, renderer, side, params: track === undefined ? {} : { track } });

test('the feature density rows equal the rows Python gives', () => {
  assert.ok(CASES.length > 0);
  for (const { id, mode, kind, options, rows, axisIndex, trackIds, expected } of CASES) {
    if (kind === 'default') {
      const build = mode === 'circular' ? createDefaultCircularTrackSlots : createDefaultLinearTrackSlots;
      const slots = build({ ...options, featureDensityTrackIds: trackIds });
      assert.deepEqual(slots.map(asRow).map(({ id: slotId, renderer, track, side }) => ({ id: slotId, renderer, track, side })),
        expected.rows.map((row) => ({ track: undefined, side: undefined, ...row })), id);
      continue;
    }
    const stack = rows.map(givenSlot);
    const axis = insertFeatureDensitySlots(stack, axisIndex, trackIds, makeSlot(mode));
    assert.deepEqual(stack.map(asRow), expected.rows, id);
    assert.equal(axis, expected.axisIndex, id);
  }
});
