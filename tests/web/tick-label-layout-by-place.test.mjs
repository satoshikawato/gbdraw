// OV-3189: the label layout of a ticks row that names none follows its place by the features row, and Python owns
// that rule (gbdraw/diagrams/circular/presets.py `_TICK_LABEL_LAYOUT_BY_FEATURE_PLACE`). The Web orients its ticks
// rows by the same table, carried by mode-profiles.generated.js, so a CLI Session that Load reads draws what the CLI
// drew, and the Web keeps no copy of the two layouts.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import test from 'node:test';
import { MODE_PROFILE_DATA } from '../../gbdraw/web/js/mode-profiles.generated.js';
import {
  createDefaultCircularTrackSlots, normalizeCircularTrackSlots, syncSlotsFromAxisIndex
} from '../../gbdraw/web/js/services/circular-track-slot-model.js';

const BY_PLACE = MODE_PROFILE_DATA.circularTrackPresets.tickLabelLayoutByFeaturePlace;
const ROWS = {
  features: { id: 'features', renderer: 'features', params: {} },
  ticks: { id: 'ticks', renderer: 'ticks', params: {} },
  gc_content: { id: 'gc_content', renderer: 'dinucleotide_content', params: {} }
};
/** A stack as the CLI writes it into a Session: no side, no tick_label_layout. */
const loaded = (order) => {
  const slots = normalizeCircularTrackSlots(order.map((id) => structuredClone(ROWS[id])));
  syncSlotsFromAxisIndex(slots, null);
  return slots.find((slot) => slot.renderer === 'ticks')?.params.tick_label_layout;
};

test('the table comes from Python and names two layouts', () => {
  assert.deepEqual({ ...BY_PLACE }, { before: 'label_out_tick_in', after: 'label_in_tick_out' });
});

test('a loaded ticks row without a layout takes the layout of its place by the features row', () => {
  assert.equal(loaded(['ticks', 'features', 'gc_content']), BY_PLACE.before);
  assert.equal(loaded(['ticks', 'gc_content', 'features']), BY_PLACE.before);
  assert.equal(loaded(['features', 'ticks', 'gc_content']), BY_PLACE.after);
  assert.equal(loaded(['gc_content', 'features', 'ticks']), BY_PLACE.after);
});

test('the Web\'s own stack puts its ticks row after the features with the "after" layout', () => {
  const ticks = createDefaultCircularTrackSlots().find((slot) => slot.renderer === 'ticks');
  assert.equal(ticks?.params.tick_label_layout, BY_PLACE.after);
});

test('the slot model holds no copy of the layouts by place', () => {
  const source = readFileSync(new URL('../../gbdraw/web/js/services/circular-track-slot-model.js', import.meta.url), 'utf8');
  assert.match(source, /MODE_PROFILE_DATA\.circularTrackPresets\.tickLabelLayoutByFeaturePlace/);
  assert.doesNotMatch(source, /DEFAULT_INNER_TICK_LABEL_LAYOUT|defaultPresetTickLabelLayout/);
});
