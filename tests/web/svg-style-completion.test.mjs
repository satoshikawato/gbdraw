import assert from 'node:assert/strict';
import { test } from 'node:test';
import { fixture } from './helpers/svg-style-fixture.mjs';

test('depth OFF and ON publish parent visibility at the style owner', () => {
  const { actions, state } = fixture({ default: '#d3d3d3' }, []);
  for (const enabled of [false, true]) {
    state.form.show_depth = enabled;
    actions.applyTrackVisibility();
    assert.equal(JSON.parse(state.results.value[0].content).depth.display, enabled ? undefined : 'none');
  }
});
// OV-108 (PD-OI-086): each mode has its own drawing, so the palette styles a
// Result with the skew slot colors of the drawing of that Result's mode.
test('a palette keeps the skew slot colors of the drawing of the displayed Result (OV-108)', () => {
  const { actions, state, skewFills } = fixture({ default: '#d3d3d3', skew_high: '#abcdef', skew_low: '#fedcba' }, []);
  const circular = { ...state.drawings.circular, adv: {
    circular_track_slots_enabled: true,
    circular_track_slots: [{ id: 'gc_skew', renderer: 'dinucleotide_skew', params: { positive_color: '#112233', negative_color: '#445566' } }]
  } };
  const linear = { ...state.drawings.circular, adv: { linear_track_slots_enabled: false, linear_track_slots: [] } };
  Object.defineProperties(state, {
    drawings: { value: Object.freeze({ circular, linear }), configurable: true },
    activeDrawing: { value: () => (state.mode.value === 'linear' ? linear : circular), configurable: true }
  });
  state.generatedMode = { value: 'circular' };
  state.mode.value = 'linear';
  actions.applyPaletteToSvg();
  assert.deepEqual(skewFills, ['#112233', '#445566']);
});
