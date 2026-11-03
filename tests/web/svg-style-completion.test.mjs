import assert from 'node:assert/strict';
import { test } from 'node:test';
import { fixture } from './helpers/svg-style-fixture.mjs';
import { modeDrawing, withModeDrawings } from './helpers/drawing-state.mjs';

// OV-108 (PD-OI-086), OV-488: the palette styles a Result with the skew slot
// colors of the request that drew it, whichever mode is shown; the track slots
// of neither drawing reach it before Generate (R1).
test('a palette keeps the skew slot colors of the request that drew the displayed Result (OV-108)', () => {
  const { actions, state, skewFills } = fixture({ default: '#d3d3d3', skew_high: '#abcdef', skew_low: '#fedcba' }, [], { slots: [
    { id: 'gc_skew', renderer: 'dinucleotide_skew', enabled: true, params: { positive_color: '#112233', negative_color: '#445566' } }
  ] });
  const draft = (mode, colors) => ({ ...modeDrawing(state, 'circular'), adv: {
    [`${mode}_track_slots_enabled`]: true,
    [`${mode}_track_slots`]: [{ id: 'gc_skew', renderer: 'dinucleotide_skew', params: colors }]
  } });
  const circular = draft('circular', { positive_color: '#00aa00' });
  const linear = draft('linear', { negative_color: '#0000aa' });
  withModeDrawings(state, { circular, linear });
  state.generatedDrawingId = { value: 'circular' };
  state.mode.value = 'linear';
  actions.applyPaletteToSvg();
  assert.deepEqual(skewFills, ['#112233', '#445566']);
});
