// One circle (PD-OI-019 scenario revision 2): the record layout owner (`circularRecordLayoutOf`), the choice
// that enters and leaves one circle (`chooseRecordLayout`): leaving returns to the layout one circle left with
// the GC content mode and Hide GC Skew it left, and the choice is one History step also with a custom track
// stack (no Hide GC Skew question); Undo and Duplicate keep the one-circle settings.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { modeDrawing } from './helpers/drawing-state.mjs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = { createElement: () => ({ click: () => {} }), querySelector: () => null, documentElement: { setAttribute: () => {} } };
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { duplicateDrawingState } = await import('../../gbdraw/web/js/services/config.js');
const { createDefaultCircularTrackSlots } = await import('../../gbdraw/web/js/services/circular-track-slot-model.js');
const {
  CIRCULAR_RECORD_LAYOUTS, ONE_CIRCLE_LIMITS, chooseRecordLayout, circularRecordLayoutOf, oneCircleLeftOfRequest,
  readOneCircleLeft
} = await import('../../gbdraw/web/js/services/one-circle-layout.js');
const { readGcContentMode, GC_CONTENT_MODES } = await import('../../gbdraw/web/js/mode-profiles.js');
const { MODE_PROFILE_DATA } = await import('../../gbdraw/web/js/mode-profiles.generated.js');
const { defaultSettingsContainers, installStoredSettings } = await import('../../gbdraw/web/js/services/config-field-readers.js');
const { CANONICAL_REQUEST_DESCRIPTION } = await import('../../gbdraw/web/js/contracts/request-shape.generated.js');
const { ONE_CIRCLE_PROFILE } = await import('../../gbdraw/web/js/mode-profiles.js');

const draft = (layout, gc = 'deviation', suppressSkew = false) => ({
  form: { multi_record_canvas: layout !== 'separate', suppress_skew: suppressSkew },
  adv: { circular_record_layout: layout === 'one_circle' ? 'one_circle' : 'grid', gc_content_mode: gc, one_circle_left: null }
});

test('the layouts and limits are Python\'s, read from the generated request contract', () => {
  const fields = CANONICAL_REQUEST_DESCRIPTION.types.CanonicalCircularLayout.fields;
  assert.deepEqual(CIRCULAR_RECORD_LAYOUTS, fields.recordLayout.values);
  assert.deepEqual(ONE_CIRCLE_LIMITS, {
    gapDegrees: { minimum: fields.oneCircleGapDegrees.minimum, maximum: fields.oneCircleGapDegrees.maximum },
    minLength: { minimum: fields.oneCircleMinLength.minimum }
  });
});

test('leaving one circle returns to the layout it left with the GC values it left', () => {
  for (const [layout, gc, hidden] of [['grid', 'deviation', false], ['separate', 'deviation', false], ['grid', 'percent', true]]) {
    const drawing = draft(layout, gc, hidden);
    const skew = chooseRecordLayout(drawing, 'one_circle');
    drawing.form.suppress_skew = skew;
    assert.equal(circularRecordLayoutOf(drawing), 'one_circle');
    assert.equal(drawing.adv.gc_content_mode, ONE_CIRCLE_PROFILE.gcContentMode);
    assert.equal(skew, !ONE_CIRCLE_PROFILE.showSkew);
    // A value the user changes while one circle is drawn is theirs there, and the layout left keeps its own.
    drawing.adv.gc_content_mode = 'deviation';
    drawing.form.suppress_skew = chooseRecordLayout(drawing, 'left');
    assert.equal(circularRecordLayoutOf(drawing), layout, layout);
    assert.deepEqual([drawing.adv.gc_content_mode, drawing.form.suppress_skew], [gc, hidden]);
    assert.equal(drawing.adv.one_circle_left, null);
  }
  // The grid checkbox leaves one circle for the grid and also puts the values back.
  const drawing = draft('separate');
  drawing.form.suppress_skew = chooseRecordLayout(drawing, 'one_circle');
  assert.equal(chooseRecordLayout(drawing, 'grid'), false);
  assert.equal(circularRecordLayoutOf(drawing), 'grid');
  assert.equal(drawing.adv.gc_content_mode, 'deviation');
});

test('a loaded one-circle request left the fresh values where it draws the profile (Python\'s unset case)', () => {
  // The rule covers every value of Python's ONE_CIRCLE_CONFIG_OVERRIDES, read through the generated mode profile
  // (tools/generate_mode_profiles.py --check): a key Python adds there fails here until the Web reads it too.
  assert.deepEqual(Object.keys(MODE_PROFILE_DATA.oneCircleConfigOverrides).sort(), ['canvas.show_skew', 'objects.gc_content.mode']);
  assert.deepEqual(ONE_CIRCLE_PROFILE, {
    showSkew: MODE_PROFILE_DATA.oneCircleConfigOverrides['canvas.show_skew'],
    gcContentMode: MODE_PROFILE_DATA.oneCircleConfigOverrides['objects.gc_content.mode']
  });
  const otherMode = GC_CONTENT_MODES.find((mode) => mode !== ONE_CIRCLE_PROFILE.gcContentMode);
  const profile = { gc_content_mode: ONE_CIRCLE_PROFILE.gcContentMode, suppress_skew: !ONE_CIRCLE_PROFILE.showSkew };
  const typed = { gc_content_mode: otherMode, suppress_skew: ONE_CIRCLE_PROFILE.showSkew };
  for (const fresh of [profile, typed]) {
    // A value equal to the profile's reads as unset: leaving one circle draws the fresh value.
    assert.deepEqual(oneCircleLeftOfRequest(profile, fresh), { layout: 'grid', ...fresh });
    // Another value was typed: leaving one circle keeps it.
    assert.deepEqual(oneCircleLeftOfRequest(typed, fresh), { layout: 'grid', ...typed });
  }
});

test('a stored GC content mode reads the same in the Advanced setting and in what one circle left (OV-3170)', () => {
  assert.deepEqual(GC_CONTENT_MODES, MODE_PROFILE_DATA.gcContentModes);
  const left = (gc) => readOneCircleLeft({ layout: 'grid', gc_content_mode: gc, suppress_skew: false });
  for (const gc of [...GC_CONTENT_MODES, 'PERCENT', 'bogus', '']) {
    const reading = readGcContentMode(gc);
    const stored = left(gc);
    // A mode the setting would not hold is never written back when one circle is left.
    assert.deepEqual(stored, reading.lost ? { value: null, lost: true }
      : { value: { layout: 'grid', gc_content_mode: reading.value, suppress_skew: false }, lost: false }, gc);
  }
  assert.equal(left('bogus').lost, true);
  // Load reads both stored fields with that reader: a stored "bogus" is lost in both places.
  const targets = defaultSettingsContainers('circular');
  const { losses } = installStoredSettings(targets, { adv: {
    gc_content_mode: 'bogus', circular_record_layout: 'one_circle',
    one_circle_left: { layout: 'grid', gc_content_mode: 'bogus', suppress_skew: false }
  } }, 'circular');
  assert.equal(targets.adv.gc_content_mode, readGcContentMode('bogus').value);
  assert.equal(targets.adv.one_circle_left, null);
  assert.deepEqual(losses.map(({ at }) => at).filter((at) => /gc_content_mode|one_circle_left/.test(at)).sort(),
    ['adv.gc_content_mode', 'adv.one_circle_left']);
});

test('choosing one circle with a custom stack is one History step without the Hide GC Skew question; Undo and Duplicate keep the settings', async () => {
  const app = createAppSetup();
  const history = window.__GBDRAW_HISTORY__;
  const drawing = modeDrawing(state, 'circular');
  drawing.adv.circular_track_slots = createDefaultCircularTrackSlots({ nt: 'GC', showGc: true, showSkew: true, showTicks: true });
  drawing.adv.circular_track_slots_enabled = true;
  const skewRows = () => drawing.adv.circular_track_slots.filter((slot) => slot.renderer === 'dinucleotide_skew');
  assert.ok(skewRows().length > 0 && skewRows().every((slot) => slot.enabled !== false));
  const steps = history.getUndoCount();

  await history.runUndoable('One circle', () => app.chooseCircularRecordLayout('one_circle'));
  // Applied at once: no question was left open, and the skew rows are hidden by the layout's mark.
  assert.equal(drawing.form.suppress_skew, true);
  assert.ok(skewRows().every((slot) => slot.enabled === false));
  assert.equal(history.getUndoCount(), steps + 1, 'one action, one History step');
  // The readers that name the layout drawing every record name one circle, not Multi-Record Canvas.
  assert.equal(app.circularRecordLayout.value, 'one_circle');
  drawing.adv.one_circle_gap_degrees = 2.5;
  drawing.adv.one_circle_min_length = 500;
  drawing.adv.one_circle_record_band = false;
  const settings = (target) => ({
    layout: circularRecordLayoutOf(target),
    gap: target.adv.one_circle_gap_degrees,
    min: target.adv.one_circle_min_length,
    band: target.adv.one_circle_record_band,
    left: target.adv.one_circle_left
  });
  const onCircle = settings(drawing);
  assert.deepEqual(settings(duplicateDrawingState(drawing, 'circular-2')), onCircle, 'Duplicate copies them');

  await history.runUndoable('Gap', () => { drawing.adv.one_circle_gap_degrees = 9; });
  await app.undoHistory();
  assert.deepEqual(settings(drawing), onCircle, 'Undo of a later edit keeps the one-circle settings');

  await history.runUndoable('Leave one circle', () => app.chooseCircularRecordLayout('left'));
  assert.equal(circularRecordLayoutOf(drawing), 'grid');
  assert.deepEqual([drawing.adv.gc_content_mode, drawing.form.suppress_skew], ['deviation', false]);
  assert.equal(app.circularRecordLayout.value, 'grid');
  assert.ok(skewRows().every((slot) => slot.enabled !== false), 'the marked skew rows come back');
  await app.undoHistory();
  assert.deepEqual(settings(drawing), onCircle, 'Undo of leaving puts one circle and its settings back');
  assert.equal(drawing.form.suppress_skew, true);
});
