// PD-OI-126 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md: the Feature density section's editor (#199,
// app/feature-density-tracks.js): Add, Duplicate, Move and Remove
// find a track by its id; a track's stack row follows it into a stack that draws or that a Session or the user saved,
// and not into a stack no Session saved (Use custom stack builds that one again with every track); Remove asks first
// when another row starts at or extends to the track's row (PD-OI-104); a field writes the request's shape and
// leaves the value checks to Python (R7).
import assert from 'node:assert/strict';
import test from 'node:test';
import { FEATURE_DENSITY_CHOICE_WORDS, createFeatureDensityTrackEditor } from '../../gbdraw/web/js/app/feature-density-tracks.js';
import { createAskingStackPorts } from '../../gbdraw/web/js/app/track-slot-edits.js';
import { createDefaultAdv } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import { createDefaultCircularTrackSlots } from '../../gbdraw/web/js/services/circular-track-slot-model.js';
import { createDefaultLinearTrackSlots } from '../../gbdraw/web/js/services/linear-track-slot-model.js';
import { MODE_PROFILE_DATA } from '../../gbdraw/web/js/mode-profiles.generated.js';

const SLOTS_KEY = { circular: 'circular_track_slots', linear: 'linear_track_slots' };

/** A drawing of `mode` with its defaults: Use custom stack off, the stack no Session saved. */
const drawingOf = (mode) => ({
  mode,
  adv: createDefaultAdv(mode),
  form: { track_type: 'tuckin', linear_track_layout: 'middle' },
  annotationSets: []
});

const setup = (drawing, { answer = 'hide', busy = null } = {}) => {
  const asked = [];
  const steps = [];
  const normalized = [];
  const { askingStackPorts } = createAskingStackPorts({
    askChoice: async (question) => { asked.push(question); return answer; },
    runUndoable: (label, fn) => { steps.push(label); return fn(); },
    busy: () => busy,
    shownDrawing: () => drawing
  });
  const editor = createFeatureDensityTrackEditor({
    shownDrawing: () => drawing, busy: () => busy, askingStackPorts, normalizeStack: (each) => normalized.push(each)
  });
  return { editor, asked, steps, normalized };
};

const densityRows = (drawing) => drawing.adv[SLOTS_KEY[drawing.mode]]
  .filter(({ renderer }) => renderer === 'feature_density').map(({ id, params }) => `${id}:${params.track}`);

for (const mode of ['circular', 'linear']) {
  test(`${mode}: Add starts a track with Python's defaults and the drawing's first feature type`, () => {
    const drawing = drawingOf(mode);
    const { editor } = setup(drawing);
    editor.addFeatureDensityTrack();
    editor.addFeatureDensityTrack();
    const { featureDensityTrack: choices } = MODE_PROFILE_DATA;
    assert.deepEqual(drawing.adv.feature_density_tracks.map(({ id }) => id), ['track', 'track_2']);
    assert.deepEqual(drawing.adv.feature_density_tracks[0], {
      id: 'track',
      query: { featureType: drawing.adv.features[0], qualifier: null, value: null, minLength: null, maxLength: null },
      metric: choices.defaultMetric, style: choices.defaultStyle, window: null, step: null, color: null, caption: null
    });
    // A drawing that draws no feature type counts every type.
    drawing.adv.features = [];
    editor.addFeatureDensityTrack();
    assert.equal(drawing.adv.feature_density_tracks[2].query.featureType, choices.everyFeatureType);
  });

  test(`${mode}: a stack no Session saved gets no row; a drawn or saved stack gets one after the last density row`, () => {
    const drawing = drawingOf(mode);
    const { editor, normalized } = setup(drawing);
    editor.addFeatureDensityTrack();
    assert.deepEqual(densityRows(drawing), [], 'Use custom stack builds the stack no Session saved again with every track');
    assert.equal(normalized.length, 0);
    drawing.adv[`${SLOTS_KEY[mode]}_enabled`] = true;
    editor.addFeatureDensityTrack();
    editor.duplicateFeatureDensityTrack('track');
    assert.deepEqual(densityRows(drawing), ['density_track_2:track_2', 'density_track_3:track_3']);
    const rows = drawing.adv[SLOTS_KEY[mode]].map(({ id }) => id);
    assert.equal(rows[rows.indexOf('features') + 1], 'density_track_2', JSON.stringify(rows));
    assert.equal(normalized.length, 2);
    // With Use custom stack off, a stack the user changed (it differs from every stack no Session saved) gets the row.
    const saved = drawingOf(mode);
    saved.adv[SLOTS_KEY[mode]].pop();
    setup(saved).editor.addFeatureDensityTrack();
    assert.deepEqual(densityRows(saved), ['density_track:track']);
  });

  test(`${mode}: Duplicate copies a track after it under the next copy id; Move changes the order`, () => {
    const drawing = drawingOf(mode);
    const { editor } = setup(drawing);
    editor.addFeatureDensityTrack();
    editor.addFeatureDensityTrack();
    editor.setFeatureDensityTrackField('track', 'value', ' LTR ');
    editor.duplicateFeatureDensityTrack('track');
    assert.deepEqual(drawing.adv.feature_density_tracks.map(({ id }) => id), ['track', 'track_3', 'track_2']);
    assert.equal(drawing.adv.feature_density_tracks[1].query.value, 'LTR');
    assert.notEqual(drawing.adv.feature_density_tracks[1].query, drawing.adv.feature_density_tracks[0].query, 'a copy');
    editor.moveFeatureDensityTrack('track_2', -1);
    editor.moveFeatureDensityTrack('track', -1);
    assert.deepEqual(drawing.adv.feature_density_tracks.map(({ id }) => id), ['track', 'track_2', 'track_3']);
  });

  test(`${mode}: Remove takes the track and its rows out, and asks when another row reaches them (PD-OI-104)`, async () => {
    const drawing = drawingOf(mode);
    drawing.adv[`${SLOTS_KEY[mode]}_enabled`] = true;
    const { editor, asked, steps } = setup(drawing);
    editor.addFeatureDensityTrack();
    editor.addFeatureDensityTrack();
    const axisBefore = drawing.adv[`${SLOTS_KEY[mode]}_axis_index`];
    assert.equal(editor.removeFeatureDensityTrack('track'), undefined, 'nothing asks');
    assert.deepEqual(asked, []);
    assert.deepEqual(drawing.adv.feature_density_tracks.map(({ id }) => id), ['track_2']);
    assert.deepEqual(densityRows(drawing), ['density_track_2:track_2']);
    if (axisBefore !== null) assert.ok(drawing.adv[`${SLOTS_KEY[mode]}_axis_index`] <= axisBefore);
    // A row placed over the track's row asks; the chosen removal is its own step and turns that row off.
    drawing.adv[SLOTS_KEY[mode]].push({ id: 'regions', renderer: 'annotations', enabled: true, side: 'overlay', z: -1,
      params: { set_id: 'review', anchor_slot: 'density_track_2', layer: 'underlay' } });
    drawing.annotationSets = [{ id: 'review' }];
    await editor.removeFeatureDensityTrack('track_2');
    assert.equal(asked.length, 1);
    assert.match(asked[0].title, /Remove feature density track track_2 and turn off/);
    assert.deepEqual(steps, ['Remove feature density track']);
    assert.deepEqual(drawing.adv.feature_density_tracks, []);
    assert.deepEqual(densityRows(drawing), []);
    assert.equal(drawing.adv[SLOTS_KEY[mode]].find(({ id }) => id === 'regions').enabled, false);
  });

  test(`${mode}: Cancel keeps the track and its rows`, async () => {
    const drawing = drawingOf(mode);
    drawing.adv[`${SLOTS_KEY[mode]}_enabled`] = true;
    const { editor, steps } = setup(drawing, { answer: null });
    editor.addFeatureDensityTrack();
    drawing.adv[SLOTS_KEY[mode]].push({ id: 'regions', renderer: 'annotations', enabled: true, side: 'overlay', z: -1,
      params: { set_id: 'review', anchor_slot: 'density_track', layer: 'underlay' } });
    await editor.removeFeatureDensityTrack('track');
    assert.deepEqual(steps, []);
    assert.deepEqual(drawing.adv.feature_density_tracks.map(({ id }) => id), ['track']);
    assert.deepEqual(densityRows(drawing), ['density_track:track']);
  });
}

test('a field writes the request\'s shape and leaves the value checks to Python (R7)', () => {
  const drawing = drawingOf('circular');
  const { editor } = setup(drawing);
  editor.addFeatureDensityTrack();
  const [track] = drawing.adv.feature_density_tracks;
  const write = (field, value) => editor.setFeatureDensityTrackField('track', field, value);
  assert.equal(write('qualifier', ' rpt_type '), true);
  assert.equal(track.query.qualifier, 'rpt_type');
  assert.equal(write('qualifier', '  '), true);
  assert.equal(track.query.qualifier, null, 'a blank is null');
  assert.equal(write('featureType', ' '), true);
  assert.equal(track.query.featureType, '*', 'a blank type is every type');
  assert.equal(write('minLength', '1000'), true);
  assert.equal(track.query.minLength, 1000);
  assert.equal(write('window', '0'), true);
  assert.equal(track.window, 0, 'Python refuses 0 at Generate, in its words');
  assert.equal(write('window', ''), true);
  assert.equal(track.window, null, 'blank: Auto');
  assert.equal(write('step', 'abc'), false, 'a text that is no number changes nothing');
  assert.equal(write('metric', 'coverage'), true);
  assert.equal(write('metric', 'mean'), false, 'only Python\'s choices');
  assert.equal(write('style', 'histogram'), true);
  assert.equal(write('color', '#d9944a'), true);
  assert.equal(write('color', null), true);
  assert.equal(track.color, null, 'Auto: the Depth color');
  assert.equal(write('caption', 'LTR coverage'), true);
  assert.equal(write('caption', 'LTR coverage'), false, 'unchanged');
  assert.deepEqual({ metric: track.metric, style: track.style, caption: track.caption },
    { metric: 'coverage', style: 'histogram', caption: 'LTR coverage' });
  assert.equal(editor.setFeatureDensityTrackField('missing', 'caption', 'x'), false);
});

test('the editor has words for exactly Python\'s metric and style choices', () => {
  const { metrics, styles } = MODE_PROFILE_DATA.featureDensityTrack;
  assert.deepEqual(Object.keys(FEATURE_DENSITY_CHOICE_WORDS).sort(), [...metrics, ...styles].sort());
});

test('a busy Session stops every edit', () => {
  const drawing = drawingOf('linear');
  const { editor } = setup(drawing, { busy: 'Saving' });
  assert.equal(editor.addFeatureDensityTrack(), 'Saving');
  assert.deepEqual(drawing.adv.feature_density_tracks, []);
});

test('the stacks no Session saved include a new drawing\'s default stack in both modes (OV-3276)', () => {
  assert.deepEqual(createDefaultAdv('circular').circular_track_slots, createDefaultCircularTrackSlots());
  assert.deepEqual(createDefaultAdv('linear').linear_track_slots, createDefaultLinearTrackSlots());
});
