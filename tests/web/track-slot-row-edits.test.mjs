// Custom stack row edits through the Circular and Linear editors (TK-04,
// TK-07, TK-09, TK-10, GX-01): a renderer change keeps only the params the new
// renderer accepts, an invalid Depth track index is not replaced silently,
// Hide GC Skew reaches only the GC rows, Reset follows the loaded Depth
// series, and Duplicate is unavailable while a Session operation runs.
import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-track-slot-row-edits-'));
for (const directory of ['app', 'utils', 'contracts']) {
  await cp(join(repoRoot, 'gbdraw', 'web', 'js', directory), join(tempRoot, directory), { recursive: true });
}
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// track-slot-display.js reads the generated size-class boundaries.
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.generated.js'), join(tempRoot, 'mode-profiles.generated.js'));
// The services the slot editors import.
for (const leaf of ['conservation-series.js', 'current-option-values.js', 'depth-track-state.js', 'linear-label-visibility.js', 'losat-normalization.js',
  'track-slot-display.js', 'track-slot-validation.js', 'circular-track-slot-model.js', 'linear-track-slot-model.js']) {
  await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services', leaf), join(tempRoot, 'services', leaf));
}

const { createCircularTrackSlotEditor } = await import(pathToFileURL(join(tempRoot, 'app', 'circular-track-slots.js')));
const { createLinearTrackSlotEditor } = await import(pathToFileURL(join(tempRoot, 'app', 'linear-track-slots.js')));

const depthFile = { name: 'a.depth.tsv', size: 1, lastModified: 1 };

const circularState = (slots, { showDepth = false, depthFiles = [] } = {}) => withDrawings({
  mode: { value: 'circular' },
  form: { track_type: 'tuckin', show_depth: showDepth, suppress_gc: false, suppress_skew: false, show_scale: true, separate_strands: true },
  adv: {
    circular_track_slots_enabled: true,
    circular_track_slots_axis_index: 1,
    circular_track_slots: structuredClone(slots),
    nt: 'GC',
    features: ['CDS'],
    feature_shapes: { CDS: 'arrow' },
    depth_tracks: [],
    feature_width_circular: null,
    depth_width_circular: null,
    gc_content_width_circular: null,
    gc_content_radius_circular: null,
    gc_skew_width_circular: null,
    gc_skew_radius_circular: null
  },
  files: { c_depth: depthFiles, c_conservation_blasts: [], c_conservation_fastas: [] },
  circularConservation: { enabled: false, source: 'upload', labels: '', series: [] },
  annotationSets: [],
  circularRecordList: { value: [] }
});

const linearState = (slots, { showDepth = false, depthFiles = [] } = {}) => withDrawings({
  form: { linear_track_layout: 'middle', show_depth: showDepth, show_gc: true, show_skew: true },
  adv: {
    nt: 'GC',
    linear_track_slots_enabled: true,
    linear_track_slots_axis_index: 0,
    linear_track_slots: structuredClone(slots),
    depth_tracks: [],
    depth_height: null,
    gc_height: null
  },
  linearSeqs: [{ depth: depthFiles }],
  annotationSets: []
});

const circularRow = (id, renderer, params = {}, side = 'inside') => ({ id, renderer, enabled: true, side, z: 0, params });
const linearRow = (id, renderer, params = {}, side = 'below') => ({ id, renderer, enabled: true, side, height: '', spacing: '', z: 0, params });

const MODES = {
  circular: (slots, options, askChoice = async () => 'hide') => {
    const state = circularState(slots, options);
    const editor = createCircularTrackSlotEditor({ state, askChoice });
    return {
      state,
      editor,
      slots: () => state.adv.circular_track_slots,
      setRenderer: (slot, renderer) => editor.updateCircularTrackSlotRenderer(slot, renderer),
      issue: (index) => editor.circularTrackSlotIssue(state.adv.circular_track_slots[index], index),
      reset: () => editor.resetCircularTrackSlotsFromSimpleControls(),
      canDuplicate: (slot) => editor.canDuplicateCircularTrackSlot(slot),
      duplicate: (index) => editor.duplicateCircularTrackSlot(index)
    };
  },
  linear: (slots, options) => {
    const state = linearState(slots, options);
    const editor = createLinearTrackSlotEditor({ state });
    return {
      state,
      editor,
      slots: () => state.adv.linear_track_slots,
      setRenderer: (slot, renderer) => editor.updateLinearTrackSlotRenderer(slot, renderer),
      issue: (index) => editor.linearTrackSlotIssue(state.adv.linear_track_slots[index], index),
      reset: () => editor.resetLinearTrackSlotsFromSimpleControls(),
      canDuplicate: (slot) => editor.canDuplicateLinearTrackSlot(slot),
      duplicate: (index) => editor.duplicateLinearTrackSlot(index)
    };
  }
};

const ROWS = {
  circular: () => [
    circularRow('features', 'features', { lane_direction: 'outside' }, 'outside'),
    circularRow('ticks', 'ticks', { tick_label_layout: 'label_out_tick_in' }),
    circularRow('gc_skew', 'dinucleotide_skew', { nt: 'GC', positive_color: '#ff0000' })
  ],
  linear: () => [
    linearRow('features', 'features', {}, 'overlay'),
    linearRow('gc_skew', 'dinucleotide_skew', { nt: 'GC', positive_color: '#ff0000' })
  ]
};

for (const [mode, create] of Object.entries(MODES)) {
  // TK-04: the params of the previous renderer made the row invalid, also
  // after switching back.
  test(`${mode}: a renderer change keeps only the params the new renderer accepts`, () => {
    const editor = create(ROWS[mode](), { depthFiles: [depthFile] });
    const index = editor.slots().findIndex((slot) => slot.id === 'gc_skew');
    const slot = editor.slots()[index];

    editor.setRenderer(slot, 'depth');
    // Circular leaves an unset index (series 0); Linear writes 0.
    assert.equal(editor.slots()[index].params.track_index ?? 0, 0);
    assert.deepEqual(Object.keys(editor.slots()[index].params).filter((key) => key !== 'track_index'), []);
    assert.equal(editor.issue(index), '');

    // A Depth row's legend label is its series title (syncDepthSlotLabels).
    editor.slots()[index].params.legend_label = 'hmmt_depth_a';
    editor.setRenderer(editor.slots()[index], 'dinucleotide_content');
    assert.deepEqual({ ...editor.slots()[index].params }, { nt: 'GC' });
    assert.equal(editor.issue(index), '');

    editor.slots()[index].params.nt = 'AT';
    editor.setRenderer(editor.slots()[index], 'dinucleotide_skew');
    assert.deepEqual({ ...editor.slots()[index].params }, { nt: 'AT' }, 'nt belongs to both dinucleotide renderers');
    assert.equal(editor.issue(index), '');
  });

  // GX-01: Duplicate uses the availability term of the other row controls.
  test(`${mode}: Duplicate is unavailable while a Session operation runs`, () => {
    const editor = create(ROWS[mode]());
    const slot = editor.slots().find((candidate) => candidate.id === 'gc_skew');
    assert.equal(editor.canDuplicate(slot), true);
    editor.state.sessionOperationAvailability = () => ({ code: 'SESSION_BUSY' });
    assert.equal(editor.canDuplicate(slot), false);
  });

  // OV-1099: both modes name a copy by its renderer's default ID and the first
  // unused number, so a copy of a copy is numbered on, not suffixed again.
  test(`${mode}: Duplicate names the copy with the next unused number`, () => {
    const editor = create(ROWS[mode]());
    const ids = () => editor.slots().map((slot) => slot.id);
    editor.duplicate(ids().indexOf('gc_skew'));
    assert.deepEqual(ids().slice(-2), ['gc_skew', 'gc_skew_2']);
    editor.duplicate(ids().indexOf('gc_skew_2'));
    assert.deepEqual(ids().slice(-3), ['gc_skew', 'gc_skew_2', 'gc_skew_3']);
    const copy = editor.slots().at(-1);
    assert.equal(copy.params.positive_color, '#ff0000', 'the copy keeps the source params');
  });

  // TK-10: Reset rebuilds from the loaded Depth series; Show Depth is not an
  // input (docs/REFERENCE/web-app.md, Custom Track Slots).
  test(`${mode}: Reset keeps one Depth row per loaded series while Show Depth is off`, () => {
    for (const showDepth of [false, true]) {
      const editor = create(ROWS[mode](), { showDepth, depthFiles: [depthFile] });
      editor.reset();
      const depthRows = editor.slots().filter((slot) => slot.renderer === 'depth');
      assert.deepEqual(depthRows.map((slot) => slot.params.track_index ?? 0), [0], `Show Depth ${showDepth}`);
    }
    const withoutDepth = create(ROWS[mode](), { showDepth: false, depthFiles: [] });
    withoutDepth.reset();
    assert.equal(withoutDepth.slots().some((slot) => slot.renderer === 'depth'), false);
  });
}

test('circular: Reset to preset keeps the Depth row while Show Depth is off', () => {
  const editor = MODES.circular(ROWS.circular(), { showDepth: false, depthFiles: [depthFile] });
  editor.editor.resetCircularTrackSlotsToPreset('middle');
  assert.equal(editor.slots().filter((slot) => slot.renderer === 'depth').length, 1);
});

test('circular: a renderer change back to Ticks leaves a valid Ticks row', () => {
  const editor = MODES.circular(ROWS.circular());
  const index = editor.slots().findIndex((slot) => slot.id === 'ticks');
  editor.setRenderer(editor.slots()[index], 'dinucleotide_content');
  assert.deepEqual({ ...editor.slots()[index].params }, { nt: 'GC' });
  assert.equal(editor.issue(index), '');
  editor.setRenderer(editor.slots()[index], 'ticks');
  assert.equal(editor.slots()[index].params.nt, undefined);
  assert.equal(editor.issue(index), '');
});

// TK-09: Hide GC Skew draws no GC skew in the simple path; in the custom stack
// it disables only the rows of the drawing's dinucleotide.
test('circular: Hide GC Skew disables only the GC skew rows and names them in the right number', async () => {
  // The app's choice dialog asks; the stub records the question and hides.
  const questions = [];
  const editor = MODES.circular([
    ...ROWS.circular(),
    circularRow('gc_skew_2', 'dinucleotide_skew', { nt: 'AT' })
  ], undefined, async (request) => { questions.push(request); return 'hide'; });
  await editor.editor.setCircularSkewSuppressed(true);
  const rows = Object.fromEntries(editor.slots().map((slot) => [slot.id, slot.enabled]));
  assert.equal(rows.gc_skew, false);
  assert.equal(rows.gc_skew_2, true, 'the AT skew row stays enabled');
  assert.equal(editor.editor.circularTrackSlotSuppressMessage(editor.slots().find((slot) => slot.id === 'gc_skew_2')), '');
  assert.equal(questions.length, 1);
  assert.equal(questions[0].title, 'Hide GC skew?');
  assert.match(questions[0].description, /has 1 GC skew track\. Hiding GC skew turns that track off\./);
  assert.deepEqual(questions[0].actions.map((action) => action.label), ['Hide GC skew']);

  await editor.editor.setCircularSkewSuppressed(false);
  assert.equal(editor.slots().every((slot) => slot.enabled), true);

  editor.editor.duplicateCircularTrackSlot(editor.slots().findIndex((slot) => slot.id === 'gc_skew'));
  await editor.editor.setCircularSkewSuppressed(true);
  assert.match(questions[1].description, /has 2 GC skew tracks\. Hiding GC skew turns those tracks off\./);

  // Keep leaves the setting and the rows unchanged.
  const kept = MODES.circular(ROWS.circular(), undefined, async () => 'cancel');
  assert.equal(await kept.editor.setCircularSkewSuppressed(true), false);
  assert.equal(kept.state.form.suppress_skew, false);
  assert.equal(kept.slots().every((slot) => slot.enabled), true);
});

// The question is asked before any History step (app-setup.js gives the editor
// History's `runUndoable` and the app's one choice prompt): the chosen hide is
// one step recorded at the answer, and Keep records none. A hide chosen while
// another action is still running applies at once as its own step, and what
// that action writes afterwards is a step under its name (OV-460).
test('circular: a hide that asks is one History step recorded at the answer, also while another action runs', async () => {
  const { createChoicePrompt } = await import(pathToFileURL(join(tempRoot, 'app', 'app-messages.js')));
  const { createHistoryManager } = await import('../../gbdraw/web/js/services/history.js');
  const state = Object.assign(circularState(ROWS.circular()), { loaded: 0 });
  const { ask, resolve } = createChoicePrompt({ reactive: (value) => value });
  const history = createHistoryManager({
    buildIntent: () => structuredClone({ hide: state.form.suppress_skew, slots: state.adv.circular_track_slots, loaded: state.loaded }),
    applyIntent: (intent) => {
      state.form.suppress_skew = intent.hide;
      state.adv.circular_track_slots.splice(0, Infinity, ...structuredClone(intent.slots));
      state.loaded = intent.loaded;
    },
    buildCheckpoint: () => ({}), applyCheckpoint: () => {}
  });
  const editor = createCircularTrackSlotEditor({ state, askChoice: ask, runUndoable: history.runUndoable });
  const read = () => [state.form.suppress_skew, state.adv.circular_track_slots.find((slot) => slot.id === 'gc_skew').enabled, state.loaded];
  const steps = () => [history.getUndoCount(), history.undoLabel()];
  const tick = () => new Promise((settle) => setImmediate(settle));
  await history.initializeIntentBaseline();
  /** @param {string} choice */
  const hide = async (choice) => {
    const asked = editor.setCircularSkewSuppressed(true);
    await tick();
    assert.deepEqual(read().slice(0, 2), [false, true], 'nothing is written before the answer');
    resolve(choice);
    return asked;
  };

  const asking = editor.setCircularSkewSuppressed(true);
  await tick();
  assert.equal(history.mutationPending(), false, 'no step is open while the question is asked');
  resolve('cancel');
  assert.equal(await asking, false);
  assert.deepEqual([read(), history.getUndoCount()], [[false, true, 0], 0]);

  await hide('hide');
  assert.deepEqual([read(), steps()], [[true, false, 0], [1, 'Hide GC Skew']]);
  await history.undo();
  assert.deepEqual(read(), [false, true, 0]);
  await history.redo();
  assert.deepEqual(read(), [true, false, 0]);
  await history.undo();

  let release;
  const held = new Promise((settle) => { release = settle; });
  const load = history.runUndoable('Load label edits', async () => { await held; state.loaded = 1; });
  await tick();
  await hide('hide');
  assert.deepEqual([read(), steps(), history.mutationPending()], [[true, false, 0], [1, 'Hide GC Skew'], false]);
  release();
  await load;
  assert.deepEqual([read(), steps()], [[true, false, 1], [2, 'Load label edits']]);
  await history.undo();
  await history.undo();
  assert.deepEqual(read(), [false, true, 0]);
  await history.redo();
  await history.redo();
  assert.deepEqual(read(), [true, false, 1]);
});

// TK-07: `1.5` or `-1` became track index 1 without a message.
test('circular: an invalid Depth track index stays visible with a field error and the row keeps its index', async () => {
  const Vue = loadPageVue();
  globalThis.window = { Vue };
  const { reactive, nextTick } = Vue;
  const { DepthTrackIndexInput } = await import(pathToFileURL(join(tempRoot, 'app', 'circular-track-slots', 'track-index-input.js')));
  const props = reactive({ modelValue: 1, slotId: 'depth_2', controlId: 'slot-track-index', helpId: 'help-index', disabled: false });
  const emissions = [];
  const field = DepthTrackIndexInput.setup(props, { emit: (event, value) => {
    emissions.push(value);
    props.modelValue = value;
  } });

  for (const entry of ['1.5', '-1', '']) {
    field.commit(entry);
    assert.equal(field.text.value, entry);
    assert.equal(field.error.value, 'Enter a whole number of 0 or more. The row keeps track index 1.');
    assert.equal(field.describedBy.value, 'help-index slot-track-index-error');
  }
  assert.deepEqual(emissions, [], 'the row keeps its last valid index');

  field.commit('2');
  assert.deepEqual(emissions, [2]);
  assert.equal(field.text.value, '2');
  assert.equal(field.error.value, '');
  assert.equal(field.describedBy.value, 'help-index');

  field.commit('-1');
  props.modelValue = 0; // Undo restores another index and drops the rejected entry.
  await nextTick();
  assert.equal(field.text.value, '0');
  assert.equal(field.error.value, '');
});

// TK-09 review: the Dinucleotide setting and a row's dinucleotide are inputs of
// the Hide GC projection, so their edits run it; the editor state equals what
// Generate projects (run-analysis applies the same projection to the request).
test('circular: Dinucleotide and row dinucleotide edits re-run the Hide GC Skew projection', async () => {
  const { applyCircularSuppressControlsToSlots } = await import(pathToFileURL(join(tempRoot, 'app', 'circular-track-slots.js')));
  const { applyCircularTrackOrderPlacements } = await import(pathToFileURL(join(tempRoot, 'services', 'circular-track-slot-model.js')));
  {
    const editor = MODES.circular(ROWS.circular());
    const { state } = editor;
    const requestSlots = () => applyCircularSuppressControlsToSlots(
      applyCircularTrackOrderPlacements(state.adv.circular_track_slots, state.adv.nt, state.form.track_type, state.adv.circular_track_slots_axis_index),
      state.form,
      state.adv.nt
    );
    const enabled = () => editor.slots().find((slot) => slot.id === 'gc_skew').enabled;
    await editor.editor.setCircularSkewSuppressed(true);
    assert.equal(enabled(), false);

    editor.editor.setCircularDinucleotide('AT');
    assert.equal(state.adv.nt, 'AT');
    assert.equal(enabled(), true, 'the GC row no longer uses the Dinucleotide setting');
    assert.deepEqual(JSON.parse(JSON.stringify(editor.slots())), JSON.parse(JSON.stringify(requestSlots())));

    const row = editor.slots().find((slot) => slot.id === 'gc_skew');
    editor.editor.setCircularTrackSlotNt(row, ' AT ');
    assert.equal(row.params.nt, 'AT');
    assert.equal(enabled(), false, 'the row now uses the Dinucleotide setting');
    assert.deepEqual(JSON.parse(JSON.stringify(editor.slots())), JSON.parse(JSON.stringify(requestSlots())));
  }
});
