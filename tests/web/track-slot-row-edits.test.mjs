// Custom stack row edits through the Circular and Linear editors (TK-04,
// TK-07, TK-09, TK-10, GX-01): a renderer change keeps only the params the new
// renderer accepts, an invalid Depth track index is not replaced silently,
// Hide GC Skew reaches only the GC rows, Reset follows the loaded Depth
// series, and Duplicate is unavailable while a Session operation runs.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import vm from 'node:vm';
import { withDrawings } from './helpers/drawing-state.mjs';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-track-slot-row-edits-'));
for (const directory of ['app', 'utils']) {
  await cp(join(repoRoot, 'gbdraw', 'web', 'js', directory), join(tempRoot, directory), { recursive: true });
}
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// The services the slot editors import.
for (const leaf of ['conservation-series.js', 'current-option-values.js', 'depth-track-state.js', 'track-slot-display.js',
  'track-slot-validation.js', 'circular-track-slot-model.js', 'linear-track-slot-model.js']) {
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
  circular: (slots, options) => {
    const state = circularState(slots, options);
    const editor = createCircularTrackSlotEditor({ state });
    return {
      state,
      editor,
      slots: () => state.adv.circular_track_slots,
      setRenderer: (slot, renderer) => editor.updateCircularTrackSlotRenderer(slot, renderer),
      issue: (index) => editor.circularTrackSlotIssue(state.adv.circular_track_slots[index], index),
      reset: () => editor.resetCircularTrackSlotsFromSimpleControls(),
      canDuplicate: (slot) => editor.canDuplicateCircularTrackSlot(slot)
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
      canDuplicate: (slot) => editor.canDuplicateLinearTrackSlot(slot)
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
test('circular: Hide GC Skew disables only the GC skew rows and names them in the right number', () => {
  const messages = [];
  const previousConfirm = globalThis.confirm;
  globalThis.confirm = (message) => { messages.push(message); return true; };
  try {
    const editor = MODES.circular([
      ...ROWS.circular(),
      circularRow('gc_skew_2', 'dinucleotide_skew', { nt: 'AT' })
    ]);
    editor.editor.setCircularSkewSuppressed(true);
    const rows = Object.fromEntries(editor.slots().map((slot) => [slot.id, slot.enabled]));
    assert.equal(rows.gc_skew, false);
    assert.equal(rows.gc_skew_2, true, 'the AT skew row stays enabled');
    assert.equal(editor.editor.circularTrackSlotSuppressMessage(editor.slots().find((slot) => slot.id === 'gc_skew_2')), '');
    assert.equal(messages.length, 1);
    assert.match(messages[0], /include 1 enabled GC skew track\./);
    assert.match(messages[0], /will disable that custom track slot and/);

    editor.editor.setCircularSkewSuppressed(false);
    assert.equal(editor.slots().every((slot) => slot.enabled), true);

    editor.editor.duplicateCircularTrackSlot(editor.slots().findIndex((slot) => slot.id === 'gc_skew'));
    editor.editor.setCircularSkewSuppressed(true);
    assert.match(messages[1], /include 2 enabled GC skew tracks\./);
    assert.match(messages[1], /will disable those custom track slots and/);
  } finally {
    globalThis.confirm = previousConfirm;
  }
});

// TK-07: `1.5` or `-1` became track index 1 without a message.
test('circular: an invalid Depth track index stays visible with a field error and the row keeps its index', async () => {
  const context = vm.createContext({ console });
  vm.runInContext(readFileSync('gbdraw/web/vendor/vue/vue.global.js', 'utf8'), context);
  globalThis.window = { Vue: context.Vue };
  const { reactive, nextTick } = context.Vue;
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
