import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-circular-track-slots-'));
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app'),
  join(tempRoot, 'app'),
  { recursive: true }
);
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'utils'),
  join(tempRoot, 'utils'),
  { recursive: true }
);
// circular-track-slots.js reads the conservation-series row helpers.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'conservation-series.js'),
  join(tempRoot, 'services', 'conservation-series.js')
);
// The slice editors and svg-serialization resolve a drawing's ID by its mode.
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'drawing-view.js'), join(tempRoot, 'services', 'drawing-view.js'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// track-slot-display.js reads the generated size-class boundaries.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'mode-profiles.generated.js'),
  join(tempRoot, 'mode-profiles.generated.js')
);
// color-utils.js reads the generated palette name aliases.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'contracts', 'schema.generated.js'),
  join(tempRoot, 'contracts', 'schema.generated.js')
);
// The track-slot leaves the slot editors import.
for (const leaf of ['depth-track-state.js', 'track-slot-display.js', 'track-slot-validation.js', 'circular-track-slot-model.js']) {
  await cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'services', leaf),
    join(tempRoot, 'services', leaf)
  );
}

const { createCircularTrackSlotEditor } = await import(
  pathToFileURL(join(tempRoot, 'app', 'circular-track-slots.js'))
);

const featureSlot = () => ({
  id: 'features',
  renderer: 'features',
  enabled: true,
  side: 'outside',
  z: 0,
  params: { lane_direction: 'outside' }
});

const spacerSlot = () => ({
  id: 'outer_gap',
  renderer: 'spacer',
  enabled: true,
  side: 'outside',
  z: 0,
  params: {}
});

const gcSlot = () => ({
  id: 'gc_content',
  renderer: 'dinucleotide_content',
  enabled: true,
  side: 'inside',
  z: 0,
  params: { nt: 'GC' }
});

const tickSlot = (side = 'overlay') => ({
  id: 'ticks',
  renderer: 'ticks',
  enabled: true,
  side,
  z: 0,
  params: { tick_label_layout: 'label_out_tick_in' }
});

const depthSlot = (side = 'inside') => ({
  id: 'depth',
  renderer: 'depth',
  enabled: true,
  side,
  z: 0,
  params: { track_index: 0 }
});

const conservationSlot = (side = 'inside') => ({
  id: 'conservation_pair',
  renderer: 'sequence_conservation',
  enabled: true,
  side,
  z: 0,
  params: {
    managed: 'circular_conservation',
    series_key: 'pair.tsv|1|0|0',
    source_index: 0,
    track_index: 1,
    label: 'Pair',
    color: '#4e79a7',
    fileName: 'pair.tsv'
  }
});

const createState = ({
  slots = [featureSlot(), spacerSlot(), gcSlot()],
  axisIndex = 2,
  showDepth = false,
  depthFiles = [],
  conservationEnabled = false,
  conservationFiles = []
} = {}) => ({
  mode: { value: 'circular' },
  form: {
    track_type: 'tuckin',
    show_depth: showDepth,
    suppress_gc: false,
    suppress_skew: true,
    show_scale: false,
    separate_strands: true
  },
  adv: {
    circular_track_slots_enabled: true,
    circular_track_slots_axis_index: axisIndex,
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
  files: {
    c_depth: depthFiles,
    c_conservation_blasts: conservationFiles,
    c_conservation_fastas: []
  },
  circularConservation: {
    enabled: conservationEnabled,
    source: 'upload',
    labels: '',
    series: []
  },
  annotationSets: [],
  circularRecordList: { value: [] }
});

const slotSides = (state) => Object.fromEntries(
  state.adv.circular_track_slots.map((slot) => [slot.id, slot.side])
);

// The preview estimates a row with the size class Generate will use: 60,000 bp
// is short by the packaged Circular boundary and long by a boundary of 50 kb
// that the drawing keeps from its Session.
test('the track preview uses the short/long boundary the Circular drawing keeps', () => {
  const featureWidthNote = (unmanagedConfigOverrides) => {
    const state = { ...createState(), unmanagedConfigOverrides };
    const editor = createCircularTrackSlotEditor({ state: withDrawings(state), estimateRecordLengths: () => [60_000] });
    return editor.circularTrackSlotGeometryAutoText(state.adv.circular_track_slots[0], 0, 'width');
  };
  assert.equal(featureWidthNote({}), '≈ 37 px (estimate)');
  assert.equal(featureWidthNote({ 'labels.length_threshold.circular': 50_000 }), '≈ 18.5 px (estimate)');
});

// OV-494: the estimates take the size class of the records they describe
// (`circularEstimateLengths`), and the long class only while no length is known.
test('the track preview estimates the records it is given', () => {
  const featureWidthNote = (/** @type {number[]} */ recordLengths) => {
    const state = { ...createState(), unmanagedConfigOverrides: {} };
    const editor = createCircularTrackSlotEditor({ state: withDrawings(state), estimateRecordLengths: () => recordLengths });
    return editor.circularTrackSlotGeometryAutoText(state.adv.circular_track_slots[0], 0, 'width');
  };
  assert.equal(featureWidthNote([30_000]), '≈ 37 px (estimate)');
  assert.equal(featureWidthNote([30_000, 300_000]), '≈ 18.5 px (estimate)');
  assert.equal(featureWidthNote([]), '≈ 18.5 px (estimate)', 'long while no length is known');
});

// OV-1409: the global issue list keys each issue by its code, not its shown text.
test('the global issues keep their codes', () => {
  const editor = createCircularTrackSlotEditor({ state: withDrawings(createState({ axisIndex: 9 })) });
  assert.deepEqual(editor.circularTrackGlobalIssues().map((issue) => [issue.code, issue.parts]), [
    ['axis_out_of_range', [{ text: 'The Axis boundary must be an integer from 0 to {max}.', params: { max: 3 } }]]
  ]);
  assert.deepEqual(createCircularTrackSlotEditor({ state: withDrawings(createState()) }).circularTrackGlobalIssues(), []);
});

// OV-487: Python draws a content or skew slot without a dinucleotide with the
// drawing's, so the closed disclosure names that one, not GC.
test('the slot options summary names the dinucleotide the slot is drawn with', () => {
  const state = createState();
  state.adv.nt = 'AT';
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  const summary = (renderer, params) => editor.circularTrackSlotOptionsSummary({ id: 'slot', renderer, params });
  assert.equal(summary('dinucleotide_skew', {}), 'AT · legend label, gaps, colors');
  assert.equal(summary('dinucleotide_content', { nt: ' ' }), 'AT · legend label, gaps');
  assert.equal(summary('dinucleotide_content', { nt: 'gc' }), 'GC · legend label, gaps');
  assert.equal(summary('features', {}), 'Gaps');
});

// The default tracks draw the drawing's dinucleotide, so the preset summary
// names it too. The tracks are a list (`listText`) of translated names.
test('the preset summary names the default tracks by the dinucleotide of the drawing', () => {
  const state = createState();
  state.form.suppress_skew = false;
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  assert.match(editor.circularTrackPresetSummary().detail, /GC content and GC skew inside$/);
  state.adv.nt = ' at ';
  assert.match(editor.circularTrackPresetSummary().detail, /AT content and AT skew inside$/);
  state.form.show_depth = true;
  assert.match(editor.circularTrackPresetSummary().detail, /Depth, AT content, and AT skew inside$/);
});

test('Circular placement actions preserve the Axis boundary and radial sides', () => {
  const state = createState({
    slots: [
      { ...gcSlot(), side: 'outside' },
      { ...featureSlot(), side: 'inside', params: { lane_direction: 'inside' } },
      tickSlot('inside'),
      {
        id: 'gc_skew',
        renderer: 'dinucleotide_skew',
        enabled: true,
        side: 'inside',
        z: 0,
        params: { nt: 'GC' }
      }
    ],
    axisIndex: 1
  });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });

  assert.equal(editor.canMoveCircularTrackSlot(1, -1), false);
  editor.moveCircularTrackSlot(1, 0);
  assert.deepEqual(
    state.adv.circular_track_slots.map((slot) => slot.id),
    ['gc_content', 'features', 'ticks', 'gc_skew']
  );

  editor.moveCircularTrackSlotOutside(1);
  assert.deepEqual(slotSides(state), {
    gc_content: 'outside',
    features: 'outside',
    ticks: 'inside',
    gc_skew: 'inside'
  });
  assert.equal(state.adv.circular_track_slots[1].params.lane_direction, 'outside');
});

const depthFile = (name) => ({ name, size: 1, lastModified: 0 });

test('managed Depth insertion preserves the existing Axis sides', () => {
  const state = createState({ showDepth: true });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });

  editor.changeCircularDepthSources(() => {
    state.files.c_depth = [depthFile('depth.tsv')];
  });

  assert.deepEqual(
    state.adv.circular_track_slots.map((slot) => slot.id),
    ['features', 'outer_gap', 'depth', 'gc_content']
  );
  assert.equal(state.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(state), {
    features: 'outside',
    outer_gap: 'outside',
    depth: 'inside',
    gc_content: 'inside'
  });
});

test('managed Depth insertion stays after an on-Axis Tick', () => {
  const existingDepth = {
    ...depthSlot('outside'),
    id: 'depth_1'
  };
  const state = createState({
    slots: [featureSlot(), existingDepth, tickSlot(), gcSlot()],
    axisIndex: 2,
    showDepth: true,
    depthFiles: [depthFile('depth-a.tsv')]
  });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });

  editor.changeCircularDepthSources(() => {
    state.files.c_depth = [depthFile('depth-a.tsv'), depthFile('depth-b.tsv')];
  });

  assert.deepEqual(
    state.adv.circular_track_slots.map((slot) => slot.id),
    ['features', 'depth_1', 'ticks', 'depth_2', 'gc_content']
  );
  assert.equal(state.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(state), {
    features: 'outside',
    depth_1: 'outside',
    ticks: 'overlay',
    depth_2: 'inside',
    gc_content: 'inside'
  });
  assert.equal(editor.circularTrackSlotIssue(state.adv.circular_track_slots[3], 3), '');
});

test('managed Depth removal rebases the Axis only when the row is before it', () => {
  const removeSource = (state) => createCircularTrackSlotEditor({ state: withDrawings(state) }).changeCircularDepthSources(() => {
    state.files.c_depth = [];
  });
  const beforeState = createState({
    slots: [depthSlot('outside'), featureSlot(), spacerSlot(), gcSlot()],
    axisIndex: 3,
    depthFiles: [depthFile('depth.tsv')]
  });
  removeSource(beforeState);
  assert.equal(beforeState.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(beforeState), {
    features: 'outside',
    outer_gap: 'outside',
    gc_content: 'inside'
  });

  const afterState = createState({
    slots: [featureSlot(), spacerSlot(), gcSlot(), depthSlot()],
    axisIndex: 2,
    depthFiles: [depthFile('depth.tsv')]
  });
  removeSource(afterState);
  assert.equal(afterState.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(afterState), {
    features: 'outside',
    outer_gap: 'outside',
    gc_content: 'inside'
  });
});

test('managed Conservation insertion preserves the existing Axis sides', () => {
  const state = createState({
    conservationEnabled: true,
    conservationFiles: [{ name: 'pair.tsv', size: 1, lastModified: 0 }]
  });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });

  editor.syncCircularConservationSlots();

  const managed = state.adv.circular_track_slots.find(
    (slot) => slot.renderer === 'sequence_conservation'
  );
  assert.ok(managed);
  assert.deepEqual(
    state.adv.circular_track_slots.map((slot) => slot.id),
    ['features', 'outer_gap', managed.id, 'gc_content']
  );
  assert.equal(state.adv.circular_track_slots_axis_index, 2);
  assert.equal(managed.side, 'inside');
  assert.equal(slotSides(state).features, 'outside');
  assert.equal(slotSides(state).outer_gap, 'outside');
  assert.equal(slotSides(state).gc_content, 'inside');
});

test('managed Conservation insertion stays after an on-Axis Tick', () => {
  const state = createState({
    slots: [featureSlot(), tickSlot(), gcSlot()],
    axisIndex: 1,
    conservationEnabled: true,
    conservationFiles: [{ name: 'pair.tsv', size: 1, lastModified: 0 }]
  });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });

  editor.syncCircularConservationSlots();

  const managed = state.adv.circular_track_slots.find(
    (slot) => slot.renderer === 'sequence_conservation'
  );
  assert.ok(managed);
  assert.deepEqual(
    state.adv.circular_track_slots.map((slot) => slot.id),
    ['features', 'ticks', managed.id, 'gc_content']
  );
  assert.equal(state.adv.circular_track_slots_axis_index, 1);
  assert.deepEqual(slotSides(state), {
    features: 'outside',
    ticks: 'overlay',
    [managed.id]: 'inside',
    gc_content: 'inside'
  });
  assert.equal(editor.circularTrackSlotIssue(managed, 2), '');
});

test('managed Conservation removal rebases the Axis only when the row is before it', () => {
  const beforeState = createState({
    slots: [conservationSlot('outside'), featureSlot(), spacerSlot(), gcSlot()],
    axisIndex: 3
  });
  createCircularTrackSlotEditor({ state: withDrawings(beforeState) }).syncCircularConservationSlots();
  assert.equal(beforeState.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(beforeState), {
    features: 'outside',
    outer_gap: 'outside',
    gc_content: 'inside'
  });

  const afterState = createState({
    slots: [featureSlot(), spacerSlot(), gcSlot(), conservationSlot()],
    axisIndex: 2
  });
  createCircularTrackSlotEditor({ state: withDrawings(afterState) }).syncCircularConservationSlots();
  assert.equal(afterState.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(afterState), {
    features: 'outside',
    outer_gap: 'outside',
    gc_content: 'inside'
  });
});

test('Circular measure action updates only the owned scalar and preserves slot identity', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  const slot = state.adv.circular_track_slots[0];
  const key = editor.circularTrackSlotEditorKey(slot);
  const scalar = { value: '1e', unit: 'px' };
  const radius = slot.radius;
  editor.updateCircularTrackSlotMeasure(slot, 'width', scalar);
  assert.strictEqual(slot.width, scalar);
  assert.strictEqual(slot.radius, radius);
  assert.equal(editor.circularTrackSlotEditorKey(slot), key);
  editor.updateCircularTrackSlotMeasure(slot, 'width', scalar);
  assert.strictEqual(slot.width, scalar);
  editor.updateCircularTrackSlotMeasure(slot, 'renderer', scalar);
  assert.equal(slot.renderer, 'features');
  const detached = { width: 'saved' };
  editor.updateCircularTrackSlotMeasure(detached, 'width', null);
  assert.equal(detached.width, 'saved');
});

test('Circular measure action refuses nonfinite numeric leaves before changing a slot', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  const slot = state.adv.circular_track_slots[0];
  for (const field of ['width', 'radius']) {
    const before = slot[field];
    for (const scalar of [NaN, Infinity, -Infinity,
      { value: NaN, unit: 'px' }, { value: Infinity, unit: 'factor' },
      { value: -Infinity, unit: 'px' }]) {
      assert.throws(() => editor.updateCircularTrackSlotMeasure(slot, field, scalar), /must be a number greater than 0, in px or ×R/);
      assert.strictEqual(slot[field], before);
    }
  }
  for (const scalar of [0.5, { value: 20, unit: 'px' }, null,
    { value: '1e-3', unit: 'factor' }, { value: 'Infinity', unit: 'px' },
    { value: '1e', unit: 'px' }, { value: '20em', unit: 'px' }]) {
    editor.updateCircularTrackSlotMeasure(slot, 'width', scalar);
    assert.strictEqual(slot.width, scalar);
    editor.updateCircularTrackSlotMeasure(slot, 'width', scalar);
    assert.strictEqual(slot.width, scalar);
  }
});

const gcRows = (state) => state.adv.circular_track_slots
  .filter((slot) => slot.renderer === 'dinucleotide_content')
  .map((slot) => ({ id: slot.id, enabled: slot.enabled, marked: '_suppressed_by_global' in slot.params }));

// A2 W4: turning a non-empty custom stack on adds the managed comparison ring
// rows in that call (no watcher follows the toggle).
test('turning the custom stack on adds the managed comparison ring rows', () => {
  const state = createState({
    conservationEnabled: true,
    conservationFiles: [{ name: 'pair.tsv', size: 1, lastModified: 0 }]
  });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  editor.setCircularTrackSlotsEnabled(false);
  assert.equal(state.adv.circular_track_slots.some((slot) => slot.renderer === 'sequence_conservation'), false);

  editor.setCircularTrackSlotsEnabled(true);

  assert.deepEqual(state.adv.circular_track_slots.filter((slot) => slot.renderer === 'sequence_conservation')
    .map((slot) => slot.params.fileName), ['pair.tsv']);
});

// TR-06: Hide GC Content owns the visibility of its rows whether or not the
// custom stack is in use.
test('un-hiding GC while the custom stack is off restores the row the hide disabled', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });

  editor.setCircularGcSuppressed(true);
  editor.setCircularTrackSlotsEnabled(false);
  editor.setCircularGcSuppressed(false);
  editor.setCircularTrackSlotsEnabled(true);

  assert.deepEqual(gcRows(state), [{ id: 'gc_content', enabled: true, marked: false }]);
});

test('a GC row the user disabled stays disabled across Hide and un-hide', async () => {
  const state = createState({
    slots: [featureSlot(), spacerSlot(), { ...gcSlot(), enabled: false }, { ...gcSlot(), id: 'gc_content_2' }]
  });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state), askChoice: async () => 'hide' });

  await editor.setCircularGcSuppressed(true);
  assert.deepEqual(gcRows(state), [
    { id: 'gc_content', enabled: false, marked: false },
    { id: 'gc_content_2', enabled: false, marked: true }
  ]);
  editor.duplicateCircularTrackSlot(2);
  editor.setCircularGcSuppressed(false);

  assert.deepEqual(gcRows(state), [
    { id: 'gc_content', enabled: false, marked: false },
    { id: 'gc_content_3', enabled: false, marked: false },
    { id: 'gc_content_2', enabled: true, marked: false }
  ]);
});

// TR-09: every Reset honors Show Coordinate Scale; the first Use custom stack
// still uses the saved stack (D-34).
test('Reset and the preset resets omit Ticks while Show Coordinate Scale is off', () => {
  for (const preset of ['tuckin', 'middle', 'spreadout']) {
    const state = createState();
    createCircularTrackSlotEditor({ state: withDrawings(state) }).resetCircularTrackSlotsToPreset(preset);
    assert.equal(state.adv.circular_track_slots.some((slot) => slot.renderer === 'ticks'), false, preset);
  }

  const plain = createState();
  const preset = createState();
  createCircularTrackSlotEditor({ state: withDrawings(plain) }).resetCircularTrackSlotsFromSimpleControls();
  createCircularTrackSlotEditor({ state: withDrawings(preset) }).resetCircularTrackSlotsToPreset(preset.form.track_type);
  assert.deepEqual(plain.adv.circular_track_slots, preset.adv.circular_track_slots);
  assert.equal(plain.adv.circular_track_slots_axis_index, preset.adv.circular_track_slots_axis_index);
});

// OV-529: a renderer name from a hand-edited Session reads no inherited member of the renderer tables.
test('a renderer named like an inherited member keeps its own name as the label', () => {
  const editor = createCircularTrackSlotEditor({ state: withDrawings(createState()) });

  for (const renderer of ['constructor', '__proto__', 'toString']) {
    assert.equal(editor.circularTrackRendererLabel(renderer), renderer);
    assert.equal(editor.circularTrackSlotDisplayLabel({ ...spacerSlot(), renderer }), renderer);
  }
  assert.equal(editor.circularTrackRendererLabel('depth'), 'Depth');
  assert.equal(editor.circularTrackSlotDisplayLabel(spacerSlot()), 'Spacer');
});

const annotationSlot = (side = 'inside') => ({
  id: 'annotations',
  renderer: 'annotations',
  enabled: true,
  side,
  z: 0,
  params: { set_id: 'plastome_regions' }
});

// OV-1263 (audit G5b-03): the default stack with an Annotations row added at
// its end and set to Overlay, every other row inside the Axis. Outside axis
// moves that row across the Axis (before, the select snapped back to Inside);
// the other rows keep their sides, and Inside axis brings the row back. The
// flip of the other rows came from Generate's Axis index; the Playwright case
// in track-layout-transition.playwright.spec.js pins it.
test('an Annotations row placement moves only that row across the Axis (OV-1263)', () => {
  const inside = (slot) => ({ ...slot, side: 'inside' });
  const state = createState({
    slots: [
      inside({ ...featureSlot(), params: { lane_direction: 'inside' } }),
      inside(tickSlot()),
      inside(gcSlot()),
      { id: 'gc_skew', renderer: 'dinucleotide_skew', enabled: true, side: 'inside', z: 0, params: {} },
      { ...annotationSlot('overlay'), params: { set_id: 'plastome_regions', to_slot: 'gc_skew' } }
    ],
    axisIndex: 0
  });
  state.annotationSets = [{ id: 'plastome_regions' }];
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  const others = () => state.adv.circular_track_slots
    .filter((slot) => slot.id !== 'annotations').map((slot) => [slot.id, slot.side, slot.params.lane_direction]);
  const before = others();
  const row = state.adv.circular_track_slots[4];

  editor.updateCircularTrackSlotPlacement(row, 'outside');
  assert.deepEqual(others(), before);
  assert.equal(state.adv.circular_track_slots.find((slot) => slot.id === 'annotations').side, 'outside');

  editor.updateCircularTrackSlotPlacement(row, 'inside');
  assert.deepEqual(others(), before);
  assert.equal(state.adv.circular_track_slots.find((slot) => slot.id === 'annotations').side, 'inside');

  // A refused change (a Session operation runs) writes nothing.
  const settled = structuredClone(state.adv.circular_track_slots);
  const settledAxis = state.adv.circular_track_slots_axis_index;
  state.sessionOperationAvailability = () => 'busy';
  assert.equal(editor.updateCircularTrackSlotPlacement(row, 'outside'), 'busy');
  assert.deepEqual(state.adv.circular_track_slots, settled);
  assert.equal(state.adv.circular_track_slots_axis_index, settledAxis);
});

// OV-1254: the disabled Add track names the missing Region annotation set,
// not a Features-row conflict that Circular does not have.
test('Add track names the missing source of the chosen renderer (OV-1254)', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  assert.equal(editor.circularTrackAddBlockReason('annotations'), 'Import or add a Region annotation set in Region Annotations first');
  assert.equal(editor.circularTrackAddBlockReason('depth'), 'Load a Depth TSV first');
  assert.equal(editor.circularTrackAddBlockReason('spacer'), '');
  assert.equal(editor.canAddCircularTrackRenderer('annotations'), false);
  state.annotationSets = [{ id: 'plastome_regions' }];
  assert.equal(editor.circularTrackAddBlockReason('annotations'), '');
  assert.equal(editor.canAddCircularTrackRenderer('annotations'), true);
});

// D-51: the slot editor shows each row issue as one sentence in the UI language.
test('Circular editor shows a row issue in the UI language (D-51)', async () => {
  const { setUiLocale } = await import(pathToFileURL(join(tempRoot, 'utils', 'i18n.js')));
  const state = createState({ slots: [featureSlot(), { ...spacerSlot(), id: '' }, gcSlot()] });
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  assert.equal(editor.circularTrackSlotIssue(state.adv.circular_track_slots[1], 1), 'This row needs a non-empty ID.');
  await setUiLocale('ja');
  try {
    assert.equal(editor.circularTrackSlotIssue(state.adv.circular_track_slots[1], 1), 'この行には空でないIDが必要です。');
  } finally {
    await setUiLocale('en');
  }
});

test('OV-1656: the Circular gap writer stores one text, so `1px` is `1` before Save and after Load', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state: withDrawings(state) });
  const slot = state.adv.circular_track_slots[0];
  for (const field of /** @type {const} */ (['inner_gap_px', 'outer_gap_px'])) {
    editor.updateCircularTrackSlotGap(slot, field, ' 1px ');
    assert.strictEqual(slot[field], '1', field);
    editor.updateCircularTrackSlotGap(slot, field, '2.5 PX');
    assert.strictEqual(slot[field], '2.5', field);
    editor.updateCircularTrackSlotGap(slot, field, '10em');
    assert.strictEqual(slot[field], '10em', 'an invalid draft stays visible');
    editor.updateCircularTrackSlotGap(slot, field, '');
    assert.strictEqual(slot[field], null, 'blank is Auto');
  }
  editor.updateCircularTrackSlotGap(slot, 'width', '3px');
  assert.notEqual(slot.width, '3', 'only the gap fields');
});
