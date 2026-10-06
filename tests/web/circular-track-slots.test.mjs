import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';

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
// circular-track-slots.js reports through the dependency-free wording owner.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'error-normalization.js'),
  join(tempRoot, 'services', 'error-normalization.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// The track-slot leaves the slot editors import.
for (const leaf of ['depth-track-state.js', 'track-slot-display.js', 'track-slot-validation.js']) {
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
  const editor = createCircularTrackSlotEditor({ state });

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
  const editor = createCircularTrackSlotEditor({ state });

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
  const editor = createCircularTrackSlotEditor({ state });

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
  const removeSource = (state) => createCircularTrackSlotEditor({ state }).changeCircularDepthSources(() => {
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
  const editor = createCircularTrackSlotEditor({ state });

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
  const editor = createCircularTrackSlotEditor({ state });

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
  createCircularTrackSlotEditor({ state: beforeState }).syncCircularConservationSlots();
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
  createCircularTrackSlotEditor({ state: afterState }).syncCircularConservationSlots();
  assert.equal(afterState.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(afterState), {
    features: 'outside',
    outer_gap: 'outside',
    gc_content: 'inside'
  });
});

test('Circular measure action updates only the owned scalar and preserves slot identity', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state });
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
  const editor = createCircularTrackSlotEditor({ state });
  const slot = state.adv.circular_track_slots[0];
  for (const field of ['width', 'radius']) {
    const before = slot[field];
    for (const scalar of [NaN, Infinity, -Infinity,
      { value: NaN, unit: 'px' }, { value: Infinity, unit: 'factor' },
      { value: -Infinity, unit: 'px' }]) {
      assert.throws(() => editor.updateCircularTrackSlotMeasure(slot, field, scalar), /positive finite px or factor scalar/);
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

// TR-06: Hide GC Content owns the visibility of its rows whether or not the
// custom stack is in use.
test('un-hiding GC while the custom stack is off restores the row the hide disabled', () => {
  const state = createState();
  const editor = createCircularTrackSlotEditor({ state });

  editor.setCircularGcSuppressed(true);
  editor.setCircularTrackSlotsEnabled(false);
  editor.setCircularGcSuppressed(false);
  editor.setCircularTrackSlotsEnabled(true);

  assert.deepEqual(gcRows(state), [{ id: 'gc_content', enabled: true, marked: false }]);
});

test('a GC row the user disabled stays disabled across Hide and un-hide', () => {
  const state = createState({
    slots: [featureSlot(), spacerSlot(), { ...gcSlot(), enabled: false }, { ...gcSlot(), id: 'gc_content_2' }]
  });
  const editor = createCircularTrackSlotEditor({ state });

  editor.setCircularGcSuppressed(true);
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
    createCircularTrackSlotEditor({ state }).resetCircularTrackSlotsToPreset(preset);
    assert.equal(state.adv.circular_track_slots.some((slot) => slot.renderer === 'ticks'), false, preset);
  }

  const plain = createState();
  const preset = createState();
  createCircularTrackSlotEditor({ state: plain }).resetCircularTrackSlotsFromSimpleControls();
  createCircularTrackSlotEditor({ state: preset }).resetCircularTrackSlotsToPreset(preset.form.track_type);
  assert.deepEqual(plain.adv.circular_track_slots, preset.adv.circular_track_slots);
  assert.equal(plain.adv.circular_track_slots_axis_index, preset.adv.circular_track_slots_axis_index);
});
