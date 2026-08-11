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
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

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

test('managed Depth insertion preserves the existing Axis sides', () => {
  const state = createState({
    showDepth: true,
    depthFiles: [{ name: 'depth.tsv', size: 1, lastModified: 0 }]
  });
  const editor = createCircularTrackSlotEditor({ state });

  editor.ensureCircularTrackDepthSlot();

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
    depthFiles: [
      { name: 'depth-a.tsv', size: 1, lastModified: 0 },
      { name: 'depth-b.tsv', size: 1, lastModified: 0 }
    ]
  });
  const editor = createCircularTrackSlotEditor({ state });

  editor.ensureCircularTrackDepthSlot();

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
  const beforeState = createState({
    slots: [depthSlot('outside'), featureSlot(), spacerSlot(), gcSlot()],
    axisIndex: 3
  });
  createCircularTrackSlotEditor({ state: beforeState }).ensureCircularTrackDepthSlot();
  assert.equal(beforeState.adv.circular_track_slots_axis_index, 2);
  assert.deepEqual(slotSides(beforeState), {
    features: 'outside',
    outer_gap: 'outside',
    gc_content: 'inside'
  });

  const afterState = createState({
    slots: [featureSlot(), spacerSlot(), gcSlot(), depthSlot()],
    axisIndex: 2
  });
  createCircularTrackSlotEditor({ state: afterState }).ensureCircularTrackDepthSlot();
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
