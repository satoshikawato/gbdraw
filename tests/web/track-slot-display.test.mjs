import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-track-slot-display-'));
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

const { findTrackSlotGeometry } = await import(
  pathToFileURL(join(tempRoot, 'app', 'track-slot-display.js'))
);
const { createLinearTrackSlotEditor } = await import(
  pathToFileURL(join(tempRoot, 'app', 'linear-track-slots.js'))
);
const { createCircularTrackSlotEditor } = await import(
  pathToFileURL(join(tempRoot, 'app', 'circular-track-slots.js'))
);

test('resolves shifted Linear public geometry by slot ID', () => {
  const geometry = {
    records: [{
      resultIndex: 0,
      recordIndex: 0,
      slots: [
        {
          slotIndex: 0,
          slotId: '__gbdraw_auto_feature_underlay_slot__',
          renderer: 'annotations',
          height: 14
        },
        { slotIndex: 1, slotId: 'depth_3', renderer: 'depth', height: 32 },
        { slotIndex: 2, slotId: 'features', renderer: 'features', height: 44 }
      ]
    }]
  };

  assert.deepEqual(
    findTrackSlotGeometry({ geometry, slotIndex: 0, slotId: 'depth_3' }),
    geometry.records[0].slots[1]
  );
  assert.deepEqual(
    findTrackSlotGeometry({ geometry, slotIndex: 1, slotId: 'features' }),
    geometry.records[0].slots[2]
  );
});

test('resolves shifted Circular public geometry by slot ID', () => {
  const geometry = {
    records: [{
      resultIndex: 1,
      recordIndex: 2,
      slots: [
        {
          slotIndex: 0,
          slotId: '__gbdraw_auto_feature_underlay_slot__',
          renderer: 'annotations',
          width: 11,
          radius: 1
        },
        { slotIndex: 1, slotId: 'features', renderer: 'features', width: 18, radius: 0.9 },
        { slotIndex: 2, slotId: 'ticks', renderer: 'ticks', width: 0, radius: 1 },
        { slotIndex: 3, slotId: 'depth', renderer: 'depth', width: 20, radius: 0.7 },
        { slotIndex: 4, slotId: 'gc_content', renderer: 'dinucleotide_content', width: 24, radius: 0.5 }
      ]
    }]
  };

  for (const [publicIndex, slotId, resolvedIndex] of [
    [0, 'features', 1],
    [1, 'ticks', 2],
    [2, 'depth', 3],
    [3, 'gc_content', 4]
  ]) {
    assert.deepEqual(
      findTrackSlotGeometry({
        geometry,
        resultIndex: 1,
        recordIndex: 2,
        slotIndex: publicIndex,
        slotId
      }),
      geometry.records[0].slots[resolvedIndex]
    );
  }
});

test('falls back to slot index for older geometry without slot IDs', () => {
  const geometry = {
    records: [{
      resultIndex: 0,
      recordIndex: 0,
      slots: [
        { slotIndex: 0, renderer: 'features', height: 41 },
        { slotIndex: 1, renderer: 'depth', height: 17 }
      ]
    }]
  };

  assert.deepEqual(
    findTrackSlotGeometry({ geometry, slotIndex: 1, slotId: 'depth' }),
    geometry.records[0].slots[1]
  );
  assert.deepEqual(
    findTrackSlotGeometry({ geometry, slotIndex: 0 }),
    geometry.records[0].slots[0]
  );
});

test('Linear placement actions preserve the Axis boundary and row sides', () => {
  const state = {
    adv: {
      nt: 'GC',
      linear_track_slots_enabled: true,
      linear_track_slots_axis_index: 1,
      linear_track_slots: [
        { id: 'gc_content', renderer: 'dinucleotide_content', side: 'above', params: { nt: 'GC' } },
        { id: 'features', renderer: 'features', side: 'overlay', params: {} },
        { id: 'gc_skew', renderer: 'dinucleotide_skew', side: 'below', params: { nt: 'GC' } }
      ]
    },
    form: {
      linear_track_layout: 'middle',
      show_depth: false,
      show_gc: true,
      show_skew: true
    },
    linearSeqs: []
  };
  const editor = createLinearTrackSlotEditor({ state });
  const order = () => state.adv.linear_track_slots.map((slot) => slot.id);
  const sides = () => state.adv.linear_track_slots.map((slot) => slot.side);

  assert.equal(editor.canMoveLinearTrackSlot(0, 1), false);
  editor.moveLinearTrackSlot(0, 1);
  assert.deepEqual(order(), ['gc_content', 'features', 'gc_skew']);

  editor.moveLinearTrackSlotAbove(2);
  assert.deepEqual(order(), ['gc_content', 'gc_skew', 'features']);
  assert.equal(state.adv.linear_track_slots_axis_index, 2);
  assert.deepEqual(sides(), ['above', 'above', 'overlay']);

  editor.moveLinearTrackSlotBelow(0);
  assert.deepEqual(order(), ['gc_skew', 'features', 'gc_content']);
  assert.equal(state.adv.linear_track_slots_axis_index, 1);
  assert.deepEqual(sides(), ['above', 'overlay', 'below']);

  let feature = state.adv.linear_track_slots.find((slot) => slot.id === 'features');
  editor.updateLinearTrackSlotPlacement(feature, 'below');
  feature = state.adv.linear_track_slots.find((slot) => slot.id === 'features');
  assert.equal(feature.side, 'below');
  editor.updateLinearTrackSlotPlacement(feature, 'overlay');
  assert.equal(
    state.adv.linear_track_slots.find((slot) => slot.id === 'features').side,
    'overlay'
  );
});

test('Linear editor requests resolved geometry with each public slot ID', () => {
  const linearSlots = [
    {
      id: 'depth_3', renderer: 'depth', enabled: true, side: 'above',
      height: null, spacing: null, z: 0, params: { track_index: 0 }
    },
    {
      id: 'features', renderer: 'features', enabled: true, side: 'overlay',
      height: null, spacing: null, z: 0, params: {}
    }
  ];
  const linearGeometry = {
    mode: 'linear',
    records: [{
      resultIndex: 0,
      recordIndex: 0,
      slots: [
        {
          slotIndex: 0,
          slotId: '__gbdraw_auto_feature_underlay_slot__',
          heightPx: 14
        },
        { slotIndex: 1, slotId: 'depth_3', heightPx: 32 },
        { slotIndex: 2, slotId: 'features', heightPx: 44 }
      ]
    }]
  };
  const linearEditor = createLinearTrackSlotEditor({
    state: {
      form: { linear_track_layout: 'middle', show_depth: true },
      adv: {
        linear_track_slots: linearSlots,
        linear_track_slots_axis_index: 1,
        depth_tracks: [{ height: null }],
        depth_color: '#4A90E2',
        depth_height: null,
        gc_height: null,
        nt: 'GC'
      },
      files: { linearSeqs: [] },
      annotationSets: [],
      selectedResultIndex: { value: 0 },
      trackSlotResolvedGeometry: { value: linearGeometry }
    }
  });
  assert.equal(
    linearEditor.linearTrackSlotGeometryAutoText(linearSlots[0], 0, 'height'),
    '32 px (auto)'
  );
  assert.equal(
    linearEditor.linearTrackSlotGeometryAutoText(linearSlots[1], 1, 'height'),
    '44 px (auto; varies by record)'
  );
});

test('Circular editor requests resolved geometry with each public slot ID', () => {
  const circularSlots = [
    { id: 'features', renderer: 'features', enabled: true, side: 'inside', params: {} },
    { id: 'ticks', renderer: 'ticks', enabled: true, side: 'inside', params: {} },
    { id: 'depth', renderer: 'depth', enabled: true, side: 'inside', params: { track_index: 0 } },
    { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'inside', params: { nt: 'GC' } }
  ];
  const circularGeometry = {
    mode: 'circular',
    records: [{
      resultIndex: 0,
      recordIndex: 0,
      slots: [
        {
          slotIndex: 0,
          slotId: '__gbdraw_auto_feature_underlay_slot__',
          widthPx: 11,
          radiusFactor: 1
        },
        { slotIndex: 1, slotId: 'features', widthPx: 18, radiusFactor: 0.9 },
        { slotIndex: 2, slotId: 'ticks', widthPx: 0, radiusFactor: 1 },
        { slotIndex: 3, slotId: 'depth', widthPx: 20, radiusFactor: 0.7 },
        { slotIndex: 4, slotId: 'gc_content', widthPx: 24, radiusFactor: 0.5 }
      ]
    }]
  };
  const circularEditor = createCircularTrackSlotEditor({
    state: {
      mode: { value: 'circular' },
      form: {
        track_type: 'tuckin', show_depth: true, suppress_gc: false,
        suppress_skew: true, show_scale: true, separate_strands: true
      },
      adv: {
        circular_track_slots: circularSlots,
        circular_track_slots_axis_index: 0,
        nt: 'GC',
        features: ['repeat_region'],
        feature_shapes: { repeat_region: 'underlay' },
        depth_tracks: [],
        feature_width_circular: null,
        depth_width_circular: null,
        gc_content_width_circular: null,
        gc_content_radius_circular: null,
        gc_skew_width_circular: null,
        gc_skew_radius_circular: null
      },
      files: { c_depth: [] },
      circularConservation: { enabled: false, series: [] },
      circularRecordList: { value: [] },
      annotationSets: [],
      selectedResultIndex: { value: 0 },
      trackSlotResolvedGeometry: { value: circularGeometry }
    }
  });
  for (const [slotIndex, expected] of [
    [0, '18 px (auto)'],
    [1, '0 px (auto)'],
    [2, '20 px (auto)'],
    [3, '24 px (auto)']
  ]) {
    assert.equal(
      circularEditor.circularTrackSlotGeometryAutoText(
        circularSlots[slotIndex],
        slotIndex,
        'width'
      ),
      expected
    );
  }
});
