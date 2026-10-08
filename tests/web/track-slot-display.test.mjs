import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';
import { withDrawings } from './helpers/drawing-state.mjs';

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
// circular-track-slots.js reads the conservation-series row helpers.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'conservation-series.js'),
  join(tempRoot, 'services', 'conservation-series.js')
);
// linear-track-slots.js reads the current option values.
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'current-option-values.js'),
  join(tempRoot, 'services', 'current-option-values.js')
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');
// The track-slot leaves the slot editors import.
for (const leaf of ['depth-track-state.js', 'track-slot-display.js', 'track-slot-validation.js',
  'circular-track-slot-model.js', 'linear-track-slot-model.js']) {
  await cp(
    join(repoRoot, 'gbdraw', 'web', 'js', 'services', leaf),
    join(tempRoot, 'services', leaf)
  );
}

const { findTrackSlotGeometry, formatRadiusFactorAuto, tickAnchorRadiusFactor } = await import(
  pathToFileURL(join(tempRoot, 'services', 'track-slot-display.js'))
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
  const editor = createLinearTrackSlotEditor({ state: withDrawings(state) });
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
    state: withDrawings({
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
    })
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
    state: withDrawings({
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
    })
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

// TR-08: Python emits geometry only for rendered rows, so a slot index does not
// identify a row. Geometry is found by slot ID alone.
test('a row without rendered geometry does not resolve to another row', () => {
  const geometry = {
    records: [{
      resultIndex: 0,
      recordIndex: 0,
      slots: [
        { slotIndex: 0, slotId: 'features', renderer: 'features', widthPx: 60 },
        { slotIndex: 1, slotId: 'gc_content', renderer: 'gc_content', widthPx: 74.1 }
      ]
    }]
  };
  assert.equal(findTrackSlotGeometry({ geometry, slotId: 'ticks' }), null);
  assert.equal(findTrackSlotGeometry({ geometry, slotIndex: 1, slotId: 'ticks' }), null);
  assert.equal(findTrackSlotGeometry({ geometry, slotIndex: 1 }), null);
  assert.deepEqual(findTrackSlotGeometry({ geometry, slotId: 'gc_content' }), geometry.records[0].slots[1]);
});

// TR-08: a disabled row, including one that shares an ID with a rendered row,
// shows its estimate instead of resolved geometry.
test('disabled rows show their estimate in both editors', () => {
  const circularGeometry = {
    mode: 'circular',
    records: [{
      resultIndex: 0,
      recordIndex: 0,
      slots: [
        { slotIndex: 0, slotId: 'features', widthPx: 60, radiusFactor: 0.9 },
        { slotIndex: 1, slotId: 'gc_content', widthPx: 74.1, radiusFactor: 0.69 }
      ]
    }]
  };
  const circularState = (geometry) => ({
    mode: { value: 'circular' },
    form: {
      track_type: 'tuckin', show_depth: false, suppress_gc: false,
      suppress_skew: true, show_scale: true, separate_strands: true
    },
    adv: {
      circular_track_slots: [
        { id: 'features', renderer: 'features', enabled: true, side: 'inside', params: {} },
        { id: 'ticks', renderer: 'ticks', enabled: false, side: 'inside', params: {} },
        { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'inside', params: { nt: 'GC' } },
        { id: 'gc_content', renderer: 'dinucleotide_content', enabled: false, side: 'inside', params: { nt: 'GC' } }
      ],
      circular_track_slots_axis_index: 0,
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
    files: { c_depth: [] },
    circularConservation: { enabled: false, series: [] },
    circularRecordList: { value: [] },
    annotationSets: [],
    selectedResultIndex: { value: 0 },
    trackSlotResolvedGeometry: { value: geometry }
  });
  const resolved = createCircularTrackSlotEditor({ state: withDrawings(circularState(circularGeometry)) });
  const estimated = createCircularTrackSlotEditor({ state: withDrawings(circularState(null)) });
  // circularTrackSlots() lists row entries; the note reads the entry's slot.
  const circularText = (editor, index, field) => (
    editor.circularTrackSlotGeometryAutoText(editor.circularTrackSlots()[index].slot, index, field)
  );
  assert.equal(circularText(resolved, 2, 'width'), '74.1 px (auto)');
  for (const index of [1, 3]) {
    for (const field of ['width', 'radius']) {
      assert.equal(circularText(resolved, index, field), circularText(estimated, index, field), `${index} ${field}`);
    }
  }

  const linearSlots = [
    { id: 'features', renderer: 'features', enabled: true, side: 'overlay', height: null, spacing: null, z: 0, params: {} },
    { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'below', height: null, spacing: null, z: 0, params: { nt: 'GC' } },
    { id: 'gc_content', renderer: 'dinucleotide_content', enabled: false, side: 'below', height: null, spacing: null, z: 0, params: { nt: 'GC' } }
  ];
  const linearEditor = createLinearTrackSlotEditor({
    state: withDrawings({
      form: { linear_track_layout: 'middle', show_depth: false },
      adv: {
        linear_track_slots: linearSlots,
        linear_track_slots_axis_index: 0,
        depth_tracks: [],
        depth_height: null,
        gc_height: null,
        nt: 'GC'
      },
      linearSeqs: [],
      annotationSets: [],
      selectedResultIndex: { value: 0 },
      trackSlotResolvedGeometry: {
        value: {
          mode: 'linear',
          records: [{ resultIndex: 0, recordIndex: 0, slots: [{ slotIndex: 1, slotId: 'gc_content', heightPx: 77 }] }]
        }
      }
    })
  });
  assert.equal(linearEditor.linearTrackSlotGeometryAutoText(linearSlots[1], 1, 'height'), '77 px (auto)');
  assert.notEqual(linearEditor.linearTrackSlotGeometryAutoText(linearSlots[2], 2, 'height'), '77 px (auto)');
});

// GX-18: Python's resolver (HmmtDNA, Axis radius 390 px; ids `ticks` with the
// preset's 7.8 px width and `scale` with an Auto width) gives these band
// centres and anchors: [tick_label_layout, side, widthPx, radiusFactor, anchor].
const TICK_ANCHOR_VECTORS = [
  ['label_out_tick_in', 'inside', 7.8, 0.9274365651709401, 0.9374365651709402],
  ['label_out_tick_in', 'outside', 7.8, 1.0199999999999998, 1.03],
  ['label_in_tick_out', 'inside', 7.8, 0.9800000000000001, 0.9700000000000001],
  ['label_in_tick_out', 'outside', 7.8, 1.0199999999999998, 1.01],
  ['tick_only', 'inside', 7.8, 0.9700000000000002, 0.9800000000000001],
  ['tick_only', 'outside', 7.8, 1.0199999999999998, 1.01],
  ['label_only', 'inside', 7.8, 0.9800000000000001, 0.9800000000000001],
  ['label_only', 'outside', 7.8, 1.01, 1.01],
  ['label_out_tick_in', 'inside', 0.0, 0.9249365651709403, 0.9374365651709402],
  ['label_out_tick_in', 'outside', 0.0, 1.0225, 1.035],
  ['label_in_tick_out', 'inside', 0.0, 0.9775, 0.9650000000000001],
  ['label_in_tick_out', 'outside', 0.0, 1.0225, 1.01],
  ['tick_only', 'inside', 0.0, 0.9775, 0.9900000000000001],
  ['tick_only', 'outside', 0.0, 1.0225, 1.01],
  ['label_only', 'inside', 0.0, 0.9900000000000001, 0.9900000000000001],
  ['label_only', 'outside', 0.0, 1.01, 1.01]
];

test('a ticks row note converts the band centre to the anchor that r pins (GX-18)', () => {
  for (const [layout, side, widthPx, radiusFactor, anchor] of TICK_ANCHOR_VECTORS) {
    const actual = tickAnchorRadiusFactor({ radiusFactor, widthPx, side }, 390, layout);
    assert.ok(Math.abs(actual - anchor) < 1e-12, `${layout} ${side} ${widthPx}: ${actual} != ${anchor}`);
  }
  assert.equal(tickAnchorRadiusFactor(null, 390, 'tick_only'), null);
  assert.equal(tickAnchorRadiusFactor({ radiusFactor: 0.9, widthPx: 0 }, 0, 'tick_only'), null);
  // Two decimals, rounded away from the Axis: typing the note back keeps the
  // ticks off the neighbouring row (MG1655 Tuckin 0.855, Middle 0.9175).
  for (const [anchor, text] of [[0.855, '0.85 R (auto)'], [0.9175, '0.91 R (auto)'], [0.76, '0.76 R (auto)'],
    [0.7599999999999999, '0.76 R (auto)'], [1.0125, '1.02 R (auto)'], [1.01, '1.01 R (auto)'], [1, '1 R (auto)']]) {
    assert.equal(formatRadiusFactorAuto(anchor, { awayFromAxis: true }), text, String(anchor));
  }
});

const circularNoteEditor = (geometry) => {
  const slots = [
    { id: 'features', renderer: 'features', enabled: true, side: 'inside', params: { lane_direction: 'inside' } },
    { id: 'ticks', renderer: 'ticks', enabled: true, side: 'inside', params: { tick_label_layout: 'label_in_tick_out' } },
    { id: 'gc_content', renderer: 'dinucleotide_content', enabled: true, side: 'inside', params: { nt: 'GC' } }
  ];
  const editor = createCircularTrackSlotEditor({
    state: withDrawings({
      mode: { value: 'circular' },
      form: { track_type: 'tuckin', show_depth: false, suppress_gc: false, suppress_skew: true, show_scale: true, separate_strands: true },
      adv: {
        circular_track_slots: slots, circular_track_slots_axis_index: 0, nt: 'GC',
        features: ['CDS'], feature_shapes: { CDS: 'arrow' }, depth_tracks: [],
        feature_width_circular: null, depth_width_circular: null, gc_content_width_circular: null,
        gc_content_radius_circular: null, gc_skew_width_circular: null, gc_skew_radius_circular: null
      },
      files: { c_depth: [] },
      circularConservation: { enabled: false, series: [] },
      circularRecordList: { value: [] },
      annotationSets: [],
      selectedResultIndex: { value: 0 },
      trackSlotResolvedGeometry: { value: geometry }
    })
  });
  const note = (index, field) => editor.circularTrackSlotGeometryAutoText(slots[index], index, field);
  return Object.assign(note, { editor, slots });
};

// TK-15: before a render the note is an estimate and says so; a rendered row
// shows the measured value. An Auto ticks row is never 0 px wide.
test('Circular notes mark estimates and show measured values after a render (TK-15, GX-18)', () => {
  const estimated = circularNoteEditor(null);
  assert.equal(estimated(1, 'width'), '≈ 9.8 px (estimate)');
  for (const [index, field] of [[0, 'width'], [0, 'radius'], [1, 'radius'], [2, 'width'], [2, 'radius'], [2, 'inner_gap_px']]) {
    assert.match(estimated(index, field), /^≈ [0-9.]+ (?:px|R) \(estimate\)$/, `${index} ${field}`);
  }
  const rendered = circularNoteEditor({
    mode: 'circular',
    records: [{ resultIndex: 0, recordIndex: 0, axisRadiusPx: 390, slots: [
      { slotIndex: 0, slotId: 'features', renderer: 'features', side: 'inside', widthPx: 37, radiusFactor: 0.89, innerGapPx: 3.9, outerGapPx: 3.9 },
      { slotIndex: 1, slotId: 'ticks', renderer: 'ticks', side: 'inside', widthPx: 7.8, radiusFactor: 0.77, innerGapPx: 3.9, outerGapPx: 3.9 },
      { slotIndex: 2, slotId: 'gc_content', renderer: 'dinucleotide_content', side: 'inside', widthPx: 54.1, radiusFactor: 0.54, innerGapPx: 3.9, outerGapPx: 3.9 }
    ] }]
  });
  assert.equal(rendered(0, 'radius'), '0.89 R (auto)');
  assert.equal(rendered(1, 'width'), '7.8 px (auto)');
  // The ticks grow outward (label_in_tick_out) from the anchor 0.77 - 7.8 / 390 / 2.
  assert.equal(rendered(1, 'radius'), '0.76 R (auto)');
  assert.equal(rendered(2, 'width'), '54.1 px (auto)');
  assert.equal(rendered(2, 'inner_gap_px'), '3.9 px (auto)');
});

// TK-12: an invalid Width or Radius is reported once, by its own field; the
// row alert does not repeat it.
test('the row alert leaves an invalid Width or Radius to its field (TK-12)', () => {
  const { editor, slots } = circularNoteEditor(null);
  slots[2].width = { value: '0', unit: 'factor' };
  slots[0].radius = '0x10';
  assert.equal(editor.circularTrackSlotIssue(slots[2], 2), '');
  assert.equal(editor.circularTrackSlotIssue(slots[0], 0), '');
  slots[2].z = 'x';
  assert.match(editor.circularTrackSlotIssue(slots[2], 2), /z must be an integer/);
});
