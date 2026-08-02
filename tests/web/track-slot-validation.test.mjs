import assert from 'node:assert/strict';
import { cp, mkdtemp, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import test from 'node:test';

const repoRoot = process.cwd();
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-track-slot-validation-'));
await cp(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app'),
  join(tempRoot, 'app'),
  { recursive: true }
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}', 'utf8');

const {
  CustomTrackPlanValidationError,
  assertValidCustomTrackPlan,
  validateCustomTrackPlan
} = await import(
  pathToFileURL(join(tempRoot, 'app', 'track-slot-validation.js'))
);

const feature = (id = 'features', overrides = {}) => ({
  id,
  renderer: 'features',
  enabled: true,
  side: 'overlay',
  z: 0,
  params: {},
  ...overrides
});

const depth = (id = 'depth', trackIndex = 0, overrides = {}) => ({
  id,
  renderer: 'depth',
  enabled: true,
  side: 'below',
  z: 0,
  params: { track_index: trackIndex },
  ...overrides
});

const annotation = (id = 'annotations', overrides = {}) => ({
  id,
  renderer: 'annotations',
  enabled: true,
  side: 'above',
  z: 1,
  params: { set_id: 'review', layer: 'foreground' },
  ...overrides
});

const validate = (overrides = {}) => validateCustomTrackPlan({
  mode: 'linear',
  slots: [feature()],
  axisIndex: 0,
  trackType: 'middle',
  depthTrackCount: 0,
  annotationSetIds: [],
  visibleFeatureUnderlays: false,
  conservationSeries: [],
  ...overrides
});

const rowCodes = (plan, index) => (
  (plan.rowIssues.get(index) || []).map((issue) => issue.code)
);

test('returns enabled projection and rebases the draft Axis boundary', () => {
  const slots = [
    depth('above', 0, { side: 'above' }),
    depth('disabled-before-axis', 0, { enabled: false, side: 'above' }),
    feature(),
    depth('disabled-after-axis', 0, { enabled: false }),
    depth('below', 0)
  ];
  const plan = validate({
    slots,
    axisIndex: 2,
    depthTrackCount: 1
  });

  assert.deepEqual(plan.enabledSlots.map((slot) => slot.id), ['above', 'features', 'below']);
  assert.equal(plan.emittedAxisIndex, 1);
  assert.equal(plan.rowIssues.size, 0);
  assert.deepEqual(plan.globalIssues, []);
  assert.notEqual(plan.draftSlots[0], slots[0], 'validation projection must be detached');
});

test('handles all Axis boundary edge cases after disabled-row filtering', () => {
  const slots = [
    depth('a', 0, { enabled: false, side: 'above' }),
    depth('b', 0, { enabled: false, side: 'above' }),
    depth('c', 0)
  ];
  assert.equal(validate({ slots, axisIndex: 0, depthTrackCount: 1 }).emittedAxisIndex, 0);
  assert.equal(validate({ slots, axisIndex: 2, depthTrackCount: 1 }).emittedAxisIndex, 0);
  assert.equal(validate({ slots, axisIndex: 3, depthTrackCount: 1 }).emittedAxisIndex, 1);
  assert.equal(validate({
    slots: slots.map((slot) => ({ ...slot, enabled: false })),
    axisIndex: 3,
    depthTrackCount: 1
  }).emittedAxisIndex, 0);

  const invalid = validate({ slots, axisIndex: 4, depthTrackCount: 1 });
  assert.deepEqual(invalid.globalIssues.map((issue) => issue.code), ['axis_out_of_range']);
  assert.equal(invalid.emittedAxisIndex, null);
});

test('reports duplicate and empty enabled IDs without rewriting the draft', () => {
  const slots = [feature(''), depth('', 0), feature('same'), depth('same', 0)];
  const plan = validate({ slots, depthTrackCount: 1 });
  assert.ok(rowCodes(plan, 0).includes('id_required'));
  assert.ok(rowCodes(plan, 1).includes('id_required'));
  assert.ok(rowCodes(plan, 2).includes('id_duplicate'));
  assert.ok(rowCodes(plan, 3).includes('id_duplicate'));
  assert.equal(plan.draftSlots[0].id, '');
});

test('limits Linear to one enabled Features row and ignores a disabled duplicate', () => {
  const invalid = validate({ slots: [feature('a'), feature('b')] });
  assert.ok(rowCodes(invalid, 1).includes('features_multiple'));

  const valid = validate({
    slots: [feature('a'), feature('b', { enabled: false })]
  });
  assert.equal(valid.rowIssues.size, 0);
});

test('requires exactly one Circular Features row only for visible feature underlays', () => {
  const duplicate = validate({
    mode: 'circular',
    slots: [feature('a', { side: 'inside' }), feature('b', { side: 'outside' })],
    axisIndex: 1,
    visibleFeatureUnderlays: 1
  });
  assert.deepEqual(
    duplicate.globalIssues.map((issue) => issue.code),
    ['feature_underlay_features_count']
  );

  const missing = validate({
    mode: 'circular',
    slots: [{ id: 'ticks', renderer: 'ticks', enabled: true, side: 'inside', z: 0, params: {} }],
    axisIndex: 0,
    visibleFeatureUnderlays: true
  });
  assert.deepEqual(
    missing.globalIssues.map((issue) => issue.code),
    ['feature_underlay_features_count']
  );

  const noUnderlay = validate({
    mode: 'circular',
    slots: [feature('a', { side: 'inside' }), feature('b', { side: 'outside' })],
    axisIndex: 1,
    visibleFeatureUnderlays: false
  });
  assert.equal(noUnderlay.globalIssues.length, 0);
});

test('rejects nonpositive Circular radius and width before request serialization', () => {
  for (const geometry of [
    { radius: 0 },
    { radius: '-0.5' },
    { width: 0 },
    { width: '-2px' },
    { width: { value: 0, unit: 'px' } }
  ]) {
    const plan = validate({
      mode: 'circular',
      slots: [feature('features', { side: 'inside', ...geometry })],
      axisIndex: 0
    });
    assert.ok(
      rowCodes(plan, 0).includes('geometry_invalid'),
      JSON.stringify(geometry)
    );
  }

  const valid = validate({
    mode: 'circular',
    slots: [feature('features', {
      side: 'inside',
      radius: { value: 0.75, unit: 'factor' },
      width: '12px'
    })],
    axisIndex: 0
  });
  assert.equal(valid.rowIssues.size, 0);
});

test('validates Depth binding and allows disabling an invalid row', () => {
  const missing = validate({ slots: [depth()], depthTrackCount: 0 });
  assert.ok(rowCodes(missing, 0).includes('depth_source_missing'));

  const outOfRange = validate({ slots: [depth('depth', 2)], depthTrackCount: 1 });
  assert.ok(rowCodes(outOfRange, 0).includes('depth_track_index_range'));

  const malformed = validate({
    slots: [depth('depth', 0, { params: { track_index: 1.5 } })],
    depthTrackCount: 2
  });
  assert.ok(rowCodes(malformed, 0).includes('depth_track_index'));

  const disabled = validate({
    slots: [depth('depth', 2, { enabled: false })],
    depthTrackCount: 0
  });
  assert.equal(disabled.rowIssues.size, 0);
});

test('validates Annotation set inventory and advanced params', () => {
  const plan = validate({
    slots: [
      annotation('blank', { params: { set_id: '' } }),
      annotation('unknown', { params: { set_id: 'missing' } }),
      annotation('options', {
        params: {
          set_id: 'review',
          marks: ['line', 'mystery'],
          lane_gap_px: -1,
          padding_px: 'nan',
          cover_anchor: 'yes'
        }
      })
    ],
    annotationSetIds: ['review']
  });
  assert.ok(rowCodes(plan, 0).includes('annotation_set_required'));
  assert.ok(rowCodes(plan, 1).includes('annotation_set_unknown'));
  assert.ok(rowCodes(plan, 2).includes('annotation_marks'));
  assert.ok(rowCodes(plan, 2).includes('annotation_lane_gap'));
  assert.ok(rowCodes(plan, 2).includes('annotation_padding'));
  assert.ok(rowCodes(plan, 2).includes('annotation_cover_anchor'));
});

test('validates overlay anchors and does not treat Spacer, Ticks, disabled, or Annotation rows as drawable', () => {
  const anchorCases = [
    ['missing', [], 'annotation_anchor_unknown'],
    ['disabled', [feature('disabled', { enabled: false })], 'annotation_anchor_unknown'],
    ['annotation', [annotation('annotation')], 'annotation_anchor_ineligible'],
    ['ticks', [{ id: 'ticks', renderer: 'ticks', enabled: true, side: 'below', z: 0, params: {} }], 'annotation_anchor_ineligible'],
    ['spacer', [{ id: 'spacer', renderer: 'spacer', enabled: true, side: 'below', z: 0, params: {} }], 'annotation_anchor_ineligible']
  ];

  for (const [anchorId, anchors, code] of anchorCases) {
    const plan = validate({
      slots: [
        ...anchors,
        annotation('overlay', {
          side: 'overlay',
          z: 2,
          params: {
            set_id: 'review',
            anchor_slot: anchorId,
            layer: 'foreground'
          }
        })
      ],
      annotationSetIds: ['review']
    });
    assert.ok(
      rowCodes(plan, anchors.length).includes(code),
      `${anchorId} should report ${code}`
    );
  }
});

test('accepts only source-bound managed Circular conservation rows', () => {
  const unmanaged = validate({
    mode: 'circular',
    slots: [{
      id: 'comparison',
      renderer: 'sequence_conservation',
      enabled: true,
      side: 'inside',
      z: 0,
      params: { source_index: 0, series_key: 'source-a' }
    }],
    axisIndex: 0,
    conservationSeries: [{ sourceKey: 'source-a', sourceIndex: 0 }]
  });
  assert.ok(rowCodes(unmanaged, 0).includes('conservation_unmanaged'));

  const missing = validate({
    mode: 'circular',
    slots: [{
      id: 'comparison',
      renderer: 'sequence_conservation',
      enabled: true,
      side: 'inside',
      z: 0,
      params: {
        managed: 'circular_conservation',
        source_index: 2,
        series_key: 'missing'
      }
    }],
    axisIndex: 0,
    conservationSeries: [{ sourceKey: 'source-a', sourceIndex: 0 }]
  });
  assert.ok(rowCodes(missing, 0).includes('conservation_source_missing'));

  const mismatchedBinding = validate({
    mode: 'circular',
    slots: [{
      id: 'comparison',
      renderer: 'sequence_conservation',
      enabled: true,
      side: 'inside',
      z: 0,
      params: {
        managed: 'circular_conservation',
        source_index: 1,
        series_key: 'source-a'
      }
    }],
    axisIndex: 0,
    conservationSeries: [
      { sourceKey: 'source-a', sourceIndex: 0, orderIndex: 0 },
      { sourceKey: 'source-b', sourceIndex: 1, orderIndex: 1 }
    ]
  });
  assert.ok(rowCodes(mismatchedBinding, 0).includes('conservation_source_missing'));

  const valid = validate({
    mode: 'circular',
    slots: [{
      id: 'comparison',
      renderer: 'sequence_conservation',
      enabled: true,
      side: 'inside',
      z: 0,
      params: {
        managed: 'circular_conservation',
        source_index: 0,
        series_key: 'source-a'
      }
    }],
    axisIndex: 0,
    conservationSeries: [{ sourceKey: 'source-a', sourceIndex: 0 }]
  });
  assert.equal(valid.rowIssues.size, 0);
});

test('reports renderer-specific params invalidated by a renderer change', () => {
  const plan = validate({
    mode: 'circular',
    slots: [
      {
        id: 'changed',
        renderer: 'ticks',
        enabled: true,
        side: 'inside',
        z: 0,
        params: { anchor_slot: 'features', axis: true }
      },
      {
        id: 'changed_skew',
        renderer: 'spacer',
        enabled: true,
        side: 'inside',
        z: 0,
        params: { nt: 'GC', positive_color: '#112233' }
      }
    ],
    axisIndex: 0
  });
  assert.ok(rowCodes(plan, 0).includes('ticks_obsolete_param'));
  assert.ok(rowCodes(plan, 0).includes('renderer_param_mismatch'));
  assert.ok(rowCodes(plan, 1).includes('renderer_param_mismatch'));
});

test('throws the typed summary only at an explicit boundary', () => {
  const plan = validate({ slots: [depth()], depthTrackCount: 0 });
  assert.throws(
    () => assertValidCustomTrackPlan(plan),
    (error) => (
      error instanceof CustomTrackPlanValidationError &&
      error.issues.some((issue) => issue.field === 'params.track_index') &&
      /Depth/.test(error.message)
    )
  );
});
