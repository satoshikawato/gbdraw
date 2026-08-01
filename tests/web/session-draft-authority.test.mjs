import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};
globalThis.File = class File extends Blob {
  constructor(parts, name, options = {}) {
    super(parts, options);
    this.name = String(name || 'file');
    this.lastModified = options.lastModified ?? Date.now();
  }
};
globalThis.alert = () => {};

const {
  importSession,
  overlayCurrentWriterDraftConfig,
  validateCurrentWriterActiveDraft
} = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');

const canonicalFeature = {
  id: 'features',
  renderer: 'features',
  enabled: true,
  side: 'inside',
  width: null,
  radius: null,
  inner_gap_px: null,
  outer_gap_px: null,
  z: 0,
  params: { lane_direction: 'split' }
};
const disabledDraft = {
  id: 'disabled-draft',
  renderer: 'depth',
  enabled: false,
  side: 'outside',
  width: '27px',
  radius: '1.2',
  inner_gap_px: '4px',
  outer_gap_px: '5px',
  z: 3,
  params: { track_index: 99, nested: { keep: true } }
};
const projectedConfig = {
  form: { track_type: 'tuckin' },
  adv: {
    nt: 'GC',
    circular_track_slots_enabled: true,
    circular_track_slots: [canonicalFeature],
    circular_track_slots_axis_index: 1,
    linear_track_slots_enabled: false,
    linear_track_slots: [],
    linear_track_slots_axis_index: null
  }
};
const storedConfig = {
  form: { track_type: 'tuckin' },
  adv: {
    ...projectedConfig.adv,
    circular_track_slots: [disabledDraft, canonicalFeature],
    circular_track_slots_axis_index: 2,
    linear_track_slots_enabled: false,
    linear_track_slots: [{
      id: 'inactive-linear',
      renderer: 'spacer',
      enabled: false,
      side: 'below',
      height: '8px',
      spacing: '2px',
      z: 0,
      params: {}
    }],
    linear_track_slots_axis_index: 1,
    feature_width_circular: 19,
    depth_width_circular: 23
  }
};

assert.doesNotThrow(() => validateCurrentWriterActiveDraft({
  mode: 'circular',
  projectedConfig,
  storedConfig
}));
const restored = overlayCurrentWriterDraftConfig(projectedConfig, storedConfig);
assert.deepEqual(restored.adv.circular_track_slots, storedConfig.adv.circular_track_slots);
assert.deepEqual(restored.adv.linear_track_slots, storedConfig.adv.linear_track_slots);
assert.equal(restored.adv.circular_track_slots_axis_index, 2);
assert.equal(restored.adv.feature_width_circular, 19);
assert.equal(restored.adv.depth_width_circular, 23);

const mismatched = structuredClone(storedConfig);
mismatched.adv.circular_track_slots[1].width = '12px';
assert.throws(
  () => validateCurrentWriterActiveDraft({
    mode: 'circular',
    projectedConfig,
    storedConfig: mismatched
  }),
  /does not match the committed render request/
);

const inactive = structuredClone(mismatched);
inactive.adv.circular_track_slots_enabled = false;
assert.doesNotThrow(() => validateCurrentWriterActiveDraft({
  mode: 'circular',
  projectedConfig,
  storedConfig: inactive
}));

const divergentSession = JSON.parse(await readFile(
  'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json',
  'utf8'
));
const committedComparisonCount = divergentSession.renderRequest.comparisons.length;
assert.ok(committedComparisonCount > 0);
divergentSession.config.linearComparisonPlan = {
  mode: 'none',
  defaultSource: 'losat',
  edges: []
};
const importEvent = {
  target: {
    files: [new Blob([JSON.stringify(divergentSession)], { type: 'application/json' })],
    value: 'selected'
  }
};

const imported = await importSession(importEvent);

assert.equal(imported.status, 'ok');
assert.equal(state.linearComparisonPlan.mode, 'none');
assert.deepEqual(state.linearComparisonPlan.edges, []);
assert.equal(
  imported.data.renderRequest.comparisons.length,
  committedComparisonCount,
  'the editable comparison draft must not replace the last committed render request'
);
