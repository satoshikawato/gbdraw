// Feature density slots in the Web's track stacks (#199 E4): the renderer names come from Python (D-29), a
// `feature_density` row keeps its renderer and names its track by id, the stacks the Web builds from the simple
// settings hold one row per track as Python's do (OV-3174), and a row that names no track of the drawing is a row
// issue (OV-3173 is the Load that refused such a row).
import assert from 'node:assert/strict';
import test from 'node:test';
import { MODE_PROFILE_DATA } from '../../gbdraw/web/js/mode-profiles.generated.js';
import {
  CIRCULAR_TRACK_RENDERERS, buildCircularTrackSlotPayload, createDefaultCircularTrackSlots, normalizeCircularTrackSlot,
  parseCircularTrackSlotSpec
} from '../../gbdraw/web/js/services/circular-track-slot-model.js';
import {
  LINEAR_TRACK_RENDERERS, buildLinearTrackSlotPayload, createDefaultLinearTrackSlots, normalizeLinearTrackSlots
} from '../../gbdraw/web/js/services/linear-track-slot-model.js';
import { validateCustomTrackPlan, validateTrackSlotBindingInvariants } from '../../gbdraw/web/js/services/track-slot-validation.js';
import { insertFeatureDensitySlots } from '../../gbdraw/web/js/services/feature-density-tracks.js';

test('the slot models and the slot checks take each mode\'s renderers from Python (D-29)', () => {
  assert.deepEqual([...CIRCULAR_TRACK_RENDERERS], MODE_PROFILE_DATA.modes.circular.trackRenderers);
  assert.deepEqual([...LINEAR_TRACK_RENDERERS], MODE_PROFILE_DATA.modes.linear.trackRenderers);
  for (const mode of ['circular', 'linear']) {
    assert.ok(MODE_PROFILE_DATA.modes[mode].trackRenderers.includes('feature_density'), mode);
    const plan = validateCustomTrackPlan({
      mode, slots: [{ id: 'features', renderer: 'features', enabled: true, side: mode === 'circular' ? null : 'overlay', z: 0, params: {} },
        { id: 'density_genes', renderer: 'feature_density', enabled: true, side: null, z: 0, params: { track: 'genes' } }],
      featureDensityTracks: [{ id: 'genes' }]
    });
    assert.deepEqual([...plan.rowIssues.values()].flat().map(({ code }) => code), [], mode);
  }
});

test('a feature density row keeps its renderer and its track id in both modes', () => {
  const circular = normalizeCircularTrackSlot({ id: 'density_genes', renderer: 'feature_density', enabled: true, params: { track: ' genes ' } });
  assert.equal(circular.renderer, 'feature_density');
  assert.deepEqual(circular.params, { track: 'genes' });
  assert.deepEqual(buildCircularTrackSlotPayload(circular).params, { track: 'genes' });
  assert.deepEqual(parseCircularTrackSlotSpec('density_ltr:feature_density@track=ltr,legend_label=LTR').params,
    { track: 'ltr', legend_label: 'LTR' });
  const [linear] = normalizeLinearTrackSlots([{ id: 'density_genes', renderer: 'feature_density', enabled: true, side: 'below', params: { track: 'genes' } }]);
  assert.equal(linear.renderer, 'feature_density');
  assert.deepEqual(buildLinearTrackSlotPayload(linear).params, { track: 'genes' });
});

test('a Load of a stack with feature density rows passes the binding check (OV-3173)', () => {
  for (const [mode, sides] of [['circular', ['inside', 'outside', 'overlay']], ['linear', ['above', 'below', 'overlay']]]) {
    assert.doesNotThrow(() => validateTrackSlotBindingInvariants([
      { id: 'features', renderer: 'features', enabled: true, side: 'overlay', z: 0, params: {} },
      { id: 'density_genes', renderer: 'feature_density', enabled: true, side: null, z: 0, params: { track: 'genes' } }
    ], {
      modeLabel: mode, layoutKind: mode, supportedRenderers: MODE_PROFILE_DATA.modes[mode].trackRenderers,
      supportedSides: sides, anchorlessRenderers: ['spacer'], depthTrackCount: 0
    }), mode);
  }
});

test('the stacks the Web builds hold one feature density row per track after the Features row (OV-3174)', () => {
  const circular = createDefaultCircularTrackSlots({ featureDensityTrackIds: ['genes', 'trna'] });
  assert.deepEqual(circular.map(({ id, renderer }) => `${id}:${renderer}`),
    ['features:features', 'density_genes:feature_density', 'density_trna:feature_density', 'ticks:ticks',
      'gc_content:dinucleotide_content', 'gc_skew:dinucleotide_skew']);
  assert.deepEqual(circular[1].params, { track: 'genes' });
  const linear = createDefaultLinearTrackSlots({ showGc: true, featureDensityTrackIds: ['genes'] });
  assert.deepEqual(linear.map(({ id, renderer, side }) => `${id}:${renderer}:${side}`),
    ['features:features:overlay', 'density_genes:feature_density:below', 'gc_content:dinucleotide_content:below']);
  assert.deepEqual(createDefaultCircularTrackSlots().map(({ id }) => id), ['features', 'ticks', 'gc_content', 'gc_skew']);
});

test('a feature density row must name one of the drawing\'s tracks', () => {
  const issues = (track, tracks) => [...validateCustomTrackPlan({
    mode: 'circular',
    slots: [{ id: 'features', renderer: 'features', enabled: true, side: null, z: 0, params: {} },
      { id: 'density', renderer: 'feature_density', enabled: true, side: null, z: 0, params: track === null ? {} : { track } }],
    featureDensityTracks: tracks
  }).rowIssues.values()].flat().map(({ code }) => code);
  assert.deepEqual(issues('genes', [{ id: 'genes' }]), []);
  assert.deepEqual(issues('repeats', [{ id: 'genes' }]), ['feature_density_track_unknown']);
  assert.deepEqual(issues(null, [{ id: 'genes' }]), ['feature_density_track_required']);
});

test('Use custom stack keeps the feature density tracks the simple stack drew (OV-3174)', () => {
  const rows = (...ids) => ids.map((id) => ({ id, renderer: id === 'features' ? 'features' : (id.startsWith('density') ? 'feature_density' : 'spacer') }));
  const make = (track) => ({ id: `density_${track}`, renderer: 'feature_density' });
  // After the Features row; a row outside the Axis keeps its side.
  const stack = rows('ticks', 'features', 'gc');
  assert.equal(insertFeatureDensitySlots(stack, 1, ['genes', 'ltr'], make), 1);
  assert.deepEqual(stack.map(({ id }) => id), ['ticks', 'features', 'density_genes', 'density_ltr', 'gc']);
  const outside = rows('features', 'ticks', 'gc');
  assert.equal(insertFeatureDensitySlots(outside, 2, ['genes'], make), 3);
  // A stack that draws a feature density track already is kept; no tracks add nothing; no Features row: the top.
  const kept = rows('features', 'density_genes');
  assert.equal(insertFeatureDensitySlots(kept, null, ['genes', 'ltr'], make), null);
  assert.deepEqual(kept.map(({ id }) => id), ['features', 'density_genes']);
  const top = rows('gc');
  insertFeatureDensitySlots(top, null, ['genes'], make);
  assert.deepEqual(top.map(({ id }) => id), ['density_genes', 'gc']);
});
