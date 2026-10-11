import assert from 'node:assert/strict';

assert.equal(globalThis.window, undefined);
assert.equal(globalThis.document, undefined);

const {
  CURRENT_WRITER_ADV_FIELDS,
  CURRENT_WRITER_FORM_FIELDS,
  createDefaultAdv,
  createDefaultForm,
  createDefaultLosat,
  createDefaultLosatExecution,
  holdsCliWriterConfig,
  validateCurrentWriterActiveConfig,
  validateImportedCircularTrackSlots
} = await import('../../gbdraw/web/js/services/session-active-config-contract.js');
const { normalizeUserFacingError } = await import('../../gbdraw/web/js/utils/error-normalization.js');
const { normalizeCurrentPairwiseMatchStyle } = await import(
  '../../gbdraw/web/js/services/current-option-values.js'
);

const storedConfig = {
  form: createDefaultForm(),
  adv: createDefaultAdv('circular')
};

assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
  mode: 'circular',
  storedConfig
}));
assert.deepEqual(CURRENT_WRITER_FORM_FIELDS, [
  ...Object.keys(createDefaultForm()),
  'legend'
]);
// `rich_feature_popup` is read from a Session 40-44 draft only; Session 46
// keeps it as the app-level `ui.richFeaturePopup` (PR-1).
assert.deepEqual(CURRENT_WRITER_ADV_FIELDS, [
  ...Object.keys(createDefaultAdv()),
  'plot_title_position',
  'losatProgram',
  'rich_feature_popup'
]);
assert.equal(Object.hasOwn(createDefaultAdv(), 'rich_feature_popup'), false);
assert.equal(createDefaultLosat().blastp.candidateLimit, null);
// How LOSAT runs is one app-level setting (`ui.losatExecution`).
assert.equal(Object.hasOwn(createDefaultLosat(), 'executionMode'), false);
assert.equal(createDefaultLosatExecution().executionMode, 'threaded');
assert.equal(createDefaultLosat().blastp.collinearSearchScope, 'adjacent');
assert.equal(createDefaultLosat().blastp.collinearMergeOrientation, 'either');
assert.equal(createDefaultAdv('circular').pairwise_match_style, 'ribbon');
assert.equal(createDefaultAdv('linear').pairwise_match_style, 'curve');
assert.equal(normalizeCurrentPairwiseMatchStyle('curve', 'ribbon'), 'curve');
assert.equal(normalizeCurrentPairwiseMatchStyle(undefined, 'ribbon'), 'ribbon');
assert.equal(normalizeCurrentPairwiseMatchStyle('invalid', 'curve'), 'curve');
assert.throws(
  () => normalizeCurrentPairwiseMatchStyle(undefined, 'invalid'),
  /fallback must be one of/
);

assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
  mode: 'circular',
  storedConfig: {
    ...storedConfig,
    unmanagedConfigOverrides: {
      'objects.gc_content.percent_background_opacity': 0.42
    }
  }
}));
assert.throws(
  () => validateCurrentWriterActiveConfig({
    mode: 'circular',
    storedConfig: { ...storedConfig, unmanagedConfigOverrides: [] }
  }),
  /config\.unmanagedConfigOverrides must be object/
);
const unsafeStoredConfig = JSON.parse(JSON.stringify({
  ...storedConfig,
  unmanagedConfigOverrides: {
    'labels.filtering.raw': JSON.parse('{"__proto__":{"polluted":true}}')
  }
}));
assert.throws(
  () => validateCurrentWriterActiveConfig({
    mode: 'circular',
    storedConfig: unsafeStoredConfig
  }),
  /unsafe key __proto__/
);

for (const [candidateLimit, collinearSearchScope] of [[null, 'adjacent'], [9, 'all']]) {
  assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
    mode: 'linear',
    storedConfig: {
      ...storedConfig,
      losat: {
        ...createDefaultLosat(),
        blastp: {
          ...createDefaultLosat().blastp,
          mode: 'collinear',
          candidateLimit,
          collinearSearchScope
        }
      }
    }
  }));
}
for (const blastp of [
  { mode: 'unsupported' },
  { candidateLimit: 0 },
  { maxHits: 0 },
  { orthogroupMembershipMode: 'legacy' },
  { orthogroupMemberMaxHits: 0 },
  { collinearMinAnchors: 0 },
  { collinearMaxUnitGap: -1 },
  { collinearMaxDiagonalDrift: -1 },
  { collinearMaxConflictsInMergeGap: -1 },
  { collinearMaxParalogLinksPerOrthogroup: 0 },
  { collinearUnitMode: 'gene' },
  { collinearAnchorMode: 'top1' },
  { collinearMergeOrientation: 'both' },
  { collinearColorMode: 'score' },
  { collinearSearchScope: 'global' }
]) {
  assert.throws(() => validateCurrentWriterActiveConfig({
    mode: 'linear',
    storedConfig: {
      ...storedConfig,
      losat: { blastp }
    }
  }));
}

const retiredTrackOrder = structuredClone(storedConfig);
retiredTrackOrder.adv.cli_circular_track_order = ['features'];
assert.throws(
  () => validateCurrentWriterActiveConfig({
    mode: 'circular',
    storedConfig: retiredTrackOrder
  }),
  /config\.adv.*cli_circular_track_order/
);

const retiredTrackSlots = structuredClone(storedConfig);
retiredTrackSlots.adv.cli_circular_track_slots = [];
assert.throws(
  () => validateCurrentWriterActiveConfig({
    mode: 'circular',
    storedConfig: retiredTrackSlots
  }),
  /config\.adv.*cli_circular_track_slots/
);

assert.equal(createDefaultForm().keep_definition_left_aligned, true);
assert.equal(createDefaultAdv('linear').linear_show_replicon, false);
for (const locked of [false, true]) {
  assert.doesNotThrow(() => validateCurrentWriterActiveConfig({
    mode: 'linear', storedConfig: {
      ...storedConfig, form: { ...storedConfig.form, keep_definition_left_aligned: locked }
    }
  }));
}
for (const malformed of [null, 'false', 'true', 0, 1, [], {}, undefined]) {
  assert.throws(() => validateCurrentWriterActiveConfig({
    mode: 'linear', storedConfig: {
      ...storedConfig, form: { ...storedConfig.form, keep_definition_left_aligned: malformed }
    }
  }), /keep_definition_left_aligned must be a boolean/);
}

for (const field of ['width', 'radius']) {
  for (const number of [NaN, Infinity, -Infinity]) {
    for (const value of [number, { value: number, unit: 'px' }]) {
      const invalid = structuredClone(storedConfig);
      invalid.adv.circular_track_slots = [{
        id: 'features', renderer: 'features', enabled: true, side: 'outside', z: 0,
        params: { lane_direction: 'outside' }, width: null, radius: null,
        inner_gap_px: null, outer_gap_px: null,
        [field]: value
      }];
      assert.throws(() => validateCurrentWriterActiveConfig({
        mode: 'circular', storedConfig: invalid
      }), /must be a number greater than 0, in px or ×R/);
    }
  }
}

// OV-38 (R6): an obsolete Circular track slot field is a classified failure that
// names the field and the row, whatever its value and whether the slots are on;
// readers drop only the lossless legacy null (config.js).
const obsoleteSlotFailure = (adv) => {
  try {
    validateImportedCircularTrackSlots({ adv });
  } catch (error) {
    return normalizeUserFacingError(error);
  }
  return null;
};
const featureSlot = {
  id: 'features', renderer: 'features', enabled: true, side: 'inside', z: 0,
  params: { lane_direction: 'inside' }, width: null, radius: null, inner_gap_px: null, outer_gap_px: null
};
const ticksSlot = {
  id: 'ticks', renderer: 'ticks', enabled: true, side: 'inside', z: 0,
  params: { tick_label_layout: 'label_in_tick_out' }, width: null, radius: null, inner_gap_px: null, outer_gap_px: null
};
for (const [enabled, slots, field, slotIndex] of [
  [false, [featureSlot, { ...ticksSlot, spacing: null }], 'spacing', 1],
  [true, [{ ...featureSlot, spacing: null }, ticksSlot], 'spacing', 0],
  [false, [{ ...featureSlot, spacing: '4px' }, ticksSlot], 'spacing', 0],
  [false, [featureSlot, { ...ticksSlot, inner_radius: 0.5 }], 'inner_radius', 1],
  [false, [featureSlot, { ...ticksSlot, params: { ...ticksSlot.params, spacing: null } }], 'spacing', 1]
]) {
  const failure = obsoleteSlotFailure({
    ...createDefaultAdv('circular'), circular_track_slots_enabled: enabled, circular_track_slots: slots
  });
  assert.equal(failure?.code, 'TRACK_INVALID', field);
  assert.deepEqual(failure.context, { field, reason: 'OBSOLETE_TRACK_FIELD', slotIndex });
  assert.equal(failure.summary, `The track settings are invalid. Track row ${slotIndex + 1}. Field: ${field}. `
    + 'Custom Track Slots no longer read this field. Use slot-level radius, width, inner_gap_px, '
    + 'outer_gap_px, side, and z fields.');
}
assert.equal(obsoleteSlotFailure({
  ...createDefaultAdv('circular'), circular_track_slots: [featureSlot, ticksSlot]
}), null);

// OV-269: shared with tests/test_draft_value_migrations.py.
// A CLI-written Session 40-41 config holds no Web draft.
const { readFile } = await import('node:fs/promises');
const cliWriterConfigVectors = JSON.parse(await readFile(
  new URL('../fixtures/cli-writer-config-vectors.json', import.meta.url), 'utf8'
)).cases;
for (const { name, session, holdsCliWriterConfig: expected } of cliWriterConfigVectors) {
  assert.equal(holdsCliWriterConfig(session), expected, name);
}
