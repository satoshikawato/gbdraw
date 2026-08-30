import assert from 'node:assert/strict';

assert.equal(globalThis.window, undefined);
assert.equal(globalThis.document, undefined);

const {
  CURRENT_WRITER_ADV_FIELDS,
  CURRENT_WRITER_FORM_FIELDS,
  createDefaultAdv,
  createDefaultForm,
  createDefaultLosat,
  validateCurrentWriterActiveConfig
} = await import('../../gbdraw/web/js/services/session-active-config-contract.js');

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
assert.deepEqual(CURRENT_WRITER_ADV_FIELDS, [
  ...Object.keys(createDefaultAdv()),
  'plot_title_position',
  'losatProgram'
]);
assert.equal(createDefaultLosat().blastp.candidateLimit, null);
assert.equal(createDefaultLosat().blastp.collinearSearchScope, 'adjacent');
assert.equal(createDefaultLosat().blastp.collinearMergeOrientation, 'either');

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
