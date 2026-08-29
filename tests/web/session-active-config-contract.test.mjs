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
