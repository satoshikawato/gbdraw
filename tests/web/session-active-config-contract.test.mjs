import assert from 'node:assert/strict';

assert.equal(globalThis.window, undefined);
assert.equal(globalThis.document, undefined);

const {
  CURRENT_WRITER_ADV_FIELDS,
  CURRENT_WRITER_FORM_FIELDS,
  createDefaultAdv,
  createDefaultForm,
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
