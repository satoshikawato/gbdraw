import assert from 'node:assert/strict';

import {
  describeLinearLabelVisibility,
  migrateLegacyLinearLabelVisibility,
  requireLinearLabelVisibilityMode,
  resolveLinearLabelVisibility
} from '../../gbdraw/web/js/app/linear-label-visibility.js';

for (const recordRows of [[1], [1, 2, 3, 4]]) {
  const hasSharedRow = new Set(recordRows).size !== recordRows.length;
  assert.equal(resolveLinearLabelVisibility('auto', { hasSharedRow }), true);
}
assert.equal(resolveLinearLabelVisibility('auto', { hasSharedRow: true }), false);
assert.equal(resolveLinearLabelVisibility('show', { hasSharedRow: true }), true);
assert.equal(resolveLinearLabelVisibility('hide', { hasSharedRow: false }), false);
assert.deepEqual([
  resolveLinearLabelVisibility('show', { hasSharedRow: true }),
  resolveLinearLabelVisibility('auto', { hasSharedRow: true })
], [true, false], 'Accession=Show and Length=Auto resolve independently');
assert.deepEqual([
  resolveLinearLabelVisibility('hide', { hasSharedRow: false }),
  resolveLinearLabelVisibility('show', { hasSharedRow: false })
], [false, true], 'Accession=Hide and Length=Show resolve independently');
const selectedModes = ['auto', 'show'];
assert.deepEqual(
  [true, false, true].map((hasSharedRow) => selectedModes.map((mode) => (
    resolveLinearLabelVisibility(mode, { hasSharedRow })
  ))),
  [[false, true], [true, true], [false, true]],
  'shared-row changes only the effective Auto result, not the selected modes'
);
assert.deepEqual(selectedModes, ['auto', 'show']);
assert.equal(describeLinearLabelVisibility('auto'), 'Auto · Shown');
assert.equal(describeLinearLabelVisibility('auto', { hasSharedRow: true }), 'Auto · Hidden');
assert.throws(() => requireLinearLabelVisibilityMode('sometimes'), /must be one of/);

assert.deepEqual(migrateLegacyLinearLabelVisibility({
  linear_show_accession: true,
  linear_show_length: false
}), {
  linear_accession_visibility: 'show',
  linear_length_visibility: 'hide'
});
assert.deepEqual(migrateLegacyLinearLabelVisibility({}), {
  linear_accession_visibility: 'show',
  linear_length_visibility: 'show'
});
assert.deepEqual(migrateLegacyLinearLabelVisibility({
  linear_accession_visibility: 'auto',
  linear_length_visibility: 'hide',
  linear_show_accession: false,
  linear_show_length: true
}), {
  linear_accession_visibility: 'auto',
  linear_length_visibility: 'hide'
});
assert.throws(
  () => migrateLegacyLinearLabelVisibility({ linear_show_accession: 'yes' }),
  /legacy value must be a boolean/
);
