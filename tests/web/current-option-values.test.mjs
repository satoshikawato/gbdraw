import assert from 'node:assert/strict';

import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames,
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentWebStateFieldNames
} from '../../gbdraw/web/js/app/current-option-values.js';

assert.equal(requireCurrentCircularMultiRecordSizeMode(), 'auto');
assert.equal(requireCurrentCircularMultiRecordSizeMode('linear'), 'linear');
assert.throws(
  () => requireCurrentCircularMultiRecordSizeMode('sqrt'),
  /Circular multi-record size mode/
);

assert.equal(requireCurrentLinearTrackLayout(), 'middle');
assert.equal(requireCurrentLinearTrackLayout('below'), 'below');
assert.throws(
  () => requireCurrentLinearTrackLayout('spreadout'),
  /Linear track layout/
);
assert.throws(
  () => requireCurrentLinearTrackLayout('tuckin'),
  /Linear track layout/
);

assert.equal(requireCurrentLinearLabelPlacement(), 'auto');
assert.equal(requireCurrentLinearLabelPlacement('above_feature'), 'above_feature');
assert.throws(
  () => requireCurrentLinearLabelPlacement('on_feature'),
  /Linear label placement/
);

assert.equal(migratePersistedCircularMultiRecordSizeMode('sqrt'), 'auto');
assert.equal(migratePersistedLinearTrackLayout('spreadout'), 'above');
assert.equal(migratePersistedLinearTrackLayout('tuckin'), 'below');
assert.equal(migratePersistedLinearLabelPlacement('on_feature'), 'above_feature');

const persisted = {
  adv: {
    depth_tick_interval: 25,
    depth_tracks: [
      { tick_interval: 10 },
      { tick_interval: 20, large_tick_interval: 30 }
    ]
  },
  losat: {
    blastp: {
      collinearMaxGeneGap: 4
    }
  }
};
const migrated = migratePersistedWebStateFieldNames(persisted);
assert.equal(migrated.adv.depth_large_tick_interval, 25);
assert.equal(migrated.adv.depth_tracks[0].large_tick_interval, 10);
assert.equal(migrated.adv.depth_tracks[1].large_tick_interval, 30);
assert.equal(migrated.losat.blastp.collinearMaxUnitGap, 4);
assert.equal(Object.hasOwn(migrated.adv, 'depth_tick_interval'), false);
assert.equal(Object.hasOwn(migrated.adv.depth_tracks[0], 'tick_interval'), false);
assert.equal(Object.hasOwn(migrated.losat.blastp, 'collinearMaxGeneGap'), false);
assert.equal(persisted.adv.depth_tick_interval, 25, 'migration must not mutate persisted data');
assert.equal(persisted.adv.depth_tracks[0].tick_interval, 10);
assert.equal(persisted.losat.blastp.collinearMaxGeneGap, 4);
assert.equal(requireCurrentWebStateFieldNames(migrated), migrated);

assert.throws(
  () => requireCurrentWebStateFieldNames({ adv: { depth_tick_interval: 10 } }),
  /depth_tick_interval is obsolete/
);
assert.throws(
  () => requireCurrentWebStateFieldNames({
    adv: { depth_tracks: [{ tick_interval: 10 }] }
  }),
  /tick_interval is obsolete/
);
assert.throws(
  () => requireCurrentWebStateFieldNames({
    losat: { blastp: { collinearMaxGeneGap: 2 } }
  }),
  /collinearMaxGeneGap is obsolete/
);
