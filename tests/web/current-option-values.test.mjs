import assert from 'node:assert/strict';

import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames,
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentCollinearAnchorMode,
  requireCurrentCollinearColorMode,
  requireCurrentCollinearMaxConflicts,
  requireCurrentCollinearMaxDiagonalDrift,
  requireCurrentCollinearMaxParalogLinks,
  requireCurrentCollinearMaxUnitGap,
  requireCurrentCollinearMergeOrientation,
  requireCurrentCollinearMinAnchors,
  requireCurrentCollinearSearchScope,
  requireCurrentCollinearUnitMode,
  requireCurrentLinearLabelPlacement,
  requireCurrentLinearTrackLayout,
  requireCurrentOrthogroupMemberMaxHits,
  requireCurrentOrthogroupMembershipMode,
  requireCurrentProteinBlastpCandidateLimit,
  requireCurrentProteinBlastpMaxHits,
  requireCurrentProteinBlastpMode,
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

assert.equal(requireCurrentProteinBlastpMode(), 'orthogroup');
for (const mode of ['pairwise', 'orthogroup', 'collinear']) {
  assert.equal(requireCurrentProteinBlastpMode(mode), mode);
}
assert.throws(() => requireCurrentProteinBlastpMode('similarity'), /Protein BLASTP mode/);

assert.equal(requireCurrentCollinearSearchScope(), 'adjacent');
assert.equal(requireCurrentCollinearSearchScope('adjacent'), 'adjacent');
assert.equal(requireCurrentCollinearSearchScope('all'), 'all');
assert.throws(() => requireCurrentCollinearSearchScope('all-records'), /search scope/i);

for (const [validator, validValues, invalidValue] of [
  [requireCurrentCollinearUnitMode, ['auto', 'cds', 'locus'], 'gene'],
  [requireCurrentCollinearAnchorMode, ['all', 'one_to_one', 'rbh'], 'top1'],
  [requireCurrentCollinearMergeOrientation, ['strand', 'order', 'either'], 'both'],
  [requireCurrentCollinearColorMode,
    ['average_identity', 'orientation', 'orientation_identity'], 'score'],
  [requireCurrentOrthogroupMembershipMode, ['anchor_core_v1'], 'legacy']
]) {
  for (const value of validValues) assert.equal(validator(value), value);
  assert.throws(() => validator(invalidValue));
}
assert.equal(requireCurrentCollinearColorMode('identity'), 'average_identity');

for (const [validator, defaultValue, minimum] of [
  [requireCurrentProteinBlastpMaxHits, 5, 1],
  [requireCurrentOrthogroupMemberMaxHits, 5, 1],
  [requireCurrentCollinearMinAnchors, 1, 1],
  [requireCurrentCollinearMaxUnitGap, 0, 0],
  [requireCurrentCollinearMaxDiagonalDrift, 0, 0],
  [requireCurrentCollinearMaxConflicts, 1, 0],
  [requireCurrentCollinearMaxParalogLinks, 2, 1]
]) {
  assert.equal(validator(), defaultValue);
  assert.equal(validator(minimum), minimum);
  assert.throws(() => validator(minimum - 1));
  assert.throws(() => validator(1.5));
}

for (const omitted of [undefined, null, '']) {
  assert.equal(requireCurrentProteinBlastpCandidateLimit(omitted), null);
}
assert.equal(requireCurrentProteinBlastpCandidateLimit(7), 7);
assert.equal(requireCurrentProteinBlastpCandidateLimit('7'), 7);
for (const invalid of [0, -1, 1.5, 'unbounded']) {
  assert.throws(
    () => requireCurrentProteinBlastpCandidateLimit(invalid),
    /Candidate limit/
  );
}

assert.equal(migratePersistedCircularMultiRecordSizeMode('sqrt'), 'auto');
assert.equal(migratePersistedLinearTrackLayout('spreadout'), 'above');
assert.equal(migratePersistedLinearTrackLayout('tuckin'), 'below');
assert.equal(migratePersistedLinearLabelPlacement('on_feature'), 'above_feature');

const persisted = {
  adv: {
    arrow_shaft_width_ratio: 0.75,
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
assert.equal(migrated.adv.arrow_shaft_width_ratio, 0.75);
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
