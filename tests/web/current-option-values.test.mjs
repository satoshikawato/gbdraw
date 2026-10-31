import assert from 'node:assert/strict';

import {
  migratePersistedCircularMultiRecordSizeMode,
  migratePersistedLinearLabelPlacement,
  migratePersistedLinearTrackLayout,
  migratePersistedWebStateFieldNames,
  requireCurrentCircularMultiRecordSizeMode,
  requireCurrentCollinearAnchorMode,
  requireCurrentCollinearColorMode,
  requireCurrentCollinearInferOrthogroups,
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
} from '../../gbdraw/web/js/services/current-option-values.js';
import {
  ABSENT_COLLINEAR_INFER_ORTHOGROUPS,
  readCollinearAnchorMode,
  readCollinearColorMode,
  readCollinearSearchScope,
  readCollinearUnitMode,
  readOrthogroupMembershipMode
} from '../../gbdraw/web/js/services/losat-normalization.js';
import { CANONICAL_REQUEST_DESCRIPTION } from '../../gbdraw/web/js/contracts/request-shape.generated.js';
import { readSpelling } from '../../gbdraw/web/js/utils/contract-shape.js';

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
assert.throws(() => requireCurrentCollinearSearchScope('all-records'), /COLLINEAR_SEARCH_SCOPE/);

for (const [validator, validValues, invalidValue] of [
  [requireCurrentCollinearUnitMode, ['auto', 'cds', 'locus'], 'gene'],
  [requireCurrentCollinearAnchorMode, ['all', 'rbh'], 'one_to_one'],
  [requireCurrentCollinearMergeOrientation, ['strand', 'order'], 'both'],
  [requireCurrentCollinearColorMode,
    ['average_identity', 'orientation', 'orientation_identity'], 'score'],
  [requireCurrentOrthogroupMembershipMode, ['anchor_core_v1'], 'banana']
]) {
  for (const value of validValues) assert.equal(validator(value), value);
  assert.throws(() => validator(invalidValue));
}
assert.equal(requireCurrentCollinearColorMode('identity'), 'average_identity');
// "either" is the old spelling of "strand" (OV-2301).
assert.equal(requireCurrentCollinearMergeOrientation('either'), 'strand');

for (const [validator, defaultValue, minimum] of [
  [requireCurrentProteinBlastpMaxHits, 5, 1],
  [requireCurrentOrthogroupMemberMaxHits, null, 1],
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
  assert.equal(requireCurrentOrthogroupMemberMaxHits(omitted), null);
  assert.equal(requireCurrentProteinBlastpCandidateLimit(omitted), null);
}
assert.equal(requireCurrentProteinBlastpCandidateLimit(7), 7);
assert.equal(requireCurrentProteinBlastpCandidateLimit('7'), 7);
for (const invalid of [0, -1, 1.5, 'unbounded']) {
  assert.throws(
    () => requireCurrentProteinBlastpCandidateLimit(invalid),
    /Max target seqs/
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

// A stored draft without the inference flag inferred its Collinear groups: the
// absent reading of the generated description, not the fresh default (off).
assert.equal(ABSENT_COLLINEAR_INFER_ORTHOGROUPS, true);
assert.equal(requireCurrentCollinearInferOrthogroups(undefined), ABSENT_COLLINEAR_INFER_ORTHOGROUPS);
assert.equal(requireCurrentCollinearInferOrthogroups(false), false);

// A comparison setting refuses with a diagnostic, not a message that repeats
// the setting's value list for the error table to match.
for (const [require, field] of /** @type {const} */ ([
  [requireCurrentCollinearSearchScope, 'collinear_search_scope'],
  [requireCurrentCollinearUnitMode, 'collinear_unit_mode'],
  [requireCurrentCollinearAnchorMode, 'collinear_anchor_mode'],
  [requireCurrentCollinearColorMode, 'collinear_color_mode'],
  [requireCurrentOrthogroupMembershipMode, 'orthogroup_membership_mode']
])) {
  assert.throws(() => require('banana'), (error) => {
    assert.equal(error.code, 'INPUT_INVALID', field);
    assert.deepEqual(error.context, { field, reason: field.toUpperCase() });
    return true;
  });
}

// A draft that holds the removed anchor mode one_to_one, or a spelling of it,
// loads it as rbh with a Load row (OV-2303).
for (const retired of ['one_to_one', 'top1', 'mutual-best']) {
  assert.deepEqual(readCollinearAnchorMode(retired), { value: 'rbh', lost: true }, retired);
}

// The five comparison settings: the Web reads each spelling exactly as Python's
// read_spelling does, from the generated description (OV-776). Load gives the
// name, or the setting's default as a lost value; Generate gives the name or
// refuses. A spelling only the Web once read (`all_pairs`, `global`) is lost.
const comparisonSettings = CANONICAL_REQUEST_DESCRIPTION.types.CanonicalProteinComparisonSettings.fields;
for (const [name, read, require] of /** @type {const} */ ([
  ['collinearityUnitMode', readCollinearUnitMode, requireCurrentCollinearUnitMode],
  ['collinearityAnchorMode', readCollinearAnchorMode, requireCurrentCollinearAnchorMode],
  ['collinearitySearchScope', readCollinearSearchScope, requireCurrentCollinearSearchScope],
  ['collinearityColorMode', readCollinearColorMode, requireCurrentCollinearColorMode],
  ['orthogroupMembershipMode', readOrthogroupMembershipMode, requireCurrentOrthogroupMembershipMode]
])) {
  const shape = comparisonSettings[name];
  const aliases = Object.keys(shape.spellings.aliases ?? {});
  const spellings = [
    ...shape.values, ...aliases, ...[...shape.values, ...aliases].flatMap((value) => [
      ` ${value.toUpperCase()} `, value.replaceAll('_', '-')
    ]),
    '', null, undefined, 'banana', 5, 'all_pairs', 'adjacent_pairs', 'global', 'all-records'
  ];
  for (const spelling of spellings) {
    const expected = readSpelling(spelling, shape);
    assert.deepEqual(read(spelling), expected === null
      ? { value: shape.spellings.blank, lost: true } : { value: expected, lost: false }, `${name} reads ${String(spelling)}`);
    if (expected === null) assert.throws(() => require(spelling), Error, `${name} refuses ${String(spelling)}`);
    else assert.equal(require(spelling), expected, `${name} requires ${String(spelling)}`);
  }
}
assert.deepEqual(readCollinearSearchScope('all_pairs'), { value: 'adjacent', lost: true });
assert.deepEqual(readCollinearSearchScope('global'), { value: 'adjacent', lost: true });
