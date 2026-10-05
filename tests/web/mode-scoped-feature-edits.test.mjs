// R2 (residual R-1 of the override-precedence audit): Circular and Linear can
// use the same record key (`record-1` in a Gallery or Python Session) for the
// same feature. Each draft row names its mode (`scope`), so a request, the
// live projection, and every reconcile reach only the rows of their own mode,
// and a mode change keeps the other mode's rows in the draft.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync, statSync } from 'node:fs';
import { join, relative } from 'node:path';
import test from 'node:test';
import {
  countUnresolvedFeatureEdits,
  pruneUnmatchedFeatureOverrides,
  removeUnresolvedFeatureEdits
} from '../../gbdraw/web/js/app/feature-visibility.js';
import { replaceFeatureEdits } from '../../gbdraw/web/js/app/feature-editor/feature-edit-table.js';
import { migrateSessionFeaturePlacements } from '../../gbdraw/web/js/services/feature-edit-migration.js';
import {
  canonicalFeatureOverrides,
  canonicalFeaturePlacements,
  featureIdentityKeyOf,
  requestFeatureOverrides,
  requestFeaturePlacements,
  updateFeatureOverride
} from '../../gbdraw/web/js/services/feature-placement.js';

const MODES = ['circular', 'linear'];
const SIDES = { circular: ['outward', 'inward'], linear: ['above', 'below'] };
const other = (mode) => (mode === 'circular' ? 'linear' : 'circular');
// A catalog feature: an admitted catalog gives each feature its Result's mode.
const feature = (scope, recordKey, biologicalFeatureId) => ({
  scope, record_key: recordKey, biological_feature_id: biologicalFeatureId
});
// The records of both modes share keys; `all-input` is an ALL record.
const records = [{ recordKey: 'record-1' }, { recordKey: 'all-input', cardinality: 'all' }];
const recordKeys = ['record-1', 'all-input:2'];

// Both draft maps hold, for each mode, every edit kind on the same identities.
const drafts = () => {
  const featureOverrides = {};
  const featurePlacementOverrides = {};
  const bulkLabelText = {};
  MODES.forEach((mode) => recordKeys.forEach((recordKey) => ['f1', 'f2', 'f3'].forEach((id, index) => {
    const item = feature(mode, recordKey, id);
    updateFeatureOverride(featureOverrides, item, [
      { featureVisibility: 'off' }, { labelVisibility: 'on', labelText: `${mode} text` }, { labelSourceText: 'source' }
    ][index]);
    const row = { scope: mode, recordKey, biologicalFeatureId: id, placement: index === 0
      ? { kind: 'main' } : { kind: 'lane', side: SIDES[mode][index - 1], level: 1 } };
    featurePlacementOverrides[featureIdentityKeyOf(row)] = row;
    bulkLabelText[featureIdentityKeyOf(item)] = `${mode} bulk`;
  })));
  return { featureOverrides, featurePlacementOverrides, bulkLabelText };
};

test('a draft row reaches only the requests of its own mode (both draft maps)', () => {
  const { featureOverrides, featurePlacementOverrides, bulkLabelText } = drafts();
  // The same identity in the two modes is two draft rows.
  assert.equal(Object.keys(featureOverrides).length, 2 * recordKeys.length * 3);
  assert.equal(Object.keys(featurePlacementOverrides).length, 2 * recordKeys.length * 3);
  MODES.forEach((mode) => {
    const placements = requestFeaturePlacements(featurePlacementOverrides, mode, records);
    assert.equal(placements.length, recordKeys.length * 3);
    placements.forEach((row) => assert.ok(row.placement.kind === 'main' || SIDES[mode].includes(row.placement.side)));
    // Python accepts the rows as the request's (no draft-only field).
    assert.deepEqual(canonicalFeaturePlacements(placements, mode), placements);
    const overrides = requestFeatureOverrides(featureOverrides, mode, records, { bulkLabelText });
    assert.equal(overrides.length, recordKeys.length * 3);
    overrides.forEach((row) => {
      assert.ok(!Object.hasOwn(row, 'scope'));
      assert.ok(row.labelText === null || row.labelText.startsWith(mode), JSON.stringify(row));
    });
    assert.deepEqual(canonicalFeatureOverrides(overrides), overrides);
  });
});

test('the draft key encodes the row mode, and a row without one is invalid', () => {
  const row = { scope: 'circular', recordKey: 'record-1', biologicalFeatureId: 'f1', placement: { kind: 'main' } };
  assert.equal(featureIdentityKeyOf(row), JSON.stringify(['circular', 'record-1', 'f1']));
  assert.equal(featureIdentityKeyOf({ recordKey: 'record-1', biologicalFeatureId: 'f1' }), '');
  assert.throws(() => canonicalFeaturePlacements({ [JSON.stringify(['linear', 'record-1', 'f1'])]: row }));
  assert.throws(() => canonicalFeaturePlacements({ [JSON.stringify(['record-1', 'f1'])]: row }));
  // A lane names a side of its own mode.
  assert.throws(() => canonicalFeaturePlacements({ [featureIdentityKeyOf(row)]: {
    ...row, placement: { kind: 'lane', side: 'above', level: 1 } } }));
  const { scope: _scope, ...modeless } = row;
  assert.throws(() => canonicalFeaturePlacements({ [JSON.stringify(['record-1', 'f1'])]: modeless }));
});

test('notices, the source-replacing reconcile, and Remove unmatched touch only their mode', () => {
  MODES.forEach((mode) => {
    const notices = recordKeys.map((recordKey) => ({
      recordKey, biologicalFeatureId: 'f1', status: 'unresolved', kinds: ['placement', 'feature_visibility'], resultIndex: 0
    }));
    const { featureOverrides, featurePlacementOverrides } = drafts();
    assert.equal(countUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices, scope: mode }),
      2 * recordKeys.length);
    assert.equal(removeUnresolvedFeatureEdits({ featureOverrides, featurePlacementOverrides, notices, scope: mode }),
      2 * recordKeys.length);
    const kept = drafts();
    Object.entries(kept.featurePlacementOverrides).forEach(([key, row]) => {
      if (row.scope === other(mode)) assert.deepEqual(featurePlacementOverrides[key], row);
    });
    Object.entries(kept.featureOverrides).forEach(([key, row]) => {
      if (row.scope === other(mode)) assert.deepEqual(featureOverrides[key], row);
    });

    const pruned = drafts();
    const removed = pruneUnmatchedFeatureOverrides({
      ...pruned, notices, scope: mode, replacedRecordKeys: recordKeys,
      previousRecords: [...records, { recordKey: 'dropped' }], currentRecords: records, biologicalFeatures: []
    });
    // The unresolved placements and visibility edits, and the source texts of
    // features the replaced source lost, of this mode only.
    assert.equal(removed, 2 * recordKeys.length);
    Object.entries(kept.featureOverrides).forEach(([key, row]) => {
      if (row.scope === other(mode)) assert.deepEqual(pruned.featureOverrides[key], row);
      else assert.equal(pruned.featureOverrides[key]?.labelSourceText ?? null, null);
    });
  });
});

test('Load Feature Edits TSV replaces only the edits of the committed mode', () => {
  const { featureOverrides } = drafts();
  const before = structuredClone(featureOverrides);
  replaceFeatureEdits(featureOverrides, [{
    recordKey: 'record-1', biologicalFeatureId: 'f9', featureVisibility: 'on', labelVisibility: null, labelText: null
  }], 'linear', records);
  Object.entries(before).forEach(([key, row]) => {
    if (row.scope === 'circular') assert.deepEqual(featureOverrides[key], row);
  });
  assert.deepEqual(requestFeatureOverrides(featureOverrides, 'linear', records), [{
    recordKey: 'record-1', biologicalFeatureId: 'f9', featureVisibility: 'on', labelVisibility: null, labelText: null
  }]);
});

// A Session 44 draft placement reached every request with its record key (a
// lane failed the other mode's Generate). Session 45 keeps that reach.
test('Session 44 placement drafts take a mode: a lane its side, a Main row both', () => {
  const main = { recordKey: 'record-1', biologicalFeatureId: 'f1', placement: { kind: 'main' } };
  const lane = { recordKey: 'record-1', biologicalFeatureId: 'f2', placement: { kind: 'lane', side: 'below', level: 1 } };
  const migrated = migrateSessionFeaturePlacements({
    [JSON.stringify(['record-1', 'f1'])]: main,
    [JSON.stringify(['record-1', 'f2'])]: lane
  });
  assert.deepEqual(migrated, {
    [JSON.stringify(['circular', 'record-1', 'f1'])]: { scope: 'circular', ...main },
    [JSON.stringify(['linear', 'record-1', 'f1'])]: { scope: 'linear', ...main },
    [JSON.stringify(['linear', 'record-1', 'f2'])]: { scope: 'linear', ...lane }
  });
  canonicalFeaturePlacements(migrated);
});

// The next path that indexes a draft map by a hand-built [recordKey, feature]
// key would read or write the other mode's row: draft keys come only from
// services/feature-placement.js. The Session 44 migration reads the old key,
// and the exported SVG's runtime indexes the catalog of its one Result.
test('only the draft owner builds a draft identity key', () => {
  const root = new URL('../../gbdraw/web/js/', import.meta.url).pathname;
  const files = [];
  const walk = (dir) => readdirSync(dir).forEach((name) => {
    const path = join(dir, name);
    if (statSync(path).isDirectory()) walk(path);
    else if (name.endsWith('.js')) files.push(path);
  });
  walk(root);
  const owners = new Set([
    'services/feature-placement.js',
    'services/feature-edit-migration.js',
    'services/standalone-interactivity-assets.js'
  ]);
  const offenders = files.map((path) => relative(root, path)).filter((path) => !owners.has(path))
    .flatMap((path) => {
      const source = readFileSync(join(root, path), 'utf8').replace(/^\s*\/\/.*$/gm, '');
      return [...source.matchAll(/JSON\.stringify\(\[([^[\]]*)\]/g)]
        .filter(([, items]) => /record_?[kK]ey/.test(items) && /biological_?[fF]eature_?[iI]d/.test(items))
        .map((match) => `${path}:${source.slice(0, match.index).split('\n').length}`);
    });
  assert.deepEqual(offenders, []);
});
