// R2 (residual R-1 of the override-precedence audit): Circular and Linear can
// use the same record key (`record-1` in a Gallery or Python Session) for the
// same feature. PR-1: each mode has its own drawing, so the same identity has
// one draft row in each drawing, keyed by the identity pair without a mode; a
// request, the live projection, and every reconcile read the drawing of their
// own mode, and the other drawing keeps its rows.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync, statSync } from 'node:fs';
import { join, relative } from 'node:path';
import test from 'node:test';
import {
  countUnresolvedFeatureEdits,
  pruneUnmatchedFeatureOverrides,
  removeUnresolvedFeatureEdits
} from '../../gbdraw/web/js/services/feature-visibility.js';
import { replaceFeatureEdits } from '../../gbdraw/web/js/app/feature-editor/feature-edit-table.js';
import { migrateSessionFeaturePlacements } from '../../gbdraw/web/js/services/feature-edit-migration.js';
import { unscopedDraftRows } from '../../gbdraw/web/js/services/mode-scoped-migration.js';
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
// A catalog feature: an admitted catalog gives each feature its Result's mode.
const feature = (scope, recordKey, biologicalFeatureId) => ({
  scope, record_key: recordKey, biological_feature_id: biologicalFeatureId
});
// The records of both modes share keys; `all-input` is an ALL record.
const records = [{ recordKey: 'record-1' }, { recordKey: 'all-input', cardinality: 'all' }];
const recordKeys = ['record-1', 'all-input:2'];

// Each mode's drawing holds every edit kind on the same identities.
const drawingDrafts = (mode) => {
  const featureOverrides = {};
  const featurePlacementOverrides = {};
  const bulkLabelText = {};
  recordKeys.forEach((recordKey) => ['f1', 'f2', 'f3'].forEach((id, index) => {
    const item = feature(mode, recordKey, id);
    updateFeatureOverride(featureOverrides, item, [
      { featureVisibility: 'off' }, { labelVisibility: 'on', labelText: `${mode} text` }, { labelSourceText: 'source' }
    ][index]);
    const row = { recordKey, biologicalFeatureId: id, placement: index === 0
      ? { kind: 'main' } : { kind: 'lane', side: SIDES[mode][index - 1], level: 1 } };
    featurePlacementOverrides[featureIdentityKeyOf(row)] = row;
    bulkLabelText[featureIdentityKeyOf(item)] = `${mode} bulk`;
  }));
  return { featureOverrides, featurePlacementOverrides, bulkLabelText };
};
const drawings = () => ({ circular: drawingDrafts('circular'), linear: drawingDrafts('linear') });

test('each drawing reaches only the requests of its own mode (both draft maps)', () => {
  const drafts = drawings();
  MODES.forEach((mode) => {
    const { featureOverrides, featurePlacementOverrides, bulkLabelText } = drafts[mode];
    // One row per identity in each drawing.
    assert.equal(Object.keys(featureOverrides).length, recordKeys.length * 3);
    assert.equal(Object.keys(featurePlacementOverrides).length, recordKeys.length * 3);
    const placements = requestFeaturePlacements(featurePlacementOverrides, mode, records);
    assert.equal(placements.length, recordKeys.length * 3);
    placements.forEach((row) => assert.ok(row.placement.kind === 'main' || SIDES[mode].includes(row.placement.side)));
    // Python accepts the rows as the request's (no draft-only field).
    assert.deepEqual(canonicalFeaturePlacements(placements, mode), placements);
    const overrides = requestFeatureOverrides(featureOverrides, records, { bulkLabelText });
    assert.equal(overrides.length, recordKeys.length * 3);
    overrides.forEach((row) => {
      assert.ok(!Object.hasOwn(row, 'scope'));
      assert.ok(row.labelText === null || row.labelText.startsWith(mode), JSON.stringify(row));
    });
    assert.deepEqual(canonicalFeatureOverrides(overrides), overrides);
  });
});

test('the draft key is the identity pair, and a row that names a mode is invalid', () => {
  const row = { recordKey: 'record-1', biologicalFeatureId: 'f1', placement: { kind: 'main' } };
  assert.equal(featureIdentityKeyOf(row), JSON.stringify(['record-1', 'f1']));
  assert.equal(featureIdentityKeyOf({ recordKey: '', biologicalFeatureId: 'f1' }), '');
  assert.deepEqual(canonicalFeaturePlacements({ [featureIdentityKeyOf(row)]: row }, 'circular'), [row]);
  // The mode is the drawing's, never a row field.
  assert.throws(() => canonicalFeaturePlacements({ [featureIdentityKeyOf(row)]: { scope: 'circular', ...row } }, 'circular'));
  assert.throws(() => canonicalFeaturePlacements({ [JSON.stringify(['circular', 'record-1', 'f1'])]: row }, 'circular'));
  // A lane names a side of its drawing's mode.
  assert.throws(() => canonicalFeaturePlacements({ [featureIdentityKeyOf(row)]: {
    ...row, placement: { kind: 'lane', side: 'above', level: 1 } } }, 'circular'));
});

test('notices, the source-replacing reconcile, and Remove unmatched touch only the drawing they are given', () => {
  MODES.forEach((mode) => {
    const other = mode === 'circular' ? 'linear' : 'circular';
    const notices = recordKeys.map((recordKey) => ({
      recordKey, biologicalFeatureId: 'f1', status: 'unresolved', kinds: ['placement', 'feature_visibility'], resultIndex: 0
    }));
    const drafts = drawings();
    const otherBefore = structuredClone(drafts[other]);
    assert.equal(countUnresolvedFeatureEdits({ ...drafts[mode], notices }), 2 * recordKeys.length);
    assert.equal(removeUnresolvedFeatureEdits({ ...drafts[mode], notices }), 2 * recordKeys.length);
    assert.deepEqual(drafts[other], otherBefore);

    const pruned = drawings();
    const removed = pruneUnmatchedFeatureOverrides({
      ...pruned[mode], notices, replacedRecordKeys: recordKeys,
      previousRecords: [...records, { recordKey: 'dropped' }], currentRecords: records, biologicalFeatures: []
    });
    // The unresolved placements and visibility edits, and the source texts of
    // features the replaced source lost.
    assert.equal(removed, 2 * recordKeys.length);
    Object.values(pruned[mode].featureOverrides).forEach((row) => assert.equal(row.labelSourceText ?? null, null));
    assert.deepEqual(pruned[other], otherBefore);
  });
});

test('Load Feature Edits TSV replaces the edits of its drawing only', () => {
  const drafts = drawings();
  const circularBefore = structuredClone(drafts.circular.featureOverrides);
  replaceFeatureEdits(drafts.linear.featureOverrides, [{
    recordKey: 'record-1', biologicalFeatureId: 'f9', featureVisibility: 'on', labelVisibility: null, labelText: null
  }], records);
  assert.deepEqual(drafts.circular.featureOverrides, circularBefore);
  assert.deepEqual(requestFeatureOverrides(drafts.linear.featureOverrides, records), [{
    recordKey: 'record-1', biologicalFeatureId: 'f9', featureVisibility: 'on', labelVisibility: null, labelText: null
  }]);
});

// A Session 44 draft placement reached every request with its record key (a
// lane failed the other mode's Generate). The migration gives each row the
// modes it reaches, and the split puts it into those drawings (plan 4.2).
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
  const circular = unscopedDraftRows(migrated, 'circular');
  const linear = unscopedDraftRows(migrated, 'linear');
  assert.deepEqual(circular, { [JSON.stringify(['record-1', 'f1'])]: main });
  assert.deepEqual(linear, { [JSON.stringify(['record-1', 'f1'])]: main, [JSON.stringify(['record-1', 'f2'])]: lane });
  canonicalFeaturePlacements(circular, 'circular');
  canonicalFeaturePlacements(linear, 'linear');
});

// A path that indexes a draft map by a hand-built [recordKey, feature] key
// would drift from the draft owner's key: draft keys come only from
// services/feature-placement.js. The Session 44 migration reads and writes the
// old keys, and the exported SVG's runtime indexes the catalog of its one Result.
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
