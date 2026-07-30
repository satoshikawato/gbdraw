import assert from 'node:assert/strict';
import test from 'node:test';

import {
  featureOverrideKey,
  getFeatureOverride,
  migrateLegacyFeatureOverrides
} from '../../gbdraw/web/js/services/feature-override-identity.js';

const feature = (recordKey, biologicalFeatureId, aliases = {}) => ({
  recordKey,
  biologicalFeatureId,
  ...aliases
});

test('feature overrides use the compound biological identity', () => {
  const first = feature('record-a', 'feature-1', {
    id: 'f0',
    svg_id: 'rendered-a'
  });
  const second = feature('record-b', 'feature-1', {
    id: 'f0',
    svg_id: 'rendered-b'
  });

  assert.equal(featureOverrideKey(first), 'record-a\u0000feature-1');
  assert.equal(featureOverrideKey(second), 'record-b\u0000feature-1');
  assert.notEqual(featureOverrideKey(first), featureOverrideKey(second));
});

test('feature overrides follow stable record keys when records are reordered', () => {
  const first = feature('record-a', 'feature-1', {
    fileIdx: 0,
    svg_id: 'feature-1_record_1'
  });
  const second = feature('record-b', 'feature-2', {
    fileIdx: 1,
    svg_id: 'feature-2_record_2'
  });
  const overrides = {
    [featureOverrideKey(first)]: { color: '#ff00ff', strokeWidth: 3 }
  };
  const reordered = [
    { ...second, fileIdx: 0, svg_id: 'feature-2_record_1' },
    { ...first, fileIdx: 1, svg_id: 'feature-1_record_2' }
  ];

  assert.equal(getFeatureOverride(overrides, reordered[0]), undefined);
  assert.deepEqual(getFeatureOverride(overrides, reordered[1]), {
    color: '#ff00ff',
    strokeWidth: 3
  });
});

test('a unique version-39 alias migrates to the compound key', () => {
  const overrides = {
    f0: { color: '#ff00ff' }
  };
  const current = feature('record-a', 'feature-1', {
    svg_id: 'rendered-a',
    stable_svg_id: 'stable-a',
    record_id: 'accession-a'
  });
  const legacy = {
    id: 'f0',
    stable_svg_id: 'stable-a',
    record_id: 'accession-a'
  };

  const result = migrateLegacyFeatureOverrides(overrides, [current], {
    legacyFeatures: [legacy],
    warn: () => assert.fail('a unique migration must not warn')
  });

  assert.deepEqual(overrides, {
    'record-a\u0000feature-1': { color: '#ff00ff' }
  });
  assert.deepEqual(result, {
    migrated: 1,
    unresolved: 0,
    ambiguous: 0,
    collisions: 0
  });
  assert.deepEqual(getFeatureOverride(overrides, current), {
    color: '#ff00ff'
  });
});

test('a positional alias without a matching legacy catalog stays unresolved', () => {
  const overrides = {
    f0: { color: '#ff00ff' }
  };
  const unrelatedCurrent = feature('record-b', 'feature-2', {
    svg_id: 'rendered-b',
    stable_svg_id: 'stable-b'
  });
  const warnings = [];

  const result = migrateLegacyFeatureOverrides(overrides, [unrelatedCurrent], {
    warn: (message) => warnings.push(message)
  });

  assert.equal(result.migrated, 0);
  assert.equal(result.unresolved, 1);
  assert.deepEqual(overrides, {
    f0: { color: '#ff00ff' }
  });
  assert.equal(getFeatureOverride(overrides, unrelatedCurrent), undefined);
  assert.equal(warnings.length, 1);
});

test('unmatched, ambiguous, and colliding legacy keys remain unresolved', () => {
  const first = feature('record-a', 'feature-1', {
    svg_id: 'shared-rendered',
    stable_svg_id: 'stable-a'
  });
  const second = feature('record-b', 'feature-2', {
    svg_id: 'shared-rendered',
    stable_svg_id: 'stable-b'
  });
  const firstKey = featureOverrideKey(first);
  const overrides = {
    missing: { color: '#111111' },
    'shared-rendered': { color: '#222222' },
    'stable-a': { color: '#333333' },
    [firstKey]: { color: '#444444' }
  };
  const warnings = [];

  const result = migrateLegacyFeatureOverrides(overrides, [first, second], {
    warn: (message) => warnings.push(message)
  });

  assert.deepEqual(result, {
    migrated: 0,
    unresolved: 1,
    ambiguous: 1,
    collisions: 1
  });
  assert.equal(overrides.missing.color, '#111111');
  assert.equal(overrides['shared-rendered'].color, '#222222');
  assert.equal(overrides['stable-a'].color, '#333333');
  assert.equal(overrides[firstKey].color, '#444444');
  assert.equal(getFeatureOverride(overrides, first).color, '#444444');
  assert.equal(getFeatureOverride(overrides, second), undefined);
  assert.equal(warnings.length, 1);
  assert.match(warnings[0], /3 legacy feature overrides remain unresolved/);
});

test('legacy-only feature state remains readable without inventing an identity', () => {
  const legacy = { id: 'f7', svg_id: 'rendered-7' };
  const overrides = { f7: { strokeWidth: 5 } };

  assert.equal(featureOverrideKey(legacy), 'f7');
  assert.deepEqual(getFeatureOverride(overrides, legacy), { strokeWidth: 5 });
});

test('a version-39 linear file key correlates through the prior stable feature', () => {
  const overrides = {
    file0_f2: { color: '#abcdef' }
  };
  const current = feature('persistent-record-a', 'biological-a', {
    stableFeatureId: 'stable-a',
    record_id: 'accession-a',
    svg_id: 'rendered-a'
  });
  const legacy = {
    id: 'file0_f2',
    fileIdx: 0,
    stable_svg_id: 'stable-a',
    record_id: 'accession-a',
    svg_id: 'stable-a_record_1'
  };

  const result = migrateLegacyFeatureOverrides(overrides, [current], {
    legacyFeatures: [legacy],
    warn: () => assert.fail('a unique stable correlation must not warn')
  });

  assert.equal(result.migrated, 1);
  assert.deepEqual(overrides[featureOverrideKey(current)], {
    color: '#abcdef'
  });
  assert.equal(overrides.file0_f2, undefined);
});

test('duplicate-accession legacy correlations remain ambiguous without record keys', () => {
  const overrides = {
    file0_f2: { color: '#abcdef' }
  };
  const current = [
    feature('persistent-record-a', 'biological-a', {
      stableFeatureId: 'shared-stable',
      record_id: 'duplicated-accession'
    }),
    feature('persistent-record-b', 'biological-b', {
      stableFeatureId: 'shared-stable',
      record_id: 'duplicated-accession'
    })
  ];
  const warnings = [];

  const result = migrateLegacyFeatureOverrides(overrides, current, {
    legacyFeatures: [{
      id: 'file0_f2',
      fileIdx: 0,
      stable_svg_id: 'shared-stable',
      record_id: 'duplicated-accession'
    }],
    warn: (message) => warnings.push(message)
  });

  assert.equal(result.ambiguous, 1);
  assert.equal(overrides.file0_f2.color, '#abcdef');
  assert.equal(warnings.length, 1);
});
