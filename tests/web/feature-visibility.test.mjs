import assert from 'node:assert/strict';
import {
  buildExactQualifierFeatureVisibilityRule,
  applyFeatureVisibilityOverrideChanges,
  buildFeatureVisibilityChanges,
  exactRegexValue,
  getFeatureVisibilityOverride,
  normalizeVisibilityMode,
  parseFeatureVisibilityRules,
  serializeFeatureVisibilityRules,
  setFeatureVisibilityOverride,
  splitLegacyVisibilityRules,
  upsertEditorQualifierFeatureVisibilityRule
} from '../../gbdraw/web/js/app/feature-visibility.js';

assert.deepEqual(
  parseFeatureVisibilityRules('*\tCDS\tgene\t^geneA$\toff\n').rules.map((rule) => ({
    recordId: rule.recordId,
    featureType: rule.featureType,
    qualifier: rule.qualifier,
    value: rule.value,
    action: rule.action
  })),
  [{ recordId: '*', featureType: 'CDS', qualifier: 'gene', value: '^geneA$', action: 'off' }]
);

assert.deepEqual(
  parseFeatureVisibilityRules(
    '# comment\n\nrecord_id\tfeature_type\tqualifier\tvalue\taction\n' +
      '*\t*\tproduct\ttransposase\ton\n' +
      'rec1\tmisc_feature\tnote\tpseudo\tsuppress\n' +
      'rec2\tCDS\tprotein_id\t^YP_009725295\\.1$\texclude_matching\n'
  ).rules.map((rule) => rule.action),
  ['show', 'exclude_matching', 'exclude_matching']
);

assert.throws(
  () => parseFeatureVisibilityRules('*\tCDS\tgene\t^geneA$\n'),
  /Missing feature visibility columns/
);
assert.throws(
  () => parseFeatureVisibilityRules('*\tCDS\tgene\t^geneA$\thide\textra\n'),
  /Malformed feature visibility row/
);
assert.throws(
  () => parseFeatureVisibilityRules('*\tCDS\tgene\t^geneA$\tmaybe\n'),
  /Invalid feature visibility action/
);
assert.throws(
  () => parseFeatureVisibilityRules('*\tCDS\tgene\t^geneA$\texclude\n'),
  /Invalid feature visibility action/
);

assert.equal(
  serializeFeatureVisibilityRules([
    { recordId: '*', featureType: 'CDS', qualifier: 'gene', value: '^b$', action: 'show' },
    { recordId: '*', featureType: 'CDS', qualifier: 'gene', value: '^a$', action: 'exclude_matching' }
  ]),
  '*\tCDS\tgene\t^b$\tshow\n*\tCDS\tgene\t^a$\texclude_matching\n'
);

assert.equal(exactRegexValue('YP_009725295.1'), '^YP_009725295\\.1$');

assert.equal(normalizeVisibilityMode('suppress'), 'exclude_matching');
assert.equal(normalizeVisibilityMode('default'), 'default');
assert.equal(normalizeVisibilityMode('bad'), 'default');

// Per-feature visibility is the identity row's `featureVisibility` (design Q4),
// in the mode of the feature's Result (R2).
{
  const feature = (svgId, biologicalFeatureId) => ({
    svg_id: svgId, scope: 'circular', record_key: 'rec', biological_feature_id: biologicalFeatureId
  });
  const key = (biologicalFeatureId) => JSON.stringify(['circular', 'rec', biologicalFeatureId]);
  const overrides = {};
  assert.equal(getFeatureVisibilityOverride(overrides, feature('f.1_record_1', 'b1')), 'default');
  setFeatureVisibilityOverride(overrides, feature('f.1_record_1', 'b1'), 'off');
  assert.deepEqual(overrides, { [key('b1')]: {
    scope: 'circular', recordKey: 'rec', biologicalFeatureId: 'b1', featureVisibility: 'off',
    labelVisibility: null, labelText: null, labelSourceText: null
  } });
  // The same identity drawn with another rendered ID reads the same edit (OV-12).
  assert.equal(getFeatureVisibilityOverride(overrides, feature('f.1_record_2', 'b1')), 'off');
  // The same identity in a Result of the other mode is another draft row.
  assert.equal(getFeatureVisibilityOverride(overrides, { ...feature('f.1', 'b1'), scope: 'linear' }), 'default');
  setFeatureVisibilityOverride(overrides, feature('f.1_record_1', 'b1'), 'default');
  assert.deepEqual(overrides, {});

  setFeatureVisibilityOverride(overrides, feature('f.1', 'b1'), 'off');
  const changes = buildFeatureVisibilityChanges(
    [feature('f.1', 'b1'), feature('f.2', 'b2'), feature('f.2', 'b2')],
    'exclude_matching',
    overrides
  );
  assert.deepEqual(changes, [
    { scope: 'circular', recordKey: 'rec', biologicalFeatureId: 'b1', featureId: 'f.1', before: 'off', after: 'exclude_matching' },
    { scope: 'circular', recordKey: 'rec', biologicalFeatureId: 'b2', featureId: 'f.2', before: 'default',
      after: 'exclude_matching' }
  ]);
  applyFeatureVisibilityOverrideChanges(overrides, changes.map((change) => ({ ...change, mode: change.after })));
  assert.deepEqual(Object.values(overrides).map((row) => [row.biologicalFeatureId, row.featureVisibility]), [
    ['b1', 'exclude_matching'], ['b2', 'exclude_matching']
  ]);
}

{
  const rule = buildExactQualifierFeatureVisibilityRule({
    featureType: 'CDS',
    qualifier: 'protein_id',
    value: 'YP_009725295.1',
    action: 'off',
    label: 'Exact protein'
  });
  assert.equal(rule.source, 'editor');
  assert.equal(rule.featureType, 'CDS');
  assert.equal(rule.qualifier, 'protein_id');
  assert.equal(rule.value, '^YP_009725295\\.1$');
  assert.equal(rule.action, 'off');
}

{
  const rules = [
    { source: 'manual', recordId: '*', featureType: 'CDS', qualifier: 'product', value: '.*', action: 'off' }
  ];
  upsertEditorQualifierFeatureVisibilityRule(
    rules,
    { featureType: 'CDS', qualifier: 'product', value: 'ORF1a polyprotein' },
    'off'
  );
  assert.equal(rules[0].qualifier, 'product');
  assert.equal(rules[1].source, 'manual');
}

{
  const split = splitLegacyVisibilityRules([
    {
      source: 'editor',
      featureId: 'f.1',
      label: 'Gene A',
      recordId: 'rec1',
      featureType: 'CDS',
      qualifier: 'protein_id',
      value: '^P1$',
      action: 'off'
    },
    { source: 'manual', recordId: '*', featureType: 'CDS', qualifier: 'product', value: '.*', action: 'show' },
    { source: 'editor', recordId: '*', featureType: 'CDS', qualifier: 'product', value: '^ORF1$', action: 'off' }
  ]);
  assert.deepEqual(split.overrides, { 'f.1': 'off' });
  assert.equal(split.manualRules.length, 2);
  assert.deepEqual(split.manualRules.map((rule) => [rule.source, rule.featureId, rule.qualifier]), [
    ['manual', '', 'product'],
    ['editor', '', 'product']
  ]);
}

console.log('feature visibility tests passed');
