import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-visibility.js');
const selectorSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-selector.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-visibility-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'feature-visibility.js'),
  await readFile(sourcePath, 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'feature-selector.js'),
  await readFile(selectorSourcePath, 'utf8'),
  'utf8'
);

const {
  buildEditorFeatureVisibilityRule,
  buildExactHashFeatureVisibilityRule,
  buildExactQualifierFeatureVisibilityRule,
  applyFeatureVisibilityOverrideChanges,
  buildFeatureVisibilityChanges,
  buildFeatureVisibilityOverrideCache,
  buildFeatureVisibilitySelectorCache,
  deriveFeatureVisibilityRulesForBoundary,
  exactRegexValue,
  featureVisibilityOverridesToRules,
  featureVisibilityRulesFromOverrideCache,
  getEditorFeatureVisibilityMode,
  getFeatureVisibilityOverride,
  normalizeVisibilityMode,
  parseFeatureVisibilityRules,
  removeEditorFeatureVisibilityRule,
  resolveEffectiveFeatureVisibility,
  serializeFeatureVisibilityRules,
  setFeatureVisibilityOverride,
  splitLegacyVisibilityRules,
  upsertEditorQualifierFeatureVisibilityRule,
  upsertEditorFeatureVisibilityRule
} = await import(pathToFileURL(join(tempDir, 'feature-visibility.js')));

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

{
  const overrides = {};
  assert.equal(getFeatureVisibilityOverride(overrides, 'f.1'), 'default');
  setFeatureVisibilityOverride(overrides, 'f.1', 'off');
  assert.deepEqual(overrides, { 'f.1': 'off' });
  setFeatureVisibilityOverride(overrides, 'f.1', 'default');
  assert.deepEqual(overrides, {});

  const changes = buildFeatureVisibilityChanges(
    [{ svg_id: 'f.1' }, { svg_id: 'f.2' }, { svg_id: 'f.2' }],
    'exclude_matching',
    { 'f.1': 'off' }
  );
  assert.deepEqual(changes, [
    { featureId: 'f.1', before: 'off', after: 'exclude_matching' },
    { featureId: 'f.2', before: 'default', after: 'exclude_matching' }
  ]);
  applyFeatureVisibilityOverrideChanges(overrides, changes.map((change) => ({
    featureId: change.featureId,
    mode: change.after
  })));
  assert.deepEqual(overrides, { 'f.1': 'exclude_matching', 'f.2': 'exclude_matching' });
}

{
  const hashRule = buildExactHashFeatureVisibilityRule(
    { svg_id: 'f.1', label: 'Gene A', type: 'CDS', start: 10, end: 20 },
    'exclude_matching'
  );
  assert.deepEqual(hashRule, {
    id: hashRule.id,
    source: 'editor',
    featureId: 'f.1',
    label: 'Gene A',
    recordId: '*',
    featureType: '*',
    qualifier: 'hash',
    value: '^f\\.1$',
    action: 'exclude_matching'
  });
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
  const feat = {
    svg_id: 'f.1',
    label: 'Gene A',
    record_id: 'rec1',
    type: 'CDS',
    selector: {
      hash: 'f.1',
      record_location: 'rec1:10..20:+',
      qualifiers: {
        protein_id: ['P1'],
        locus_tag: ['L1']
      }
    }
  };
  const rule = buildEditorFeatureVisibilityRule(
    feat,
    {
      selectorSafetyScope: [
        { record_id: 'rec1', feature_type: 'CDS', selector: feat.selector }
      ]
    },
    'exclude_matching'
  );
  assert.equal(rule.source, 'editor');
  assert.equal(rule.featureId, 'f.1');
  assert.equal(rule.recordId, 'rec1');
  assert.equal(rule.featureType, 'CDS');
  assert.equal(rule.qualifier, 'protein_id');
  assert.equal(rule.value, '^P1$');
  assert.equal(rule.action, 'exclude_matching');
}

{
  const feat = {
    svg_id: 'f.1',
    label: 'Gene A',
    record_id: 'rec1',
    type: 'CDS',
    selector: {
      hash: 'f.1',
      record_location: 'rec1:10..20:+',
      qualifiers: { protein_id: ['P1'] }
    }
  };
  const rule = buildEditorFeatureVisibilityRule(feat, {}, 'off');
  assert.equal(rule.qualifier, 'hash');
  assert.equal(rule.value, '^f\\.1$');
}

{
  const feat = {
    svg_id: 'f.1',
    label: 'Gene A',
    record_id: 'rec1',
    type: 'CDS',
    selector: {
      hash: 'f.1',
      record_location: 'rec1:10..20:+',
      qualifiers: {
        protein_id: ['P1']
      }
    }
  };
  const selectorCache = buildFeatureVisibilitySelectorCache(
    [feat],
    [{ record_id: 'rec1', feature_type: 'CDS', selector: feat.selector }]
  );
  assert.deepEqual(selectorCache['f.1'], {
    recordId: 'rec1',
    featureType: 'CDS',
    qualifier: 'protein_id',
    value: 'P1',
    label: 'Gene A'
  });

  const manualRules = [
    { source: 'manual', recordId: '*', featureType: 'CDS', qualifier: 'product', value: 'transposase', action: 'off' }
  ];
  const rules = deriveFeatureVisibilityRulesForBoundary(manualRules, { 'f.1': 'on' }, selectorCache);
  assert.deepEqual(
    rules.map((rule) => [rule.source, rule.featureId, rule.recordId, rule.featureType, rule.qualifier, rule.value, rule.action]),
    [
      ['editor', 'f.1', 'rec1', 'CDS', 'protein_id', '^P1$', 'show'],
      ['manual', '', '*', 'CDS', 'product', 'transposase', 'off']
    ]
  );
  assert.equal(
    serializeFeatureVisibilityRules(rules),
    'rec1\tCDS\tprotein_id\t^P1$\tshow\n*\tCDS\tproduct\ttransposase\toff\n'
  );
  assert.equal(featureVisibilityOverridesToRules({ 'f.missing': 'off' }, selectorCache)[0].qualifier, 'hash');
}

{
  const rules = [
    { source: 'manual', recordId: '*', featureType: 'CDS', qualifier: 'product', value: '.*', action: 'off' }
  ];
  const feat = { svg_id: 'f.1', label: 'Gene A', type: 'CDS', start: 10, end: 20 };
  upsertEditorFeatureVisibilityRule(rules, feat, 'exclude_matching');
  assert.equal(rules[0].source, 'editor');
  assert.equal(rules[0].value, '^f\\.1$');
  assert.equal(rules[0].label, 'Gene A');
  assert.equal(getEditorFeatureVisibilityMode(rules, 'f.1'), 'exclude_matching');
  assert.deepEqual(buildFeatureVisibilityOverrideCache(rules), { 'f.1': 'exclude_matching' });

  upsertEditorFeatureVisibilityRule(rules, feat, 'on');
  assert.equal(rules.length, 2);
  assert.equal(rules[0].action, 'show');

  upsertEditorQualifierFeatureVisibilityRule(
    rules,
    { featureType: 'CDS', qualifier: 'product', value: 'ORF1a polyprotein' },
    'off'
  );
  assert.equal(rules[0].qualifier, 'hash');
  assert.equal(rules[1].qualifier, 'product');
  assert.equal(rules[2].source, 'manual');

  removeEditorFeatureVisibilityRule(rules, 'f.1');
  assert.equal(rules.length, 2);
  assert.equal(getEditorFeatureVisibilityMode(rules, 'f.1'), 'default');
}

{
  const rules = [{
    source: 'editor',
    featureId: 'f.1',
    label: 'Gene A',
    recordId: 'rec1',
    featureType: 'CDS',
    qualifier: 'locus_tag',
    value: '^L1$',
    action: 'off'
  }];
  assert.equal(getEditorFeatureVisibilityMode(rules, 'f.1'), 'off');
  assert.deepEqual(buildFeatureVisibilityOverrideCache(rules), { 'f.1': 'off' });
  removeEditorFeatureVisibilityRule(rules, 'f.1');
  assert.equal(rules.length, 0);
}

assert.deepEqual(
  featureVisibilityRulesFromOverrideCache({ abc: 'off', def: 'suppress', ignored: 'default' })
    .map((rule) => [rule.featureId, rule.action]),
  [['abc', 'off'], ['def', 'exclude_matching']]
);

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

{
  assert.equal(
    resolveEffectiveFeatureVisibility(
      'f.1',
      {},
      null,
      [{ recordId: '*', featureType: '*', qualifier: 'hash', value: '^f\\.1$', action: 'off' }]
    ),
    'off'
  );
  assert.equal(
    resolveEffectiveFeatureVisibility(
      'f.1',
      { 'f.1': 'on' },
      null,
      [{ recordId: '*', featureType: '*', qualifier: 'hash', value: '^f\\.1$', action: 'off' }]
    ),
    'on'
  );
}

console.log('feature visibility tests passed');
