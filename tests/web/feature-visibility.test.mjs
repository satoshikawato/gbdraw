import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-visibility.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-visibility-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempDir, 'feature-visibility.js'),
  await readFile(sourcePath, 'utf8'),
  'utf8'
);

const {
  buildExactHashFeatureVisibilityRule,
  buildExactQualifierFeatureVisibilityRule,
  buildFeatureVisibilityOverrideCache,
  exactRegexValue,
  featureVisibilityRulesFromOverrideCache,
  getEditorFeatureVisibilityMode,
  parseFeatureVisibilityRules,
  removeEditorFeatureVisibilityRule,
  serializeFeatureVisibilityRules,
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

assert.deepEqual(
  featureVisibilityRulesFromOverrideCache({ abc: 'off', def: 'suppress', ignored: 'default' })
    .map((rule) => [rule.featureId, rule.action]),
  [['abc', 'off'], ['def', 'exclude_matching']]
);

console.log('feature visibility tests passed');
