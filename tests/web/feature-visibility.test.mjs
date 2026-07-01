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
  buildFeatureVisibilityOverrideCache,
  featureVisibilityRulesFromOverrideCache,
  getEditorFeatureVisibilityMode,
  parseFeatureVisibilityRules,
  removeEditorFeatureVisibilityRule,
  serializeFeatureVisibilityRules,
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
  [{ recordId: '*', featureType: 'CDS', qualifier: 'gene', value: '^geneA$', action: 'hide' }]
);

assert.deepEqual(
  parseFeatureVisibilityRules(
    '# comment\n\nrecord_id\tfeature_type\tqualifier\tvalue\taction\n' +
      '*\t*\tproduct\ttransposase\tinclude\n' +
      'rec1\tmisc_feature\tnote\tpseudo\texclude\n'
  ).rules.map((rule) => rule.action),
  ['show', 'suppress']
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

assert.equal(
  serializeFeatureVisibilityRules([
    { recordId: '*', featureType: 'CDS', qualifier: 'gene', value: '^b$', action: 'show' },
    { recordId: '*', featureType: 'CDS', qualifier: 'gene', value: '^a$', action: 'suppress' }
  ]),
  '*\tCDS\tgene\t^b$\tshow\n*\tCDS\tgene\t^a$\tsuppress\n'
);

{
  const rules = [
    { source: 'manual', recordId: '*', featureType: 'CDS', qualifier: 'product', value: '.*', action: 'hide' }
  ];
  const feat = { svg_id: 'f.1', label: 'Gene A', type: 'CDS', start: 10, end: 20 };
  upsertEditorFeatureVisibilityRule(rules, feat, 'suppress');
  assert.equal(rules[0].source, 'editor');
  assert.equal(rules[0].value, '^f\\.1$');
  assert.equal(rules[0].label, 'Gene A');
  assert.equal(getEditorFeatureVisibilityMode(rules, 'f.1'), 'suppress');
  assert.deepEqual(buildFeatureVisibilityOverrideCache(rules), { 'f.1': 'suppress' });

  upsertEditorFeatureVisibilityRule(rules, feat, 'on');
  assert.equal(rules.length, 2);
  assert.equal(rules[0].action, 'show');

  removeEditorFeatureVisibilityRule(rules, 'f.1');
  assert.equal(rules.length, 1);
  assert.equal(getEditorFeatureVisibilityMode(rules, 'f.1'), 'default');
}

assert.deepEqual(
  featureVisibilityRulesFromOverrideCache({ abc: 'off', def: 'suppress', ignored: 'default' })
    .map((rule) => [rule.featureId, rule.action]),
  [['abc', 'hide'], ['def', 'suppress']]
);

console.log('feature visibility tests passed');
