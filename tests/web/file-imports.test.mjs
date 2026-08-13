import assert from 'node:assert/strict';
import { mkdir, readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const sourceServicesDir = join(repoRoot, 'gbdraw', 'web', 'js', 'services');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-file-imports-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'));
await mkdir(join(tempDir, 'app', 'feature-editor'));
await mkdir(join(tempDir, 'services'));
await writeFile(
  join(tempDir, 'app', 'file-imports.js'),
  await readFile(join(sourceDir, 'file-imports.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'color-utils.js'),
  await readFile(join(sourceDir, 'color-utils.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'specific-color-rules.js'),
  await readFile(join(sourceDir, 'specific-color-rules.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'feature-instance-identity.js'),
  await readFile(join(sourceServicesDir, 'feature-instance-identity.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'json-clone.js'),
  await readFile(join(sourceServicesDir, 'json-clone.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'feature-editor', 'semantic-fill-selectors.js'),
  await readFile(join(sourceDir, 'feature-editor', 'semantic-fill-selectors.js'), 'utf8'),
  'utf8'
);

const { parseColorTable, parsePriorityRules, parseSpecificRules, serializeSpecificRules } = await import(
  pathToFileURL(join(tempDir, 'app', 'file-imports.js'))
);
const {
  applySpecificRuleProvenance,
  buildLegendIntents,
  diffLegendIntents,
  prepareSpecificColorImport
} = await import(pathToFileURL(join(tempDir, 'app', 'specific-color-rules.js')));

const rules = [
  { feat: 'CDS', qual: 'custom_annotation', val: '^alpha$', color: '#111111', cap: 'Alpha' },
  { feat: 'tRNA', qual: 'color', val: '^beta$', color: '#222222', cap: 'Beta' }
];

assert.equal(
  serializeSpecificRules(rules),
  'CDS\tcustom_annotation\t^alpha$\t#111111\tAlpha\n' +
    'tRNA\tcolor\t^beta$\t#222222\tBeta\n'
);

const instanceHash = 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa';
const literalRule = parseSpecificRules(
  `CDS\t__gbdraw_instance_hash__\t${instanceHash}\tnone\tHidden\n`,
  { schema: 6 }
).rules[0];
assert.equal(literalRule.match, 'literal');
assert.equal(literalRule.color, 'none');
assert.equal(
  serializeSpecificRules([literalRule], { schema: 6 }),
  `CDS\t__gbdraw_instance_hash__\t${instanceHash}\tnone\tHidden\n`
);
assert.throws(
  () => serializeSpecificRules([literalRule], { schema: 5 }),
  /require canonical request schema 6/
);
assert.throws(
  () => parseSpecificRules(
    'CDS\t__gbdraw_instance_hash__\tfi1_bad\t#112233\tBad\n',
    { schema: 6 }
  ),
  /Invalid feature-instance hash/
);
assert.equal(
  parseSpecificRules(
    'CDS\t__gbdraw_instance_hash__\tfi1_.*\t#112233\tLegacy biological qualifier\n',
    { schema: 5 }
  ).rules[0].match,
  'regex'
);

const semanticLiteralRule = parseSpecificRules(
  '*\t__gbdraw_semantic_scope__\tfs1:source-annotation-label:Shared%20protein\t#abcdef\tShared protein\n',
  { schema: 6 }
).rules[0];
assert.equal(semanticLiteralRule.match, 'literal');
assert.equal(
  serializeSpecificRules([semanticLiteralRule], { schema: 6 }),
  '*\t__gbdraw_semantic_scope__\tfs1:source-annotation-label:Shared%20protein\t#abcdef\tShared protein\n'
);
assert.throws(
  () => serializeSpecificRules([semanticLiteralRule], { schema: 5 }),
  /require canonical request schema 6/
);
assert.throws(
  () => parseSpecificRules(
    '*\t__gbdraw_semantic_scope__\tfs1:source-annotation-label:%zz\t#abcdef\tBad\n',
    { schema: 6 }
  ),
  /Invalid Feature semantic selector/
);
assert.equal(
  parseSpecificRules(
    '*\t__gbdraw_semantic_scope__\t^bio.*$\t#abcdef\tLegacy biological qualifier\n',
    { schema: 5 }
  ).rules[0].match,
  'regex'
);

assert.deepEqual(
  parseSpecificRules(serializeSpecificRules(rules)).rules.map(({ feat, qual, val, color, cap }) => ({
    feat,
    qual,
    val,
    color,
    cap
  })),
  rules
);

assert.equal(
  serializeSpecificRules([
    { feat: 'CDS', qual: 'product', val: 'line\nbreak', color: '#333333', cap: 'has\ttab' }
  ]),
  'CDS\tproduct\tline break\t#333333\thas tab\n'
);

assert.equal(
  serializeSpecificRules([{ feat: 'CDS', qual: 'product', val: '', color: '#333333', cap: 'Skipped' }]),
  ''
);

assert.deepEqual(
  parseSpecificRules(
    'feature_type\tqualifier_key\tvalue\tcolor\tcaption\nCDS\tgene\tpsaA\t#00662c\tphotosystem I\n'
  ).rules.map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })),
  [{ feat: 'CDS', qual: 'gene', val: 'psaA', color: '#00662c', cap: 'photosystem I' }]
);
assert.deepEqual(
  parsePriorityRules('feature_type\tpriorities\nCDS\tgene,old_locus_tag\n').rules,
  [{ feat: 'CDS', order: 'gene,old_locus_tag' }]
);
assert.deepEqual(
  parseColorTable('feature_type\tcolor\nCDS\t#54bcf8\n').colors,
  { CDS: '#54bcf8' }
);

assert.throws(
  () => parseSpecificRules('CDS\tgene\talpha\n'),
  /line 1: expected 4 or 5 columns/
);
assert.throws(
  () => parseSpecificRules('feature_type\n'),
  /line 1: expected 4 or 5 columns/
);
assert.throws(
  () => parseSpecificRules('CDS\tgene\t[\t#112233\tAlpha\n'),
  /Invalid specific-color regex at line 1/
);
assert.throws(
  () => parseSpecificRules('CDS\tgene\talpha\tnot-a-color\tAlpha\n'),
  /Invalid specific-color value at line 1/
);
assert.equal(
  parseSpecificRules(
    'CDS\tgene\talpha\t#112233\tAlpha\nCDS\tgene\talpha\t#112233\tAlpha\n'
  ).count,
  1
);

const prepared = prepareSpecificColorImport(
  'CDS\tgene\talpha\t#112233\tAlpha\nCDS\tproduct\tbeta\t#112233\tAlpha\n',
  [
    { feat: 'CDS', qual: 'gene', val: 'old', color: '#999999', cap: 'Old', fromFile: true },
    { feat: 'tRNA', qual: 'gene', val: 'manual', color: '#445566', cap: 'Manual' }
  ]
);
assert.equal(prepared.nextRules.length, 3);
assert.equal(prepared.nextRules.filter((rule) => rule.fromFile).length, 2);
assert.deepEqual(prepared.intents, [{ caption: 'Alpha', color: '#112233' }]);
assert.throws(
  () => prepareSpecificColorImport(
    'CDS\tgene\talpha\t#112233\tAlpha\nCDS\tgene\tbeta\t#445566\tAlpha\n',
    []
  ),
  /caption "Alpha" uses multiple colors/
);

assert.deepEqual(diffLegendIntents(
  [
    { caption: 'Keep', color: '#111111' },
    { caption: 'Recolor', color: '#222222' },
    { caption: 'Remove', color: '#333333' }
  ],
  [
    { caption: 'Keep', color: '#111111' },
    { caption: 'Recolor', color: '#abcdef' },
    { caption: 'Add', color: '#444444' }
  ]
), {
  add: [{ caption: 'Add', color: '#444444' }],
  update: [{ caption: 'Recolor', color: '#abcdef' }],
  remove: [{ caption: 'Remove', color: '#333333' }],
  unchanged: [{ caption: 'Keep', color: '#111111' }]
});

const canonicalRules = [{ feat: 'CDS', qual: 'gene', val: 'alpha', color: '#112233', cap: 'Alpha' }];
assert.equal(applySpecificRuleProvenance(canonicalRules, [
  { ...canonicalRules[0], fromFile: true }
])[0].fromFile, true);
assert.equal(applySpecificRuleProvenance(canonicalRules, [
  { ...canonicalRules[0], color: '#ffffff', fromFile: true }
])[0].fromFile, undefined);
assert.deepEqual(buildLegendIntents(canonicalRules).intents, [
  { caption: 'Alpha', color: '#112233' }
]);
