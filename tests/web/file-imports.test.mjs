import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp, mkdir } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-file-imports-'));
const tempDir = join(tempRoot, 'app');
await mkdir(tempDir);
await mkdir(join(tempRoot, 'services'));
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempRoot, 'services', 'error-normalization.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'error-normalization.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'file-imports.js'),
  await readFile(join(sourceDir, 'file-imports.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'color-utils.js'),
  await readFile(join(sourceDir, 'color-utils.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'specific-color-rules.js'),
  await readFile(join(sourceDir, 'specific-color-rules.js'), 'utf8'),
  'utf8'
);

const { parseColorTable, parsePriorityRules, parseSpecificRules, serializeSpecificRules } = await import(
  pathToFileURL(join(tempDir, 'file-imports.js'))
);
const {
  applySpecificRuleProvenance,
  buildLegendIntents,
  diffLegendIntents,
  prepareSpecificColorImport
} = await import(pathToFileURL(join(tempDir, 'specific-color-rules.js')));

const rules = [
  { feat: 'CDS', qual: 'custom_annotation', val: '^alpha$', color: '#111111', cap: 'Alpha' },
  { feat: 'tRNA', qual: 'color', val: '^beta$', color: '#222222', cap: 'Beta' }
];

assert.equal(
  serializeSpecificRules(rules),
  'CDS\tcustom_annotation\t^alpha$\t#111111\tAlpha\n' +
    'tRNA\tcolor\t^beta$\t#222222\tBeta\n'
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
  { code: 'TABLE_INVALID', context: { row: 1, reason: 'SPECIFIC_COLUMNS' } }
);
assert.throws(
  () => parseSpecificRules('feature_type\n'),
  { code: 'TABLE_INVALID', context: { row: 1, reason: 'SPECIFIC_COLUMNS' } }
);
// The TSV codec preserves Python patterns; validation belongs to the Worker.
assert.equal(parseSpecificRules('CDS\tgene\t(?i)alpha\t#112233\tAlpha\n').rules[0].val, '(?i)alpha');
assert.throws(
  () => parseSpecificRules('CDS\tgene\talpha\tnot-a-color\tAlpha\n'),
  { code: 'TABLE_INVALID', context: { row: 1, field: 'color', reason: 'COLOR' } }
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
// Import assembles the full draft; allocation belongs to Python preparation.
assert.deepEqual(prepareSpecificColorImport(
  'CDS\tgene\talpha\t#112233\tAlpha\nCDS\tgene\tbeta\t#445566\tAlpha\n', []
).nextRules.map(rule => rule.cap), ['Alpha', 'Alpha']);

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

{
  // G-D: the Specific-color domain is shared with Python's table reader (FE-12, D-39).
  const domain = JSON.parse(await readFile(
    join(repoRoot, 'tests', 'fixtures', 'specific_color_domain.json'),
    'utf8'
  ));
  const parseColor = (value) => parseSpecificRules(`CDS\tproduct\tx\t${value}\tcap\n`).rules[0].color;
  const rejects = (value) => assert.throws(
    () => parseColor(value),
    { code: 'TABLE_INVALID', context: { row: 1, field: 'color', reason: 'COLOR' } },
    value
  );
  // DOM-free JavaScript cannot resolve a color name; Python validates it.
  [...domain.valid, ...domain.invalid].forEach((entry) => {
    if (entry.word) assert.equal(parseColor(entry.value), entry.value.toLowerCase(), entry.value);
    else if (entry.normalized) assert.equal(parseColor(entry.value), entry.normalized, entry.value);
    else rejects(entry.value);
  });
  const noneRule = { feat: 'CDS', qual: 'hash', val: 'fx', color: 'none', cap: 'hollow' };
  assert.deepEqual(
    parseSpecificRules(serializeSpecificRules([noneRule])).rules
      .map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })),
    [noneRule]
  );

  const browserColors = new Map(domain.valid
    .filter((entry) => entry.browser)
    .map((entry) => [entry.value.toLowerCase(), entry.browser]));
  let fillStyle = '#000000';
  const context = {
    get fillStyle() { return fillStyle; },
    set fillStyle(value) {
      const text = String(value).toLowerCase();
      if (browserColors.has(text)) fillStyle = browserColors.get(text);
      else if (/^#[0-9a-f]{6}$/.test(text)) fillStyle = text;
    }
  };
  globalThis.document = { createElement: () => ({ getContext: () => context }) };
  try {
    domain.valid.forEach((entry) => {
      assert.equal(parseColor(entry.value), entry.browser || entry.normalized, entry.value);
    });
    domain.invalid.forEach((entry) => rejects(entry.value));
  } finally {
    delete globalThis.document;
  }
}
