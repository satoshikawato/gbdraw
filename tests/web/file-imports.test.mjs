import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp, mkdir } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'services');
const tempRoot = await mkdtemp(join(tmpdir(), 'gbdraw-file-imports-'));
const tempDir = join(tempRoot, 'services');
await mkdir(tempDir);
await mkdir(join(tempRoot, 'utils'));
await writeFile(
  join(tempRoot, 'utils', 'tsv-cell.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'tsv-cell.js'), 'utf8'),
  'utf8'
);
await writeFile(join(tempRoot, 'package.json'), '{"type":"module"}\n', 'utf8');
await writeFile(
  join(tempRoot, 'utils', 'error-normalization.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'error-normalization.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'file-imports.js'),
  await readFile(join(sourceDir, 'file-imports.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempRoot, 'utils', 'color-utils.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'color-utils.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempRoot, 'utils', 'named-colors.js'),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'named-colors.js'), 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'specific-color-rules.js'),
  await readFile(join(sourceDir, 'specific-color-rules.js'), 'utf8'),
  'utf8'
);

const { parseColorTable, parsePriorityRules, parseSpecificRules, parseWhitelistRules, serializeSpecificRules } = await import(
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
// The CLI and Python store the label whitelist with a header row.
assert.deepEqual(
  parseWhitelistRules('feature_type\tqualifier\tkeyword\nCDS\tproduct\tterminase\n').rules,
  [{ feat: 'CDS', qual: 'product', key: 'terminase' }]
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
  const parsesTheDomain = () => {
    domain.valid.forEach((entry) => assert.equal(parseColor(entry.value), entry.normalized, entry.value));
    domain.invalid.forEach((entry) => rejects(entry.value));
  };
  parsesTheDomain();
  const noneRule = { feat: 'CDS', qual: 'hash', val: 'fx', color: 'none', cap: 'hollow' };
  assert.deepEqual(
    parseSpecificRules(serializeSpecificRules([noneRule])).rules
      .map(({ feat, qual, val, color, cap }) => ({ feat, qual, val, color, cap })),
    [noneRule]
  );

  // OV-271: a browser canvas, which also reads currentColor and system colors
  // as Chromium does, does not widen the domain.
  const browserColors = new Map([
    ['red', '#ff0000'], ['grey', '#808080'], ['darkgrey', '#a9a9a9'],
    ['currentcolor', '#000000'], ['buttonface', '#efefef'], ['canvas', '#ffffff']
  ]);
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
    parsesTheDomain();
  } finally {
    delete globalThis.document;
  }
}

// OV-28: the import parsers reject a row of the wrong width like Python does
// (OV-25), with the diagnostic the Label override parser uses.
const wrongWidthDiagnostic = (row, columnCount) => ({
  code: 'TABLE_INVALID',
  context: { row, columnCount, reason: 'FIELDS' }
});
for (const [name, parse, width, valid] of [
  ['Label whitelist', parseWhitelistRules, 3, 'CDS\tproduct\tkinase'],
  ['Qualifier priority', parsePriorityRules, 2, 'CDS\tgene,product'],
  ['Default colors', parseColorTable, 2, 'CDS\t#54bcf8']
]) {
  const cells = valid.split('\t');
  const longRow = [...cells, 'extra'].join('\t');
  const shortRow = cells.slice(0, -1).join('\t');
  assert.throws(() => parse(`${valid}\n${longRow}\n`), wrongWidthDiagnostic(2, width), `${name}: extra column`);
  assert.throws(() => parse(`${longRow}\n${valid}\n`), wrongWidthDiagnostic(1, width), `${name}: extra column on the first row`);
  assert.throws(() => parse(`${valid}\n\n${shortRow}\n`), wrongWidthDiagnostic(3, width), `${name}: short row`);
  assert.throws(() => parse(`${valid}\r\n${longRow}\r\n`), wrongWidthDiagnostic(2, width), `${name}: CRLF`);
  assert.doesNotThrow(() => parse(`\n${valid}\n\n${valid}\r\n`), `${name}: valid rows and blank lines`);
}
// A whitelist keyword may be empty; the row still has three cells.
assert.deepEqual(parseWhitelistRules('CDS\tproduct\t\n').rules, [{ feat: 'CDS', qual: 'product', key: '' }]);

// OV-40: Sessions 31–39 stored these tables as their Web writer wrote them,
// before cell values were normalized, so a tab in a value made extra cells and a
// line break a short row. With legacyRows a parser reads such a row as the
// current writer writes it: the extra cells join the last field with one space,
// and a row without a required field is dropped. It lists each such row; every
// other caller keeps the strict reader above.
for (const [name, parse, text, read, expected, repairs, valid] of [
  ['Label whitelist', parseWhitelistRules,
    'CDS\tproduct\tkinase\nCDS\tproduct\tDNA\t\tpolymerase \npolymerase\n\tproduct\tATP\nCDS\tproduct\t\n',
    (parsed) => parsed.rules,
    [{ feat: 'CDS', qual: 'product', key: 'kinase' }, { feat: 'CDS', qual: 'product', key: 'DNA polymerase' },
      { feat: 'CDS', qual: 'product', key: '' }],
    [{ row: 2, repair: 'joined' }, { row: 3, repair: 'dropped' }, { row: 4, repair: 'dropped' }], 'CDS\tproduct\tkinase'],
  ['Qualifier priority', parsePriorityRules,
    'feature_type\tpriorities\nCDS\tgene,\tproduct\nrRNA\n\tproduct\ntRNA\t\n',
    (parsed) => parsed.rules,
    [{ feat: 'CDS', order: 'gene, product' }],
    [{ row: 2, repair: 'joined' }, { row: 3, repair: 'dropped' }, { row: 4, repair: 'dropped' }, { row: 5, repair: 'dropped' }],
    'CDS\tgene,product'],
  ['Default colors', parseColorTable,
    'feature_type\tcolor\nCDS\t#54bcf8\ntRNA\tred\t\nrRNA\n\tblue\n',
    (parsed) => parsed.colors,
    // A color name reads as its CSS hex, as in the browser (OV-160).
    { CDS: '#54bcf8', tRNA: '#FF0000' },
    [{ row: 3, repair: 'joined' }, { row: 4, repair: 'dropped' }, { row: 5, repair: 'dropped' }], 'CDS\t#54bcf8']
]) {
  const parsed = parse(text, { legacyRows: true });
  assert.deepEqual(read(parsed), expected, `${name}: legacy rows`);
  assert.deepEqual(parsed.repairs, repairs, `${name}: repaired rows`);
  assert.throws(() => parse(text), { code: 'TABLE_INVALID' }, `${name}: strict reader`);
  assert.deepEqual(parse(`${valid}\n`, { legacyRows: true }), parse(`${valid}\n`), `${name}: a valid row reads the same`);
}

// OV-26/OV-27: Python skips whole-line comments in these tables too, so the Web
// parsers and Python agree: an indented or plain `#` line is skipped, a `#` or `"`
// inside a value is data.
assert.deepEqual(
  parseSpecificRules('# note\n  # indented\n\t# tab\nCDS\tproduct\tGene #1\t#ff0000\tCaption #2\n').rules.map((rule) => [rule.val, rule.cap]),
  [['Gene #1', 'Caption #2']]
);
assert.deepEqual(parseColorTable('# note\n  # indented\n\t# tab\nCDS\t#ff0000\n').colors, { CDS: '#ff0000' });
assert.deepEqual(
  parsePriorityRules('# note\n  # indented\n\t# tab\nCDS\t"note#1\n').rules,
  [{ feat: 'CDS', order: '"note#1' }]
);
