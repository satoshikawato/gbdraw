import assert from 'node:assert/strict';
import { readFile, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourceDir = join(repoRoot, 'gbdraw', 'web', 'js', 'app');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-file-imports-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
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

const { parseColorTable, parsePriorityRules, parseSpecificRules, serializeSpecificRules } = await import(
  pathToFileURL(join(tempDir, 'file-imports.js'))
);

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
