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

const { parseSpecificRules, serializeSpecificRules } = await import(
  pathToFileURL(join(tempDir, 'file-imports.js'))
);

const rules = [
  { feat: 'CDS', qual: 'product', val: '^alpha$', color: '#111111', cap: 'Alpha' },
  { feat: 'tRNA', qual: 'gene', val: '^beta$', color: '#222222', cap: 'Beta' }
];

assert.equal(
  serializeSpecificRules(rules),
  'CDS\tproduct\t^alpha$\t#111111\tAlpha\n' +
    'tRNA\tgene\t^beta$\t#222222\tBeta\n'
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
