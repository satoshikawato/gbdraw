import assert from 'node:assert/strict';
import { readdir, readFile } from 'node:fs/promises';
import { join } from 'node:path';
import test from 'node:test';
import { pathToFileURL } from 'node:url';

const jsRoot = join(process.cwd(), 'gbdraw', 'web', 'js');
const { normalizeTsvCell } = await import(pathToFileURL(join(jsRoot, 'utils', 'tsv-cell.js')));
const { serializeLabelWhitelistRules, serializeQualifierPriorityRules } = await import(
  pathToFileURL(join(jsRoot, 'services', 'file-imports.js'))
);

test('a TSV cell never carries a tab or line break and is trimmed', () => {
  assert.equal(normalizeTsvCell('  a\tb\r\nc  '), 'a b c');
  assert.equal(normalizeTsvCell(null), '');
  assert.equal(normalizeTsvCell(undefined), '');
});

test('the Label whitelist writer keeps one row of three cells per rule', () => {
  const text = serializeLabelWhitelistRules([
    { feat: 'CDS', qual: 'product', key: 'two\twords' },
    { feat: ' CDS ', qual: 'gene', key: 'a\nb' },
    { feat: '', qual: 'gene', key: 'dropped' },
    { feat: 'CDS', qual: '\t', key: 'dropped' },
    null
  ]);
  assert.equal(text, 'CDS\tproduct\ttwo words\nCDS\tgene\ta b\n');
  text.trimEnd().split('\n').forEach((line) => assert.equal(line.split('\t').length, 3));
  assert.equal(serializeLabelWhitelistRules([]), '');
});

test('the Qualifier priority writer keeps one row of two cells per rule', () => {
  const text = serializeQualifierPriorityRules([
    { feat: 'CDS', order: 'gene,\tproduct' },
    { feat: 'tRNA', order: '' },
    { feat: '', order: 'gene' }
  ]);
  assert.equal(text, 'CDS\tgene, product\n');
  assert.equal(serializeQualifierPriorityRules([]), '');
});

// Guard: every writer of these tables uses the shared serializers, so a new
// writer cannot interpolate a raw cell (the Generate and the request paths had
// each their own copy that skipped the normalization).
test('no module interpolates a raw whitelist or priority cell into a TSV row', async () => {
  const offenders = [];
  const walk = async (directory) => {
    for (const entry of await readdir(directory, { withFileTypes: true })) {
      const path = join(directory, entry.name);
      if (entry.isDirectory()) await walk(path);
      else if (entry.name.endsWith('.js')) {
        const text = await readFile(path, 'utf8');
        if (/\$\{\s*r(?:ule)?\.(?:feat|qual|key|order)\s*(?:\|\|\s*'')?\}\\t/.test(text)) offenders.push(path);
      }
    }
  };
  await walk(jsRoot);
  assert.deepEqual(offenders, []);
});
