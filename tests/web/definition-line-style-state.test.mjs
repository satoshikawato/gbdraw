import assert from 'node:assert/strict';
import { cp, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-definition-line-style-state-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempDir, 'js'), { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');
const modulePath = join(tempDir, 'js', 'app', 'definition-line-style-state.js');

const {
  createDefaultLinearDefinitionLineStyles,
  normalizeDefinitionLineStyleState
} = await import(pathToFileURL(modulePath));

assert.deepEqual(createDefaultLinearDefinitionLineStyles(), {
  name: { font_size: null, font_weight: null, fill: null },
  subtitle: { font_size: null, font_weight: null, fill: null },
  replicon: { font_size: null, font_weight: null, fill: null },
  accession: { font_size: null, font_weight: null, fill: null },
  length: { font_size: null, font_weight: null, fill: null }
});

assert.deepEqual(
  normalizeDefinitionLineStyleState({
    name: { font_size: '12', font_weight: 'Bold', fill: '#111111' },
    accession: { font_size: '', font_weight: 'Default', fill: '' },
    replicon: { font_weight: 'normal' },
    length: { font_size: '9', font_weight: '700', fill: 'rgb(1,2,3)' }
  }),
  {
    name: { font_size: 12, font_weight: 'bold', fill: '#111111' },
    subtitle: { font_size: null, font_weight: null, fill: null },
    replicon: { font_size: null, font_weight: null, fill: null },
    accession: { font_size: null, font_weight: null, fill: null },
    length: { font_size: 9, font_weight: '700', fill: 'rgb(1,2,3)' }
  }
);
