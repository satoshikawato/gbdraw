import assert from 'node:assert/strict';
import { cp, writeFile, mkdtemp } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-definition-line-style-state-'));
await cp(join(repoRoot, 'gbdraw', 'web', 'js'), join(tempDir, 'js'), { recursive: true });
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}', 'utf8');
const modulePath = join(tempDir, 'js', 'services', 'definition-line-style-state.js');

const {
  createDefaultLinearDefinitionLineStyles,
  normalizeDefinitionLineStyleState,
  refusedDefinitionLineFontSizePath
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

// OV-1270 (PD-OI-081): Generate names a font size of zero or less by its
// setting path instead of drawing it as Auto; Auto and positive sizes pass.
assert.equal(refusedDefinitionLineFontSizePath({ name: { font_size: '-1' } }), 'objects.definition.linear.line_styles.name.font_size');
assert.equal(refusedDefinitionLineFontSizePath({ name: { font_size: '12' }, accession: { font_size: 0 } }), 'objects.definition.linear.line_styles.accession.font_size');
assert.equal(refusedDefinitionLineFontSizePath({ name: { font_size: '12', font_weight: 'heavy' }, length: { font_size: '' } }), '');
assert.equal(refusedDefinitionLineFontSizePath(null), '');
