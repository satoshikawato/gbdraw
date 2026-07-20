import assert from 'node:assert/strict';
import { mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-feature-shapes-'));
const modulePath = join(tempDir, 'feature-rendering.mjs');
await writeFile(
  modulePath,
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'utils', 'feature-rendering.js'), 'utf8'),
  'utf8'
);

const {
  DEFAULT_FEATURE_RENDERINGS,
  FEATURE_RENDERING_VALUES,
  createDefaultFeatureRenderings,
  defaultFeatureRendering,
  normalizeFeatureRendering,
  normalizeFeatureRenderingMap
} = await import(pathToFileURL(modulePath));

assert.deepEqual(FEATURE_RENDERING_VALUES, ['arrow', 'rectangle', 'underlay']);
assert.deepEqual(DEFAULT_FEATURE_RENDERINGS, {
  CDS: 'arrow',
  rRNA: 'arrow',
  tRNA: 'arrow',
  tmRNA: 'arrow',
  ncRNA: 'arrow',
  misc_RNA: 'arrow',
  repeat_region: 'underlay'
});
assert.equal(defaultFeatureRendering('repeat_region'), 'underlay');
assert.equal(defaultFeatureRendering('misc_feature'), 'rectangle');
assert.equal(normalizeFeatureRendering(' UNDERLAY '), 'underlay');
assert.deepEqual(
  normalizeFeatureRenderingMap({ repeat_region: 'RECTANGLE', CDS: 'underlay' }),
  { repeat_region: 'rectangle', CDS: 'underlay' }
);
assert.throws(() => normalizeFeatureRendering('triangle'), /Unsupported feature rendering/);
assert.throws(() => normalizeFeatureRenderingMap({ ' ': 'arrow' }), /must not be empty/);

const first = createDefaultFeatureRenderings();
const second = createDefaultFeatureRenderings();
first.repeat_region = 'rectangle';
assert.equal(second.repeat_region, 'underlay');
