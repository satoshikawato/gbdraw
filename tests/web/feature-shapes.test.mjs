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
  DEFAULT_ARROW_HEAD_LENGTH_RATIO,
  DEFAULT_ARROWHEAD_SHAFT_WIDTH_RATIO,
  DEFAULT_FEATURE_RENDERINGS,
  FEATURE_RENDERING_VALUES,
  arrowHeadLengthRatioForState,
  createDefaultFeatureRenderings,
  defaultFeatureRendering,
  normalizeArrowHeadLengthRatio,
  normalizeArrowheadShaftWidthRatio,
  normalizeFeatureRendering,
  normalizeFeatureRenderingMap
} = await import(pathToFileURL(modulePath));

assert.deepEqual(FEATURE_RENDERING_VALUES, ['arrow', 'arrowhead', 'rectangle', 'underlay']);
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
assert.equal(normalizeFeatureRendering(' ARROWHEAD '), 'arrowhead');
assert.equal(normalizeFeatureRendering(' UNDERLAY '), 'underlay');
assert.deepEqual(
  normalizeFeatureRenderingMap({ repeat_region: 'RECTANGLE', CDS: 'underlay' }),
  { repeat_region: 'rectangle', CDS: 'underlay' }
);
assert.throws(() => normalizeFeatureRendering('triangle'), /Unsupported feature rendering/);
assert.throws(() => normalizeFeatureRenderingMap({ ' ': 'arrow' }), /must not be empty/);

assert.equal(DEFAULT_ARROW_HEAD_LENGTH_RATIO, 'auto');
assert.equal(DEFAULT_ARROWHEAD_SHAFT_WIDTH_RATIO, 0.5);
assert.equal(normalizeArrowHeadLengthRatio(undefined), 'auto');
assert.equal(normalizeArrowHeadLengthRatio(' AUTO '), 'auto');
assert.equal(normalizeArrowHeadLengthRatio(1.25), 1.25);
assert.equal(arrowHeadLengthRatioForState('auto'), null);
assert.equal(arrowHeadLengthRatioForState(0.75), 0.75);
for (const invalid of [true, false, 0, -1, Number.NaN, Number.POSITIVE_INFINITY, '1.25', [], [1], {}]) {
  assert.throws(() => normalizeArrowHeadLengthRatio(invalid), /positive finite number/);
}
assert.equal(normalizeArrowheadShaftWidthRatio(undefined), 0.5);
assert.equal(normalizeArrowheadShaftWidthRatio(0.25), 0.25);
assert.equal(normalizeArrowheadShaftWidthRatio(1), 1);
for (const invalid of [true, false, 0, -1, 1.01, Number.NaN, Number.POSITIVE_INFINITY, '0.25', 'auto', [], [0.25], {}]) {
  assert.throws(() => normalizeArrowheadShaftWidthRatio(invalid), /at most 1/);
}

const first = createDefaultFeatureRenderings();
const second = createDefaultFeatureRenderings();
first.repeat_region = 'rectangle';
assert.equal(second.repeat_region, 'underlay');
