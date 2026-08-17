import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'preview-runtime.js');
const featureDomSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-dom.js');
const serviceNames = [
  'runtime-test-hooks.js',
  'session-feature-metadata.js',
  'svg-result-normalization.js',
  'svg-result-ingestion.js',
  'svg-sanitization.js',
  'svg-serialization.js'
];
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-preview-runtime-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(join(tempDir, 'app', 'preview-runtime.js'), await readFile(sourcePath, 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-dom.js'), await readFile(featureDomSourcePath, 'utf8'), 'utf8');
await Promise.all(serviceNames.map(async (name) => {
  const servicePath = join(repoRoot, 'gbdraw', 'web', 'js', 'services', name);
  await writeFile(join(tempDir, 'services', name), await readFile(servicePath, 'utf8'), 'utf8');
}));

const { createPreviewRuntime } = await import(pathToFileURL(join(tempDir, 'app', 'preview-runtime.js')));

const ref = (value) => ({ value });

class FakeFeatureElement {
  constructor(id, { featureId = id, part = 'block', fill = null } = {}) {
    this.id = id;
    this.attributes = {
      'data-gbdraw-feature-id': featureId,
      'data-gbdraw-feature-part': part
    };
    if (fill !== null) this.attributes.fill = fill;
  }

  getAttribute(name) {
    return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null;
  }

  setAttribute(name, value) {
    this.attributes[name] = String(value);
  }

  removeAttribute(name) {
    delete this.attributes[name];
  }
}

const makeSvg = (elements) => ({
  elements,
  querySelectorAll: () => elements,
  getElementById: (id) => elements.find((element) => element.id === id) || null
});

const featureA = new FakeFeatureElement('feature-a');
const featureBBlock = new FakeFeatureElement('feature-b__part1', {
  featureId: 'feature-b',
  fill: '#111111'
});
const featureBConnector = new FakeFeatureElement('feature-b__line1', {
  featureId: 'feature-b',
  part: 'connector',
  fill: 'none'
});
const svg = makeSvg([featureA, featureBConnector, featureBBlock]);
let resultReplacementCount = 0;
const trackedResults = new Proxy([
  { name: 'one.svg', content: '<svg id="old-one"></svg>' },
  { name: 'two.svg', content: '<svg id="old-two"></svg>' }
], {
  set(target, property, value) {
    if (/^\d+$/.test(String(property))) resultReplacementCount += 1;
    return Reflect.set(target, property, value);
  }
});
const state = {
  results: ref(trackedResults),
  selectedResultIndex: ref(0),
  skipCaptureBaseConfig: ref(false),
  svgContainer: ref({
    querySelector: (selector) => (selector === 'svg' ? svg : null)
  })
};
let serializeCount = 0;
const runtime = createPreviewRuntime({
  state,
  serializeSvg: (targetSvg) => {
    serializeCount += 1;
    return `<svg data-count="${serializeCount}" data-elements="${targetSvg.elements.length}"></svg>`;
  }
});

runtime.mountResultSvg(0, svg);
assert.equal(runtime.applyFeatureVisibilityChanges([{ featureId: 'feature-a', mode: 'off' }]), true);
assert.equal(featureA.getAttribute('display'), 'none');
assert.equal(serializeCount, 1);
assert.equal(resultReplacementCount, 1);
assert.equal(state.results.value[0].content, '<svg data-count="1" data-elements="3"></svg>');
assert.equal(state.skipCaptureBaseConfig.value, true);

state.skipCaptureBaseConfig.value = false;
const resultAfterSingleEdit = state.results.value[0];
assert.equal(runtime.applyFeatureVisibilityChanges([{ featureId: 'feature-a', mode: 'off' }]), false);
assert.equal(serializeCount, 1);
assert.equal(resultReplacementCount, 1);
assert.equal(state.results.value[0], resultAfterSingleEdit);
assert.equal(state.skipCaptureBaseConfig.value, false);

assert.equal(runtime.applyFeatureVisibilityChanges([
  { featureId: 'feature-a', mode: 'on' },
  { featureId: 'feature-b', mode: 'off' }
]), true);
assert.equal(featureA.getAttribute('display'), null);
assert.equal(featureBBlock.getAttribute('display'), 'none');
assert.equal(featureBConnector.getAttribute('display'), 'none');
assert.equal(serializeCount, 2);
assert.equal(resultReplacementCount, 2);

state.skipCaptureBaseConfig.value = false;
const resultAfterBulkEdit = state.results.value[0];
assert.equal(runtime.applyFeatureVisibilityChanges([
  { featureId: 'feature-a', mode: 'on' },
  { featureId: 'feature-b', mode: 'off' }
]), false);
assert.equal(serializeCount, 2);
assert.equal(resultReplacementCount, 2);
assert.equal(state.results.value[0], resultAfterBulkEdit);
assert.equal(state.skipCaptureBaseConfig.value, false);

assert.equal(runtime.applyFeatureFillChanges([{ featureId: 'feature-b', color: '#111111' }]), false);
featureBBlock.setAttribute('stroke', '#222222');
featureBBlock.setAttribute('stroke-width', '3');
featureBConnector.setAttribute('stroke', '#222222');
featureBConnector.setAttribute('stroke-width', '3');
assert.equal(runtime.applyFeatureStrokeChanges([{
  featureId: 'feature-b',
  strokeColor: '#222222',
  strokeWidth: 3
}]), false);
assert.equal(runtime.getActiveRuntime().dirty, false);
assert.equal(runtime.flushActiveResult(), false);
assert.equal(serializeCount, 2);
assert.equal(state.skipCaptureBaseConfig.value, false);
state.skipCaptureBaseConfig.value = false;
assert.equal(runtime.flushActiveResult(), false);
assert.equal(serializeCount, 2);

runtime.applyFeatureFillChanges([{ featureId: 'feature-b', color: '#abcdef' }]);
assert.equal(featureBBlock.getAttribute('fill'), '#abcdef');
assert.equal(featureBConnector.getAttribute('fill'), 'none');
runtime.selectResult(1);
assert.equal(serializeCount, 3);
assert.equal(resultReplacementCount, 3);
assert.equal(state.skipCaptureBaseConfig.value, false);
assert.equal(state.selectedResultIndex.value, 1);
assert.equal(state.results.value[0].content, '<svg data-count="3" data-elements="3"></svg>');

const legacyConnector = new FakeFeatureElement('feature-c__line1', {
  featureId: '',
  part: '',
  fill: 'none'
});
delete legacyConnector.attributes['data-gbdraw-feature-id'];
delete legacyConnector.attributes['data-gbdraw-feature-part'];
const legacyBlock = new FakeFeatureElement('feature-c__part1', {
  featureId: '',
  part: '',
  fill: '#222222'
});
delete legacyBlock.attributes['data-gbdraw-feature-id'];
delete legacyBlock.attributes['data-gbdraw-feature-part'];
const legacySvg = makeSvg([legacyConnector, legacyBlock]);
state.svgContainer.value = { querySelector: (selector) => (selector === 'svg' ? legacySvg : null) };
runtime.mountResultSvg(1, legacySvg);
runtime.applyFeatureFillChanges([{ featureId: 'feature-c', color: '#fedcba' }]);
assert.equal(legacyBlock.getAttribute('fill'), '#fedcba');
assert.equal(legacyConnector.getAttribute('fill'), 'none');

console.log('preview runtime tests passed');
