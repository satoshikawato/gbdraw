import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'preview-runtime.js');
const featureDomSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-dom.js');
const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-preview-runtime-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app'), { recursive: true });
await writeFile(join(tempDir, 'app', 'preview-runtime.js'), await readFile(sourcePath, 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-dom.js'), await readFile(featureDomSourcePath, 'utf8'), 'utf8');

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
const state = {
  results: ref([
    { name: 'one.svg', content: '<svg id="old-one"></svg>' },
    { name: 'two.svg', content: '<svg id="old-two"></svg>' }
  ]),
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
assert.equal(state.results.value[0].content, '<svg id="old-one"></svg>');
assert.equal(serializeCount, 0);
assert.equal(state.skipCaptureBaseConfig.value, false);

assert.equal(runtime.flushActiveResult(), true);
assert.equal(serializeCount, 1);
assert.equal(state.results.value[0].content, '<svg data-count="1" data-elements="3"></svg>');
assert.equal(state.skipCaptureBaseConfig.value, true);
state.skipCaptureBaseConfig.value = false;
assert.equal(runtime.flushActiveResult(), false);
assert.equal(serializeCount, 1);

runtime.applyFeatureFillChanges([{ featureId: 'feature-b', color: '#abcdef' }]);
assert.equal(featureBBlock.getAttribute('fill'), '#abcdef');
assert.equal(featureBConnector.getAttribute('fill'), 'none');
runtime.applyFeatureVisibilityChanges([{ featureId: 'feature-b', mode: 'off' }]);
assert.equal(featureBBlock.getAttribute('display'), 'none');
assert.equal(featureBConnector.getAttribute('display'), 'none');
runtime.selectResult(1);
assert.equal(serializeCount, 2);
assert.equal(state.skipCaptureBaseConfig.value, false);
assert.equal(state.selectedResultIndex.value, 1);
assert.equal(state.results.value[0].content, '<svg data-count="2" data-elements="3"></svg>');

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
