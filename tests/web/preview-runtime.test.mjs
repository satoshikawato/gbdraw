import assert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const sourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'preview-runtime.js');
const featureDomSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-dom.js');
const legendUtilsSourcePath = join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'legend', 'utils.js');
const legendTransformSourcePath = join(
  repoRoot,
  'gbdraw',
  'web',
  'js',
  'app',
  'legend-layout',
  'transform-utils.js'
);
const serviceNames = [
  'current-worker-result-source.js',
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
await mkdir(join(tempDir, 'app', 'legend'), { recursive: true });
await mkdir(join(tempDir, 'app', 'legend-layout'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });
await writeFile(join(tempDir, 'app', 'preview-runtime.js'), await readFile(sourcePath, 'utf8'), 'utf8');
await writeFile(join(tempDir, 'app', 'feature-dom.js'), await readFile(featureDomSourcePath, 'utf8'), 'utf8');
await writeFile(
  join(tempDir, 'app', 'legend', 'utils.js'),
  await readFile(legendUtilsSourcePath, 'utf8'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'legend-layout', 'transform-utils.js'),
  await readFile(legendTransformSourcePath, 'utf8'),
  'utf8'
);
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
  isConnected: true,
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
let trackedResults = [
  { name: 'one.svg', content: '<svg id="old-one"></svg>' },
  { name: 'two.svg', content: '<svg id="old-two"></svg>' }
];
const resultsRef = {
  get value() { return trackedResults; },
  set value(nextResults) {
    resultReplacementCount += 1;
    trackedResults = nextResults;
  }
};
const state = {
  results: resultsRef,
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

const structuralMetrics = [];
const lifecycleEvents = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric(metric) {
    structuralMetrics.push(metric);
  },
  onSessionLifecycleEvent(event) {
    lifecycleEvents.push(event);
  }
};
const readinessResult = { name: 'ready.svg', content: '<svg></svg>' };
const readinessRoot = makeSvg([]);
const readinessState = {
  results: ref([readinessResult]),
  selectedResultIndex: ref(0),
  featureCatalog: ref({ schema: 3, items: [] }),
  svgContainer: ref({
    querySelector: (selector) => (selector === 'svg' ? readinessRoot : null)
  })
};
const readinessRuntime = createPreviewRuntime({
  state: readinessState,
  serializeSvg: () => '<svg></svg>'
});
const binderCalls = new Map();
const countStep = (name) => {
  binderCalls.set(name, Number(binderCalls.get(name) || 0) + 1);
};
readinessRuntime.configureMountedResultBinder({
  adoptLegend: () => countStep('legend'),
  bindComposition: () => countStep('composition'),
  setupDragAffordances: () => countStep('drag'),
  installDelegatedInteractions: () => countStep('interactions'),
  synchronizeLabelEditor: () => countStep('label'),
  initializeStrokeAndCanvas: () => countStep('stroke-canvas'),
  reconcileSelection: () => countStep('selection')
});
const readyExpectation = readinessRuntime.registerReadinessExpectation({
  result: readinessResult,
  resultIndex: 0,
  artifactIdentity: { fingerprint: 'artifact-ready' },
  generationToken: 'generation-1',
  catalogState: readinessState.featureCatalog.value,
  phase: 'generate',
  isCurrent: () => true
});
const readyContext = readinessRuntime.createMountedResultContext({
  root: readinessRoot,
  result: readinessResult,
  resultIndex: 0
});
const readyReceipt = await readinessRuntime.bindMountedResult(readyContext);
assert.equal(await readyExpectation.promise, readyReceipt);
assert.deepEqual(Object.keys(readyReceipt).sort(), [
  'artifactIdentity',
  'bindSequence',
  'generationToken',
  'phase',
  'readyTimestamp',
  'requiredBindingFlags',
  'resultIdentity',
  'resultIndex',
  'rootGeneration',
  'rootIdentity'
]);
assert.deepEqual(readyReceipt.requiredBindingFlags, {
  rootAdopted: true,
  legendReady: true,
  compositionReady: true,
  dragReady: true,
  interactionsReady: true,
  labelEditorReady: true,
  strokeCanvasReady: true,
  selectionReady: true
});
assert.equal(readinessRuntime.isActiveResultReady(), true);
for (const step of [
  'legend',
  'composition',
  'drag',
  'interactions',
  'label',
  'stroke-canvas',
  'selection'
]) {
  assert.equal(binderCalls.get(step), 1, step);
}

const repeatedContext = readinessRuntime.createMountedResultContext({
  root: readinessRoot,
  result: readinessResult,
  resultIndex: 0
});
assert.equal(await readinessRuntime.bindMountedResult(repeatedContext), readyReceipt);
assert.equal(binderCalls.get('legend'), 1);

const incompleteResult = { name: 'incomplete.svg', content: '<svg></svg>' };
readinessState.results.value = [incompleteResult];
const incompleteRoot = makeSvg([]);
readinessState.svgContainer.value = {
  querySelector: (selector) => (selector === 'svg' ? incompleteRoot : null)
};
const incompleteExpectation = readinessRuntime.registerReadinessExpectation({
  result: incompleteResult,
  resultIndex: 0,
  artifactIdentity: 'artifact-incomplete',
  generationToken: 'generation-incomplete',
  phase: 'generate-incomplete',
  isCurrent: () => true
});
const incompleteContext = readinessRuntime.createMountedResultContext({
  root: incompleteRoot,
  result: incompleteResult,
  resultIndex: 0
});
assert.equal(readinessRuntime.acceptReadyReceipt(Object.freeze({
  artifactIdentity: incompleteContext.artifactIdentity,
  resultIdentity: incompleteContext.resultIdentity,
  resultIndex: incompleteContext.resultIndex,
  generationToken: incompleteContext.generationToken,
  rootGeneration: incompleteContext.rootGeneration,
  bindSequence: incompleteContext.bindSequence,
  requiredBindingFlags: Object.freeze({ rootAdopted: true }),
  phase: incompleteContext.phase
})).accepted, false);
await assert.rejects(incompleteExpectation.promise, /did not complete/);

const canceledExpectation = readinessRuntime.registerReadinessExpectation({
  result: incompleteResult,
  resultIndex: 0,
  artifactIdentity: 'artifact-canceled',
  generationToken: 'generation-canceled',
  phase: 'generate-canceled',
  isCurrent: () => true
});
const canceledContext = readinessRuntime.createMountedResultContext({
  root: incompleteRoot,
  result: incompleteResult,
  resultIndex: 0
});
assert.equal(
  readinessRuntime.invalidateReadinessExpectation(
    'generation-canceled',
    new Error('canceled generation')
  ),
  true
);
await assert.rejects(canceledExpectation.promise, /canceled generation/);
assert.equal(readinessRuntime.acceptReadyReceipt(Object.freeze({
  artifactIdentity: canceledContext.artifactIdentity,
  resultIdentity: canceledContext.resultIdentity,
  resultIndex: canceledContext.resultIndex,
  generationToken: canceledContext.generationToken,
  rootIdentity: canceledContext.rootGeneration,
  rootGeneration: canceledContext.rootGeneration,
  bindSequence: canceledContext.bindSequence,
  requiredBindingFlags: Object.freeze({
    rootAdopted: true,
    legendReady: true,
    compositionReady: true,
    dragReady: true,
    interactionsReady: true,
    labelEditorReady: true,
    strokeCanvasReady: true,
    selectionReady: true
  }),
  readyTimestamp: Date.now(),
  phase: canceledContext.phase
})).accepted, false);

const staleResult = { name: 'stale.svg', content: '<svg></svg>' };
const staleRoot = makeSvg([]);
readinessState.results.value = [staleResult];
readinessState.svgContainer.value = {
  querySelector: (selector) => (selector === 'svg' ? staleRoot : null)
};
const staleExpectation = readinessRuntime.registerReadinessExpectation({
  result: staleResult,
  resultIndex: 0,
  artifactIdentity: 'artifact-stale',
  generationToken: 'generation-stale',
  phase: 'generate-stale',
  isCurrent: () => false
});
const staleContext = readinessRuntime.createMountedResultContext({
  root: staleRoot,
  result: staleResult,
  resultIndex: 0
});
await assert.rejects(readinessRuntime.bindMountedResult(staleContext), /stale or canceled/);
await assert.rejects(staleExpectation.promise, /stale or canceled/);

const failedResult = { name: 'failed.svg', content: '<svg></svg>' };
const failedRoot = makeSvg([]);
readinessState.results.value = [failedResult];
readinessState.svgContainer.value = {
  querySelector: (selector) => (selector === 'svg' ? failedRoot : null)
};
readinessRuntime.configureMountedResultBinder({
  adoptLegend: () => { throw new Error('forced binder failure'); }
});
const failedExpectation = readinessRuntime.registerReadinessExpectation({
  result: failedResult,
  resultIndex: 0,
  artifactIdentity: 'artifact-failed',
  generationToken: 'generation-failed',
  phase: 'generate-failed',
  isCurrent: () => true
});
const failedContext = readinessRuntime.createMountedResultContext({
  root: failedRoot,
  result: failedResult,
  resultIndex: 0
});
await assert.rejects(
  readinessRuntime.bindMountedResult(failedContext),
  /forced binder failure/
);
await assert.rejects(failedExpectation.promise, /forced binder failure/);

const selectedResult = { name: 'selected.svg', content: '<svg></svg>' };
const selectedRoot = makeSvg([]);
readinessState.results.value = [failedResult, selectedResult];
readinessRuntime.configureMountedResultBinder({
  adoptLegend: () => countStep('selection-flow-legend')
});
assert.equal(readinessRuntime.selectResult(1), true);
readinessState.svgContainer.value = {
  querySelector: (selector) => (selector === 'svg' ? selectedRoot : null)
};
const selectedContext = readinessRuntime.createMountedResultContext({
  root: selectedRoot,
  result: selectedResult,
  resultIndex: 1
});
const selectedReceipt = await readinessRuntime.bindMountedResult(selectedContext);
assert.equal(selectedReceipt.resultIndex, 1);
assert.equal(binderCalls.get('selection-flow-legend'), 1);

const metricTotal = (name, phase = null) => structuralMetrics
  .filter((metric) => metric.name === name && (phase === null || metric.phase === phase))
  .reduce((total, metric) => total + Number(metric.value || 0), 0);
assert.equal(metricTotal('previewMaterializationObservedCount', 'generate'), 1);
assert.equal(metricTotal('previewMountAdoptionCount', 'generate'), 1);
assert.equal(metricTotal('previewBinderInvocationCount', 'generate'), 1);
assert.equal(metricTotal('previewReadyReceiptEmittedCount', 'generate'), 1);
assert.equal(metricTotal('previewReadyReceiptAcceptedCount', 'generate'), 1);
assert.equal(metricTotal('previewDuplicateBindRejectedCount', 'generate'), 1);
assert.equal(metricTotal('featureDomFullScanCount', 'generate'), 0);
assert.equal(metricTotal('perFeatureListenerRegistrationCount', 'generate'), 0);
assert.equal(metricTotal('featureSearchIndexBuildCount', 'generate'), 0);
assert.ok(lifecycleEvents.some(({ name }) => name === 'preview.ready-receipt-accepted'));

console.log('preview runtime tests passed');
