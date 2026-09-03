import nodeAssert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = process.cwd();
const ref = (value) => ({ value });
let assertionCount = 0;
const caseNames = [];
const assert = {
  equal(...args) {
    assertionCount += 1;
    nodeAssert.equal(...args);
  },
  match(...args) {
    assertionCount += 1;
    nodeAssert.match(...args);
  },
  doesNotMatch(...args) {
    assertionCount += 1;
    nodeAssert.doesNotMatch(...args);
  }
};
const completeCase = (name) => caseNames.push(name);

const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = {
  onStructuralMetric(metric) {
    metrics.push({ ...metric });
  }
};

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-preview-transform-'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
await mkdir(join(tempDir, 'app', 'feature-editor'), { recursive: true });
await mkdir(join(tempDir, 'app', 'legend'), { recursive: true });
await mkdir(join(tempDir, 'services'), { recursive: true });

const copyModule = async (source, destination) => {
  await writeFile(
    join(tempDir, destination),
    await readFile(join(repoRoot, 'gbdraw', 'web', 'js', source), 'utf8'),
    'utf8'
  );
};

await copyModule('app/ui.js', 'app/ui.js');
await copyModule('app/feature-dom.js', 'app/feature-dom.js');
await copyModule('app/feature-editor/svg-actions.js', 'app/feature-editor/svg-actions.js');
await copyModule('services/runtime-test-hooks.js', 'services/runtime-test-hooks.js');
await writeFile(
  join(tempDir, 'app', 'color-utils.js'),
  'export const resolveColorToHex = (value) => value || "#94a3b8";\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'feature-utils.js'),
  [
    'export const getFeatureCaption = (feature) => feature?.label || feature?.id || "Feature";',
    'export const normalizeStringArray = (value) => Array.isArray(value) ? value : (value ? [String(value)] : []);',
    'export const resolveDisplayProteinId = () => "";'
  ].join('\n'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'pairwise-match-popup.js'),
  [
    'export const PAIRWISE_MATCH_SELECTOR = "path[data-gbdraw-pairwise-match-id]";',
    'export const buildPairwiseMatchHoverSummary = (element) => ({ id: element.getAttribute("data-gbdraw-pairwise-match-id"), title: "Pairwise match", subtitle: "", fill: "#94a3b8", rows: [] });',
    'export const buildPairwiseMatchPayload = () => null;'
  ].join('\n'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'feature-sequence-fasta.js'),
  'export const buildFeatureSequenceFastas = () => ({ nucleotideFasta: "", aminoAcidFasta: "" });\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'svg-serialization.js'),
  'export const serializeCleanSvg = () => "<svg></svg>";\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'feature-override-identity.js'),
  'export const getFeatureOverride = () => null;\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'legend', 'utils.js'),
  'export const COMPARISON_LEGEND_SELECTOR = "[data-gbdraw-role=comparison-legend]";\n',
  'utf8'
);

let nextTimerId = 1;
const timers = new Map();
let nextFrameId = 1;
const frames = new Map();
const runTimersAtDelay = (delay) => {
  [...timers.entries()]
    .filter(([, timer]) => timer.delay === delay)
    .forEach(([id, timer]) => {
      timers.delete(id);
      timer.callback();
    });
};
const flushFrames = () => {
  const pending = [...frames.entries()];
  frames.clear();
  pending.forEach(([, callback]) => callback());
};

globalThis.window = {
  innerWidth: 1440,
  innerHeight: 1000,
  setTimeout(callback, delay) {
    const id = nextTimerId++;
    timers.set(id, { callback, delay });
    return id;
  },
  clearTimeout(id) {
    timers.delete(id);
  },
  requestAnimationFrame(callback) {
    const id = nextFrameId++;
    frames.set(id, callback);
    return id;
  },
  cancelAnimationFrame(id) {
    frames.delete(id);
  },
  matchMedia: () => ({ matches: false })
};
globalThis.requestAnimationFrame = window.requestAnimationFrame;
globalThis.cancelAnimationFrame = window.cancelAnimationFrame;
globalThis.CSS = { escape: (value) => String(value) };

class ListenerTarget {
  constructor() {
    this.listeners = new Map();
    this.listenerAdds = new Map();
    this.listenerRemoves = new Map();
  }

  addEventListener(type, listener) {
    if (!this.listeners.has(type)) this.listeners.set(type, new Set());
    this.listeners.get(type).add(listener);
    this.listenerAdds.set(type, Number(this.listenerAdds.get(type) || 0) + 1);
  }

  removeEventListener(type, listener) {
    this.listeners.get(type)?.delete(listener);
    this.listenerRemoves.set(type, Number(this.listenerRemoves.get(type) || 0) + 1);
  }

  dispatch(type, event = {}) {
    for (const listener of this.listeners.get(type) || []) listener({ target: this, ...event });
  }
}

const wrapper = new ListenerTarget();
wrapper.style = {};
const canvas = { style: {} };
const uiState = {
  zoom: ref(1),
  layoutRepositionMode: ref(false),
  isPanning: ref(false),
  panStart: { x: 0, y: 0, panX: 0, panY: 0 },
  canvasPan: { x: 0, y: 0 },
  canvasContainerRef: ref(canvas),
  svgContainer: ref(wrapper)
};
const { createPanZoom } = await import(pathToFileURL(join(tempDir, 'app', 'ui.js')));
const panZoom = createPanZoom(uiState);
const interactionChanges = [];
panZoom.previewTransformInteraction.subscribe((change) => interactionChanges.push(change));

const wheelEvent = { deltaY: -100, clientX: 500, clientY: 400 };
panZoom.handleWheel(wheelEvent);
panZoom.handleWheel(wheelEvent);
panZoom.handleWheel(wheelEvent);
assert.equal(panZoom.previewTransformInteraction.isActive(), true);
assert.equal(interactionChanges.filter((change) => change.active).length, 1);
assert.equal(wrapper.listenerAdds.get('transitionend'), 1);
assert.equal(wrapper.listeners.get('transitionend').size, 1);
assert.equal([...timers.values()].filter((timer) => timer.delay === 220).length, 1);
assert.equal([...timers.values()].filter((timer) => timer.delay === 260).length, 1);
wrapper.dispatch('transitionend', { propertyName: 'opacity' });
wrapper.dispatch('transitionend', { target: canvas, propertyName: 'transform' });
assert.equal(panZoom.previewTransformInteraction.isActive(), true);
runTimersAtDelay(220);
assert.equal(panZoom.previewTransformInteraction.isActive(), true);
wrapper.dispatch('transitionend', { propertyName: 'transform' });
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
assert.equal(wrapper.listeners.get('transitionend').size, 0);
assert.equal(interactionChanges.filter((change) => !change.active).length, 1);
completeCase('repeated wheel uses one source, timer pair, and transition listener');

panZoom.handleWheel({ ...wheelEvent, deltaY: 100 });
runTimersAtDelay(260);
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
assert.equal(wrapper.listenerAdds.get('transitionend'), 2);
assert.equal(wrapper.listenerRemoves.get('transitionend'), 2);
completeCase('wheel fallback ends a clamped or missing-transition window');

const backgroundTarget = {
  tagName: 'rect',
  closest: () => null
};
panZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 100,
  clientY: 120,
  target: backgroundTarget
});
assert.equal(uiState.isPanning.value, true);
assert.equal(panZoom.previewTransformInteraction.isActive(), true);
assert.equal(canvas.style.cursor, 'grabbing');
panZoom.doPan({ clientX: 130, clientY: 150 });
flushFrames();
assert.equal(wrapper.style.transform, 'translate(30px, 30px) scale(1.2)');
panZoom.endPan({ clientX: 130, clientY: 150 });
assert.equal(uiState.isPanning.value, false);
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
assert.equal(canvas.style.cursor, 'grab');
assert.equal(wrapper.style.transition, 'transform 0.2s');
completeCase('pan updates the wrapper transform and restores the existing transition');

const changesBeforePanThenWheel = interactionChanges.length;
panZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 130,
  clientY: 150,
  target: backgroundTarget
});
panZoom.handleWheel(wheelEvent);
panZoom.endPan({ clientX: 130, clientY: 150 });
assert.equal(panZoom.previewTransformInteraction.isActive(), true);
assert.equal(
  interactionChanges.slice(changesBeforePanThenWheel).filter((change) => change.active).length,
  1
);
runTimersAtDelay(220);
wrapper.dispatch('transitionend', { propertyName: 'transform' });
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
completeCase('pan then wheel keeps wheel active after pan ends');

const changesBeforeWheelThenPan = interactionChanges.length;
panZoom.handleWheel(wheelEvent);
panZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 130,
  clientY: 150,
  target: backgroundTarget
});
runTimersAtDelay(260);
assert.equal(panZoom.previewTransformInteraction.isActive(), true);
assert.equal(
  interactionChanges.slice(changesBeforeWheelThenPan).filter((change) => !change.active).length,
  0
);
panZoom.endPan({ clientX: 130, clientY: 150 });
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
completeCase('wheel then pan keeps pan active after wheel ends');

panZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 130,
  clientY: 150,
  target: backgroundTarget
});
panZoom.endPan({ type: 'pointercancel', clientX: 130, clientY: 150 });
assert.equal(uiState.isPanning.value, false);
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
completeCase('pointercancel routes to complete pan cleanup');

panZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 130,
  clientY: 150,
  target: backgroundTarget
});
panZoom.endPan({ type: 'lostpointercapture', clientX: 130, clientY: 150 });
assert.equal(uiState.isPanning.value, false);
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
completeCase('lostpointercapture routes to complete pan cleanup');

panZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 130,
  clientY: 150,
  target: backgroundTarget
});
panZoom.handleWheel(wheelEvent);
panZoom.cancelPreviewTransformInteraction({ reconcile: false });
assert.equal(uiState.isPanning.value, false);
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
assert.equal(timers.size, 0);
assert.equal(wrapper.listeners.get('transitionend').size, 0);
completeCase('explicit canonical cancel empties every active source and resource');

panZoom.handleWheel(wheelEvent);
panZoom.resetPreviewViewport();
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
assert.equal(timers.size, 0);
completeCase('reset cancels the canonical interaction without reconciliation');

panZoom.handleWheel(wheelEvent);
panZoom.disposePanZoom();
assert.equal(panZoom.previewTransformInteraction.isActive(), false);
assert.equal(timers.size, 0);
assert.equal(wrapper.listeners.get('transitionend').size, 0);
completeCase('pan/zoom dispose removes sources, timers, listeners, and subscriptions');

class FakeStyle {
  constructor() {
    this.values = Object.create(null);
    this.writeCount = 0;
  }

  get opacity() { return this.values.opacity || ''; }
  set opacity(value) { this.values.opacity = String(value); this.writeCount += 1; }
  get filter() { return this.values.filter || ''; }
  set filter(value) { this.values.filter = String(value); this.writeCount += 1; }
  get cursor() { return this.values.cursor || ''; }
  set cursor(value) { this.values.cursor = String(value); this.writeCount += 1; }
  removeProperty(name) { delete this.values[name]; this.writeCount += 1; }
}

class FakeElement {
  constructor(kind, attributes = {}) {
    this.kind = kind;
    this.attributes = { ...attributes };
    this.id = attributes.id || '';
    this.style = new FakeStyle();
    this.attributeMutationCount = 0;
  }

  getAttribute(name) {
    return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null;
  }

  setAttribute(name, value) {
    this.attributes[name] = String(value);
    this.attributeMutationCount += 1;
  }

  removeAttribute(name) {
    if (Object.hasOwn(this.attributes, name)) {
      delete this.attributes[name];
      this.attributeMutationCount += 1;
    }
  }

  hasAttribute(name) {
    return Object.hasOwn(this.attributes, name);
  }

  matches(selector) {
    if (this.kind === 'feature') {
      return selector.includes('data-gbdraw-feature-id') || selector.includes('path[id^="f"]');
    }
    if (this.kind === 'match') return selector.includes('data-gbdraw-pairwise-match-id');
    return false;
  }

  closest(selector) {
    if (selector === 'g[id]' || selector === 'svg' || selector.includes('button')) return null;
    return this.matches(selector) ? this : null;
  }
}

class FakeSvg extends ListenerTarget {
  constructor(elements) {
    super();
    this.elements = elements;
    this.queryCounts = new Map();
  }

  contains(element) {
    return element === this || this.elements.includes(element);
  }

  querySelectorAll(selector) {
    this.queryCounts.set(selector, Number(this.queryCounts.get(selector) || 0) + 1);
    if (selector === '[data-orthogroup-id]') {
      return this.elements.filter((element) => element.hasAttribute('data-orthogroup-id'));
    }
    if (selector.includes('data-gbdraw-pairwise-match-id')) {
      return this.elements.filter((element) => element.kind === 'match');
    }
    return this.elements.filter((element) => element.kind === 'feature');
  }

  getElementById(id) {
    return this.elements.find((element) => element.id === id) || null;
  }
}

const featureOne = new FakeElement('feature', {
  id: 'f1',
  'data-gbdraw-feature-id': 'f1',
  'data-gbdraw-feature-part': 'block',
  fill: '#123456',
  'data-orthogroup-id': 'OG1'
});
const featureTwo = new FakeElement('feature', {
  id: 'f2',
  'data-gbdraw-feature-id': 'f2',
  'data-gbdraw-feature-part': 'block',
  fill: '#654321',
  'data-orthogroup-id': 'OG1'
});
const match = new FakeElement('match', {
  id: 'm1',
  'data-gbdraw-pairwise-match-id': 'm1',
  'data-orthogroup-id': 'OG1'
});
const unrelatedPath = new FakeElement('generic', {
  id: 'axis-path',
  'data-orthogroup-id': 'OG2'
});
const svg = new FakeSvg([featureOne, featureTwo, match, unrelatedPath]);
globalThis.document = {
  elementsFromPoint: () => [featureOne],
  createElement: () => { throw new Error('Hover summary should stay disabled in this Contract.'); },
  body: { appendChild: () => {} }
};

const featureInteractionWrapper = new ListenerTarget();
featureInteractionWrapper.style = {};
const featureInteractionCanvas = { style: {} };
const featureInteractionState = {
  zoom: ref(1),
  layoutRepositionMode: ref(false),
  isPanning: ref(false),
  panStart: { x: 0, y: 0, panX: 0, panY: 0 },
  canvasPan: { x: 0, y: 0 },
  canvasContainerRef: ref(featureInteractionCanvas),
  svgContainer: ref(featureInteractionWrapper)
};
const featurePanZoom = createPanZoom(featureInteractionState);

const featureMap = new Map([
  ['f1', { id: 'f1', rendered_svg_id: 'f1', label: 'One', orthogroupId: 'OG1', start: 0, end: 10 }],
  ['f2', { id: 'f2', rendered_svg_id: 'f2', label: 'Two', orthogroupId: 'OG1', start: 10, end: 20 }]
]);
let selectedFeatureId = '';
const featureSelection = {
  consumeSuppressNextClick: () => false,
  getSelectableFeatureTarget: (eventLike, root) => (
    eventLike?.target?.kind === 'feature' && root.contains(eventLike.target)
      ? eventLike.target
      : null
  ),
  toggleFeatureSelection(featureId) {
    selectedFeatureId = featureId;
  },
  markPlainFeatureClick() {}
};
const state = {
  results: ref([{ content: '<svg></svg>' }]),
  selectedResultIndex: ref(0),
  isPanning: ref(false),
  orthogroups: ref([]),
  collinearGroups: ref([]),
  orthogroupNameOverrides: {},
  orthogroupDescriptionOverrides: {},
  extractedFeatures: ref([...featureMap.values()]),
  biologicalFeatures: ref([]),
  featuresBySvgId: ref(featureMap),
  featureColorOverrides: {},
  featureVisibilityOverrides: {},
  svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svg : null }),
  clickedFeature: ref(null),
  clickedFeaturePos: { x: 0, y: 0 },
  clickedPairwiseMatch: ref(null),
  clickedPairwiseMatchPos: { x: 0, y: 0 },
  matchSequenceRegistry: null,
  selectedAnnotation: ref(null),
  featurePopupSize: { width: 0, height: 0 },
  featureSelectionDrag: { active: false },
  skipCaptureBaseConfig: ref(false),
  adv: { rich_feature_popup: false }
};
const { createFeatureSvgActions } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'svg-actions.js'))
);
const featureActions = createFeatureSvgActions({
  state,
  getFeatureColor: () => '#123456',
  getEffectiveLegendCaption: () => '',
  featureSelection,
  previewTransformInteraction: featurePanZoom.previewTransformInteraction
});
assert.equal(featureActions.attachSvgFeatureHandlers({ root: svg }), true);
assert.equal(featureActions.attachSvgFeatureHandlers({ root: svg }), false);
for (const type of ['mouseover', 'mousemove', 'mouseout', 'click', 'keydown']) {
  assert.equal(svg.listeners.get(type).size, 1, type);
}
completeCase('delegated Feature listeners install once per SVG root');

const eventAtFeature = {
  target: featureOne,
  relatedTarget: null,
  clientX: 500,
  clientY: 400,
  stopPropagation() {},
  preventDefault() {}
};
featurePanZoom.startPan({
  button: 0,
  shiftKey: false,
  clientX: 500,
  clientY: 400,
  target: backgroundTarget
});
assert.equal([...svg.queryCounts.values()].reduce((total, count) => total + count, 0), 0);
assert.equal(
  [featureOne, featureTwo, match, unrelatedPath]
    .reduce((total, element) => total + element.style.writeCount, 0),
  0
);
featurePanZoom.cancelPreviewTransformInteraction({ reconcile: false });
completeCase('absent-hover gesture start performs no scan or presentation mutation');

assert.equal(featureActions.preparePairwiseInteractionAffordances({ root: svg }), true);
assert.equal(match.style.cursor, 'pointer');
assert.equal(featureOne.style.cursor, '');
completeCase('pairwise keyboard affordance remains separate from Feature cursor ownership');

svg.dispatch('mouseover', eventAtFeature);
assert.equal(metrics.filter((metric) => metric.name === 'featureDomFullScanCount').length, 1);
assert.equal(featureOne.style.opacity, '0.7');
assert.equal(featureTwo.style.opacity, '0.7');
assert.equal(match.style.opacity, '0.7');
assert.equal(featureOne.style.cursor, '');
assert.equal(featureTwo.style.cursor, '');
assert.equal(unrelatedPath.style.opacity, '');
completeCase('ordinary Feature hover builds one cold index and tracks only its fan-out');

const queryCountBeforeCleanup = [...svg.queryCounts.values()]
  .reduce((total, count) => total + count, 0);
featurePanZoom.startPan({ ...eventAtFeature, button: 0, shiftKey: false, target: backgroundTarget });
featurePanZoom.handleWheel(wheelEvent);
assert.equal(metrics.filter((metric) => metric.name === 'previewTransformHoverCleanupCount').length, 2);
const styleWritesAfterCleanup = [featureOne, featureTwo, match]
  .reduce((total, element) => total + element.style.writeCount, 0);
const attributeWritesAfterCleanup = [featureOne, featureTwo, match]
  .reduce((total, element) => total + element.attributeMutationCount, 0);
const queryCountAfterCleanup = [...svg.queryCounts.values()]
  .reduce((total, count) => total + count, 0);
assert.equal(queryCountAfterCleanup, queryCountBeforeCleanup);
assert.equal(unrelatedPath.style.writeCount, 0);
completeCase('existing-hover gesture clear uses tracked elements with zero scan delta');

for (let index = 0; index < 20; index += 1) {
  svg.dispatch('mouseover', { ...eventAtFeature, target: index % 2 ? featureOne : featureTwo });
  svg.dispatch('mouseout', { ...eventAtFeature, target: featureOne, relatedTarget: featureTwo });
  svg.dispatch('mousemove', { ...eventAtFeature, target: featureTwo });
}
assert.equal(
  [featureOne, featureTwo, match].reduce((total, element) => total + element.style.writeCount, 0),
  styleWritesAfterCleanup
);
assert.equal(
  [featureOne, featureTwo, match].reduce((total, element) => total + element.attributeMutationCount, 0),
  attributeWritesAfterCleanup
);
assert.equal(
  [...svg.queryCounts.values()].reduce((total, count) => total + count, 0),
  queryCountAfterCleanup
);
assert.equal(metrics.filter((metric) => metric.name === 'previewTransformHoverCleanupCount').length, 2);
assert.equal(unrelatedPath.style.writeCount, 0);
completeCase('active delegated hover events cause zero scan, style, and attribute deltas');

featurePanZoom.endPan(eventAtFeature);
assert.equal(featurePanZoom.previewTransformInteraction.isActive(), true);
runTimersAtDelay(220);
featureInteractionWrapper.dispatch('transitionend', { propertyName: 'transform' });
assert.equal(featurePanZoom.previewTransformInteraction.isActive(), false);
flushFrames();
assert.equal(featureOne.style.opacity, '0.7');
assert.equal(metrics.filter((metric) => metric.name === 'previewTransformHoverReconcileCount').length, 1);
assert.equal(metrics.filter((metric) => metric.name === 'featureDomFullScanCount').length, 1);
svg.dispatch('mouseout', eventAtFeature);
assert.equal(featureOne.style.opacity, '');
completeCase('ordinary Feature hover resumes after the final source settles');

svg.dispatch('click', eventAtFeature);
assert.equal(state.clickedFeature.value?.id, 'f1');
svg.dispatch('click', { ...eventAtFeature, ctrlKey: true });
assert.equal(selectedFeatureId, 'f1');
completeCase('Feature click and modifier selection remain available after settle');

featurePanZoom.startPan({ ...eventAtFeature, button: 0, shiftKey: false, target: backgroundTarget });
featurePanZoom.endPan(eventAtFeature);
assert.equal(frames.size, 1);
const replacementFeature = new FakeElement('feature', {
  id: 'f1',
  'data-gbdraw-feature-id': 'f1',
  'data-gbdraw-feature-part': 'block'
});
const replacementSvg = new FakeSvg([replacementFeature]);
assert.equal(featureActions.attachSvgFeatureHandlers({ root: replacementSvg }), true);
assert.equal(frames.size, 0);
for (const listeners of svg.listeners.values()) assert.equal(listeners.size, 0);
for (const type of ['mouseover', 'mousemove', 'mouseout', 'click', 'keydown']) {
  assert.equal(replacementSvg.listeners.get(type).size, 1, type);
}
completeCase('SVG root replacement cancels reconciliation and removes old listeners');

featureActions.dispose();
for (const listeners of replacementSvg.listeners.values()) assert.equal(listeners.size, 0);
const cleanupCountBeforeDisposedInteraction = metrics
  .filter((metric) => metric.name === 'previewTransformHoverCleanupCount')
  .length;
featurePanZoom.startPan({ ...eventAtFeature, button: 0, shiftKey: false, target: backgroundTarget });
featurePanZoom.cancelPreviewTransformInteraction({ reconcile: false });
assert.equal(
  metrics.filter((metric) => metric.name === 'previewTransformHoverCleanupCount').length,
  cleanupCountBeforeDisposedInteraction
);
featurePanZoom.disposePanZoom();
completeCase('Feature action dispose removes delegated listeners and lifecycle subscription');

const featureDomSource = await readFile(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-dom.js'),
  'utf8'
);
assert.doesNotMatch(featureDomSource, /markCursor|style\.cursor\s*=\s*['"]pointer['"]/);
const indexHtml = await readFile(join(repoRoot, 'gbdraw', 'web', 'index.html'), 'utf8');
const cursorRule = indexHtml.match(
  /\.gbdraw-preview-surface\s+:is\(([\s\S]*?)\)\s*\{\s*cursor:\s*pointer;\s*\}/
);
assert.equal(Boolean(cursorRule), true);
const cssFeatureSelector = cursorRule[1]
  .split(',')
  .map((selector) => selector.trim())
  .join(', ');
const { FEATURE_SELECTOR } = await import(pathToFileURL(join(tempDir, 'app', 'feature-dom.js')));
assert.equal(cssFeatureSelector, FEATURE_SELECTOR);
for (const excludedSelector of [
  'data-gbdraw-pairwise-match-id',
  'data-match-kind',
  'data-orthogroup-id',
  'legend',
  'text',
  'circle',
  'line'
]) {
  assert.equal(cursorRule[1].includes(excludedSelector), false, excludedSelector);
}
completeCase('declarative Feature cursor exactly matches the production Feature selector');

const svgActionsSource = await readFile(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'feature-editor', 'svg-actions.js'),
  'utf8'
);
const trackedClearSource = svgActionsSource.match(
  /const clearTrackedHoverStyles = \(\) => \{[\s\S]*?const clearActiveMatchHover = \(\) => \{[\s\S]*?\n\s*\};/
)?.[0] || '';
assert.match(trackedClearSource, /activeHoverElements/);
assert.doesNotMatch(
  trackedClearSource,
  /ensureFeaturePaths|ensureFeatureLookup|ensureFeatureOrthogroupIndex|ensureComparisonIndexes|querySelectorAll|cursor/
);
completeCase('gesture-start clear source has no ensure, scan, or cursor call path');

const appSetupSource = await readFile(
  join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'app-setup.js'),
  'utf8'
);
assert.match(
  appSetupSource,
  /installDelegatedInteractions\(context\)\s*\{\s*cancelPreviewTransformInteraction\(\{ reconcile: false \}\);\s*featureActions\.attachSvgFeatureHandlers/
);
completeCase('preview replacement cancels the canonical lifecycle before rebinding');

assert.match(indexHtml, /transition:\s*isPanning\s*\?\s*'none'\s*:\s*'transform 0\.2s'/);
assert.match(indexHtml, /@pointercancel="endPan"/);
assert.match(indexHtml, /@lostpointercapture="endPan"/);
completeCase('template preserves transition and routes interrupted pan cleanup');

delete globalThis.__GBDRAW_TEST_HOOKS__;
console.log(JSON.stringify({
  status: 'preview transform interaction tests passed',
  assertionCount,
  cases: caseNames
}, null, 2));
