// The Result preview's hover, click, popup, alignment and search paths do no
// work the user does not see (band-hover-perf-20261015, PR-A): each case pins
// one removed waste with a counter and checks that what is shown is unchanged.
import assert from 'node:assert/strict';
import { cp, mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { test } from 'node:test';
import { pathToFileURL } from 'node:url';
import { withDrawings } from './helpers/drawing-state.mjs';

const repoRoot = process.cwd();
const ref = (value) => ({ value });
const counts = {
  hoverSummaryBuilds: 0,
  renderedIdentity: 0,
  matchMedia: 0
};
globalThis.__PREVIEW_WASTE_COUNTS__ = counts;

const tempDir = await mkdtemp(join(tmpdir(), 'gbdraw-preview-waste-'));
// The slice editors and svg-serialization resolve a drawing's ID by its mode.
await cp(join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'drawing-view.js'), join(tempDir, 'services', 'drawing-view.js'));
await writeFile(join(tempDir, 'package.json'), '{"type":"module"}\n', 'utf8');
for (const folder of ['app/feature-editor', 'services/drawing', 'utils']) {
  await mkdir(join(tempDir, folder), { recursive: true });
}
const copyModule = async (source, destination = source) => writeFile(
  join(tempDir, destination),
  await readFile(join(repoRoot, 'gbdraw', 'web', 'js', source), 'utf8'),
  'utf8'
);
for (const source of [
  'utils/error-normalization.js', 'app/ui.js', 'services/feature-dom.js',
  'app/feature-editor/svg-actions.js', 'utils/popup-bounds.js', 'services/runtime-test-hooks.js',
  'services/feature-placement.js', 'services/json-clone.js', 'app/feature-selection.js',
  'services/svg-serialization.js', 'services/settings-sections.js', 'services/drawn-groups.js',
  'services/drawing/feature-edits.js', 'services/drawing/paint.js', 'services/drawing/slice-handle.js', 'utils/reactivity.js'
]) await copyModule(source);
// The identity owner, with its rendered-identity reads counted.
await copyModule('services/feature-identity.js', 'services/feature-identity-owner.js');
const stubs = {
  'services/feature-identity.js': [
    "import * as owner from './feature-identity-owner.js';",
    "export * from './feature-identity-owner.js';",
    'export const renderedFeatureIdentity = (feature) => {',
    '  globalThis.__PREVIEW_WASTE_COUNTS__.renderedIdentity += 1;',
    '  return owner.renderedFeatureIdentity(feature);',
    '};'
  ],
  'services/feature-catalog.js': [
    'export const resultCatalogFeatures = () => null;',
    'export const stableFeatureOverrideKey = () => "";'
  ],
  'utils/color-utils.js': ['export const resolveColorToHex = (value) => value || "#94a3b8";'],
  'services/feature-utils.js': [
    'export const formatFeatureLength = (feature) => `${feature.end - feature.start} bp`;',
    'export const formatFeatureLocation = (feature) => `${feature.start + 1}..${feature.end}`;',
    'export const getFeatureCaption = (feature) => feature?.label || feature?.id || "Feature";',
    'export const normalizeStringArray = (value) => Array.isArray(value) ? value : (value ? [String(value)] : []);',
    'export const resolveDisplayProteinId = () => "";',
    'export const resolveFeatureProteinId = () => "";'
  ],
  // The hover summary of a band: its builder is counted.
  'app/pairwise-match-popup.js': [
    'export const PAIRWISE_MATCH_SELECTOR = "path[data-gbdraw-pairwise-match-id]";',
    'export const buildPairwiseMatchHoverSummary = (element) => {',
    '  globalThis.__PREVIEW_WASTE_COUNTS__.hoverSummaryBuilds += 1;',
    '  const id = element.getAttribute("data-gbdraw-pairwise-match-id");',
    '  return { id, title: "Pairwise match", subtitle: id, fill: "#94a3b8", rows: [',
    '    { label: "Kind", value: "pairwise" },',
    '    { label: "Identity", value: element.getAttribute("data-identity") },',
    '    { label: "Query", value: element.getAttribute("data-qstart") + ".." + element.getAttribute("data-qend") },',
    '    { label: "Subject", value: "" }',
    '  ] };',
    '};',
    'export const buildMatchPopupPayload = () => null;'
  ],
  'app/rule-matching.js': [
    'export const reportRuleRunFailure = (_state, _operation, run) => run();',
    'export const runWhenPrepared = (_state, _preparations, commit) => commit();'
  ],
  'services/feature-sequence-fasta.js': ['export const buildFeatureSequenceFastas = () => ({ nucleotideFasta: "", aminoAcidFasta: "" });'],
  'services/feature-override-identity.js': ['export const getFeatureOverride = () => null;'],
  'services/legend-svg.js': [
    'export const COMPARISON_LEGEND_SELECTOR = "[data-gbdraw-role=comparison-legend]";',
    'export const legendEntryKey = (entry) => String(entry?.caption || "");'
  ],
  'app/preview-item-popup.js': [
    'export const PREVIEW_ITEM_SELECTOR = "g[data-legend-key]";',
    'export const previewItemTarget = () => null;'
  ]
};
for (const [path, lines] of Object.entries(stubs)) {
  await writeFile(join(tempDir, path), `${lines.join('\n')}\n`, 'utf8');
}

let nextTimerId = 1;
const timers = new Map();
let nextFrameId = 1;
const frames = new Map();
const runTimersAtDelay = (delay) => {
  [...timers.entries()].filter(([, timer]) => timer.delay === delay).forEach(([id, timer]) => {
    timers.delete(id);
    timer.callback();
  });
};
const flushFrames = () => {
  const pending = [...frames.values()];
  frames.clear();
  pending.forEach((callback) => callback());
};

globalThis.window = {
  innerWidth: 1440,
  innerHeight: 1000,
  Vue: { computed: (getter) => ({ get value() { return getter(); } }), reactive: (value) => value },
  setTimeout(callback, delay) {
    const id = nextTimerId++;
    timers.set(id, { callback, delay });
    return id;
  },
  clearTimeout(id) { timers.delete(id); },
  requestAnimationFrame(callback) {
    const id = nextFrameId++;
    frames.set(id, callback);
    return id;
  },
  cancelAnimationFrame(id) { frames.delete(id); },
  matchMedia() {
    counts.matchMedia += 1;
    return { matches: true };
  },
  addEventListener() {},
  removeEventListener() {}
};
globalThis.requestAnimationFrame = window.requestAnimationFrame;
globalThis.cancelAnimationFrame = window.cancelAnimationFrame;
globalThis.CSS = { escape: (value) => String(value) };
const mutationObservers = [];
globalThis.MutationObserver = class {
  constructor(callback) {
    this.callback = callback;
    mutationObservers.push(this);
  }

  observe() {}

  disconnect() {}
};

// The HTML the tooltip and the alignment overlay build, counting the DOM
// changes that make the browser restyle (node insertions and class changes).
const domLog = { childList: 0, classChanges: 0 };
class FakeText {
  constructor(data) {
    this.data = String(data);
    this.parentNode = null;
  }

  get textContent() { return this.data; }
}

const classListOf = (node) => ({
  toggle(token, force) {
    const tokens = new Set(String(node.className || '').split(/\s+/).filter(Boolean));
    const has = tokens.has(token);
    const wanted = force === undefined ? !has : Boolean(force);
    if (wanted === has) return wanted;
    if (wanted) tokens.add(token);
    else tokens.delete(token);
    node.className = [...tokens].join(' ');
    domLog.classChanges += 1;
    return wanted;
  },
  add(token) { this.toggle(token, true); },
  remove(token) { this.toggle(token, false); },
  contains(token) { return String(node.className || '').split(/\s+/).includes(token); }
});

class FakeNode {
  constructor(tagName) {
    this.tagName = tagName;
    this.className = '';
    this.children = [];
    this.parentNode = null;
    this.hidden = false;
    this.style = {};
    this.attributes = {};
    this.dataset = {};
    this.rectCalls = 0;
    this.listeners = new Map();
    this.innerHTML = '';
    this.offsetWidth = 40;
  }

  get classList() { return classListOf(this); }

  get parentElement() { return this.parentNode; }

  get isConnected() {
    for (let node = this; node; node = node.parentNode) if (node === fakeBody || node.connectedRoot) return true;
    return false;
  }

  get textContent() { return this.children.map((child) => child.textContent).join(''); }

  set textContent(value) {
    this.children = [];
    this.appendChild(new FakeText(value));
  }

  appendChild(child) {
    child.parentNode = this;
    this.children.push(child);
    if (this.isConnected) domLog.childList += 1;
    return child;
  }

  append(...children) {
    children.forEach((child) => this.appendChild(typeof child === 'string' ? new FakeText(child) : child));
  }

  replaceChildren(...children) {
    this.children = [];
    if (this.isConnected) domLog.childList += 1;
    this.append(...children);
  }

  remove() {
    if (!this.parentNode) return;
    this.parentNode.children = this.parentNode.children.filter((child) => child !== this);
    this.parentNode = null;
  }

  contains(node) {
    for (let current = node; current; current = current.parentNode) if (current === this) return true;
    return false;
  }

  setAttribute(name, value) { this.attributes[name] = String(value); }

  getAttribute(name) { return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null; }

  addEventListener(type, listener) { this.listeners.set(type, listener); }

  removeEventListener(type) { this.listeners.delete(type); }

  focus() {}

  getBoundingClientRect() {
    this.rectCalls += 1;
    return { left: 0, top: 0, right: 300, bottom: 120, width: 300, height: 120 };
  }
}
const fakeBody = new FakeNode('body');
globalThis.HTMLElement = FakeNode;
globalThis.document = {
  body: fakeBody,
  activeElement: null,
  elementsFromPoint: () => [],
  querySelector: () => null,
  createElement: (tagName) => new FakeNode(tagName),
  createElementNS: (_namespace, tagName) => new FakeNode(tagName),
  createTextNode: (data) => new FakeText(data),
  addEventListener() {},
  removeEventListener() {}
};

// The Result SVG: features and bands answer the selectors the handlers use.
class FakeElement {
  constructor(kind, attributes = {}) {
    this.kind = kind;
    this.attributes = { ...attributes };
    this.id = attributes.id || '';
    this.className = '';
    this.style = {};
  }

  get classList() { return classListOf(this); }

  getAttribute(name) { return Object.hasOwn(this.attributes, name) ? this.attributes[name] : null; }

  setAttribute(name, value) { this.attributes[name] = String(value); }

  removeAttribute(name) { delete this.attributes[name]; }

  hasAttribute(name) { return Object.hasOwn(this.attributes, name); }

  matches(selector) {
    if (this.kind === 'feature') return selector.includes('data-gbdraw-feature-id') || selector.includes('path[id^="f"]');
    if (this.kind === 'match') return selector.includes('data-gbdraw-pairwise-match-id');
    return false;
  }

  closest(selector) {
    if (selector === 'g[id]' || selector === 'svg' || selector.includes('button')) return null;
    return this.matches(selector) ? this : null;
  }
}

class FakeSvg {
  constructor(elements) {
    this.elements = elements;
    this.listeners = new Map();
    this.queryCount = 0;
    this.isConnected = true;
  }

  addEventListener(type, listener) { this.listeners.set(type, listener); }

  removeEventListener(type) { this.listeners.delete(type); }

  dispatch(type, event) { this.listeners.get(type)?.({ stopPropagation() {}, preventDefault() {}, ...event }); }

  contains(element) { return element === this || this.elements.includes(element); }

  querySelectorAll(selector) {
    this.queryCount += 1;
    if (selector === '[data-orthogroup-id]') return this.elements.filter((element) => element.hasAttribute('data-orthogroup-id'));
    if (selector.includes('data-gbdraw-pairwise-match-id')) return this.elements.filter((element) => element.kind === 'match');
    return this.elements.filter((element) => element.kind === 'feature');
  }

  getElementById(id) { return this.elements.find((element) => element.id === id) || null; }
}

const { createFeatureSvgActions, HOVER_HIGHLIGHT_DWELL_MS } = await import(pathToFileURL(join(tempDir, 'app', 'feature-editor', 'svg-actions.js')));
const { createPanZoom } = await import(pathToFileURL(join(tempDir, 'app', 'ui.js')));

const featureRow = (index, extra = {}) => ({
  id: `f${index}`, rendered_svg_id: `f${index}`, record_index: 0, feature_index: index, type: 'CDS',
  start: index * 100, end: index * 100 + 90, record_id: 'REC', ...extra
});
const band = (id, attributes = {}) => new FakeElement('match', {
  id, 'data-gbdraw-pairwise-match-id': id, 'data-identity': '91.5', 'data-qstart': '1', 'data-qend': '90', ...attributes
});

// One preview: an SVG with features and bands, the SVG owner, and a pan/zoom owner.
const createPreview = ({ features, bands }) => {
  const featureElements = features.map((feature) => new FakeElement('feature', {
    id: feature.rendered_svg_id, 'data-gbdraw-feature-id': feature.rendered_svg_id, 'data-gbdraw-feature-part': 'block',
    fill: '#123456', ...(feature.orthogroupId ? { 'data-orthogroup-id': feature.orthogroupId } : {})
  }));
  const svg = new FakeSvg([...featureElements, ...bands]);
  const viewport = new FakeNode('div');
  viewport.connectedRoot = true;
  const scrollRoot = new FakeNode('div');
  viewport.appendChild(scrollRoot);
  const surface = new FakeNode('div');
  scrollRoot.appendChild(surface);
  Object.assign(surface, { querySelector: (selector) => (selector === 'svg' ? svg : null), getAnimations: () => [] });
  const canvas = { style: {} };
  const panZoom = createPanZoom({
    zoom: ref(1), layoutRepositionMode: ref(false), isPanning: ref(false),
    panStart: { x: 0, y: 0, panX: 0, panY: 0 }, canvasPan: { x: 0, y: 0 },
    canvasContainerRef: ref(canvas), svgContainer: ref(surface)
  });
  const featureMap = new Map(features.map((feature) => [feature.rendered_svg_id, feature]));
  const state = withDrawings({
    orthogroups: ref([]), orthogroupNameOverrides: {}, orthogroupDescriptionOverrides: {},
    extractedFeatures: ref(features), biologicalFeatures: ref([]), featuresBySvgId: ref(featureMap),
    svgContainer: ref(surface), clickedFeature: ref(null), clickedFeaturePos: { x: 0, y: 0 },
    clickedPairwiseMatch: ref(null), clickedPairwiseMatchPos: { x: 0, y: 0 }, matchSequenceRegistry: null,
    selectedAnnotation: ref(null), featurePopupSize: { width: 0, height: 0 }, featureSelectionDrag: { active: false }
  });
  const actions = createFeatureSvgActions({
    state,
    getFeatureColor: () => '#123456',
    getEffectiveLegendCaption: (feature) => feature.legend || '',
    previewTransformInteraction: panZoom.previewTransformInteraction
  });
  assert.equal(actions.attachSvgFeatureHandlers({ root: svg }), true);
  const hover = (target) => svg.dispatch('mouseover', { target, relatedTarget: null, clientX: 400, clientY: 300 });
  const leave = (target) => svg.dispatch('mouseout', { target, relatedTarget: null, clientX: 0, clientY: 0 });
  return { svg, viewport, surface, panZoom, actions, featureElements, hover, leave };
};

const tooltip = () => fakeBody.children.find((node) => node.attributes.role === 'tooltip');
// What the tooltip shows: its heading, subtitle, and the rows it does not hide.
const shownSummary = () => {
  const element = tooltip();
  const [title, ...rows] = element.children;
  const [heading, subtitle] = title.children[1].children;
  return {
    hidden: element.hidden,
    heading: heading.textContent,
    subtitle: subtitle.textContent,
    rows: rows.filter((row) => row.style.display !== 'none').map((row) => [
      row.children[0].textContent, row.children[1].textContent, row.children[1].classList.contains('is-clamped')
    ])
  };
};

const features = [
  featureRow(0, { gene: 'abcA', locus_tag: 'L0', product: 'ABC transporter', note: 'membrane', orthogroupId: 'OG1', legend: 'Transport' }),
  featureRow(1, { locus_tag: 'L1', product: 'kinase' }),
  ...Array.from({ length: 38 }, (_, index) => featureRow(index + 2))
];
const bands = [band('m1'), band('m2', { 'data-qstart': '5', 'data-qend': '60' }), band('m3')];
const preview = createPreview({ features, bands });

test('OV-915: a band hover builds its summary only when the summary shows', () => {
  for (let round = 0; round < 3; round += 1) {
    bands.forEach((element) => {
      preview.hover(element);
      preview.leave(element);
    });
  }
  assert.equal(counts.hoverSummaryBuilds, 0, 'a passing pointer builds no summary');
  preview.hover(bands[0]);
  assert.equal(counts.hoverSummaryBuilds, 0);
  runTimersAtDelay(180);
  assert.equal(counts.hoverSummaryBuilds, 1);
  assert.deepEqual(shownSummary(), {
    hidden: false, heading: 'Pairwise match', subtitle: 'm1',
    rows: [['Kind', 'pairwise', false], ['Identity', '91.5', false], ['Query', '1..90', true]]
  });
  // While the summary shows, the next band's summary replaces it at once.
  preview.hover(bands[1]);
  assert.equal(counts.hoverSummaryBuilds, 2);
  assert.deepEqual(shownSummary().rows[2], ['Query', '5..60', true]);
  preview.leave(bands[1]);
  assert.equal(tooltip().hidden, true);
});

test('OV-933: a shown summary rewrites its text, adding no node and changing no class between like summaries', () => {
  preview.hover(preview.featureElements[0]);
  runTimersAtDelay(180);
  const featureSummary = {
    hidden: false, heading: 'CDS: abcA', subtitle: '1..90',
    rows: [
      ['Locus', 'L0', false], ['Product', 'ABC transporter', true], ['Note', 'membrane', true],
      ['Length', '90 bp', false], ['Location', '1..90', false], ['Record', 'REC', false],
      ['Similarity group', 'OG1', false], ['Legend', 'Transport', false]
    ]
  };
  assert.deepEqual(shownSummary(), featureSummary);
  const before = { ...domLog };
  preview.hover(bands[0]);
  const betweenBands = domLog.classChanges;
  preview.hover(bands[1]);
  preview.hover(bands[2]);
  assert.equal(domLog.classChanges - betweenBands, 0, 'like summaries change no class');
  assert.deepEqual(shownSummary().rows, [['Kind', 'pairwise', false], ['Identity', '91.5', false], ['Query', '1..90', true]]);
  preview.hover(preview.featureElements[1]);
  assert.deepEqual(shownSummary(), {
    hidden: false, heading: 'CDS: L1', subtitle: '101..190',
    rows: [['Product', 'kinase', true], ['Length', '90 bp', false], ['Location', '101..190', false], ['Record', 'REC', false]]
  });
  preview.hover(preview.featureElements[0]);
  assert.deepEqual(shownSummary(), featureSummary);
  assert.equal(domLog.childList - before.childList, 0, 'no node is added or removed');
  // Only a row whose clamping changes changes its class: one from the feature
  // to the bands, two to the second feature, three back to the first.
  assert.equal(domLog.classChanges - before.classChanges, 6);
  assert.equal(tooltip().textContent.includes('membrane'), true);
  preview.leave(preview.featureElements[0]);
});

test('OV-938: the hover reads the pointer query once and places the summary in the next frame', () => {
  counts.matchMedia = 0;
  for (let index = 0; index < 10; index += 1) preview.hover(preview.featureElements[index % 3]);
  preview.leave(preview.featureElements[0]);
  assert.ok(counts.matchMedia <= 1, `matchMedia ran ${counts.matchMedia} times`);
  flushFrames();
  const element = tooltip();
  preview.hover(bands[0]);
  runTimersAtDelay(180);
  const rectCallsAtShow = element.rectCalls;
  assert.equal(element.hidden, false);
  flushFrames();
  assert.equal(element.rectCalls - rectCallsAtShow, 1);
  assert.equal(rectCallsAtShow, 0, 'the show itself measures nothing');
  assert.deepEqual([element.style.left, element.style.top], ['414px', '314px']);
  preview.leave(bands[0]);
});

test('OV-937: the open match popup marks its bands from the band index', () => {
  const ribbonOne = band('og-a', { 'data-match-kind': 'orthogroup', 'data-orthogroup-id': 'OG7' });
  const ribbonTwo = band('og-b', { 'data-match-kind': 'orthogroup', 'data-orthogroup-id': 'OG8;OG7' });
  const otherRibbon = band('og-c', { 'data-match-kind': 'orthogroup', 'data-orthogroup-id': 'OG8' });
  const pairwise = band('pw-1', { 'data-orthogroup-id': 'OG7', 'data-match-kind': 'pairwise' });
  const pairwiseFragment = band('pw-1-b', { 'data-gbdraw-pairwise-match-id': 'pw-1', 'data-match-kind': 'pairwise' });
  const local = createPreview({ features: [featureRow(0)], bands: [ribbonOne, ribbonTwo, otherRibbon, pairwise, pairwiseFragment] });
  const selected = () => [ribbonOne, ribbonTwo, otherRibbon, pairwise, pairwiseFragment]
    .filter((element) => element.classList.contains('gbdraw-match-selected'))
    .map((element) => element.id);
  local.actions.markSelectedPairwiseMatch({ id: 'og-a', matchKind: 'orthogroup', orthogroupId: 'OG7' });
  assert.deepEqual(selected(), ['og-a', 'og-b']);
  const scans = local.svg.queryCount;
  local.actions.markSelectedPairwiseMatch({ id: 'pw-1', matchKind: 'pairwise', orthogroupId: 'OG7' });
  assert.deepEqual(selected(), ['pw-1', 'pw-1-b']);
  local.actions.markSelectedPairwiseMatch({ id: 'og-c', matchKind: 'orthogroup', orthogroupId: 'OG8' });
  assert.deepEqual(selected(), ['og-b', 'og-c']);
  local.actions.markSelectedPairwiseMatch(null);
  assert.deepEqual(selected(), []);
  assert.equal(local.svg.queryCount, scans, 'no band scan once the index exists');
  local.actions.dispose();
});

test('OV-937: the popup watch of the composition root leaves the bands to the SVG owner', async () => {
  const source = await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'app', 'app-setup.js'), 'utf8');
  const start = source.lastIndexOf('watch(clickedPairwiseMatch,', source.indexOf("selectedPairwiseBlockOrthogroupId.value = '';"));
  assert.ok(start > 0);
  const watchSource = source.slice(start, source.indexOf('\n  });', start));
  assert.match(watchSource, /featureActions\.markSelectedPairwiseMatch\(match\)/);
  assert.doesNotMatch(watchSource, /querySelectorAll/);
});

test('OV-935: the Review alignment finds each candidate through one identity index', () => {
  counts.renderedIdentity = 0;
  for (let index = 0; index < 5; index += 1) {
    assert.equal(preview.actions.previewAlignmentCandidate({ recordIndex: 0, featureIndex: index * 7 }), true);
    assert.equal(preview.featureElements[index * 7].style.opacity, '0.7');
    preview.actions.clearAlignmentCandidatePreview();
  }
  assert.equal(preview.actions.previewAlignmentCandidate({ recordIndex: 0, featureIndex: 99 }), false);
  assert.equal(preview.actions.previewAlignmentCandidate({ recordIndex: 0 }), false);
  assert.ok(counts.renderedIdentity <= features.length, `${counts.renderedIdentity} rendered identities for ${features.length} features`);
});

test('OV-936: the alignment overlay is not measured while the preview pans', () => {
  const local = createPreview({ features: [featureRow(0)], bands: [] });
  local.actions.showAlignmentOverlay({
    reference: { key: 'r', recordKey: 'REC', anchor: {}, name: 'Reference' },
    records: [], onSelect() {}, onHover() {}
  });
  flushFrames();
  flushFrames();
  const settled = local.viewport.rectCalls;
  assert.ok(settled > 0);
  local.panZoom.startPan({ button: 0, shiftKey: false, clientX: 100, clientY: 120, target: { tagName: 'rect', closest: () => null } });
  assert.equal(local.panZoom.previewTransformInteraction.isActive(), true);
  for (let step = 0; step < 5; step += 1) {
    mutationObservers.at(-1).callback();
    flushFrames();
  }
  assert.equal(local.viewport.rectCalls, settled, 'no overlay pass while panning');
  local.panZoom.endPan({ clientX: 130, clientY: 150 });
  flushFrames();
  assert.equal(local.viewport.rectCalls, settled + 1, 'the end of the pan places the overlay once');
  local.actions.dispose();
});

// OV-900 (D-01): the hover highlight waits until the pointer rests on a target
// for HOVER_HIGHLIGHT_DWELL_MS, so a pointer passing over bands and features
// writes no highlight; resting writes the same highlight as before.
const dwellFeatures = [
  featureRow(0, { orthogroupId: 'OG5' }), featureRow(1, { orthogroupId: 'OG5' }),
  featureRow(2, { orthogroupId: 'OG6' }), featureRow(3)
];
const ribbonFive = band('og5', { 'data-match-kind': 'orthogroup', 'data-orthogroup-id': 'OG5' });
const ribbonSix = band('og6', { 'data-match-kind': 'orthogroup', 'data-orthogroup-id': 'OG6' });
const dwellPreview = createPreview({ features: dwellFeatures, bands: [ribbonFive, ribbonSix] });
const [groupFiveA, groupFiveB, groupSix, ungrouped] = dwellPreview.featureElements;
const dwellElements = [...dwellPreview.featureElements, ribbonFive, ribbonSix];
const highlightWrites = { count: 0 };
dwellElements.forEach((element) => {
  element.style = new Proxy({}, {
    set(target, key, value) {
      if (key === 'opacity' || key === 'filter') highlightWrites.count += 1;
      target[key] = value;
      return true;
    }
  });
});
const pointerEvent = (target, relatedTarget) => ({ target, relatedTarget, clientX: 400, clientY: 300 });
const enter = (target, from = null) => dwellPreview.svg.dispatch('mouseover', pointerEvent(target, from));
const exit = (target, to = null) => dwellPreview.svg.dispatch('mouseout', pointerEvent(target, to));
// The pointer crosses the targets in order and leaves the Result, all within one dwell.
const sweep = (targets) => targets.forEach((target, index) => {
  enter(target, targets[index - 1] || null);
  exit(target, targets[index + 1] || null);
});
const lit = () => dwellElements.filter((element) => element.style.opacity === '0.7').map((element) => element.id);
const pendingDwells = () => [...timers.values()].filter((timer) => timer.delay === HOVER_HIGHLIGHT_DWELL_MS).length;
// The highlight as it was written before the dwell: inline opacity and filter, with the earlier values kept.
const highlightOf = (element) => ({
  opacity: element.style.opacity,
  filter: element.style.filter,
  savedOpacity: element.getAttribute('data-gbdraw-hover-opacity'),
  savedFilter: element.getAttribute('data-gbdraw-hover-filter')
});
const HIGHLIGHT = { opacity: '0.7', filter: 'brightness(1.2)', savedOpacity: '', savedFilter: '' };
const CLEARED = { opacity: '', filter: '', savedOpacity: null, savedFilter: null };

test('OV-900: a pointer sweep faster than the dwell writes no highlight', () => {
  const targets = [groupFiveA, ribbonSix, ungrouped, groupSix, ribbonFive, groupFiveB];
  for (let round = 0; round < 3; round += 1) sweep(targets);
  assert.equal(highlightWrites.count, 0, 'no opacity or filter is written');
  assert.deepEqual(lit(), []);
  assert.equal(pendingDwells(), 0, 'leaving the Result cancels the pending highlight');
});

test('OV-900: resting for the dwell writes the same highlight as before', () => {
  assert.equal(typeof HOVER_HIGHLIGHT_DWELL_MS, 'number');
  assert.ok(HOVER_HIGHLIGHT_DWELL_MS >= 60 && HOVER_HIGHLIGHT_DWELL_MS <= 80, `${HOVER_HIGHLIGHT_DWELL_MS} ms`);
  assert.ok(HOVER_HIGHLIGHT_DWELL_MS < 180, 'the highlight comes before the summary');
  enter(groupFiveA);
  assert.deepEqual(lit(), [], 'nothing is written while the pointer has not rested');
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  assert.deepEqual(lit(), ['f0', 'f1', 'og5'], 'the whole similarity group');
  [groupFiveA, groupFiveB, ribbonFive].forEach((element) => assert.deepEqual(highlightOf(element), HIGHLIGHT));
  exit(groupFiveA);
  assert.deepEqual(lit(), [], 'leaving clears at once');
  [groupFiveA, groupFiveB, ribbonFive].forEach((element) => assert.deepEqual(highlightOf(element), CLEARED));
  enter(ribbonSix);
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  assert.deepEqual(lit(), ['f2', 'og6']);
  exit(ribbonSix);
  enter(ungrouped);
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  assert.deepEqual(lit(), ['f3'], 'a feature without a group lights alone');
  exit(ungrouped);
  assert.deepEqual(lit(), []);
});

test('OV-900: moving from a lit target straight to another swaps the highlight once, without a gap', () => {
  enter(groupFiveA);
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  // To the group's own band: the group stays lit throughout.
  exit(groupFiveA, ribbonFive);
  enter(ribbonFive, groupFiveA);
  assert.deepEqual(lit(), ['f0', 'f1', 'og5']);
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  assert.deepEqual(lit(), ['f0', 'f1', 'og5']);
  // To another group's band: the first group stays lit until the second one replaces it.
  exit(ribbonFive, ribbonSix);
  enter(ribbonSix, ribbonFive);
  assert.deepEqual(lit(), ['f0', 'f1', 'og5']);
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  assert.deepEqual(lit(), ['f2', 'og6']);
  exit(ribbonSix);
  assert.deepEqual(lit(), []);
});

test('OV-900: the Review alignment preview cancels a pending highlight and restores the rested one', () => {
  enter(groupFiveA);
  assert.equal(dwellPreview.actions.previewAlignmentCandidate({ recordIndex: 0, featureIndex: 3 }), true);
  runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
  assert.deepEqual(lit(), ['f3'], 'only the candidate is lit');
  dwellPreview.actions.clearAlignmentCandidatePreview();
  assert.deepEqual(lit(), ['f0', 'f1', 'og5'], 'the hover under the pointer comes back at once');
  exit(groupFiveA);
  assert.deepEqual(lit(), []);
  dwellPreview.actions.dispose();
});

test('OV-913: clearing an empty selection writes no new selection and restyles nothing', async () => {
  const { createFeatureSelection, FEATURE_SELECTED_CLASS } = await import(pathToFileURL(join(tempDir, 'app', 'feature-selection.js')));
  let selectionWrites = 0;
  let selection = new Set();
  const selectedFeatureIds = {
    get value() { return selection; },
    set value(next) { selectionWrites += 1; selection = next; }
  };
  const selectedElement = new FakeElement('feature', { 'data-gbdraw-feature-id': 'f1' });
  let styledQueries = 0;
  const svg = {
    querySelectorAll(selector) {
      styledQueries += 1;
      return selector.startsWith('.') ? [] : [selectedElement];
    },
    querySelector: (selector) => (selector.split(', ').some((part) => selectedElement.classList.contains(part.slice(1))) ? selectedElement : null)
  };
  const state = {
    selectedFeatureIds, selectedFeatureAnchorId: ref(''), featureSelectionStatus: ref('x'),
    featureSelectionSuppressNextClick: ref(false), featureSelectionDrag: { active: true, committed: true },
    extractedFeatures: ref([]), featuresBySvgId: ref(new Map()), selectedFeatures: ref([]),
    svgContainer: ref({ querySelector: () => svg }), clickedFeature: ref(null), clickedPairwiseMatch: ref(null),
    clickedLabel: ref(null), mode: ref('linear')
  };
  const selectionOwner = createFeatureSelection({ state });
  selectionOwner.clearFeatureSelection({ clearStatus: true });
  selectionOwner.markPlainFeatureClick('f1');
  selectionOwner.markPlainFeatureClick('f2');
  assert.equal(selectionWrites, 0);
  assert.equal(styledQueries, 0);
  assert.equal(state.selectedFeatureAnchorId.value, 'f2', 'a plain click still sets the anchor');
  assert.equal(state.featureSelectionStatus.value, '');
  assert.deepEqual([state.featureSelectionDrag.active, state.featureSelectionDrag.committed], [false, false]);
  // A selection is cleared and its classes removed as before.
  selection = new Set(['f1']);
  selectedElement.classList.toggle(FEATURE_SELECTED_CLASS, true);
  selectionOwner.clearFeatureSelection({ clearStatus: true });
  assert.equal(selectionWrites, 1);
  assert.equal(selection.size, 0);
  assert.ok(styledQueries > 0);
  // A class left on the Result by a restore is still removed.
  const queriesBefore = styledQueries;
  selectedElement.classList.toggle(FEATURE_SELECTED_CLASS, true);
  selectionOwner.clearFeatureSelection();
  assert.ok(styledQueries > queriesBefore);
  assert.equal(selectionWrites, 1);
});
