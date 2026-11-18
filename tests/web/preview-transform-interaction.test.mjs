import nodeAssert from 'node:assert/strict';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { pathToFileURL } from 'node:url';
import { withDrawings } from './helpers/drawing-state.mjs';
import { loadPageVue } from './helpers/vue-runtime.mjs';

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
  },
  async rejects(...args) {
    assertionCount += 1;
    await nodeAssert.rejects(...args);
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
await mkdir(join(tempDir, 'services', 'drawing'), { recursive: true });
await mkdir(join(tempDir, 'utils'), { recursive: true });
await mkdir(join(tempDir, 'contracts'), { recursive: true });

const copyModule = async (source, destination) => {
  await writeFile(
    join(tempDir, destination),
    await readFile(join(repoRoot, 'gbdraw', 'web', 'js', source), 'utf8'),
    'utf8'
  );
};

await copyModule('utils/error-normalization.js', 'utils/error-normalization.js');
// error-normalization.js reads Python's figure length units from the generated mode profiles.
await copyModule('mode-profiles.generated.js', 'mode-profiles.generated.js');
await copyModule('utils/i18n.js', 'utils/i18n.js');
await copyModule('app/ui.js', 'app/ui.js');
await copyModule('services/feature-dom.js', 'services/feature-dom.js');
await copyModule('app/feature-editor/svg-actions.js', 'app/feature-editor/svg-actions.js');
await copyModule('services/feature-popup.js', 'services/feature-popup.js');
await copyModule('services/preview-affordances.js', 'services/preview-affordances.js');
await copyModule('utils/popup-bounds.js', 'utils/popup-bounds.js');
await copyModule('services/runtime-test-hooks.js', 'services/runtime-test-hooks.js');
await copyModule('services/drawn-groups.js', 'services/drawn-groups.js');
await copyModule('services/drawn-attributes.generated.js', 'services/drawn-attributes.generated.js');
await copyModule('services/feature-identity.js', 'services/feature-identity.js');
await copyModule('services/feature-placement.js', 'services/feature-placement.js');
await copyModule('services/json-clone.js', 'services/json-clone.js');
await copyModule('utils/vue-raw.js', 'utils/vue-raw.js');
// A match popup names its spans' records from the shown Result's catalog item (D5); this test opens none.
await writeFile(join(tempDir, 'services', 'catalog-sequences.js'), 'export const matchSequenceRecordsOf = () => null;\n', 'utf8');
for (const path of ['services/drawing/feature-edits.js', 'services/drawing/paint.js', 'services/drawing/slice-handle.js', 'utils/reactivity.js',
  // feature-placement.js reads the record keys of a drawing from their owner.
  'services/record-draw-selection.js', 'utils/download-names.js', 'contracts/schema.generated.js',
  // feature-identity.js reads a record source's files from the input-type table.
  'services/input-types.js',
  // svg-actions.js finds the feature a modifier click selects with feature-selection.js, which reads a Shift+click
  // range on the record axis and holds its gesture.
  'app/feature-selection.js', 'services/svg-serialization.js', 'services/region-axis.js', 'services/record-options.js',
  'utils/pointer-gesture.js',
  // svg-actions.js gives a similarity ring hit the ring hover (ring-hover.js); its Legend name and colour reads are
  // stubbed below.
  'app/feature-editor/ring-hover.js', 'services/legend-layout.js', 'services/font-families.js',
  'utils/font-family-rules.generated.js']) await copyModule(path, path);
// No admitted catalog: features open as the Result draws them.
await writeFile(
  join(tempDir, 'services', 'feature-catalog.js'),
  'export const resultCatalogFeatures = () => null;\nexport const stableFeatureOverrideKey = () => "";\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'utils', 'color-utils.js'),
  'export const resolveColorToHex = (value) => value || "#94a3b8";\n'
    + 'export const hexToRgb = () => ({ r: 148, g: 163, b: 184 });\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'feature-utils.js'),
  [
    'export const formatFeatureLength = () => "";',
    'export const formatFeatureLocation = () => "";',
    'export const getFeatureCaption = (feature) => feature?.label || feature?.id || "Feature";',
    'export const normalizeStringArray = (value) => Array.isArray(value) ? value : (value ? [String(value)] : []);',
    'export const resolveDisplayProteinId = () => "";',
    'export const resolveFeatureProteinId = () => "";',
    "export const getQualifierFirstValue = (feature, key) => { const first = (v) => (Array.isArray(v) ? v : [v]).map((x) => String(x ?? '').trim()).find(Boolean) || ''; const d = feature?.[key]; if (d !== null && d !== undefined && d !== '') return first(d); return first(feature?.qualifiers?.[key]); };"
  ].join('\n'),
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'pairwise-match-popup.js'),
  [
    'export const PAIRWISE_MATCH_SELECTOR = "path[data-gbdraw-match-id]";',
    'export const buildPairwiseMatchHoverSummary = (element) => ({ id: element.getAttribute("data-gbdraw-match-id"), title: "Pairwise match", subtitle: "", fill: "#94a3b8", rows: [] });',
    'export const buildMatchPopupPayload = () => null;'
  ].join('\n'),
  'utf8'
);
await writeFile(
  join(tempDir, 'app', 'rule-matching.js'),
  'export const reportRuleRunFailure = (_state, _operation, run) => run();\nexport const runWhenPrepared = (_state, _preparations, commit) => commit();\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'feature-sequence-fasta.js'),
  'export const buildFeatureSequenceFastas = () => ({ nucleotideFasta: "", aminoAcidFasta: "" });\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'feature-override-identity.js'),
  'export const getFeatureOverride = () => null;\n',
  'utf8'
);
await writeFile(
  join(tempDir, 'services', 'legend-svg.js'),
  'export const COMPARISON_LEGEND_SELECTOR = "[data-gbdraw-role=comparison-legend]";\n' +
    'export const legendEntryKey = (entry) => String(entry?.caption || "");\n' +
    'export const drawnLegendRowName = () => null;\n',
  'utf8'
);
// No Legend row or record definition in these fixtures opens a popup.
await writeFile(
  join(tempDir, 'services', 'match-frame-popup.js'),
  'export const PREVIEW_ITEM_SELECTOR = "g[data-legend-key]";\nexport const previewItemTarget = () => null;\n'
    + "export const previewItemPopupOf = (opened) => (opened && 'previewItem' in opened ? opened : null);\nexport const matchPopupOf = (opened) => (opened && !('previewItem' in opened) ? opened : null);\n",
  'utf8'
);
// No drawn element in these fixtures opens a settings section.
await writeFile(
  join(tempDir, 'services', 'settings-sections.js'),
  'export const DRAWN_ELEMENT_SELECTOR = "g[data-gbdraw-slot-renderer]";\nexport const drawnElementSection = () => null;\n',
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

const Vue = loadPageVue();
globalThis.window = {
  Vue,
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
// The fake roots below have no preview host element, so the ring hover binds nothing.
globalThis.HTMLElement = class {};

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
// The view is reactive, as in state.js, and the surface shows it through the
// same projection `v-preview-transform` binds: the pan and zoom owner writes the view only.
const uiState = {
  zoom: Vue.ref(1),
  layoutRepositionMode: ref(false),
  isPanning: Vue.ref(false),
  panStart: { x: 0, y: 0, panX: 0, panY: 0 },
  canvasPan: Vue.reactive({ x: 0, y: 0 }),
  previewHiddenUntilShown: Vue.ref(false),
  previewTransitionHeld: Vue.ref(false),
  canvasContainerRef: ref(canvas),
  svgContainer: ref(wrapper)
};
const { createPanZoom, projectPreviewSurface } = await import(pathToFileURL(join(tempDir, 'app', 'ui.js')));
projectPreviewSurface(wrapper, uiState);
const panZoom = createPanZoom(uiState);
const interactionChanges = [];
panZoom.subscribePreviewTransformChanges((change) => interactionChanges.push(change));

const wheelEvent = { deltaY: -100, clientX: 500, clientY: 400 };
panZoom.handleWheel(wheelEvent);
panZoom.handleWheel(wheelEvent);
panZoom.handleWheel(wheelEvent);
assert.equal(panZoom.isPreviewTransformActive(), true);
assert.equal(interactionChanges.filter((change) => change.active).length, 1);
assert.equal(wrapper.listenerAdds.get('transitionend'), 1);
assert.equal(wrapper.listeners.get('transitionend').size, 1);
assert.equal([...timers.values()].filter((timer) => timer.delay === 220).length, 1);
assert.equal([...timers.values()].filter((timer) => timer.delay === 260).length, 1);
wrapper.dispatch('transitionend', { propertyName: 'opacity' });
wrapper.dispatch('transitionend', { target: canvas, propertyName: 'transform' });
assert.equal(panZoom.isPreviewTransformActive(), true);
runTimersAtDelay(220);
assert.equal(panZoom.isPreviewTransformActive(), true);
wrapper.dispatch('transitionend', { propertyName: 'transform' });
assert.equal(panZoom.isPreviewTransformActive(), false);
assert.equal(wrapper.listeners.get('transitionend').size, 0);
assert.equal(interactionChanges.filter((change) => !change.active).length, 1);
completeCase('repeated wheel uses one source, timer pair, and transition listener');

panZoom.handleWheel({ ...wheelEvent, deltaY: 100 });
runTimersAtDelay(260);
assert.equal(panZoom.isPreviewTransformActive(), false);
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
assert.equal(panZoom.isPreviewTransformActive(), true);
assert.equal(canvas.style.cursor, 'grabbing');
panZoom.doPan({ clientX: 130, clientY: 150 });
flushFrames();
// The surface shows the view at Vue's flush, before the frame is drawn.
await Vue.nextTick();
assert.equal(wrapper.style.transform, 'translate(30px, 30px) scale(1.2)');
// The pan's live position is the view's (OV-2722), shown at once as a promoted layer.
assert.equal(`${uiState.canvasPan.x},${uiState.canvasPan.y}`, '30,30');
assert.equal(`${wrapper.style.transition}|${wrapper.style.willChange}`, 'none|transform');
panZoom.endPan({ clientX: 130, clientY: 150 });
assert.equal(uiState.isPanning.value, false);
assert.equal(panZoom.isPreviewTransformActive(), false);
assert.equal(canvas.style.cursor, 'grab');
await Vue.nextTick();
assert.equal(wrapper.style.transition, 'transform 0.2s');
completeCase('pan updates the wrapper transform and restores the existing transition');

const legendTarget = {
  tagName: 'text',
  closest: (selector) => selector === 'g[id]' ? { id: 'feature_legend' }
    : selector === '.gbdraw-preview-layout-target' ? { id: 'legend' } : null
};
panZoom.startPan({ button: 0, shiftKey: false, clientX: 100, clientY: 120, target: legendTarget });
assert.equal(uiState.isPanning.value, true, 'a bound legend must pan while Layout edit is OFF');
panZoom.endPan({ clientX: 100, clientY: 120 });
panZoom.startPan({ button: 0, shiftKey: false, clientX: 100, clientY: 120,
  target: { tagName: 'path', closest: (selector) => selector === 'g[id]' ? { id: 'f1' } : null } });
assert.equal(uiState.isPanning.value, false, 'unbound feature groups retain gesture priority');
completeCase('eligible legend bypasses the feature-prefix pan exclusion without changing feature priority');


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
assert.equal(panZoom.isPreviewTransformActive(), true);
assert.equal(
  interactionChanges.slice(changesBeforePanThenWheel).filter((change) => change.active).length,
  1
);
runTimersAtDelay(220);
wrapper.dispatch('transitionend', { propertyName: 'transform' });
assert.equal(panZoom.isPreviewTransformActive(), false);
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
assert.equal(panZoom.isPreviewTransformActive(), true);
assert.equal(
  interactionChanges.slice(changesBeforeWheelThenPan).filter((change) => !change.active).length,
  0
);
panZoom.endPan({ clientX: 130, clientY: 150 });
assert.equal(panZoom.isPreviewTransformActive(), false);
completeCase('wheel then pan keeps pan active after wheel ends');

for (const type of ['pointercancel', 'lostpointercapture']) {
  panZoom.resetPreviewViewport();
  let defaultPrevented = false;
  panZoom.startPan({
    button: 0, shiftKey: false, clientX: 130, clientY: 150,
    target: backgroundTarget,
    preventDefault() { defaultPrevented = true; }
  });
  panZoom.doPan({ clientX: 160, clientY: 190 });
  // Native text dragging can cancel a mouse pointer with zero coordinates.
  // Flush the pending real move even if its animation frame has not run.
  panZoom.endPan({ type, clientX: 0, clientY: 0 });
  assert.equal(uiState.canvasPan.x, 30);
  assert.equal(uiState.canvasPan.y, 40);
  assert.equal(defaultPrevented, true);
  assert.equal(uiState.isPanning.value, false);
  assert.equal(panZoom.isPreviewTransformActive(), false);
  completeCase(`${type} preserves the last real pointer and completes cleanup`);
}

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
assert.equal(panZoom.isPreviewTransformActive(), false);
assert.equal(timers.size, 0);
assert.equal(wrapper.listeners.get('transitionend').size, 0);
completeCase('explicit canonical cancel empties every active source and resource');

panZoom.resetPreviewViewport();
panZoom.startPan({ button: 0, shiftKey: false, clientX: 130, clientY: 150, target: backgroundTarget });
panZoom.doPan({ clientX: 160, clientY: 190 });
assert.equal(frames.size, 1);
assert.equal(canvas.style.cursor, 'grabbing');
panZoom.cancelPreviewTransformInteraction({ reconcile: false });
assert.equal(frames.size, 0, 'the cancel stops the pending pan frame');
assert.equal(canvas.style.cursor, 'grab');
assert.equal(uiState.isPanning.value, false);
panZoom.endPan({ type: 'pointercancel', clientX: 0, clientY: 0 });
assert.equal(`${uiState.canvasPan.x},${uiState.canvasPan.y}`, '0,0', 'a cancelled pan moves the figure no more');
completeCase('a cancel during a pan stops its frame and resets the cursor');

panZoom.handleWheel(wheelEvent);
panZoom.resetPreviewViewport();
assert.equal(panZoom.isPreviewTransformActive(), false);
assert.equal(timers.size, 0);
completeCase('reset cancels the canonical interaction without reconciliation');

// UI-11 Fit: a fake layout whose surface (and SVG) has its top-center origin at
// (606.5, 58) in a 896x386 frame at (100, 50), scaled and panned like the CSS.
// The stylesheet's Editor reserve resolves through a probe element.
const layout = { width: 997, height: 817 };
let editorReserve = 0;
/** The box of the figure under a pan and zoom; `measured` rounds it to single precision, as Chromium measures. */
const boxUnder = ({ x, y, zoom }, measured = false) => {
  const round = measured ? Math.fround : (value) => value;
  const width = zoom * layout.width;
  const height = zoom * layout.height;
  const left = 606.5 + x - width / 2;
  const top = 58 + y;
  return { left: round(left), top: round(top), width: round(width), height: round(height),
    right: round(left + width), bottom: round(top + height) };
};
// The view the figure shows once Vue has projected it.
const visualBox = () => boxUnder({ x: uiState.canvasPan.x, y: uiState.canvasPan.y, zoom: uiState.zoom.value });
// What a measure reads: the transform the surface's style holds now.
const measuredBox = () => {
  const [, x, y, zoom] = wrapper.style.transform.match(/^translate\((.+)px, (.+)px\) scale\((.+)\)$/).map(Number);
  return boxUnder({ x, y, zoom }, true);
};
const fakeSvg = { getBoundingClientRect: measuredBox };
globalThis.getComputedStyle = (element) => ({ transform: element.style.transform });
let probes = 0;
// The canvas's computed padding: where it lays out its content (p-2, 0.5rem).
let canvasPaddingLeft = '8px';
const computedStyle = (element) => ({
  transform: element.style?.transform, paddingLeft: element === canvas ? canvasPaddingLeft : '0px'
});
globalThis.getComputedStyle = computedStyle;
// OV-3588: the control over the canvas edge (the Editor toggle) beside the canvas in its workspace,
// measured through a copy of it placed beside it.
/** @type {{ left: number, top: number, width: number, height: number } | null} */
let edgeControl = null;
const controlCopies = [];
const control = {
  style: { transition: 'transform 0.3s' },
  cloneNode: (deep) => {
    const copy = {
      deep, placed: false, style: {},
      getBoundingClientRect: () => {
        assert.equal(copy.placed, true, 'the copy is measured in place');
        return { ...edgeControl, right: edgeControl.left + edgeControl.width, bottom: edgeControl.top + edgeControl.height };
      },
      remove: () => { copy.placed = false; }
    };
    controlCopies.push(copy);
    return copy;
  },
  after: (copy) => { copy.placed = true; }
};
const workspace = {
  querySelector: (selector) => (selector === ':scope > [data-preview-edge-control]' ? control : null)
};
Object.assign(wrapper, { querySelector: (selector) => (selector === 'svg' ? fakeSvg : null), getBoundingClientRect: measuredBox });
Object.defineProperty(canvas, 'parentElement', { get: () => (edgeControl ? workspace : null), configurable: true });
// The canvas's classic scrollbars (8 px, index.html) narrow its client box while the shown figure overflows it.
let scrollbars = false;
const clientInset = () => {
  if (!scrollbars) return 0;
  const box = measuredBox();
  return box.left < 100 || box.right > 996 || box.top < 50 || box.bottom > 436 ? 8 : 0;
};
Object.defineProperty(canvas, 'clientWidth', { get: () => 896 - clientInset(), configurable: true });
Object.defineProperty(canvas, 'clientHeight', { get: () => 386 - clientInset(), configurable: true });
Object.assign(canvas, {
  clientLeft: 0, clientTop: 0,
  getBoundingClientRect: () => ({ left: 100, top: 50, width: 896, height: 386 }),
  appendChild: () => { probes += 1; },
  ownerDocument: {
    createElement: () => {
      const probe = {
        style: { cssText: '' },
        getBoundingClientRect: () => ({
          width: probe.style.cssText.includes('width:var(--preview-editor-reserve,0px)') ? editorReserve : 0
        }),
        remove: () => { probes -= 1; }
      };
      return probe;
    }
  }
});
const expectFitted = async (zoom) => {
  await Vue.nextTick();
  const box = visualBox();
  const right = 996 - editorReserve;
  const margins = [box.left - 100, right - box.right, box.top - 50, 436 - box.bottom];
  assert.equal(uiState.zoom.value, zoom);
  assert.equal(margins.every((margin) => margin >= 8 - 1e-9), true);
  assert.equal(Math.abs(margins[0] - margins[1]) < 1e-9, true);
  assert.equal(Math.abs(margins[2] - margins[3]) < 1e-9, true);
  assert.equal(probes, 0);
  assert.equal(wrapper.style.transform, `translate(${uiState.canvasPan.x}px, ${uiState.canvasPan.y}px) scale(${zoom})`);
  assert.equal(wrapper.style.transition, 'transform 0.2s');
};
uiState.zoom.value = 2;
Object.assign(uiState.canvasPan, { x: 30, y: -20 });
panZoom.handleWheel(wheelEvent);
panZoom.fitPreviewToViewport();
// 370 / 817 = 0.4529 rounds down to the whole percent, so the diagram fills the frame.
await expectFitted(0.45);
assert.equal(panZoom.isPreviewTransformActive(), false);
assert.equal(timers.size, 0);
// OV-3697: the Fit reads the layout, not the single-precision box of the pan and zoom it starts from, so a
// second Fit, and a Fit from another view, give the same view to the last bit.
const fittedPan = { ...uiState.canvasPan };
panZoom.fitPreviewToViewport();
await expectFitted(0.45);
assert.equal(JSON.stringify(uiState.canvasPan), JSON.stringify(fittedPan));
for (const [zoom, x, y] of [[0.3, -754.671875, 0.83625], [2.7, 123.456, -98.7]]) {
  uiState.zoom.value = zoom;
  Object.assign(uiState.canvasPan, { x, y });
  await Vue.nextTick();
  panZoom.fitPreviewToViewport();
  await expectFitted(0.45);
  assert.equal(JSON.stringify(uiState.canvasPan), JSON.stringify(fittedPan), `a Fit from ${zoom} at (${x}, ${y})`);
}
// The frame is measured in the same layout read: a 300x200 figure zoomed to 3 overflows the canvas and its
// scrollbars take 8 px off the client box, but the Fit reads the room the fitted view has, (386 - 16) / 200.
scrollbars = true;
Object.assign(layout, { width: 300, height: 200 });
uiState.zoom.value = 3;
Object.assign(uiState.canvasPan, { x: 0, y: 0 });
await Vue.nextTick();
assert.equal(canvas.clientWidth, 888);
panZoom.fitPreviewToViewport();
await expectFitted(1.85);
const fittedSmallPan = { ...uiState.canvasPan };
panZoom.fitPreviewToViewport();
await expectFitted(1.85);
assert.equal(JSON.stringify(uiState.canvasPan), JSON.stringify(fittedSmallPan));
scrollbars = false;
Object.assign(layout, { width: 997, height: 817 });
// An open Editor drawer reserves 360 px: Fit uses the 536 px left of it,
// (536 - 16) / 997 = 0.5215 for a wide diagram.
editorReserve = 360;
Object.assign(layout, { width: 997, height: 300 });
panZoom.fitPreviewToViewport();
await expectFitted(0.52);
editorReserve = 0;
Object.assign(layout, { width: 50, height: 40 });
panZoom.fitPreviewToViewport();
await expectFitted(5);
// A 40x40 control in the middle of the right edge, from x 956. The fitted 50x40 diagram (250x200 px,
// centred at x 548) stays clear of it and keeps its Fit.
edgeControl = { left: 956, top: 223, width: 40, height: 40 };
panZoom.fitPreviewToViewport();
await expectFitted(5);
// A wide diagram fitted to the whole frame, 0.88 of 997 px, would end at x 987 over the control: it is
// fitted to the 856 px left of the control instead, (856 - 16) / 997 = 0.8425, with the same margins.
Object.assign(layout, { width: 997, height: 300 });
panZoom.fitPreviewToViewport();
assert.equal(uiState.zoom.value, 0.84);
const besideControl = visualBox();
assert.equal(Math.abs((besideControl.left - 100) - (956 - besideControl.right)) < 1e-9, true);
assert.equal(besideControl.left - 100 >= 8, true);
assert.equal(Math.abs((besideControl.top - 50) - (436 - besideControl.bottom)) < 1e-9, true);
// Only the copies were measured, each one removed, and the control keeps its transition.
assert.equal(JSON.stringify(controlCopies.map(({ deep, placed }) => ({ deep, placed }))),
  JSON.stringify([{ deep: true, placed: false }, { deep: true, placed: false }]));
assert.equal(control.style.transition, 'transform 0.3s');
edgeControl = { left: 0, top: 0, width: 0, height: 0 };
panZoom.fitPreviewToViewport();
await expectFitted(0.88);
edgeControl = null;
Object.assign(layout, { width: 50, height: 40 });
panZoom.fitPreviewToViewport();
await expectFitted(5);
// The Fit measures the layout under the identity transform with the transition held off, and the browser
// takes the shown view back before the transition returns; the Fit shown at once reads its fitted view with
// no transition before the transition comes back.
const transitionsSeen = [];
const svgBox = fakeSvg.getBoundingClientRect;
fakeSvg.getBoundingClientRect = () => {
  transitionsSeen.push(`measure ${wrapper.style.transform} ${wrapper.style.transition}`);
  return svgBox();
};
globalThis.getComputedStyle = (element) => {
  if (element !== canvas) transitionsSeen.push(`read ${element.style.transform} ${element.style.transition}`);
  return computedStyle(element);
};
Object.assign(layout, { width: 997, height: 817 });
uiState.zoom.value = 2;
Object.assign(uiState.canvasPan, { x: 30, y: -20 });
panZoom.showHiddenPreview({ fit: true });
const fittedTransform = `translate(${uiState.canvasPan.x}px, ${uiState.canvasPan.y}px) scale(0.45)`;
assert.equal(transitionsSeen.join('\n'), [
  'measure translate(0px, 0px) scale(1) none', 'read translate(30px, -20px) scale(2) none', `read ${fittedTransform} none`
].join('\n'));
await expectFitted(0.45);
assert.equal(uiState.previewTransitionHeld.value, false);
fakeSvg.getBoundingClientRect = svgBox;
globalThis.getComputedStyle = computedStyle;
completeCase('fit centres the whole diagram left of the Editor reserve at the largest whole percent and cancels the interaction');

// D-57 (OV-3620, OV-3621): the figure a drawing switch shows keeps the zoom, its top stays where pan y 0
// puts it, and it is centred in the free room Fit measures (left of the Editor reserve, and left of the edge
// control when centred in the whole room it would meet it). It is hidden until placed, placed at once.
/** @param {number} right The free room's right edge. */
const expectCentred = (zoom, right) => {
  const box = visualBox();
  assert.equal(uiState.zoom.value, zoom);
  assert.equal(uiState.canvasPan.y, 0);
  // Centred to the whole pixel: the pan is a whole pixel.
  assert.equal(Number.isInteger(uiState.canvasPan.x), true);
  assert.equal(Math.abs((box.left - 100) - (right - box.right)) <= 1 + 1e-9, true);
  assert.equal(uiState.previewHiddenUntilShown.value, false);
  assert.equal(wrapper.style.visibility, '');
  assert.equal(wrapper.style.transition, 'transform 0.2s');
};
const placements = [];
globalThis.getComputedStyle = (element) => {
  if (element !== canvas) placements.push(`read ${element.style.transform} ${element.style.transition} ${element.style.visibility}`);
  return computedStyle(element);
};
Object.assign(layout, { width: 997, height: 300 });
uiState.zoom.value = 0.5;
Object.assign(uiState.canvasPan, { x: 30, y: -20 });
const arriving = panZoom.showArrivingFigure();
assert.equal(uiState.previewHiddenUntilShown.value, true);
await arriving;
expectCentred(0.5, 996);
// The layout read takes the departing view back with no transition; the placed view then shows at once.
assert.equal(placements.join('\n'), [
  'read translate(30px, -20px) scale(0.5) none hidden',
  `read translate(${uiState.canvasPan.x}px, 0px) scale(0.5) none hidden`
].join('\n'));
// An open Editor drawer: the room ends at its reserve.
editorReserve = 360;
await panZoom.showArrivingFigure();
expectCentred(0.5, 636);
editorReserve = 0;
// A figure centred in the whole room at 0.84 would end at x 967 over the control: it is centred left of it.
edgeControl = { left: 956, top: 223, width: 40, height: 40 };
uiState.zoom.value = 0.84;
await panZoom.showArrivingFigure();
expectCentred(0.84, 956);
// Clear of the control, a small figure stays centred in the whole room, at its own zoom.
Object.assign(layout, { width: 50, height: 40 });
uiState.zoom.value = 1;
await panZoom.showArrivingFigure();
expectCentred(1, 996);
edgeControl = null;
// A figure wider than the room at the kept zoom shows its start where the canvas lays out its content, inside
// its padding, as Generate shows it: 8 px in, and 10 px in with a 20 px browser font (p-2 is 0.5rem).
Object.assign(layout, { width: 997, height: 300 });
uiState.zoom.value = 2;
await panZoom.showArrivingFigure();
assert.equal(uiState.zoom.value, 2);
assert.equal(uiState.canvasPan.y, 0);
assert.equal(Math.abs(visualBox().left - 108) <= 0.5, true);
canvasPaddingLeft = '10px';
await panZoom.showArrivingFigure();
assert.equal(Math.abs(visualBox().left - 110) <= 0.5, true, `the start at the canvas padding: ${visualBox().left}`);
canvasPaddingLeft = '8px';
// N-c: a placement that throws never leaves the figure hidden.
const canvasBox = canvas.getBoundingClientRect;
canvas.getBoundingClientRect = () => { throw new Error('no box'); };
await assert.rejects(panZoom.showArrivingFigure(), /no box/);
assert.equal(uiState.previewHiddenUntilShown.value, false);
canvas.getBoundingClientRect = canvasBox;
completeCase('a drawing switch keeps the zoom and centres the arriving figure in the free room, shown at once');

wrapper.querySelector = () => null;
panZoom.fitPreviewToViewport();
assert.equal(uiState.zoom.value, 2, 'with no figure the Fit changes nothing');
for (const key of ['querySelector', 'getBoundingClientRect']) delete wrapper[key];
for (const key of ['clientLeft', 'clientTop', 'clientWidth', 'clientHeight', 'getBoundingClientRect', 'appendChild', 'parentElement', 'ownerDocument']) delete canvas[key];
delete globalThis.getComputedStyle;

panZoom.handleWheel(wheelEvent);
panZoom.disposePanZoom();
assert.equal(panZoom.isPreviewTransformActive(), false);
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
    if (this.kind === 'match') return selector.includes('data-gbdraw-match-id');
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
    if (selector.includes('data-gbdraw-match-id')) {
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
  'data-gbdraw-match-id': 'm1',
  'data-orthogroup-id': 'OG1'
});
const unrelatedPath = new FakeElement('generic', {
  id: 'axis-path',
  'data-orthogroup-id': 'OG2'
});
const svg = new FakeSvg([featureOne, featureTwo, match, unrelatedPath]);
globalThis.document = {
  elementsFromPoint: () => [featureOne],
  querySelector: () => null,
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
  previewHiddenUntilShown: ref(false),
  previewTransitionHeld: ref(false),
  canvasContainerRef: ref(featureInteractionCanvas),
  svgContainer: ref(featureInteractionWrapper)
};
const featurePanZoom = createPanZoom(featureInteractionState);

const featureMap = new Map([
  ['f1', { id: 'f1', rendered_svg_id: 'f1', label: 'One', orthogroupId: 'OG1', start: 0, end: 10 }],
  ['f2', { id: 'f2', rendered_svg_id: 'f2', label: 'Two', orthogroupId: 'OG1', start: 10, end: 20 }]
]);
let selectedFeatureId = '';
// The selection owner's ports: a modifier click records the feature it selects.
const selectionPorts = {
  consumeSuppressNextClick: () => false,
  selectFeatureByClick(featureId) {
    selectedFeatureId = featureId;
  },
  markPlainFeatureClick() {},
  clearFeatureSelection() {},
  startMarqueePointer: () => false,
  moveMarqueePointer: () => false,
  commitMarqueePointer: () => false,
  cancelMarqueePointer() {}
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
  selectedAnnotation: ref(null),
  featurePopupSize: { width: 0, height: 0 },
  featureSelectionDrag: { active: false },
  skipCaptureBaseConfig: ref(false),
  richFeaturePopup: ref(false)
};
const { createFeatureSvgActions, HOVER_HIGHLIGHT_DWELL_MS } = await import(
  pathToFileURL(join(tempDir, 'app', 'feature-editor', 'svg-actions.js'))
);
const featureActions = createFeatureSvgActions({
  state: withDrawings(state),
  getFeatureColor: () => '#123456',
  getEffectiveLegendCaption: () => '',
  ...selectionPorts,
  isPreviewTransformActive: featurePanZoom.isPreviewTransformActive,
  subscribePreviewTransformChanges: featurePanZoom.subscribePreviewTransformChanges
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
// The highlight is written once the pointer has rested (OV-900).
assert.equal(featureOne.style.opacity, '');
runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
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
assert.equal(featurePanZoom.isPreviewTransformActive(), true);
runTimersAtDelay(220);
featureInteractionWrapper.dispatch('transitionend', { propertyName: 'transform' });
assert.equal(featurePanZoom.isPreviewTransformActive(), false);
flushFrames();
runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
assert.equal(featureOne.style.opacity, '0.7');
assert.equal(metrics.filter((metric) => metric.name === 'previewTransformHoverReconcileCount').length, 1);
assert.equal(metrics.filter((metric) => metric.name === 'featureDomFullScanCount').length, 1);
svg.dispatch('mouseout', eventAtFeature);
assert.equal(featureOne.style.opacity, '');
completeCase('ordinary Feature hover resumes after the final source settles');

// OV-2281: `elementsFromPoint` here finds Feature One at any point, as it
// finds the drawing under a popup or toolbar the pointer moved onto.
const eventOffSvg = { ...eventAtFeature, relatedTarget: null, clientX: 2000 };
featurePanZoom.startPan({ ...eventAtFeature, button: 0, shiftKey: false, target: backgroundTarget });
featurePanZoom.endPan(eventAtFeature);
assert.equal(frames.size, 1);
svg.dispatch('mouseout', eventOffSvg);
flushFrames();
runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
assert.equal(featureOne.style.opacity, '');
featurePanZoom.startPan({ ...eventAtFeature, button: 0, shiftKey: false, target: backgroundTarget });
svg.dispatch('mouseout', eventOffSvg);
featurePanZoom.endPan(eventOffSvg);
flushFrames();
runTimersAtDelay(HOVER_HIGHLIGHT_DWELL_MS);
assert.equal(featureOne.style.opacity, '');
completeCase('a pointer that leaves the svg during a pan or before its hover replay leaves no hover');

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
  join(repoRoot, 'gbdraw', 'web', 'js', 'services', 'feature-dom.js'),
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
const { FEATURE_SELECTOR } = await import(pathToFileURL(join(tempDir, 'services', 'feature-dom.js')));
assert.equal(cssFeatureSelector, FEATURE_SELECTOR);
for (const excludedSelector of [
  'data-gbdraw-match-id',
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
  /const clearTrackedHoverStyles = \(\) => \{[\s\S]*?const clearHover = \(\) => \{[\s\S]*?\n\s*\};/
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

// OV-912: the surface's `v-preview-transform` (components.js) binds the projection of the view, which writes the transition.
const componentsSource = await readFile(join(repoRoot, 'gbdraw', 'web', 'js', 'components.js'), 'utf8');
assert.match(indexHtml, /ref="svgContainer"[^>]*\sv-preview-transform[\s>]/);
assert.match(componentsSource, /projectPreviewSurface\(element, state\)/);
assert.match(indexHtml, /@pointercancel="endPan"/);
assert.match(indexHtml, /@lostpointercapture="endPan"/);
completeCase('template preserves transition and routes interrupted pan cleanup');

delete globalThis.__GBDRAW_TEST_HOOKS__;
console.log(JSON.stringify({
  status: 'preview transform interaction tests passed',
  assertionCount,
  cases: caseNames
}, null, 2));
