// OV-912: a wheel tick or a zoom button changes the preview zoom. The app
// template does not read the zoom, the pan, or the panning flag, so these
// render only the zoom label and write only the surface transform; History
// and the Session still read `state.zoom` and `state.canvasPan`.
import assert from 'node:assert/strict';
import { readdirSync, readFileSync } from 'node:fs';
import test from 'node:test';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue, setTimeout, clearTimeout };
/** @type {Map<number, FrameRequestCallback>} */
const frames = new Map();
let nextFrame = 1;
globalThis.requestAnimationFrame = (callback) => { frames.set(nextFrame, callback); return nextFrame++; };
globalThis.cancelAnimationFrame = (id) => { frames.delete(id); };
const flushFrames = () => { const pending = [...frames.values()]; frames.clear(); pending.forEach((callback) => callback(0)); };
const { nextTick } = Vue;
const { state } = await import('../../gbdraw/web/js/state.js');
const { previewTransform, PreviewZoomPercent, resetLayoutState } = await import('../../gbdraw/web/js/components.js');
const { createPanZoom } = await import('../../gbdraw/web/js/app/ui.js');

test('the app template reads no preview zoom, pan, panning flag, or Reset layout availability (OV-912, OV-2723)', () => {
  const html = readFileSync(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const surface = html.slice(html.lastIndexOf('<div', html.indexOf('ref="svgContainer"')), html.indexOf('>', html.indexOf('ref="svgContainer"')) + 1);
  assert.match(surface, /\sv-preview-transform[\s>]/);
  assert.doesNotMatch(surface, /:style=/);
  // Outside event handlers, no binding reads them.
  const bindings = [...html.matchAll(/(?:\s:[\w.-]+|\sv-(?:if|show|for|bind|html|text|model)(?::[\w.-]+)?)="([^"]*)"|\{\{([^}]*)\}\}/g)]
    // String literals in an expression are text, not reads.
    .map((match) => (match[1] ?? match[2]).replace(/'[^']*'|`[^`$]*`/g, "''"));
  // Reset layout's availability changes at the first pan frame or wheel tick: the button reads it in a function
  // that `v-reset-layout-state` calls, which the app template creates but does not call (OV-2723).
  for (const name of ['zoom', 'canvasPan', 'isPanning', 'layoutResetAvailable']) {
    const readers = bindings.filter((expression) => new RegExp(`(^|[^\\w.])${name}\\b`).test(expression));
    assert.deepEqual(readers, [], `${name} is read by ${readers.join(' | ')}`);
  }
  const resetLayoutReads = [...html.matchAll(/(?<![\w.])layoutResetAvailable\b/g)].length;
  const deferred = html.match(/\sv-reset-layout-state="\(\) => \(\{[^"]*"/g) ?? [];
  assert.equal(deferred.length, 1);
  assert.equal(deferred[0].match(/layoutResetAvailable\b/g)?.length, resetLayoutReads);
  assert.match(html, /<preview-zoom-percent><\/preview-zoom-percent>/);
});

test('the surface follows the zoom, the pan, panning, and a held transition until it unmounts (OV-912)', async () => {
  const element = { style: {} };
  state.zoom.value = 1;
  state.canvasPan.x = 0;
  state.canvasPan.y = 0;
  state.isPanning.value = false;
  previewTransform.beforeMount(element);
  assert.deepEqual({ ...element.style }, {
    transform: 'translate(0px, 0px) scale(1)', transformOrigin: 'top center', transition: 'transform 0.2s', willChange: '',
    visibility: ''
  });
  // A change of the view shows at Vue's flush, once for all the members it changed.
  state.zoom.value = 1.5;
  state.canvasPan.x = 12;
  await nextTick();
  assert.equal(element.style.transform, 'translate(12px, 0px) scale(1.5)');
  state.isPanning.value = true;
  await nextTick();
  assert.deepEqual([element.style.transition, element.style.willChange], ['none', 'transform']);
  state.isPanning.value = false;
  await nextTick();
  assert.deepEqual([element.style.transition, element.style.willChange], ['transform 0.2s', '']);
  // A held transition shows the view at once, with no transition, so the Fit reads the style its write set.
  state.zoom.value = 2;
  state.previewTransitionHeld.value = true;
  assert.deepEqual([element.style.transform, element.style.transition, element.style.willChange],
    ['translate(12px, 0px) scale(2)', 'none', 'transform']);
  state.zoom.value = 1.5;
  state.previewTransitionHeld.value = false;
  assert.deepEqual([element.style.transform, element.style.transition], ['translate(12px, 0px) scale(1.5)', 'transform 0.2s']);
  // OV-2702: the figure a file Session Load mounted stays hidden until the Load shows it.
  state.previewHiddenUntilShown.value = true;
  await nextTick();
  assert.equal(element.style.visibility, 'hidden');
  state.previewHiddenUntilShown.value = false;
  await nextTick();
  assert.equal(element.style.visibility, '');
  previewTransform.unmounted(element);
  state.zoom.value = 2;
  await nextTick();
  assert.equal(element.style.transform, 'translate(12px, 0px) scale(1.5)', 'an unmounted surface is not written');
  state.zoom.value = 1;
  state.canvasPan.x = 0;
  assert.equal(PreviewZoomPercent.template.includes('zoom'), true);
  assert.equal(PreviewZoomPercent.setup().zoom, state.zoom);
});

// The surface's style has one writer: the projection of the view in app/ui.js, which
// `v-preview-transform` binds. No other module writes the surface's style or knows its format.
// It reads every module for a style write through the surface's ref (`svgContainer`) or a
// variable assigned from it (`const surface = svgContainer.value`), and for the transform format.
test('only the projection of the view writes the Preview surface style', () => {
  const root = new URL('../../gbdraw/web/js/', import.meta.url);
  const modules = readdirSync(root, { recursive: true })
    .filter((path) => path.endsWith('.js') && !path.startsWith('vendor'));
  const styleWrite = String.raw`(?:\.value)?\??\.style(?:\.(?:transform|transformOrigin|transition|willChange|visibility|cssText)\s*=(?!=)|\.setProperty\(|\s*=(?!=))`;
  const writers = [];
  for (const path of modules) {
    const source = readFileSync(new URL(path, root), 'utf8');
    const aliases = [...source.matchAll(/\b(\w+)\s*=\s*(?:\/\*\*.*?\*\/\s*\()?(?:state\.)?svgContainer\.value\b/g)].map((match) => match[1]);
    const names = ['svgContainer', ...aliases].join('|');
    if (new RegExp(String.raw`\b(?:${names})${styleWrite}`).test(source)) writers.push(`${path} (${names})`);
    if (/translate\(\$\{[^}]+\}px, \$\{[^}]+\}px\) scale\(/.test(source)) writers.push(`${path} (transform format)`);
  }
  assert.deepEqual(writers, ['app/ui.js (transform format)']);
  // In app/ui.js the five style members are written once each, by the projection.
  const owner = readFileSync(new URL('app/ui.js', root), 'utf8');
  const projection = owner.slice(owner.indexOf('export const projectPreviewSurface'), owner.indexOf('\n};\n', owner.indexOf('export const projectPreviewSurface')));
  const memberWrites = /\.style\.(?:transform|transformOrigin|transition|willChange|visibility)\s*=(?!=)/g;
  assert.equal(owner.match(memberWrites)?.length, 5);
  assert.equal(projection.match(memberWrites)?.length, 5);
});

test('a wheel zoom during a pan keeps the figure where the pan holds it (OV-2722)', async () => {
  const element = { style: {}, addEventListener() {}, removeEventListener() {} };
  const canvas = { style: {} };
  state.svgContainer.value = element;
  state.canvasContainerRef.value = canvas;
  const panZoom = createPanZoom(state);
  previewTransform.beforeMount(element);
  const background = { tagName: 'rect', closest: () => null };
  panZoom.startPan({ button: 0, shiftKey: false, clientX: 100, clientY: 120, target: background, preventDefault() {} });
  // A browser runs the microtasks queued by the press before the next animation frame.
  await nextTick();
  panZoom.doPan({ clientX: 160, clientY: 200 });
  flushFrames();
  await nextTick();
  assert.equal(element.style.transform, 'translate(60px, 80px) scale(1)');
  panZoom.handleWheel({ deltaY: -100, clientX: 160, clientY: 200 });
  await nextTick();
  assert.equal(element.style.transform, 'translate(60px, 80px) scale(1.1)');
  assert.equal(element.style.transition, 'none');
  // Every reader of the view sees the live pan, not the pan's start.
  assert.deepEqual({ ...state.canvasPan }, { x: 60, y: 80 });
  panZoom.endPan({ clientX: 160, clientY: 200 });
  await nextTick();
  assert.equal(element.style.transform, 'translate(60px, 80px) scale(1.1)');
  assert.equal(element.style.transition, 'transform 0.2s');
  panZoom.disposePanZoom();
  previewTransform.unmounted(element);
  state.svgContainer.value = null;
  state.canvasContainerRef.value = null;
  state.zoom.value = 1;
  Object.assign(state.canvasPan, { x: 0, y: 0 });
});

test('the title and disabled state of the Reset layout button follow its availability (OV-2723)', async () => {
  const element = { style: {}, addEventListener() {}, removeEventListener() {} };
  state.svgContainer.value = element;
  state.canvasContainerRef.value = { style: {} };
  const panZoom = createPanZoom(state);
  // The view part of app-setup.js's `layoutResetAvailable`, which the button's function reads.
  const available = Vue.computed(() => Boolean(Number(state.canvasPan.x) || Number(state.canvasPan.y)) || state.zoom.value !== 1);
  const button = { title: '', disabled: false };
  let calls = 0;
  resetLayoutState.beforeMount(button, { value: () => {
    calls += 1;
    return { title: available.value ? 'Return' : 'Nothing moved yet', disabled: !available.value };
  } });
  assert.deepEqual(button, { title: 'Nothing moved yet', disabled: true });
  const background = { tagName: 'rect', closest: () => null };
  panZoom.startPan({ button: 0, shiftKey: false, clientX: 100, clientY: 120, target: background, preventDefault() {} });
  await nextTick();
  for (const step of [1, 2, 3]) {
    panZoom.doPan({ clientX: 100 + 10 * step, clientY: 120 });
    flushFrames();
    await nextTick();
  }
  assert.deepEqual(button, { title: 'Return', disabled: false });
  panZoom.endPan({ clientX: 130, clientY: 120 });
  panZoom.resetPreviewViewport({ resetZoom: true });
  await nextTick();
  assert.deepEqual(button, { title: 'Nothing moved yet', disabled: true });
  panZoom.handleWheel({ deltaY: -100, clientX: 0, clientY: 0 });
  await nextTick();
  assert.deepEqual(button, { title: 'Return', disabled: false });
  const before = calls;
  resetLayoutState.unmounted(button);
  panZoom.resetPreviewViewport({ resetZoom: true });
  await nextTick();
  assert.equal(calls, before, 'an unmounted button is not written');
  panZoom.disposePanZoom();
  state.svgContainer.value = null;
  state.canvasContainerRef.value = null;
});
