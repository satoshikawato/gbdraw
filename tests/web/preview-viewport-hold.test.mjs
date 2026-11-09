import assert from 'node:assert/strict';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';
import { callerSession } from './helpers/session-caller.mjs';
import { savedDrawingView } from './helpers/saved-drawing-view.mjs';
import { sessionRollbackPort } from './helpers/session-rollback-port.mjs';

// While Review alignment is open the figure is fitted, and what is saved or captured is the view
// from before the review. The pan/zoom owner holds it; the Session writer and the UI capture with
// navigation read it through the reader the composition root registers.
globalThis.window = {
  Vue: {
    ref: value => ({ value }), reactive: value => value,
    computed: getter => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: value => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null })
};
const sessionConfig = await import('../../gbdraw/web/js/services/config.js');
const { buildUiStateData, setHeldPreviewViewportReader } = sessionConfig;
const { state } = await import('../../gbdraw/web/js/state.js');
const { exportSession } = callerSession(sessionConfig, state, await sessionRollbackPort(state));
const { createPanZoom } = await import('../../gbdraw/web/js/app/ui.js');

const panZoom = createPanZoom(state);
const setLive = (zoom, x, y) => { state.zoom.value = zoom; state.canvasPan.x = x; state.canvasPan.y = y; };
const live = () => ({ zoom: state.zoom.value, canvasPan: { x: state.canvasPan.x, y: state.canvasPan.y } });
const before = { zoom: 0.46, canvasPan: { x: 12, y: -7 } };
const fitted = { zoom: 0.31, canvasPan: { x: -40, y: 3 } };

test('the pan/zoom owner holds the view from before the review and restores or drops it', () => {
  setLive(0.46, 12, -7);
  assert.equal(panZoom.heldPreviewViewport(), null);
  panZoom.holdPreviewViewport(() => true);
  setLive(0.31, -40, 3);
  assert.deepEqual(panZoom.heldPreviewViewport(), before);
  // A second hold keeps the first view.
  panZoom.holdPreviewViewport(() => true);
  assert.deepEqual(panZoom.heldPreviewViewport(), before);

  panZoom.releasePreviewViewport({ restore: true });
  assert.deepEqual(live(), before);
  assert.equal(panZoom.heldPreviewViewport(), null);

  panZoom.holdPreviewViewport(() => true);
  setLive(0.31, -40, 3);
  panZoom.releasePreviewViewport({ restore: false });
  assert.deepEqual(live(), fitted);
  assert.equal(panZoom.heldPreviewViewport(), null);
  // Without a hold, a release changes nothing.
  panZoom.releasePreviewViewport({ restore: true });
  assert.deepEqual(live(), fitted);
});

test('the UI capture with navigation reads the held view, and the capture without it has none', () => {
  setHeldPreviewViewportReader(panZoom.heldPreviewViewport);
  setLive(0.46, 12, -7);
  panZoom.holdPreviewViewport(() => true);
  setLive(0.31, -40, 3);
  const ui = buildUiStateData(state.activeDrawing(), { includePreviewNavigation: true });
  assert.deepEqual({ zoom: ui.zoom, canvasPan: ui.canvasPan }, before);
  const intent = buildUiStateData(state.activeDrawing(), { includePreviewNavigation: false });
  assert.equal('zoom' in intent, false);
  assert.equal('canvasPan' in intent, false);
  panZoom.releasePreviewViewport({ restore: false });
  const afterDrop = buildUiStateData(state.activeDrawing(), { includePreviewNavigation: true });
  assert.deepEqual({ zoom: afterDrop.zoom, canvasPan: afterDrop.canvasPan }, fitted);
});

test('a Session saved during the hold stores the view from before it', async () => {
  const savedUi = async (title) => {
    const result = await exportSession(title);
    assert.equal(result.status, 'saved');
    // The view is the shown drawing's (Session 49 keeps it in that drawing's `ui`).
    const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
    return savedDrawingView(saved, saved.ui.activeDrawingId).ui;
  };
  setHeldPreviewViewportReader(panZoom.heldPreviewViewport);
  setLive(0.46, 12, -7);
  panZoom.holdPreviewViewport(() => true);
  setLive(0.31, -40, 3);
  const held = await savedUi('during review');
  assert.deepEqual({ zoom: held.zoom, canvasPan: held.canvasPan }, before);
  panZoom.releasePreviewViewport({ restore: true });
  const after = await savedUi('after review');
  assert.deepEqual({ zoom: after.zoom, canvasPan: after.canvasPan }, before);

  // Without a registered reader the writers read the live view, as before.
  setHeldPreviewViewportReader(null);
  setLive(0.31, -40, 3);
  const unregistered = await savedUi('no reader');
  assert.deepEqual({ zoom: unregistered.zoom, canvasPan: unregistered.canvasPan }, fitted);
});

test('the held view refits the shown figure when its room changes, but not in the room of the last Fit', () => {
  const observers = [];
  globalThis.ResizeObserver = class {
    constructor(callback) { this.callback = callback; this.observed = null; observers.push(this); }
    observe(element) { this.observed = element; }
    disconnect() { this.observed = null; }
  };
  const canvas = {
    offsetWidth: 800, offsetHeight: 400, clientWidth: 800, clientHeight: 400, clientLeft: 0, clientTop: 0,
    getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 400 }), style: {}
  };
  let fits = 0;
  const svg = { getBoundingClientRect: () => { fits += 1; return { left: 300, top: 0, width: 200, height: 100 }; } };
  state.canvasContainerRef.value = canvas;
  state.svgContainer.value = { querySelector: () => svg, getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 100 }) };
  const owner = createPanZoom(state);
  try {
    let shown = true;
    const resize = (width) => {
      canvas.offsetWidth = width;
      observers.filter(({ observed }) => observed === canvas).forEach(({ callback }) => callback());
    };
    owner.holdPreviewViewport(() => shown);
    resize(800);
    assert.equal(fits, 1, 'the first room under the hold is fitted');
    resize(800);
    assert.equal(fits, 1, 'the room of the last Fit is no change');
    resize(600);
    assert.equal(fits, 2, 'a new room is fitted');
    canvas.offsetWidth = 650;
    owner.fitPreviewToViewport();
    resize(650);
    assert.equal(fits, 3, 'the room of a Fit made under the hold is no change either');
    resize(0);
    assert.equal(fits, 3, 'a hidden canvas has no room');
    shown = false;
    resize(700);
    assert.equal(fits, 3, 'another figure shown under the hold is not fitted');
    owner.releasePreviewViewport({ restore: false });
    shown = true;
    resize(500);
    assert.equal(fits, 3, 'a released hold observes no room');
    owner.fitPreviewToViewport();
    owner.holdPreviewViewport(() => shown);
    resize(500);
    assert.equal(fits, 5, 'a new hold forgets the room of the Fit before it');
  } finally {
    owner.disposePanZoom();
    state.canvasContainerRef.value = null;
    state.svgContainer.value = null;
    delete globalThis.ResizeObserver;
  }
});
