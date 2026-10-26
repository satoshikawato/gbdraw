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
