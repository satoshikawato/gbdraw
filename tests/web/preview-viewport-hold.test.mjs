import assert from 'node:assert/strict';
import { test } from 'node:test';
import { gunzipSync } from 'node:zlib';

// D-11: while Review alignment is open the figure is fitted, and what is saved or captured is the
// view from before the review. The pan/zoom owner holds it; the Session writer and the UI capture
// with navigation read it through the one reader the composition root registers.
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
const { exportSession, buildUiStateData, setPreviewViewportReader } = await import('../../gbdraw/web/js/services/config.js');
const { state } = await import('../../gbdraw/web/js/state.js');
const { createPanZoom } = await import('../../gbdraw/web/js/app/ui.js');

const panZoom = createPanZoom(state);
const setLive = (zoom, x, y) => { state.zoom.value = zoom; state.canvasPan.x = x; state.canvasPan.y = y; };
const live = () => ({ zoom: state.zoom.value, canvasPan: { x: state.canvasPan.x, y: state.canvasPan.y } });
const before = { zoom: 0.46, canvasPan: { x: 12, y: -7 } };
const fitted = { zoom: 0.31, canvasPan: { x: -40, y: 3 } };

test('the pan/zoom owner holds the view from before the review and restores or drops it', () => {
  setLive(0.46, 12, -7);
  assert.deepEqual(panZoom.previewViewportToPersist(), before);
  panZoom.holdPreviewViewport();
  setLive(0.31, -40, 3);
  assert.deepEqual(panZoom.previewViewportToPersist(), before);
  // A second hold keeps the first view.
  panZoom.holdPreviewViewport();
  assert.deepEqual(panZoom.previewViewportToPersist(), before);

  panZoom.releasePreviewViewport({ restore: true });
  assert.deepEqual(live(), before);
  assert.deepEqual(panZoom.previewViewportToPersist(), before);

  panZoom.holdPreviewViewport();
  setLive(0.31, -40, 3);
  panZoom.releasePreviewViewport({ restore: false });
  assert.deepEqual(live(), fitted);
  assert.deepEqual(panZoom.previewViewportToPersist(), fitted);
  // Without a hold, a release changes nothing.
  panZoom.releasePreviewViewport({ restore: true });
  assert.deepEqual(live(), fitted);
});

test('the UI capture with navigation reads the held view, and the capture without it has none', () => {
  setPreviewViewportReader(panZoom.previewViewportToPersist);
  setLive(0.46, 12, -7);
  panZoom.holdPreviewViewport();
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
    return JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer()))).ui;
  };
  setPreviewViewportReader(panZoom.previewViewportToPersist);
  setLive(0.46, 12, -7);
  panZoom.holdPreviewViewport();
  setLive(0.31, -40, 3);
  const held = await savedUi('during review');
  assert.deepEqual({ zoom: held.zoom, canvasPan: held.canvasPan }, before);
  panZoom.releasePreviewViewport({ restore: true });
  const after = await savedUi('after review');
  assert.deepEqual({ zoom: after.zoom, canvasPan: after.canvasPan }, before);

  // Without a registered reader the writers read the live view, as before.
  setPreviewViewportReader(null);
  setLive(0.31, -40, 3);
  const unregistered = await savedUi('no reader');
  assert.deepEqual({ zoom: unregistered.zoom, canvasPan: unregistered.canvasPan }, fitted);
});
