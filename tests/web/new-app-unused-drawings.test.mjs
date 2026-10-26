// A new app starts with one drawing per mode (state.js). Setup leaves each at
// its mode's new-drawing defaults, the drawing `withNewDrawingDefaults`
// (services/reset.js) builds and `tools/generate_mode_slice_defaults.mjs`
// reads, so Save does not write a drawing the person never touched (R16-8,
// PD-OI-090 〔Save と Load〕, OV-1901).
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { gunzipSync } from 'node:zlib';
import test from 'node:test';
import { installFakeSvgDom } from './fake-svg-dom.mjs';
import { readFileSync } from 'node:fs';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }), shallowRef: (value) => ({ value }), reactive: (value) => value,
    isRef: (value) => Boolean(value && typeof value === 'object' && 'value' in value),
    computed: (getter) => ({ get value() { return getter(); } }), nextTick: async () => {},
    watch: () => () => {}, onMounted: () => {}, onUnmounted: () => {}, provide: () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {
  body: { appendChild: () => {} },
  createElement: () => ({ addEventListener: () => {}, click: () => {}, parentNode: null }),
  querySelector: () => null
};
installFakeSvgDom();
globalThis.File = File;

const { state } = await import('../../gbdraw/web/js/state.js');
const { createAppSetup } = await import('../../gbdraw/web/js/app/app-setup.js');
const { exportSession } = await import('../../gbdraw/web/js/services/config.js');
const { createPaletteLoader } = await import('../../gbdraw/web/js/app/palettes.js');

test("a new app's untouched Circular drawing is not saved", async () => {
  const app = createAppSetup();
  // The boot palette install, as app.js runs it after setup.
  const palettes = JSON.parse(readFileSync(new URL('../../gbdraw/web/gallery/palettes/palettes.json', import.meta.url), 'utf8'));
  assert.ok(createPaletteLoader({ state }).applyPalettes(palettes.palettes));
  app.selectDrawing('linear');
  state.activeDrawing().form.plot_title = 'Linear only';
  const result = await exportSession('ov-1901', { shownDrawing: state.activeDrawing });
  assert.equal(result.status, 'saved', result.error?.message);
  const saved = JSON.parse(gunzipSync(Buffer.from(await result.blob.arrayBuffer())));
  assert.deepEqual(saved.drawings.map((entry) => entry.id), ['linear']);
});
